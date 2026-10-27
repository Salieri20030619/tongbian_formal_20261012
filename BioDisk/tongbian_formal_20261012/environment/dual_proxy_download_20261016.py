from pathlib import Path
import os, sys, json, datetime, hashlib, subprocess, fcntl, time, signal, traceback
import threading, concurrent.futures

P=Path('/mnt/BioDisk/wangmingyuan/tongbian_formal_20261012')
D=Path('/mnt/BioData/bio/wangmingyuan/tongbian_formal_20261012/raw_reads/PRJNA1232352')
OLD=P/'logs/jobs/PRJNA1232352_direct_dual_20261014_v2'
LOCK=P/'logs/jobs/PRJNA1232352_direct_download_v1/run.lock'
J=P/'logs/jobs/PRJNA1232352_proxy_dual_20261016_v1'
MAN=P/'01_inventory/resume_20261013/PRJNA1232352_download_manifest.json'
PROXY='http://127.0.0.1:7897'
now=lambda:datetime.datetime.now(datetime.timezone.utc).isoformat()

def write(path,obj):
    tmp=path.with_suffix(path.suffix+'.tmp')
    tmp.write_text(json.dumps(obj,indent=2,ensure_ascii=False),encoding='utf-8');os.replace(tmp,path)

entries={e['filename']:e for e in json.loads(MAN.read_text(encoding='utf-8'))}
if '--launch' in sys.argv:
    assert not J.exists(),'Job already exists; inspect before any restart'
    old=json.loads((OLD/'status.json').read_text(encoding='utf-8'))
    assert not Path('/proc',str(old['runner_pid'])).exists(),'Old runner present'
    lock=LOCK.open('a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    oldcfg=json.loads((OLD/'config.json').read_text(encoding='utf-8'))
    done={n:v for n,v in old['files'].items() if v['status'] in ('VERIFIED','REUSED_VERIFIED')}
    for n,v in done.items():assert (D/n).is_file() and (D/n).stat().st_size==v['bytes'],n
    pending=[n for n in entries if n not in done]
    priority=[n for n in oldcfg['priority_failed'] if n in pending]
    priority = [n for n,v in old['files'].items() if v['status']=='INTERRUPTED' and n in pending] + priority
    priority = list(dict.fromkeys(priority))
    priority += [n for n,v in old['files'].items() if v['status']=='FAILED' and n in pending and n not in priority]
    assert len(pending)==len(set(pending)) and set(pending)|set(done)==set(entries)
    cfg={'created_utc':now(),'workers':2,'proxy':PROXY,'direct_fallback':False,
         'user_instruction':'User requested proxy downloads on current intranet and throughput comparison',
         'scheduling':'large lane: previous failures first, then largest; small lane: smallest nonpriority, then failures',
         'manifest_sha256':hashlib.sha256(MAN.read_bytes()).hexdigest(),'pending':pending,
         'priority_failed':priority,'previous_verified':done,'previous_job':str(OLD),
         'validation':'exact size + ENA MD5 + gzip -t before rename; corrupt partials retained',
         'tunnel_dependency':'Windows sgs10-codex SSH tunnel and local proxy 7897 must stay available'}
    J.mkdir();(J/'files').mkdir();write(J/'config.json',cfg)
    lock.close()
    with (J/'launcher.log').open('ab') as log:
        proc=subprocess.Popen(['python3',str(Path(__file__).resolve()),'--run'],stdin=subprocess.DEVNULL,stdout=log,stderr=log,start_new_session=True,cwd=str(P))
    (J/'runner_pid.txt').write_text(str(proc.pid),encoding='utf-8')
    print(json.dumps({'pid':proc.pid,'job_dir':str(J),'workers':2,'pending':len(pending)}));sys.exit(0)

assert '--run' in sys.argv
lock=LOCK.open('a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
cfg=json.loads((J/'config.json').read_text(encoding='utf-8'))
assert hashlib.sha256(MAN.read_bytes()).hexdigest()==cfg['manifest_sha256']
guard=threading.RLock();stopping=threading.Event();children={};pending=list(cfg['pending'])
state={'status':'RUNNING','started_utc':now(),'runner_pid':os.getpid(),'workers':2,'proxy':PROXY,
       'total_files':len(entries),'total_expected_bytes':sum(e['bytes'] for e in entries.values()),
       'files':{n:dict(v,status='REUSED_VERIFIED',verification_source=str(OLD/'status.json')) for n,v in cfg['previous_verified'].items()},'lanes':{}}

def save_locked():
    state.update(updated_utc=now(),completed_files=sum(v['status'] in ('VERIFIED','REUSED_VERIFIED') for v in state['files'].values()),
                 failed_files=sum(v['status']=='FAILED' for v in state['files'].values()),pending_files=len(pending))
    write(J/'status.json',state)

def update(n=None,**kw):
    with guard:
        if n:state['files'].setdefault(n,{}).update(kw)
        else:state.update(kw)
        save_locked()

def take(lane):
    with guard:
        if stopping.is_set() or not pending:return None
        priority=[n for n in cfg['priority_failed'] if n in pending]
        ordinary=[n for n in pending if n not in cfg['priority_failed']]
        if lane=='large':n=priority[0] if priority else max(pending,key=lambda n:entries[n]['bytes'])
        else:n=min(ordinary or pending,key=lambda n:entries[n]['bytes'])
        pending.remove(n);state['lanes'][lane]=n;save_locked();return n

def run_process(lane,cmd,log,n):
    with guard:
        if stopping.is_set():raise InterruptedError('Stop requested')
        child=subprocess.Popen(cmd,stdin=subprocess.DEVNULL,stdout=log,stderr=log)
        children[lane]=child
        state['files'][n]['child_pid']=child.pid;save_locked()
    rc=child.wait()
    with guard:children.pop(lane,None)
    if stopping.is_set():raise InterruptedError('Stop requested')
    return rc

def on_signal(sig,frame):
    stopping.set()
    with guard:
        for child in children.values():
            if child.poll() is None:child.terminate()

signal.signal(signal.SIGTERM,on_signal);signal.signal(signal.SIGINT,on_signal)

def md5(path):
    h=hashlib.md5()
    with path.open('rb') as f:
        for b in iter(lambda:f.read(8*1024*1024),b''):
            if stopping.is_set():raise InterruptedError('Stop requested')
            h.update(b)
    return h.hexdigest()

def quarantine(part,reason):
    assert part.resolve().parent==D.resolve()
    q=D/'quarantine_proxy_dual_20261016_v1';q.mkdir(exist_ok=True)
    target=q/(part.name+'.'+str(time.time_ns()))
    assert target.resolve().parent==q.resolve() and not target.exists()
    part.rename(target)
    with guard:
        with (J/'quarantine.jsonl').open('a',encoding='utf-8') as f:
            f.write(json.dumps({'utc':now(),'original':str(part),'retained':str(target),'reason':reason})+'\n')

def one(n,lane):
    e=entries[n];dest=D/n;part=D/(n+'.partial')
    update(n,status='PREPARING',lane=lane,expected_bytes=e['bytes'],started_utc=now())
    assert not dest.exists(),'Unexpected destination preserved: '+str(dest)
    for attempt in range(1,6):
        if stopping.is_set():raise InterruptedError('Stop requested')
        update(n,attempt=attempt,status='PREPARING')
        if part.exists() and part.stat().st_size>e['bytes']:quarantine(part,'Oversized partial')
        with (J/'files'/(n+'.log')).open('a',encoding='utf-8') as log:
            if not part.exists() or part.stat().st_size<e['bytes']:
                assert os.statvfs(D).f_bavail*os.statvfs(D).f_frsize>1000000000000+e['bytes'],'Disk free floor'
                cmd=['curl','--proxy',PROXY,'--noproxy','','--fail','--location','--silent','--show-error','--connect-timeout','30',
                     '--speed-time','600','--speed-limit','1024','--continue-at','-','--output',str(part),e['url']]
                log.write('\nSTART '+now()+' '+json.dumps(cmd)+'\n');log.flush()
                update(n,status='DOWNLOADING',starting_bytes=part.stat().st_size if part.exists() else 0)
                rc=run_process(lane,cmd,log,n);log.write('CURL_EXIT '+str(rc)+' '+now()+'\n');log.flush()
                update(n,curl_exit_code=rc)
                if rc!=0:
                    update(n,status='RETRY_WAIT',error='curl exit '+str(rc))
                    if attempt<5:stopping.wait(30)
                    continue
            update(n,status='VALIDATING')
            if part.stat().st_size!=e['bytes']:
                update(n,error='Size mismatch');continue
            actual=md5(part)
            if actual!=e['md5']:
                quarantine(part,'MD5 mismatch: '+actual+' expected '+e['md5']);update(n,error='MD5 mismatch; preserved');continue
            gz=run_process(lane,['gzip','-t',str(part)],log,n)
            if gz!=0:quarantine(part,'gzip test failure');update(n,error='gzip test failure');continue
            assert not dest.exists()
            part.rename(dest)
            update(n,status='VERIFIED',bytes=e['bytes'],md5=actual,gzip_exit_code=0,error=None,completed_utc=now());return
    update(n,status='FAILED',finished_utc=now(),partial_bytes=part.stat().st_size if part.exists() else None)

def worker(lane):
    while not stopping.is_set():
        n=take(lane)
        if n is None:break
        try:one(n,lane)
        except InterruptedError as ex:update(n,status='INTERRUPTED',error=str(ex));break
        except Exception as ex:update(n,status='FAILED',error=repr(ex),traceback=traceback.format_exc())
    with guard:state['lanes'][lane]=None;save_locked()

update()
try:
    with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
        futures=[pool.submit(worker,lane) for lane in ['large','small']]
        for future in futures:future.result()
    rc=130 if stopping.is_set() else (0 if state['completed_files']==len(entries) else 1)
    update(status='INTERRUPTED' if rc==130 else ('COMPLETED_VERIFIED' if rc==0 else 'FINISHED_WITH_FAILURES'),exit_code=rc,finished_utc=now())
    if rc==0:write(J/'DOWNLOAD_COMPLETE.json',{'status':'COMPLETED_VERIFIED','files':len(entries),'bytes':state['total_expected_bytes'],'utc':now()})
except BaseException as ex:
    on_signal(None,None);rc=2;update(status='FAILED_RUNNER',exit_code=rc,error=repr(ex),traceback=traceback.format_exc())
(J/'exit_code.txt').write_text(str(rc),encoding='utf-8')
with (P/'HANDOFF.md').open('a',encoding='utf-8') as h:
    h.write('\n## '+now()+' / Dual proxy job ended\n'+state['status']+'; exit '+str(rc)+'; verified '+str(state['completed_files'])+'/120. Job: '+str(J)+'\n')
sys.exit(rc)
