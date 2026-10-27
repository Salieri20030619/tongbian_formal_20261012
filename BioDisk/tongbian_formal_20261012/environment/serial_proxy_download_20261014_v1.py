from pathlib import Path
import os, sys, json, datetime, hashlib, subprocess, fcntl, time, signal, traceback

P = Path('/mnt/BioDisk/wangmingyuan/tongbian_formal_20261012')
D = Path('/mnt/BioData/bio/wangmingyuan/tongbian_formal_20261012/raw_reads/PRJNA1232352')
OLD = P/'logs/jobs/PRJNA1232352_direct_download_v1'
J = P/'logs/jobs/PRJNA1232352_proxy_serial_20261014_v1'
MAN = P/'01_inventory/resume_20261013/PRJNA1232352_download_manifest.json'
PROXY = 'http://127.0.0.1:7897'
now = lambda: datetime.datetime.now(datetime.timezone.utc).isoformat()

def write_json(path, obj):
    tmp = path.with_suffix(path.suffix+'.tmp')
    tmp.write_text(json.dumps(obj, indent=2, ensure_ascii=False), encoding='utf-8')
    os.replace(tmp, path)

def md5(path):
    h = hashlib.md5()
    with path.open('rb') as f:
        for b in iter(lambda: f.read(8*1024*1024), b''):
            h.update(b)
    return h.hexdigest()

if '--launch' in sys.argv:
    assert not J.exists(), 'Existing job: inspect instead of launching twice'
    old = json.loads((OLD/'status.json').read_text(encoding='utf-8'))
    assert not Path('/proc', str(old['runner_pid'])).exists(), 'Old runner still present'
    lock = (OLD/'run.lock').open('a')
    fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    entries = json.loads(MAN.read_text(encoding='utf-8'))
    J.mkdir(); (J/'files').mkdir()
    failed = [n for n, v in old['files'].items() if v['status']=='FAILED']
    active = [n for n, v in old['files'].items() if v['status']=='DOWNLOADING']
    done = {n:v for n,v in old['files'].items() if v['status'] in ('VERIFIED','REUSED_VERIFIED')}
    for n,v in done.items():
        assert (D/n).is_file() and (D/n).stat().st_size == v['bytes'], n
    names = failed + active + [e['filename'] for e in entries if e['filename'] not in failed+active and e['filename'] not in done]
    cfg = {'created_utc':now(), 'workers':1, 'proxy':PROXY, 'direct_fallback':False,
           'manifest':str(MAN), 'manifest_sha256':hashlib.sha256(MAN.read_bytes()).hexdigest(),
           'priority_failed':failed, 'queue':names, 'previous_verified':done,
           'user_instruction':'Prioritize failed downloads, use proxy, one file at a time',
           'integrity':'expected size + ENA MD5 + gzip -t before atomic rename',
           'retry':'fresh curl process resumes current file; no curl internal retry; 5 attempts',
           'corrupt_policy':'retain full-size MD5-failed partials under quarantine; redownload',
           'tunnel_dependency':'Windows SSH sgs10-codex forwards server 7897 to Windows 7897'}
    write_json(J/'config.json',cfg)
    fcntl.flock(lock,fcntl.LOCK_UN); lock.close()
    with (J/'launcher.log').open('ab') as log:
        proc = subprocess.Popen(['python3',str(Path(__file__).resolve()),'--run'],stdin=subprocess.DEVNULL,stdout=log,stderr=log,start_new_session=True,cwd=str(P))
    (J/'runner_pid.txt').write_text(str(proc.pid),encoding='utf-8')
    print(json.dumps({'pid':proc.pid,'priority':failed,'queued':len(names),'job':str(J)}))
    sys.exit(0)

assert '--run' in sys.argv
lock = (OLD/'run.lock').open('a')
fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
cfg = json.loads((J/'config.json').read_text(encoding='utf-8'))
assert hashlib.sha256(MAN.read_bytes()).hexdigest()==cfg['manifest_sha256']
entries = {e['filename']:e for e in json.loads(MAN.read_text(encoding='utf-8'))}
state = {'status':'RUNNING','started_utc':now(),'runner_pid':os.getpid(),'workers':1,'proxy':PROXY,
         'total_files':len(entries),'total_expected_bytes':sum(e['bytes'] for e in entries.values()),
         'files':{n:dict(v,status='REUSED_VERIFIED',verification_source=str(OLD/'status.json')) for n,v in cfg['previous_verified'].items()}}
child = None

def save():
    state['updated_utc']=now()
    state['completed_files']=sum(v['status'] in ('VERIFIED','REUSED_VERIFIED') for v in state['files'].values())
    state['failed_files']=sum(v['status']=='FAILED' for v in state['files'].values())
    write_json(J/'status.json',state)

def interrupted(sig, frame):
    global child
    if child and child.poll() is None:
        child.terminate(); child.wait(timeout=30)
    raise InterruptedError('Received signal '+str(sig))

signal.signal(signal.SIGTERM,interrupted)
signal.signal(signal.SIGINT,interrupted)

def quarantine(part, reason):
    assert part.resolve().parent == D.resolve()
    q=D/'quarantine_proxy_serial_20261014_v1'; q.mkdir(exist_ok=True)
    target=q/(part.name+'.'+str(time.time_ns()))
    assert target.resolve().parent==q.resolve() and not target.exists()
    part.rename(target)
    with (J/'quarantine.jsonl').open('a',encoding='utf-8') as log:
        log.write(json.dumps({'time_utc':now(),'original':str(part),'retained':str(target),'reason':reason})+'\n')

def one(n):
    global child
    e=entries[n]; dest=D/n; part=D/(n+'.partial')
    rec={'status':'PREPARING','expected_bytes':e['bytes'],'started_utc':now()}
    state['files'][n]=rec; state['current_file']=n; save()
    if dest.exists():
        raise RuntimeError('Unexpected existing destination; preserved: '+str(dest))
    for attempt in range(1,6):
        rec.update(attempt=attempt,status='PREPARING');save()
        if part.exists() and part.stat().st_size>e['bytes']:
            quarantine(part,'Oversized partial retained')
        if not part.exists() or part.stat().st_size<e['bytes']:
            free=os.statvfs(D).f_bavail*os.statvfs(D).f_frsize
            if free<1000000000000+e['bytes']:raise RuntimeError('Disk free floor reached')
            cmd=['curl','--proxy',PROXY,'--noproxy','','--fail','--location','--show-error','--silent',
                 '--connect-timeout','30','--speed-time','600','--speed-limit','1024',
                 '--continue-at','-','--output',str(part),e['url']]
            with (J/'files'/(n+'.log')).open('a',encoding='utf-8') as log:
                log.write('\nSTART '+now()+' '+json.dumps(cmd)+'\n');log.flush()
                child=subprocess.Popen(cmd,stdin=subprocess.DEVNULL,stdout=log,stderr=log)
                rec.update(status='DOWNLOADING',child_pid=child.pid,starting_bytes=part.stat().st_size if part.exists() else 0);save()
                rc=child.wait();child=None
                rec['curl_exit_code']=rc; log.write('CURL_EXIT '+str(rc)+' '+now()+'\n')
            if rc!=0:
                rec.update(status='RETRY_WAIT',error='curl exit '+str(rc));save()
                if attempt<5:time.sleep(30)
                continue
        rec.update(status='VALIDATING');save()
        if part.stat().st_size!=e['bytes']:
            rec['error']='Size mismatch';continue
        actual=md5(part)
        if actual!=e['md5']:
            quarantine(part,'MD5 mismatch: '+actual+' expected '+e['md5']);rec['error']='MD5 mismatch; bad partial retained';continue
        with (J/'files'/(n+'.log')).open('a',encoding='utf-8') as log:
            gz=subprocess.run(['gzip','-t',str(part)],stdout=log,stderr=log)
        if gz.returncode!=0:
            quarantine(part,'gzip -t failed');rec['error']='gzip failed';continue
        assert not dest.exists()
        part.rename(dest)
        rec.update(status='VERIFIED',bytes=e['bytes'],md5=actual,gzip_exit_code=0,completed_utc=now())
        rec.pop('error',None);save();return
    rec.update(status='FAILED',finished_utc=now(),partial_bytes=part.stat().st_size if part.exists() else None);save()

save()
try:
    for n in cfg['queue']:one(n)
    rc=0 if state['completed_files']==len(entries) else 1
    state.update(status='COMPLETED_VERIFIED' if rc==0 else 'FINISHED_WITH_FAILURES',exit_code=rc,completed_utc=now(),current_file=None);save()
    if rc==0:write_json(J/'DOWNLOAD_COMPLETE.json',{'status':'COMPLETED_VERIFIED','files':len(entries),'bytes':state['total_expected_bytes'],'completed_utc':now()})
except BaseException as ex:
    rc=130 if isinstance(ex,InterruptedError) else 2
    state.update(status='INTERRUPTED' if rc==130 else 'FAILED_RUNNER',exit_code=rc,error=repr(ex),traceback=traceback.format_exc());save()
(J/'exit_code.txt').write_text(str(rc),encoding='utf-8')
with (P/'HANDOFF.md').open('a',encoding='utf-8') as h:
    h.write('\n## '+now()+' / Serial proxy download ended\n'+state['status']+'; exit '+str(rc)+'; verified '+str(state['completed_files'])+'/120. Job: '+str(J)+'\n')
sys.exit(rc)
