from pathlib import Path
import datetime, json, time, subprocess, os, signal, traceback

P=Path('/mnt/BioDisk/wangmingyuan/tongbian_formal_20261012')
D=Path('/mnt/BioData/bio/wangmingyuan/tongbian_formal_20261012/raw_reads/PRJNA1232352')
J=P/'logs/jobs/PRJNA1232352_proxy_serial_20261014_v1'
stamp=datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%dT%H%M%SZ')
O=P/'logs/jobs'/('proxy_speed_comparison_'+stamp)
O.mkdir()
now=lambda:datetime.datetime.now(datetime.timezone.utc).isoformat()
s=json.loads((J/'status.json').read_text(encoding='utf-8'))
primary=s['current_file'];primary_pid=s['files'][primary]['child_pid']
small='SRR32583445_1.fastq.gz'
manifest=json.loads((P/'01_inventory/resume_20261013/PRJNA1232352_download_manifest.json').read_text(encoding='utf-8'))
e=next(x for x in manifest if x['filename']==small)
assert primary!=small and not (D/small).exists() and not (D/(small+'.partial')).exists()
assert '--proxy' in (Path('/proc')/str(primary_pid)/'cmdline').read_bytes().decode()
trial_dir=D/('speed_trial_'+stamp);trial_dir.mkdir()
secondary=trial_dir/(small+'.partial')
child=None;rows=[];phases=[];rc=0
def save(name,obj):
    (O/name).write_text(json.dumps(obj,indent=2),encoding='utf-8')
def sample():
    s=json.loads((J/'status.json').read_text(encoding='utf-8'))
    assert s['current_file']==primary,'Primary changed: stop experiment rather than mix transfers'
    f=D/(primary+'.partial')
    assert f.exists(),'Primary finished or moved'
    return {'time_utc':now(),'clock':time.monotonic(),'primary_bytes':f.stat().st_size,'secondary_bytes':secondary.stat().st_size if secondary.exists() else 0}
def phase(name,seconds):
    a=sample();prev=a
    for _ in range(seconds//30):
        time.sleep(30)
        b=sample();elapsed=b['clock']-prev['clock']
        r={'phase':name,'time_utc':b['time_utc'],'seconds':elapsed,
           'primary_MB_s':(b['primary_bytes']-prev['primary_bytes'])/elapsed/1e6,
           'secondary_MB_s':(b['secondary_bytes']-prev['secondary_bytes'])/elapsed/1e6}
        r['total_MB_s']=r['primary_MB_s']+r['secondary_MB_s'];rows.append(r)
        with (O/'intervals.jsonl').open('a',encoding='utf-8') as h:h.write(json.dumps(r)+'\n')
        print(json.dumps(r),flush=True);prev=b
        if name=='dual' and child.poll() is not None:
            raise RuntimeError('Secondary ended early, curl exit '+str(child.returncode))
    elapsed=b['clock']-a['clock']
    result={'phase':name,'seconds':elapsed,'primary_MB_s':(b['primary_bytes']-a['primary_bytes'])/elapsed/1e6,
            'secondary_MB_s':(b['secondary_bytes']-a['secondary_bytes'])/elapsed/1e6}
    result['total_MB_s']=result['primary_MB_s']+result['secondary_MB_s'];phases.append(result)
    save('phase_summary.json',phases)
def terminate_secondary():
    global child
    if child and child.poll() is None:
        child.terminate()
        try:child.wait(timeout=15)
        except subprocess.TimeoutExpired:child.kill();child.wait()
def interrupted(sig,frame):raise InterruptedError('Signal '+str(sig))
signal.signal(signal.SIGTERM,interrupted);signal.signal(signal.SIGINT,interrupted)
save('config.json',{'created_utc':now(),'primary':primary,'primary_pid':primary_pid,'secondary':e,
     'design':'single 60s; dual 180s; single 60s; sample every 30s','proxy':'http://127.0.0.1:7897',
     'secondary_partial':str(secondary),'existing_runner_unchanged':True})
print(json.dumps({'experiment_dir':str(O),'primary':primary,'secondary':small}),flush=True)
try:
    phase('single_before',60)
    cmd=['curl','--proxy','http://127.0.0.1:7897','--noproxy','','--fail','--location','--silent','--show-error',
         '--connect-timeout','30','--speed-time','600','--speed-limit','1024','--continue-at','-',
         '--output',str(secondary),e['url']]
    with (O/'secondary_curl.log').open('ab') as log:
        child=subprocess.Popen(cmd,stdin=subprocess.DEVNULL,stdout=log,stderr=log)
    save('secondary_process.json',{'pid':child.pid,'command':cmd,'started_utc':now()})
    phase('dual',180)
    terminate_secondary()
    phase('single_after',60)
except BaseException as ex:
    rc=1;save('error.json',{'error':repr(ex),'traceback':traceback.format_exc()})
finally:
    terminate_secondary()
    # Only move the retained trial partial into the normal queue when no process writes that exact path.
    canonical=D/(small+'.partial');conflicts=[]
    for q in Path('/proc').iterdir():
        if q.name.isdigit():
            try:
                args=(q/'cmdline').read_bytes().split(b'\0')
                if str(canonical).encode() in args:conflicts.append(int(q.name))
            except (FileNotFoundError,ProcessLookupError,PermissionError):pass
    retained=str(secondary)
    if secondary.exists() and not canonical.exists() and not (D/small).exists() and not conflicts:
        assert secondary.resolve().parent==trial_dir.resolve() and canonical.resolve().parent==D.resolve()
        secondary.rename(canonical);retained=str(canonical)
    out={'status':'COMPLETED' if rc==0 else 'FAILED','exit_code':rc,'finished_utc':now(),'phases':phases,
         'secondary_curl_exit':child.returncode if child else None,'secondary_stopped_intentionally':True,
         'secondary_partial_retained':retained,'retained_bytes':Path(retained).stat().st_size if Path(retained).exists() else 0,
         'original_serial_runner_unchanged':True,'integrity_not_yet_verified':True}
    if len(phases)==3:
        baseline=(phases[0]['total_MB_s']*phases[0]['seconds']+phases[2]['total_MB_s']*phases[2]['seconds'])/(phases[0]['seconds']+phases[2]['seconds'])
        out['single_baseline_MB_s']=baseline
        out['dual_vs_single_percent']=(phases[1]['total_MB_s']/baseline-1)*100 if baseline else None
    save('result.json',out);(O/'exit_code.txt').write_text(str(rc),encoding='utf-8')
    with (P/'HANDOFF.md').open('a',encoding='utf-8') as h:
        h.write('\n## '+now()+' / Proxy concurrency speed trial\n'+json.dumps(out)+'\nEvidence: '+str(O)+'\n')
    print('RESULT '+json.dumps(out),flush=True)
raise SystemExit(rc)
