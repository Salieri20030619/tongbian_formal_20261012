"""Concurrency-only migration: preserve scientific commands and pilot baseline."""
import argparse, concurrent.futures, fcntl, json, os, select, signal, threading, time, traceback
from pathlib import Path
import profile_pipeline as base
J=base.J
def procstat(pid):
 text=Path('/proc',str(pid),'stat').read_text()
 fields=text[text.rfind(')')+2:].split()
 return fields[0],int(fields[19]),int(fields[49])
def validate(row):
 out=base.OUT/'full'/row['run']; n=2*row['retained_pairs']
 receipt=json.loads((out/'input_receipt.json').read_text())
 assert receipt['reads']==n
 for tag in ['metaphlan','humann']:
  assert (out/(tag+'_exit_code.txt')).read_text().strip()=='0'
 mpa=base.validate_mpa(out/'metaphlan.tsv',n)
 tables={name:base.validate_humann_table(out/'humann'/(row['run']+'_'+name+'.tsv'),name=='pathcoverage') for name in ['genefamilies','pathabundance','pathcoverage']}
 result={'run':row['run'],'group':row['group'],'reads':n,'pilot':False,'metaphlan':mpa,'humann':tables,'completed_utc':base.now(),'profile_sha256':base.sha(out/'metaphlan.tsv')}
 old=out/'VALIDATED.json'
 if old.exists():
  prior=json.loads(old.read_text())
  assert result['profile_sha256']==prior['profile_sha256']
  for name in tables: assert tables[name]['sha256']==prior['humann'][name]['sha256']
  return prior
 base.write(old,result)
 return result

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--workers',type=int,required=True);ap.add_argument('--self-test',action='store_true')
 args=ap.parse_args();assert 2<=args.workers<=8
 if args.self_test:
  status,start,code=procstat(os.getpid());assert start>0 and code==0
  import subprocess
  child=subprocess.Popen(['sh','-c','exit 7'])
  fd=os.pidfd_open(child.pid);poll=select.poll();poll.register(fd,select.POLLIN);poll.poll()
  st,_,waitstatus=procstat(child.pid)
  assert st=='Z' and os.waitstatus_to_exitcode(waitstatus)==7
  assert child.wait()==7;os.close(fd)
  print('PASS: pidfd completion and /proc zombie exit status');return
 lock=(J/'worker_upgrade.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
 assert not (J/'worker_upgrade_started.json').exists(),'Migration already started'
 gate=json.loads((J/'PILOT_PASS.json').read_text())
 assert gate['status']=='PASS' and gate['pipeline_sha256']==base.sha(J/'profile_pipeline.py')
 assert gate['config_sha256']==base.sha(J/'config.json') and gate['samples_sha256']==base.sha(J/'samples.json')
 mem={x.split(':')[0]:int(x.split()[1]) for x in Path('/proc/meminfo').read_text().splitlines() if len(x.split())>=2}
 assert mem['MemAvailable']>100*1024**2
 oldpid=int((J/'full_pid.txt').read_text());oldstart=procstat(oldpid)[1]
 assert b'profile_pipeline.py\x00full' in Path('/proc',str(oldpid),'cmdline').read_bytes()
 backup=Path((J/'worker_upgrade_backup_path.txt').read_text())
 rows=base.allrows;byrun={r['run']:r for r in rows}
 # Prevalidate completed results before pausing the original scheduler.
 initial=json.loads((J/'full_status.json').read_text())
 for run,v in initial['samples'].items():
  if v['status']=='VERIFIED':validate(byrun[run])
 os.kill(oldpid,signal.SIGSTOP)
 committed=False
 try:
  for _ in range(100):
   if all(procstat(int(t.name))[0] in ('T','t') for t in Path('/proc',str(oldpid),'task').iterdir()):break
   time.sleep(.05)
  else:raise RuntimeError('Old scheduler did not stop')
  state=json.loads((J/'full_status.json').read_text())
  assert all(v['status'] in ('VERIFIED','RUNNING') for v in state['samples'].values())
  active={k:v.copy() for k,v in state['samples'].items() if v['status']=='RUNNING'}
  assert len(active)<=args.workers and active
  handles={}
  for run,v in active.items():
   assert v['phase']=='humann',('Unsupported adoption phase',run,v)
   pid=v['child_pid'];st,start,rc=procstat(pid)
   cmd=Path('/proc',str(pid),'cmdline').read_bytes()
   assert (st=='Z' or (b'humann' in cmd and run.encode() in cmd))
   # Parent remains stopped until exit statuses of all original children are captured.
   assert int(Path('/proc',str(pid),'status').read_text().split('PPid:')[1].split()[0])==oldpid
   handles[run]=(pid,start,os.pidfd_open(pid))
  pending=[r for r in rows if r['run'] not in state['samples']]
  for r in pending:assert not (base.OUT/'full'/r['run']).exists()
  base.write(backup/'full_status_at_pause.json',state)
  state.update(pid=os.getpid(),status='RUNNING',concurrent_samples=args.workers,threads_per_sample=8,previous_pid=oldpid,upgrade_utc=base.now())
  base.write(J/'worker_upgrade_started.json',{'utc':base.now(),'old_pid':oldpid,'old_start':oldstart,'new_pid':os.getpid(),'workers':args.workers,'adopted':active,'backup':str(backup),'original_pipeline_sha256':base.sha(J/'profile_pipeline.py'),'new_scheduler_sha256':base.sha(Path(__file__)),'scientific_parameters_changed':False})
  guard=threading.Lock();adopt_guard=threading.Lock();remaining=set(active);full_lock=[None]
  def update(sample_key=None,**kw):
   with guard:
    if sample_key:state['samples'].setdefault(sample_key,{}).update(kw)
    else:state.update(kw)
    state['completed_samples']=sum(v['status']=='VERIFIED' for v in state['samples'].values())
    state['updated_utc']=base.now();base.write(J/'full_status.json',state)
  def adopt(run):
   pid,start,fd=handles[run]
   try:
    update(run,adopted=True,phase='humann',child_pid=pid)
    poll=select.poll();poll.register(fd,select.POLLIN);poll.poll()
    st,actual_start,waitstatus=procstat(pid)
    assert actual_start==start and st=='Z'
    rc=os.waitstatus_to_exitcode(waitstatus)
    out=base.OUT/'full'/run
    (out/'humann_exit_code.txt').write_text(str(rc),encoding='utf-8')
    base.write(out/'adopted_exit_receipt.json',{'pid':pid,'start_ticks':start,'raw_wait_status':waitstatus,'exit_code':rc,'utc':base.now(),'source':'pidfd readiness then unreaped child /proc stat field 52; original parent stopped'})
    # Release original scheduler only after every child exit code was safely captured.
    with adopt_guard:
     remaining.remove(run)
     if not remaining:
      assert procstat(oldpid)[1]==oldstart
      os.kill(oldpid,signal.SIGKILL)
      f=(J/'full.lock').open('a');fcntl.flock(f,fcntl.LOCK_EX);full_lock[0]=f
      base.write(J/'old_scheduler_retired.json',{'utc':base.now(),'pid':oldpid,'reason':'All original child exits captured; replaced scheduling only'})
    assert rc==0,('Adopted HUMAnN failed',run,rc)
    update(run,phase='VALIDATING')
    result=validate(byrun[run]);update(run,status='VERIFIED',phase='COMPLETE',result=result)
    return result
   except BaseException as e:
    update(run,status='FAILED',error=repr(e),traceback=traceback.format_exc());raise
   finally:os.close(fd)
  update()
  (J/'full_pid.txt').write_text(str(os.getpid()),encoding='utf-8')
  base.write(J/'worker_runtime_config.json',{'concurrent_samples':args.workers,'threads_per_sample':8,'baseline_config_sha256':gate['config_sha256']})
  committed=True
  with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as pool:
   futures=[pool.submit(adopt,run) for run in active]
   futures += [pool.submit(base.one,row,'full',update) for row in pending]
   failures=[]
   for f in concurrent.futures.as_completed(futures):
    try:f.result()
    except BaseException as e:failures.append(repr(e))
  assert not failures,failures
  assert state['completed_samples']==len(rows)
  update(status='COMPLETED_VERIFIED',exit_code=0,finished_utc=base.now())
  base.write(J/'FULL_COMPLETE.json',{'status':'COMPLETED_VERIFIED','samples':len(rows),'results':[state['samples'][r['run']]['result'] for r in rows],'finished_utc':base.now()})
  (J/'full_exit_code.txt').write_text('0',encoding='utf-8')
 except BaseException:
  if not committed:
   if procstat(oldpid)[1]==oldstart:os.kill(oldpid,signal.SIGCONT)
  else:
   update(status='FAILED',error=traceback.format_exc(),exit_code=1)
   (J/'full_exit_code.txt').write_text('1',encoding='utf-8')
  raise
if __name__=='__main__':main()