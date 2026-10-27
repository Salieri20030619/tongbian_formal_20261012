import pathlib,json,subprocess,hashlib,datetime,os,select,fcntl,shutil,time,traceback,sys,csv,signal
N=pathlib.Path(__file__).parent;P=json.loads((N/'cohort_plan.json').read_text());G=pathlib.Path(P['first_sample_gate']);L=pathlib.Path(P['large_output_root']);M=json.loads((N/'worker8_migration.json').read_text());lock=(N/'controller8.lock').open('w');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB);complete=json.loads((N/'completed_samples.json').read_text());active={};failed={};old_retired=False;guard=False;bad_streak=0;good_streak=0;last_cpu=None;max_simultaneous=0

def now():return datetime.datetime.now(datetime.timezone.utc).isoformat()
def save(p,x):
 t=p.with_suffix(p.suffix+'.tmp');t.write_text(json.dumps(x,indent=2)+'\n',encoding='utf-8');t.replace(p)
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(8388608),b''):h.update(b)
 return h.hexdigest()
def proc(pid):
 p=pathlib.Path('/proc')/str(pid);f=(p/'stat').read_text().rsplit(')',1)[1].split();return {'state':f[0],'ppid':int(f[1]),'ticks':int(f[11])+int(f[12]),'starttime':int(f[19]),'exit_raw':int(f[49]),'cmd':(p/'cmdline').read_bytes().replace(b'\0',b' ').decode()}
def ownership(pid,needle):
 x=proc(pid);assert needle in x['cmd'];return x

def verify(job):
 assert (job/'exit_code.txt').read_text().strip()=='0';r=json.loads((job/'FULL_SAMPLE_COMPLETE.json').read_text());assert r['status']=='COMPLETED_VERIFIED'
 for stage in ['megahit','prodigal','dbcan']:assert (job/(stage+'_exit_code.txt')).read_text().strip()=='0'
 for name,x in r['outputs'].items():
  f=pathlib.Path(name);assert f.stat().st_size==x['bytes'] and sha(f)==x['sha256']
 return r

def metrics():
 global last_cpu
 mem=int(next(l.split()[1] for l in pathlib.Path('/proc/meminfo').read_text().splitlines() if l.startswith('MemAvailable:')))*1024
 psi={k:{l.split()[0]:float(l.split()[1].split('=')[1]) for l in pathlib.Path('/proc/pressure',k).read_text().splitlines()} for k in ['cpu','memory','io']}
 cur=[int(v) for v in pathlib.Path('/proc/stat').read_text().splitlines()[0].split()[1:9]];busy=None;iowait=None
 if last_cpu:
  d=[v-w for v,w in zip(cur,last_cpu)];total=sum(d);busy=100*(total-d[3]-d[4])/total if total else 0;iowait=100*d[4]/total if total else 0
 last_cpu=cur
 phases={};groups={};running={}
 for run,a in active.items():
  job=pathlib.Path(a['job']);s=json.loads((job/'status.json').read_text()) if (job/'status.json').exists() else {'phase':'starting','status':'RUNNING'};phase=s['phase'];phases[phase]=phases.get(phase,0)+1;running[run]={'pid':a['pid'],'phase':phase,'worker_status':s['status'],'adopted':a['adopted']}
  try:
   cg=pathlib.Path('/sys/fs/cgroup')/pathlib.Path('/proc',str(a['pid']),'cgroup').read_text().strip().split('::')[1].lstrip('/');groups[str(cg)]={'current':int((cg/'memory.current').read_text()),'events':(cg/'memory.events').read_text()}
  except FileNotFoundError:pass
 return {'utc':now(),'MemAvailable_GiB':mem/1024**3,'disk_free_GiB':shutil.disk_usage(L).free/1024**3,'cpu_busy_percent':busy,'iowait_percent':iowait,'cpu_some_avg10':psi['cpu']['some'],'memory_full_avg10':psi['memory']['full'],'io_full_avg10':psi['io']['full'],'active_count':len(active),'phases':phases,'active':running,'distinct_cgroups':groups,'cgroup_current_sum_GiB':sum(x['current'] for x in groups.values())/1024**3}

def status(st,m=None):
 save(N/'status.json',{'status':st,'utc':now(),'controller_pid':os.getpid(),'completed_samples':len(complete),'total_samples':60,'workers':8,'threads_per_sample':16,'active_count':len(active),'active':m['active'] if m else {r:{k:v for k,v in a.items() if k in ['pid','job','adopted']} for r,a in active.items()},'failed':failed,'admission_paused':guard or bool(failed),'old_scheduler_retired':old_retired,'max_simultaneous_workers':max_simultaneous})

def launch(sample):
 run=sample['run'];job=N/'jobs'/run;large=L/run;assert not job.exists() and not large.exists(),('pre-existing unadopted job',run);job.mkdir();large.mkdir();plan=json.loads((G/'plan_frozen.json').read_text());plan.update(scope='eight-worker cohort after full sample gate',sample=sample,full_input_pairs=sample['retained_pairs'],large_output=str(large),cohort_directory=str(N),cgroup_scope='shared worker8 session: cgroup peak cumulative across simultaneous workers; per-stage GNU time RSS sample specific')
 save(job/'plan_frozen.json',plan);(job/'run_full_sample.py').write_text((N/'worker_template.py').read_text(),encoding='utf-8')
 with (job/'runner.log').open('w',encoding='utf-8') as f:p=subprocess.Popen([sys.executable,'-u',str(job/'run_full_sample.py')],stdout=f,stderr=subprocess.STDOUT,start_new_session=True)
 (job/'worker_pid.txt').write_text(str(p.pid),encoding='utf-8');fd=os.pidfd_open(p.pid);active[run]={'pid':p.pid,'job':str(job),'fd':fd,'popen':p,'adopted':False};print(now(),'LAUNCHED',run,p.pid,flush=True)

def retire_old():
 global old_retired
 pid=M['old_controller_pid'];x=ownership(pid,str(N/'run_cohort.py'));assert x['starttime']==M['old_controller_starttime'] and x['state'] in ['T','t'];os.kill(pid,signal.SIGKILL);save(N/'worker8_old_scheduler_retired.json',{'utc':now(),'old_pid':pid,'reason':'adopted worker actual exit status recorded before retiring stopped scheduler','signal':'SIGKILL to scheduler PID only','compute_signal_sent':False});old_retired=True

try:
 for item in complete:
  job=pathlib.Path(item['job']);assert sha(job/'FULL_SAMPLE_COMPLETE.json')==item['receipt_sha256'];verify(job)
 old=M['active'];x=ownership(old['pid'],str(pathlib.Path(old['job'])/'run_full_sample.py'));assert x['starttime']==old['starttime'];fd=os.pidfd_open(old['pid']);active[old['run']]={'pid':old['pid'],'job':old['job'],'fd':fd,'popen':None,'adopted':True};print(now(),'ADOPTED',old['run'],old['pid'],flush=True)
 save(N/'worker8_started.json',{'utc':now(),'pid':os.getpid(),'workers':8,'adopted_run':old['run'],'adopted_pid':old['pid'],'old_scheduler_stopped':True,'scientific_worker_sha256':sha(N/'worker_template.py')})
 pending=[s for s in P['pending_samples'] if s['run'] not in {x['run'] for x in complete} and s['run'] not in active]
 cfg=json.loads((N/'worker8_runtime.json').read_text());baseline=(N/'worker_template.py').read_bytes();assert hashlib.sha256(baseline).hexdigest()==M['worker_template_sha256']
 with (N/'worker8_resources.jsonl').open('a',encoding='utf-8',buffering=1) as logfile:
  while active or pending:
   cfg=json.loads((N/'worker8_runtime.json').read_text());m=metrics();logfile.write(json.dumps(m)+'\n');bad=m['MemAvailable_GiB']<cfg['min_available_GiB'] or m['disk_free_GiB']<cfg['min_disk_free_GiB'] or m['memory_full_avg10']>cfg['memory_full_psi_pause'] or m['io_full_avg10']>cfg['io_full_psi_pause'];bad_streak=bad_streak+1 if bad else 0;good_streak=good_streak+1 if not bad and m['MemAvailable_GiB']>cfg['resume_available_GiB'] else 0
   if bad_streak>=3:guard=True
   if guard and good_streak>=3:guard=False
   # Do not refill under low immediate headroom even before the sustained-pressure guard fires.
   while pending and len(active)<int(cfg['target_workers']) and not guard and not failed and m['MemAvailable_GiB']>=cfg['min_available_GiB'] and m['disk_free_GiB']>=cfg['min_disk_free_GiB']:
    launch(pending.pop(0))
   max_simultaneous=max(max_simultaneous,len(active));m=metrics();status('RUNNING' if active else 'WAITING_RESOURCE',m)
   if not active:
    if failed:break
    time.sleep(30);continue
   ready=select.select([a['fd'] for a in active.values()],[],[],10)[0]
   for run,a in list(active.items()):
    if a['fd'] not in ready:continue
    job=pathlib.Path(a['job'])
    if a['adopted']:
     x=proc(a['pid']);assert x['state']=='Z';rc=os.waitstatus_to_exitcode(x['exit_raw']);save(N/'worker8_adopted_exit.json',{'utc':now(),'run':run,'pid':a['pid'],'raw_wait_status':x['exit_raw'],'actual_exit_code':rc,'method':'pidfd readiness plus unreaped /proc zombie exit_code; original parent stopped'});retire_old()
    else:rc=a['popen'].wait()
    os.close(a['fd']);del active[run]
    try:
     assert rc==0,('worker exit',rc);receipt=verify(job);complete.append({'run':run,'job':str(job),'receipt_sha256':sha(job/'FULL_SAMPLE_COMPLETE.json'),'adopted':a['adopted']});save(N/'completed_samples.json',complete);print(now(),'VERIFIED',run,flush=True)
    except Exception as e:failed[run]=str(e);save(N/'worker8_failures.json',failed);print(now(),'FAILED',run,str(e),flush=True)
  if failed:raise RuntimeError('sample failures; preserved outputs; no refill')
 assert len(complete)==60 and len({x['run'] for x in complete})==60;save(N/'COHORT_COMPLETE.json',{'status':'COMPLETED_VERIFIED','samples':60,'completed_utc':now(),'results':complete,'workers':8,'scope':'assembly, Prodigal and CAZyme family annotation; no specific substrate or whole-project validation'});(N/'exit_code.txt').write_text('0\n',encoding='utf-8');status('COMPLETED_VERIFIED')
except BaseException as e:
 save(N/'worker8_failure.json',{'utc':now(),'error':str(e),'active_workers_preserved':[{k:v for k,v in a.items() if k in ['pid','job','adopted']} for a in active.values()]});status('FAILED_REVIEW_REQUIRED');traceback.print_exc();raise
