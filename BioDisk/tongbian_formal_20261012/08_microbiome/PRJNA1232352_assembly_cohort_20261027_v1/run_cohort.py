import pathlib,json,subprocess,hashlib,datetime,os,select,fcntl,shutil,time,traceback,sys
N=pathlib.Path(__file__).parent;R=N.parent.parent;P=json.loads((N/'cohort_plan.json').read_text());G=pathlib.Path(P['first_sample_gate']);L=pathlib.Path(P['large_output_root']);lock=(N/'controller.lock').open('w');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB);complete=[]
def now():return datetime.datetime.now(datetime.timezone.utc).isoformat()
def save(p,x):
 t=p.with_suffix(p.suffix+'.tmp');t.write_text(json.dumps(x,indent=2)+'\n',encoding='utf-8');t.replace(p)
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(8388608),b''):h.update(b)
 return h.hexdigest()
def status(st,**kw):
 d={'status':st,'utc':now(),'controller_pid':os.getpid(),'completed_samples':len(complete),'total_samples':60,'workers':1,'threads_per_sample':16,**kw};save(N/'status.json',d);print(now(),st,kw,flush=True)
def wait_pid(pid,expected):
 p=pathlib.Path('/proc')/str(pid)/'cmdline'
 if not p.exists():return
 assert expected in p.read_bytes().decode().replace('\x00',' '),'pid ownership mismatch'
 try:fd=os.pidfd_open(pid);select.select([fd],[],[]);os.close(fd)
 except ProcessLookupError:pass
def verify(job):
 assert (job/'exit_code.txt').read_text().strip()=='0';receipt=json.loads((job/'FULL_SAMPLE_COMPLETE.json').read_text());assert receipt['status']=='COMPLETED_VERIFIED'
 for stage in ['megahit','prodigal','dbcan']:assert (job/(stage+'_exit_code.txt')).read_text().strip()=='0'
 for name,item in receipt['outputs'].items():
  p=pathlib.Path(name);assert p.exists() and p.stat().st_size==item['bytes'] and sha(p)==item['sha256'],name
 return receipt
try:
 status('WAITING_FOR_FIRST_FULL_SAMPLE')
 if not (G/'exit_code.txt').exists():wait_pid(int((G/'pid.txt').read_text()),'run_full_sample.py')
 first=verify(G);assert first['run']=='SRR32583461';complete.append({'run':first['run'],'job':str(G),'receipt_sha256':sha(G/'FULL_SAMPLE_COMPLETE.json'),'adopted':True})
 save(N/'FIRST_SAMPLE_ADOPTED.json',{'status':'PASS','utc':now(),'run':first['run'],'receipt':str(G/'FULL_SAMPLE_COMPLETE.json'),'outputs_rehashed':5,'no_recomputation':True});save(N/'completed_samples.json',complete)
 baseline=json.loads((G/'plan_frozen.json').read_text());template=(N/'worker_template.py').read_text()
 for sample in P['pending_samples']:
  run=sample['run'];job=N/'jobs'/run;large=L/run
  if (job/'FULL_SAMPLE_COMPLETE.json').exists() and (job/'exit_code.txt').read_text().strip()=='0':receipt=verify(job)
  else:
   job.mkdir(exist_ok=True);large.mkdir(exist_ok=True)
   if (job/'worker_pid.txt').exists() and not (job/'exit_code.txt').exists():
    pid=int((job/'worker_pid.txt').read_text());status('ADOPTING_ACTIVE_CHILD',run=run,worker_pid=pid);wait_pid(pid,str(job/'run_full_sample.py'));receipt=verify(job)
   else:
    assert not (large/'megahit').exists(),'existing incomplete assembly requires explicit recovery; preserving '+str(large)
    mem=lambda:int(next(l.split()[1] for l in pathlib.Path('/proc/meminfo').read_text().splitlines() if l.startswith('MemAvailable:')))*1024
    while mem()<100*1024**3 or shutil.disk_usage(L).free<300*1024**3:
     status('WAITING_RESOURCE',next_run=run);time.sleep(P['resource_wait_interval_seconds'])
    plan=dict(baseline);plan.update(scope='full-sample cohort execution after full-sample gate',sample=sample,full_input_pairs=sample['retained_pairs'],large_output=str(large),cohort_directory=str(N),cgroup_scope='shared cohort session: cumulative cgroup peak; per-stage GNU time RSS is sample-specific')
    save(job/'plan_frozen.json',plan);(job/'run_full_sample.py').write_text(template,encoding='utf-8')
    with (job/'runner.log').open('w',encoding='utf-8') as log:
     child=subprocess.Popen([sys.executable,'-u',str(job/'run_full_sample.py')],stdout=log,stderr=subprocess.STDOUT,start_new_session=True);(job/'worker_pid.txt').write_text(str(child.pid),encoding='utf-8');status('RUNNING',run=run,worker_pid=child.pid,sample_status=str(job/'status.json'));rc=child.wait()
    assert rc==0,'worker exit '+str(rc)+'; inspect '+str(job);receipt=verify(job)
  assert receipt['run']==run and receipt['full_input_pairs']==sample['retained_pairs'];complete.append({'run':run,'job':str(job),'receipt_sha256':sha(job/'FULL_SAMPLE_COMPLETE.json'),'adopted':False});save(N/'completed_samples.json',complete);status('SAMPLE_VERIFIED',run=run)
 assert len(complete)==60 and len({x['run'] for x in complete})==60
 save(N/'COHORT_COMPLETE.json',{'status':'COMPLETED_VERIFIED','samples':60,'completed_utc':now(),'results':complete,'scope':'assembly, Prodigal and CAZyme family annotation; not validated substrate reactions or whole-project completion'});(N/'exit_code.txt').write_text('0\n',encoding='utf-8');status('COMPLETED_VERIFIED')
except BaseException as e:
 save(N/'failure.json',{'error':str(e),'utc':now(),'completed':len(complete)});(N/'exit_code.txt').write_text('1\n',encoding='utf-8');status('FAILED_REVIEW_REQUIRED',error=str(e));traceback.print_exc();raise
