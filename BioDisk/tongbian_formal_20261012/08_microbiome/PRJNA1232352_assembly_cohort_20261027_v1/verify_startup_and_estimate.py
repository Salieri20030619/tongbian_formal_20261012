import pathlib,json,datetime,os,select,time,hashlib
N=pathlib.Path(__file__).parent;P=json.loads((N/'cohort_plan.json').read_text());G=pathlib.Path(P['first_sample_gate'])
if not (G/'exit_code.txt').exists():
 try:fd=os.pidfd_open(int((G/'pid.txt').read_text()));select.select([fd],[],[]);os.close(fd)
 except ProcessLookupError:pass
assert (G/'exit_code.txt').read_text().strip()=='0','first sample failed'
start=time.monotonic();status=None;job=None
while time.monotonic()-start<1800:
 status=json.loads((N/'status.json').read_text())
 if status['status']=='FAILED_REVIEW_REQUIRED':raise RuntimeError(str(status))
 if status['status']=='RUNNING':
  job=N/'jobs'/status['run'];sfile=job/'status.json'
  if sfile.exists():
   s=json.loads(sfile.read_text());childfile=job/'megahit_child.json';inputfile=job/'INPUT_ACCEPTANCE.json'
   if s['phase']=='megahit' and s['status']=='RUNNING' and childfile.exists() and inputfile.exists():
    child=json.loads(childfile.read_text());proc=pathlib.Path('/proc')/str(child['pid'])/'cmdline'
    if proc.exists() and 'megahit' in proc.read_bytes().decode():break
 time.sleep(5)
else:raise RuntimeError('startup check timeout; preserve jobs and inspect status; do not stop compute')
assert job is not None
first=json.loads((G/'FULL_SAMPLE_COMPLETE.json').read_text());assert first['status']=='COMPLETED_VERIFIED';first_start=datetime.datetime.fromisoformat(json.loads((G/'resource_baseline.json').read_text())['utc']);first_end=datetime.datetime.fromisoformat(first['finished_utc']);hours=(first_end-first_start).total_seconds()/3600
remaining=sum(x['retained_pairs'] for x in P['pending_samples']);scaled=hours*remaining/first['full_input_pairs'];equal=59*hours
estimate={'created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'first_sample':first['run'],'first_sample_pairs':first['full_input_pairs'],'first_elapsed_hours_including_input_validation':hours,'remaining_samples':59,'remaining_pairs':remaining,'workers':1,'threads_per_sample':16,'remaining_estimate_hours_reads_scaled':scaled,'remaining_estimate_hours_equal_sample_time':equal,'interpretation':'two planning assumptions based on one successful sample, not statistical confidence bounds; assumes uninterrupted sequential operation; complexity/resource contention/guard waits may change runtime','no_parallel_speedup_assumed':True,'background_computation_continues':True}
(N/'TIME_ESTIMATE.json').write_text(json.dumps(estimate,indent=2)+'\n',encoding='utf-8')
receipt={'status':'COHORT_STARTUP_VERIFIED','utc':estimate['created_utc'],'first_sample_adopted':json.loads((N/'FIRST_SAMPLE_ADOPTED.json').read_text()),'completed':status['completed_samples'],'total':60,'running_sample':status['run'],'worker_pid':status['worker_pid'],'megahit_child_pid':child['pid'],'input_acceptance':json.loads(inputfile.read_text()),'workers':1,'threads':16,'time_estimate':str(N/'TIME_ESTIMATE.json'),'monitoring_stop_scope':'conversation only; background controller and intrinsic resource recording stay active'}
(N/'STARTUP_VERIFIED.json').write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf-8');print(json.dumps({'status':receipt['status'],'running_sample':status['run'],'completed':status['completed_samples'],'time_estimate':estimate}),flush=True)
