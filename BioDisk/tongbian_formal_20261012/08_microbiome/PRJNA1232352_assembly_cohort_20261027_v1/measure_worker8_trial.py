import pathlib,json,time,datetime,hashlib,os,subprocess
N=pathlib.Path(__file__).parent;M=json.loads((N/'worker8_migration.json').read_text());start=time.monotonic();window=None;observations=[];cpu_acc={};prev={};first8=None;seen_oom={};window_start_cpu={}
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def descendants(pid):
 out={};queue=[pid];seen=set()
 while queue:
  p=queue.pop()
  if p in seen:continue
  seen.add(p)
  try:
   stat=pathlib.Path('/proc',str(p),'stat').read_text().rsplit(')',1)[1].split();kids=pathlib.Path('/proc',str(p),'task',str(p),'children').read_text().split();queue.extend(map(int,kids))
   if p!=pid:out[str(p)+':'+stat[19]]=int(stat[11])+int(stat[12])
  except FileNotFoundError:pass
 return out
with (N/'worker8_trial_observations.jsonl').open('w',encoding='utf-8',buffering=1) as output:
 while time.monotonic()-start<3600:
  state=json.loads((N/'status.json').read_text());assert state['status']!='FAILED_REVIEW_REQUIRED',state;assert not state.get('failed'),state
  raw=(N/'worker8_resources.jsonl').read_text();lines=raw.splitlines();lines=lines if raw.endswith('\n') else lines[:-1];m=json.loads(lines[-1]);active=state.get('active',{});row={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'active_count':len(active),'phases':{r:a.get('phase') for r,a in active.items()},'metric':m}
  for r,a in active.items():
   cur=descendants(a['pid']);old=prev.get(r,{});cpu_acc[r]=cpu_acc.get(r,0)+sum(max(0,v-old[k]) for k,v in cur.items() if k in old);prev[r]=cur
  for cg,x in m['distinct_cgroups'].items():
   ev={l.split()[0]:int(l.split()[1]) for l in x['events'].splitlines()};initial=seen_oom.setdefault(cg,ev)
   assert all(ev.get(k,0)==initial.get(k,0) for k in ['oom','oom_kill']),(cg,ev)
  row['descendant_CPU_seconds_observed']={r:v/os.sysconf('SC_CLK_TCK') for r,v in cpu_acc.items()};output.write(json.dumps(row)+'\n');observations.append(row)
  compute=len(active)==8 and all(a.get('phase') in ['megahit','prodigal','dbcan'] and a.get('worker_status')=='RUNNING' for a in active.values())
  if compute:
   if window is None:window=time.monotonic();first8=list(active);window_start_cpu=dict(cpu_acc);print('EIGHT_COMPUTE_STARTED',row['utc'],row['phases'],flush=True)
   if time.monotonic()-window>=300 and all(cpu_acc.get(r,0)-window_start_cpu.get(r,0)>os.sysconf('SC_CLK_TCK') for r in first8):break
  else:window=None
  time.sleep(10)
 else:raise RuntimeError('No sustained eight-compute window within 60 minutes; inspect evidence, do not infer pass')
assert sha(N/'worker_template.py')==M['worker_template_sha256'];assert sha(pathlib.Path(M['active']['job'])/'run_full_sample.py')==M['active_worker_script_sha256'];assert sha(pathlib.Path(M['active']['job'])/'plan_frozen.json')==M['active_plan_sha256']
for run in first8:
 job=N/'jobs'/run;p=json.loads((job/'plan_frozen.json').read_text());assert p['threads']==16 and p['megahit_memory_bytes']==68719476736 and p['min_contig_length']==500 and p['megahit_k_list']=='21,29,39,59,79,99,119,141'
 assert sha(job/'run_full_sample.py')==M['worker_template_sha256'];assert json.loads((job/'INPUT_ACCEPTANCE.json').read_text())['status']=='PASS'
metrics=[r['metric'] for r in observations];window_metrics=[r['metric'] for r in observations[-31:]]
summary={'status':'EIGHT_WORKER_LOAD_TRIAL_PASS','utc':observations[-1]['utc'],'workers':8,'threads_per_sample_max':16,'original_worker_preserved':M['active'],'all_eight_inputs_verified':True,'scientific_worker_and_active_plan_hashes_unchanged':True,'sustained_eight_compute_seconds':time.monotonic()-window,'window_samples':len(window_metrics),'phases_at_end':observations[-1]['phases'],'min_MemAvailable_GiB_all_observations':min(x['MemAvailable_GiB'] for x in metrics),'peak_distinct_cgroup_current_sum_GiB':max(x['cgroup_current_sum_GiB'] for x in metrics),'cpu_busy_percent_window_mean':sum(x['cpu_busy_percent'] for x in window_metrics if x['cpu_busy_percent'] is not None)/len([x for x in window_metrics if x['cpu_busy_percent'] is not None]),'iowait_percent_window_max':max(x['iowait_percent'] or 0 for x in window_metrics),'cpu_some_PSI_avg10_window_max':max(x['cpu_some_avg10'] for x in window_metrics),'memory_full_PSI_avg10_window_max':max(x['memory_full_avg10'] for x in window_metrics),'io_full_PSI_avg10_window_max':max(x['io_full_avg10'] for x in window_metrics),'no_observed_cgroup_OOM_increment':True,'compute_window_tool_CPU_seconds':{r:(cpu_acc[r]-window_start_cpu.get(r,0))/os.sysconf('SC_CLK_TCK') for r in first8},'scope':'sustained simultaneous compute load and input/provenance acceptance; not full-cohort completion or eightfold throughput proof','old_scheduler_retirement':'read worker8_old_scheduler_retired.json; may remain stopped until adopted child exits','runtime_guards':json.loads((N/'worker8_runtime.json').read_text())}
(N/'WORKER8_LOAD_ACCEPTANCE.json').write_text(json.dumps(summary,indent=2)+'\n',encoding='utf-8');print(json.dumps(summary),flush=True)
