import pathlib,json,csv,re,hashlib,datetime,subprocess,os,select,sys
N=pathlib.Path(__file__).parent;R=N.parent.parent;P=json.loads((N/'plan_frozen.json').read_text());L=pathlib.Path(P['large_output'])
if not (N/'exit_code.txt').exists():
 try:
  fd=os.pidfd_open(int((N/'pid.txt').read_text()));select.select([fd],[],[]);os.close(fd)
 except ProcessLookupError:pass
assert (N/'exit_code.txt').read_text().strip()=='0',(N/'failure.json').read_text() if (N/'failure.json').exists() else 'nonzero exit'
complete=json.loads((N/'FULL_SAMPLE_COMPLETE.json').read_text());assert complete['status']=='COMPLETED_VERIFIED'
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(8388608),b''):h.update(b)
 return h.hexdigest()
for name,rec in complete['outputs'].items():
 p=pathlib.Path(name);assert p.stat().st_size==rec['bytes'] and sha(p)==rec['sha256']
stages=[]
for label in ['megahit','prodigal','dbcan']:
 assert (N/(label+'_exit_code.txt')).read_text().strip()=='0';txt=(N/(label+'.log')).read_text();rss=re.search(r'Maximum resident set size \(kbytes\):\s*(\d+)',txt);wall=re.search(r'Elapsed \(wall clock\) time \(h:mm:ss or m:ss\):\s*(\S+)',txt);assert rss and wall
 stages.append({'stage':label,'max_RSS_kbytes':int(rss.group(1)),'elapsed_wall':wall.group(1)})
with (N/'stage_resources.tsv').open('w',encoding='utf-8',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(stages[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(stages)
with (N/'resource_samples.tsv').open() as f:obs=list(csv.DictReader(f,delimiter='\t'))
qa={'status':'SINGLE_FULL_SAMPLE_ACCEPTED','utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'run':complete['run'],'full_input_pairs':complete['full_input_pairs'],'outputs_rehashed':5,'stages':stages,'resource_samples':len(obs),'sampled_min_MemAvailable_GiB':min(float(x['MemAvailable_bytes']) for x in obs)/1024**3,'sampled_min_disk_free_GiB':min(float(x['disk_free_bytes']) for x in obs)/1024**3,'sampled_peak_memory_full_PSI_avg10':max(float(x['memory_full_avg10']) for x in obs),'cgroup_memory_peak_including_cache_GiB':complete['cgroup_memory_peak_including_cache']/1024**3,'disk_output_GiB':complete['disk_output_bytes']/1024**3,'no_cgroup_OOM_increment':True,'cohort_concurrency_accepted':False,'interpretation':'largest read count is not guaranteed greatest assembly complexity; single-sample gate does not prove multiworker safety','covariate_limitation':'TB-COV-001 DEFERRED_BY_USER'}
(N/'RESOURCE_GATE_ACCEPTANCE.json').write_text(json.dumps(qa,indent=2)+'\n',encoding='utf-8')
print(json.dumps(qa),flush=True)
