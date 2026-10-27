from pathlib import Path
import sys,os,json,datetime,socket,hashlib,zipfile,subprocess
root=Path('/mnt/BioDisk/wangmingyuan/tongbian_formal_20261012');assert socket.gethostname()=='osv01'
staging=root/'06_human_selected_metabolites/limma_staging_v1'
p=staging/'tongbian_human_PMC9366892_deidentified_inputs_20261012.zip'
assert p.stat().st_size==40716 and hashlib.sha256(p.read_bytes()).hexdigest()=='634eba9e4c718d6051b9bce01695f306a574d4c2184cf3dbdd71dcc42b324bb0'
inp=root/'03_inputs/PMC9366892_deidentified_v1';inp.mkdir(parents=True,exist_ok=True)
with zipfile.ZipFile(p) as z:
 for n in z.infolist():
  target=(inp/n.filename).resolve();assert target.is_relative_to(inp.resolve())
  if n.is_dir():target.mkdir(parents=True,exist_ok=True);continue
  data=z.read(n);target.parent.mkdir(parents=True,exist_ok=True)
  if target.exists():assert target.read_bytes()==data
  else:target.write_bytes(data)
base=inp/'tongbian_human_PMC9366892_deidentified_inputs_20261012'
verified=[]
for line in (base/'SHA256SUMS.txt').read_text().splitlines():
 if not line.strip():continue
 h,n=line.split(maxsplit=1);n=n.lstrip('*');q=base/n;assert q.resolve().is_relative_to(base.resolve());assert hashlib.sha256(q.read_bytes()).hexdigest()==h;verified.append(n)
(inp/p.name).write_bytes(p.read_bytes())
script=root/'environment/PMC9366892_limma.R';script.write_bytes((staging/script.name).read_bytes())
cfg=root/'00_protocol/PMC9366892_limma_config_frozen_v1.json';cfg.write_bytes((staging/'analysis_config_frozen.json').read_bytes())
(inp/'source_identity.json').write_text(json.dumps({'library_file_id':'libfile_9025b7a56c008191971c58e1dd17241e','version':0,'zip_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'verified_inner_files':verified,'no_public_network_downloads':True},indent=2))
R='/mnt/BioDisk/wangmingyuan/miniconda3/envs/R_4.5.2/bin/Rscript'
env=os.environ.copy();env.update(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1')
z=subprocess.run([R,'--vanilla',str(script),'--check'],capture_output=True,text=True,timeout=45,env=env)
(staging/'R_preflight.log').write_text(z.stdout+'\n'+z.stderr);print(z.stdout,flush=True);print(z.stderr,flush=True);assert z.returncode==0
job=root/'logs/jobs/PMC9366892_limma_v1';job.mkdir(parents=True,exist_ok=True)
if (job/'launch.lock').exists():print('EXISTING_JOB', (job/'status.json').read_text());raise SystemExit(0)
spec={'root':str(root),'job_dir':str(job),'command':[R,'--vanilla',str(script),'--run'],'script_sha256':hashlib.sha256(script.read_bytes()).hexdigest()}
fd=os.open(job/'launch.lock',os.O_CREAT|os.O_EXCL|os.O_WRONLY);os.write(fd,b'launch once');os.close(fd)
(job/'spec.json').write_text(json.dumps(spec,indent=2));(job/'status.json').write_text(json.dumps({'status':'QUEUED','created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat()}))
with (job/'launcher.log').open('a') as f:
 child=subprocess.Popen([sys.executable,str(root/'environment/durable_R_job_runner.py'),str(job/'spec.json')],stdout=f,stderr=subprocess.STDOUT,stdin=subprocess.DEVNULL,start_new_session=True,close_fds=True)
(job/'runner_pid.txt').write_text(str(child.pid));print('LAUNCHED',child.pid,flush=True)
