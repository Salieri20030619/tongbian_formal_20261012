from pathlib import Path
import os,json,sys,subprocess,hashlib,socket,datetime
r=Path('/mnt/BioDisk/wangmingyuan/tongbian_formal_20261012');assert socket.gethostname()=='osv01'
R='/mnt/BioDisk/wangmingyuan/miniconda3/envs/R_4.5.2/bin/Rscript'
assert json.loads((r/'logs/jobs/GSE245885_R_DESeq2_bridge_v1/status.json').read_text())['status']=='COMPLETED'
assert json.loads((r/'logs/jobs/MSigDB_package_resource_freeze/status.json').read_text())['status']=='COMPLETED'
refs=json.loads((r/'03_inputs/references/msigdbr_26.1.0_20261012/resource_summary.json').read_text())
assert refs['C2_CP']['gs_count']==4115 and refs['GO_BP']['gs_count']==7538
assert refs['C2_CP']['db_version']=='2026.1.Hs' and refs['GO_BP']['db_version']=='2026.1.Hs'
scripts=['prepare_functional_inputs.R','run_host_functional_scores.R','run_fgsea_pathways.R','sample_quality_audit.R']
env=os.environ.copy();env['OMP_NUM_THREADS']='1';env['OPENBLAS_NUM_THREADS']='1'
for s in scripts:
 p=r/'environment'/s
 ck=subprocess.run([R,'--vanilla','-e',f"invisible(parse(file='{p}'));cat('PARSE_OK {s}\\n')"],capture_output=True,text=True,timeout=30,env=env)
 print(ck.stdout,flush=True);assert ck.returncode==0,ck.stderr
frozen={s:hashlib.sha256((r/'environment'/s).read_bytes()).hexdigest() for s in scripts}
frozen['functional_plan.json']=hashlib.sha256((r/'00_protocol/functional_plan.json').read_bytes()).hexdigest()
manifest=r/'00_protocol/first_batch_script_freeze.json'
if manifest.exists():assert json.loads(manifest.read_text())['sha256']==frozen
else:manifest.write_text(json.dumps({'frozen_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'sha256':frozen},indent=2))
def launch(name,script):
 job=r/'logs/jobs'/name;job.mkdir(parents=True,exist_ok=True)
 if (job/'launch.lock').exists():
  print('EXISTING_NO_RELAUNCH',name,(job/'status.json').read_text() if (job/'status.json').exists() else 'locked',flush=True);return
 p=r/'environment'/script
 spec={'job_dir':str(job),'root':str(r),'command':[R,'--vanilla',str(p)],'script_sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
 fd=os.open(job/'launch.lock',os.O_WRONLY|os.O_CREAT|os.O_EXCL);os.write(fd,b'launch once\n');os.close(fd)
 (job/'spec.json').write_text(json.dumps(spec,indent=2));(job/'status.json').write_text(json.dumps({'status':'QUEUED','time_utc':datetime.datetime.now(datetime.timezone.utc).isoformat()}))
 with (job/'launcher.log').open('a') as f:
  p=subprocess.Popen([sys.executable,str(r/'environment/durable_R_job_runner.py'),str(job/'spec.json')],stdout=f,stderr=subprocess.STDOUT,stdin=subprocess.DEVNULL,start_new_session=True,close_fds=True)
 (job/'runner_pid.txt').write_text(str(p.pid));print('LAUNCHED',name,p.pid,flush=True)
launch('GSE245885_host_functional_scores_v1','run_host_functional_scores.R')
launch('GSE245885_sample_QC_v1','sample_quality_audit.R')
# fgsea is launched in a separate explicit check only after frozen functional inputs exist.
