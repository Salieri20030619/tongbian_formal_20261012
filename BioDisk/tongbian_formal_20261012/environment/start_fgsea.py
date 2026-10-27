from pathlib import Path
import json,os,subprocess,sys,hashlib,datetime,socket
r=Path('/mnt/BioDisk/wangmingyuan/tongbian_formal_20261012');assert socket.gethostname()=='osv01'
assert (r/'05_human_RNA/host_functional_v1/functional_inputs_frozen.rds').exists(),'Functional input preparation not complete; no launch'
summary=json.loads((r/'05_human_RNA/host_functional_v1/functional_input_summary.json').read_text())
assert summary['P_screen_applied']==False and summary['samples']==10
script=r/'environment/run_fgsea_pathways.R'
freeze=json.loads((r/'00_protocol/first_batch_script_freeze.json').read_text())
assert hashlib.sha256(script.read_bytes()).hexdigest()==freeze['sha256']['run_fgsea_pathways.R']
job=r/'logs/jobs/GSE245885_fgsea_v1';job.mkdir(parents=True,exist_ok=True)
if (job/'launch.lock').exists():print('EXISTING_NO_RELAUNCH',(job/'status.json').read_text() if (job/'status.json').exists() else 'locked');raise SystemExit(0)
spec={'job_dir':str(job),'root':str(r),'command':['/mnt/BioDisk/wangmingyuan/miniconda3/envs/R_4.5.2/bin/Rscript','--vanilla',str(script)],'script_sha256':hashlib.sha256(script.read_bytes()).hexdigest()}
fd=os.open(job/'launch.lock',os.O_CREAT|os.O_EXCL|os.O_WRONLY);os.write(fd,b'launch once\n');os.close(fd)
(job/'spec.json').write_text(json.dumps(spec,indent=2));(job/'status.json').write_text(json.dumps({'status':'QUEUED','time_utc':datetime.datetime.now(datetime.timezone.utc).isoformat()}))
with (job/'launcher.log').open('a') as f:
 p=subprocess.Popen([sys.executable,str(r/'environment/durable_R_job_runner.py'),str(job/'spec.json')],stdout=f,stderr=subprocess.STDOUT,stdin=subprocess.DEVNULL,start_new_session=True,close_fds=True)
(job/'runner_pid.txt').write_text(str(p.pid));print('LAUNCHED_FGSEA',p.pid,str(job),flush=True)
