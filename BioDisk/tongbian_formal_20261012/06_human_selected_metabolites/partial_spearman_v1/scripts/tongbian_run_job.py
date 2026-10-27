from pathlib import Path
import subprocess,os,sys,json,datetime,hashlib
P=Path('/mnt/BioDisk/wangmingyuan/tongbian_formal_20261012');name=sys.argv[1];script=P/'environment'/sys.argv[2]
J=P/'logs/jobs'/name;J.mkdir(exist_ok=False)
cmd=['/mnt/BioDisk/wangmingyuan/miniconda3/envs/R_4.5.2/bin/Rscript','--vanilla',str(script)]
spec={'root':str(P),'job_dir':str(J),'command':cmd,'script_sha256':hashlib.sha256(script.read_bytes()).hexdigest()}
(J/'spec.json').write_text(json.dumps(spec,indent=2))
with (J/'launcher.log').open('w') as log:
 p=subprocess.Popen(['python3',str(P/'environment/durable_R_job_runner.py'),str(J/'spec.json')],stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
(J/'runner_pid.txt').write_text(str(p.pid));print('STARTED',name,'PID',p.pid)
rc=p.wait();print('RUNNER_EXIT',rc);print((J/'status.json').read_text());print((J/'run.log').read_text()[-8000:])
sys.exit(0 if json.loads((J/'status.json').read_text()).get('exit_code')==0 else 1)
