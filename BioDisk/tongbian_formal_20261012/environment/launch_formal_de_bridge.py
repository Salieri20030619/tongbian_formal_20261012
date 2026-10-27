from pathlib import Path
import os,json,sys,subprocess,datetime,hashlib,socket
root=Path('/mnt/BioDisk/wangmingyuan/tongbian_formal_20261012');assert socket.gethostname()=='osv01'
R='/mnt/BioDisk/wangmingyuan/miniconda3/envs/R_4.5.2/bin/Rscript'
script=root/'environment/GSE245885_R_DESeq2_bridge.R'
assert hashlib.sha256(script.read_bytes()).hexdigest()=='70bf0e9586da074f36643f128b796e78902183c05f79a068896400eb3a5c6bee'
job=root/'logs/jobs/GSE245885_R_DESeq2_bridge_v1';job.mkdir(parents=True,exist_ok=True)
status=job/'status.json'
if (job/'launch.lock').exists():
 print('EXISTING_JOB_NO_RELAUNCH',status.read_text() if status.exists() else 'launch lock exists; inspect',flush=True);raise SystemExit(0)
assert not (root/'05_human_RNA/R_DESeq2_bridge_v1').exists(),'Existing formal output: inspect before starting'
runner=root/'environment/durable_R_job_runner.py'
code='''from pathlib import Path
import os,json,sys,subprocess,datetime,traceback
spec=json.loads(Path(sys.argv[1]).read_text());job=Path(spec['job_dir'])
def update(**kw):
 kw['updated_utc']=datetime.datetime.now(datetime.timezone.utc).isoformat()
 tmp=job/'status.tmp';tmp.write_text(json.dumps(kw,indent=2));os.replace(tmp,job/'status.json')
env=os.environ.copy()
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS']:env[k]='1'
env['TZ']='UTC';env['TMPDIR']=str(job/'tmp');Path(env['TMPDIR']).mkdir(exist_ok=True)
update(status='RUNNING',runner_pid=os.getpid(),command=spec['command'],script_sha256=spec['script_sha256'])
try:
 with (job/'run.log').open('a') as log:
  p=subprocess.Popen(spec['command'],cwd=spec['root'],stdout=log,stderr=subprocess.STDOUT,env=env)
  update(status='RUNNING',runner_pid=os.getpid(),child_pid=p.pid,command=spec['command'],script_sha256=spec['script_sha256'])
  rc=p.wait()
 update(status='COMPLETED' if rc==0 else 'FAILED',exit_code=rc,runner_pid=os.getpid(),command=spec['command'],script_sha256=spec['script_sha256'])
except BaseException as e:
 update(status='FAILED',exception=repr(e),traceback=traceback.format_exc(),runner_pid=os.getpid())
 raise
'''
if runner.exists():assert runner.read_text()==code
else:runner.write_text(code)
spec={'job_dir':str(job),'root':str(root),'command':[R,'--vanilla',str(script),'--run'],'script_sha256':hashlib.sha256(script.read_bytes()).hexdigest()}
lock=os.open(job/'launch.lock',os.O_CREAT|os.O_EXCL|os.O_WRONLY);os.write(lock,b'one launch only; inspect status before recovery\n');os.close(lock)
(job/'spec.json').write_text(json.dumps(spec,indent=2));status.write_text(json.dumps({'status':'QUEUED','created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat()},indent=2))
cfgpath=root/'00_protocol/project_config.json';cfg=json.loads(cfgpath.read_text());cfg['execution'].update(approved_for_compute=True,formal_analysis_started=True,current_authorized_scope='GSE245885 small-cohort DESeq2 bridge, fgsea, real ssGSEA, VIPER and limma scores; no macro/MD')
cfgpath.write_text(json.dumps(cfg,ensure_ascii=False,indent=2))
with (job/'launcher.log').open('a') as log:
 p=subprocess.Popen([sys.executable,str(runner),str(job/'spec.json')],stdin=subprocess.DEVNULL,stdout=log,stderr=subprocess.STDOUT,start_new_session=True,close_fds=True)
(job/'runner_pid.txt').write_text(str(p.pid))
print(json.dumps({'launched':True,'runner_pid':p.pid,'job_dir':str(job),'status_file':str(status),'result_dir':str(root/'05_human_RNA/R_DESeq2_bridge_v1'),'repeat_launch_protected':True},indent=2),flush=True)
