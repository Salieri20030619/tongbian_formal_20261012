from pathlib import Path
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
