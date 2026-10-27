from pathlib import Path
import os,sys,json,subprocess,datetime,concurrent.futures,threading,fcntl,traceback
P=Path('/mnt/BioDisk/wangmingyuan/tongbian_formal_20261012')
J=P/'environment/metagenome_prepare_20261020_v1'
D=Path('/mnt/BioData/bio/wangmingyuan/tongbian_formal_20261012')
E=D/'software/envs';CACHE=D/'software/conda_pkgs';TMP=D/'software/tmp'
for p in [E,CACHE,TMP,J/'logs']:p.mkdir(parents=True,exist_ok=True)
CONDA='/mnt/BioDisk/wangmingyuan/miniconda3/bin/conda'
env=dict(os.environ,CONDA_PKGS_DIRS=str(CACHE),CONDA_ENVS_PATH=str(E),TMPDIR=str(TMP),CONDA_REMOTE_CONNECT_TIMEOUT_SECS='15',CONDA_REMOTE_READ_TIMEOUT_SECS='90',CONDA_REMOTE_MAX_RETRIES='3',OMP_NUM_THREADS='2',OPENBLAS_NUM_THREADS='2')
for k in ['http_proxy','https_proxy','all_proxy','HTTP_PROXY','HTTPS_PROXY','ALL_PROXY']:env.pop(k,None)
now=lambda:datetime.datetime.now(datetime.timezone.utc).isoformat()
guard=threading.Lock()
state={'status':'RUNNING','pid':os.getpid(),'started_utc':now(),'scope':'software installation and reference preparation only; no sample profiling','jobs':{}}
def write(p,o):
 t=p.with_suffix(p.suffix+'.tmp');t.write_text(json.dumps(o,indent=2),encoding='utf-8');os.replace(t,p)
def update(n,**kw):
 with guard:
  state['jobs'].setdefault(n,{}).update(kw);state['updated_utc']=now();write(J/'software_status.json',state)
specs={'core':['bowtie2=2.5.5','samtools','megahit=1.2.9','prodigal=2.6.3','diamond','hmmer','seqkit'],'biobakery':['python=3.10','humann=3.9','metaphlan=4.1.1'],'dbcan':['python=3.11','dbcan=5.2.9'],'maaslin2':['bioconductor-maaslin2=1.22.0','r-base=4.5']}
tests={'core':[['bowtie2','--version'],['samtools','--version'],['megahit','--version'],['prodigal','-v'],['diamond','version'],['hmmsearch','-h'],['seqkit','version']],'biobakery':[['humann','--version'],['metaphlan','--version'],['humann_databases','--available']],'dbcan':[['run_dbcan','--help']],'maaslin2':[['Rscript','--vanilla',str(J/'check_maaslin2.R')]]}
(J/'check_maaslin2.R').write_text('library(Maaslin2)\ncat(as.character(packageVersion("Maaslin2")), "\\n")\nsessionInfo()\n',encoding='utf-8')
write(J/'software_plan.json',{'env_root':str(E),'cache':str(CACHE),'specifications':specs,'channels':['https://conda.anaconda.org/conda-forge','https://conda.anaconda.org/bioconda'],'max_parallel_installs':2,'system_environments_modified':False})
def install(n,packages):
 prefix=E/(n+'_20261020_v1')
 try:
  assert not prefix.exists(),'Prefix already exists'
  cmd=[CONDA,'create','--yes','--prefix',str(prefix),'--override-channels','--strict-channel-priority','--solver','libmamba','-c','https://conda.anaconda.org/conda-forge','-c','https://conda.anaconda.org/bioconda']+packages
  update(n,status='INSTALLING',prefix=str(prefix),command=cmd,started_utc=now())
  with (J/'logs'/(n+'_install.log')).open('wb') as log:
   proc=subprocess.Popen(cmd,env=env,stdin=subprocess.DEVNULL,stdout=log,stderr=log)
   update(n,child_pid=proc.pid);rc=proc.wait()
  update(n,install_exit_code=rc)
  assert rc==0,'Conda exited '+str(rc)
  checkenv=dict(env,PATH=str(prefix/'bin')+':'+env.get('PATH',''))
  checks=[]
  for args in tests[n]:
   cmd=[str(prefix/'bin'/args[0])]+args[1:]
   r=subprocess.run(cmd,env=checkenv,capture_output=True,timeout=180)
   (J/'logs'/(n+'_'+args[0]+'_check.log')).write_bytes(r.stdout+r.stderr)
   checks.append({'command':cmd,'exit_code':r.returncode})
   assert r.returncode==0,str(checks[-1])
  for flag,name in [('--explicit','explicit.txt'),('--json','packages.json')]:
   r=subprocess.run([CONDA,'list','--prefix',str(prefix),flag],env=env,capture_output=True,check=True)
   (J/(n+'_'+name)).write_bytes(r.stdout)
  update(n,status='VERIFIED',checks=checks,finished_utc=now())
 except BaseException as ex:update(n,status='FAILED',error=repr(ex),traceback=traceback.format_exc(),finished_utc=now())
lock=(J/'software.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
write(J/'software_status.json',state)
with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
 for f in [pool.submit(install,n,p) for n,p in specs.items()]:f.result()
rc=0 if all(v['status']=='VERIFIED' for v in state['jobs'].values()) else 1
state.update(status='COMPLETED_VERIFIED' if rc==0 else 'FINISHED_WITH_FAILURES',exit_code=rc,finished_utc=now());write(J/'software_status.json',state)
(J/'software_exit_code.txt').write_text(str(rc),encoding='utf-8')
if rc==0:write(J/'SOFTWARE_COMPLETE.json',state)
sys.exit(rc)
