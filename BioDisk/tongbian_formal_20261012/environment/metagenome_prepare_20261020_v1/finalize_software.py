from pathlib import Path
import os,sys,json,subprocess,datetime,hashlib,traceback
P=Path('/mnt/BioDisk/wangmingyuan/tongbian_formal_20261012')
J=P/'environment/metagenome_prepare_20261020_v1'
D=Path('/mnt/BioData/bio/wangmingyuan/tongbian_formal_20261012')
now=lambda:datetime.datetime.now(datetime.timezone.utc).isoformat()
def write(p,o):
 t=p.with_suffix(p.suffix+'.tmp');t.write_text(json.dumps(o,indent=2),encoding='utf-8');os.replace(t,p)
def waitpid(pid,expected):
 p=Path('/proc',str(pid),'cmdline')
 if p.exists():
  assert expected.encode() in p.read_bytes()
  subprocess.run(['tail','--pid='+str(pid),'-f','/dev/null'],check=True)
waitpid(2867265,'prepare_software.py')
state=json.loads((J/'software_status.json').read_text())
write(J/'software_status_before_entrypoint_recheck.json',state)
envbase=dict(os.environ,CONDA_PKGS_DIRS=str(D/'software/conda_pkgs'),TMPDIR=str(D/'software/tmp'),OMP_NUM_THREADS='2',OPENBLAS_NUM_THREADS='2')
n='dbcan';job=state['jobs'][n];prefix=Path(job['prefix'])
if job.get('install_exit_code')==0:
 try:
  cmd=[str(prefix/'bin/run_dbcan'),'--help']
  r=subprocess.run(cmd,env=dict(envbase,PATH=str(prefix/'bin')+':'+os.environ.get('PATH','')),capture_output=True,timeout=180)
  (J/'logs/dbcan_run_dbcan_check.log').write_bytes(r.stdout+r.stderr)
  assert r.returncode==0,r.returncode
  for flag,name in [('--explicit','explicit.txt'),('--json','packages.json')]:
   r=subprocess.run(['/mnt/BioDisk/wangmingyuan/miniconda3/bin/conda','list','--prefix',str(prefix),flag],env=envbase,capture_output=True,check=True)
   (J/(n+'_'+name)).write_bytes(r.stdout)
  job.update(status='VERIFIED',validation_correction='dbCAN5.2.9 official executable is run_dbcan; old wrong-entrypoint result retained in snapshot',checks=[{'command':cmd,'exit_code':0}],finished_utc=now())
  job.pop('error',None);job.pop('traceback',None)
 except BaseException as e:job.update(status='FAILED',error=repr(e),traceback=traceback.format_exc())
rc=0 if all(v['status']=='VERIFIED' for v in state['jobs'].values()) else 1
state.update(status='COMPLETED_VERIFIED' if rc==0 else 'FINISHED_WITH_FAILURES',exit_code=rc,finished_utc=now())
write(J/'software_status.json',state);(J/'software_exit_code.txt').write_text(str(rc),encoding='utf-8')
if rc==0:write(J/'SOFTWARE_COMPLETE.json',state)
ref=Path('/mnt/BioData/bio/Databases/STAR_Indexes/genome/GRCh38.p14.genome.fa')
h=hashlib.sha256()
with ref.open('rb') as f:
 for x in iter(lambda:f.read(8*1024*1024),b''):h.update(x)
write(J/'host_reference_reuse.json',{'path':str(ref),'bytes':ref.stat().st_size,'sha256':h.hexdigest(),'reference':'existing shared GRCh38.p14 FASTA','copied':False,'bowtie2_index':'NOT_BUILT','reference_source_metadata_review':'required before host-depletion run','utc':now()})
with (P/'HANDOFF.md').open('a',encoding='utf-8') as f:f.write('\n## '+now()+' / Software preparation validated\n'+state['status']+'; exit '+str(rc)+'. See '+str(J)+'/software_status.json. Databases downloaded separately; do not equate installed software with completed analysis.\n')
(J/'validation_exit_code.txt').write_text(str(rc),encoding='utf-8')
sys.exit(rc)
