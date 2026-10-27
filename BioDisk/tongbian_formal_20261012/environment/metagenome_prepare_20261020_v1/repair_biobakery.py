import pathlib,os,subprocess,json,datetime,traceback
P=pathlib.Path('/mnt/BioDisk/wangmingyuan/tongbian_formal_20261012')
J=P/'environment/metagenome_prepare_20261020_v1'
D=pathlib.Path('/mnt/BioData/bio/wangmingyuan/tongbian_formal_20261012')
E=D/'software/envs/biobakery_20261020_v1'
C='/mnt/BioDisk/wangmingyuan/miniconda3/bin/conda'
env=dict(os.environ,CONDA_PKGS_DIRS=str(D/'software/conda_pkgs'),TMPDIR=str(D/'software/tmp'),CONDA_NUMBER_CHANNEL_NOTICES='0',OMP_NUM_THREADS='2',OPENBLAS_NUM_THREADS='2')
now=lambda:datetime.datetime.now(datetime.timezone.utc).isoformat()
record={'status':'RUNNING','pid':os.getpid(),'started_utc':now(),'reason':'Conda HUMAnN3.9 generic build contains python3.12/site-packages; interpreter must match'}
(J/'biobakery_repair_status.json').write_text(json.dumps(record,indent=2),encoding='utf-8')
rc=1
try:
 cmd=[C,'install','--yes','--prefix',str(E),'--override-channels','--strict-channel-priority','--solver','libmamba','-c','https://conda.anaconda.org/conda-forge','-c','https://conda.anaconda.org/bioconda','python=3.12','humann=3.9','metaphlan=4.1.1']
 with (J/'logs/biobakery_python312_repair.log').open('wb') as f:subprocess.run(cmd,env=env,stdout=f,stderr=f,check=True)
 checks=[]
 for args in [['humann','--version'],['metaphlan','--version'],['humann_databases','--available']]:
  cmd=[str(E/'bin'/args[0])]+args[1:]
  r=subprocess.run(cmd,env=dict(env,PATH=str(E/'bin')+':'+env['PATH']),capture_output=True,timeout=180)
  (J/'logs'/('biobakery_'+args[0]+'_recheck.log')).write_bytes(r.stdout+r.stderr)
  assert r.returncode==0,(args,r.stderr.decode())
  checks.append({'command':cmd,'exit_code':r.returncode})
 for flag,name in [('--explicit','explicit.txt'),('--json','packages.json')]:
  r=subprocess.run([C,'list','--prefix',str(E),flag],env=env,capture_output=True,check=True);(J/('biobakery_'+name)).write_bytes(r.stdout)
 record.update(status='VERIFIED',checks=checks);rc=0
except BaseException as e:record.update(status='FAILED',error=repr(e),traceback=traceback.format_exc())
record.update(exit_code=rc,finished_utc=now())
(J/'biobakery_repair_status.json').write_text(json.dumps(record,indent=2),encoding='utf-8')
s=json.loads((J/'software_status.json').read_text())
(J/'software_status_before_python312_repair.json').write_text(json.dumps(s,indent=2),encoding='utf-8')
s['jobs']['biobakery'].update(record)
if rc==0:
 s['jobs']['biobakery'].pop('error',None);s['jobs']['biobakery'].pop('traceback',None)
allrc=0 if all(x['status']=='VERIFIED' for x in s['jobs'].values()) else 1
s.update(status='COMPLETED_VERIFIED' if allrc==0 else 'FINISHED_WITH_FAILURES',exit_code=allrc,finished_utc=now())
(J/'software_status.json').write_text(json.dumps(s,indent=2),encoding='utf-8')
(J/'software_exit_code.txt').write_text(str(allrc),encoding='utf-8')
if allrc==0:(J/'SOFTWARE_COMPLETE.json').write_text(json.dumps(s,indent=2),encoding='utf-8')
with (P/'HANDOFF.md').open('a',encoding='utf-8') as f:f.write('\n## '+now()+' / HUMAnN interpreter correction\n'+record['status']+'; all software status '+s['status']+'. See environment/metagenome_prepare_20261020_v1.\n')
