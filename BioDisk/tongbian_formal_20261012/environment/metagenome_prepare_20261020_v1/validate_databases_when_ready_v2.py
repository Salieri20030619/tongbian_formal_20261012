from pathlib import Path
import os,json,subprocess,datetime,traceback
P=Path('/mnt/BioDisk/wangmingyuan/tongbian_formal_20261012')
J=P/'environment/metagenome_prepare_20261020_v1'
K=J/'database_v2'
D=Path('/mnt/BioData/bio/wangmingyuan/tongbian_formal_20261012')
R=D/'databases/preparation_20261020_v1/prepared'
now=lambda:datetime.datetime.now(datetime.timezone.utc).isoformat()
def write(p,o):
 t=p.with_suffix(p.suffix+'.tmp');t.write_text(json.dumps(o,indent=2),encoding='utf-8');os.replace(t,p)
state={'status':'WAITING_FOR_DOWNLOADS','pid':os.getpid(),'started_utc':now(),'checks':[]}
write(J/'database_runtime_validation.json',state)
rc=1
try:
 for pid,expected in [(2877159,'database_v2/prepare_databases.py'),(2877826,'repair_biobakery.py')]:
  p=Path('/proc',str(pid),'cmdline')
  if p.exists():
   assert expected.encode() in p.read_bytes()
   subprocess.run(['tail','--pid='+str(pid),'-f','/dev/null'],check=True)
 assert (K/'database_exit_code.txt').read_text().strip()=='0','Some database downloads failed'
 assert (J/'software_exit_code.txt').read_text().strip()=='0','Software not fully verified'
 core=D/'software/envs/core_20261020_v1/bin'
 env=dict(os.environ,PATH=str(core)+':'+os.environ['PATH'],OMP_NUM_THREADS='2')
 def check(name,cmd):
  r=subprocess.run([str(x) for x in cmd],env=env,capture_output=True,timeout=600)
  (J/'database_logs'/(name+'_runtime.log')).write_bytes(r.stdout+r.stderr)
  state['checks'].append({'name':name,'command':[str(x) for x in cmd],'exit_code':r.returncode})
  assert r.returncode==0,name
 def link(src,target):
  if target.is_symlink():assert target.resolve()==src.resolve()
  else:assert not target.exists();target.symlink_to(src)
 state['status']='VALIDATING';write(J/'database_runtime_validation.json',state)
 index='mpa_vJun23_CHOCOPhlAnSGB_202403'
 combined=R/('metaphlan_'+index);combined.mkdir(exist_ok=True)
 for src in list((R/'mpa_metadata').rglob('*'))+list((R/'mpa_bt2').rglob('*')):
  if src.is_file():link(src,combined/src.name)
 assert (combined/(index+'.pkl')).is_file()
 assert len(list(combined.glob(index+'*.bt2l')))==6
 check('metaphlan_bowtie2_index',[core/'bowtie2-inspect','-s',combined/index])
 protein=list((R/'uniref90').rglob('*.dmnd'));assert protein
 for i,db in enumerate(protein):check('uniref90_'+str(i),[core/'diamond','dbinfo','--db',db])
 nucleotide=list((R/'chocophlan').rglob('*.ffn.gz'));assert len(nucleotide)>100
 dbcan=R/'dbcan_5_2_9'
 for f in dbcan.glob('*.dmnd'):check('dbcan_'+f.stem,[core/'diamond','dbinfo','--db',f])
 for f in dbcan.glob('*.hmm'):check('dbcan_'+f.stem,[core/'hmmstat',f])
 for src in (R/'dbcan_dbCAN_PUL_tar_gz').iterdir():link(src,dbcan/src.name)
 state.update(status='REFERENCE_PACKAGES_READY',paths={'metaphlan':str(combined),'chocophlan':str(nucleotide[0].parent),'uniref90':str(protein[0].parent),'utility_mapping':str(R/'humann_mapping'),'dbcan':str(dbcan)},chocophlan_genome_files=len(nucleotide),human_host_index='NOT_BUILT',sample_analysis_started=False)
 rc=0
except BaseException as e:state.update(status='FAILED_OR_BLOCKED',error=repr(e),traceback=traceback.format_exc())
state.update(exit_code=rc,finished_utc=now());write(J/'database_runtime_validation.json',state)
(J/'database_validation_exit_code.txt').write_text(str(rc),encoding='utf-8')
if rc==0:write(J/'REFERENCE_PACKAGES_READY.json',state)
with (P/'HANDOFF.md').open('a',encoding='utf-8') as f:f.write('\n## '+now()+' / Database runtime validation\n'+state['status']+'; exit '+str(rc)+'. See '+str(J)+'/database_runtime_validation.json. Human host index is separate; no sample profiling launched.\n')
