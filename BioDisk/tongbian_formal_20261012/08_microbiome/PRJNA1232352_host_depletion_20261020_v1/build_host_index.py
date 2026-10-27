from pathlib import Path
import os,sys,json,hashlib,datetime,subprocess,traceback,fcntl
P=Path('/mnt/BioDisk/wangmingyuan/tongbian_formal_20261012')
J=P/'08_microbiome/PRJNA1232352_host_depletion_20261020_v1'
D=Path('/mnt/BioData/bio/wangmingyuan/tongbian_formal_20261012')
REF=Path('/mnt/BioData/bio/Databases/STAR_Indexes/genome/GRCh38.p14.genome.fa')
OUT=D/'databases/host_GRCh38p14_20261020_v1'
BIN=D/'software/envs/core_20261020_v1/bin'
now=lambda:datetime.datetime.now(datetime.timezone.utc).isoformat()
def write(p,o):
 t=p.with_suffix(p.suffix+'.tmp');t.write_text(json.dumps(o,indent=2),encoding='utf-8');os.replace(t,p)
lock=(J/'index.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
state={'status':'VALIDATING_REFERENCE','pid':os.getpid(),'started_utc':now()}
write(J/'index_status.json',state)
rc=1
try:
 expected=json.loads((P/'environment/metagenome_prepare_20261020_v1/host_reference_reuse.json').read_text())
 h=hashlib.sha256();counts={};name=None
 with REF.open('rb') as f:
  for line in f:
   h.update(line)
   if line.startswith(b'>'):
    name=line[1:].split()[0].decode('ascii');assert name not in counts;counts[name]=0
   else:
    assert name is not None
    counts[name]+=len(line.strip())
 assert h.hexdigest()==expected['sha256']
 assert counts['chr1']==248956422 and counts['chrM']==16569 and all('chr'+str(i) in counts for i in range(1,23))
 receipt={'path':str(REF),'sha256':h.hexdigest(),'bytes':REF.stat().st_size,'source_readme':'/mnt/BioData/bio/Databases/STAR_Indexes/README','documented_source':'https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/latest_release/GRCh38.p14.genome.fa.gz','source_identity_scope':'shared-file provenance plus fixed local SHA256 and canonical contig-length checks; no claimed upstream byte-for-byte comparison','contig_lengths':counts,'total_bases':sum(counts.values())}
 write(J/'reference_frozen.json',receipt)
 OUT.mkdir(parents=True,exist_ok=False)
 cmd=[str(BIN/'bowtie2-build'),'--threads','8','--seed','20261020',str(REF),str(OUT/'GRCh38p14')]
 state.update(status='BUILDING_INDEX',command=cmd);write(J/'index_status.json',state)
 with (J/'index_build.log').open('wb') as log:
  proc=subprocess.Popen(cmd,stdout=log,stderr=log,stdin=subprocess.DEVNULL,env=dict(os.environ,PATH=str(BIN)+':'+os.environ['PATH'],OMP_NUM_THREADS='8'))
  state['child_pid']=proc.pid;write(J/'index_status.json',state)
  rc=proc.wait()
 assert rc==0,'Index builder exit '+str(rc)
 files=list(OUT.glob('*.bt2'))+list(OUT.glob('*.bt2l'));assert len(files)==6
 with (J/'index_inspect.txt').open('wb') as log:subprocess.run([str(BIN/'bowtie2-inspect'),'-s',str(OUT/'GRCh38p14')],stdout=log,stderr=subprocess.STDOUT,check=True)
 hashes={}
 for f in files:
  hh=hashlib.sha256()
  with f.open('rb') as stream:
   for x in iter(lambda:stream.read(8*1024*1024),b''):hh.update(x)
  hashes[f.name]={'bytes':f.stat().st_size,'sha256':hh.hexdigest()}
 state.update(status='COMPLETED_VERIFIED',index=str(OUT/'GRCh38p14'),files=hashes,exit_code=0,finished_utc=now());rc=0
 write(J/'INDEX_COMPLETE.json',state)
except BaseException as e:
 state.update(status='FAILED',error=repr(e),traceback=traceback.format_exc(),finished_utc=now());rc=1
state['exit_code']=rc;write(J/'index_status.json',state);(J/'index_exit_code.txt').write_text(str(rc),encoding='utf-8')
sys.exit(rc)
