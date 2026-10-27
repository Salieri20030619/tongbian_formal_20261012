from pathlib import Path
import os,sys,json,gzip,hashlib,datetime,subprocess,csv,math,re,shutil,threading,concurrent.futures,fcntl,traceback
P=Path('/mnt/BioDisk/wangmingyuan/tongbian_formal_20261012')
J=P/'08_microbiome/PRJNA1232352_profiling_20261021_v1'
D=Path('/mnt/BioData/bio/wangmingyuan/tongbian_formal_20261012')
IN=D/'large_intermediates/PRJNA1232352_host_depletion_20261020_v1/samples'
OUT=D/'large_intermediates/PRJNA1232352_profiling_20261021_v1'
E=D/'software/envs/biobakery_20261020_v1'
CORE=D/'software/envs/core_20261020_v1/bin'
DB=json.loads((P/'environment/metagenome_prepare_20261020_v1/database_runtime_validation.json').read_text())['paths']
INDEX='mpa_vJun23_CHOCOPhlAnSGB_202403'
env=dict(os.environ,PATH=str(E/'bin')+':'+str(CORE)+':'+os.environ['PATH'],TMPDIR=str(D/'software/tmp'),OMP_NUM_THREADS='8',OPENBLAS_NUM_THREADS='1')
now=lambda:datetime.datetime.now(datetime.timezone.utc).isoformat()
def write(p,v):
 t=p.with_suffix(p.suffix+'.tmp');t.write_text(json.dumps(v,indent=2),encoding='utf-8');os.replace(t,p)
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
 return h.hexdigest()
cfg=json.loads((J/'config.json').read_text());allrows=json.loads((J/'samples.json').read_text())
def prepare_input(row,out,pilot):
 run=row['run'];count=100000 if pilot else row['retained_pairs']
 source=[IN/run/(run+'_'+str(i)+'.fastq.gz') for i in [1,2]]
 target=out/'combined.fastq.gz'
 assert not target.exists()
 if pilot:
  members=[out/'pilot_R1.fastq.gz',out/'pilot_R2.fastq.gz']
  with gzip.open(source[0],'rb') as a,gzip.open(source[1],'rb') as b,gzip.open(members[0],'wb',compresslevel=3) as x,gzip.open(members[1],'wb',compresslevel=3) as y:
   for _ in range(count):
    r1=[a.readline() for _ in range(4)];r2=[b.readline() for _ in range(4)]
    assert r1[0].startswith(b'@') and r2[0].startswith(b'@') and r1[0].rstrip().removesuffix(b'/1')==r2[0].rstrip().removesuffix(b'/2')
    assert r1[2].startswith(b'+') and r2[2].startswith(b'+')
    x.writelines(r1);y.writelines(r2)
 else:members=source
 h=hashlib.sha256()
 with target.open('xb') as dest:
  for f in members:
   with f.open('rb') as src:
    for data in iter(lambda:src.read(8*1024*1024),b''):dest.write(data);h.update(data)
 r=subprocess.run([str(CORE/'seqkit'),'stats','-T','-j','2',str(target)],env=env,capture_output=True,text=True)
 (out/'input_stats.tsv').write_text(r.stdout,encoding='utf-8');assert r.returncode==0,r.stderr
 stats=list(csv.DictReader(r.stdout.splitlines(),delimiter='\t'));assert len(stats)==1 and int(stats[0]['num_seqs'])==2*count
 write(out/'input_receipt.json',{'source':[str(x) for x in source],'reads':2*count,'pairs':count,'pilot_prefix':pilot,'combined_gzip_members':2,'sha256':h.hexdigest(),'bytes':target.stat().st_size})
 return target,2*count
def validate_mpa(path,expected):
 header=[];rows=[]
 for line in path.read_text().splitlines():
  if not line:continue
  if line.startswith('#'):header.append(line);continue
  f=line.split('\t');assert len(f)>=3
  v=float(f[2]);assert math.isfinite(v) and 0<=v<=100.0001
  rows.append((f[0],v))
 assert any(INDEX in h for h in header),'Wrong marker database'
 counts=[int(m.group(1)) for h in header if (m:=re.search(r'(\d+) reads processed',h))]
 assert counts and counts[-1]==expected,('MetaPhlAn read count',counts,expected)
 assert len({k for k,v in rows})==len(rows)
 species=[(k,v) for k,v in rows if k.split('|')[-1].startswith('s__')]
 assert sum(v for k,v in species)<=100.1
 return {'taxonomic_rows':len(rows),'species_rows':len(species),'species_abundance_sum':sum(v for k,v in species),'reads_processed':counts[-1],'header':header}
def validate_humann_table(path,coverage=False):
 assert path.is_file() and path.stat().st_size>0
 rows=[];seen=set()
 with path.open(encoding='utf-8') as f:
  header=f.readline().strip().split('\t');assert len(header)==2
  for line in f:
   if not line.strip():continue
   fields=line.rstrip().split('\t');assert len(fields)==2 and fields[0] not in seen
   seen.add(fields[0]);v=float(fields[1]);assert math.isfinite(v) and v>=0
   if coverage:assert v<=1.0000001
   rows.append((fields[0],v))
 assert rows,'Missing numeric rows'
 special={'UNMAPPED','UNINTEGRATED','UNCLASSIFIED'}
 return {'rows':len(rows),'unstratified_features':sum('|' not in k and k not in special for k,v in rows),'stratified_features':sum('|' in k for k,v in rows),'positive_features':sum(v>0 and k.split('|')[0] not in special for k,v in rows),'sha256':sha(path),'bytes':path.stat().st_size}
def run_stage(cmd,log,out,tag,update,run):
 write(out/(tag+'_command.json'),cmd)
 with log.open('wb') as f:
  proc=subprocess.Popen(cmd,stdin=subprocess.DEVNULL,stdout=f,stderr=f,env=env)
  update(run,phase=tag,child_pid=proc.pid);rc=proc.wait()
 (out/(tag+'_exit_code.txt')).write_text(str(rc),encoding='utf-8')
 assert rc==0,(tag,rc,str(log))
def one(row,mode,update):
 run=row['run'];pilot=mode=='pilot';out=OUT/mode/run
 try:
  assert os.statvfs(OUT).f_bavail*os.statvfs(OUT).f_frsize>1024**4,'Require 1 TiB free'
  out.mkdir(parents=True,exist_ok=False)
  update(run,status='RUNNING',phase='PREPARE_INPUT',group=row['group'],started_utc=now())
  inputfile,n=prepare_input(row,out,pilot)
  profile=out/'metaphlan.tsv'
  cmd=[str(E/'bin/metaphlan'),str(inputfile),'--input_type','fastq','--bowtie2db',DB['metaphlan'],'--index',INDEX,'--offline','--nproc','8','--read_min_len','50','--bowtie2out',str(out/'metaphlan.bowtie2.bz2'),'--sample_id',run,'-o',str(profile)]
  run_stage(cmd,out/'metaphlan.log',out,'metaphlan',update,run)
  mpa=validate_mpa(profile,n)
  write(out/'METAPHLAN_VALIDATED.json',mpa)
  hdir=out/'humann'
  cmd=[str(E/'bin/humann'),'--input',str(inputfile),'--input-format','fastq.gz','--output',str(hdir),'--output-basename',run,'--threads','8','--taxonomic-profile',str(profile),'--nucleotide-database',DB['chocophlan'],'--protein-database',DB['uniref90'],'--search-mode','uniref90','--prescreen-threshold','0.01','--translated-alignment','diamond','--memory-use','minimum','--output-format','tsv','--o-log',str(out/'humann_internal.log')]
  run_stage(cmd,out/'humann_console.log',out,'humann',update,run)
  update(run,phase='VALIDATING')
  tables={}
  for name in ['genefamilies','pathabundance','pathcoverage']:
   tables[name]=validate_humann_table(hdir/(run+'_'+name+'.tsv'),name=='pathcoverage')
  result={'run':run,'group':row['group'],'reads':n,'pilot':pilot,'metaphlan':mpa,'humann':tables,'completed_utc':now(),'profile_sha256':sha(profile)}
  write(out/'VALIDATED.json',result)
  update(run,status='VERIFIED',phase='COMPLETE',result=result)
  return result
 except BaseException as e:
  update(run,status='FAILED',error=repr(e),traceback=traceback.format_exc());raise
if __name__=='__main__':
 mode=sys.argv[1];assert mode in ['pilot','full']
 lock=(J/(mode+'.lock')).open('a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
 assert not (J/(mode+'_status.json')).exists(),'Do not overwrite previous job'
 if mode=='full':
  gate=json.loads((J/'PILOT_PASS.json').read_text());assert gate['status']=='PASS' and gate['pipeline_sha256']==sha(Path(__file__)) and gate['config_sha256']==sha(J/'config.json')
 rows=[min([r for r in allrows if r['group']==g],key=lambda r:r['retained_pairs']) for g in ['DCTT','NCTT','H']] if mode=='pilot' else allrows
 OUT.mkdir(parents=True,exist_ok=True)
 guard=threading.Lock();state={'status':'RUNNING','mode':mode,'pid':os.getpid(),'started_utc':now(),'total_samples':len(rows),'completed_samples':0,'samples':{}}
 def update(sample_key=None,**kw):
  with guard:
   if sample_key:state['samples'].setdefault(sample_key,{}).update(kw)
   else:state.update(kw)
   state['completed_samples']=sum(v['status']=='VERIFIED' for v in state['samples'].values());state['updated_utc']=now();write(J/(mode+'_status.json'),state)
 update();rc=1
 try:
  with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
   results=list(pool.map(lambda r:one(r,mode,update),rows))
  assert state['completed_samples']==len(rows)
  update(status='COMPLETED_VERIFIED',exit_code=0,finished_utc=now())
  write(J/(mode.upper()+'_COMPLETE.json'),{'status':'COMPLETED_VERIFIED','samples':len(rows),'results':results,'finished_utc':now()})
  rc=0
 except BaseException as e:update(status='FAILED',error=repr(e),traceback=traceback.format_exc(),exit_code=1,finished_utc=now())
 (J/(mode+'_exit_code.txt')).write_text(str(rc),encoding='utf-8')
 sys.exit(rc)
