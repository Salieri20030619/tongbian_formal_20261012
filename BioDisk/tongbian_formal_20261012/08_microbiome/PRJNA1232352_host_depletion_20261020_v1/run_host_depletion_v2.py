from pathlib import Path
import os,sys,json,csv,subprocess,datetime,hashlib,gzip,threading,concurrent.futures,signal,fcntl,traceback,random,re
P=Path('/mnt/BioDisk/wangmingyuan/tongbian_formal_20261012')
J=P/'08_microbiome/PRJNA1232352_host_depletion_20261020_v1'
D=Path('/mnt/BioData/bio/wangmingyuan/tongbian_formal_20261012')
QC=P/'08_microbiome/PRJNA1232352_fastp_20261020_v1'
IN=D/'large_intermediates/PRJNA1232352_fastp_20261020_v1/full'
OUT=D/'large_intermediates/PRJNA1232352_host_depletion_20261020_v1'
BIN=D/'software/envs/core_20261020_v1/bin'
INDEX=D/'databases/host_GRCh38p14_20261020_v1/GRCh38p14'
now=lambda:datetime.datetime.now(datetime.timezone.utc).isoformat()
env=dict(os.environ,PATH=str(BIN)+':'+os.environ['PATH'],OMP_NUM_THREADS='8',OPENBLAS_NUM_THREADS='1',TMPDIR=str(D/'software/tmp'))
lock=(J/'cohort.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
guard=threading.RLock();stop=threading.Event();children={}
state={'status':'WAITING_FOR_INDEX','pid':os.getpid(),'started_utc':now(),'total_samples':60,'completed_samples':0,'samples':{},'max_parallel_samples':4,'alignment_threads_per_sample':8}
def write(p,o):
 t=p.with_suffix(p.suffix+'.tmp');t.write_text(json.dumps(o,indent=2),encoding='utf-8');os.replace(t,p)
def save():
 state['updated_utc']=now();state['completed_samples']=sum(v.get('status')=='VERIFIED' for v in state['samples'].values());write(J/'status.json',state)
def update(sample_key=None,**kw):
 with guard:
  if sample_key:state['samples'].setdefault(sample_key,{}).update(kw)
  else:state.update(kw)
  save()
def sig(sig,frame):
 stop.set()
 with guard:
  for p in children.values():
   if p.poll() is None:p.terminate()
signal.signal(signal.SIGTERM,sig);signal.signal(signal.SIGINT,sig)
def command(args,logkey):
 if stop.is_set():raise InterruptedError('Stop requested')
 with (J/(logkey+'.log')).open('ab') as log:
  p=subprocess.Popen([str(x) for x in args],stdin=subprocess.DEVNULL,stdout=log,stderr=log,env=env)
  with guard:children[logkey]=p
  rc=p.wait()
  with guard:children.pop(logkey,None)
 assert rc==0,(logkey,rc)
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
 return h.hexdigest()
def align_filter(name,r1,r2,out):
 cmd=[str(BIN/'bowtie2'),'--very-sensitive','--end-to-end','--reorder','--seed','20261020','-p','8','-x',str(INDEX),'-1',str(r1),'-2',str(r2)]
 files=[out/(name+'_1.fastq.gz.partial'),out/(name+'_2.fastq.gz.partial')]
 filt=['python3',str(J/'filter_host_sam.py'),str(files[0]),str(files[1]),str(out/'pair_counts.json')]
 write(out/'commands.json',{'bowtie2':cmd,'pair_filter':filt})
 if stop.is_set():raise InterruptedError('Stop requested')
 with (out/'bowtie2.log').open('wb') as log,(out/'filter.log').open('wb') as ferr:
  a=subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=log,env=env,stdin=subprocess.DEVNULL)
  b=subprocess.Popen(filt,stdin=a.stdout,stdout=ferr,stderr=ferr,env=env)
  a.stdout.close()
  with guard:children[name+'_align']=a;children[name+'_filter']=b
  brc=b.wait()
  if brc!=0 and a.poll() is None:a.terminate()
  arc=a.wait()
  with guard:children.pop(name+'_align',None);children.pop(name+'_filter',None)
 write(out/'exit_codes.json',{'aligner':arc,'filter':brc})
 if stop.is_set():raise InterruptedError('Stop requested')
 assert arc==0 and brc==0,(name,arc,brc)
 return files
def validate(name,out,expected):
 counts=json.loads((out/'pair_counts.json').read_text())
 assert counts['total_pairs']==expected,(name,counts['total_pairs'],expected)
 m=re.search(r'(?m)^(\d+) reads; of these:',(out/'bowtie2.log').read_text())
 assert m and int(m.group(1))==expected,'Bowtie2 summary count mismatch'
 files=[out/(name+'_1.fastq.gz.partial'),out/(name+'_2.fastq.gz.partial')]
 for f in files:command([BIN/'gzip' if (BIN/'gzip').exists() else '/usr/bin/gzip','-t',f],name+'_gzipcheck')
 r=subprocess.run([str(BIN/'seqkit'),'stats','-T','-j','2']+[str(x) for x in files],capture_output=True,text=True,env=env)
 (out/'seqkit_stats.tsv').write_text(r.stdout,encoding='utf-8');(out/'seqkit_stats.stderr').write_text(r.stderr,encoding='utf-8')
 assert r.returncode==0,r.stderr
 table=list(csv.DictReader(r.stdout.splitlines(),delimiter='\t'));assert len(table)==2
 assert all(int(x['num_seqs'])==counts['both_unmapped_pairs'] for x in table)
 with gzip.open(files[0],'rb') as a,gzip.open(files[1],'rb') as b:
  for _ in range(min(1000,counts['both_unmapped_pairs'])):
   x=[a.readline() for _ in range(4)];y=[b.readline() for _ in range(4)]
   assert x[0].rstrip().removesuffix(b'/1')==y[0].rstrip().removesuffix(b'/2')
   assert len(x[1].rstrip())==len(x[3].rstrip()) and len(y[1].rstrip())==len(y[3].rstrip())
 hashes={}
 for f in files:
  final=f.with_name(f.name.removesuffix('.partial'));assert not final.exists()
  digest=sha(f);f.rename(final);hashes[final.name]={'bytes':final.stat().st_size,'sha256':digest}
 result={**counts,'run':name,'host_removed_pairs':counts['one_mapped_pairs']+counts['both_mapped_pairs'],'host_removed_fraction':1-counts['both_unmapped_pairs']/expected,'output':hashes,'gzip_check':'PASS','full_read_count_check':'PASS','paired_prefix_check':'PASS','completed_utc':now()}
 write(out/'VALIDATED.json',result);return result
def synthetic():
 out=OUT/'synthetic_control';out.mkdir(parents=True,exist_ok=False)
 ref=json.loads((J/'reference_frozen.json').read_text())['path']
 chunks=[];length=0
 with open(ref,'rb') as f:
  assert f.readline().split()[0]==b'>chr1'
  for line in f:
   if line.startswith(b'>'):break
   chunks.append(line.strip());length+=len(line.strip())
   if length>=1200000:break
 seq=b''.join(chunks)[1000000:]
 rng=random.Random(20261020)
 def randseq():return bytes(rng.choice(b'ACGT') for _ in range(150))
 comp=bytes.maketrans(b'ACGT',b'TGCA')
 reads=[]
 for i in range(10):
  pos=i*1000;left=seq[pos:pos+150];right=seq[pos+300:pos+450].translate(comp)[::-1]
  assert len(left)==150 and len(right)==150 and set(left+right)<=set(b'ACGT')
  reads.extend([(f'human_{i}',left,right),(f'onehost_{i}',left,randseq()),(f'nonhost_{i}',randseq(),randseq())])
 a=out/'input_1.fastq.gz';b=out/'input_2.fastq.gz'
 with gzip.open(a,'wb') as x,gzip.open(b,'wb') as y:
  for name,l,r in reads:
   for f,s in [(x,l),(y,r)]:f.write(b'@'+name.encode()+b'\n'+s+b'\n+\n'+b'I'*len(s)+b'\n')
 align_filter('synthetic',a,b,out);v=validate('synthetic',out,30)
 assert v['both_unmapped_pairs']==10 and v['one_mapped_pairs']==10 and v['both_mapped_pairs']==10,v
 write(J/'SYNTHETIC_CONTROL_PASS.json',v)
def one(row):
 run=row['run_accession']
 if stop.is_set():return
 out=OUT/'samples'/run
 try:
  out.mkdir(parents=True,exist_ok=False)
  expected=qcresults[run]['retained_pairs']
  assert os.statvfs(OUT).f_bavail*os.statvfs(OUT).f_frsize>500*1024**3
  for i in [1,2]:
   f=IN/run/(run+'_'+str(i)+'.fastq.gz')
   assert f.is_file()
  update(run,status='ALIGNING',group=row['group'],input_pairs=expected,started_utc=now())
  align_filter(run,IN/run/(run+'_1.fastq.gz'),IN/run/(run+'_2.fastq.gz'),out)
  update(run,status='VALIDATING')
  result=validate(run,out,expected);result['group']=row['group'];write(out/'VALIDATED.json',result)
  update(run,status='VERIFIED',**result)
 except BaseException as ex:
  update(run,status='INTERRUPTED' if stop.is_set() else 'FAILED',error=repr(ex),traceback=traceback.format_exc())
  stop.set();raise
rc=1
write(J/'status.json',state)
try:
 assert (J/'FILTER_TESTS_PASS.json').exists()
 src=json.loads((QC/'status.json').read_text());assert src['status']=='COMPLETED_VERIFIED_FASTP_ONLY' and src['completed_samples']==60
 qcresults={r['run']:r for r in src['results']}
 rows=json.loads((QC/'samples.json').read_text())
 assert len(rows)==60
 pid=int((J/'index_pid.txt').read_text())
 p=Path('/proc',str(pid),'cmdline')
 if p.exists():
  assert b'build_host_index.py' in p.read_bytes()
  command(['tail','--pid='+str(pid),'-f','/dev/null'],'wait_for_index')
 assert (J/'index_exit_code.txt').read_text().strip()=='0' and (J/'INDEX_COMPLETE.json').exists()
 OUT.mkdir(parents=True,exist_ok=False)
 update(status='SYNTHETIC_VALIDATION');synthetic()
 update(status='FULL_SAMPLE_PILOT');one(rows[0])
 assert state['samples'][rows[0]['run_accession']]['status']=='VERIFIED'
 write(J/'FULL_SAMPLE_PILOT_PASS.json',state['samples'][rows[0]['run_accession']])
 update(status='COHORT_RUNNING')
 with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
  for f in [pool.submit(one,r) for r in rows[1:]]:f.result()
 assert state['completed_samples']==60
 with (J/'host_depletion_summary.tsv').open('w',encoding='utf-8',newline='') as f:
  fields=['run','group','total_pairs','both_unmapped_pairs','one_mapped_pairs','both_mapped_pairs','host_removed_pairs','host_removed_fraction']
  w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',extrasaction='ignore');w.writeheader();w.writerows(state['samples'][r['run_accession']] for r in rows)
 update(status='COMPLETED_VERIFIED',finished_utc=now(),exit_code=0)
 write(J/'HOST_DEPLETION_COMPLETE.json',{'status':'COMPLETED_VERIFIED','samples':60,'finished_utc':now(),'index':str(INDEX),'outputs':str(OUT/'samples'),'downstream_profiling':'NOT_STARTED'})
 rc=0
except BaseException as e:
 sig(None,None)
 update(status='FAILED_OR_INTERRUPTED',error=repr(e),traceback=traceback.format_exc(),finished_utc=now(),exit_code=1)
(J/'exit_code.txt').write_text(str(rc),encoding='utf-8')
with (P/'HANDOFF.md').open('a',encoding='utf-8') as f:f.write('\n## '+now()+' / Host depletion stage ended\n'+state['status']+'; '+str(state['completed_samples'])+'/60; exit '+str(rc)+'. See '+str(J)+'. Sample profiling not started.\n')
sys.exit(rc)
