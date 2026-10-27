from pathlib import Path
import os,sys,json,csv,subprocess,hashlib,datetime,gzip,fcntl,signal,traceback
P=Path('/mnt/BioDisk/wangmingyuan/tongbian_formal_20261012')
J=P/'08_microbiome/PRJNA1232352_fastp_20261020_v1'
D=Path('/mnt/BioData/bio/wangmingyuan/tongbian_formal_20261012')
RAW=D/'raw_reads/PRJNA1232352'
OUT=D/'large_intermediates/PRJNA1232352_fastp_20261020_v1'
now=lambda:datetime.datetime.now(datetime.timezone.utc).isoformat()
def write(p,v):
 t=p.with_suffix(p.suffix+'.tmp');t.write_text(json.dumps(v,indent=2),encoding='utf-8');os.replace(t,p)
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
 return h.hexdigest()
cfg=json.loads((J/'config.json').read_text())
samples=json.loads((J/'samples.json').read_text())
def headcheck(a,b,n=1000):
 with gzip.open(a,'rb') as f,gzip.open(b,'rb') as g:
  count=0
  for i in range(n):
   x=[f.readline() for _ in range(4)];y=[g.readline() for _ in range(4)]
   if not x[0] or not y[0]:
    assert not x[0] and not y[0];break
   assert x[0].startswith(b'@') and y[0].startswith(b'@')
   assert x[0].split()[0].removesuffix(b'/1')==y[0].split()[0].removesuffix(b'/2')
   assert x[2].startswith(b'+') and y[2].startswith(b'+')
   assert len(x[1].rstrip())==len(x[3].rstrip()) and len(y[1].rstrip())==len(y[3].rstrip())
   count+=1
  return count
def command(row,dest,pilot):
 run=row['run_accession']
 return ['/usr/bin/fastp','-i',str(RAW/(run+'_1.fastq.gz')),'-I',str(RAW/(run+'_2.fastq.gz')),'-o',str(dest/(run+'_1.fastq.gz')),'-O',str(dest/(run+'_2.fastq.gz')),'--json',str(dest/'fastp.json'),'--html',str(dest/'fastp.html'),'--thread','4','--compression','4','--detect_adapter_for_pe','--trim_poly_g','--qualified_quality_phred','20','--unqualified_percent_limit','40','--n_base_limit','5','--length_required','50','--dont_overwrite']+(['--reads_to_process','100000'] if pilot else [])
def validate(row,dest,pilot):
 q=json.loads((dest/'fastp.json').read_text())
 before=q['summary']['before_filtering']['total_reads'];after=q['summary']['after_filtering']['total_reads']
 assert before>0 and after>0 and before%2==0 and after%2==0 and after<=before
 if pilot:assert before==200000
 else:assert before==2*int(row['read_count']),(before,row['read_count'])
 for suffix in ['before_filtering','after_filtering']:
  assert q['read1_'+suffix]['total_reads']==q['read2_'+suffix]['total_reads']
  assert 2*q['read1_'+suffix]['total_reads']==q['summary'][suffix]['total_reads']
 run=row['run_accession'];a=dest/(run+'_1.fastq.gz');b=dest/(run+'_2.fastq.gz')
 assert headcheck(a,b)>0
 for f in [a,b]:subprocess.run(['gzip','-t',str(f)],check=True,stdout=subprocess.DEVNULL,stderr=subprocess.PIPE)
 result={'run':run,'group':row['group'],'input_pairs':before//2,'retained_pairs':after//2,'retention':after/before,'q30_before':q['summary']['before_filtering']['q30_rate'],'q30_after':q['summary']['after_filtering']['q30_rate'],'gzip_check':'PASS','paired_prefix_check':'PASS','host_filtering':'NOT_YET_DONE','completed_utc':now()}
 write(dest/'VALIDATED.json',result);return result
if '--pilot' in sys.argv:
 results=[]
 for row in [min([r for r in samples if r['group']==g],key=lambda r:sum(map(int,r['fastq_bytes'].split(';')))) for g in ['DCTT','NCTT','H']]:
  dest=OUT/'pilot'/row['run_accession'];dest.mkdir(parents=True,exist_ok=False)
  cmd=command(row,dest,True);write(dest/'command.json',cmd)
  with (dest/'run.log').open('wb') as f:r=subprocess.run(cmd,stdout=f,stderr=f)
  assert r.returncode==0,r.returncode
  results.append(validate(row,dest,True))
 write(J/'PILOT_COMPLETE.json',{'status':'PASS_TECHNICAL_VALIDATION','first_100000_pairs_per_group_example':True,'representative_QC_claim':False,'results':results})
 print(json.dumps(results));sys.exit(0)
assert '--run' in sys.argv
lock=(J/'run.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
assert (J/'PILOT_COMPLETE.json').exists()
assert sha(Path('/usr/bin/fastp'))==cfg['fastp_sha256']
assert sha(J/'samples.json')==cfg['samples_sha256']
state={'status':'RUNNING','pid':os.getpid(),'started_utc':now(),'completed_samples':0,'total_samples':60,'active':None,'results':[]}
child=None;stopped=False
def stop(sig,frame):
 global stopped
 stopped=True
 if child is not None and child.poll() is None:child.terminate()
signal.signal(signal.SIGTERM,stop);signal.signal(signal.SIGINT,stop)
write(J/'status.json',state)
rc=1
try:
 for row in samples:
  if stopped:raise InterruptedError('User stop')
  assert os.statvfs(OUT).f_bavail*os.statvfs(OUT).f_frsize>500*1024**3,'Free space floor 500 GiB'
  run=row['run_accession'];dest=OUT/'full'/run;dest.mkdir(parents=True,exist_ok=False)
  state.update(active=run,updated_utc=now());write(J/'status.json',state)
  cmd=command(row,dest,False);write(dest/'command.json',cmd)
  with (dest/'run.log').open('wb') as f:
   child=subprocess.Popen(cmd,stdout=f,stderr=f,stdin=subprocess.DEVNULL)
   state['child_pid']=child.pid;write(J/'status.json',state)
   resultcode=child.wait()
  (dest/'exit_code.txt').write_text(str(resultcode),encoding='utf-8')
  if stopped:raise InterruptedError('User stop')
  assert resultcode==0,(run,resultcode)
  state.update(phase='VALIDATING');write(J/'status.json',state)
  result=validate(row,dest,False)
  result['output_sha256']={f.name:sha(f) for f in [dest/(run+'_1.fastq.gz'),dest/(run+'_2.fastq.gz')]}
  write(dest/'VALIDATED.json',result)
  state['results'].append(result);state.update(completed_samples=len(state['results']),phase='FASTP',updated_utc=now());write(J/'status.json',state)
 state.update(status='COMPLETED_VERIFIED_FASTP_ONLY',active=None,finished_utc=now());rc=0
 write(J/'FASTP_COMPLETE.json',{'samples':60,'status':state['status'],'host_filtered':False,'finished_utc':now()})
except BaseException as e:
 state.update(status='INTERRUPTED' if stopped else 'FAILED',error=repr(e),traceback=traceback.format_exc(),finished_utc=now());rc=130 if stopped else 1
state['exit_code']=rc;write(J/'status.json',state);(J/'exit_code.txt').write_text(str(rc),encoding='utf-8')
with (P/'HANDOFF.md').open('a',encoding='utf-8') as f:f.write('\n## '+now()+' / Metagenome fastp stage ended\n'+state['status']+'; completed '+str(state['completed_samples'])+'/60; exit '+str(rc)+'. Job '+str(J)+'. Host depletion and profiling remain pending.\n')
sys.exit(rc)
