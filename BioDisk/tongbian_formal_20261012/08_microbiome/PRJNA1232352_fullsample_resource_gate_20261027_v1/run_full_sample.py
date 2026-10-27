import pathlib,json,subprocess,hashlib,datetime,os,csv,traceback,shutil,threading,time,fcntl,io
N=pathlib.Path(__file__).parent;R=N.parent.parent;P=json.loads((N/'plan_frozen.json').read_text());B=pathlib.Path('/mnt/BioData/bio/wangmingyuan/tongbian_formal_20261012');L=pathlib.Path(P['large_output']);C=B/'software/envs/core_20261020_v1/bin';D=B/'software/envs/dbcan_20261020_v1/bin';DB=B/'databases/preparation_20261020_v1/prepared/dbcan_5_2_9';sample=P['sample'];run=sample['run'];phase='initializing';stop=threading.Event();cg=pathlib.Path('/sys/fs/cgroup')/pathlib.Path('/proc/self/cgroup').read_text().strip().split('::')[1].lstrip('/')
lock=(N/'run.lock').open('w');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
def now():return datetime.datetime.now(datetime.timezone.utc).isoformat()
def save(p,x):
 t=p.with_suffix(p.suffix+'.tmp');t.write_text(json.dumps(x,indent=2)+'\n',encoding='utf-8');t.replace(p)
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(8388608),b''):h.update(b)
 return h.hexdigest()
def memavail():return int(next(l.split()[1] for l in pathlib.Path('/proc/meminfo').read_text().splitlines() if l.startswith('MemAvailable:')))*1024
def cgread(name):
 p=cg/name
 return p.read_text().strip() if p.exists() else None
def cgnum(name):
 s=cgread(name);return int(s) if s and s.isdigit() else None
def psival(kind):
 lines=pathlib.Path('/proc/pressure'/pathlib.Path(kind)).read_text().splitlines();return {l.split()[0]:float(l.split()[1].split('=')[1]) for l in lines}
def monitor():
 with (N/'resource_samples.tsv').open('w',encoding='utf-8',newline='') as f:
  w=csv.writer(f,delimiter='\t',lineterminator='\n');w.writerow(['utc','phase','cgroup_memory_current','cgroup_memory_peak_including_cache','MemAvailable_bytes','disk_free_bytes','memory_some_avg10','memory_full_avg10','io_some_avg10','io_full_avg10'])
  while not stop.is_set():
   try:
    m=psival('memory');i=psival('io');w.writerow([now(),phase,cgnum('memory.current'),cgnum('memory.peak'),memavail(),shutil.disk_usage(L).free,m['some'],m['full'],i['some'],i['full']]);f.flush()
   except Exception as e:print('MONITOR_WARNING',str(e),flush=True)
   stop.wait(P['resource_sampling_seconds'])
def state(value,status='RUNNING'):
 global phase
 phase=value;save(N/'status.json',{'status':status,'phase':phase,'run':run,'pid':os.getpid(),'utc':now()});print(now(),status,phase,flush=True)
def execute(cmd,label,guard=True):
 state(label)
 if guard:assert memavail()>=100*1024**3 and shutil.disk_usage(L).free>=300*1024**3,'resource guard before '+label
 env=dict(os.environ,PATH=str(D)+':'+str(C)+':'+os.environ['PATH'],OMP_NUM_THREADS='16',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',LC_ALL='C.UTF-8')
 before={'utc':now(),'memory_events':cgread('memory.events'),'MemAvailable':memavail(),'disk_free':shutil.disk_usage(L).free};save(N/(label+'_command.json'),{'command':cmd,'before':before})
 with (N/(label+'.log')).open('w',encoding='utf-8') as f:
  child=subprocess.Popen(['/usr/bin/time','-v','nice','-n','10']+cmd,stdout=f,stderr=subprocess.STDOUT,env=env);save(N/(label+'_child.json'),{'pid':child.pid,'started_utc':now()});rc=child.wait()
 (N/(label+'_exit_code.txt')).write_text(str(rc)+'\n',encoding='utf-8');assert rc==0,(label,rc)
 save(N/(label+'_resource_end.json'),{'utc':now(),'memory_events':cgread('memory.events'),'cgroup_memory_peak':cgnum('memory.peak'),'output_directory_bytes':int(subprocess.check_output(['du','-sb',str(L)]).split()[0]),'disk_free':shutil.disk_usage(L).free})
def fasta(p):
 ids=set();lens=[];key=None;length=0
 with p.open() as f:
  for line in f:
   if line.startswith('>'):
    if key is not None:assert length>0;lens.append(length)
    key=line[1:].split()[0];assert key not in ids;ids.add(key);length=0
   else:length+=len(line.strip())
 if key is not None:assert length>0;lens.append(length)
 return ids,lens
initial_events=cgread('memory.events');save(N/'resource_baseline.json',{'utc':now(),'cgroup':str(cg),'memory_events':initial_events,'memory_current':cgnum('memory.current'),'memory_peak':cgnum('memory.peak'),'MemAvailable':memavail(),'disk_free':shutil.disk_usage(L).free,'note':'cgroup includes dedicated SSH session and cache; not pure process RSS'})
monitor_thread=threading.Thread(target=monitor,daemon=True);monitor_thread.start()
def main():
 state('input_hash_verification');source=B/'large_intermediates/PRJNA1232352_host_depletion_20261020_v1/samples'/run;files=[source/(run+'_'+str(i)+'.fastq.gz') for i in [1,2]];checks=[]
 for f in files:
  expected=sample['host_output'][f.name];assert f.stat().st_size==expected['bytes'];h=sha(f);assert h==expected['sha256'];checks.append({'file':str(f),'bytes':f.stat().st_size,'sha256':h})
 state('input_read_count');p=subprocess.run([str(C/'seqkit'),'stats','-T','-j','4']+[str(f) for f in files],capture_output=True,text=True);(N/'input_seqkit_stats.tsv').write_text(p.stdout,encoding='utf-8');(N/'input_seqkit_stderr.log').write_text(p.stderr,encoding='utf-8');assert p.returncode==0,p.stderr;stats=list(csv.DictReader(io.StringIO(p.stdout),delimiter='\t'));assert len(stats)==2
 for x in stats:assert int(x['num_seqs'].replace(',',''))==sample['retained_pairs']
 save(N/'INPUT_ACCEPTANCE.json',{'status':'PASS','full_input_pairs':sample['retained_pairs'],'files':checks,'seqkit_stats':stats,'pair_ID_integrity':'preserved upstream acceptance via full-file matching SHA256','prefix_subsample':False})
 old=R/'08_microbiome/PRJNA1232352_functional_pilot_20261027_v1/database_receipt.json';db=json.loads(old.read_text());
 for name,record in db['files'].items():assert (DB/name).stat().st_size==record['bytes'] and sha(DB/name)==record['sha256']
 save(N/'database_receipt.json',db)
 execute([str(C/'megahit'),'-1',str(files[0]),'-2',str(files[1]),'-o',str(L/'megahit'),'-t','16','-m',str(P['megahit_memory_bytes']),'--min-contig-len','500','--k-list',P['megahit_k_list']],'megahit')
 contig=L/'megahit/final.contigs.fa';ids,lens=fasta(contig);assert len(ids)>0 and min(lens)>=500
 lengths=sorted(lens,reverse=True);total=sum(lengths);cum=0;n50=0
 for length in lengths:
  cum+=length
  if cum>=total/2:n50=length;break
 save(N/'ASSEMBLY_ACCEPTANCE.json',{'status':'PASS','contigs':len(ids),'bases':total,'N50':n50,'max_length':max(lens),'min_length':min(lens),'sha256':sha(contig),'utc':now()})
 execute([str(C/'prodigal'),'-i',str(contig),'-a',str(L/'proteins.faa'),'-d',str(L/'genes.fna'),'-o',str(L/'genes.gff'),'-f','gff','-p','meta'],'prodigal')
 proteins,pl=fasta(L/'proteins.faa');genes,gl=fasta(L/'genes.fna');assert proteins==genes and len(proteins)>0
 execute([str(D/'run_dbcan'),'CAZyme_annotation','--mode','protein','--input_raw_data',str(L/'proteins.faa'),'--output_dir',str(L/'dbcan'),'--db_dir',str(DB),'--methods','diamond,hmm','--threads','16','--large','--batch_size','5000','--no-run_signalp','--no-run_deeptmhmm'],'dbcan')
 overview=L/'dbcan/overview.tsv';rows=list(csv.DictReader(overview.open(),delimiter='\t'));assert all(x['Gene ID'] in proteins for x in rows)
 events=lambda s:{line.split()[0]:int(line.split()[1]) for line in s.splitlines()} if s else {}
 before=events(initial_events);after=events(cgread('memory.events'));assert all(after.get(k,0)==before.get(k,0) for k in ['oom','oom_kill'])
 outputs=[contig,L/'proteins.faa',L/'genes.fna',L/'genes.gff',overview];receipt={'status':'COMPLETED_VERIFIED','finished_utc':now(),'run':run,'group':sample['group'],'full_input_pairs':sample['retained_pairs'],'assembly':json.loads((N/'ASSEMBLY_ACCEPTANCE.json').read_text()),'proteins':len(proteins),'dbcan_any_method':len(rows),'dbcan_two_methods':sum(int(x['#ofTools'])>=2 for x in rows),'outputs':{str(f):{'bytes':f.stat().st_size,'sha256':sha(f)} for f in outputs},'cgroup_memory_peak_including_cache':cgnum('memory.peak'),'memory_events_before':before,'memory_events_after':after,'disk_output_bytes':int(subprocess.check_output(['du','-sb',str(L)]).split()[0]),'cohort_expanded':False,'next':'review full-sample resources before cohort scaling; single sample does not prove all-sample resource upper bound','covariate_issue':'TB-COV-001 remains deferred; no covariate inference'}
 save(N/'FULL_SAMPLE_COMPLETE.json',receipt);state('complete','COMPLETED_VERIFIED');(N/'exit_code.txt').write_text('0\n',encoding='utf-8')
if __name__=='__main__':
 try:main()
 except BaseException as e:
  save(N/'failure.json',{'phase':phase,'error':str(e),'utc':now()});state(phase,'FAILED_OR_BLOCKED');(N/'exit_code.txt').write_text('1\n',encoding='utf-8');traceback.print_exc();raise
 finally:stop.set();monitor_thread.join(timeout=20)
