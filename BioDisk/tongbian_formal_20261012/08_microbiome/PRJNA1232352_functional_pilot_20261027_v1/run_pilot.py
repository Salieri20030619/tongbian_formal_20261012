import pathlib,json,subprocess,gzip,hashlib,datetime,os,traceback,csv
N=pathlib.Path(__file__).parent; R=N.parent.parent; B=pathlib.Path('/mnt/BioData/bio/wangmingyuan/tongbian_formal_20261012'); C=B/'software/envs/core_20261020_v1/bin'; D=B/'software/envs/dbcan_20261020_v1/bin'; DB=B/'databases/preparation_20261020_v1/prepared/dbcan_5_2_9'; P=json.loads((N/'plan_frozen.json').read_text()); L=pathlib.Path(P['large_output'])
def now():return datetime.datetime.now(datetime.timezone.utc).isoformat()
def save(p,x):
 t=p.with_suffix(p.suffix+'.tmp');t.write_text(json.dumps(x,indent=2)+'\n');t.replace(p)
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(8388608),b''):h.update(b)
 return h.hexdigest()
def fasta(p):
 ids=set();lengths=[];seq='';key=None
 with p.open() as f:
  for line in f:
   if line.startswith('>'):
    if key is not None:assert seq;lengths.append(len(seq))
    key=line[1:].split()[0];assert key not in ids;ids.add(key);seq=''
   else:seq+=line.strip()
 if key is not None:assert seq;lengths.append(len(seq))
 return ids,lengths
results=[]
def execute(cmd,label,o):
 env=dict(os.environ,PATH=str(D)+':'+str(C)+':'+os.environ['PATH'],OMP_NUM_THREADS='4',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')
 save(o/(label+'_command.json'),{'command':cmd,'started_utc':now(),'env_thread_limit':4})
 with (o/(label+'.log')).open('w') as f:p=subprocess.run(['/usr/bin/time','-v','nice','-n','10']+cmd,stdout=f,stderr=subprocess.STDOUT,env=env)
 (o/(label+'_exit_code.txt')).write_text(str(p.returncode)+'\n');assert p.returncode==0,(label,p.returncode)
def main():
 helptext=subprocess.check_output([str(D/'run_dbcan'),'CAZyme_annotation','--help'],text=True,env=dict(os.environ,COLUMNS='200',NO_COLOR='1'));(N/'dbcan_help.txt').write_text(helptext)
 save(N/'database_receipt.json',{'files':{f:{'bytes':(DB/f).stat().st_size,'sha256':sha(DB/f)} for f in ['CAZy.dmnd','dbCAN.hmm']},'release':'dbcan 5.2.9 prepared reference','prior_runtime_receipt':str(R/'environment/metagenome_prepare_20261020_v1/database_runtime_validation.json')})
 for sample in P['samples']:
  run=sample['run'];o=L/run;o.mkdir(exist_ok=False);save(N/'status.json',{'status':'RUNNING','run':run,'phase':'subset','completed':len(results),'utc':now()})
  source=B/'large_intermediates/PRJNA1232352_host_depletion_20261020_v1/samples'/run
  sourcefiles=[source/(run+'_'+str(i)+'.fastq.gz') for i in [1,2]]
  for f in sourcefiles:assert f.stat().st_size==sample['host_output'][f.name]['bytes']
  outs=[o/'pilot_1.fastq.gz',o/'pilot_2.fastq.gz'];count=0
  with gzip.open(sourcefiles[0],'rb') as a,gzip.open(sourcefiles[1],'rb') as b,gzip.open(outs[0],'wb',compresslevel=1) as x,gzip.open(outs[1],'wb',compresslevel=1) as y:
   for j in range(P['pairs_per_sample']):
    reads=[[z.readline() for k in range(4)] for z in [a,b]]
    for z in reads:assert z[0].startswith(b'@') and z[2].startswith(b'+') and len(z[1].strip())==len(z[3].strip()) and z[3]
    key=lambda z:z[0].split()[0].removesuffix(b'/1').removesuffix(b'/2')
    assert key(reads[0])==key(reads[1]);x.write(b''.join(reads[0]));y.write(b''.join(reads[1]));count+=1
  for f in outs:
   with gzip.open(f,'rb') as z:assert sum(1 for _ in z)==4*count
  save(o/'INPUT_ACCEPTANCE.json',{'pairs':count,'paired_IDs_match':True,'gzip_full_read_pass':True,'source_files':[str(f) for f in sourcefiles],'source_receipt':sample,'pilot_files':{f.name:sha(f) for f in outs},'scope':'prefix technical subset only'})
  save(N/'status.json',{'status':'RUNNING','run':run,'phase':'megahit','completed':len(results),'utc':now()})
  execute([str(C/'megahit'),'-1',str(outs[0]),'-2',str(outs[1]),'-o',str(o/'megahit'),'-t','4','-m',str(P['megahit_memory_bytes']),'--min-contig-len','500','--k-list',P['megahit_k_list']],'megahit',o)
  contig=o/'megahit/final.contigs.fa';ids,lengths=fasta(contig);assert len(ids)>0 and min(lengths)>=500
  save(N/'status.json',{'status':'RUNNING','run':run,'phase':'prodigal','completed':len(results),'utc':now()})
  execute([str(C/'prodigal'),'-i',str(contig),'-a',str(o/'proteins.faa'),'-d',str(o/'genes.fna'),'-o',str(o/'genes.gff'),'-f','gff','-p','meta'],'prodigal',o)
  proteins,pl=fasta(o/'proteins.faa');genes,gl=fasta(o/'genes.fna');assert proteins==genes and len(proteins)>0
  save(N/'status.json',{'status':'RUNNING','run':run,'phase':'dbcan','completed':len(results),'utc':now()})
  execute([str(D/'run_dbcan'),'CAZyme_annotation','--mode','protein','--input_raw_data',str(o/'proteins.faa'),'--output_dir',str(o/'dbcan'),'--db_dir',str(DB),'--methods','diamond,hmm','--threads','4','--no-run_signalp','--no-run_deeptmhmm'],'dbcan',o)
  overview=o/'dbcan/overview.tsv';assert overview.exists();rows=list(csv.DictReader(overview.open(),delimiter='\t'))
  if rows:
   gene_col=next(k for k in rows[0] if k.lower().replace(' ','').replace('_','') in ['geneid','gene']);assert all(x[gene_col] in proteins for x in rows)
  files=[contig,o/'proteins.faa',o/'genes.fna',o/'genes.gff',overview]
  out={'run':run,'group':sample['group'],'status':'TECHNICAL_PASS','pairs':count,'contigs':len(ids),'assembled_bp':sum(lengths),'proteins':len(proteins),'dbcan_overview_rows':len(rows),'outputs':{str(f.relative_to(o)):{'bytes':f.stat().st_size,'sha256':sha(f)} for f in files},'finished_utc':now(),'biological_validation':False}
  save(o/'VALIDATED.json',out);results.append(out);save(N/'partial_results.json',results);print('VERIFIED',run,len(ids),len(proteins),len(rows),flush=True)
 save(N/'PILOT_COMPLETE.json',{'status':'TECHNICAL_PILOT_PASS','samples':3,'results':results,'finished_utc':now(),'full_cohort_started':False,'biological_conclusions':'not assessed from shallow prefix pilot','next':'full-sample assembly resource gate before cohort scale'})
 save(N/'status.json',{'status':'COMPLETED_VERIFIED','completed':3,'utc':now()});(N/'exit_code.txt').write_text('0\n')
if __name__=='__main__':
 try:main()
 except BaseException as e:
  save(N/'status.json',{'status':'FAILED','completed':len(results),'error':str(e),'utc':now()});(N/'exit_code.txt').write_text('1\n');traceback.print_exc();raise
