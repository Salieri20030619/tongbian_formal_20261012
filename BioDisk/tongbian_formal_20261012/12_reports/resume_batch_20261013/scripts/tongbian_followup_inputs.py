from pathlib import Path
import subprocess,json,csv,re,hashlib,datetime,concurrent.futures,os,gzip
P=Path('/mnt/BioDisk/wangmingyuan/tongbian_formal_20261012');O=P/'01_inventory/resume_20261013';I=P/'03_inputs/MTBLS13907_audit_v1';I.mkdir(exist_ok=True)
base='https://ftp.ebi.ac.uk/pub/databases/metabolights/studies/public/MTBLS13907/'
names=re.findall(r'href="([^"]+)"',(O/'MTBLS13907_ftp_index.html').read_text());tasks={n:base+n for n in names if n.endswith(('.txt','.tsv'))}
tasks['FILES_index.html']=base+'FILES/'
tasks['GSE319612_gene_count.txt.gz']='https://ftp.ncbi.nlm.nih.gov/geo/series/GSE319nnn/GSE319612/suppl/GSE319612_gene_count.txt.gz'
def get(item):
 n,u=item;dest=(P/'03_inputs/GSE319612_verified_v1' if n.startswith('GSE') else I)/n;part=Path(str(dest)+'.partial')
 if dest.exists():return {'file':str(dest),'status':'reused'}
 with (O/(n+'.download.log')).open('a') as log:
  r=subprocess.run(['curl','--noproxy','*','-fL','--connect-timeout','20','--max-time','600','--retry','2','--retry-delay','3','--max-filesize','30000000','-C','-','-o',str(part),u],stdout=log,stderr=log)
 if r.returncode==0:part.rename(dest)
 return {'file':str(dest),'url':u,'rc':r.returncode,'bytes':dest.stat().st_size if dest.exists() else None,'sha256':hashlib.sha256(dest.read_bytes()).hexdigest() if dest.exists() else None}
with concurrent.futures.ThreadPoolExecutor(max_workers=3) as ex:
 results=[]
 for r in ex.map(get,tasks.items()):results.append(r);print(json.dumps(r),flush=True)
(O/'followup_input_receipt.json').write_text(json.dumps(results,indent=2))
a=list(csv.DictReader((O/'PRJNA1048321_runs.tsv').open(),delimiter='\t'));b=list(csv.DictReader((O/'PRJNA1031121_runs.tsv').open(),delimiter='\t'))
overlap=[]
for x in a:
 for y in b:
  md=set(x['fastq_md5'].split(';'))&set(y['fastq_md5'].split(';'))
  if md:overlap.append({'study1_run':x['run_accession'],'study1_alias':x['sample_alias'],'study2_run':y['run_accession'],'study2_alias':y['sample_alias'],'shared_md5':';'.join(sorted(md))})
(O/'16S_exact_file_overlap.json').write_text(json.dumps(overlap,indent=2));print('16S_IDENTICAL_FILE_OVERLAP',len(overlap))
for p in I.iterdir():
 if p.suffix in ['.txt','.tsv']:
  with p.open() as f:rows=list(csv.reader(f,delimiter='\t'))
  print('MT_TABLE',p.name,'ROWS',len(rows),'HEADER',rows[0]);print('FIRST',rows[1:3])
 elif p.suffix=='.html':print(p.read_text()[:6000])
p=P/'03_inputs/GSE319612_verified_v1/GSE319612_gene_count.txt.gz'
if p.exists():
 with gzip.open(p,'rt') as f:
  for i,line in enumerate(f):
   if i<3:print('COUNT',line[:4000])
   else:break
(O/'capacity_current_available.json').write_text(json.dumps({str(p):{'available_bytes':os.statvfs(p).f_bavail*os.statvfs(p).f_frsize,'snapshot_utc':datetime.datetime.now(datetime.timezone.utc).isoformat()} for p in [P,Path('/mnt/BioData/bio/wangmingyuan/tongbian_formal_20261012')]},indent=2))
