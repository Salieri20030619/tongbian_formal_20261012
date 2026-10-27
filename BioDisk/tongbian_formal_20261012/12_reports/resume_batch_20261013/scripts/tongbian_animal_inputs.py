from pathlib import Path
import subprocess,json,gzip,csv,re,hashlib,datetime,concurrent.futures
P=Path('/mnt/BioDisk/wangmingyuan/tongbian_formal_20261012');O=P/'01_inventory/resume_20261013';A=P/'03_inputs/GSE319612_verified_v1';A.mkdir(exist_ok=True)
tasks={'GSE319612_gene_count.txt.gz':'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE319nnn/GSE319612/suppl/GSE319612_gene_count.txt.gz','MTBLS13907_ftp_index.html':'https://ftp.ebi.ac.uk/pub/databases/metabolights/studies/public/MTBLS13907/','MTBLS13907_study.json':'https://www.ebi.ac.uk/metabolights/ws/studies/MTBLS13907','PRJNA1232352.xml':'https://www.ebi.ac.uk/ena/browser/api/xml/PRJNA1232352','PRJNA1048321.xml':'https://www.ebi.ac.uk/ena/browser/api/xml/PRJNA1048321','PRJNA1031121.xml':'https://www.ebi.ac.uk/ena/browser/api/xml/PRJNA1031121'}
def get(item):
 n,u=item;p=(A if n.startswith('GSE') else O)/n
 if p.exists():return {'file':str(p),'status':'reused'}
 r=subprocess.run(['curl','--noproxy','*','-fLsS','--connect-timeout','10','--max-time','60','--max-filesize','30000000',u],capture_output=True)
 if r.returncode==0:p.write_bytes(r.stdout)
 return {'file':str(p),'url':u,'rc':r.returncode,'bytes':len(r.stdout),'error':r.stderr.decode(errors='replace'),'sha256':hashlib.sha256(r.stdout).hexdigest() if not r.returncode else None}
with concurrent.futures.ThreadPoolExecutor(max_workers=3) as ex:r=list(ex.map(get,tasks.items()))
(O/'animal_followup_retrieval.json').write_text(json.dumps(r,indent=2));print(json.dumps(r,indent=2))
soft=gzip.decompress((O/'GSE319612_family.soft.gz').read_bytes()).decode();samples=[]
for block in soft.split('^SAMPLE = ')[1:]:
 gsm=block.splitlines()[0];lib=re.search(r'!Sample_description = Library name: (.+)',block).group(1).strip();group=re.search(r'!Sample_description = Group: (.+)',block).group(1).strip();title=re.search(r'!Sample_title = (.+)',block).group(1).strip();treat=re.search(r'!Sample_characteristics_ch1 = treatment: (.+)',block).group(1).strip();assert lib.startswith(group)
 assert (group=='NC' and 'normal control' in title) or (group=='MC' and 'model control' in title) or (group=='RG' and 'Rothia' in title)
 samples.append({'GSM':gsm,'sample_id':lib,'group':group,'title':title,'treatment':treat,'group_verified':True,'crossomics_pairing_verified':False})
assert len(samples)==26 and len(set(r['sample_id'] for r in samples))==26
with (A/'sample_metadata.tsv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=list(samples[0]),delimiter='\t');w.writeheader();w.writerows(samples)
(A/'GSE319612_family.soft.gz').write_bytes((O/'GSE319612_family.soft.gz').read_bytes())
p=A/'GSE319612_gene_count.txt.gz'
if p.exists():
 with gzip.open(p,'rt') as f:
  for i,line in enumerate(f):
   if i<4: print('COUNT',line[:3500])
   else:break
print('SAMPLE_MAP',samples)
for n in ['MTBLS13907_ftp_index.html','PRJNA1048321.xml','PRJNA1031121.xml']:
 p=O/n
 if p.exists():print(n,p.read_text()[:10000])
