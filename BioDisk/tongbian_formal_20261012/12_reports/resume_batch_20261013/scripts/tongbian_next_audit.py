from pathlib import Path
import os,json,csv,hashlib,shutil,subprocess,datetime,urllib.parse,concurrent.futures,socket
P=Path('/mnt/BioDisk/wangmingyuan/tongbian_formal_20261012');O=P/'01_inventory/resume_20261013';O.mkdir(exist_ok=True)
ref=[]
for a in json.loads((P/'03_inputs/references/first_batch_reference_manifest.json').read_text()):
 p=Path(a['path']);actual=hashlib.sha256(p.read_bytes()).hexdigest();assert actual==a['sha256'];ref.append({'file':str(p),'sha256':actual,'passed':True})
p=P/'03_inputs/references/dorothea_mm_ABC_v1.22.0.tsv'
if p.exists():
 h=hashlib.sha256(p.read_bytes()).hexdigest();assert h=='3ce829885fea67f8216aa2e33bbabb430607d94eebe818d8f1b9d3b21f158bd7';ref.append({'file':str(p),'sha256':h,'passed':True})
(O/'reference_hash_verification.json').write_text(json.dumps(ref,indent=2))
toolnames=['curl','wget','aria2c','prefetch','fasterq-dump','metaphlan','humann','megahit','prodigal','diamond','run_dbcan','conda']
tools={n:shutil.which(n) for n in toolnames}
envroot=Path('/mnt/BioDisk/wangmingyuan/miniconda3/envs')
tools['environment_names']=[p.name for p in envroot.iterdir() if p.is_dir()]
tools['environment_binaries']={n:[str(p) for p in envroot.glob('*/bin/'+n)] for n in toolnames}
tools['disk']={str(p):shutil.disk_usage(p)._asdict() for p in [P,Path('/mnt/BioData/bio/wangmingyuan/tongbian_formal_20261012')]}
tools['proxy_env']={k:urllib.parse.urlsplit(v)._replace(netloc=urllib.parse.urlsplit(v).hostname or '').geturl() for k,v in os.environ.items() if k.lower() in ['http_proxy','https_proxy','all_proxy']}
try:
 s=socket.create_connection(('127.0.0.1',7897),2);s.close();proxy_available=True
except OSError:proxy_available=False
tools['loopback_7897_listening']=proxy_available
(O/'capacity_tools_network.json').write_text(json.dumps(tools,indent=2))
print(json.dumps(tools,indent=2),flush=True)
tasks={}
fields='run_accession,sample_accession,secondary_sample_accession,study_accession,sample_alias,sample_title,scientific_name,library_strategy,library_layout,library_source,instrument_platform,read_count,base_count,fastq_bytes,fastq_md5,fastq_ftp'
for acc in ['PRJNA1232352','PRJNA1048321','PRJNA1031121','PRJNA1348869']:
 tasks[acc+'_runs.tsv']='https://www.ebi.ac.uk/ena/portal/api/filereport?'+urllib.parse.urlencode({'accession':acc,'result':'read_run','fields':fields,'format':'tsv','download':'true'})
tasks['GSE319612_family.soft.gz']='https://ftp.ncbi.nlm.nih.gov/geo/series/GSE319nnn/GSE319612/soft/GSE319612_family.soft.gz'
tasks['GSE319612_suppl_index.html']='https://ftp.ncbi.nlm.nih.gov/geo/series/GSE319nnn/GSE319612/suppl/'
tasks['MTBLS13907_files.json']='https://www.ebi.ac.uk/metabolights/ws/studies/MTBLS13907/files'
def fetch(item):
 name,url=item;p=O/name
 if p.exists() and p.stat().st_size>0:return {'file':name,'status':'reused','bytes':p.stat().st_size}
 attempts=[]
 for route in ['direct']+(['existing_loopback_proxy'] if proxy_available else []):
  args=['curl','--fail','--location','--silent','--show-error','--connect-timeout','10','--max-time','45','--max-filesize','20000000']
  args+=['--noproxy','*'] if route=='direct' else ['--proxy','http://127.0.0.1:7897','--noproxy','']
  r=subprocess.run(args+[url],capture_output=True);attempts.append({'route':route,'exit_code':r.returncode,'error':r.stderr.decode(errors='replace')[:500]})
  if r.returncode==0:
   p.write_bytes(r.stdout);return {'file':name,'url':url,'status':'downloaded','bytes':len(r.stdout),'sha256':hashlib.sha256(r.stdout).hexdigest(),'attempts':attempts}
 return {'file':name,'url':url,'status':'FAILED','attempts':attempts}
with concurrent.futures.ThreadPoolExecutor(max_workers=3) as ex:results=list(ex.map(fetch,tasks.items()))
(O/'retrieval_receipt.json').write_text(json.dumps(results,indent=2))
for x in results: print(json.dumps(x),flush=True)
summary=[]
for p in O.glob('*_runs.tsv'):
 rows=list(csv.DictReader(p.open(),delimiter='\t'))
 if not rows or 'run_accession' not in rows[0]:continue
 summary.append({'accession':p.name.split('_')[0],'runs':len(rows),'samples':len(set(r['sample_accession'] for r in rows)),'strategies':sorted(set(r['library_strategy'] for r in rows)),'fastq_bytes':sum(sum(int(x) for x in r['fastq_bytes'].split(';') if x) for r in rows),'missing_fastq_rows':sum(not r['fastq_ftp'] for r in rows)})
(O/'run_manifest_summary.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary,indent=2))
