from pathlib import Path
import csv,json,subprocess,os,datetime,hashlib,time
P=Path('/mnt/BioDisk/wangmingyuan/tongbian_formal_20261012');O=P/'01_inventory/resume_20261013';B=Path('/mnt/BioData/bio/wangmingyuan/tongbian_formal_20261012');D=B/'raw_reads/PRJNA1232352';D.mkdir(exist_ok=True)
runs=list(csv.DictReader((O/'PRJNA1232352_runs.tsv').open(),delimiter='\t'));files=[]
for r in runs:
 urls=r['fastq_ftp'].split(';');md=r['fastq_md5'].split(';');sizes=r['fastq_bytes'].split(';');assert len(urls)==len(md)==len(sizes)==2
 for u,h,n in zip(urls,md,sizes):files.append({'run':r['run_accession'],'sample':r['sample_alias'],'url':'https://'+u,'filename':u.split('/')[-1],'md5':h,'bytes':int(n)})
need=sum(r['bytes'] for r in files);avail=os.statvfs(B).f_bavail*os.statvfs(B).f_frsize;assert avail>need+1_000_000_000_000
(O/'PRJNA1232352_download_manifest.json').write_text(json.dumps(files,indent=2))
first=files[0];part=D/(first['filename']+'.partial');receipt={'time_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'available_bytes':avail,'expected_compressed_bytes':need,'files':len(files),'samples':len(runs),'reserved_free_floor_bytes':1_000_000_000_000,'route':'explicit direct --noproxy *; no proxy variables/listener; no global changes','curl_version':subprocess.check_output(['curl','--version']).decode(),'probe_file':str(part),'attempts':[]}
for i in [1,2]:
 before=part.stat().st_size if part.exists() else 0;start=time.monotonic();headers=O/f'FASTQ_resume_probe_{i}.headers.txt'
 with (O/f'FASTQ_resume_probe_{i}.log').open('w') as log:r=subprocess.run(['curl','--noproxy','*','-fL','--connect-timeout','15','--max-time','30','--limit-rate','4M','-C','-','-D',str(headers),'-o',str(part),first['url']],stdout=log,stderr=log)
 after=part.stat().st_size if part.exists() else 0;hh=headers.read_text() if headers.exists() else '';receipt['attempts'].append({'before_bytes':before,'after_bytes':after,'elapsed_seconds':time.monotonic()-start,'exit_code':r.returncode,'HTTP_status_lines':[s for s in hh.splitlines() if s.startswith('HTTP/')],'content_range':[s for s in hh.splitlines() if s.lower().startswith('content-range')],'bytes_per_second':(after-before)/(time.monotonic()-start)})
 if after==0:break
receipt['resume_confirmed']=len(receipt['attempts'])==2 and receipt['attempts'][1]['before_bytes']>0 and receipt['attempts'][1]['after_bytes']>receipt['attempts'][1]['before_bytes'] and any('206' in s for s in receipt['attempts'][1]['HTTP_status_lines'])
receipt['bulk_gate']=receipt['resume_confirmed'] and receipt['attempts'][-1]['bytes_per_second']>=1_000_000
receipt['bulk_gate_reason']='requires actual resumed large-file response plus >=1MB/s observed to avoid days-long blind transfer on this shared machine; retain partial for next attempt'
(O/'FASTQ_download_preflight.json').write_text(json.dumps(receipt,indent=2));print(json.dumps(receipt,indent=2))
