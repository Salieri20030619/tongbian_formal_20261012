from pathlib import Path
import os,sys,json,subprocess,hashlib,re,datetime,concurrent.futures,threading,tarfile,fcntl,traceback
P=Path('/mnt/BioDisk/wangmingyuan/tongbian_formal_20261012')
J=P/'environment/metagenome_prepare_20261020_v1'
D=Path('/mnt/BioData/bio/wangmingyuan/tongbian_formal_20261012/databases/preparation_20261020_v1')
for p in [D/'archives',D/'prepared',J/'database_logs']:p.mkdir(parents=True,exist_ok=True)
now=lambda:datetime.datetime.now(datetime.timezone.utc).isoformat()
guard=threading.Lock()
state={'status':'RUNNING','pid':os.getpid(),'started_utc':now(),'max_downloads':4,'artifacts':{}}
def write(p,v):
 t=p.with_suffix(p.suffix+'.tmp');t.write_text(json.dumps(v,indent=2),encoding='utf-8');os.replace(t,p)
def update(n,**kw):
 with guard:state['artifacts'].setdefault(n,{}).update(kw);state['updated_utc']=now();write(J/'database_status.json',state)
def digests(p):
 a=hashlib.sha256();b=hashlib.md5()
 with p.open('rb') as f:
  for x in iter(lambda:f.read(8*1024*1024),b''):a.update(x);b.update(x)
 return a.hexdigest(),b.hexdigest()
def one(e):
 n=e['id'];part=D/'archives'/(e['filename']+'.partial');dest=D/'archives'/e['filename']
 try:
  assert not dest.exists(),'Archive already exists'
  assert os.statvfs(D).f_bavail*os.statvfs(D).f_frsize>1024**4,'Require 1 TiB free'
  update(n,status='HEAD',url=e['url'],partial=str(part))
  h=subprocess.run(['curl','--noproxy','*','-fsSIL','--connect-timeout','15','--max-time','90','--retry','2',e['url']],capture_output=True,text=True,check=True)
  (J/'database_logs'/(n+'_headers.txt')).write_text(h.stdout,encoding='utf-8')
  lengths=re.findall(r'(?im)^content-length:\s*(\d+)',h.stdout)
  expected=int(lengths[-1]) if lengths else e.get('bytes')
  if e.get('bytes') and expected:assert expected==e['bytes'],(expected,e['bytes'])
  tags=re.findall(r'(?im)^etag:\s*"([0-9a-f]{32})"',h.stdout)
  md5_expected=e.get('md5') or (tags[-1] if e['url'].startswith('https://dbcan.s3.') and tags else None)
  update(n,status='DOWNLOADING',expected_bytes=expected,md5_expected=md5_expected)
  rc=1
  for attempt in range(1,6):
   cmd=['curl','--noproxy','*','--fail','--location','--silent','--show-error','--connect-timeout','20','--speed-limit','1024','--speed-time','600','--continue-at','-','--output',str(part),e['url']]
   with (J/'database_logs'/(n+'.log')).open('ab') as log:
    proc=subprocess.Popen(cmd,stdin=subprocess.DEVNULL,stdout=log,stderr=log)
    update(n,attempt=attempt,child_pid=proc.pid);rc=proc.wait()
   update(n,curl_exit_code=rc,bytes=part.stat().st_size if part.exists() else 0)
   if rc==0:break
  assert rc==0,'curl exit '+str(rc)
  if expected is not None:assert part.stat().st_size==expected,'Size mismatch'
  update(n,status='VALIDATING')
  sha,md5=digests(part)
  if md5_expected:assert md5==md5_expected,'MD5 mismatch; partial retained'
  if e.get('archive'):
   with tarfile.open(part,'r:*') as tar:
    members=tar.getmembers()
    for m in members:
     p=Path(m.name)
     assert not p.is_absolute() and '..' not in p.parts and (m.isfile() or m.isdir()),'Unsafe archive member '+m.name
   if e['filename'].endswith('.gz'):subprocess.run(['gzip','-t',str(part)],check=True,capture_output=True)
  part.rename(dest)
  result={'url':e['url'],'archive':str(dest),'bytes':dest.stat().st_size,'sha256':sha,'md5':md5,'upstream_md5_checked':bool(md5_expected),'verification_scope':'HTTPS + size + archive integrity + local SHA256; upstream MD5 only when available'}
  write(J/'database_logs'/(n+'_receipt.json'),result)
  if e.get('archive'):
   out=D/'prepared'/n
   out.mkdir(exist_ok=False)
   update(n,status='EXTRACTING',output=str(out))
   subprocess.run(['tar','--extract','--file',str(dest),'--directory',str(out),'--no-same-owner','--no-same-permissions'],check=True,capture_output=True)
   count=sum(1 for p in out.rglob('*') if p.is_file());assert count>0
   result.update(extracted=str(out),extracted_files=count)
  else:
   out=D/'prepared'/e.get('group',n);out.mkdir(exist_ok=True)
   link=out/e.get('output_name',e['filename'])
   assert not link.exists();link.symlink_to(dest)
   result['prepared_path']=str(link)
  update(n,status='VERIFIED_PREPARED',finished_utc=now(),**result)
 except BaseException as ex:update(n,status='FAILED',error=repr(ex),traceback=traceback.format_exc(),finished_utc=now())
lock=(J/'databases.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
entries=json.loads((J/'database_manifest.json').read_text())
write(J/'database_status.json',state)
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
 for f in [pool.submit(one,e) for e in entries]:f.result()
rc=0 if all(v['status']=='VERIFIED_PREPARED' for v in state['artifacts'].values()) else 1
state.update(status='DOWNLOADS_PREPARED_PENDING_RUNTIME_VALIDATION' if rc==0 else 'FINISHED_WITH_FAILURES',exit_code=rc,finished_utc=now())
write(J/'database_status.json',state);(J/'database_exit_code.txt').write_text(str(rc),encoding='utf-8')
if rc==0:write(J/'DATABASE_DOWNLOADS_COMPLETE.json',state)
with (P/'HANDOFF.md').open('a',encoding='utf-8') as f:f.write('\n## '+now()+' / Reference package preparation ended\n'+state['status']+'; exit '+str(rc)+'. See '+str(J)+'. Runtime compatibility validation is separate.\n')
sys.exit(rc)
