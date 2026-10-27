from pathlib import Path
import os,sys,json,subprocess,hashlib,datetime,threading,concurrent.futures,fcntl,traceback
P=Path('/mnt/BioDisk/wangmingyuan/tongbian_formal_20261012');B=Path('/mnt/BioData/bio/wangmingyuan/tongbian_formal_20261012');D=B/'raw_reads/PRJNA1232352';J=P/'logs/jobs/PRJNA1232352_direct_download_v1';manifest=P/'01_inventory/resume_20261013/PRJNA1232352_download_manifest.json'
now=lambda:datetime.datetime.now(datetime.timezone.utc).isoformat()
if '--run' not in sys.argv:
 if J.exists():
  print('EXISTING_JOB_DO_NOT_RELAUNCH', (J/'status.json').read_text() if (J/'status.json').exists() else str(J));sys.exit(0)
 J.mkdir(parents=True);(J/'files').mkdir();D.mkdir(parents=True,exist_ok=True)
 cfg={'created_utc':now(),'manifest':str(manifest),'manifest_sha256':hashlib.sha256(manifest.read_bytes()).hexdigest(),'data_root':str(D),'workers':4,'route':'explicit direct --noproxy *','user_instruction':'2026-10-13: background direct download now, slow speed accepted; reconnect proxy tomorrow; no proxy connection started','file_completion':'expected bytes AND ENA MD5 AND gzip -t exit0; atomic rename from .partial','existing_partial':'resume with curl -C -','network_settings_changed':False,'disk_free_floor_bytes':1000000000000,'monitoring':'persistent Python runner waits natively on curl subprocesses; no agent-level polling needed','tomorrow':'inspect status/PIDs before any intervention; completed files reusable; preserve partials; do not run a second writer to same paths'}
 (J/'config.json').write_text(json.dumps(cfg,indent=2));(J/'status.json').write_text(json.dumps({'status':'STARTING','updated_utc':now()}))
 with (J/'launcher.log').open('w') as log:p=subprocess.Popen(['python3',str(Path(__file__).resolve()),'--run'],stdout=log,stderr=subprocess.STDOUT,start_new_session=True,cwd=str(P))
 (J/'runner_pid.txt').write_text(str(p.pid));print(json.dumps({'started':True,'runner_pid':p.pid,'job_dir':str(J),'route':'direct','workers':4}));sys.exit(0)
lock=(J/'run.lock').open('w');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
entries=json.loads(manifest.read_text());guard=threading.Lock();state={'status':'RUNNING','started_utc':now(),'runner_pid':os.getpid(),'total_files':len(entries),'total_expected_bytes':sum(x['bytes'] for x in entries),'route':'direct','files':{}}
def update(name=None,**kw):
 with guard:
  if name:state['files'][name]=dict(state['files'].get(name,{}),**kw)
  else:state.update(kw)
  state['updated_utc']=now();state['completed_files']=sum(v.get('status') in ['VERIFIED','REUSED_VERIFIED'] for v in state['files'].values());state['failed_files']=sum(v.get('status')=='FAILED' for v in state['files'].values());tmp=J/'status.tmp';tmp.write_text(json.dumps(state,indent=2));os.replace(tmp,J/'status.json')
def md5(p):
 h=hashlib.md5()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
 return h.hexdigest()
def one(e):
 n=e['filename'];dest=D/n;part=D/(n+'.partial');logp=J/'files'/(n+'.log')
 try:
  if dest.exists():
   assert dest.stat().st_size==e['bytes'] and md5(dest)==e['md5'],'Existing completed file does not match manifest; preserved'
   update(n,status='REUSED_VERIFIED',bytes=e['bytes'],md5=e['md5']);return True
  assert os.statvfs(B).f_bavail*os.statvfs(B).f_frsize>1000000000000+e['bytes'],'Free-space floor reached'
  cmd=['curl','--noproxy','*','--fail','--location','--show-error','--connect-timeout','30','--retry','5','--retry-delay','30','--retry-all-errors','--speed-time','600','--speed-limit','1024','--continue-at','-','--output',str(part),e['url']]
  with logp.open('a') as log:
   log.write('\nSTART '+now()+'\n'+json.dumps(cmd)+'\n');log.flush();p=subprocess.Popen(cmd,stdout=log,stderr=subprocess.STDOUT);update(n,status='DOWNLOADING',child_pid=p.pid,partial_path=str(part),starting_bytes=part.stat().st_size if part.exists() else 0,expected_bytes=e['bytes'],started_utc=now());rc=p.wait();log.write('\nCURL_EXIT '+str(rc)+' '+now()+'\n');log.flush()
   if rc!=0:raise RuntimeError('curl exit '+str(rc)+'; partial retained')
   assert part.stat().st_size==e['bytes'],'Size mismatch; partial retained'
   actual=md5(part);assert actual==e['md5'],'MD5 mismatch; partial retained'
   gz=subprocess.run(['gzip','-t',str(part)],stdout=log,stderr=log);assert gz.returncode==0,'gzip integrity failure; partial retained'
   part.rename(dest);update(n,status='VERIFIED',bytes=e['bytes'],md5=actual,curl_exit_code=rc,gzip_exit_code=gz.returncode,completed_utc=now());return True
 except BaseException as ex:
  update(n,status='FAILED',error=repr(ex),partial_bytes=part.stat().st_size if part.exists() else None,finished_utc=now());return False
update()
try:
 with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:ok=list(pool.map(one,entries))
 rc=0 if all(ok) else 1;update(status='COMPLETED_VERIFIED' if rc==0 else 'FINISHED_WITH_FAILURES',exit_code=rc,completed_utc=now());(J/'exit_code.txt').write_text(str(rc))
 if rc==0:(J/'DOWNLOAD_COMPLETE.json').write_text(json.dumps({'status':'COMPLETED_VERIFIED','files':len(entries),'bytes':sum(x['bytes'] for x in entries),'completed_utc':now()},indent=2))
 with (P/'HANDOFF.md').open('a') as f:f.write('\n## '+now()+' / 宏基因组直连下载任务结束\n状态：'+state['status']+'；退出码：'+str(rc)+'；已校验文件：'+str(state['completed_files'])+'/'+str(len(entries))+'。详见logs/jobs/PRJNA1232352_direct_download_v1/status.json。\n')
 sys.exit(rc)
except BaseException as ex:
 if isinstance(ex,SystemExit):raise
 update(status='FAILED_RUNNER',exit_code=2,error=repr(ex),traceback=traceback.format_exc());(J/'exit_code.txt').write_text('2');raise
