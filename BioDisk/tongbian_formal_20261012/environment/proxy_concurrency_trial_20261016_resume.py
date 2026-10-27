
from pathlib import Path
import os,signal,time,json,subprocess,datetime,fcntl,traceback
P=Path('/mnt/BioDisk/wangmingyuan/tongbian_formal_20261012');D=Path('/mnt/BioData/bio/wangmingyuan/tongbian_formal_20261012/raw_reads/PRJNA1232352')
OLD=P/'logs/jobs/PRJNA1232352_proxy_dual_20261016_v1';J=P/'logs/jobs/proxy_concurrency_trial_20261016_v1';J.mkdir(exist_ok=True)
def save(n,x):(J/n).write_text(json.dumps(x,indent=2),encoding='utf-8')
def snap():
 out={}
 for f in D.rglob('*'):
  if f.is_file():
   s=f.stat();out[s.st_ino]=s.st_size
 return out
def observe(job,label,count):
 rows=[];a=snap();t=time.monotonic()
 for i in range(count):
  time.sleep(30);b=snap();nt=time.monotonic();delta=sum(max(0,v-a.get(k,0)) for k,v in b.items());s=json.loads((job/'status.json').read_text());row={'phase':label,'seconds':nt-t,'bytes':delta,'MB_s':delta/(nt-t)/1e6,'completed':s['completed_files'],'failed':s['failed_files'],'status':s['status']};rows.append(row);print(json.dumps(row),flush=True);save(label+'.json',rows);a=b;t=nt
 return sum(r['bytes'] for r in rows)/sum(r['seconds'] for r in rows)/1e6
def stop(job,script):
 s=json.loads((job/'status.json').read_text());pid=s['runner_pid'];assert script.encode() in Path('/proc',str(pid),'cmdline').read_bytes();os.kill(pid,signal.SIGTERM)
 for i in range(45):
  if (job/'exit_code.txt').exists() and not Path('/proc',str(pid)).exists():break
  time.sleep(1)
 assert (job/'exit_code.txt').exists() and not Path('/proc',str(pid)).exists(),'Old runner did not stop'
 for proc in Path('/proc').iterdir():
  if not proc.name.isdigit():continue
  try:c=(proc/'cmdline').read_bytes()
  except OSError:continue
  assert not(c.startswith(b'curl\x00') and str(D).encode() in c),'Another writer remains'
 lock=(P/'logs/jobs/PRJNA1232352_direct_download_v1/run.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB);lock.close()
def launch(previous,jobname,workers):
 source=(P/'environment/dual_proxy_download_20261016.py').read_text()
 source=source.replace("OLD=P/'logs/jobs/PRJNA1232352_direct_dual_20261014_v2'","OLD=P/'logs/jobs/"+previous.name+"'")
 source=source.replace("J=P/'logs/jobs/PRJNA1232352_proxy_dual_20261016_v1'","J=P/'logs/jobs/"+jobname+"'").replace("'workers':2","'workers':"+str(workers)).replace('max_workers=2','max_workers='+str(workers))
 source=source.replace("['large','small']",repr(['large','small']+['extra'+str(i) for i in range(workers-2)]))
 source=source.replace('User requested proxy downloads on current intranet and throughput comparison','User requested increasing concurrent downloads; controlled throughput trial')
 dest=P/'environment'/(jobname+'.py');compile(source,str(dest),'exec');dest.write_text(source,encoding='utf-8')
 r=subprocess.run(['python3',str(dest),'--launch'],capture_output=True,text=True);print(r.stdout,r.stderr,flush=True);assert r.returncode==0
 return P/'logs/jobs'/jobname
try:
 rows=json.loads((J/'two_before.json').read_text());baseline=sum(r['bytes'] for r in rows)/sum(r['seconds'] for r in rows)/1e6
 assert (OLD/'exit_code.txt').read_text()=='130'
 four=launch(OLD,'PRJNA1232352_proxy_four_20261016_v1',4)
 speed=observe(four,'four',4)
 chosen=four;workers=4
 if speed < baseline*1.15:
  stop(four,'PRJNA1232352_proxy_four_20261016_v1.py');chosen=launch(four,'PRJNA1232352_proxy_dual_20261016_v2',2);workers=2
 result={'status':'COMPLETED','two_MB_s':baseline,'four_MB_s':speed,'ratio':speed/baseline if baseline else None,'selected_workers':workers,'active_job':str(chosen),'selection_rule':'Retain four only when aggregate throughput improves at least 15%; short sequential trial, network variability remains','utc':datetime.datetime.now(datetime.timezone.utc).isoformat()}
 save('result.json',result)
 with (P/'HANDOFF.md').open('a',encoding='utf-8') as f:f.write('\n## '+result['utc']+' / Proxy concurrency trial\n'+json.dumps(result)+'\nAll prior verified files and partials retained; exact size, ENA MD5 and gzip checks remain required.\n')
 (J/'resume_exit_code.txt').write_text('0',encoding='utf-8');print(json.dumps(result),flush=True)
except BaseException:
 (J/'resume_exit_code.txt').write_text('1',encoding='utf-8');(J/'error.txt').write_text(traceback.format_exc(),encoding='utf-8');raise
