from pathlib import Path
import os,json,subprocess,hashlib,datetime,fcntl,traceback
P=Path('/mnt/BioDisk/wangmingyuan/tongbian_formal_20261012')
J=P/'08_microbiome/PRJNA1232352_profiling_20261021_v1'
D=Path('/mnt/BioData/bio/wangmingyuan/tongbian_formal_20261012/large_intermediates/PRJNA1232352_profiling_20261021_v1')
now=lambda:datetime.datetime.now(datetime.timezone.utc).isoformat()
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
 return h.hexdigest()
def write(p,v):
 t=p.with_suffix(p.suffix+'.tmp');t.write_text(json.dumps(v,indent=2),encoding='utf-8');os.replace(t,p)
lock=(J/'launch_full.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
freeze=json.loads((J/'launch_freeze.json').read_text())
state={'status':'WAITING_FOR_PILOT','pid':os.getpid(),'started_utc':now()}
write(J/'gate_status.json',state);rc=1
try:
 pid=int((J/'pilot_pid.txt').read_text());p=Path('/proc',str(pid),'cmdline')
 if p.exists():
  assert b'profile_pipeline.py' in p.read_bytes() and b'pilot' in p.read_bytes()
  subprocess.run(['tail','--pid='+str(pid),'-f','/dev/null'],check=True)
 assert (J/'pilot_exit_code.txt').read_text().strip()=='0','Pilot failed; full launch withheld'
 s=json.loads((J/'pilot_status.json').read_text());assert s['status']=='COMPLETED_VERIFIED' and s['completed_samples']==3
 results=json.loads((J/'PILOT_COMPLETE.json').read_text())['results']
 assert len(results)==3 and {x['group'] for x in results}=={'DCTT','NCTT','H'}
 assert sha(J/'profile_pipeline.py')==freeze['pipeline_sha256'] and sha(J/'config.json')==freeze['config_sha256']
 assert sha(J/'samples.json')==freeze['samples_sha256']
 checked=[]
 for r in results:
  out=D/'pilot'/r['run']
  assert r['reads']==200000 and r['metaphlan']['reads_processed']==200000
  assert json.loads((out/'input_receipt.json').read_text())['reads']==200000
  assert sha(out/'metaphlan.tsv')==r['profile_sha256']
  for name,v in r['humann'].items():assert sha(out/'humann'/(r['run']+'_'+name+'.tsv'))==v['sha256']
  for name in ['metaphlan','humann']:assert (out/(name+'_exit_code.txt')).read_text().strip()=='0'
  checked.append({'run':r['run'],'group':r['group'],'reads':r['reads'],'species_rows':r['metaphlan']['species_rows'],'gene_features':r['humann']['genefamilies']['unstratified_features'],'pathway_features':r['humann']['pathabundance']['unstratified_features'],'all_output_hashes_match':True})
 write(J/'PILOT_PASS.json',{'status':'PASS','utc':now(),**freeze,'results':checked,'biological_significance_required':False})
 assert not (J/'full_status.json').exists() and not (J/'full_pid.txt').exists(),'Full job already exists'
 with (J/'full_launcher.log').open('ab') as log:
  proc=subprocess.Popen(['nice','-n','10','python3',str(J/'profile_pipeline.py'),'full'],stdin=subprocess.DEVNULL,stdout=log,stderr=log,start_new_session=True)
 (J/'full_pid.txt').write_text(str(proc.pid),encoding='utf-8')
 state.update(status='PILOT_PASSED_FULL_LAUNCHED',full_pid=proc.pid,verified_pilot=checked,finished_utc=now(),exit_code=0);rc=0
except BaseException as e:state.update(status='FULL_LAUNCH_WITHHELD',error=repr(e),traceback=traceback.format_exc(),finished_utc=now(),exit_code=1)
write(J/'gate_status.json',state);(J/'gate_exit_code.txt').write_text(str(rc),encoding='utf-8')
(P/'00_protocol/profiling_state_20261021.json').write_text(json.dumps({'updated_utc':now(),'job':str(J),**state},indent=2),encoding='utf-8')
with (P/'HANDOFF.md').open('a',encoding='utf-8') as f:f.write('\n## '+now()+' / Profiling pilot gate\n'+json.dumps(state)+'. Full profiling completion remains separate; use full_status.json and FULL_COMPLETE.json with exit code.\n')
