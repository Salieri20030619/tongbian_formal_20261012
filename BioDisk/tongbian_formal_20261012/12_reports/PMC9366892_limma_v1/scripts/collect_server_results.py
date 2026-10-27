from pathlib import Path
import json,hashlib,tarfile,socket,shutil,datetime
r=Path('/mnt/BioDisk/wangmingyuan/tongbian_formal_20261012');assert socket.gethostname()=='osv01'
job=r/'logs/jobs/PMC9366892_limma_v1';status=json.loads((job/'status.json').read_text());print(json.dumps(status),flush=True)
if status['status']!='COMPLETED':
 print((job/'run.log').read_text()[-8000:] if (job/'run.log').exists() else 'NO_LOG');raise SystemExit(2)
out=r/'06_human_selected_metabolites/limma_v1';print((out/'audit/result_summary.txt').read_text(),flush=True)
files=[]
for folder in [out,job,r/'03_inputs/PMC9366892_deidentified_v1']:
 files.extend(p for p in folder.rglob('*') if p.is_file() and '/tmp/' not in str(p))
files.extend([r/'environment/PMC9366892_limma.R',r/'environment/durable_R_job_runner.py',r/'00_protocol/PMC9366892_limma_config_frozen_v1.json',r/'06_human_selected_metabolites/limma_staging_v1/R_preflight.log'])
manifest=[{'path':str(p.relative_to(r)),'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in sorted(files)]
mp=out/'audit/server_output_manifest.json';mp.write_text(json.dumps(manifest,indent=2));files.append(mp)
bundle=r/'12_reports/PMC9366892_limma_v1_server_bundle.tar.gz'
with tarfile.open(bundle,'w:gz') as t:
 for p in files:t.add(p,arcname=str(p.relative_to(r)))
print('BUNDLE',json.dumps({'path':str(bundle),'bytes':bundle.stat().st_size,'sha256':hashlib.sha256(bundle.read_bytes()).hexdigest(),'manifest_files':len(manifest)}),flush=True)
