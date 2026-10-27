from pathlib import Path
import json,os
r=Path('/mnt/BioDisk/wangmingyuan/tongbian_formal_20261012')
for name in ['GSE245885_R_DESeq2_bridge_v1','MSigDB_package_resource_freeze','GSE245885_host_functional_scores_v1','GSE245885_sample_QC_v1','GSE245885_fgsea_v1']:
 job=r/'logs/jobs'/name
 print('\nJOB',name)
 if not (job/'status.json').exists():print('NOT_STARTED');continue
 s=json.loads((job/'status.json').read_text());print(json.dumps(s))
 if s['status'] in ['RUNNING','FAILED']:
  for k in ['runner_pid','child_pid']:
   if k in s:print(k,s[k],'proc_exists',Path('/proc',str(s[k])).exists())
  if (job/'run.log').exists():print('LOG_TAIL\n'+'\n'.join((job/'run.log').read_text(errors='replace').splitlines()[-16:]))
for rel in ['05_human_RNA/host_functional_v1/functional_input_summary.json','05_human_RNA/host_functional_v1/scores_status.json','05_human_RNA/sample_QC_v1/sample_QC_summary.json','05_human_RNA/fgsea_v1/collection_execution_status.json']:
 p=r/rel
 if p.exists():print('\nOUTPUT',rel,p.read_text())
