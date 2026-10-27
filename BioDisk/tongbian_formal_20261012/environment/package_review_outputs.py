from pathlib import Path
import zipfile,json,csv
r=Path('/mnt/BioDisk/wangmingyuan/tongbian_formal_20261012')
include=[]
for p in (r/'05_human_RNA').rglob('*'):
 if p.suffix in ['.png','.json'] or (p.suffix=='.tsv' and any(s in p.name for s in ['limma_ALL','sample_QC_all','fgsea_ALL','coverage','PCA_','bridge_testability'])):include.append(p)
include.extend([r/'00_protocol/functional_plan.json',r/'03_inputs/references/first_batch_reference_manifest.json'])
z=r/'12_reports/first_batch_review_outputs.zip'
with zipfile.ZipFile(z,'w',zipfile.ZIP_DEFLATED) as f:
 for p in include:
  if p.is_file():f.write(p,str(p.relative_to(r)))
print('BUNDLE',str(z),'bytes',z.stat().st_size)
print('FGSEA_SUMMARY',(r/'05_human_RNA/fgsea_v1/fgsea_summary.json').read_text())
for c in ['C2_CP','GO_BP']:
 with (r/'05_human_RNA/fgsea_v1'/f'{c}_fgsea_ALL_results.tsv').open() as f:a=list(csv.DictReader(f,delimiter='\t'))
 a=sorted(a,key=lambda x:float(x['pval']))
 print('TOP',c,json.dumps([{k:z[k] for k in ['pathway','pval','BH_within_collection','NES','size','log2err']} for z in a[:12]]))
