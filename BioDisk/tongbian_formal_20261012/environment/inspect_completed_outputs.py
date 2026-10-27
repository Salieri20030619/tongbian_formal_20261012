from pathlib import Path
import csv,json,hashlib,math
r=Path('/mnt/BioDisk/wangmingyuan/tongbian_formal_20261012')
def rows(p):
 with (r/p).open() as f:return list(csv.DictReader(f,delimiter='\t'))
for p in ['05_human_RNA/host_functional_v1/ssGSEA_four_prespecified_modules_limma_ALL_results.tsv','05_human_RNA/host_functional_v1/ssGSEA_module_coverage.tsv','05_human_RNA/sample_QC_v1/sample_QC_all_samples.tsv','05_human_RNA/sample_QC_v1/PCA_top500_variable_VST_scores.tsv']:
 print('\nTABLE',p,json.dumps(rows(p),ensure_ascii=False))
t=rows('05_human_RNA/host_functional_v1/VIPER_all_TFs_limma_ALL_results.tsv')
print('\nTF_TOP20',json.dumps(sorted(t,key=lambda x:float(x['P_value']))[:20]))
print('\nTF_NOT_ELIGIBLE',json.dumps([z for z in rows('05_human_RNA/host_functional_v1/VIPER_regulon_coverage_and_return_status.tsv') if z['eligible']!='TRUE']))
for n in ['GSE245885_host_functional_scores_v1','GSE245885_sample_QC_v1','GSE245885_fgsea_v1']:
 j=r/'logs/jobs'/n
 print('\nJOB',n,(j/'status.json').read_text())
 print('LOG',(j/'run.log').read_text(errors='replace')[-7000:])
print('\nRESOURCE_LOADER',(r/'03_inputs/references/msigdbr_26.1.0_20261012/package_resource_loader_provenance.txt').read_text())
manifest=[]
for p in sorted((r/'03_inputs/references/msigdbr_26.1.0_20261012').rglob('*')):
 if p.is_file():manifest.append(dict(path=str(p),bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest()))
for rel in ['03_inputs/references/dorothea_hs_ABC_v1.22.0.tsv','03_inputs/GSE245885_pilot/genesets/host_modules_frozen.gmt']:
 p=r/rel;manifest.append(dict(path=str(p),bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest()))
out=r/'03_inputs/references/first_batch_reference_manifest.json';out.write_text(json.dumps(manifest,indent=2))
print('\nRESOURCE_HASHES',json.dumps(manifest))
