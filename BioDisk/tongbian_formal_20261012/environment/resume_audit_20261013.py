from pathlib import Path
import json,hashlib,csv,os,shutil,subprocess
P=Path('/mnt/BioDisk/wangmingyuan/tongbian_formal_20261012')
def read(p):
 print('\nFILE',str(p.relative_to(P))); print(p.read_text())
for p in [P/'01_inventory/server_staged_input_verification.json',*sorted((P/'logs/jobs').glob('*/status.json'))]:read(p)
for p in sorted(P.glob('03_inputs/references/**/*manifest*.json')):read(p)
for p in sorted(P.glob('05_human_RNA/R_DESeq2_bridge_v1/*json')):read(p)
for p in sorted(P.glob('06_human_selected_metabolites/limma_v1/*json')):read(p)
I=P/'03_inputs/PMC9366892_deidentified_v1/tongbian_human_PMC9366892_deidentified_inputs_20261012'
for p in I.glob('*.csv'):
 with p.open() as f: r=list(csv.reader(f))
 print('\nTABLE',p.name,'rows',len(r)-1,'columns',r[0]); print('FIRST',r[1:3])
 if 'conflict' in p.name: print(r)
for p in [P/'environment/PMC9366892_limma.R',P/'environment/durable_R_job_runner.py']:read(p)
print('\nNEXT MODULE FILES')
for name in ['02_metadata','07_animal_multiomics','08_microbiome','09_modelB']:
 for p in sorted((P/name).rglob('*')):
  if p.is_file():print(p.relative_to(P),p.stat().st_size)
print('\nLARGE INPUTS')
B=Path('/mnt/BioData/bio/wangmingyuan/tongbian_formal_20261012')
for p in B.rglob('*'):
 if p.is_file():print(p,p.stat().st_size)
