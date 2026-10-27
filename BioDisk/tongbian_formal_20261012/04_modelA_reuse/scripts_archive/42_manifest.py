"""Write a SHA256 manifest of deliverables, raw responses and reproducible source code."""
from pathlib import Path
import hashlib,pandas as pd
ROOT=Path(__file__).resolve().parents[1]
rows=[]
for dirname in ['raw','tables','figures','scripts','logs','results','output','wheelhouse']:
 for p in sorted((ROOT/dirname).rglob('*')):
  if not p.is_file() or '__pycache__' in p.parts or p.name=='final_sha256_manifest.tsv':continue
  h=hashlib.sha256()
  with p.open('rb') as f:
   while x:=f.read(1024*1024):h.update(x)
  rows.append({'relative_path':p.relative_to(ROOT).as_posix(),'bytes':p.stat().st_size,'sha256':h.hexdigest()})
pd.DataFrame(rows).to_csv(ROOT/'logs/final_sha256_manifest.tsv',sep='\t',index=False,encoding='utf-8-sig')
print('Manifest files:',len(rows),'bytes:',sum(x['bytes'] for x in rows))
