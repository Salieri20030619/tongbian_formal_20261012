from direct_http import *
import pandas as pd,re
cache=ROOT/'raw'/'pubchem_properties';cache.mkdir(exist_ok=True)
src=pd.read_csv(ROOT/'tables'/'17_herb_ingredient_details.tsv',sep='\t',dtype=str).fillna('')
props='InChIKey,IsomericSMILES,CanonicalSMILES,MolecularFormula,InChI'
cids=sorted({x for value in src.get('PubChem id',pd.Series(dtype=str)) for x in re.findall(r'\b\d+\b',value)},key=int)
todo=[c for c in cids if not (cache/f'{c}.json').exists()];failures=[]
def fetch(vals):
 label='pubchem_cid_'+hashlib.sha256(','.join(vals).encode()).hexdigest()[:16]
 url='https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/'+','.join(vals)+'/property/'+props+'/JSON'
 r,resp=get(url,label,timeout=(12,40))
 if r.get('status') in [401,403,429]:raise SystemExit('Access/rate block: stop')
 data=resp.json() if resp is not None and 'json' in resp.headers.get('Content-Type','') else {}
 if 'PropertyTable' in data:
  for rec in data['PropertyTable']['Properties']:
   rec.update(retrieval_raw=f'raw/{label}.bin',retrieval_url=url)
   (cache/f"{rec['CID']}.json").write_text(json.dumps(rec,ensure_ascii=False),encoding='utf-8')
 elif len(vals)>1 and r.get('status') in [400,404]:
  half=len(vals)//2;fetch(vals[:half]);fetch(vals[half:])
 else:failures.append({'cids':','.join(vals),'status':r.get('status'),'raw':f'raw/{label}.bin'})
 time.sleep(.3)
for i in range(0,len(todo),50):fetch(todo[i:i+50]);print('new HERB CID',min(i+50,len(todo)),len(todo),flush=True)
cached=[json.loads(p.read_text(encoding='utf-8')) for p in cache.glob('*.json')]
pd.DataFrame(cached).to_csv(ROOT/'tables'/'22_pubchem_structure_properties.tsv',sep='\t',index=False,encoding='utf-8-sig')
with (ROOT/'logs'/'herb_pubchem_cid_batches.jsonl').open('a',encoding='utf-8') as f:f.write(json.dumps({'time_utc':datetime.now(timezone.utc).isoformat(),'source_rows':len(src),'cids':len(cids),'queried_new_cids':len(todo),'cache_records':len(cached),'failures':failures})+'\n')
