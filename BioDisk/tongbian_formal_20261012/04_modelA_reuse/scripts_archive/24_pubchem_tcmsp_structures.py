from direct_http import *
import pandas as pd,re
from urllib.parse import quote
cache=ROOT/'raw'/'pubchem_properties';cache.mkdir(exist_ok=True)
src=pd.read_csv(ROOT/'tables'/'11_tcmsp_structure_identifiers.tsv',sep='\t',dtype=str).fillna('')
properties='InChIKey,IsomericSMILES,CanonicalSMILES,MolecularFormula,InChI'
collected=[];failures=[]
def fetch(kind,values):
 label='pubchem_'+kind+'_'+hashlib.sha256(','.join(values).encode()).hexdigest()[:16]
 url='https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/'+kind+'/'+quote(','.join(values),safe=',-')+'/property/'+properties+'/JSON'
 p=ROOT/'raw'/f'{label}.bin'
 if p.exists():
  try:data=json.loads(p.read_text(encoding='utf-8'));status=200 if 'PropertyTable' in data else 400
  except ValueError:data={};status=400
 else:
  r,resp=get(url,label,timeout=(12,40));status=r.get('status');data={}
  if resp is not None:
   try:data=resp.json()
   except ValueError:pass
  if status in [401,403,429]:raise SystemExit('PubChem access/rate block; stop without retry')
  time.sleep(.3)
 if 'PropertyTable' in data:
  for rec in data['PropertyTable']['Properties']:
   rec['retrieval_raw']=f'raw/{label}.bin';rec['retrieval_url']=url
   (cache/f"{rec['CID']}.json").write_text(json.dumps(rec,ensure_ascii=False),encoding='utf-8');collected.append(rec)
  return
 if len(values)>1 and status in [400,404]:
  half=len(values)//2;fetch(kind,values[:half]);fetch(kind,values[half:])
 else:failures.append({'kind':kind,'query':','.join(values),'status':status,'raw':f'raw/{label}.bin','fault':json.dumps(data)[:1000]})
cids=sorted({x for value in src.pubchem_cids for x in value.split(';') if x.isdigit()},key=int)
todo=[c for c in cids if not (cache/f'{c}.json').exists()]
for i in range(0,len(todo),50):
 fetch('cid',todo[i:i+50]);print('CID batch',min(i+50,len(todo)),len(todo),flush=True)
cached=[json.loads(p.read_text(encoding='utf-8')) for p in cache.glob('*.json')]
known_keys={x['InChIKey'] for x in cached}
keys=sorted({k for k in src.InChIKey if re.fullmatch(r'[A-Z]{14}-[A-Z]{10}-[A-Z]',k) and k not in known_keys})
for i in range(0,len(keys),20):
 fetch('inchikey',keys[i:i+20]);print('InChIKey batch',min(i+20,len(keys)),len(keys),flush=True)
cached=[json.loads(p.read_text(encoding='utf-8')) for p in cache.glob('*.json')]
pd.DataFrame(cached).to_csv(ROOT/'tables'/'22_pubchem_structure_properties.tsv',sep='\t',index=False,encoding='utf-8-sig')
pd.DataFrame(failures,columns=['kind','query','status','raw','fault']).to_csv(ROOT/'tables'/'23_pubchem_lookup_failures.tsv',sep='\t',index=False,encoding='utf-8-sig')
print(json.dumps({'source_cids':len(cids),'property_records_cached':len(cached),'source_keys_needed_lookup':len(keys),'failed_queries':len(failures)}),flush=True)
