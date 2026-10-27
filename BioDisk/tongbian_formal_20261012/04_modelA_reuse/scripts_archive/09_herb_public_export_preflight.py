from direct_http import *
from concurrent.futures import ThreadPoolExecutor,as_completed
# These paths and file_path parameter are copied from the public Download menu
# and its onClick handler in raw/herb_js.bin; no invented path or hidden endpoint.
SOURCES=[('herb_herb_info','HERB_herb_info.txt'),('herb_ingredient_info','HERB_ingredient_info.txt'),('herb_target_info','HERB_target_info.txt')]
def fetch_export(label,filename):
 return get('http://herb.ac.cn/download/file/',label,params={'file_path':'/data/Web_server/HERB_web/static/download_data/'+filename},timeout=(12,45))
with ThreadPoolExecutor(max_workers=2) as pool:
 futures=[pool.submit(fetch_export,n,f) for n,f in SOURCES]
 futures.append(pool.submit(get,'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/6321405/property/InChIKey,CanonicalSMILES/JSON','pubchem_formula_sample'))
 for fut in as_completed(futures):
  r,resp=fut.result();print(json.dumps({k:r.get(k) for k in ['label','status','bytes','content_type','error']},ensure_ascii=False),flush=True)
