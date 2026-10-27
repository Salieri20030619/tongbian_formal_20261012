"""Collect every source molecule detail in the predeclared complete 16-herb cohort."""
from direct_http import *
from bs4 import BeautifulSoup
from concurrent.futures import ThreadPoolExecutor,as_completed
from urllib.parse import urljoin
import pandas as pd,re
mol=pd.read_csv(ROOT/'tables'/'05_tcmsp_ingredient_rows.tsv',sep='\t',dtype=str).drop_duplicates('MOL_ID').sort_values('MOL_ID')
LIMITER=threading.Lock();NEXT=[0.0];STOP=threading.Event()
def parse(data,mid,mol_id):
 soup=BeautifulSoup(data,'html.parser');out={'MOL_ID':mol_id,'molecule_ID':mid,'source_url':f'https://old.tcmsp-e.com/molecule.php?qn={mid}','source_raw':f'raw/mol_{mol_id}.bin'}
 for tr in soup.find_all('tr'):
  cells=tr.find_all(['td','th'],recursive=False)
  if len(cells)<2:continue
  key=cells[0].get_text(' ',strip=True);val=cells[1].get_text(' ',strip=True)
  if key in ['Molecule ID','Molecule name','InChIKey','Pubchem Cid','CAS']:out[key]=val
 cids=[]
 for a in soup.find_all('a',href=True):
  url=a['href']
  if 'pubchem' in url.lower():
   cids+=re.findall(r'(?:cid=|/compound/)(\d+)',url)
  if re.search(r'\.mol2?(?:\?|$)|\.sdf(?:\?|$)',url):out['structure_download_url']=urljoin(out['source_url'],url)
 out['pubchem_cids']=';'.join(sorted(set(cids),key=int))
 out['parse_status']='ok' if out.get('Molecule ID')==mol_id else 'missing_or_mismatch'
 return out
def one(row):
 mid=row['molecule_ID'];mol_id=row['MOL_ID'];path=ROOT/'raw'/f'mol_{mol_id}.bin'
 if path.exists():return parse(path.read_bytes(),mid,mol_id)
 if mol_id=='MOL000018':
  data=(ROOT/'raw'/'tcmsp_molecule_sample.bin').read_bytes();path.write_bytes(data);return parse(data,mid,mol_id)
 if STOP.is_set():return {'MOL_ID':mol_id,'molecule_ID':mid,'parse_status':'not_requested_after_block'}
 with LIMITER:
  delay=max(0,NEXT[0]-time.monotonic())
  if delay:time.sleep(delay)
  NEXT[0]=time.monotonic()+.6
 r,resp=get(f'https://old.tcmsp-e.com/molecule.php?qn={mid}',f'mol_{mol_id}',timeout=(12,30))
 if resp is not None and (resp.status_code in [401,403,429] or 'captcha' in resp.text.lower()):STOP.set()
 if resp is None or resp.status_code!=200:return {'MOL_ID':mol_id,'molecule_ID':mid,'parse_status':'request_failed','error':r.get('error'),'http_status':r.get('status')}
 return parse(resp.content,mid,mol_id)
records=[]
with ThreadPoolExecutor(max_workers=2) as pool:
 for fut in as_completed([pool.submit(one,row) for row in mol.to_dict('records')]):
  records.append(fut.result())
  if len(records)%50==0 or len(records)==len(mol):
   pd.DataFrame(records).sort_values('MOL_ID').to_csv(ROOT/'tables'/'11_tcmsp_structure_identifiers.tsv',sep='\t',index=False,encoding='utf-8-sig')
   print(json.dumps({'completed':len(records),'total':len(mol),'ok':sum(r.get('parse_status')=='ok' for r in records),'blocked':STOP.is_set()}),flush=True)
