from direct_http import *
from concurrent.futures import ThreadPoolExecutor,as_completed
import pandas as pd
src=pd.read_csv(ROOT/'tables'/'26_tcmsp_structure_qc.tsv',sep='\t',dtype=str).fillna('')
rows=src[(src.structure_status=='unresolved')&src.structure_download_url.str.startswith('https://old.tcmsp-e.com/tcmspmol/')].to_dict('records')
RATE=threading.Lock();NEXT=[0.0];STOP=threading.Event()
def one(s):
 label='mol2_'+s['MOL_ID']
 if (ROOT/'raw'/f'{label}.bin').exists():return 'cached'
 if STOP.is_set():return 'not_requested_after_block'
 with RATE:
  delay=max(0,NEXT[0]-time.monotonic())
  if delay:time.sleep(delay)
  NEXT[0]=time.monotonic()+.6
 r,resp=get(s['structure_download_url'],label)
 if r.get('status') in [401,403,429]:STOP.set()
 return r.get('status')
counts={}
with ThreadPoolExecutor(max_workers=2) as pool:
 for f in as_completed([pool.submit(one,r) for r in rows]):
  k=str(f.result());counts[k]=counts.get(k,0)+1
print(json.dumps({'requested_cohort':len(rows),'results':counts}),flush=True)
