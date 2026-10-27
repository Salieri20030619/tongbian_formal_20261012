from direct_http import *
from concurrent.futures import ThreadPoolExecutor,as_completed
import pandas as pd,ast,re
edges_path=ROOT/'tables'/'18_herb_ingredient_target_edges_UNCLASSIFIED.tsv'
edges=pd.read_csv(edges_path,sep='\t',dtype=str).fillna('')
allowed=set(pd.read_csv(ROOT/'tables'/'06_tcmsp_target_rows_UNCLASSIFIED.tsv',sep='\t',dtype=str).target_ID)
cross=[]
for row in edges.to_dict('records'):
 if not row.get('Database sources'):continue
 try:sources=ast.literal_eval(row['Database sources'])
 except (ValueError,SyntaxError):continue
 for source in sources:
  m=re.fullmatch(r'TCMSP:TAR(\d+)',source.get('title',''))
  if m and str(int(m.group(1))) in allowed:
   cross.append({'tcmsp_target_id':str(int(m.group(1))),'herb_target_id':row['Target id'],'herb_target_symbol':row['Target name'],'crosswalk_evidence_raw':row['source_raw'],'observed_ingredient_id':row['herb_ingredient_id']})
df=pd.DataFrame(cross);df.to_csv(ROOT/'tables'/'20_tcmsp_herb_target_crosswalk_evidence.tsv',sep='\t',index=False,encoding='utf-8-sig')
ids=sorted(df.herb_target_id.unique());RL=threading.Lock();NEXT=[0.0];STOP=threading.Event()
def flatten(v):return v.get('title','') if isinstance(v,dict) else v
def parse(dbid,data):
 d=json.loads(data);s=d.get('summary',[]);out={'herb_target_id':dbid,'source_raw':f'raw/htarget_{dbid}.bin','source_url':f'http://herb.ac.cn/Detail/?v={dbid}&label=Target'}
 if len(s)>=2:out.update(dict(zip([flatten(x) for x in s[0]],[flatten(x) for x in s[1]])))
 out['parse_status']='ok' if out.get('Target id')==dbid else 'missing_or_mismatch';return out
def one(dbid):
 p=ROOT/'raw'/f'htarget_{dbid}.bin'
 if p.exists():return parse(dbid,p.read_text(encoding='utf-8'))
 if dbid=='HBTAR003295':p.write_bytes((ROOT/'raw'/'herb_target_detail_sample.bin').read_bytes());return parse(dbid,p.read_text(encoding='utf-8'))
 if STOP.is_set():return {'herb_target_id':dbid,'parse_status':'not_requested_after_block'}
 with RL:
  delay=max(0,NEXT[0]-time.monotonic())
  if delay:time.sleep(delay)
  NEXT[0]=time.monotonic()+1.0
 r,resp=get('http://herb.ac.cn/chedi/api/',f'htarget_{dbid}',method='POST',json_body={'func_name':'detail_api','label':'Target','v':dbid,'key_id':dbid},timeout=(12,45))
 if resp is not None and (resp.status_code in [401,403,429] or 'captcha' in resp.text.lower()):STOP.set()
 if resp is None or resp.status_code!=200:return {'herb_target_id':dbid,'parse_status':'request_failed','error':r.get('error')}
 try:return parse(dbid,resp.text)
 except Exception as ex:return {'herb_target_id':dbid,'parse_status':'parse_failed','error':str(ex)}
results=[]
with ThreadPoolExecutor(max_workers=2) as pool:
 for f in as_completed([pool.submit(one,x) for x in ids]):
  results.append(f.result())
  if len(results)%50==0 or len(results)==len(ids):
   pd.DataFrame(results).sort_values('herb_target_id').to_csv(ROOT/'tables'/'21_herb_target_taxonomy_verified.tsv',sep='\t',index=False,encoding='utf-8-sig')
   print(json.dumps({'target_details':len(results),'total_snapshot_targets':len(ids),'source_taxon9606':sum(x.get('Taxonomy id')=='9606' for x in results),'blocked':STOP.is_set()}),flush=True)
