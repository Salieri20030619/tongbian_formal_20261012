from direct_http import *
from concurrent.futures import ThreadPoolExecutor,as_completed
import pandas as pd
rows=pd.read_csv(ROOT/'tables'/'14_herb_ingredients_CANDIDATE.tsv',sep='\t',dtype=str)
scope=pd.read_csv(ROOT/'tables'/'13_herb_details_scope_review.tsv',sep='\t',dtype=str)
decisions=[]
for r in scope.to_dict('records'):
 db=r['herb_db_id'];status='medicinal_drug_record';note='数据库药材成分关联；实际供应基原和炮制未核验'
 if db=='HERB001952':status='exclude_part_mixed';note='整瓜蒌记录混合果实/皮/种子，不替代瓜蒌仁'
 elif db=='HERB000742':status='exclude_species_unknown';note='同名当归缺基原和部位，不能排除其他Angelica种'
 elif db=='HERB003708':status='conditional_metadata_conflict';note='中文麦冬+拉丁Ophiopogon japonicus+块根匹配；英文写Liriope，保留冲突'
 elif db=='HERB005012':status='conditional_processing_conflict';note='中文生地黄，英文Fresh root，部位字段填物种；未把它等同干燥生地'
 elif db=='HERB001783':status='conditional_raw_dried_candidate';note='干地黄+Rehmannia glutinosa+dried root；为生地黄候选同义材料，保留与原药名映射待核验'
 elif db=='HERB000309':note+='；Use part写root，但拉丁/英文明确rhizome，保留字段冲突'
 decisions.append(dict(r,material_status=status,material_note=note))
dec=pd.DataFrame(decisions);dec.to_csv(ROOT/'tables'/'15_herb_material_scope_decisions.tsv',sep='\t',index=False,encoding='utf-8-sig')
rows=rows.merge(dec[['herb_db_id','material_status','material_note']],on='herb_db_id',validate='many_to_one')
eligible=rows[~rows.material_status.str.startswith('exclude_')].copy()
eligible.to_csv(ROOT/'tables'/'16_herb_memberships_eligible.tsv',sep='\t',index=False,encoding='utf-8-sig')
ids=sorted(eligible['Ingredient id'].unique());LOCK_RATE=threading.Lock();NEXT=[0.0];STOP=threading.Event()
def flatten(v):return v.get('title','') if isinstance(v,dict) else v
def parse(dbid,data):
 d=json.loads(data);summary=d.get('summary',[]);out={'herb_ingredient_id':dbid,'source_raw':f'raw/hing_{dbid}.bin','source_url':f'http://herb.ac.cn/Detail/?v={dbid}&label=Ingredient'}
 if len(summary)>=2:out.update(dict(zip([flatten(x) for x in summary[0]],[flatten(x) for x in summary[1]])))
 out['parse_status']='ok' if out.get('Ingredient id')==dbid else 'missing_or_mismatch'
 edges=[]
 for key,cat in [('ingredient_target','HERB_database_association_UNCLASSIFIED'),('drug_paper_target','HERB_reference_mining_NOT_binding_validation')]:
  arr=d.get(key,[])
  for row in arr[1:]:
   edge=dict(zip([flatten(x) for x in arr[0]],[flatten(x) for x in row]));edge.update(herb_ingredient_id=dbid,evidence_class=cat,source_raw=out['source_raw'],source_url=out['source_url']);edges.append(edge)
 return out,edges
def one(dbid):
 p=ROOT/'raw'/f'hing_{dbid}.bin';sample=ROOT/'raw'/f'herb_ing_sample_{dbid}.bin'
 if p.exists():return parse(dbid,p.read_text(encoding='utf-8'))
 if sample.exists():p.write_bytes(sample.read_bytes());return parse(dbid,p.read_text(encoding='utf-8'))
 if STOP.is_set():return {'herb_ingredient_id':dbid,'parse_status':'not_requested_after_block'},[]
 with LOCK_RATE:
  delay=max(0,NEXT[0]-time.monotonic())
  if delay:time.sleep(delay)
  NEXT[0]=time.monotonic()+.4
 r,resp=get('http://herb.ac.cn/chedi/api/',f'hing_{dbid}',method='POST',json_body={'func_name':'detail_api','label':'Ingredient','v':dbid,'key_id':dbid},timeout=(12,35))
 if resp is not None and (resp.status_code in [401,403,429] or 'captcha' in resp.text.lower()):STOP.set()
 if resp is None or resp.status_code!=200:return {'herb_ingredient_id':dbid,'parse_status':'request_failed','http_status':r.get('status'),'error':r.get('error')},[]
 try:return parse(dbid,resp.text)
 except Exception as ex:return {'herb_ingredient_id':dbid,'parse_status':'parse_failed','error':str(ex)},[]
records=[];edges=[]
with ThreadPoolExecutor(max_workers=2) as pool:
 for f in as_completed([pool.submit(one,x) for x in ids]):
  rec,ed=f.result();records.append(rec);edges.extend(ed)
  if len(records)%100==0 or len(records)==len(ids):
   pd.DataFrame(records).sort_values('herb_ingredient_id').to_csv(ROOT/'tables'/'17_herb_ingredient_details.tsv',sep='\t',index=False,encoding='utf-8-sig')
   pd.DataFrame(edges).to_csv(ROOT/'tables'/'18_herb_ingredient_target_edges_UNCLASSIFIED.tsv',sep='\t',index=False,encoding='utf-8-sig')
   print(json.dumps({'completed':len(records),'total':len(ids),'parsed_ok':sum(x.get('parse_status')=='ok' for x in records),'edges':len(edges),'blocked':STOP.is_set()}),flush=True)
