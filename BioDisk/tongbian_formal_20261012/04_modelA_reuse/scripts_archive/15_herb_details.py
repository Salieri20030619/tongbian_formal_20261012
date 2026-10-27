from direct_http import *
import pandas as pd
MAPPINGS={'H01':['HERB000309'],'H02':['HERB002560'],'H03':['HERB004692'],'H04':['HERB001210','HERB000742'],'H05':['HERB000775'],'H06':['HERB001952'],'H07':['HERB003174'],'H08':['HERB002640'],'H09':['HERB007034'],'H10':['HERB000497'],'H11':['HERB002342','HERB002339'],'H12':['HERB002876'],'H13':['HERB003708'],'H14':['HERB006280'],'H15':['HERB005012','HERB001783'],'H16':['HERB004938'],'H17':['HERB005502'],'H18':['HERB002534'],'H19':['HERB001211'],'H20':['HERB001915','HERB004184']}
def text(v):return v.get('title','') if isinstance(v,dict) else v
out=[];ing=[]
for hid,ids in MAPPINGS.items():
 for dbid in ids:
  label=f'herb_detail_{dbid}';p=ROOT/'raw'/f'{label}.bin'
  if p.exists():d=json.loads(p.read_text(encoding='utf-8'))
  else:
   r,resp=get('http://herb.ac.cn/chedi/api/',label,method='POST',json_body={'func_name':'detail_api','label':'Herb','v':dbid,'key_id':dbid},timeout=(12,40))
   if resp is None or resp.status_code!=200:continue
   d=resp.json()
  summary=d.get('summary');s={}
  if summary and len(summary)>=2:s=dict(zip([text(x) for x in summary[0]],[text(x) for x in summary[1]]))
  ingred=d.get('herb_ingredient',[])
  out.append({'herb_id':hid,'herb_db_id':dbid,'summary':json.dumps(s,ensure_ascii=False),'ingredient_rows':max(0,len(ingred)-1),'source_raw':f'raw/{label}.bin'})
  if ingred:
   for row in ingred[1:]:
    ing.append(dict(zip([text(x) for x in ingred[0]],[text(x) for x in row]),herb_id=hid,herb_db_id=dbid,source_raw=f'raw/{label}.bin',source_url=f'http://herb.ac.cn/Detail/?v={dbid}&label=Herb'))
  print(hid,dbid,'summary',json.dumps(s,ensure_ascii=False),'ingredients',len(ingred)-1,'keys',list(d.keys()),flush=True)
  time.sleep(.35)
pd.DataFrame(out).to_csv(ROOT/'tables'/'13_herb_details_scope_review.tsv',sep='\t',index=False,encoding='utf-8-sig')
pd.DataFrame(ing).to_csv(ROOT/'tables'/'14_herb_ingredients_CANDIDATE.tsv',sep='\t',index=False,encoding='utf-8-sig')
