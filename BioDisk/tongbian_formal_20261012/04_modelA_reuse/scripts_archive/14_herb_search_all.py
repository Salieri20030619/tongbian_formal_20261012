from direct_http import *
import pandas as pd
def val(x):return x.get('title','') if isinstance(x,dict) else x
formula=pd.read_csv(ROOT/'tables'/'02_formula_standardization.tsv',sep='\t',dtype=str)
queries=[(r.herb_id,('生地黄' if r.herb_id=='H15' else r.search_drug_name)) for r in formula.itertuples()]
queries += [('H06','栝楼子'),('H06','栝楼仁'),('H06','Trichosanthis Semen'),('H06','Trichosanthes')]
out=[]
for ix,(hid,q) in enumerate(queries):
 label=f'herb_all_search_{hid}_{ix:02d}';p=ROOT/'raw'/f'{label}.bin'
 if p.exists():data=json.loads(p.read_text(encoding='utf-8'))
 else:
  rec,resp=get('http://herb.ac.cn/chedi/api/',label,method='POST',json_body={'func_name':'search_api','label':'Herb','keyword':q})
  if resp is None or resp.status_code!=200:continue
  data=resp.json()
 records=data.get('res_data',[])
 for row in records[1:]:
  out.append(dict(herb_id=hid,query=q,herb_db_id=val(row[0]),pinyin=val(row[1]),chinese=val(row[2]),english=val(row[3]),latin=val(row[4]),source_raw=f'raw/{label}.bin',source_url='http://herb.ac.cn'+row[0].get('link','')))
 print(hid,q,[(val(r[0]),val(r[2]),val(r[4])) for r in records[1:]],flush=True)
 time.sleep(.3)
pd.DataFrame(out).to_csv(ROOT/'tables'/'12_herb_all_search_candidates.tsv',sep='\t',index=False,encoding='utf-8-sig')
