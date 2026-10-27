from direct_http import *
import pandas as pd
rows=pd.read_csv(ROOT/'tables'/'14_herb_ingredients_CANDIDATE.tsv',sep='\t',dtype=str)
ids=sorted(set(rows[rows.herb_id.isin(['H13','H15','H17'])]['Ingredient id']))[:2]+['HBIN020499']
for dbid in ids:
 r,resp=get('http://herb.ac.cn/chedi/api/',f'herb_ing_sample_{dbid}',method='POST',json_body={'func_name':'detail_api','label':'Ingredient','v':dbid,'key_id':dbid},timeout=(12,40))
 if resp is None:print(r);continue
 d=resp.json();print(dbid,json.dumps({k:(v[:3] if isinstance(v,list) else str(v)[:400]) for k,v in d.items() if k!='Merge_diff_gene'},ensure_ascii=False),flush=True)
