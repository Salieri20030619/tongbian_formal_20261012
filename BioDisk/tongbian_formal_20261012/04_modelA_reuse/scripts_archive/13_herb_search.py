from direct_http import *
from concurrent.futures import ThreadPoolExecutor,as_completed
queries=[('gualouren','瓜蒌仁'),('gualouzi','瓜蒌子'),('maidong','麦冬'),('dihuang','地黄'),('tiepi','铁皮石斛')]
def query(label,term):
 r,resp=get('http://herb.ac.cn/chedi/api/',f'herb_search_{label}',method='POST',json_body={'func_name':'search_api','label':'Herb','keyword':term})
 print(json.dumps({'label':label,'status':r.get('status'),'bytes':r.get('bytes'),'body_preview':resp.text[:2500] if resp is not None else r.get('error')},ensure_ascii=False),flush=True)
with ThreadPoolExecutor(max_workers=2) as pool:
 for f in as_completed([pool.submit(query,n,q) for n,q in queries]):f.result()
r,resp=get('http://47.92.70.12','herb_official_alt_home')
print(json.dumps({'label':'official_alternate','status':r.get('status'),'bytes':r.get('bytes'),'body_preview':resp.text[:1100] if resp is not None else r.get('error')},ensure_ascii=False),flush=True)
