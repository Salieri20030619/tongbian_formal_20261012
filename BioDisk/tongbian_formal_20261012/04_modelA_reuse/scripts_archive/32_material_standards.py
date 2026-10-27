from direct_http import *
from bs4 import BeautifulSoup
r,resp=get('https://www.cmro.gov.hk/html/eng/useful_information/hkcmms/cmmlist.html','cmro_cmm_list')
if resp is not None:
 soup=BeautifulSoup(resp.content,'html.parser')
 for tr in soup.find_all('tr'):
  text=tr.get_text(' ',strip=True)
  if any(w.lower() in text.lower() for w in ['Rehmann','Ophiopog','Trichosan','地黄','麥冬','麦冬']):print(str(tr),flush=True)
for dbid in ['HERB001010','HERB003367']:
 label='herb_detail_'+dbid;p=ROOT/'raw'/f'{label}.bin'
 if p.exists():d=json.loads(p.read_text(encoding='utf-8'))
 else:
  r,resp=get('http://herb.ac.cn/chedi/api/',label,method='POST',json_body={'func_name':'detail_api','label':'Herb','v':dbid,'key_id':dbid},timeout=(12,35))
  if resp is None:continue
  d=resp.json()
 print(dbid,json.dumps(d.get('summary'),ensure_ascii=False),'ingredients',len(d.get('herb_ingredient',[]))-1,flush=True)
