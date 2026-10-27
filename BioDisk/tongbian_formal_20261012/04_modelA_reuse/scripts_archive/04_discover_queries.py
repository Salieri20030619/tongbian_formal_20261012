from direct_http import *
from concurrent.futures import ThreadPoolExecutor,as_completed
from bs4 import BeautifulSoup
SOURCES=[('herb_js','http://herb.ac.cn/static/umi.js?version=474'),('sea_js','https://sea.bkslab.org/assets/index-D6mhr3nD.js'),('herb_v2','http://herb.ac.cn/v2'),('tcmsp_browse','https://old.tcmsp-e.com/browse.php?qc=herbs')]
with ThreadPoolExecutor(max_workers=3) as pool:
 for fut in as_completed([pool.submit(get,u,n) for n,u in SOURCES]):
    r,resp=fut.result();print(json.dumps({k:r.get(k) for k in ['label','status','bytes','error']},ensure_ascii=False),flush=True)
session=new_session()
r,resp=get('https://old.tcmsp-e.com/tcmsp.php','tcmsp_form_session',session=session)
soup=BeautifulSoup(resp.content,'html.parser')
token=soup.select_one('input[name=token]')['value']
r,resp=get('https://old.tcmsp-e.com/tcmspsearch.php','tcmsp_sample_baizhu',params={'qs':'herb_all_name','q':'白术','token':token},session=session)
print(json.dumps({k:r.get(k) for k in ['label','status','bytes','error']},ensure_ascii=False),flush=True)
print(BeautifulSoup(resp.content,'html.parser').get_text(' ',strip=True)[-15000:])
