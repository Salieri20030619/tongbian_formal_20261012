from direct_http import *
from bs4 import BeautifulSoup
queries=[('maimendong','麦门冬'),('ophiopogon','Ophiopogon japonicus'),('ophiopogon_radix','Radix Ophiopogonis')]
for label,q in queries:
 r,resp=get('http://herb.ac.cn/chedi/api/','herb_alias_'+label,method='POST',json_body={'func_name':'search_api','label':'Herb','keyword':q})
 print(label,resp.text[:4500] if resp is not None else r.get('error'),flush=True);time.sleep(.4)
checks=[('cmro_hkcmms','https://www.cmro.gov.hk/html/eng/useful_information/hkcmms/volumes.html'),('bing_rehmannia_official','https://www.bing.com/search?q=site%3Acmro.gov.hk+Rehmanniae+Radix+dried'),('bing_shengdihuang_official','https://www.bing.com/search?q=site%3Agov.cn+%E5%B9%B2%E5%9C%B0%E9%BB%84+%E7%94%9F%E5%9C%B0%E9%BB%84')]
for label,url in checks:
 r,resp=get(url,label)
 print(label,r.get('status'),BeautifulSoup(resp.content,'html.parser').get_text(' ',strip=True)[:13000] if resp is not None else r.get('error'),flush=True)
