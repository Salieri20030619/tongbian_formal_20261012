from direct_http import *
from concurrent.futures import ThreadPoolExecutor,as_completed
from bs4 import BeautifulSoup
SOURCES=[
('tcmsp_old_search','https://old.tcmsp-e.com/tcmsp.php'),
('swiss_terms','https://www.swisstargetprediction.ch/termsofuse.php'),
('diseases_knowledge','https://download.jensenlab.org/human_disease_knowledge_full.tsv'),
('diseases_experiments','https://download.jensenlab.org/human_disease_experiments_full.tsv'),
('hpo_ontology','https://purl.obolibrary.org/obo/hp.obo'),
('herb_robots','http://herb.ac.cn/robots.txt'),
('sea_robots','https://sea.bkslab.org/robots.txt'),
('tcmsp_robots','https://old.tcmsp-e.com/robots.txt')]
with ThreadPoolExecutor(max_workers=4) as pool:
 for fut in as_completed([pool.submit(get,url,name,timeout=(12,60)) for name,url in SOURCES]):
    r,resp=fut.result();print(json.dumps({k:r.get(k) for k in ['label','status','bytes','error']},ensure_ascii=False),flush=True)
for name in ['tcmsp_old_search','herb','sea','swiss_terms']:
 p=ROOT/'raw'/f'{name}.bin'
 if not p.exists():continue
 soup=BeautifulSoup(p.read_bytes(),'html.parser')
 print('\n'+name+' TEXT: '+soup.get_text(' ',strip=True)[:9000])
 print('SCRIPTS: '+str([x.get('src') for x in soup.find_all('script')]))
 print('FORMS: '+str([str(x)[:3000] for x in soup.find_all('form')]))
