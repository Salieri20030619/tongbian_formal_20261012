from direct_http import *
from concurrent.futures import ThreadPoolExecutor,as_completed
SOURCES=[('tcmsp_molecule_sample','https://old.tcmsp-e.com/molecule.php?qn=18'),('tcmsp_target_sample','https://old.tcmsp-e.com/target.php?qt=38'),('sea_search','https://sea.bkslab.org/search'),('sea_captcha_status','https://sea.bkslab.org/api/captcha/challenge')]
with ThreadPoolExecutor(max_workers=3) as pool:
 for fut in as_completed([pool.submit(get,u,n) for n,u in SOURCES]):
  r,resp=fut.result();print(json.dumps({k:r.get(k) for k in ['label','status','bytes','error']},ensure_ascii=False),flush=True)
