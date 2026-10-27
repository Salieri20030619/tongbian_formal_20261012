from direct_http import *
from concurrent.futures import ThreadPoolExecutor,as_completed
sources=[('official_ophiopogon_standard','https://www.cmro.gov.hk/files/hkcmms/vol3/pdf_e/Radix_Ophiopogonis_v3_e.pdf'),('official_rehmannia_standard','https://www.cmro.gov.hk/files/hkcmms/vol3/pdf_e/Radix_Rehmanniae_v3_e.pdf')]
with ThreadPoolExecutor(max_workers=2) as pool:
 for f in as_completed([pool.submit(get,u,n,timeout=(20,40)) for n,u in sources]):
  r,resp=f.result();print(json.dumps({k:r.get(k) for k in ['label','status','bytes','error']},ensure_ascii=False),flush=True)
