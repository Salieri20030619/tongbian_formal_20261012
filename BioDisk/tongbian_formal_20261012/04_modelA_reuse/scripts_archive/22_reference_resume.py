from direct_http import *
from concurrent.futures import ThreadPoolExecutor,as_completed
SOURCES=[('hpo_genes_to_phenotype','https://purl.obolibrary.org/obo/hp/hpoa/genes_to_phenotype.txt'),('hpo_phenotype_hpoa','https://purl.obolibrary.org/obo/hp/hpoa/phenotype.hpoa'),('hgnc_complete_set','https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt')]
with ThreadPoolExecutor(max_workers=2) as pool:
 for f in as_completed([pool.submit(get,u,n,timeout=(12,60)) for n,u in SOURCES if not (ROOT/'raw'/f'{n}.bin').exists()]):
  r,resp=f.result();print(json.dumps({k:r.get(k) for k in ['label','status','bytes','error']},ensure_ascii=False),flush=True)
