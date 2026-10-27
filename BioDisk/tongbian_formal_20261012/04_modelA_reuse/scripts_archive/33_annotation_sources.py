from direct_http import *
from concurrent.futures import ThreadPoolExecutor,as_completed
sources=[
('uniprot_human_go','https://rest.uniprot.org/uniprotkb/stream',{'query':'(organism_id:9606) AND (reviewed:true)','format':'tsv','fields':'accession,gene_primary,organism_id,xref_geneid,go_id'}),
('go_basic_ontology','https://purl.obolibrary.org/obo/go/go-basic.obo',None),
('kegg_human_gene_pathways','https://rest.kegg.jp/link/pathway/hsa',None),
('kegg_database_version','https://rest.kegg.jp/info/pathway',None),
('string_version','https://string-db.org/api/json/version',None),
('cmro_vol3_main_retry','https://www.cmro.gov.hk/files/hkcmms/vol3/main.html',None)]
with ThreadPoolExecutor(max_workers=3) as pool:
 for f in as_completed([pool.submit(get,u,n,params=p,timeout=(20,60)) for n,u,p in sources]):
  r,resp=f.result();print(json.dumps({k:r.get(k) for k in ['label','status','bytes','error']},ensure_ascii=False),flush=True)
