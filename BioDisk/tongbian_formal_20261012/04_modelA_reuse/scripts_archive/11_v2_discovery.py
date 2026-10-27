from direct_http import *
from concurrent.futures import ThreadPoolExecutor,as_completed
SOURCES=[('herb_v2_js','http://herb.ac.cn/static/umi.js?version=762'),('herb_v2_download_page','http://herb.ac.cn/v2/Download/'),('tcmsp_query_alias_sample','https://old.tcmsp-e.com/tcmspsearch.php?qs=herb_all_name&q=Rehmannia&token=2e7649d5eac783b1b3b1a5f47db819ce'),('uniprot_fields_sample','https://rest.uniprot.org/uniprotkb/search?query=%28protein_name%3A%22Muscarinic+acetylcholine+receptor+M1%22%29+AND+%28organism_id%3A9606%29&format=tsv&fields=accession,id,protein_name,gene_names,organism_id,reviewed,xref_geneid&size=100')]
with ThreadPoolExecutor(max_workers=3) as pool:
 for fut in as_completed([pool.submit(get,u,n) for n,u in SOURCES]):
  r,resp=fut.result();print(json.dumps({k:r.get(k) for k in ['label','status','bytes','error']},ensure_ascii=False),flush=True)
