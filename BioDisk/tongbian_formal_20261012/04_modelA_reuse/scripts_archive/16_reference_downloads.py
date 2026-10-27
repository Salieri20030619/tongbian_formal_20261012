from direct_http import *
from concurrent.futures import ThreadPoolExecutor,as_completed
SOURCES=[('uniprot_reviewed_human_proteome','https://rest.uniprot.org/uniprotkb/stream',{'query':'(organism_id:9606) AND (reviewed:true)','format':'tsv','fields':'accession,id,protein_name,gene_names,organism_id,reviewed,xref_geneid'}),('hpo_annotation_docs','https://raw.githubusercontent.com/obophenotype/human-phenotype-ontology/master/README.md',None),('hpo_download_docs','https://hpo.jax.org/data/annotations',None),('hpo_genes_to_phenotype','https://purl.obolibrary.org/obo/hp/hpoa/genes_to_phenotype.txt',None),('hpo_phenotype_hpoa','https://purl.obolibrary.org/obo/hp/hpoa/phenotype.hpoa',None)]
with ThreadPoolExecutor(max_workers=3) as pool:
 for f in as_completed([pool.submit(get,u,n,params=p,timeout=(12,60)) for n,u,p in SOURCES]):
  r,resp=f.result();print(json.dumps({k:r.get(k) for k in ['label','status','bytes','error']},ensure_ascii=False),flush=True)
