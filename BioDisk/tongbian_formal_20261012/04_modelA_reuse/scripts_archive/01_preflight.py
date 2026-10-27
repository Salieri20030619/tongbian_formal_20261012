from direct_http import *
from concurrent.futures import ThreadPoolExecutor,as_completed
import csv,sys,platform,importlib.metadata
SOURCES=[
 ('tcmsp_old','https://old.tcmsp-e.com/'),
 ('tcmsp','https://tcmsp-e.com/tcmsp.php'),
 ('herb','http://herb.ac.cn/'),
 ('pubchem','https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/2244/property/InChIKey,CanonicalSMILES/JSON'),
 ('sea','https://sea.bkslab.org/'),
 ('swiss_target_prediction','http://www.swisstargetprediction.ch/'),
 ('genecards','https://www.genecards.org/'),
 ('string','https://string-db.org/api/tsv/get_string_ids?identifiers=TP53&species=9606&limit=1'),
 ('ctd','https://ctdbase.org/downloads/'),
 ('disgenet','https://disgenet.com/'),
 ('diseases','https://diseases.jensenlab.org/Downloads'),
 ('ncbi_eutils','https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=medgen&term=constipation&retmode=json'),
 ('chembl','https://www.ebi.ac.uk/chembl/api/data/status.json'),
 ('uniprot','https://rest.uniprot.org/uniprotkb/P04637.json'),
 ('kegg','https://rest.kegg.jp/list/pathway/hsa'),
 ('gprofiler','https://biit.cs.ut.ee/gprofiler/api/util/organisms_list'),
 ('imppat','https://cb.imsc.res.in/imppat/'),
 ('symmap','https://www.symmap.org/'),
 ('hpo','https://purl.obolibrary.org/obo/hp/hpoa/phenotype_to_genes.txt'),
 ('opentargets','https://api.platform.opentargets.org/api/v4/graphql')
]
rows=[]
with ThreadPoolExecutor(max_workers=5) as pool:
    for future in as_completed([pool.submit(get,url,name) for name,url in SOURCES]):
        rec,resp=future.result();rows.append(rec)
        print(json.dumps({k:rec.get(k) for k in ['label','status','bytes','elapsed_seconds','error']},ensure_ascii=False),flush=True)
(ROOT/'results'/('preflight_'+('host_direct' if '--host-direct' in sys.argv else 'sandbox')+'.json')).write_text(json.dumps(rows,ensure_ascii=False,indent=2),encoding='utf-8')
fields=['label','url','status','bytes','elapsed_seconds','error','sha256','raw_file']
with (ROOT/'tables'/'01_source_preflight.tsv').open('w',encoding='utf-8-sig',newline='') as f:
    w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',extrasaction='ignore');w.writeheader();w.writerows(rows)
versions={'python':sys.version,'platform':platform.platform(),'packages':{p:importlib.metadata.version(p) for p in ['requests','pandas','numpy','scipy','statsmodels','networkx','matplotlib','reportlab','openpyxl','beautifulsoup4']}}
(ROOT/'logs'/'environment.json').write_text(json.dumps(versions,indent=2),encoding='utf-8')
