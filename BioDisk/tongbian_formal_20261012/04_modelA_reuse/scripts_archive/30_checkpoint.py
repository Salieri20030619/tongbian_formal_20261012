from pathlib import Path
import json,pandas as pd
from datetime import datetime,timezone
ROOT=Path(__file__).resolve().parents[1]
t=pd.read_csv(ROOT/'tables'/'21_herb_target_taxonomy_verified.tsv',sep='\t',dtype=str)
d=pd.read_csv(ROOT/'tables'/'17_herb_ingredient_details.tsv',sep='\t',dtype=str)
sq=json.loads((ROOT/'results'/'tcmsp_structure_qc_summary.json').read_text(encoding='utf-8'))
out={'time_utc':datetime.now(timezone.utc).isoformat(),'tcmsp_structure_qc':sq,'herb_saved_details':len(d),'herb_cached_details':len(list((ROOT/'raw').glob('hing_*.bin'))),'herb_total_cohort':2472,'human_target_details_taxon9606':int((t['Taxonomy id']=='9606').sum()),'unique_verified_gene_ids':t.loc[t['Taxonomy id']=='9606','Gene id'].nunique(),'genecards_status':'awaiting_user_manual_export; no HPO replacement','ongoing_process':'scripts/19_herb_collect_ingredients.py (unified exec session 79647)','next_actions':['finish all HERB ingredient details using existing cache','refresh only new PubChem CIDs with scripts/28_pubchem_herb_cids.py','normalize HERB structures with scripts/29_herb_structure_qc.py','refresh TCMSP/HERB crosswalk and missing taxonomy details with scripts/23_herb_tcmsp_crosswalk.py','assemble source-specific human edges and report before disease analysis']}
(ROOT/'results'/'阶段进度_续作断点.json').write_text(json.dumps(out,ensure_ascii=False,indent=2),encoding='utf-8');print(json.dumps(out,ensure_ascii=False,indent=2))
