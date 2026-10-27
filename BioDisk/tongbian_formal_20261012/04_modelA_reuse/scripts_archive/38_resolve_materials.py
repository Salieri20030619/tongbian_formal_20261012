"""Resolve two medicinal-material aliases against downloaded official HKCMMS standards."""
from pathlib import Path
import pandas as pd,json,hashlib
ROOT=Path(__file__).resolve().parents[1]
def txt(x):return x.get('title','') if isinstance(x,dict) else x
rows=[];audit=[]
for hid,dbid,pdf,page,url,note in [
 ('H13','HERB001010','official_ophiopogon_standard',168,'https://www.cmro.gov.hk/files/hkcmms/vol3/pdf_e/Radix_Ophiopogonis_v3_e.pdf','HKCMMS vol3 printed p168 defines Maidong/Radix Ophiopogonis as dried root tuber of Ophiopogon japonicus. HERB Cundong explicitly records Radix Ophiopogonis and Ophiopogon japonicus; use-part field is misfiled as species, recovered using the official medicinal-drug standard. No batch processing or chemical exposure confirmed.'),
 ('H15','HERB001783','official_rehmannia_standard',236,'https://www.cmro.gov.hk/files/hkcmms/vol3/pdf_e/Radix_Rehmanniae_v3_e.pdf','HKCMMS vol3 printed p236 defines Radix Rehmanniae as dried root tuber of Rehmannia glutinosa, commonly called Sheng Dihuang. HERB Gan Dihuang explicitly records Rehmannia glutinosa and dried root. Fresh and prepared Rehmannia records remain excluded; batch processing is not measured.')]:
  raw=f'raw/herb_detail_{dbid}.bin';d=json.loads((ROOT/raw).read_text(encoding='utf-8'))
  for row in d['herb_ingredient'][1:]:
   r=dict(zip([txt(x) for x in d['herb_ingredient'][0]],[txt(x) for x in row]))
   rows.append({'herb_id':hid,'source_compound_id':r['Ingredient id'],'compound_name':r['Ingredient name'],'source_material_id':dbid,'source_url':f'http://herb.ac.cn/Detail/?v={dbid}&label=Herb','source_raw':raw,'primary_material_eligible':True,'resolution':'official_standard_medical_drug_synonym_resolved','resolution_evidence':note,'official_standard_url':url,'official_standard_pdf_page':2,'official_standard_printed_page':page})
  audit.append({'herb_id':hid,'source_material_id':dbid,'record_count':len(d['herb_ingredient'])-1,'evidence':note,'source_url':url,'source_raw':f'raw/{pdf}.bin','source_sha256':hashlib.sha256((ROOT/f'raw/{pdf}.bin').read_bytes()).hexdigest()})
pd.DataFrame(rows).to_csv(ROOT/'tables/28_material_resolution_memberships.tsv',sep='\t',index=False,encoding='utf-8-sig')
(ROOT/'results/material_resolution_audit.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps(audit,ensure_ascii=False,indent=2))
