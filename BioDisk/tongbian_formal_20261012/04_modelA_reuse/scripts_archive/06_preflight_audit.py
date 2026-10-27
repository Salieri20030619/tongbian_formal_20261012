from direct_http import *
import re,pandas as pd
from bs4 import BeautifulSoup
catalog=pd.read_csv(ROOT/'tables'/'03_tcmsp_catalog.tsv',sep='\t')
print('RELATED CATALOG',catalog[catalog.herb_cn_name.str.contains('麦|地黄|生地|石斛|瓜蒌|栝楼')].to_dict('records'))
for name in ['H01','H02','H03','H16','H20']:
 text=(ROOT/'raw'/f'tcmsp_{name}.bin').read_text(encoding='utf-8')
 for match in re.finditer(r'\bdata:\s*(\[.*?\])\s*,',text,re.S):
  try:json.loads(match.group(1))
  except Exception as exc:
   print(name,'JSON ERROR',str(exc));print(match.group(1)[:300]);print('around error',match.group(1)[max(0,exc.pos-100):exc.pos+150])
text=(ROOT/'raw'/'tcmsp_H01.bin').read_text(encoding='utf-8')
for match in re.finditer('template:',text):
 print(text[match.start():match.start()+300])
hpo=pd.read_csv(ROOT/'raw'/'hpo.bin',sep='\t')
terms=['constipation','chronic constipation','functional constipation','chronic idiopathic constipation','slow transit constipation']
chosen=hpo[hpo.hpo_name.str.lower().isin(terms)]
chosen.to_csv(ROOT/'tables'/'07_HPO_constipation_CANDIDATE_NOT_PRIMARY.tsv',sep='\t',index=False,encoding='utf-8-sig')
out={'download_rows':len(hpo),'exact_term_match_rows':len(chosen),'unique_genes':chosen.ncbi_gene_id.nunique(),'hpo_terms':chosen[['hpo_id','hpo_name']].drop_duplicates().to_dict('records'),'per_term':chosen.groupby(['hpo_id','hpo_name']).agg(rows=('ncbi_gene_id','size'),genes=('ncbi_gene_id','nunique')).reset_index().to_dict('records'),'interpretation':'候选源覆盖预检，尚未获根审阅，不作为正式疾病基因集；与遗传综合征表型关联有偏倚'}
(ROOT/'results'/'hpo_candidate_coverage.json').write_text(json.dumps(out,ensure_ascii=False,indent=2),encoding='utf-8');print(json.dumps(out,ensure_ascii=False))
