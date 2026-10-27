from pathlib import Path
import pandas as pd,re,json
ROOT=Path(__file__).resolve().parents[1]
text=(ROOT/'raw'/'hpo_ontology.bin').read_text(encoding='utf-8');nodes={}
for block in text.split('[Term]')[1:]:
 mid=re.search(r'^id: (HP:\d+)',block,re.M);name=re.search(r'^name: (.+)',block,re.M)
 if mid:nodes[mid[1]]={'name':name[1] if name else '', 'parents':re.findall(r'^is_a: (HP:\d+)',block,re.M),'obsolete':'is_obsolete: true' in block}
desc={'HP:0002019'}
while True:
 new=desc|{k for k,v in nodes.items() if set(v['parents'])&desc and not v['obsolete']}
 if new==desc:break
 desc=new
a=pd.read_csv(ROOT/'raw'/'hpo_phenotype_hpoa.bin',sep='\t',comment='#',dtype=str).fillna('')
subset=a[a.hpo_id.isin(desc)].copy()
subset.to_csv(ROOT/'tables'/'24_HPO_constipation_direct_disease_annotations_NOT_PRIMARY.tsv',sep='\t',index=False,encoding='utf-8-sig')
pos=subset[subset.qualifier!='NOT']; disease_terms=pos.groupby('database_id').hpo_id.agg(set).to_dict()
cand=pd.read_csv(ROOT/'tables'/'07_HPO_constipation_CANDIDATE_NOT_PRIMARY.tsv',sep='\t',dtype=str)
def classify(row):
 terms=disease_terms.get(row.disease_id,set())
 if row.hpo_id in terms:return 'same_term_direct_disease_annotation_present'
 if terms&desc:return 'descendant_disease_annotation_only_in_snapshot'
 return 'not_reconciled_in_downloaded_annotation_snapshot'
cand['disease_annotation_support']=cand.apply(classify,axis=1)
cand.to_csv(ROOT/'tables'/'25_HPO_candidate_propagation_audit_NOT_PRIMARY.tsv',sep='\t',index=False,encoding='utf-8-sig')
out={'descendant_terms_including_root':[dict(hpo_id=k,name=nodes[k]['name']) for k in sorted(desc)],'hpoa_metadata_lines':(ROOT/'raw'/'hpo_phenotype_hpoa.bin').read_text(encoding='utf-8').splitlines()[:4],'annotation_rows_by_term':subset.groupby(['hpo_id','qualifier']).size().to_dict().__str__(),'candidate_row_support_counts':cand.disease_annotation_support.value_counts().to_dict(),'conclusion':'phenotype_to_genes的精确名称匹配不能解读为直接注释；这里以疾病ID回查原始HPOA，区分直接同词支持与仅后代词支持。未将其用于主疾病分析。'}
(ROOT/'results'/'hpo_annotation_audit.json').write_text(json.dumps(out,ensure_ascii=False,indent=2),encoding='utf-8');print(json.dumps(out,ensure_ascii=False,indent=2))
