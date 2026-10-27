from direct_http import *
import re,pandas as pd
coverage=pd.read_csv(ROOT/'tables'/'04_tcmsp_herb_coverage.tsv',sep='\t').to_dict('records')
ingredients=[];targets=[];disease_edges=[];arrays_by_herb={}
for cr in coverage:
 if cr['collection_status']!='public_grid_downloaded':continue
 text=(ROOT/cr['raw_file']).read_text(encoding='utf-8');arrays=[]
 for m in re.finditer(r'\bdata:\s*(?=\[)',text):
  arrays.append(json.JSONDecoder().raw_decode(text[m.end():])[0])
 arrays_by_herb[cr['herb_id']]=arrays
 for arr in arrays:
  if not arr:continue
  extra={'herb_id':cr['herb_id'],'prescription_label':cr['prescription_label'],'source_url':cr['source_url'],'source_raw':cr['raw_file']}
  if 'ob' in arr[0]:
   ingredients.extend([dict(x,**extra) for x in arr]);cr['ingredient_rows']=len(arr)
  elif 'target_name' in arr[0] and 'molecule_ID' in arr[0]:
   targets.extend([dict(x,**extra) for x in arr]);cr['target_rows']=len(arr)
  elif 'target_name' in arr[0] and 'disease_name' in arr[0]:disease_edges.extend([dict(x,**extra) for x in arr])
for name,data in [('04_tcmsp_herb_coverage',coverage),('05_tcmsp_ingredient_rows',ingredients),('06_tcmsp_target_rows_UNCLASSIFIED',targets),('08_tcmsp_disease_edges_NOT_USED',disease_edges)]:
 pd.DataFrame(data).to_csv(ROOT/'tables'/f'{name}.tsv',sep='\t',index=False,encoding='utf-8-sig')
(ROOT/'results'/'tcmsp_grid_arrays.json').write_text(json.dumps(arrays_by_herb,ensure_ascii=False,indent=2),encoding='utf-8')
i=pd.DataFrame(ingredients);t=pd.DataFrame(targets)
summary={'herb_name_matches':len(arrays_by_herb),'herb_denominator':20,'herb_ingredient_rows':len(i),'unique_tcmsp_mol_ids_NOT_structure_dedup':i.MOL_ID.nunique(),'herb_compound_target_rows_UNCLASSIFIED':len(t),'unique_compound_target_pairs_all_species_unverified':len(t[['MOL_ID','target_ID']].drop_duplicates()),'unique_target_ids_all_species_unverified':t.target_ID.nunique(),'database_validated_label_rows':int((t.validated=='validated').sum()),'species_status':'未映射；表中包含细菌等非人靶点，不能把这些计数当人靶点网络','structure_deduplicated_parent_compounds':None,'human_targets':None,'disease_intersection':None,'ppi_nodes':None,'ppi_edges':None,'go_kegg':None}
(ROOT/'results'/'tcmsp_acquisition_summary.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2),encoding='utf-8');print(json.dumps(summary,ensure_ascii=False,indent=2))
