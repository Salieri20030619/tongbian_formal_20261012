"""Assemble a traceable parent-compound baseline; no disease-directed ingredient selection."""
from pathlib import Path
import pandas as pd,json,re,hashlib
from datetime import datetime,timezone
ROOT=Path(__file__).resolve().parents[1]
def read(n):return pd.read_csv(ROOT/'tables'/n,sep='\t',dtype=str,keep_default_na=False)
formula=read('02_formula_standardization.tsv')
tc=read('26_tcmsp_structure_qc.tsv');hc=read('27_herb_structure_qc.tsv')
if len(hc)!=2472:raise SystemExit('HERB cohort incomplete; do not publish complete baseline')
tc_ok=tc[tc.network_structure_eligible=='True'];hc_ok=hc[hc.network_structure_eligible=='True']
tm=read('05_tcmsp_ingredient_rows.tsv');hm=read('16_herb_memberships_eligible.tsv')
members=[];excluded=[]
for r in tm.to_dict('records'):
 q=tc_ok[tc_ok.MOL_ID==r['MOL_ID']]
 status='accepted_database_drug_record';reason='药名匹配；源记录未区分生品炮制，非制剂实测'
 if r['herb_id'] in ['H02','H07']:
  status='exclude_botanical_label_conflict';reason='TCMSP拉丁物种与处方药材冲突；保留HERB正确药材条目的独立归属'
 if q.empty:excluded.append(dict(r,exclusion='structure_not_eligible',material_status=status));continue
 if status.startswith('exclude'):excluded.append(dict(r,exclusion=reason,material_status=status));continue
 q=q.iloc[0];members.append({'herb_id':r['herb_id'],'compound_key':q.computed_inchikey,'compound_name':r['molecule_name'],'source':'TCMSP','source_compound_id':r['MOL_ID'],'source_material_id':r['herb_id'],'source_url':r['source_url'],'source_raw':r['source_raw'],'material_status':status,'material_note':reason,'structure_status':q.structure_status,'structure_source_raw':q.get('structure_source_raw',''),'smiles':q.canonical_isomeric_smiles})
for r in hm.to_dict('records'):
 q=hc_ok[hc_ok.herb_ingredient_id==r['Ingredient id']]
 if r['material_status'].startswith('conditional_'):
  excluded.append(dict(r,exclusion=r['material_status']));continue
 if q.empty:excluded.append(dict(r,exclusion='structure_not_eligible'));continue
 q=q.iloc[0];members.append({'herb_id':r['herb_id'],'compound_key':q.computed_inchikey,'compound_name':r['Ingredient name'],'source':'HERB','source_compound_id':r['Ingredient id'],'source_material_id':r['herb_db_id'],'source_url':r['source_url'],'source_raw':r['source_raw'],'material_status':r['material_status'],'material_note':r['material_note'],'structure_status':q.structure_status,'structure_source_raw':q.get('pubchem_source_raw','') or q.source_raw,'smiles':q.canonical_isomeric_smiles})
# Add independently verified medicinal-drug synonym records only after their documented resolution.
extra=ROOT/'tables'/'28_material_resolution_memberships.tsv'
if extra.exists():
 ex=pd.read_csv(extra,sep='\t',dtype=str,keep_default_na=False)
 for r in ex[ex.primary_material_eligible=='True'].to_dict('records'):
  q=hc_ok[hc_ok.herb_ingredient_id==r['source_compound_id']]
  if q.empty:excluded.append(dict(r,exclusion='structure_not_eligible'));continue
  q=q.iloc[0];members.append({'herb_id':r['herb_id'],'compound_key':q.computed_inchikey,'compound_name':r['compound_name'],'source':'HERB','source_compound_id':r['source_compound_id'],'source_material_id':r['source_material_id'],'source_url':r['source_url'],'source_raw':r['source_raw'],'material_status':r['resolution'],'material_note':r['resolution_evidence'],'structure_status':q.structure_status,'structure_source_raw':q.get('pubchem_source_raw','') or q.source_raw,'smiles':q.canonical_isomeric_smiles})
m=pd.DataFrame(members).drop_duplicates(['herb_id','source','source_compound_id','source_material_id'])
# Final Model A audit: source derivative-tagged records cannot be assumed to be
# intact medicinal-material constituents. This identity rule is independent of
# target overlap/enrichment and never treats the derivatives as microbial input.
derivative=m.compound_name.str.contains(r'_qt\b',case=False,regex=True)
for r in m[derivative].to_dict('records'):
 excluded.append(dict(r,exclusion='source_qt_derivative_tag_no_independent_intact_parent_membership'))
m=m[~derivative].copy()
m.to_csv(ROOT/'tables'/'30_primary_herb_compound_evidence.tsv',sep='\t',index=False,encoding='utf-8-sig')
pd.DataFrame(excluded).to_csv(ROOT/'tables'/'31_primary_membership_exclusions.tsv',sep='\t',index=False,encoding='utf-8-sig')
keys=set(m.compound_key)
comp=m.groupby('compound_key',sort=True).agg(compound_names=('compound_name',lambda x:' | '.join(sorted(set(x)))),herbs=('herb_id',lambda x:';'.join(sorted(set(x)))),source_compound_ids=('source_compound_id',lambda x:';'.join(sorted(set(x)))),sources=('source',lambda x:';'.join(sorted(set(x)))),canonical_isomeric_smiles=('smiles','first')).reset_index()
comp.to_csv(ROOT/'tables'/'32_primary_unique_parent_structures.tsv',sep='\t',index=False,encoding='utf-8-sig')
cross=read('20_tcmsp_herb_target_crosswalk_evidence.tsv')[['tcmsp_target_id','herb_target_id']].drop_duplicates()
tax=read('21_herb_target_taxonomy_verified.tsv');tax=tax[(tax['Taxonomy id']=='9606')&(tax['Type of gene']=='protein-coding')]
cross=cross.merge(tax[['herb_target_id','Gene id','Gene name','Taxonomy id','source_raw','source_url']],on='herb_target_id').rename(columns={'Gene id':'ncbi_gene_id','Gene name':'source_gene_symbol','source_raw':'taxonomy_source_raw','source_url':'taxonomy_source_url'})
ambiguous=cross.groupby('tcmsp_target_id').ncbi_gene_id.nunique();bad=set(ambiguous[ambiguous>1].index)
cross[cross.tcmsp_target_id.isin(bad)].to_csv(ROOT/'tables'/'33_target_crosswalk_ambiguities.tsv',sep='\t',index=False,encoding='utf-8-sig')
cross=cross[~cross.tcmsp_target_id.isin(bad)].drop_duplicates(['tcmsp_target_id','ncbi_gene_id'])
edges=read('06_tcmsp_target_rows_UNCLASSIFIED.tsv').drop_duplicates(['MOL_ID','target_ID','validated','SVM_score','RF_score'])
edges=edges.merge(tc_ok[['MOL_ID','computed_inchikey']],on='MOL_ID').rename(columns={'computed_inchikey':'compound_key'})
edges=edges[edges.compound_key.isin(keys)]
edges=edges.merge(cross,left_on='target_ID',right_on='tcmsp_target_id')
registry=pd.read_csv(ROOT/'raw'/'hgnc_complete_set.bin',sep='\t',dtype=str,keep_default_na=False)
registry=registry[(registry.status=='Approved')&registry.entrez_id.ne('')][['hgnc_id','symbol','entrez_id','locus_group']].drop_duplicates('entrez_id')
edges=edges.merge(registry,left_on='ncbi_gene_id',right_on='entrez_id',how='left').fillna('')
edges['approved_symbol']=edges['symbol'].where(edges['symbol'].ne(''),edges.source_gene_symbol)
edges['evidence_class']=edges.validated.map(lambda x:'TCMSP_existing_validated_label_NOT_revalidated' if x=='validated' else 'TCMSP_existing_SVM_RF_prediction')
edges['new_prediction_run']=False
edges['species_evidence']='HERB explicit TCMSP target cross-reference + Taxonomy id 9606; original assay species not independently confirmed'
edges.to_csv(ROOT/'tables'/'34_primary_compound_human_target_evidence.tsv',sep='\t',index=False,encoding='utf-8-sig')
pairs=edges[['compound_key','ncbi_gene_id','approved_symbol']].drop_duplicates()
targets=edges[['ncbi_gene_id','approved_symbol','hgnc_id']].drop_duplicates()
targets.to_csv(ROOT/'tables'/'35_primary_human_targets.tsv',sep='\t',index=False,encoding='utf-8-sig')
gc=read('GeneCards_human_gene_mapping.tsv');gset=set(gc.loc[gc.mapping_status.str.startswith('unique_')&gc.ncbi_gene_id.ne(''),'ncbi_gene_id'])
inter=targets[targets.ncbi_gene_id.isin(gset)].copy()
inter=inter.merge(gc[['ncbi_gene_id','input_symbol','source_score_unfiltered','genecards_type']].drop_duplicates('ncbi_gene_id'),on='ncbi_gene_id')
inter.to_csv(ROOT/'tables'/'36_GeneCards_primary_target_intersection.tsv',sep='\t',index=False,encoding='utf-8-sig')
cov=[]
for r in formula.to_dict('records'):
 mm=m[m.herb_id==r['herb_id']];ck=set(mm.compound_key);ep=pairs[pairs.compound_key.isin(ck)];ip=ep[ep.ncbi_gene_id.isin(gset)]
 cov.append(dict(herb_id=r['herb_id'],prescription_label=r['prescription_label'],relative_input=r['relative_input'],raw_tcmsp_memberships=int((tm.herb_id==r['herb_id']).sum()),raw_herb_memberships=int((hm.herb_id==r['herb_id']).sum()),primary_structure_count=len(ck),primary_target_linked_structures=ep.compound_key.nunique(),primary_human_targets=ep.ncbi_gene_id.nunique(),disease_intersection_targets=ip.ncbi_gene_id.nunique(),primary_herb_compound_evidence_rows=len(mm)))
pd.DataFrame(cov).to_csv(ROOT/'tables'/'37_formula_primary_coverage.tsv',sep='\t',index=False,encoding='utf-8-sig')
summary={'time_utc':datetime.now(timezone.utc).isoformat(),'formula_herbs':20,'primary_herbs_with_structures':m.herb_id.nunique(),'primary_herbs_with_human_targets':sum(r['primary_human_targets']>0 for r in cov),'primary_unique_parent_structures':len(keys),'primary_herb_compound_unique_pairs':len(m[['herb_id','compound_key']].drop_duplicates()),'primary_target_linked_structures':pairs.compound_key.nunique(),'primary_human_target_genes':len(targets),'primary_compound_human_gene_unique_pairs':len(pairs),'source_evidence_rows':len(edges),'evidence_class_unique_pair_counts':edges.groupby('evidence_class').apply(lambda x:len(x[['compound_key','ncbi_gene_id']].drop_duplicates()),include_groups=False).to_dict(),'GeneCards_unique_source_symbols':len(gc),'GeneCards_unique_mapped_ncbi_genes':len(gset),'disease_intersection_genes':len(inter),'intersection_fraction_of_primary_targets':len(inter)/len(targets) if len(targets) else None,'human_identity_verified_before_structure_scope_filter':tax['Gene id'].nunique(),'unresolved_tcmsp_target_ids_excluded':805-cross.tcmsp_target_id.nunique(),'scope_note':'Primary target evidence restricted to TCMSP existing labelled validated/predicted edges with explicit human ID crosswalk. HERB unclassified/reference-mining target associations retained separately; not silently called experimental binding or new predictions.','GeneCards_filter':'all valid exported symbols; no score cutoff','go_kegg_status':'not_yet_run','ppi_status':'not_yet_run'}
(ROOT/'results'/'primary_baseline_summary.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2),encoding='utf-8');print(json.dumps(summary,ensure_ascii=False,indent=2))
