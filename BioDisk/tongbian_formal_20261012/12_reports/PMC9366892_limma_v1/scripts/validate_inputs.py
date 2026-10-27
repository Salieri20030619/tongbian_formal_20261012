from pathlib import Path
import hashlib,json,csv,datetime,collections
import numpy as np,pandas as pd
r=Path(__file__).resolve().parent
base=r/'inputs/tongbian_human_PMC9366892_deidentified_inputs_20261012'
checks=[]
for line in (base/'SHA256SUMS.txt').read_text().splitlines():
 if not line.strip():continue
 expected,name=line.split(maxsplit=1);name=name.lstrip('*');p=base/name
 assert p.resolve().is_relative_to(base.resolve());actual=hashlib.sha256(p.read_bytes()).hexdigest();assert actual==expected
 checks.append({'file':name,'bytes':p.stat().st_size,'sha256':actual})
meta=pd.read_csv(base/'metadata_deidentified.csv');ids=meta.sample_id.tolist();assert len(ids)==len(set(ids))==60
assert set(ids)=={f'{g}{n}' for g in ['C','T'] for n in range(1,31)}
assert all((sid.startswith('C') and g=='Control') or (sid.startswith('T') and g=='STC') for sid,g in zip(meta.sample_id,meta.group))
matrices={};sample_alignment=[]
for name in ['metabolites_60x21_peak_area.csv','joint_selected_60x17.csv','microbial_selected_taxa_60x14.csv','secondary_bile_acids_60x7_peak_area.csv','alpha_diversity_60x7.csv']:
 df=pd.read_csv(base/name);assert len(df)==df.sample_id.nunique()==60 and set(df.sample_id)==set(ids)
 sample_alignment.append({'file':name,'same_order_as_metadata':df.sample_id.tolist()==ids,'matched_unique_ids':60})
 matrices[name]=df.set_index('sample_id').loc[ids]
x=matrices['metabolites_60x21_peak_area.csv'];ann=pd.read_csv(base/'metabolites_21_annotation.csv').set_index('ID');assert x.shape==(60,21) and set(x.columns)==set(ann.index)
assert np.isfinite(x.to_numpy()).all() and (x.to_numpy()>0).all()
assert abs(x.loc['C6','M407T218']-22146617.1867411)<1e-6
j=matrices['joint_selected_60x17.csv'];assert abs(j.loc['C6','Cholic acid']-146617.187)<1e-6
assert ann.loc['M498T146','Name']=='Taurochenodeoxycholate' and 'M500T143' in matrices['secondary_bile_acids_60x7_peak_area.csv'].columns
assert meta[['Age','Sex','BMI']].notna().all().all();assert set(meta.Sex)=={1,2}
designs=[]
for omit in [False,True]:
 m=meta.loc[meta.sample_id!='C6'] if omit else meta
 for adjusted in [False,True]:
  a=np.column_stack([np.ones(len(m)),(m.group=='STC').astype(int)])
  if adjusted:a=np.column_stack([a,m.Age,(m.Sex==2).astype(int),m.BMI])
  rank=int(np.linalg.matrix_rank(a));assert rank==a.shape[1]
  designs.append({'excluded_C6':omit,'adjusted':adjusted,'n':len(m),'columns':a.shape[1],'rank':rank,'residual_df':len(m)-rank,'group_n':m.group.value_counts().to_dict()})
tax=matrices['microbial_selected_taxa_60x14.csv'];duplicate_taxa=[]
for i,a in enumerate(tax):
 for b in list(tax)[i+1:]:
  if np.array_equal(tax[a].values,tax[b].values):duplicate_taxa.append([a,b])
audit={'time_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'zip_sha256':'634eba9e4c718d6051b9bce01695f306a574d4c2184cf3dbdd71dcc42b324bb0','library_file_id':'libfile_9025b7a56c008191971c58e1dd17241e','version':0,'inner_files_verified':checks,'matrix_alignment':sample_alignment,'raw_peak_area':{'shape':list(x.shape),'minimum':float(x.min().min()),'maximum':float(x.max().max()),'zero_cells':0,'missing_cells':0,'negative_cells':0,'log2_without_pseudocount':True,'scale_evidence':'source README and file name label raw peak area; magnitude and annotation arithmetic FC audit corroborate'},'metadata_core_missing':meta[['Age','Sex','BMI']].isna().sum().to_dict(),'sex_coding':'categorical code1/code2; biological meanings not inferred','design_checks':designs,'exact_duplicate_taxon_vectors':duplicate_taxa,'microbe_correlation_run':False,'C6_conflict_verified':{'feature_id':'M407T218','source_A':float(x.loc['C6','M407T218']),'source_B':float(j.loc['C6','Cholic acid']),'modified_input':False},'duplicate_name_rule':'M498T146 and M500T143 kept separate; only first is in selected21 analysis','names_present_in_package':'cloud package_validation confirmed zero original name matches; original workbook absent'}
out=r/'validation';out.mkdir(exist_ok=True)
(out/'local_input_validation.json').write_text(json.dumps(audit,indent=2),encoding='utf-8')
config={'frozen_before_fit':True,'time_utc':audit['time_utc'],'seed':20261012,'input_sha256':audit['zip_sha256'],'input_library_file_id':audit['library_file_id'],'transform':'log2(original strictly positive peak area); no pseudocount, additional normalization, imputation, or FC gate','primary_20':{'exclude_feature':'M407T218','subjects':'all60','fit':'lmFit + eBayes(trend=FALSE, robust=FALSE)','design':'~ group','reference':'Control','coefficient':'groupSTC','BH_family':20},'all21_original_record':{'C6_value':22146617.1867411,'flag':'source conflict; excluded from main20 summary','BH_family':21},'prespecified_cholic_scenarios':['source_A_21_metabolites','source_B_C6_only_from_joint','exclude_C6_entire_sample'],'sensitivity_design':'~ group + Age + factor(Sex_original_code) + BMI','sensitivity_covariates_excluded':['symptom flags','disease duration','Wexner','GIQLI','Bristol'], 'source_B_rule':'change only C6 M407T218 in an in-memory copy; retain original input and other original precision','primary_P':0.05,'effect':'STC - Control log2 peak area (logFC)','BH':'separately for each frozen fit family (20 primary or 21 full); no significance prescreen','sensitivity_BH':'cholic P/BH extracted from the corresponding all21 fit; no result-driven version selection','ppcor':'NOT_INSTALLED; no partial correlation or substitute algorithm','correlation_across_RNA_cohorts':False,'evidence_label':'selected feature reanalysis in same cohort; neither independent validation nor untargeted discovery'}
(r/'analysis_config_frozen.json').write_text(json.dumps(config,indent=2),encoding='utf-8')
print(json.dumps({'status':'PASS','inner_hashes_verified':len(checks),'raw_min':audit['raw_peak_area']['minimum'],'raw_max':audit['raw_peak_area']['maximum'],'designs':designs,'C6':audit['C6_conflict_verified'],'duplicate_taxon_pairs':len(duplicate_taxa)},indent=2))
