#!/usr/bin/env python3
"""Offline audit and exploratory bulk RNA-seq analysis of GSE245885."""
import os
os.environ.setdefault('MPLCONFIGDIR','/tmp/tongbian_mpl');os.environ.setdefault('XDG_CACHE_HOME','/tmp/tongbian_cache')
os.environ.setdefault('OPENBLAS_NUM_THREADS','1');os.environ.setdefault('OMP_NUM_THREADS','1')
import pathlib,json,sys,platform,datetime,pickle,importlib.metadata
import numpy as np,pandas as pd
from scipy.stats import false_discovery_control
from sklearn.decomposition import PCA
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from pydeseq2.default_inference import DefaultInference
ROOT=pathlib.Path(__file__).resolve().parents[1]
meta=pd.read_csv(ROOT/'metadata/samples.tsv',sep='\t').set_index('sample')
counts=[];fpkm=[];annotation=None; audits=[]
for sample,r in meta.iterrows():
    d=pd.read_csv(ROOT/'raw'/r['supplement_file'],sep='\t',keep_default_na=False)
    assert d.columns.tolist()==['ensembl_id','gene_name',sample+'.count',sample+'.FPKM'], (sample,d.columns)
    assert not d.ensembl_id.duplicated().any(),sample
    assert d.ensembl_id.str.fullmatch('ENSG[0-9]+').all(),sample
    c=pd.to_numeric(d[sample+'.count'],errors='raise');f=pd.to_numeric(d[sample+'.FPKM'],errors='raise')
    assert np.isfinite(c).all() and (c>=0).all() and (c==np.floor(c)).all(),sample
    assert np.isfinite(f).all() and (f>=0).all(),sample
    a=d.set_index('ensembl_id')[['gene_name']]
    if annotation is None: annotation=a
    else: pd.testing.assert_frame_equal(annotation,a)
    counts.append(pd.Series(c.to_numpy(dtype=np.int64),index=d.ensembl_id,name=sample))
    fpkm.append(pd.Series(f.to_numpy(dtype=float),index=d.ensembl_id,name=sample))
    audits.append({'sample':sample,'gsm':r.gsm,'condition':r.condition,'rows':len(d),'count_sum':int(c.sum()),'detected_count_gt0':int((c>0).sum()),'count_ge10':int((c>=10).sum()),'zero_fraction':float((c==0).mean()),'noninteger_count_values':int((c!=np.floor(c)).sum()),'duplicate_ensembl_ids':int(d.ensembl_id.duplicated().sum()),'missing_count_values':int(c.isna().sum()),'max_count':int(c.max()),'top10_count_fraction':float(c.nlargest(10).sum()/c.sum())})
counts=pd.concat(counts,axis=1);fpkm=pd.concat(fpkm,axis=1)
assert counts.columns.equals(meta.index)
qc=pd.DataFrame(audits).set_index('sample')
# Counts and FPKM must be distinct; all FPKM columns have fractional values.
assert (fpkm.to_numpy()!=np.floor(fpkm.to_numpy())).any()
counts.to_csv(ROOT/'results/counts_all_genes.tsv.gz',sep='\t',compression='gzip')
fpkm.to_csv(ROOT/'results/fpkm_all_genes_audit_only.tsv.gz',sep='\t',compression='gzip')
annotation.to_csv(ROOT/'results/gene_annotation.tsv',sep='\t')
qc.to_csv(ROOT/'results/sample_qc.tsv',sep='\t')
keep=(counts>=10).sum(axis=1)>=5
# Original annotation is never collapsed for differential expression.
audit={'completed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'shape':[int(x) for x in counts.shape],'sample_order':counts.columns.tolist(),'all_count_values_nonnegative_integers':True,'all_gene_ids_and_symbols_identical_across_samples':True,'missing_numeric_values':0,'duplicate_ensembl_ids':0,'duplicate_gene_symbol_rows':int(annotation.gene_name.duplicated(keep=False).sum()),'empty_gene_symbols':int(annotation.gene_name.eq('').sum()),'original_genes':len(counts),'all_zero_genes':int(counts.sum(axis=1).eq(0).sum()),'eligible_genes':int(keep.sum()),'prefilter':'count >=10 in >=5 samples','min_library_size':int(qc.count_sum.min()),'max_library_size':int(qc.count_sum.max()),'no_samples_excluded':True,'raw_read_qc_available':False,'count_provenance':'GEO submission states HTSeq v0.13.5; no FASTQ re-quantification performed; integer checks validate file structure, not absolute provenance'}
(ROOT/'results/input_audit.json').write_text(json.dumps(audit,indent=2));print(json.dumps(audit,indent=2),flush=True);print(qc.to_string(),flush=True)
versions={'python':sys.version,'platform':platform.platform(),'started_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'packages':{n:importlib.metadata.version(n) for n in ['pydeseq2','anndata','numpy','scipy','pandas','scikit-learn','matplotlib','formulaic']},'parameters':{'design':'~condition','contrast':['condition','STC','Control'],'fit_type':'parametric','size_factors_fit_type':'ratio','refit_cooks':True,'min_replicates':7,'cooks_filter':True,'independent_filter':False,'n_cpus':2}}
(ROOT/'metadata/session_info.json').write_text(json.dumps(versions,indent=2))
metadata=meta[['condition']].copy();metadata['condition']=pd.Categorical(metadata.condition,categories=['Control','STC'])
inference=DefaultInference(n_cpus=2)
dds=DeseqDataSet(counts=counts.loc[keep].T,metadata=metadata,design='~condition',refit_cooks=True,fit_type='parametric',size_factors_fit_type='ratio',inference=inference)
dds.deseq2()
# Save fit before variance stabilization potentially refits dispersion trend.
with (ROOT/'results/pydeseq2_fit.pkl').open('wb') as f: pickle.dump(dds,f)
ds=DeseqStats(dds,contrast=['condition','STC','Control'],cooks_filter=True,independent_filter=False,inference=inference)
ds.summary();res=ds.results_df.copy();res['log2FC_CI95_low']=res.log2FoldChange-1.959963984540054*res.lfcSE;res['log2FC_CI95_high']=res.log2FoldChange+1.959963984540054*res.lfcSE
res=res.rename(columns={'padj':'BH_FDR_all_tested'})
res['cooks_filtered']=res.pvalue.isna()
# Retain unfiltered Wald evidence transparently as diagnostics, not as primary results.
dsraw=DeseqStats(dds,contrast=['condition','STC','Control'],cooks_filter=False,independent_filter=False,inference=inference,quiet=True);dsraw.summary()
res['pvalue_before_cooks_filter']=dsraw.results_df.pvalue
res['BH_FDR_before_cooks_filter']=dsraw.results_df.padj
for col in ['_LFC_converged','_MAP_converged','_genewise_converged','_outlier_genes','_pvalue_cooks_outlier']:
    if col in dds.var: res[col]=dds.var[col]
res['test_status']=np.where(res.pvalue.isna(),'cooks_filtered_or_not_estimable','tested')
full=annotation.join(res);full['passed_prefilter']=keep;full['test_status']=full.test_status.fillna('low_count_prefilter')
full.to_csv(ROOT/'results/differential_all_genes.tsv.gz',sep='\t',compression='gzip')
full.loc[keep].sort_values('pvalue').to_csv(ROOT/'results/differential_tested_genes.tsv',sep='\t')
qc['size_factor']=dds.obs['size_factors'];qc.to_csv(ROOT/'results/sample_qc.tsv',sep='\t')
norm=pd.DataFrame(dds.layers['normed_counts'].T,index=dds.var_names,columns=dds.obs_names)
norm.to_csv(ROOT/'results/normalized_counts_tested.tsv.gz',sep='\t',compression='gzip')
# Blind VST for QC and predetermined gene-set score analysis; no condition fitted here.
dds.vst(use_design=False)
vst=pd.DataFrame(dds.layers['vst_counts'].T,index=dds.var_names,columns=dds.obs_names)
vst.to_csv(ROOT/'results/blind_vst_tested.tsv.gz',sep='\t',compression='gzip')
for name,x in [('all',vst),('top500',vst.loc[vst.var(axis=1).nlargest(min(500,len(vst))).index])]:
    pca=PCA(n_components=5,svd_solver='full');coords=pca.fit_transform(x.T)
    pc=pd.DataFrame(coords,index=x.columns,columns=[f'PC{i+1}' for i in range(5)]).join(meta[['gsm','condition']])
    pc.to_csv(ROOT/f'results/pca_{name}_coordinates.tsv',sep='\t')
    (ROOT/f'results/pca_{name}_variance.json').write_text(json.dumps(pca.explained_variance_ratio_.tolist()))
vst.corr(method='pearson').to_csv(ROOT/'results/sample_pearson_vst.tsv',sep='\t')
counts.corr(method='spearman').to_csv(ROOT/'results/sample_spearman_rawcounts.tsv',sep='\t')
summary={'completed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'eligible_genes':len(res),'finite_primary_p':int(res.pvalue.notna().sum()),'cooks_filtered_or_na':int(res.pvalue.isna().sum()),'nominal_p_lt_0.05':int((res.pvalue<.05).sum()),'nominal_p_lt_0.05_up':int(((res.pvalue<.05)&(res.log2FoldChange>0)).sum()),'nominal_p_lt_0.05_down':int(((res.pvalue<.05)&(res.log2FoldChange<0)).sum()),'BH_lt_0.05':int((res.BH_FDR_all_tested<.05).sum()),'BH_lt_0.1':int((res.BH_FDR_all_tested<.1).sum()),'minimum_p':float(res.pvalue.min()),'minimum_BH_FDR':float(res.BH_FDR_all_tested.min()),'actual_dispersion_fit_type':dds.uns.get('disp_function_type','not available')}
(ROOT/'results/model_summary.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary,indent=2),flush=True)
print('OFFLINE MODEL COMPLETE',flush=True)
