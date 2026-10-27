#!/usr/bin/env python3
"""Reproduce offline ancillary file checks and fitting diagnostics already reported."""
import pathlib,json
import pandas as pd,numpy as np
p=pathlib.Path(__file__).resolve().parents[1]
c=pd.read_csv(p/'results/counts_all_genes.tsv.gz',sep='\t',index_col=0)
f=pd.read_csv(p/'results/fpkm_all_genes_audit_only.tsv.gz',sep='\t',index_col=0)
bins=pd.DataFrame({s:{'zeros':int((c[s]==0).sum()),'1':int((c[s]==1).sum()),'2_4':int(c[s].between(2,4).sum()),'5_9':int(c[s].between(5,9).sum()),'10_19':int(c[s].between(10,19).sum()),'20_99':int(c[s].between(20,99).sum()),'100plus':int((c[s]>=100).sum())} for s in c}).T
bins.to_csv(p/'results/low_count_distribution.tsv',sep='\t')
a=(c>=100).all(axis=1); q=f.loc[a]/c.loc[a]
logratio=np.log(q).sub(np.log(q).mean(axis=1),axis=0);resid=logratio.sub(logratio.median(axis=0),axis=1)
d={'genes_count_ge100_all_samples':int(a.sum()),'max_abs_log_relative_FPKM_count_ratio_deviation_after_sample_scaling':float(resid.abs().max().max()),'median_abs_deviation':float(resid.abs().median().median()),'count_positive_FPKM_zero':int(((c>0)&(f==0)).sum().sum()),'count_zero_FPKM_positive':int(((c==0)&(f>0)).sum().sum()),'interpretation':'Near-zero deviations support consistent gene-length factors and per-sample library scaling. This is internal consistency, not independent validation of raw counting.'}
(p/'results/fpkm_count_consistency.json').write_text(json.dumps(d,indent=2))
r=pd.read_csv(p/'results/differential_all_genes.tsv.gz',sep='\t');eligible=r[r.passed_prefilter]
flags=eligible[(eligible['_MAP_converged']==False)|(eligible['_genewise_converged']==False)]
flags.to_csv(p/'results/fit_convergence_flags.tsv',sep='\t',index=False)
qc=pd.read_csv(p/'results/sample_qc.tsv',sep='\t')
qc.groupby('condition')[['genes_count_ge1','genes_count_ge5','genes_count_ge10','zero_fraction','size_factor']].agg(['min','max','mean']).to_csv(p/'results/qc_group_descriptives.tsv',sep='\t')
print('Additional offline QC checks saved; no extra hypotheses tested.')
