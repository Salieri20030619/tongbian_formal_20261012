#!/usr/bin/env python3
"""Offline fixed-set scores, exact sample permutation, and publication-style plots."""
import os
os.environ.setdefault('MPLCONFIGDIR','/tmp/tongbian_mpl');os.environ.setdefault('XDG_CACHE_HOME','/tmp/tongbian_cache')
import pathlib,json,itertools
import numpy as np,pandas as pd
from scipy.stats import false_discovery_control,t
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
ROOT=pathlib.Path(__file__).resolve().parents[1]
meta=pd.read_csv(ROOT/'metadata/samples.tsv',sep='\t').set_index('sample')
counts=pd.read_csv(ROOT/'results/counts_all_genes.tsv.gz',sep='\t',index_col=0)
qc=pd.read_csv(ROOT/'results/sample_qc.tsv',sep='\t',index_col=0)
vst=pd.read_csv(ROOT/'results/blind_vst_tested.tsv.gz',sep='\t',index_col=0)
res=pd.read_csv(ROOT/'results/differential_all_genes.tsv.gz',sep='\t',index_col=0)
modules=json.loads((ROOT/'genesets/host_modules_frozen.json').read_text())['modules']
# Equal-weight average of VST profiles for duplicate symbols, then each symbol receives one weight.
# Neither aggregation nor standardization uses condition labels or differential results.
symbol_vst=vst.join(res[['gene_name']]).groupby('gene_name',sort=True).mean()
z=symbol_vst.sub(symbol_vst.mean(axis=1),axis=0).div(symbol_vst.std(axis=1,ddof=1),axis=0)
finite=np.isfinite(z).all(axis=1);z=z.loc[finite]
stc=meta.condition.eq('STC').to_numpy();ctrl=~stc
partitions=list(itertools.combinations(range(len(meta)),5));assert len(partitions)==252
score=pd.DataFrame(index=meta.index);rows=[];mapping=[];member_results=[];loo=[]
labels={'WATER_SALT_TRANSPORT':'Water / Na / Cl transport','SMOOTH_MUSCLE_CONTRACTION':'Smooth muscle contraction','ENTERIC_NERVOUS_SYSTEM_DEVELOPMENT':'Enteric nervous system development','BICELLULAR_TIGHT_JUNCTION_ASSEMBLY':'Bicellular tight-junction assembly'}
for key,m in modules.items():
    wanted=m['gene_symbols'];matched=sorted(set(wanted)&set(z.index));assert len(matched)>=5
    score[key]=z.loc[matched].mean(axis=0)
    for gene in wanted:
        original=res.index[res.gene_name.eq(gene)].tolist();tested=[g for g in original if g in vst.index]
        mapping.append({'module':key,'symbol':gene,'in_original_matrix':bool(original),'passed_count_filter':bool(tested),'in_score':gene in matched,'n_original_ensembl_ids':len(original),'n_tested_ensembl_ids':len(tested),'ensembl_ids':';'.join(original),'aggregation':'equal_weight_mean_VST_for_duplicate_symbols_then_zscore_all_10_samples'})
        for ens in original:
            rr=res.loc[ens].to_dict();rr.update({'module':key,'ensembl_id':ens,'symbol':gene,'in_module_score':gene in matched});member_results.append(rr)
    s=score[key].to_numpy();effect=s[stc].mean()-s[ctrl].mean()
    diffs=[]
    for inds in partitions:
        mask=np.zeros(10,dtype=bool);mask[list(inds)]=True
        diffs.append(s[mask].mean()-s[~mask].mean())
    p=float(np.mean(np.abs(diffs)>=abs(effect)-1e-12))
    v1=s[stc].var(ddof=1)/5;v0=s[ctrl].var(ddof=1)/5;se=np.sqrt(v1+v0);df=(v1+v0)**2/((v1*v1+v0*v0)/4);ci=t.ppf(.975,df)*se
    # LOO is descriptive and not refitting/transformation selection or an additional hypothesis test.
    looe=[]
    for i,sample in enumerate(meta.index):
        keep=np.arange(10)!=i;value=s[stc&keep].mean()-s[ctrl&keep].mean();looe.append(value);loo.append({'module':key,'left_out':sample,'score_difference_STC_minus_control':value})
    gm=res[res.gene_name.isin(matched)&res.passed_prefilter].groupby('gene_name').log2FoldChange.mean()
    rows.append({'module':key,'label':labels[key],'source_ids':';'.join(m['source_ids']),'n_source_symbols':len(wanted),'n_detected_symbols':len(set(wanted)&set(res.gene_name)),'n_scored_symbols':len(matched),'fraction_scored':len(matched)/len(wanted),'control_mean_score':float(s[ctrl].mean()),'STC_mean_score':float(s[stc].mean()),'difference_STC_minus_control':float(effect),'descriptive_Welch_CI95_low':effect-ci,'descriptive_Welch_CI95_high':effect+ci,'exact_two_sided_p_252_permutations':p,'median_member_log2FC':float(gm.median()),'fraction_member_log2FC_positive':float((gm>0).mean()),'LOO_min_score_difference':float(min(looe)),'LOO_max_score_difference':float(max(looe)),'all_LOO_same_direction':bool(all(np.sign(e)==np.sign(effect) for e in looe))})
summary=pd.DataFrame(rows);summary['BH_FDR_four_modules']=false_discovery_control(summary.exact_two_sided_p_252_permutations)
summary.to_csv(ROOT/'results/module_results_all_four.tsv',sep='\t',index=False)
score.join(meta[['gsm','condition']]).to_csv(ROOT/'results/module_scores_each_sample.tsv',sep='\t')
pd.DataFrame(mapping).to_csv(ROOT/'results/module_gene_mapping_audit.tsv',sep='\t',index=False)
pd.DataFrame(member_results).to_csv(ROOT/'results/module_member_all_gene_results.tsv',sep='\t',index=False)
pd.DataFrame(loo).to_csv(ROOT/'results/module_score_leave_one_out.tsv',sep='\t',index=False)
# Extra requested QC, descriptive only.
for n in [1,5,10,20,100]:qc['genes_count_ge'+str(n)]=(counts>=n).sum(axis=0)
qc.to_csv(ROOT/'results/sample_qc.tsv',sep='\t')
print(summary.to_string(index=False),flush=True)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'figure.dpi':150,'savefig.dpi':300,'pdf.fonttype':42,'ps.fonttype':42})
blue='#3575A6';orange='#D66B36';grey='#A7AEB5';colors=[orange if g else blue for g in stc]
legend=[Line2D([0],[0],marker='o',lw=0,color=blue,label='Cancer-surgery control (n=5)'),Line2D([0],[0],marker='o',lw=0,color=orange,label='Severe STC (n=5)')]
def save(fig,name):
    fig.savefig(ROOT/'figures'/f'{name}.png',bbox_inches='tight');fig.savefig(ROOT/'figures'/f'{name}.pdf',bbox_inches='tight');plt.close(fig)
# Fig 1: QC. All-sample unsupervised PCA uses all prefiltered genes.
fig,ax=plt.subplots(2,2,figsize=(11.5,8.7));x=np.arange(10)
ax[0,0].bar(x,qc.count_sum/1e6,color=colors,width=.73);ax[0,0].set_xticks(x,meta.index,rotation=45,ha='right');ax[0,0].set_ylabel('Assigned gene counts (million)');ax[0,0].set_title('A  Count library sizes',loc='left',fontweight='bold')
for n,marker,color in [(1,'o','#3A4C62'),(5,'s','#B0914C'),(10,'^','#59866F')]:ax[0,1].plot(x,qc['genes_count_ge'+str(n)]/1000,marker=marker,color=color,label=f'Count >= {n}',lw=1.2)
ax[0,1].set_xticks(x,meta.index,rotation=45,ha='right');ax[0,1].set_ylabel('Detected genes (thousand)');ax[0,1].legend(frameon=False,fontsize=9);ax[0,1].set_title('B  Genes detected across count thresholds',loc='left',fontweight='bold')
pc=pd.read_csv(ROOT/'results/pca_all_coordinates.tsv',sep='\t',index_col=0);v=json.loads((ROOT/'results/pca_all_variance.json').read_text())
ax[1,0].scatter(pc.PC1,pc.PC2,c=colors,s=65,edgecolors='white',linewidth=.7)
for sample,r in pc.iterrows():ax[1,0].annotate(sample,(r.PC1,r.PC2),xytext=(5,4),textcoords='offset points',fontsize=8)
ax[1,0].set_xlabel(f'PC1 ({v[0]*100:.1f}% variance)');ax[1,0].set_ylabel(f'PC2 ({v[1]*100:.1f}% variance)');ax[1,0].set_title('C  Blind VST PCA: 20,479 filtered genes',loc='left',fontweight='bold');ax[1,0].margins(.22)
corr=pd.read_csv(ROOT/'results/sample_pearson_vst.tsv',sep='\t',index_col=0)
im=ax[1,1].imshow(corr,cmap='viridis',vmin=float(corr.to_numpy().min()),vmax=1,aspect='auto');ax[1,1].set_xticks(x,meta.index,rotation=45,ha='right');ax[1,1].set_yticks(x,meta.index);ax[1,1].set_title('D  Sample correlation (blind VST)',loc='left',fontweight='bold');cb=fig.colorbar(im,ax=ax[1,1],fraction=.046,pad=.04);cb.set_label('Pearson r')
fig.legend(handles=legend,loc='upper center',ncol=2,frameon=False,bbox_to_anchor=(.5,1.01));fig.suptitle('GSE245885: sample quality and transcriptomic variation',y=1.05,fontsize=14,fontweight='bold');fig.tight_layout(h_pad=2.0,w_pad=2.0);save(fig,'Figure1_QC_PCA')
# Fig 2: transcriptome-wide effect/p landscape, with all tests transparent.
d=res[res.test_status.eq('tested')].copy();nom=d.pvalue<.05
fig,ax=plt.subplots(1,2,figsize=(11.7,4.7))
ax[0].scatter(np.log10(d.baseMean+1),d.log2FoldChange,s=4,c=np.where(nom,np.where(d.log2FoldChange>0,orange,blue),grey),alpha=.45,linewidth=0,rasterized=True);ax[0].axhline(0,color='#444444',lw=.8);ax[0].set_xlabel('log10(mean normalized count + 1)');ax[0].set_ylabel('STC / control log2 fold change');ax[0].set_title('A  Effect size across all tested genes',loc='left',fontweight='bold')
ax[1].scatter(d.log2FoldChange,-np.log10(d.pvalue.clip(lower=1e-300)),s=4,c=np.where(nom,np.where(d.log2FoldChange>0,orange,blue),grey),alpha=.45,linewidth=0,rasterized=True);ax[1].axhline(-np.log10(.05),ls='--',c='#444444',lw=.8);ax[1].axvline(0,color='#444444',lw=.8);ax[1].set_xlabel('STC / control log2 fold change');ax[1].set_ylabel('-log10(nominal Wald p)');ax[1].set_title('B  Genome-wide Wald test results',loc='left',fontweight='bold')
fig.suptitle('PyDESeq2 negative-binomial model | 5 STC vs 5 cancer-surgery controls',fontweight='bold',y=1.04)
fig.text(.5,-.035,'20,479 eligible genes; 20,373 finite primary tests; 106 Cook-filtered. Colors: nominal p < 0.05, without a fold-change gate.\n12,476 nominal p < 0.05; 11,747 BH FDR < 0.05 (reported for transparency). No formula intervention is present.',ha='center',fontsize=9)
fig.tight_layout(w_pad=2);save(fig,'Figure2_Transcriptome_Wide_Results')
# Fig 3: all four fixed modules, all 10 samples. No selected gene-only panel.
fig,axes=plt.subplots(2,2,figsize=(11.8,8.7),sharey=True)
for i,(key,m) in enumerate(modules.items()):
    a=axes.flat[i];r=summary.set_index('module').loc[key];s=score[key]
    for gx,mask,color in [(0,ctrl,blue),(1,stc,orange)]:
        vals=s.to_numpy()[mask];xs=gx+np.linspace(-.12,.12,len(vals));a.scatter(xs,vals,c=color,s=48,zorder=3,edgecolors='white',lw=.5)
        a.plot([gx-.23,gx+.23],[vals.mean(),vals.mean()],color=color,lw=2.4)
        for j,(xx,yy) in enumerate(zip(xs,vals)):a.annotate(meta.index[mask][j],(xx,yy),xytext=(4,3),textcoords='offset points',fontsize=7)
    a.axhline(0,color='#cccccc',lw=.7);a.set_xticks([0,1],['Cancer-surgery\ncontrol','Severe STC']);a.set_xlim(-.45,1.55);a.set_ylabel('Mean standardized gene expression');a.set_title(chr(65+i)+'  '+labels[key],loc='left',fontweight='bold',fontsize=10)
    a.text(.03,.96,f'{int(r.n_scored_symbols)}/{int(r.n_source_symbols)} symbols scored\nDifference = {r.difference_STC_minus_control:+.3f}\nExact p = {r.exact_two_sided_p_252_permutations:.4f}; BH = {r.BH_FDR_four_modules:.4f}',transform=a.transAxes,va='top',fontsize=9,bbox={'facecolor':'white','edgecolor':'none','alpha':.85})
    a.set_ylim(-1.25,1.25)
fig.suptitle('Prespecified host-expression modules',y=1.01,fontweight='bold',fontsize=14)
fig.text(.5,-.028,'Symbols standardized across all 10 samples after blind VST. Exact two-sided p: all 252 balanced label assignments.\nExpression scores do not measure water flux, muscle force, neural activity, barrier function, or treatment efficacy.',ha='center',fontsize=9)
fig.tight_layout(h_pad=2.4,w_pad=2.5);save(fig,'Figure3_Prespecified_Modules')
print('MODULES AND FIGURES COMPLETE',flush=True)
