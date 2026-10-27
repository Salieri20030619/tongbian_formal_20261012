from pathlib import Path
import csv,json,math,hashlib,datetime
r=Path('/mnt/BioDisk/wangmingyuan/tongbian_formal_20261012')
out=r/'05_human_RNA/technical_validation_v1';out.mkdir(exist_ok=True)
def rows(p):
 with (r/p).open() as f:return list(csv.DictReader(f,delimiter='\t'))
def write(name,rs):
 with (out/name).open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rs[0]),delimiter='\t');w.writeheader();w.writerows(rs)
def bh(ps):
 idx=sorted(range(len(ps)),key=lambda i:ps[i]);res=[None]*len(ps);v=1.
 for k in range(len(idx)-1,-1,-1):
  v=min(v,ps[idx[k]]*len(ps)/(k+1));res[idx[k]]=v
 return res
def gmt(p):
 d={}
 for line in p.read_text().splitlines():
  a=line.split('\t');assert a[0] not in d;d[a[0]]=set(a[2:])
 return d
checks=[]
def check(name,value):
 assert value,name;checks.append({'check':name,'pass':True})
ref=r/'03_inputs/references/msigdbr_26.1.0_20261012'
rank=rows('05_human_RNA/host_functional_v1/all_testable_signed_Wald_ranking.tsv');rankset={z['gene_symbol'] for z in rank}
check('unique unfiltered-by-P Wald symbols 20365',len(rank)==len(rankset)==20365)
covsummary={};gseaall=[]
for fam,n in [('C2_CP',4115),('GO_BP',7538)]:
 sets=gmt(ref/f'{fam}.human.symbols.gmt');check(f'{fam} full catalog count',len(sets)==n)
 a=rows(f'05_human_RNA/fgsea_v1/{fam}_fgsea_ALL_results.tsv')
 cov=rows(f'05_human_RNA/fgsea_v1/{fam}_ALL_set_coverage.tsv')
 expected={z['pathway'] for z in cov if z['tested_by_size_rule']=='TRUE'}
 check(f'{fam} complete eligible result set',expected=={z['pathway'] for z in a})
 check(f'{fam} finite positive P and no exact zero',all(0<float(z['pval'])<=1 for z in a))
 check(f'{fam} finite NES and log2err',all(math.isfinite(float(z[k])) for z in a for k in ['NES','log2err']))
 check(f'{fam} matching size and leading edge',all(int(z['size'])==len(sets[z['pathway']]&rankset) and set(z['leadingEdge'].split(';'))<=sets[z['pathway']]&rankset for z in a))
 check(f'{fam} independent BH verification',max(abs(x-float(z['BH_within_collection'])) for x,z in zip(bh([float(z['pval']) for z in a]),a))<1e-12)
 covsummary[fam]={'catalog_sets':n,'matched_size_lt15':sum(int(z['matched_ranked_genes'])<15 for z in cov),'matched_size_gt500':sum(int(z['matched_ranked_genes'])>500 for z in cov),'tested':len(a),'P_lt005_up':sum(float(z['pval'])<.05 and float(z['NES'])>0 for z in a),'P_lt005_down':sum(float(z['pval'])<.05 and float(z['NES'])<0 for z in a),'max_log2err':max(float(z['log2err']) for z in a)}
 gseaall+=a
pooled=rows('05_human_RNA/fgsea_v1/fgsea_all_collections_ALL_results.tsv')
check('pooled BH independent verification',len(pooled)==6301 and max(abs(x-float(z['BH_across_collections'])) for x,z in zip(bh([float(z['pval']) for z in pooled]),pooled))<1e-12)
for f,n in [('ssGSEA_four_prespecified_modules',4),('VIPER_all_TFs',270)]:
 a=rows(f'05_human_RNA/host_functional_v1/{f}_limma_ALL_results.tsv')
 check(f'{f} feature and sample completeness',len(a)==n and all(z['n_Control']==z['n_STC']=='5' for z in a))
 check(f'{f} effect equals group mean difference',all(abs(float(z['mean_STC'])-float(z['mean_Control'])-float(z['delta_score']))<1e-10 for z in a))
 check(f'{f} finite inference',all(math.isfinite(float(z[k])) for z in a for k in ['delta_score','P_value','BH','CI95_low','CI95_high']))
 check(f'{f} independent BH verification',max(abs(x-float(z['BH'])) for x,z in zip(bh([float(z['P_value']) for z in a]),a))<1e-12)
tf=rows('05_human_RNA/host_functional_v1/VIPER_regulon_coverage_and_return_status.tsv')
check('TF return exactly matches prespecified minsize',all(z['eligible']==z['returned_by_VIPER'] for z in tf))
original=gmt(r/'03_inputs/references/c2.cp.v2026.1.Hs.symbols.gmt');chosen=gmt(ref/'C2_CP.human.symbols.gmt')
compare=[]
for k in sorted(set(original)|set(chosen)):
 a=original.get(k,set());b=chosen.get(k,set());compare.append({'pathway':k,'original_present':k in original,'package_present':k in chosen,'original_size':len(a),'package_size':len(b),'common_symbols':len(a&b),'original_only':';'.join(sorted(a-b)),'package_only':';'.join(sorted(b-a)),'identical_membership':a==b})
write('C2_CP_package_vs_archived_GMT_membership_audit.tsv',compare)
packagecompare={'same_set_names':set(original)==set(chosen),'identical_membership_sets':sum(z['identical_membership'] for z in compare),'different_membership_sets':sum(not z['identical_membership'] for z in compare),'decision':'Frozen existing msigdbr 26.1.0 human resource used consistently for both collections; local archived GMT retained only as audit comparison.'}
cache=Path(json.loads((out/'technical_validation_summary.json').read_text())['msigdbr_source_zip'])
cacheledger={'path':str(cache),'bytes':cache.stat().st_size,'sha256':hashlib.sha256(cache.read_bytes()).hexdigest(),'md5':hashlib.md5(cache.read_bytes()).hexdigest(),'official_package_expected_md5':'512ba99c6827141a9d471972b812d4ac','source':'https://zenodo.org/records/18968178/files/msigdb.2026.1.zip?download=1'}
summary={'status':'PASS','time_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'checks':checks,'collection_coverage':covsummary,'C2_resource_comparison':packagecompare,'source_cache_integrity':cacheledger,'pooled_BH_lt005':sum(float(z['BH_across_collections'])<.05 for z in pooled),'primary_results_modified':False}
(out/'output_integrity_summary.json').write_text(json.dumps(summary,indent=2))
write('output_integrity_checks.tsv',checks)
print(json.dumps(summary,indent=2))
print('TECHNICAL_VALIDATION',(out/'technical_validation_summary.json').read_text())
print('VIPER_FORMALS',(out/'VIPER_actual_formals.txt').read_text())
print('COOKS_LOGIC',(out/'DESeq2_results_Cooks_logic_excerpt.txt').read_text())
