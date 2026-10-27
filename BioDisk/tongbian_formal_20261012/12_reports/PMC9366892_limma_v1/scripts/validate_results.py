from pathlib import Path
import hashlib,tarfile,json,datetime
import pandas as pd,numpy as np
from scipy import stats
r=Path(__file__).resolve().parent;p=r/'PMC9366892_limma_v1_server_bundle.tar.gz'
assert p.stat().st_size==178780 and hashlib.sha256(p.read_bytes()).hexdigest()=='acacaed4446f35817c6b447559e64e2e047f9575acf86881deb55773e567ea5a'
dest=r/'server_results';dest.mkdir(exist_ok=True)
with tarfile.open(p) as t:
 for m in t.getmembers():
  q=(dest/m.name).resolve();assert q.is_relative_to(dest.resolve()) and m.isfile()
  q.parent.mkdir(parents=True,exist_ok=True);q.write_bytes(t.extractfile(m).read())
out=dest/'06_human_selected_metabolites/limma_v1'
manifest=json.loads((out/'audit/server_output_manifest.json').read_text())
for x in manifest:
 q=dest/x['path'];assert q.stat().st_size==x['bytes'] and hashlib.sha256(q.read_bytes()).hexdigest()==x['sha256']
result=pd.read_csv(out/'tables/all_fits_results.tsv',sep='\t');assert len(result)==166 and result.analysis_id.nunique()==8
Y=pd.read_csv(out/'audit/log2_peak_area_aligned.tsv',sep='\t').set_index('feature_id')
checks=[]
for analysis,tab in result.groupby('analysis_id',sort=False):
 design=pd.read_csv(out/f'audit/design_{analysis}.tsv',sep='\t').set_index('sample_id')
 yy=Y.loc[tab.feature_id,design.index].copy()
 if tab.scenario.iloc[0]=='source_B':yy.loc['M407T218','C6']=np.log2(146617.187)
 estimates=np.linalg.lstsq(design.to_numpy(),yy.to_numpy().T,rcond=None)[0][list(design.columns).index('groupSTC')]
 maxdiff=float(np.max(np.abs(estimates-tab.logFC.to_numpy())));assert maxdiff<1e-9
 pcalc=2*stats.t.sf(np.abs(tab.moderated_t),tab.df_total);assert np.allclose(pcalc,tab.P,atol=1e-15,rtol=1e-10)
 ix=np.argsort(tab.P.to_numpy());sortedP=tab.P.to_numpy()[ix];bh=np.minimum.accumulate((sortedP*len(tab)/np.arange(1,len(tab)+1))[::-1])[::-1];back=np.empty(len(tab));back[ix]=np.minimum(bh,1)
 assert np.allclose(back,tab.BH,atol=1e-15,rtol=1e-12)
 ci=stats.t.ppf(.975,tab.df_total)*tab.moderated_SE
 assert np.allclose(tab.logFC-ci,tab.CI95_low) and np.allclose(tab.logFC+ci,tab.CI95_high)
 checks.append({'analysis_id':analysis,'n_features':len(tab),'design_rank':int(np.linalg.matrix_rank(design.to_numpy())),'OLS_max_abs_logFC_difference':maxdiff,'moderated_P_recomputed':True,'BH_recomputed':True,'CI_recomputed':True})
ch=result[result.feature_id=='M407T218'];assert len(ch)==6 and ch.C6_source_conflict.all()
primary=result[result.analysis_id=='primary20_group'];adj=result[result.analysis_id=='primary20_age_sex_BMI']
primary=primary.set_index('feature_id');adj=adj.set_index('feature_id').loc[primary.index]
sign_agree=int((np.sign(primary.logFC)==np.sign(adj.logFC)).sum())
sig_agree=int(((primary.P<.05)==(adj.P<.05)).sum())
receipt={'status':'PASS','utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'server_manifest_hashes_verified':len(manifest),'fit_count':len(checks),'result_rows':len(result),'independent_checks':checks,'primary_adjusted_direction_agreement':f'{sign_agree}/20','primary_adjusted_P_status_agreement':f'{sig_agree}/20','formal_partial_correlation_run':False,'original_inputs_modified':False}
(r/'validation/independent_statistical_QA.json').write_text(json.dumps(receipt,indent=2))
print(json.dumps(receipt,indent=2))
print('PRIMARY20_EFFECTS')
print(primary[['Name','logFC','P','BH']].to_string())
