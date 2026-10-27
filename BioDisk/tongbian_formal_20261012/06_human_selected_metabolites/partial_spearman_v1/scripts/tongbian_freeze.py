from pathlib import Path
import csv,json,hashlib,datetime,subprocess,os,shutil
P=Path('/mnt/BioDisk/wangmingyuan/tongbian_formal_20261012'); O=P/'06_human_selected_metabolites/partial_spearman_v1'
I=P/'03_inputs/PMC9366892_deidentified_v1/tongbian_human_PMC9366892_deidentified_inputs_20261012'
assert not (O/'analysis_config_FROZEN.json').exists(), 'Already frozen; do not repeat'
O.mkdir(exist_ok=True)
for d in ['audit','tables','figures','scripts','inputs']: (O/d).mkdir(exist_ok=True)
def rows(n):return list(csv.DictReader((I/n).open()))
joint=rows('joint_selected_60x17.csv');ann=rows('metabolites_21_annotation.csv');taxa=[x for x in joint[0] if x.startswith(('g__','f__'))]
aliases={'Arachidonic Acid (peroxide free)':'Arachidonic Acid'}
main=[a['ID'] for a in ann if aliases.get(a['Name'],a['Name']) in joint[0] and a['Name']!='Cholic acid'];ext=[a['ID'] for a in ann if a['Name']!='Cholic acid']
assert len(main)==8 and len(ext)==20 and len(taxa)==5
cfg={'frozen_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'input':str(I),'output':str(O),'taxa':taxa,'primary40':main,'extended100':ext,'cholic5':['M407T218'],'models':['unadjusted','group','group_age_sex_BMI'],'primary_model':'group','cholic_scenarios':['source_A','source_B','exclude_C6'],'BH_scope':'within each prespecified family x model x scenario; primary40 and extended100 are overlapping, not independent replication','method':'ppcor::pcor.test(method=spearman), two sided, df=n-k-2','missing':'complete cases per test; no imputation','alpha':0.05,'no_FC_threshold':True,'sex':'factor code1/code2, biological coding unknown','source_matrix':'full precision peak-area and 14-taxa matrices; joint_selected defines feature membership; verify agreement before tests'}
cfg['joint_name_aliases']=aliases
cfg['preflight_note']='Initial exact-name gate stopped before installation/tests: M303T39 annotation includes (peroxide free); explicit alias added, numeric/rank agreement required before analysis.'
(O/'analysis_config_FROZEN.json').write_text(json.dumps(cfg,indent=2))
with (O/'audit/input_sha256.tsv').open('w') as f:
 f.write('file\tbytes\tsha256\n')
 for p in sorted(I.iterdir()):
  if p.is_file():
   h=hashlib.sha256(p.read_bytes()).hexdigest();f.write(f'{p.name}\t{p.stat().st_size}\t{h}\n');shutil.copy2(p,O/'inputs'/p.name)
checks=[]
for line in (I/'SHA256SUMS.txt').read_text().splitlines():
 h,n=line.split(maxsplit=1);n=n.lstrip('*');actual=hashlib.sha256((I/n).read_bytes()).hexdigest();assert h==actual,(n,h,actual);checks.append(n)
(O/'audit/input_original_manifest_verified.json').write_text(json.dumps({'verified_files':checks,'all_passed':True},indent=2))
for n in ['tongbian_ppcor_install.R','tongbian_partial.R','tongbian_freeze.py','tongbian_run_job.py']:
 p=P/'environment'/n
 if p.exists():shutil.copy2(p,O/'scripts'/n)
print(json.dumps(cfg,indent=2))
