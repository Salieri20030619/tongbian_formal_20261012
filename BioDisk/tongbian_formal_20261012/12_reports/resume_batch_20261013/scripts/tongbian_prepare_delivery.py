from pathlib import Path
import subprocess,json,csv,hashlib,datetime,shutil,tarfile
P=Path('/mnt/BioDisk/wangmingyuan/tongbian_formal_20261012');O=P/'01_inventory/resume_20261013';I=P/'03_inputs/MTBLS13907_audit_v1'
receipts=[]
for n in ['meta_intensity_neg.xls','meta_intensity_pos.xls']:
 dest=I/n;url='https://ftp.ebi.ac.uk/pub/databases/metabolights/studies/public/MTBLS13907/FILES/DERIVED_FILES/'+n
 if not dest.exists():
  part=Path(str(dest)+'.partial')
  with (O/(n+'.download.log')).open('w') as log:r=subprocess.run(['curl','--noproxy','*','-fL','--connect-timeout','15','--max-time','300','--retry','2','-C','-','-o',str(part),url],stdout=log,stderr=log)
  if r.returncode==0:part.rename(dest)
 if dest.exists():
  b=dest.read_bytes();receipts.append({'file':n,'url':url,'bytes':len(b),'sha256':hashlib.sha256(b).hexdigest(),'magic_hex':b[:8].hex()});print('METABOLITE_FILE',n,'BYTES',len(b),'HEADER',b[:800],flush=True)
(O/'processed_metabolite_files_receipt.json').write_text(json.dumps(receipts,indent=2))
E=P/'03_inputs/references/Ensembl_human_mouse_one2one_20261013'
rr=subprocess.run(['curl','--noproxy','*','-fLsS','--max-time','30','https://www.ensembl.org/biomart/martservice?type=registry'],capture_output=True)
if rr.returncode==0:(E/'biomart_registry.xml').write_bytes(rr.stdout);print('BIOMART_REGISTRY',rr.stdout[:1800])
baseline={}
for n,expected in [('32_primary_unique_parent_structures.tsv',2059),('35_primary_human_targets.tsv',331),('36_GeneCards_primary_target_intersection.tsv',294)]:
 p=P/'04_modelA_reuse/tables'/n;rows=list(csv.DictReader(p.open(encoding='utf-8-sig'),delimiter='\t'));assert len(rows)==expected;baseline[n]={'rows':len(rows),'expected':expected,'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'passed':True}
(O/'ModelA_reuse_count_verification.json').write_text(json.dumps(baseline,indent=2));print('ModelA',baseline)
# Snapshot only this round, without touching older result batches.
D=P/'12_reports/resume_batch_20261013';D.mkdir(exist_ok=True)
for name in ['GSE319612_DESeq2_v1','GSE319612_functional_integrative_v1']:
 src=P/'07_animal_multiomics'/name
 if not (D/name).exists():shutil.copytree(src,D/name)
for name,src in [('input_audit',O),('GSE319612_inputs',P/'03_inputs/GSE319612_verified_v1'),('MTBLS13907_inputs',I),('Ensembl_mapping',E)]:
 if not (D/name).exists():shutil.copytree(src,D/name)
(D/'scripts').mkdir(exist_ok=True);(D/'jobs').mkdir(exist_ok=True)
for p in (P/'environment').glob('tongbian_*.py'):shutil.copy2(p,D/'scripts'/p.name)
for n in ['tongbian_animal_RNA.R','tongbian_animal_functional.R','tongbian_validate_animal.R']:shutil.copy2(P/'environment'/n,D/'scripts'/n)
for n in ['GSE319612_DESeq2_v1','GSE319612_functional_integrative_v1','GSE319612_validation_v1']:
 src=P/'logs/jobs'/n
 if src.exists() and not (D/'jobs'/n).exists():shutil.copytree(src,D/'jobs'/n)
with (D/'SERVER_SHA256.tsv').open('w') as f:
 f.write('file\tbytes\tsha256\n')
 for p in sorted(D.rglob('*')):
  if p.is_file() and p.name!='SERVER_SHA256.tsv':f.write(f'{p.relative_to(D).as_posix()}\t{p.stat().st_size}\t{hashlib.sha256(p.read_bytes()).hexdigest()}\n')
with tarfile.open(P/'12_reports/resume_batch_20261013.tar.gz','w:gz') as tar:tar.add(D,arcname=D.name)
print('DELIVERY_BUNDLE_READY',P/'12_reports/resume_batch_20261013.tar.gz')
