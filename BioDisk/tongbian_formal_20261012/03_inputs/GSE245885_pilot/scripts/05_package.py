#!/usr/bin/env python3
"""Build compact, reproducible deliverable, excluding environments and model pickles."""
import pathlib,hashlib,zipfile,json
p=pathlib.Path(__file__).resolve().parents[1]
allowed_audit={'GSE245885_design_audit_zh.md','module_independent_verification.json','sample_metadata_audited.tsv'}
files=[]
for f in sorted(p.rglob('*')):
 if not f.is_file():continue
 rel=f.relative_to(p)
 if any(x.startswith('.') for x in rel.parts):continue
 if f.suffix in {'.pkl','.zip'}:continue
 if str(rel).startswith('audit/') and f.name not in allowed_audit:continue
 if f.name=='DELIVERABLE_SHA256SUMS.txt':continue
 files.append(f)
checksums=''.join(hashlib.sha256(f.read_bytes()).hexdigest()+'  '+str(f.relative_to(p))+'\n' for f in files)
checksum_file=p/'DELIVERABLE_SHA256SUMS.txt';checksum_file.write_text(checksums);files.append(checksum_file)
out=p.parent/'tongbian_GSE245885_pilot_20261010.zip'
with zipfile.ZipFile(out,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=7) as z:
 for f in files:z.write(f,arcname=p.name+'/'+str(f.relative_to(p)))
with zipfile.ZipFile(out) as z:assert z.testzip() is None
result={'zip_path':str(out),'bytes':out.stat().st_size,'sha256':hashlib.sha256(out.read_bytes()).hexdigest(),'files':len(files)}
print(json.dumps(result,indent=2));(p.parent/'tongbian_GSE245885_pilot_20261010_zip_manifest.json').write_text(json.dumps(result,indent=2))
