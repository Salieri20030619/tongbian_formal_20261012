from pathlib import Path,PurePosixPath
import json,hashlib,zipfile,shutil,datetime,sys,re
base=Path(__file__).resolve().parent
z=base/'GSE245885_first_batch_complete.zip'
receipt=json.loads((base/'GSE245885_first_batch_bundle_receipt.json').read_text())
assert z.stat().st_size==receipt['bytes']
assert hashlib.sha256(z.read_bytes()).hexdigest()==receipt['sha256'],'Whole bundle SHA256 mismatch'
stage=base/'first_batch_complete';stage.mkdir(exist_ok=True)
with zipfile.ZipFile(z) as f:
 names=f.namelist();assert len(names)==len(set(names)),'Duplicate ZIP members'
 for name in names:
  p=PurePosixPath(name);assert not p.is_absolute() and '..' not in p.parts
  target=stage.joinpath(*p.parts).resolve();assert target.is_relative_to(stage.resolve())
 f.extractall(stage)
manifest=json.loads((stage/'12_reports/first_batch_20261012/deliverable_manifest.json').read_text())
for d in manifest['files']:
 p=stage.joinpath(*PurePosixPath(d['relative_path']).parts)
 assert p.stat().st_size==d['bytes'],str(p)
 assert hashlib.sha256(p.read_bytes()).hexdigest()==d['sha256'],str(p)
report=stage/'12_reports/first_batch_20261012/首批正式分析_报告.html'
markup=report.read_text(encoding='utf-8')
for rel in re.findall(r'(?:href|src)="([^"]+)"',markup):
 if not rel.startswith(('https:','http:','#')):assert (report.parent/rel).resolve().exists(),rel
print('STAGED_CHECKSUM_AND_HTML_LINKS_PASS',len(manifest['files']),flush=True)
if '--deliver' not in sys.argv:raise SystemExit(0)
parent=Path(r'C:\Users\31937\Downloads\改良通便丸_正式分析_20261012').resolve()
assert parent.exists() and parent.name=='改良通便丸_正式分析_20261012'
dest=(parent/'07_first_batch_analysis').resolve();assert dest.is_relative_to(parent)
dest.mkdir(exist_ok=True)
files=[p for p in stage.rglob('*') if p.is_file()]
for p in files:
 target=dest/p.relative_to(stage);target.parent.mkdir(parents=True,exist_ok=True)
 if target.exists():assert hashlib.sha256(target.read_bytes()).hexdigest()==hashlib.sha256(p.read_bytes()).hexdigest(),'Existing delivery differs: '+str(target)
 else:shutil.copyfile(p,target)
for d in manifest['files']:
 p=dest.joinpath(*PurePosixPath(d['relative_path']).parts);assert p.stat().st_size==d['bytes'] and hashlib.sha256(p.read_bytes()).hexdigest()==d['sha256']
sync={'status':'WINDOWS_SYNC_VERIFIED','verified_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'server_root':manifest['root'],'Windows_destination':str(dest),'bundle_sha256':receipt['sha256'],'bundle_bytes':receipt['bytes'],'delivered_payload_files':len(files),'manifest_entries_verified':len(manifest['files']),'all_file_size_and_sha256_match':True,'report_html_local_links_valid':True,'old_preparation_and_legacy_results_modified':False,'statistical_jobs_restarted':False}
for p in [base/'windows_sync_receipt.json',dest/'windows_sync_receipt.json']:p.write_text(json.dumps(sync,ensure_ascii=False,indent=2),encoding='utf-8')
(dest/'00_打开结果.txt').write_text('首批正式分析已完成并逐文件SHA256核验。\n报告：12_reports\\first_batch_20261012\\首批正式分析_报告.html\n文本：12_reports\\first_batch_20261012\\首批正式分析_结果与限制.txt\n最终图：12_reports\\first_batch_figures\n全量结果：05_human_RNA\n脚本：environment；日志：logs\\jobs\n复现与输入哈希：12_reports\\first_batch_20261012\\REPRODUCIBILITY_README.txt\n当前执行状态以本目录00_protocol和windows_sync_receipt.json为准；根目录旧准备文件保留为历史。\n',encoding='utf-8-sig')
print(json.dumps(sync,ensure_ascii=False,indent=2),flush=True)
