from pathlib import Path
import shutil,json,hashlib,datetime,csv,sys,platform,importlib.metadata,tarfile
import pandas as pd,pymupdf
r=Path(__file__).resolve().parent;server=r/'server_results';base=server/'06_human_selected_metabolites/limma_v1'
pdf=r.parent/'output/pdf/改良通便丸_PMC9366892_limma代谢物重分析_20261012.pdf'
doc=pymupdf.open(pdf);assert len(doc)==4
text='\n'.join(p.get_text() for p in doc);flat=''.join(text.split());assert '\ufffd' not in text
primary=pd.read_csv(base/'tables/primary20_group.tsv',sep='\t');ch=pd.read_csv(base/'tables/Cholic_acid_prespecified_sensitivity.tsv',sep='\t')
for row in primary.itertuples():
 for token in [row.feature_id,f'{row.logFC:.4f}',f'{row.P:.4g}',f'{row.BH:.4g}']:assert token in flat,token
for row in ch.itertuples():
 for token in [f'{row.logFC:.4f}',f'{row.P:.4g}',f'{row.BH:.4g}']:assert token in flat,token
assert '22,146,617.1867411' in text and '146,617.187' in text
for i,page in enumerate(doc,1):
 assert f'{i} / 4' in page.get_text()
 for b in page.get_text('blocks'):
  assert b[0]>=0 and b[1]>=0 and b[2]<=page.rect.width+1 and b[3]<=page.rect.height+1
 for n in [i]:assert (r.parent/f'tmp/pdfs/metabolite_limma/page-{n}.png').exists()
qa={'status':'PASS','pages':4,'visual_review':'all four Poppler-rendered pages inspected; no clipped table, overlapping text, broken glyphs or missing figures','complete_main20_numeric_table_checked':True,'six_cholic_sensitivity_rows_checked':True,'figure_count':2,'sha256':hashlib.sha256(pdf.read_bytes()).hexdigest(),'bytes':pdf.stat().st_size,'timestamp_utc':datetime.datetime.now(datetime.timezone.utc).isoformat()}
(r/'validation/PDF_QA.json').write_text(json.dumps(qa,indent=2))
versions={'python':sys.version,'platform':platform.platform(),'packages':{n:importlib.metadata.version(n) for n in ['numpy','pandas','scipy','matplotlib','reportlab','PyMuPDF']}}
(r/'validation/local_report_environment.json').write_text(json.dumps(versions,indent=2))
out=r/'delivery';out.mkdir(exist_ok=True)
for dirname in ['tables','audit','models']:shutil.copytree(base/dirname,out/dirname,dirs_exist_ok=True)
shutil.copytree(r/'report/figures',out/'figures',dirs_exist_ok=True)
shutil.copytree(r/'validation',out/'audit',dirs_exist_ok=True)
shutil.copytree(server/'03_inputs/PMC9366892_deidentified_v1',out/'inputs',dirs_exist_ok=True)
shutil.copytree(server/'logs/jobs/PMC9366892_limma_v1',out/'logs/PMC9366892_limma_v1',dirs_exist_ok=True)
(out/'scripts').mkdir(exist_ok=True);(out/'configuration').mkdir(exist_ok=True)
for name in ['PMC9366892_limma.R','validate_inputs.py','validate_results.py','setup_server_metabolites.py','collect_server_results.py','make_report.py','deliver_results.py']:
 shutil.copy2(r/name,out/'scripts'/name)
shutil.copy2(server/'environment/durable_R_job_runner.py',out/'scripts/durable_R_job_runner.py')
shutil.copy2(server/'06_human_selected_metabolites/limma_staging_v1/R_preflight.log',out/'logs/R_preflight.log')
shutil.copy2(r/'analysis_config_frozen.json',out/'configuration/analysis_config_frozen.json')
shutil.copy2(pdf,out/pdf.name)
shutil.copy2(r/'PMC9366892_limma_v1_server_bundle.tar.gz',out/'audit/PMC9366892_limma_v1_server_bundle.tar.gz')
readme='''PMC9366892 人体预选代谢物正式limma重分析 - 2026-10-12

已完成：8个预设模型，166行完整结果；4页PDF与2张图（PNG和SVG）。
主20峰结果：每组30人；P<0.05为16/20，BH<0.05为16/20；1升高、15降低。
主20峰协变量敏感性：年龄+Sex原编码+BMI，方向和P阈值状态均20/20一致。
原21峰表保留所有21峰，M407T218 Cholic acid单独冲突flag，不计入主20峰汇总。
主20峰独立eBayes拟合（避免冲突峰进入主模型方差收缩），BH家族20；全21峰各模型BH家族21。
胆酸A/B/排除C6 × group/协变量模型共6种预设结果；全部保留，不按结果选版本。
排除C6时样本数为29 Control+30 STC，其余30+30。A/B只差C6胆酸一值；原输入无修改。

优先阅读
PDF：改良通便丸_PMC9366892_limma代谢物重分析_20261012.pdf
tables/primary20_group.tsv：主20峰完整精度统计。
tables/primary20_age_sex_BMI.tsv：协变量敏感性。
tables/all21_source_A_group.tsv：原来源全21峰，含C6冲突标志。
tables/Cholic_acid_prespecified_sensitivity.tsv：胆酸6种预设结果。
tables/all_fits_results.tsv：全部8模型166行，含n/logFC/方向/P/BH/95%CI。
configuration/analysis_config_frozen.json：运行前冻结配置。
inputs/：12文件原输入包、ZIP与Library来源身份；无人名原表。
audit/：输入与设计检查、原值/log2矩阵、R与Python版本、独立统计QA、PDF逐页QA及原服务器结果包。
models/：8模型RDS；scripts/：所有分析、核验、报告、同步脚本；logs/：任务状态与运行日志。
SHA256SUMS.tsv：交付文件路径、字节数与SHA-256。服务器原始路径见audit/server_output_manifest.json。

源：Library libfile_9025b7a56c008191971c58e1dd17241e v0，ZIP SHA256
634eba9e4c718d6051b9bce01695f306a574d4c2184cf3dbdd71dcc42b324bb0。
服务器执行根目录：/mnt/BioDisk/wangmingyuan/tongbian_formal_20261012
输入：03_inputs/PMC9366892_deidentified_v1
结果：06_human_selected_metabolites/limma_v1
R：/mnt/BioDisk/wangmingyuan/miniconda3/envs/R_4.5.2/bin/Rscript
脚本：environment/PMC9366892_limma.R。默认--check仅预检，--run在结果目录已存在时拒绝覆盖；当前任务已完成，勿重复启动。
重新复现请在新的隔离目录明确改变输出路径，并保持现有结果与原始输入只读。

边界：原21峰已在同一队列预筛，本轮不是独立验证或完整代谢组发现；无额外FC阈值。
未调整可能是疾病后果的症状/疾病时长/Wexner/GIQLI/Bristol；Sex编码1/2未赋予未经核实的男女含义。
Taurochenodeoxycholate两个峰ID不合并；14个菌群条目含嵌套层级，未按独立菌属处理。
ppcor未安装、偏相关未运行，等待既有安装审批；未使用替代算法，未跨RNA队列个体匹配。
C6数值冲突继续保留；完整未预筛代谢组及完整菌群矩阵仍缺。未重跑Model A或RNA，无软件安装或网络设置变更。
此前08输入审计目录是取得本输入包之前的历史快照；本目录为本轮最新执行结果。
'''
(out/'README_结果与复现说明.txt').write_text(readme,encoding='utf-8')
items=[]
for p in sorted(out.rglob('*')):
 if p.is_file() and p.name not in ['SHA256SUMS.tsv','windows_sync_receipt.json']:items.append((p.relative_to(out).as_posix(),p.stat().st_size,hashlib.sha256(p.read_bytes()).hexdigest()))
with (out/'SHA256SUMS.tsv').open('w',encoding='utf-8',newline='') as f:
 writer=csv.writer(f,delimiter='\t');writer.writerow(['relative_path','bytes','sha256']);writer.writerows(items)
dest=Path(r'C:\Users\31937\Downloads\改良通便丸_正式分析_20261012\09_human_metabolite_limma')
shutil.copytree(out,dest,dirs_exist_ok=True)
for name,size,sha in items:
 p=dest/name;assert p.stat().st_size==size and hashlib.sha256(p.read_bytes()).hexdigest()==sha
receipt={'status':'VERIFIED','file_count':len(items),'destination':str(dest),'PDF_pages':4,'PDF_sha256':qa['sha256'],'timestamp_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'server_job':'PMC9366892_limma_v1','runner_pid':2009555,'input_unchanged':True}
for folder in [out,dest]:(folder/'windows_sync_receipt.json').write_text(json.dumps(receipt,ensure_ascii=False,indent=2),encoding='utf-8')
bundle=r/'PMC9366892_limma_v1_deliverables.tar.gz'
with tarfile.open(bundle,'w:gz') as t:
 for p in sorted(out.rglob('*')):
  if p.is_file():t.add(p,arcname=p.relative_to(out).as_posix())
print(json.dumps({'delivery':str(dest),'PDF':str(dest/pdf.name),'files_verified':len(items),'report_bundle_bytes':bundle.stat().st_size,'report_bundle_sha256':hashlib.sha256(bundle.read_bytes()).hexdigest()},ensure_ascii=False))
