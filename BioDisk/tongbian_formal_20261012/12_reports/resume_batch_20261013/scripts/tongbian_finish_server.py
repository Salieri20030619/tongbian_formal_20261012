from pathlib import Path
import json,hashlib,shutil,datetime,os
P=Path('/mnt/BioDisk/wangmingyuan/tongbian_formal_20261012');now=datetime.datetime.now(datetime.timezone.utc).isoformat();A=P/'06_human_selected_metabolites/partial_spearman_v1';D=P/'12_reports/resume_batch_20261013';J=P/'logs/jobs/PRJNA1232352_direct_download_v1'
for O in [A,D]:assert json.loads((O/'PDF_QA_PASS.json').read_text())['status']=='PDF_QA_PASS'
for n in ['ppcor_install_20261013','PMC9366892_partial_spearman_v1','GSE319612_DESeq2_v1','GSE319612_functional_integrative_v1','GSE319612_validation_v1']:assert json.loads((P/'logs/jobs'/n/'status.json').read_text())['exit_code']==0
for sub in ['audit/input_resource_hashes.json','audit/VALIDATION_COMPLETE.json','audit/result_integrity_checks.tsv']:
 src=P/'07_animal_multiomics/GSE319612_functional_integrative_v1'/sub;dest=D/'GSE319612_functional_integrative_v1'/sub;shutil.copy2(src,dest)
dst=D/'jobs/PRJNA1232352_direct_download_v1';dst.mkdir(exist_ok=True)
for n in ['config.json','status.json','runner_pid.txt','launcher.log']:
 if (J/n).exists():shutil.copy2(J/n,dst/n)
shutil.copy2(P/'environment/tongbian_background_download.py',D/'scripts/tongbian_background_download.py');shutil.copy2(P/'environment/tongbian_finish_server.py',D/'scripts/tongbian_finish_server.py')
state=json.loads((J/'status.json').read_text());pid=state['runner_pid'];alive=Path('/proc')/str(pid)/'cmdline';assert alive.exists() and b'tongbian_background_download.py' in alive.read_bytes()
probe={'time_utc':now,'runner_pid':pid,'runner_identity_verified':True,'status':state['status'],'completed_files':state['completed_files'],'partial_files':[]}
B=Path('/mnt/BioData/bio/wangmingyuan/tongbian_formal_20261012')
for f in (B/'raw_reads/PRJNA1232352').glob('*.partial'):probe['partial_files'].append({'file':str(f),'bytes':f.stat().st_size})
(dst/'delivery_runtime_verification.json').write_text(json.dumps(probe,indent=2))
for O,code in [(A,'PMC9366892_PARTIAL_SPEARMAN'),(D,'ANIMAL_HOST_AND_INPUT_AUDIT')]:
 marker={'status':'BATCH_COMPLETED_VERIFIED','batch':code,'completed_utc':now,'PDF_QA':'PASS','analysis_jobs_exit_zero':True,'project_complete':False,'ongoing_raw_download_excluded_from_completion':True,'windows_sync':'verify using FINAL_WINDOWS_SYNC_VERIFIED.json'}
 (O/'DELIVERY_COMPLETE.json').write_text(json.dumps(marker,indent=2))
 with (O/'FINAL_SERVER_SHA256.tsv').open('w') as f:
  f.write('file\tbytes\tsha256\n')
  for p in sorted(O.rglob('*')):
   if p.is_file() and p.name not in ['FINAL_SERVER_SHA256.tsv','FINAL_WINDOWS_SYNC_VERIFIED.json']:f.write(f'{p.relative_to(O).as_posix()}\t{p.stat().st_size}\t{hashlib.sha256(p.read_bytes()).hexdigest()}\n')
cp=P/'00_protocol/project_config.json';c=json.loads(cp.read_text());c['updated_utc']=now;c['network']['large_download_allowed']=True;c['network']['current_download_route']='direct by explicit user instruction; slow rate accepted; no proxy tunnel created';c['storage']['large_metagenome_download_allowed']=True;c['storage']['large_data_mount_available_bytes']=os.statvfs(B).f_bavail*os.statvfs(B).f_frsize;c['storage']['current_available_bytes']=os.statvfs(P).f_bavail*os.statvfs(P).f_frsize;c['storage']['capacity_updated_utc']=now;c['execution']['project_complete']=False;c['execution']['ongoing_download_job']=str(J);cp.write_text(json.dumps(c,ensure_ascii=False,indent=2))
sp=P/'00_protocol/resume_state_20261013.json';s=json.loads(sp.read_text());s['updated_utc']=now;s['project_complete']=False;s['delivered_batches']=[str(A),str(D)];s['raw_download']=probe;s['raw_download']['user_preference']='background direct overnight; proxy tomorrow; no new tunnel this turn';sp.write_text(json.dumps(s,ensure_ascii=False,indent=2))
hp=P/'HANDOFF.md';archive=P/'00_protocol/HANDOFF_history_through_20261013.md'
if not archive.exists():shutil.copy2(hp,archive)
text=f'''# 改良通便丸公开多组学研究：当前交接

更新时间：{now}。整个项目未完成；本轮两批计算与报告完成，核心宏基因组后台直连下载中。

## 明日首先执行
1. 使用SSH别名sgs10，实际hostname osv01、用户wangmingyuan；沿用现有认证。
2. 先读 logs/jobs/PRJNA1232352_direct_download_v1/status.json、runner_pid.txt 和 config.json；当前主PID {pid}，4路直连。核对进程身份和文件大小，不根据进程消失判定成功。
3. 用户最新指令：今晚后台直连，速度慢也继续；明日再连代理。本轮未启动代理隧道、未改网络或认证配置。此前FASTQ_download_preflight.json中的速度门槛已被此明确指令覆盖。
4. 如明日切换代理，先确认本项目原下载进程及其子进程停止写入，再在明确验证的任务代理路线续传；保留所有.partial；严禁并发双写和重复启动。逐文件大小、ENA MD5、gzip校验完成后才改名。下载器使用独立进程组和run.lock。全部完成才生成DOWNLOAD_COMPLETE.json。

## 本轮真实完成
- 项目agent.md、交接及结构化状态建立；旧状态归档保留。Model A表复核2059/331/294；7个历史任务退出0；冻结MSigDB与DoRothEA哈希复核通过。没有重跑既有人RNA、limma或Model A。
- ppcor 1.1从CRAN装入独立项目库并测试成功。60例按sample_id配对；主40对、扩展100对和胆酸5对预先冻结。465行实际计算/独立公式/行序核验通过。主group调整15对名义显著、12对BH显著；协变量敏感性14/11；扩展主模型19/5。胆酸冲突原值保留。
- GSE319612：NC8/MC8/RG10组别与count列核实，DESeq2两对比完成。各16961有效检验；疾病1247/9、干预1494/4（P/BH）。98个两对比名义显著反向基因；不称恢复正常。
- 动物真实ssGSEA、VIPER+limma、fgsea及人鼠Fisher完成，23项结果核验通过。ssGSEA疾病0/4、干预2/4名义显著；VIPER疾病16/269、干预30/269；这些评分均无BH显著且无双名义反向候选。
- fgsea C2:CP两对比各2264项，疾病P/BH=569/209、干预864/581；GO BP各3628项，疾病807/282、干预1102/537。沿用冻结2026.1.Hs并按Ensembl Genes116高可信一对一映射，不是新下载的鼠MSigDB。
- 人鼠背景12367；上调重叠190，OR1.165、P0.0608、BH0.1216；下调130，OR0.788、P0.989。当前数据不支持整体同向富集，此主要阴性结果保留。
- 两份中文PDF已逐页Poppler检查；完整表、输入、版本、脚本、日志、图和哈希在新增10/11批次，旧07/08/09批次未覆盖。

## 正在运行
- 核心PRJNA1232352：60样本、120文件、293286068526字节，DCTT/NCTT/H各20。
- 数据目录：/mnt/BioData/bio/wangmingyuan/tongbian_formal_20261012/raw_reads/PRJNA1232352。
- 运行记录：{J}；本次交付时进程身份已核实，完成文件{state['completed_files']}/120。实时进度以服务器状态和实际.partial大小为准，Windows仅为交付时快照。
- 下载器原生等待curl，不需要反复agent级轮询；不会自动启动MetaPhlAn/HUMAnN或组装。

## 未解决事项与后续
- MTBLS13907：原ISA30记录中22条组标签冲突、同组内另有编号错位；26动物+3QC+1blank。完整正3193/负2115峰强度已取到，均为处理后强度；不得重复归一化/填补。组别未解，不做差异/反向；跨组学无可靠个体交叉表。
- PRJNA1348869为26扩增子+1菌株WGS，不能整体当小鼠宏基因组；不按组内编号与RNA/代谢配对。
- PRJNA1048321和PRJNA1031121独立性未确认：共同机构及伦理2021KL-023，文件MD5无完全相同不等于人群独立；暂不作为独立验证。
- MetaPhlAn、HUMAnN、MEGAHIT、Prodigal、DIAMOND、dbCAN在已检查PATH及现有conda环境bin中未找到；正常依赖安装已授权，后续用隔离环境和可信源。大库安装前重新查容量和路线。
- 原始reads校验后继续菌群功能与转化层；模型B需实验依据酶/底物/反应方向。Model C为证据分层，不要求各层均显著。尚未启动GraphBAN/对接/MD，不将公开疾病或Rothia干预解释为本方药效。

## 交付位置
- Windows根：C:/Users/31937/Downloads/改良通便丸_正式分析_20261012。
- 10_partial_spearman_20261013/partial_spearman_v1：偏相关完整结果与7页中文PDF。
- 11_animal_and_next_modules_20261013/resume_batch_20261013：动物宿主结果、后续审计与8页中文PDF。
- 服务器对应：{A}；{D}。
- 每批FINAL_SERVER_SHA256.tsv与FINAL_WINDOWS_SYNC_VERIFIED.json核验同步；长下载完成状态独立，不混入本批完成标记。
'''
hp.write_text(text,encoding='utf-8')
with (P/'agent.md').open('a') as f:f.write('\n## 2026-10-13 最新后台下载安排\n用户允许本项目先后台直连下载，明日再连代理，低速不是停止理由。当前下载任务PRJNA1232352_direct_download_v1；先核对PID/状态和.partial后接续，不并发双写，不覆盖已校验文件。计算、报告和下载分别标记完成；详细现状以HANDOFF.md为准。\n')
print(json.dumps(probe,indent=2));print('SERVER_FINAL_PACKAGES_AND_HANDOFF_READY')
