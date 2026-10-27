from pathlib import Path
import csv,json,math,hashlib,datetime,html,shutil,zipfile,platform
r=Path('/mnt/BioDisk/wangmingyuan/tongbian_formal_20261012')
out=r/'12_reports/first_batch_20261012';out.mkdir(parents=True,exist_ok=True)
def readj(p):return json.loads((r/p).read_text())
def rows(p):
 with (r/p).open(encoding='utf-8-sig') as f:return list(csv.DictReader(f,delimiter='\t'))
def wt(path,rs):
 with path.open('w',encoding='utf-8',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rs[0]),delimiter='\t');w.writeheader();w.writerows(rs)
def fmt(v):return f'{float(v):.5g}'
jobs={n:readj('logs/jobs/'+n+'/status.json') for n in ['GSE245885_R_DESeq2_bridge_v1','MSigDB_package_resource_freeze','GSE245885_host_functional_scores_v1','GSE245885_sample_QC_v1','GSE245885_fgsea_v1','GSE245885_technical_validation_v1']}
assert all(z['status']=='COMPLETED' and z['exit_code']==0 for z in jobs.values())
validation=readj('05_human_RNA/technical_validation_v1/technical_validation_summary.json')
integrity=readj('05_human_RNA/technical_validation_v1/output_integrity_summary.json');assert integrity['status']=='PASS'
fg=readj('05_human_RNA/fgsea_v1/fgsea_summary.json');qc=readj('05_human_RNA/sample_QC_v1/sample_QC_summary.json')
ss=rows('05_human_RNA/host_functional_v1/ssGSEA_four_prespecified_modules_limma_ALL_results.tsv')
tf=sorted(rows('05_human_RNA/host_functional_v1/VIPER_all_TFs_limma_ALL_results.tsv'),key=lambda z:float(z['P_value']))
sscov={z['module']:z for z in rows('05_human_RNA/host_functional_v1/ssGSEA_module_coverage.tsv')}
tfcov={z['TF']:z for z in rows('05_human_RNA/host_functional_v1/VIPER_regulon_coverage_and_return_status.tsv')}
nice={'WATER_SALT_TRANSPORT':'水盐转运','SMOOTH_MUSCLE_CONTRACTION':'平滑肌收缩','ENTERIC_NERVOUS_SYSTEM_DEVELOPMENT':'肠神经系统发育','BICELLULAR_TIGHT_JUNCTION_ASSEMBLY':'双细胞紧密连接组装'}
ssshort=[]
for z in ss:
 c=sscov[z['feature']];ssshort.append({'module':z['feature'],'中文':nice[z['feature']],'matched_genes':c['matched_symbols'],'total_genes':c['total_members'],'delta_STC_minus_Control':z['delta_score'],'P':z['P_value'],'BH':z['BH'],'n_STC':5,'n_Control':5})
wt(out/'four_prespecified_modules_summary.tsv',ssshort)
tfshort=[{**{k:z[k] for k in ['feature','delta_score','P_value','BH','CI95_low','CI95_high','n_STC','n_Control']},'matched_targets':tfcov[z['feature']]['matched_targets'],'total_targets':tfcov[z['feature']]['total_targets']} for z in tf]
wt(out/'all_270_TFs_ranked_by_P.tsv',tfshort)
allfg=rows('05_human_RNA/fgsea_v1/fgsea_all_collections_ALL_results.tsv')
cols=['family','pathway','NES','pval','BH_within_collection','BH_across_collections','size','log2err','P_display']
top=[]
for fam in ['C2_CP','GO_BP']:
 a=sorted([z for z in allfg if z['family']==fam],key=lambda z:float(z['pval']))
 top += [{k:z[k] for k in cols} for z in a[:20]]
wt(out/'pathways_top20_each_collection_DISPLAY_ONLY.tsv',top)
summary_rows=[{'module':'R DESeq2','eligible_or_scored':20373,'P_lt_005':12484,'BH_lt_005':11756,'note':'20,479 prefiltered; 106 Cook-filtered; BH across finite gene tests'},
 {'module':'C2:CP fgsea','eligible_or_scored':2452,'P_lt_005':577,'BH_lt_005':158,'note':'complete 4,115-set catalog; size 15..500 after mapping'},
 {'module':'GO BP fgsea','eligible_or_scored':3849,'P_lt_005':689,'BH_lt_005':190,'note':'complete 7,538-set catalog; size 15..500 after mapping'},
 {'module':'four ssGSEA module limma','eligible_or_scored':4,'P_lt_005':3,'BH_lt_005':3,'note':'separate BH family of four prespecified modules'},
 {'module':'VIPER TF limma','eligible_or_scored':270,'P_lt_005':153,'BH_lt_005':139,'note':'separate BH family of all 270 scored TFs; 1 insufficient regulon'}]
wt(out/'completed_modules_summary.tsv',summary_rows)
jobrows=[{'job':k,'status':v['status'],'exit_code':v['exit_code'],'runner_pid_completed':v['runner_pid'],'completed_utc':v['updated_utc'],'log':str(r/'logs/jobs'/k/'run.log'),'script_sha256':v['script_sha256']} for k,v in jobs.items()]
wt(out/'job_completion_receipts.tsv',jobrows)
warnings=[{'source':'host_functional_v1 run.log','warning':'R stack imbalance in :: and <-','impact':'Original triggering native component not established. Base-R ssGSEA and direct VIPER cross-checks agree within floating-point tolerance. No missing scores.','disposition':'retain original logs; disclose software warning, no package installation or primary result changes'},
 {'source':'PNG drawing jobs','warning':'Fontconfig: using without calling FcInit()','impact':'PNG files produced; final plots visually checked after correcting label overlap/clipping.','disposition':'render-only warning; numeric results unchanged'},
 {'source':'fgsea C2:CP and GO BP','warning':'none from fgsea algorithm','impact':'all 6,301 returned P values finite and positive; no NA rescue needed','disposition':'eps=0; min P and log2err retained; no exact-zero claim'}]
wt(out/'warnings_and_dispositions.tsv',warnings)
sources=[{'source':'GEO GSE245885','url':'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE245885','role':'cohort and count input; local raw and SOFT frozen in prior preparation'},
 {'source':'MSigDB human collection catalog','url':'https://www.gsea-msigdb.org/gsea/msigdb/collections.jsp','role':'2026.1.Hs C2:CP 4115 and GO BP 7538 complete catalog counts'},
 {'source':'msigdbr official package documentation','url':'https://igordot.github.io/msigdbr/articles/msigdbr-intro.html','role':'existing package data access and species/ID semantics'},
 {'source':'msigdbr 26.1.0 embedded public release URL','url':integrity['source_cache_integrity']['source'],'role':'cache MD5 matches package expectation; SHA256 recorded; no authentication bypass'},
 {'source':'fgsea official vignette','url':'https://bioconductor.org/packages/release/bioc/vignettes/fgsea/inst/doc/fgsea-tutorial.html','role':'fgseaMultilevel, eps and log2err'},
 {'source':'GSVA official vignette','url':'https://bioconductor.org/packages/release/bioc/vignettes/GSVA/inst/doc/GSVA.html','role':'ssgseaParam and gsva interface'},
 {'source':'decoupleR official run_viper reference','url':'https://saezlab.github.io/decoupleR/reference/run_viper.html','role':'VIPER scoring interface; actual installed source/formals also archived'}]
wt(out/'methods_and_resource_sources.tsv',sources)
now=datetime.datetime.now(datetime.timezone.utc).isoformat()
lines=['改良通便丸项目：GSE245885 首批正式宿主分析实绩（2026-10-12）',f'汇总时间（UTC）：{now}',
 '状态：R DESeq2、fgseaMultilevel、真实 ssGSEA、decoupleR VIPER、limma 评分比较及样本质控全部完成。6 个持久任务退出码均为 0。',
 '本报告是实际输出汇总；Model A 仍直接复用，未重跑；未开展宏基因组、动物组学或分子动力学。',
 '', '一、分析口径与输入',
 'GSE245885：5 STC 与 5 癌症手术正常肠段对照，独立组，物种 Homo sapiens；没有本方干预。对照不是健康志愿者。',
 '从 10 个 *.genes.FPKM.txt.gz 的 *.count 整数列建模，未用 FPKM 拟合；35,671 行。继承 count>=10 至少5样本的低表达预过滤后20,479行；不是 FC 阈值。',
 '设计 ~condition；Control 为基线；STC-Control 为方向。年龄、性别、肠段、批次逐样本协变量缺失，未捏造调整；不删样本。名义 P<0.05 为探索性主口径，保留效应、N、置信区间及 BH，不另加 FC 门槛。',
 '', '二、已完成结果规模（BH 按各预设检验族）']
for z in summary_rows:lines.append(f"{z['module']}：有效/评分 {z['eligible_or_scored']}；P<0.05 {z['P_lt_005']}；BH<0.05 {z['BH_lt_005']}。")
lines += ['R DESeq2：名义显著上调6,386、下调6,098；106行主检验P被Cook规则过滤，无原始P=0。',
 '与旧 PyDESeq2 0.5.4 桥接：共同可检验20,373，效应方向100%一致，log2FC差值绝对值中位数2.20003e-06、最大0.002309。共同名义显著12,475；R独有9，Python独有1。原结果和边界变动逐行审计均保留。',
 '', '三、四个预设宿主模块：GSVA真实ssGSEA，再用limma比较（每组n=5）']
for z in ssshort:lines.append(f"{z['中文']}：覆盖{z['matched_genes']}/{z['total_genes']}基因；Δscore={fmt(z['delta_STC_minus_Control'])}；P={fmt(z['P'])}；BH={fmt(z['BH'])}。")
lines += ['平滑肌模块分数升高不等于收缩功能增强；肠神经发育和水盐转运分数下降是本队列转录组评分现象。紧密连接模块未达P<0.05。',
 '这些模块不能被改写为通便丸疗效或逆转证据；组织组成、癌症对照背景和未观测协变量可能影响分数。',
 '', '四、全谱通路与TF：所有有效检验均保留',
 'fgsea：20,373个可检验Ensembl行按基因符号去重后20,365个有符号Wald分数；没有P值预筛，也无Wald并列值。重复符号选合并样本baseMean最高行，平手按Ensembl ID；选择不依赖效应或P。',
 '完整C2:CP 4,115集：匹配后<15基因1,644集、>500基因19集，实际检验2,452集。完整GO BP 7,538集：<15基因3,468集、>500基因221集，实际检验3,849集。未检验属于大小规则，不是资源缺失。',
 'C2:CP名义显著正向440/负向137；GO BP正向525/负向164。两库合并6,301检验另保留pooled BH，BH<0.05共349。通路高度重叠，这些计数不等于独立机制数量。',
 'fgsea eps=0、sampleSize=101、nPermSimple=1000、minSize=15、maxSize=500、nproc=1、scoreType=std；随机种子20261013/20261014。全部P有限且>0，无算法警告，无NA重试。最小P：C2:CP 7.51307e-14；GO BP 1.40807e-15；log2err最大分别0.9545/1.0175，P为数值估计不是精确概率。',
 '通路示例：GOBP_RNA_SPLICING NES=1.9221，P=1.4081e-15，BH=5.4197e-12；REACTOME_INITIAL_TRIGGERING_OF_COMPLEMENT NES=-2.4392，P=6.2132e-09，BH=2.0096e-06。展示Top20按P固定选取；全量表未筛删。',
 'VIPER：冻结人DoRothEA1.22.0 A-C，271TF/13,223边；270TF通过至少5个靶基因覆盖。NR5A1只匹配4/15靶点，未评分。所有270TF×10样本评分完整。',
 'VIPER输入为20,471符号的VST矩阵；重复符号取各Ensembl VST行算术均值，无常量符号移除。method=scale（基因在10样本内标准化）、minsize=5、pleiotropy=TRUE、eset.filter=FALSE、bootstraps=0、cores=1、nes=TRUE。mor按TF最大绝对值归一化（本资源abs(mor)=1），likelihood全为1；pleiotropyArgs使用归档的包默认值。',
 'VIPER评分是根据靶基因推断的TF活性，不是TF表达量。limma分别对4个ssGSEA模块与270TF用lmFit/eBayes，trend=FALSE、robust=FALSE，效应为Δscore。',
 '按P排序的TF前5：']
for z in tfshort[:5]:lines.append(f"{z['feature']}：Δscore={fmt(z['delta_score'])}；P={fmt(z['P_value'])}；BH={fmt(z['BH'])}；靶点覆盖{z['matched_targets']}/{z['total_targets']}。")
lines += ['', '五、高差异比例的质控依据与解释边界',
 '原始基因计数库总量23,712,981–28,432,527，最大/最小1.199；这些是供应表的gene-assigned counts，不等于原始FASTQ读数。size factors 0.896448–1.142308，比例1.2743。',
 '最易变500基因PCA：PC1=82.5945%、PC2=4.9605%；全部20,479基因PCA：PC1=65.9713%、PC2=9.7436%；10个样本在PC1按组分离。组内Pearson范围0.9442–0.9906，组间0.8678–0.9338。',
 '对照NP2组内平均相关最低（0.95193），到自身组质心距离最大（74.143），Cook>8.6491的基因98个。已标记，但不删样本、未改变设计。NP3亦有较高组内距离（69.539），完整逐样本表同时保留。',
 '检测到非零计数的基因数：对照21,973–26,291，STC28,811–29,406；与组别相关的检测/组成差异值得核查，不可只凭库大小相近排除技术或组织混杂。',
 '119个基因任一样本Cook值超过阈值，实际主P过滤106个：另13个依DESeq2二水平因子内置规则保留（至少3个样本计数高于Cook最大点的计数）；源码逻辑与逐基因审计已保存。未关闭主分析的Cook规则，也未替换原结果。',
 '10样本计数向量无完全重复；设计秩2/2、残差自由度8。PCA/相关提示广泛组间结构，但不能辨别疾病作用与未记录临床、肠段、批次或细胞组成。较多P<0.05不构成本方机制或疗效证据。',
 '', '六、资源、覆盖与可重复性',
 'R4.5.2；DESeq2 1.50.2；fgsea1.36.2；GSVA2.4.9；decoupleR2.16.0；dorothea1.22.0；viper1.44.0；limma3.66.0；msigdbr26.1.0。实际sessionInfo逐模块保留。',
 'MSigDB2026.1.Hs来自已安装msigdbr的公开缓存，zip的MD5与包固定值512ba99c6827141a9d471972b812d4ac一致；源URL、zip SHA256、GMT、注释长表RDS和逐文件哈希齐全。未绕过登录，也未安装包。',
 '人类HS源/目标注释逐行核验；采用供应基因符号精确匹配，未猜测别名或跨物种同源转换。C2:CP资源13,957唯一符号，Wald匹配11,162；GO BP18,215符号，Wald匹配13,641。未匹配符号与每集覆盖均导出。',
 '包内C2:CP与原Windows归档GMT的4,115集名称和成员集合全部一致；字节哈希差异来自输出文本表示，未更换成员。人DoRothEA和4个预设模块成员保持冻结。',
 'R DESeq2原模型RDS、VST、全基因结果、桥接审计、fgsea全检验/leading edge/覆盖/误差、样本评分、limma拟合RDS、输入RDS、参数、随机种子、脚本哈希、日志和完成状态全部保留。',
 '', '七、警告与核验',
 f"评分原日志含R stack imbalance in :: / <-；具体触发组件尚未定位。独立base-R加权CDF ssGSEA与正式分数最大绝对差{validation['ssGSEA_independent_max_abs_difference']:.3g}；直接VIPER API与decoupleR正式分数最大差{validation['VIPER_direct_API_max_abs_difference']:.3g}。两项均通过，未更改原结果。",
 '图形Fontconfig提示保留；只重绘PCA重叠标签、热图裁切标题并补色标，统计矩阵未变。独立完整性检查23项通过（样本N、效应、BH、集合覆盖、leading edge、有限P、TF minsize等）。',
 '', '八、当前准确阻碍与下一步',
 '本批宿主模块不存在未完成计算。解释限制仍是n=5/组、癌症手术对照、缺失临床/批次/肠段协变量；原R环境警告保留且已用独立数值核验降低影响疑虑。',
 '动物GSE319612/MTBLS13907/PRJNA1348869仍为Rothia外部干预：ISA样本与assay组标签冲突未解，阻塞组间分析及3.8；跨组学样本号缺乏真实crosswalk，另行阻塞3.9配对关联。',
 'PMC9366892预筛21代谢物的实际表格、去标识化和C6 Cholic acid两表冲突仍待处理；不能宣称完整无偏代谢组或完整ORA。原有人名不输出。',
 '微生物项目仍缺经核验run/样本/组别/协变量和处理后矩阵；不启动大规模读段下载或长计算。PRJNA1048321/1031121队列与物种/分组独立性仍待核验。',
 '推荐下一项小任务：以本次全量宿主表连接已冻结Model A 331靶点/294交集，做覆盖和证据指向登记，同时继续补齐GSE245885样本来源信息。此项尚未执行，不等于重新预测或重跑Model A；也不能称为本方逆转验证。',
 '', '九、结果位置',f'服务器：sgs10（实际主机osv01）；项目根：{r}',
 '全结果：05_human_RNA/R_DESeq2_bridge_v1、host_functional_v1、fgsea_v1、sample_QC_v1、technical_validation_v1。',
 '汇总：12_reports/first_batch_20261012；修订后交付图：12_reports/first_batch_figures；日志：logs/jobs。',
 'Windows交付：C:\\Users\\31937\\Downloads\\改良通便丸_正式分析_20261012\\07_first_batch_analysis。之前准备文件与旧Model A/RNA试分析保持不变。',
 '打开本目录首批正式分析_报告.html可查看摘要、表格入口和图。文件manifest含相对路径、大小与SHA256；旧大输入不重复打包，见输入路径/哈希清单。',
 '', '方法及资源出处：']
for z in sources:lines.append(z['source']+' '+z['url'])
report='\n'.join(lines)+'\n';(out/'首批正式分析_结果与限制.txt').write_text(report,encoding='utf-8-sig')
def table(rs,keys):
 return '<table><thead><tr>'+''.join('<th>'+html.escape(k)+'</th>' for k in keys)+'</tr></thead><tbody>'+''.join('<tr>'+''.join('<td>'+html.escape(str(z[k]))+'</td>' for k in keys)+'</tr>' for z in rs)+'</tbody></table>'
sections=[]
for block in report.split('\n\n'):
 sections.append('<p>'+html.escape(block).replace('\n','<br>')+'</p>')
figs=[('Figure_sample_QC.png','库大小、size factors 与PCA'),('Figure_sample_correlation.png','样本相关性'),('Figure_ssGSEA_four_modules.png','四个预设模块：真实ssGSEA'),('Figure_fgsea_top_pathways_display_only.png','fgsea按P前列通路，展示选择不改变检验全集'),('Figure_VIPER_top20_display_only.png','VIPER按P前20TF；全270TF推断均保留'),('Figure_dispersion_and_P.png','离散度与P值分布')]
doc='<!doctype html><html lang="zh-CN"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>GSE245885首批正式分析</title><style>body{max-width:1100px;margin:36px auto;padding:0 24px;font:16px/1.75 system-ui,"Microsoft YaHei",sans-serif;color:#243744;background:#f7f9fa}h1{font-size:28px}p{background:white;padding:18px 24px;border-left:4px solid #2f6f8f;white-space:normal;overflow-wrap:anywhere}table{border-collapse:collapse;width:100%;background:white;font-size:14px}td,th{border:1px solid #d6e0e6;padding:8px;text-align:left}th{background:#e8f0f4}img{width:100%;background:white}figure{margin:28px 0}figcaption{font-weight:600}a{color:#17658c}</style><h1>GSE245885 首批正式宿主分析</h1><p>已完成并逐模块核验。5 STC vs 5癌症手术正常肠段对照；观察性疾病证据，不是本方疗效证据。</p>'+table(summary_rows,['module','eligible_or_scored','P_lt_005','BH_lt_005'])
doc+='<p><a href="four_prespecified_modules_summary.tsv">四模块表</a> · <a href="all_270_TFs_ranked_by_P.tsv">全部270 TF</a> · <a href="../../05_human_RNA/fgsea_v1/fgsea_all_collections_ALL_results.tsv">全部6,301通路结果</a> · <a href="job_completion_receipts.tsv">任务完成凭据</a> · <a href="首批正式分析_结果与限制.txt">完整文本报告</a></p>'
doc+=''.join(sections)
for filename,caption in figs:doc+='<figure><figcaption>'+html.escape(caption)+'</figcaption><img src="../first_batch_figures/'+filename+'" alt="'+html.escape(caption)+'"></figure>'
doc+='</html>';(out/'首批正式分析_报告.html').write_text(doc,encoding='utf-8')
# Current state is updated with history retained; no old statistical result changed.
statefile=r/'00_protocol/execution_state.json';cfgfile=r/'00_protocol/project_config.json'
archive=r/'00_protocol/history_before_first_batch_completion';archive.mkdir(exist_ok=True)
for p in [statefile,cfgfile,r/'00_protocol/module_status.tsv',r/'00_protocol/blockers.tsv']:
 if p.exists() and not (archive/p.name).exists():shutil.copyfile(p,archive/p.name)
state=readj('00_protocol/execution_state.json') if statefile.exists() else {}
state.update(first_batch_status='COMPLETED_WITH_DISCLOSED_ENVIRONMENT_WARNING',first_batch_completed_utc=now,server_computation_started=True,formal_analysis_started=True,current_scope='GSE245885 human host modules 3.3-3.7 plus QC, complete',first_batch_receipts=str(out/'job_completion_receipts.tsv'),ModelA_rerun=False,active_first_batch_jobs=0)
state['preparation_stage']='COMPLETE_AND_FIRST_HUMAN_BATCH_COMPLETE'
state['completed']=['previous preparation and frozen input validation','GSE245885 R DESeq2 and PyDESeq2 bridge audit','full human MSigDB and DoRothEA resource freeze','fgsea C2:CP and GO BP','real ssGSEA four fixed modules','VIPER 270 TF scores','limma 4-module and 270-TF comparisons','sample/count/design/PCA/correlation/Cook QC','independent numerical and output integrity validation']
state['not_completed']=['animal group-label conflict resolution and analysis','cross-omics pairing resolution','deidentified selected human metabolite input/C6 conflict audit','microbiome run/phenotype inputs and downstream models','Model A linkage to new host evidence (not started)']
state['recommended_next_module']='freeze host evidence and map to existing Model A targets without rerunning Model A; cohort metadata audit remains necessary'
statefile.write_text(json.dumps(state,ensure_ascii=False,indent=2))
cfg=readj('00_protocol/project_config.json');cfg['execution'].update(formal_analysis_started=True,first_batch_status='COMPLETED',current_authorized_scope='GSE245885 batch completed; no macro/animal/MD computation started')
cfgfile.write_text(json.dumps(cfg,ensure_ascii=False,indent=2))
mp=r/'00_protocol/module_status.tsv'
if mp.exists():
 ms=rows('00_protocol/module_status.tsv')
 for z in ms:
  if z.get('method_id') in ['3.3','3.4','3.5','3.6','3.7']:
   z['status']='HUMAN_GSE245885_COMPLETED_ANIMAL_NOT_RUN';z['executed_this_run']='True';z['data_condition']='GSE245885 first batch all completed, exact outputs and coverage archived';z['blocker']='No human compute blocker; n=5/group and missing covariates; animal group/resource gates remain'
 wt(mp,ms)
bp=r/'00_protocol/blockers.tsv'
if bp.exists():
 bs=rows('00_protocol/blockers.tsv')
 for z in bs:
  if z.get('id')=='B03':z.update(status='RESOLVED_HUMAN_BRIDGE',evidence='R/Python bridge completed; 9 R-only and 1 Python-only nominal calls; effect directions agree',unlock_condition='completed',blocks='none for this human batch')
  if z.get('id')=='B09':z.update(status='RESOLVED_HUMAN_ONLY',evidence='Complete human MSigDB2026.1.Hs and DoRothEA1.22.0 A-C frozen and validated',unlock_condition='human completed; animal resources and groups separate',blocks='no human resource block; animal resource gates remain')
 wt(bp,bs)
status={'status':'COMPLETED','updated_utc':now,'job_status':jobs,'statistical_summary':summary_rows,'warnings_retained':True,'all_integrity_checks_pass':True,'no_modelA_rerun':True,'no_sample_removed':True,'Windows_sync_status':'See separate windows_sync_receipt.json after checksum-verified transfer'}
(out/'first_batch_completion.json').write_text(json.dumps(status,indent=2))
print('REPORTS_AND_STATE_COMPLETE',str(out))
print('REPORT_CHARS',len(report),'HTML_BYTES',(out/'首批正式分析_报告.html').stat().st_size)
