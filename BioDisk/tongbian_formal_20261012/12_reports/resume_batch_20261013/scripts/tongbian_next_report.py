from pathlib import Path
import pandas as pd,numpy as np,json,shutil,hashlib
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from reportlab.platypus import SimpleDocTemplate,Paragraph,Spacer,Table,TableStyle,Image,PageBreak
from reportlab.lib.styles import getSampleStyleSheet,ParagraphStyle
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.lib import colors
from PIL import Image as PI
import pymupdf
P=Path(r'C:\Users\31937\Downloads\改良通便丸_正式分析_20261012');D=P/'11_animal_and_next_modules_20261013/resume_batch_20261013';A=D/'GSE319612_DESeq2_v1';F=D/'GSE319612_functional_integrative_v1';G=D/'report_figures';G.mkdir(exist_ok=True)
rd=lambda p:pd.read_csv(p,sep='\t');meta=rd(A/'audit/metadata_verified.tsv');ss=rd(F/'tables/ssGSEA_scores.tsv').set_index('feature');ss=ss[meta.sample_id]
names={'WATER_SALT_TRANSPORT':'Water / salt transport','SMOOTH_MUSCLE_CONTRACTION':'Smooth muscle contraction','ENTERIC_NERVOUS_SYSTEM_DEVELOPMENT':'Enteric nervous system','BICELLULAR_TIGHT_JUNCTION_ASSEMBLY':'Tight junction assembly'}
cn={'WATER_SALT_TRANSPORT':'水盐转运','SMOOTH_MUSCLE_CONTRACTION':'平滑肌收缩','ENTERIC_NERVOUS_SYSTEM_DEVELOPMENT':'肠神经系统发育','BICELLULAR_TIGHT_JUNCTION_ASSEMBLY':'双细胞紧密连接组装'}
fig,axs=plt.subplots(2,2,figsize=(9,6));pal={'NC':'#32788a','MC':'#c47a4b','RG':'#628c56'}
for ax,(feature,row) in zip(axs.flat,ss.iterrows()):
 for i,gr in enumerate(['NC','MC','RG']):
  y=row[meta.loc[meta.group==gr,'sample_id']].to_numpy(float);x=i+np.linspace(-.14,.14,len(y));ax.scatter(x,y,c=pal[gr],s=23);ax.hlines(y.mean(),i-.2,i+.2,color='black',lw=2)
 ax.set_xticks([0,1,2],['NC (8)','MC (8)','RG (10)']);ax.set_title(names[feature],fontsize=11);ax.set_ylabel('ssGSEA score');ax.grid(axis='y',alpha=.15)
fig.tight_layout();fig.savefig(G/'Figure_ssGSEA_all_modules.png',dpi=210);plt.close(fig)
fig,axs=plt.subplots(1,2,figsize=(11,5))
for ax,contrast in zip(axs,['MC_vs_NC','RG_vs_MC']):
 t=rd(F/'tables'/f'C2_CP_{contrast}_fgsea_ALL.tsv').nsmallest(8,'pval').sort_values('NES');labs=[x.replace('REACTOME_','').replace('KEGG_','').replace('WP_','').replace('_',' ') for x in t.pathway];labs=[x[:41]+'...' if len(x)>44 else x for x in labs];ax.barh(range(len(t)),t.NES,color=['#b66643' if v>0 else '#357d95' for v in t.NES]);ax.set_yticks(range(len(t)),labs,fontsize=8);ax.axvline(0,c='gray',lw=.7);ax.set_xlabel('NES');ax.set_title(contrast.replace('_',' '))
fig.tight_layout();fig.savefig(G/'Figure_fgsea_display_top8_C2CP.png',dpi=220);plt.close(fig)
conf=rd(D/'input_audit/MTBLS13907_sample_assay_conflicts.tsv');gr=['NC','MC','RG','QC','blank'];cross=pd.crosstab(conf.sample_table_group,conf.assay_prefix_group.replace({'QC1':'QC','QC2':'QC','QC3':'QC'})).reindex(index=gr,columns=gr,fill_value=0)
fig,ax=plt.subplots(figsize=(5.5,4.5));ax.imshow(cross,cmap='Blues',vmin=0,vmax=10);ax.set_xticks(range(5),gr);ax.set_yticks(range(5),gr);ax.set_xlabel('Assay / extract prefix');ax.set_ylabel('Sample-table group')
for i in range(5):
 for j in range(5):ax.text(j,i,str(cross.iloc[i,j]),ha='center',va='center',color='white' if cross.iloc[i,j]>=6 else 'black')
fig.tight_layout();fig.savefig(G/'Figure_metabolite_label_conflicts.png',dpi=200);plt.close(fig)
pdfmetrics.registerFont(TTFont('CN',r'C:\Windows\Fonts\simsun.ttc',subfontIndex=0));st=getSampleStyleSheet();st.add(ParagraphStyle(name='B',fontName='CN',fontSize=10.3,leading=16,spaceAfter=8,wordWrap='CJK'));st.add(ParagraphStyle(name='S',fontName='CN',fontSize=8.1,leading=11.5,spaceAfter=6,wordWrap='CJK'));st.add(ParagraphStyle(name='H',fontName='CN',fontSize=16,leading=23,spaceAfter=10));st.add(ParagraphStyle(name='T',fontName='CN',fontSize=23,leading=31,spaceAfter=14));story=[]
def p(t,k='B'):story.append(Paragraph(t.replace('−','-'),st[k]))
def h(t):p(t,'H')
def page():story.append(PageBreak())
def table(rows,widths,size=8.2):
 data=[[Paragraph(str(v),ParagraphStyle(name='c',fontName='CN',fontSize=size,leading=size+3,wordWrap='CJK')) for v in row] for row in rows];t=Table(data,colWidths=widths,repeatRows=1,hAlign='LEFT');t.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,0),colors.HexColor('#dfebef')),('ROWBACKGROUNDS',(0,1),(-1,-1),[colors.white,colors.HexColor('#f4f7f7')]),('VALIGN',(0,0),(-1,-1),'TOP'),('LINEBELOW',(0,0),(-1,0),.6,colors.HexColor('#356c7c')),('TOPPADDING',(0,0),(-1,-1),4),('BOTTOMPADDING',(0,0),(-1,-1),4)]));story.append(t);story.append(Spacer(1,7))
def pic(path,maxh=350,w=500):
 z=PI.open(path);scale=min(w/z.width,maxh/z.height);story.append(Image(str(path),z.width*scale,z.height*scale));story.append(Spacer(1,6))
fmt=lambda x:f'{x:.3g}'
p('改良通便丸公开多组学研究','T');p('动物宿主分析与后续模块<br/>正式结果、输入审计及接续报告','T');p('GSE319612 | NC 8、MC 8、RG 10 | 2026-10-13')
h('双对比识别98个反向变化基因，跨物种整体同向富集未获支持')
p('Rothia外部干预RNA数据已完成DESeq2、真实ssGSEA、VIPER、limma、fgsea及人鼠集合重叠检验。疾病对比MC-NC与干预对比RG-MC共有98个基因满足两个原始P均＜0.05且log2FC异号。该筛选使用预设名义P，不附加fold-change阈值。')
summary=rd(F/'tables/analysis_summary.tsv');rows=[['方法/检验族','疾病 MC-NC<br/>有效 / P / BH','干预 RG-MC<br/>有效 / P / BH'],['DESeq2','16961 / 1247 / 9','16961 / 1494 / 4']]
for fam,title in [('ssGSEA','4个ssGSEA模块'),('VIPER','VIPER TF活性'),('fgsea_C2_CP','fgsea C2:CP'),('fgsea_GO_BP','fgsea GO BP')]:
 rr=[]
 for c in ['MC_vs_NC','RG_vs_MC']:
  t=summary[(summary.family==fam)&(summary.contrast==c)].iloc[0];rr.append(f'{t.tested} / {t.nominal} / {t.BH}')
 rows.append([title,*rr])
table(rows,[180,160,160]);p('表1．P与BH列为各自＜0.05的数量；每个方法家族与对比独立校正。VIPER及ssGSEA干预信号为名义显著，未通过BH。ssGSEA与VIPER均无“两对比同时名义显著且反向”的候选。','S')
p('人鼠比较以12,367个共同可检验的高可信一对一同源基因为背景。上调重叠190个，OR=1.165、单侧P=0.0608、BH=0.1216；下调重叠130个，OR=0.788、P=0.9890、BH=0.9890。当前数据不支持疾病转录变化的整体同向富集。')
page();h('Methods：使用已核实的动物组别，保留预设统计规则')
p('GEO SOFT明确给出26个样本的文库名和分组。NC1-NC8为正常对照；MC9-MC16为洛哌丁胺模型；RG17-RG26为Rothia mucilaginosa干预。矩阵列名与文库名逐一连接，全部为非负整数count，54,532个唯一基因。研究对象为雄性C57BL/6小鼠，样本为结肠组织。')
p('DESeq2 1.50.2采用~group设计，26个独立动物、设计秩3；count≥10且至少8个样本满足的基因进入拟合，共17,029个。两个Wald对比为MC-NC与RG-MC。关闭独立过滤，保留默认Cook过滤；自动异常count替换预先禁用。原人RNA采用最小组5例的过滤与minReplicatesForReplace=7（每组5例不触发替换）；动物使用最小组8例并禁用替换，差异已记录，不为追求显著性调整参数。')
p('VST采用blind=TRUE。相同gene symbol的VST取均值；同名排序基因按最高baseMean及gene ID确定性选择。fgsea输入为全部可检验基因的有符号Wald，不按P预筛。沿用完整冻结MSigDB 2026.1.Hs C2:CP和GO BP，使用Ensembl Genes 116、置信度1且ortholog_one2one的人鼠映射转换为鼠符号；映射表、歧义符号排除及全部覆盖度随附。未下载或替换MSigDB版本。')
p('fgseaMultilevel采用minSize=15、maxSize=500、eps=0、nPermSimple=10000、种子20261012及单进程。真实ssGSEA采用GSVA 2.4.9，minSize=5、alpha=0.25、normalize=TRUE。4模块检测覆盖为水盐243、收缩89、肠神经13、紧密连接43个基因。')
p('VIPER使用冻结DoRothEA 1.22.0鼠A-C网络，decoupleR::run_viper设minsize=5、pleiotropy=TRUE、eset.filter=FALSE、method="scale"、bootstraps=0、cores=1。269个TF返回完整26例评分。两类评分分别使用limma 3.66.0估计两个对比、95%区间、P及BH。')
p('人鼠Fisher使用同一背景内各自P＜0.05的上调、下调集合，预设单侧greater检验；两方向作为独立于其他分析的2项检验族计算BH。方向比较仅使用人STC-对照与鼠MC-NC，不以干预方向替代疾病方向。23项数值及集合完整性检查通过。')
page();h('样本质控显示组内差异，全部26例保留')
pic(A/'figures/Figure_animal_RNA_PCA.png',420)
p('图1．17,029个预过滤基因的blind VST PCA。PC1解释25.1%，PC2解释14.3%；点为独立动物，颜色表示组别。各组存在交叠和组内离散，PCA不是分组正确性或药效验证；没有按图形或结果选择删样。','S')
p('全基因文库计数约1,397万至2,128万，检测基因22,405至24,246个，size factor约0.909至1.348。样本分组依据GEO文本及精确文库名，未凭PCA聚类反推标签。完整PCA坐标、解释率及样本QC表随附。')
h('疾病与干预差异均以名义P保留候选')
p('MC-NC有670个上调、577个下调名义显著基因；RG-MC有596个上调、898个下调名义显著基因。两个对比各16,961个基因具有有效P，原始全基因表保留未过滤与不可检验行，BH在各自有效检验家族内计算。')
page();h('98个基因符合反向模式，保留两个效应及统计证据')
pic(A/'figures/Figure_animal_reverse.png',390)
p('图2．横轴为MC-NC log2FC，纵轴为RG-MC log2FC；红点表示两对比均P＜0.05且效应异号的98个候选，灰点为其他可显示基因。图中不设fold-change门槛；完整候选与全部54,532基因双对比结果均在TSV中。','S')
r=rd(A/'tables/opposite_direction_nominal_candidates.tsv');r['display_P']=r[['MC_vs_NC_P','RG_vs_MC_P']].max(axis=1);top=r.nsmallest(8,'display_P')
table([['基因','MC-NC log2FC / P / BH','RG-MC log2FC / P / BH']]+[[str(x.gene_name),f'{x.MC_vs_NC_log2FC:.3f} / {fmt(x.MC_vs_NC_P)} / {fmt(x.MC_vs_NC_BH)}',f'{x.RG_vs_MC_log2FC:.3f} / {fmt(x.RG_vs_MC_P)} / {fmt(x.RG_vs_MC_BH)}'] for x in top.itertuples()],[115,193,192],7.8)
p('表2．仅为显示排序：按两个P的较大值取前8个；分析集合仍为全部98个。共享MC组使两个对比相关，反向效应不等同恢复正常，也不证明本方作用。','S')
page();h('干预伴随收缩评分升高与紧密连接组装评分下降')
pic(G/'Figure_ssGSEA_all_modules.png',310)
p('图3．4个预设ssGSEA模块的全部样本评分，点为独立动物，黑线为组均值。未按显著性挑选模块或删除样本；各模块的符号映射和检测覆盖度完整保留。','S')
rows=[['模块','对比','Δscore','P','BH']]
for c,label in [('MC_vs_NC','疾病'),('RG_vs_MC','干预')]:
 t=rd(F/'tables'/f'ssGSEA_{c}_ALL.tsv')
 for x in t.itertuples():rows.append([cn[x.feature],label,f'{x.delta_score:.4f}',fmt(x.P),fmt(x.BH)])
table(rows,[170,60,90,90,90],8)
p('表3．疾病对比4模块均未达P＜0.05。干预对比平滑肌收缩Δ=0.0223、P=0.0325，紧密连接组装Δ=-0.0132、P=0.0427；两者BH=0.0855，均为名义显著探索性信号。不能将紧密连接评分下降省略或描述为整体屏障改善。VIPER两对比16及30个名义显著TF均未通过BH；完整效应、区间及全部269个TF已交付。','S')
page();h('全排序通路分析提供功能线索，未以差异基因预筛')
pic(G/'Figure_fgsea_display_top8_C2CP.png',320)
p('图4．C2:CP两个对比各按P取前8条通路用于显示，横轴为NES；正值表示对比中前者富集，负值表示后者富集。图的显示集合不用于统计筛选或重新校正，完整C2:CP与GO BP结果、leading-edge及全部资源覆盖表均随附。','S')
p('C2:CP每个对比2,264项有效检验，疾病569项名义显著、209项BH显著；干预864项名义显著、581项BH显著。GO BP每个对比3,628项有效检验，疾病807项名义显著、282项BH显著；干预1,102项名义显著、537项BH显著。原冻结全集保留，实际检验数由鼠同源映射、可检验排序覆盖及预设15-500基因范围决定。')
h('人鼠集合检验未支持整体同向富集')
f=rd(F/'tables/human_mouse_directional_Fisher.tsv');table([['方向','人集合','鼠集合','重叠','OR','P / BH']]+[[x.direction,x.human_set_n,x.mouse_set_n,x.overlap_n,fmt(x.OR),f'{fmt(x.P)} / {fmt(x.BH)}'] for x in f.itertuples()],[60,65,65,65,75,170])
p('表4．背景为两队列共同可检验的12,367个一对一同源基因。单侧Fisher精确检验与独立超几何尾概率计算一致。此阴性结果限制跨物种整体一致性主张，但不否定各队列已经观察到的具体效应。','S')
page();h('动物代谢组确认22条组标签冲突，暂不进入差异分析')
pic(G/'Figure_metabolite_label_conflicts.png',290)
p('图5．ISA sample组别与同一assay行中extract/assay名称前缀的交叉表。NC/MC互换，部分QC/RG也冲突。30条记录中22条组标签冲突；同组内仍有编号错位，不能据相同组别恢复个体配对。','S')
p('已取得完整ISA、5,308行MAF注释以及正/负离子处理强度文件。负离子2,115峰、正离子3,193峰；各有26个动物列和3个QC列，峰ID在模式内唯一，强度无缺失且为正数。文件扩展名虽为.xls，实际为GB18030编码的制表文本。')
p('来源协议已做KNN填补、按总信号归一化至QC1并按QC CV筛选；本轮未重复归一化或填补。原始组标签、提取名、assay名与强度列均保留，不根据RNA方向、图形聚类或编号顺序选择哪套标签正确。差异、反向变化及样本级相关仍未执行。')
p('PRJNA1348869含27个run：26个扩增子样本和1个菌株WGS记录，不能把整个项目当作小鼠宏基因组。RNA为NC1-NC8、MC9-MC16、RG17-RG26；菌群为各组内编号。相似数字不建立跨组学个体配对。')
page();h('核心宏基因组已后台直连下载，明日从状态接续')
p('核心PRJNA1232352的60个WGS样本包含DCTT、NCTT及H前缀各20；ENA清单给出120个压缩FASTQ，共293,286,068,526字节（约293.3 GB）。BioData本轮核查可用约6.09 TB；大数据均写入既定BioData根目录。实际大文件续传返回HTTP 206。')
p('用户最新安排为先后台直连、明日再连接代理；慢速直连已获明确接受。下载器主PID 103511，4路并行，逐文件保存.partial并原生等待curl，完整文件以预期大小、ENA MD5及gzip -t共同验证后改名。后台状态：logs/jobs/PRJNA1232352_direct_download_v1/status.json。未创建代理隧道，未改SSH、代理、防火墙或密钥配置。')
p('明日先核对下载器PID、状态、日志和现有文件。只有确认本项目旧写入进程已停止，才可切换任务代理并续传同一.partial；不得并发双写或重复下载已校验文件。所有文件校验通过后才生成DOWNLOAD_COMPLETE.json；本报告不将下载、MetaPhlAn/HUMAnN、组装或Model B标记完成。')
h('16S独立性未确认，暂不作为独立重复验证')
p('PRJNA1048321为60个扩增子样本，PRJNA1031121为42个。两研究均来自成都中医药大学相关团队，并共用伦理编号2021KL-023；原始FASTQ的MD5没有完全相同文件，但这一点不能证明受试者独立。项目保留队列独立性门槛，不将两个研究直接计为独立验证，更不以16S直接证明菌酶基因。')
h('Discussion与项目状态')
p('本批动物数据属于Rothia干预外部证据，人RNA对照为癌症手术正常肠段。各队列独立分析，结果汇入Model C的不同证据层；Model C尚非已训练模型。现阶段已有候选关联、宿主效应及明确阴性证据，转化层仍需实验依据酶、底物和生成/消耗反应支持。尚未进入GraphBAN、对接、动力学及结合能评估；整个项目未完成。')
p('来源：GEO GSE319612 SOFT及count；DOI 10.1186/s12866-026-05126-0；MetaboLights MTBLS13907官方FTP；ENA四项目run清单；Ensembl Genes 116 BioMart；DOI 10.3389/fmicb.2023.1320567、10.3389/fmicb.2023.1323877。精确URL、下载日志、版本、脚本、原始标签和SHA256见本批附件。','S')
pdf=D/'改良通便丸_动物宿主分析与后续审计_完整报告_20261013.pdf'
def footer(c,doc):c.setFont('CN',8);c.setFillColor(colors.HexColor('#567078'));c.drawString(43,25,'改良通便丸公开多组学研究 | Rothia外部证据与输入审计 | 2026-10-13');c.drawRightString(553,25,str(doc.page))
SimpleDocTemplate(str(pdf),pagesize=(595.28,841.89),leftMargin=43,rightMargin=43,topMargin=42,bottomMargin=43).build(story,onFirstPage=footer,onLaterPages=footer)
Q=D/'report_QA';Q.mkdir(exist_ok=True);doc=pymupdf.open(pdf)
for i,pag in enumerate(doc):pag.get_pixmap(matrix=pymupdf.Matrix(1.2,1.2)).save(Q/f'page_{i+1:02}.png')
(Q/'text.txt').write_text('\n'.join(x.get_text() for x in doc),encoding='utf-8');shutil.copy2(__file__,D/'scripts'/Path(__file__).name);print('REPORT',pdf,'PAGES',len(doc))
