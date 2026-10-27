import pathlib,json,sys,hashlib,datetime,csv,subprocess
sys.path.append('/usr/lib/python3/dist-packages')
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.platypus import SimpleDocTemplate,Paragraph,Table,TableStyle,PageBreak,Image,Spacer
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib import colors
from xml.sax.saxutils import escape
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
N=pathlib.Path(__file__).parent; R=N.parent.parent;C=R/'08_microbiome/PRJNA1232352_covariate_recheck_20261027_v1';A=R/'08_microbiome/PRJNA1232352_association_20261027_v1'
pilot=json.loads((N/'PILOT_COMPLETE.json').read_text()); fun=json.loads((N/'FUNCTIONAL_CANDIDATES_COMPLETE.json').read_text()); search=json.loads((C/'COVARIATE_SEARCH_COMPLETE.json').read_text())
assert pilot['status']=='TECHNICAL_PILOT_PASS' and (N/'exit_code.txt').read_text().strip()=='0'
pdfmetrics.registerFont(TTFont('CJK','/usr/share/fonts/truetype/arphic/uming.ttc'))
body=ParagraphStyle('body',fontName='CJK',fontSize=10,leading=16,spaceAfter=9,wordWrap='CJK');head=ParagraphStyle('head',parent=body,fontSize=15,leading=22,spaceAfter=12);small=ParagraphStyle('small',parent=body,fontSize=8,leading=12); title=ParagraphStyle('title',parent=head,fontSize=20,leading=28)
story=[]
def p(t,st=body):story.append(Paragraph(escape(t),st))
def tab(rows,widths):
 t=Table(rows,colWidths=widths,repeatRows=1);t.setStyle(TableStyle([('FONTNAME',(0,0),(-1,-1),'CJK'),('FONTSIZE',(0,0),(-1,-1),9),('BACKGROUND',(0,0),(-1,0),colors.HexColor('#e9eff5')),('GRID',(0,0),(-1,-1),.4,colors.HexColor('#b7c4d1')),('TOPPADDING',(0,0),(-1,-1),8),('BOTTOMPADDING',(0,0),(-1,-1),8)]));story.append(t)
p('便秘项目公开资料复核与功能分支进展',title);p('内部补充记录｜PRJNA1232352｜2026-10-27')
p('已找到组级人口学资料，未找到可匹配的个体协变量',head)
p('对应论文为 Hu 等发表于 Frontiers in Microbiology 的研究，DOI 10.3389/fmicb.2025.1591697，PMID 40469730，PMCID PMC12134075。正文 Table 1 报告了年龄、女性比例和 BMI 的组级汇总。')
tab([['变量','DCTT（20例）','NCTT（20例）','HC（20例）'],['年龄（岁，均值±SD）','53.55 ± 17.03','55.40 ± 14.31','52.80 ± 9.85'],['女性比例','60%','60%','65%'],['BMI（均值±SD）','20.75 ± 1.00','21.12 ± 0.73','20.77 ± 1.01']],[150,105,105,105])
p('这些汇总值不能替代逐人的年龄、性别或 BMI。已核查 BioProject、60 份 BioSample XML、论文全文及唯一补充 DOCX。补充文件包含测序质控、组装统计、菌群–代谢物及临床指标相关性结果和菌群丰度图，未提供与样本编号对应的个体人口学表。')
p('处理决定：结束当前协变量检索；保留已验收的 group-only 模型，不用组均值填充个人数据，不虚构或跨队列拼接协变量。原关联分析、原报告与原输入保持不变。')
p('新增来源证据：60 个样本与论文测序表全部匹配',head)
p('补充表 S1 中的 CleanReads/CleanBases 与公开库 60 个样本的 reads/碱基数全部一致，同时核实 HC 与归档 H 样本别名的对应。因此，本项目下载的输入是作者已质控后提交的测序数据；项目目录 raw_reads 表示下载输入位置，不代表未经作者处理的仪器原始读段。后续本地 fastp、去宿主及 profiling 处理记录均保留。')
p('来源：论文正文 Table 1、Supplementary_file_1.docx 的 S1–S4、ENA run manifest 和 BioSample 快照。检索结论指本次检查的公开材料，不推断作者未持有个体资料。',small)
story.append(PageBreak());p('下一阶段已完成候选通路与物种贡献整理',head)
p(f"在原始目标对比 P<0.05 的规则下，保留 {fun['candidate_contrast_rows']} 条通路–对比关联，涉及 {fun['unique_candidate_pathways']} 条唯一通路。所有候选保留系数、原始 P、BH 及检验家族，未以 BH 作为下游硬门槛。")
p(f"其中 {fun['pathways_with_stratified_contributors']} 条通路具有 HUMAnN 物种分层记录；已按组汇总贡献丰度的均值、中位数、检出样本数，并保留完整贡献表和每组前 5 个贡献者表。")
p('以下 2 条候选通路没有物种分层记录，仍保留在候选表中：PWY-5088（L-glutamate degradation VIII，to propanoate）及 PWY-7723（bacterial bioluminescence）。不为其补造贡献物种。')
fig,ax=plt.subplots(figsize=(8,3.5));ax.bar(['With stratified contributors','Without stratification'],[250,2],color=['#327ba8','#b9c4cc']);ax.set_ylabel('Unique nominal candidate pathways');ax.set_ylim(0,280)
for i,v in enumerate([250,2]):ax.text(i,v+5,str(v),ha='center')
fig.tight_layout();fig.savefig(N/'candidate_coverage.png',dpi=160);plt.close(fig);story.append(Image(str(N/'candidate_coverage.png'),width=450,height=197));p('图 1：252 条名义显著候选通路的物种贡献记录可用性。图中仅统计是否存在分层记录，不是新的显著性检验。',small)
p('物种分层贡献采用已验收的通路 CPM 分母；其总和不强制等于群落层通路丰度。贡献表是描述性证据，不将贡献丰度当作独立样本，也不把物种名称或通路名称直接解释为已验证的底物转化机制。')
p('已另外保存通路名称关键词标签供人工浏览；标签不用于修改筛选集合，也不构成特定底物、酶活性或反应方向的验证。')
story.append(PageBreak());p('三组组装与 CAZyme 注释技术试跑通过',head)
p('每组按去宿主后保留读对数选择中位附近样本，选择不依赖差异结果。每个样本使用前 500,000 对读段，顺序执行，单任务 4 线程；MEGAHIT 1.2.9 以 8 GiB 作为内存参数，最短 contig 500 bp。k-mer 列表为 21,29,39,59,79,99,119,141。随后运行 Prodigal 2.6.3 的 meta 模式，以及 dbCAN 5.2.9 的 CAZy-DIAMOND 与 dbCAN-HMM 两种注释方法。采用固定版本默认阈值：DIAMOND E-value 1e-102；dbCAN-HMM E-value 1e-15、覆盖度 0.35。')
support=json.loads((N/'PILOT_OUTPUT_QA.json').read_text())['method_support']
rows=[['组别 / SRR','读对数','contigs','预测蛋白','注释记录','双方法']]
for x in pilot['results']:rows.append([x['group']+' / '+x['run'].replace('SRR',''),str(x['pairs']),str(x['contigs']),str(x['proteins']),str(x['dbcan_overview_rows']),str(support[x['run']]['two_methods'])])
tab(rows,[135,70,65,70,65,60]);p('注释记录为 overview.tsv 的记录数，不等于两种方法共同支持的数量；完整注释证据和工具支持数保留在逐样本输出中。')
p('验收包括：双端标识一致、读对数量、压缩文件完整读取、各步骤退出码、FASTA 标识与蛋白/基因对应、注释条目可回溯至预测蛋白，以及关键结果 SHA256。技术试跑只检验流程可运行及输出完整性，不代表全样本组装质量，也不用于组间生物学推断。')
p('后续进入条件',head)
p('全量组装前须用完整样本评估实际峰值内存、磁盘及运行时间，再确定队列规模和并行数。CAZyme 家族注释属于功能候选，尚不能确认具体底物或本方成分转化。Model B、Model C 和结构评价仍需各自证据门槛。')
p('记录位置',head)
p('covariate_recheck_20261027_v1/COVARIATE_SEARCH_COMPLETE.json 保存公开资料核查；functional_pilot_20261027_v1/FUNCTIONAL_CANDIDATES_COMPLETE.json 保存候选整理验收；PILOT_COMPLETE.json、exit_code.txt 与各样本 VALIDATED.json 保存技术试跑凭据。',small)
p('公开来源：https://www.frontiersin.org/journals/microbiology/articles/10.3389/fmicb.2025.1591697/full；https://www.ebi.ac.uk/europepmc/webservices/rest/PMC12134075/supplementaryFiles。',small)
def foot(c,d):c.setFont('CJK',8);c.drawString(42,25,'便秘项目 | 公开资料复核与功能分支 | 内部记录');c.drawRightString(550,25,str(d.page))
SimpleDocTemplate(str(N/'tongbian_covariate_functional_update_20261027.pdf'),pagesize=(595.28,841.89),leftMargin=55,rightMargin=55,topMargin=40,bottomMargin=48).build(story,onFirstPage=foot,onLaterPages=foot)
print('REPORT_CREATED')
