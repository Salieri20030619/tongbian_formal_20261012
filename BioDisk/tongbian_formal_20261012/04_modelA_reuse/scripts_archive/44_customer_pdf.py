"""Approved customer version; wording revision only; figures and analysis unchanged."""
from pathlib import Path
import json,html
import pandas as pd
from reportlab.pdfgen import canvas
from reportlab.platypus import SimpleDocTemplate,Paragraph,Spacer,Table,TableStyle,Image,PageBreak,KeepTogether
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib import colors
from reportlab.lib.enums import TA_LEFT
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.lib.utils import ImageReader
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'output/pdf';OUT.mkdir(parents=True,exist_ok=True)
pdfmetrics.registerFont(TTFont('CN','C:/Windows/Fonts/simhei.ttf'))
pdfmetrics.registerFontFamily('CN',normal='CN',bold='CN',italic='CN',boldItalic='CN')
W,H=595.28,841.89;CW=W-88
navy=colors.HexColor('#203C4A');teal=colors.HexColor('#147D83');orange=colors.HexColor('#CD7449');muted=colors.HexColor('#586A72');light=colors.HexColor('#EEF5F3')
styles={
 'title':ParagraphStyle('title',fontName='CN',fontSize=27,leading=38,textColor=navy,spaceAfter=16),
 'h1':ParagraphStyle('h1',fontName='CN',fontSize=18,leading=26,textColor=navy,spaceAfter=14),
 'h2':ParagraphStyle('h2',fontName='CN',fontSize=12,leading=18,textColor=teal,spaceBefore=10,spaceAfter=7),
 'body':ParagraphStyle('body',fontName='CN',fontSize=10,leading=16.5,textColor=navy,spaceAfter=9,wordWrap='CJK'),
 'small':ParagraphStyle('small',fontName='CN',fontSize=8.2,leading=12.4,textColor=muted,spaceAfter=7,wordWrap='CJK'),
 'cell':ParagraphStyle('cell',fontName='CN',fontSize=8.5,leading=12.5,textColor=navy,wordWrap='CJK'),
 'head':ParagraphStyle('head',fontName='CN',fontSize=8.5,leading=12.5,textColor=colors.white,wordWrap='CJK')}
story=[]
def p(t,s='body'):return Paragraph(t,styles[s])
def add(t,s='body'):story.append(p(t,s))
def new(title):story.append(PageBreak());add(title,'h1')
def table(rows,widths,fs=8.5):
 data=[[Paragraph(html.escape(str(v)).replace('\n','<br/>'),styles['head' if j==0 else 'cell']) for v in row] for j,row in enumerate(rows)]
 t=Table(data,colWidths=widths,repeatRows=1,hAlign='LEFT');t.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,0),navy),('VALIGN',(0,0),(-1,-1),'TOP'),('LEFTPADDING',(0,0),(-1,-1),7),('RIGHTPADDING',(0,0),(-1,-1),7),('TOPPADDING',(0,0),(-1,-1),6),('BOTTOMPADDING',(0,0),(-1,-1),6),('ROWBACKGROUNDS',(0,1),(-1,-1),[colors.white,light]),('LINEBELOW',(0,0),(-1,0),.6,navy)]));story.append(t);story.append(Spacer(1,9))
def img(name,width=CW,maxh=370):
 path=ROOT/'figures'/name;iw,ih=ImageReader(str(path)).getSize();scale=min(width/iw,maxh/ih);story.append(Image(str(path),width=iw*scale,height=ih*scale,hAlign='CENTER'));story.append(Spacer(1,9))
def read(n):return pd.read_csv(ROOT/'tables'/n,sep='\t',keep_default_na=False)
b=json.loads((ROOT/'results/primary_baseline_summary.json').read_text(encoding='utf-8'));ps=json.loads((ROOT/'results/STRING_PPI_summary.json').read_text(encoding='utf-8'));ns=json.loads((ROOT/'results/network_export_summary.json').read_text(encoding='utf-8'));es=json.loads((ROOT/'results/enrichment_summary.json').read_text(encoding='utf-8'));gc=json.loads((ROOT/'results/GeneCards_manual_import_metadata.json').read_text(encoding='utf-8'))
cov=read('37_formula_primary_coverage.tsv');er=read('40_GO_KEGG_ALL_tested_results.tsv');pr=er[er.analysis.str.startswith('PRIMARY')];bg=read('42_enrichment_effective_backgrounds.tsv');pn=read('46_STRING_PPI_nodes.tsv');ex=read('55_reviewable_formula_examples_FULL_provenance.tsv')
add('改良通便丸<br/>网络药理初步验证','title')
add('Model A · 20味方剂原形成分基线 · 客户版','h2')
add('数据获取与计算：2026年10月10日 | 完整20味处方 | 本机直连执行','small')
story.append(Spacer(1,14))
table([['处方覆盖','明确结构候选','人靶点','便秘交集'],['19 / 20味',f"{b['primary_unique_parent_structures']:,}",str(b['primary_human_target_genes']),str(b['disease_intersection_genes'])]],[CW/4]*4)
add('已完成20味方剂的公开数据库原形成分基线：20味中19味有合格结构记录、17味可关联到本次纳入的人靶点证据。2,059个唯一结构中580个具有靶点关联，合计331个人靶点；其中294个与GeneCards便秘基因集重叠，294个交集中292个进入已核验的STRING网络（含3个孤立节点）。')
add('本次主要结果','h2')
add('已形成可追溯的药味-结构-靶点网络，包含2,410个节点和7,161条唯一边。全人背景分析显示，便秘交集在神经活性配体信号、钙信号等功能条目中存在统计富集，并可回溯到厚朴、槟榔、陈皮、枳实及生地黄的具体数据库候选作用关联。')
add('294个交集占331个人靶点的88.82%。进一步以处方靶点为背景，7,786个GO/KEGG条目统一BH校正后，仅1个宽泛GO条目具有额外富集，KEGG额外富集显著数为0。这与全人背景的功能分布结果回答不同问题，尚不支持便秘特异性。')
add('阅读主线','h2')
add('这20味药 → 原形候选结构与人靶点 → 便秘交集与功能解释。第一步已获得可追溯原形基线，可用于后续A/B比较；本轮为公开数据库初步验证，尚未确认本方药效。')
add('本文“原形候选结构”指来源药材数据库中的明确小分子结构。投料、成分记录、靶点证据和疾病交集分别保存，便于后续在一致口径下比较。','small')

new('01  处方原样保留，逐味报告覆盖')
add('投料数值仅按相对量记录，总和340；单位和每剂用量未确认，不写作g/day。各药一半研粉、一半煎汤浓缩制丸，不等同于50%的成分暴露比例，也未用于加权靶点。')
rows=[['药味','相对投料','唯一结构','关联人靶点','便秘交集']]
for r in cov.to_dict('records'):rows.append([r['prescription_label'],r['relative_input'],r['primary_structure_count'],r['primary_human_targets'],r['disease_intersection_targets']])
table(rows,[125,78,90,106,CW-399])
add('计数按每味去重；同一结构或靶点可跨药味出现，故逐味数不能直接相加得到整方唯一数。0表示本次范围内未获得合格记录，不表示该药无成分或无作用。','small')
add('完整标准化、检索和排除记录：tables/02、12-16、28、31、37。','small')

new('02  19味有结构记录，17味有合格靶点')
img('01_herb_coverage.png',maxh=420)
add('药材对应关系的处理','h2')
add('麦冬以HERB“寸冬”条目核对Radix Ophiopogonis及Ophiopogon japonicus，获得9个合格结构；生地黄以“干地黄”条目核对Rehmannia glutinosa干燥块根，获得126个结构。香港中药材标准第3册第168、236页支持对应关系，且明确干燥地黄俗称生地黄。')
add('药材范围按具体物种、部位和炮制要求核对：瓜蒌仁限定种子，铁皮石斛限定Dendrobium officinale茎，生地黄限定未熟制的干燥块根。逐味实际记录数见上图，覆盖边界集中列于第9页。')
add('TCMSP黄芪、苦杏仁的物种标签与处方范围不一致，主集采用HERB匹配药材记录；存在Liriope标签冲突的麦冬记录及鲜地黄记录留作补充。上述规则及官方药材标准来源均保留在药材核验表。','small')

new('03  从数据库记录到可追溯网络')
img('02_full_formula_network.png',maxh=310)
table([['层级','实际纳入','解释'],['药味-结构','2,728条唯一边','19味有结构；20味均作为网络节点保留'],['有靶点的结构','580 / 2,059','其余1,479个结构保留在结构基线中'],['结构-人基因','4,433条唯一配对','485条源库validated标签；3,948条源库SVM/RF预测']],[110,115,CW-225])
add('结构规则与证据边界','h2')
add('按来源标识符取得PubChem结构或原库结构，使用RDKit生成并核对完整InChIKey，保留立体化学差异；混合物、多片段、通配原子、结构冲突及无法解析记录不进入明确单分子网络。无OB/DL硬筛选。47个带_qt标记、缺乏独立原形来源支持的衍生结构在Model A最终核查中排除，保留排除表。')
add('主靶点边来自TCMSP已有证据，经HERB显式TCMSP靶点交叉标识及Taxonomy 9606核对人身份，再经HGNC规范化。485条validated为源数据库状态标签，其证据可涉及表达变化或功能关联，不统一视为直接结合证据。HERB文献挖掘和未分型关联另表保留。','small')
add('本图由NetworkX/matplotlib绘制；节点表、边表和GraphML随报告提供，便于继续查看完整网络。','small')

new('04  数据库候选作用关联：6个具体实例')
img('06_formula_evidence_examples.png',maxh=300)
add('图中6条均来自完整药味-结构-靶点网络，分别展示受体、酶、炎症相关因子及转运体方向的数据库候选作用关联。')
add('厚朴的magnolol、honokiol与CHRM3，以及槟榔的homoarecoline与CHRM2，为源库已有SVM/RF预测，形成胆碱能受体方向的候选关联。该层证据用于界定候选作用关系，不区分激动或拮抗方向。')
add('陈皮hesperidin-PTGS2、枳实naringin-TNF、生地黄γ-aminobutyric acid-SLC6A1具有源库validated标签，分别提供酶、炎症相关因子及转运体方向的关联。该标签包含不同证据类型，可涉及表达或功能变化，不等同于已确认直接结合。')
add('如何核查任一配对','h2')
add('表55保存上述实例的完整InChIKey、来源成分ID、药材条目、成分-靶点原网页、人靶点交叉标识和原始缓存路径。表50保留全部14,022行来源证据；表51为8,745条去重的药味-结构-基因证据组合。')
add('对于任一PPI节点，可用表54回溯其全部药味与结构。共同hub用于描述数据库连接度；药味、成分及证据类型分别保留，可支持后续一致口径的比较。','small')

new('05  便秘交集及真实STRING功能互作')
add('用户提供的GeneCards constipation导出含9,230条有效基因记录、9,230个唯一符号；另1条版权页脚单独排除。保留所有有效记录及相关性分数，不设阈值。HGNC映射后8,877个唯一NCBI基因，与处方331个人靶点交集294个。')
img('03_STRING_PPI.png',maxh=325)
table([['交集输入','可核验STRING节点','真实互作边','孤立节点'],['294','292','4,349','3']],[CW/4]*4)
add('STRING 12.5，Homo sapiens 9606，functional网络，score≥0.4，add_nodes=0。最大连通分量289个节点。ACP3通过HGNC的UniProt标识核实为旧名ACPP；MAPK10、VEGFA返回错误名称或无明确匹配，未接入PPI，仍保留在疾病交集及功能分析。')
top=pn.head(5);add('最高连接度：'+ '；'.join(f"{r.approved_symbol}（{r.degree}）" for r in top.itertuples())+'。这些节点反映本网络的连接特征。STRING组合分数及各证据通道均保留；功能互作与疾病关联可能共享文献，按关联网络进行解释。','small')

new('06  处方内部背景：交集的额外功能富集')
add('本分析回答：相对于处方本身的人靶点功能分布，便秘交集是否具有额外的功能偏重？输入294个交集、背景331个处方靶点；每个功能体系采用有注释的有效背景，分母如下。全人背景的功能分布结果见下一页。')
rows=[['功能体系','有效背景 / 交集','全部测试条目','BH显著条目']]
for fam,name in [('GO_BP','GO生物过程'),('GO_MF','GO分子功能'),('GO_CC','GO细胞组分'),('KEGG','KEGG')]:
 rr=bg[(bg.analysis.str.startswith('PRIMARY'))&(bg.family==fam)].iloc[0];dd=pr[pr.family==fam];rows.append([name,f'{rr.effective_universe_count} / {rr.effective_query_count}',len(dd),int((dd.p_adjust_BH_all_families<.05).sum())])
table(rows,[116,146,123,CW-385])
img('04_primary_enrichment.png',maxh=235)
add('唯一具有显著额外富集的条目为GO:0032501“多细胞生物体过程”：交集174/294、背景180/330，富集倍数1.085，原始P=9.33×10^-7，统一BH q=0.00727。该条目为宽泛上位过程，尚不足以支持便秘特异性。')
add('KEGG共测试325条通路，额外富集显著数为0，最小统一BH q=0.9673。交集覆盖了处方靶点的大部分，处方内部的功能区分度有限。这不表示交集缺乏生物功能，也不否定其相对全人背景呈现的功能分布。')
add('采用单侧超几何检验；全部K≥1的条目均测试，k=0也保留；对每个分析的全部GO+KEGG条目统一BH校正，另存按功能体系BH值。GO沿is_a/part_of传播；未将regulates关系作为注释继承。','small')

new('07  全人背景：交集呈现的功能分布')
img('05_background_comparison.png',maxh=250)
add('敏感性分析将同一294个交集与已注释的人蛋白编码基因比较。HGNC提供19,297个有NCBI ID的蛋白编码基因；GO BP/MF/CC有效背景分别16,560/15,589/18,137，KEGG有效背景8,484。各项注释背景均完整留表。')
add('在全人背景下，交集呈现神经活性配体信号、钙信号、AGE-RAGE及代谢等功能条目的统计富集。该分析描述交集相对于人蛋白编码基因的功能分布；上一页的处方内部分析则检验额外功能偏重。两者结果互补，所回答的问题不同。')
ss=er[(er.analysis.str.startswith('SENSITIVITY'))&(er.family=='KEGG')].head(6)
rows=[['人背景敏感性：KEGG前6条','命中数','统一BH q']]
short={'hsa04933':'AGE-RAGE信号（糖尿病并发症）','hsa04082':'神经活性配体信号','hsa05200':'癌症通路','hsa05417':'脂质与动脉粥样硬化','hsa01100':'代谢通路','hsa05418':'流体剪切力与动脉粥样硬化'}
for rr in ss.itertuples():rows.append([rr.term_id+' '+short.get(rr.term_id,rr.term_name),rr.k_overlap,f'{rr.p_adjust_BH_all_families:.2e}'])
table(rows,[CW-150,60,90])
add('上表按统一BH q值顺序列出前6条通路，完整结果见表40。神经活性配体等功能结果可与第5页具体候选作用关联对应阅读。全部331个处方靶点在人背景下也有广泛富集；结合处方内部仅宽泛GO具有额外富集的结果，目前尚不支持便秘特异性。','small')

new('08  原方案执行状态与整体可行性')
table([['原方案环节','本次实际执行情况'],['TCMSP / HERB成分','完成公开药材检索、全量预定成分详情收集和范围审计；HERB 2,472/2,472条详情完成。'],['PubChem / 结构去重','实际下载CID/InChIKey性质，RDKit核对；明确单分子结构进入主集。'],['SEA / SwissTargetPrediction','本轮未运行新预测：SEA入口需验证码；Swiss条款禁止自动访问。'],['GeneCards疾病基因','采用用户手工导出CSV；保存原件和版权文字，全量纳入有效符号。HPO未进入主分析。'],['STRING PPI','实际调用官方API，导出真实组合分数及各通道；固定0.4阈值、不加节点。'],['Cytoscape网络','未运行软件；实际使用NetworkX/matplotlib，并提供可导入节点、边和GraphML。'],['GO / KEGG','真实注释、明确有效背景、全部测试及BH校正完成；主分析与人背景敏感性分别报告。']],[132,CW-132])
add('本轮结论','h2')
add('第一步已获得可追溯的20味方剂原形成分基线，完成明确结构、人靶点、便秘交集、功能互作及两套背景的功能分析，并提供具体数据库候选作用关联。该基线可用于后续A/B比较；本轮为公开数据库初步验证，尚未确认本方药效。')
add('本轮覆盖边界：瓜蒌仁种子记录不足；麦冬9个结构、铁皮石斛1个结构尚未获得本证据范围内的人靶点。多糖、脂肪油整体及非明确单分子作用未由当前网络代表。原数据库药材字段、衍生结构与证据类别已分层记录；源库validated可包含表达和功能证据，不等同于直接结合。原始药材、本丸剂成分含量及体内或肠腔暴露未实测。本报告范围限定为Model A原形成分基线。','small')

new('09  数据版本、许可及交付索引')
add('所有网络请求在本机进程中禁用HTTP_PROXY/HTTPS_PROXY/ALL_PROXY，requests trust_env=False，NO_PROXY=*，TLS校验开启。未修改系统或Clash设置，未回退代理，未绕过登录、验证码或下载权限。','small')
table([['数据与版本','保存与用途'],['下载日期','2026-10-10 UTC；每次请求的时间、URL、状态和SHA256见logs/http_requests.jsonl。'],['GO / KEGG / STRING','GO releases/2026-07-26；KEGG pathway更新2026-10-09；STRING API 12.5。UniProt响应版本头、HGNC原文件和全人背景均留存。'],['计算环境','Python 3.12.7；pandas 2.2.2；SciPy 1.13.1；statsmodels 0.14.2；NetworkX 3.3；matplotlib 3.9.2；RDKit 2026.9.1；ReportLab 5.0.1。'],['许可备注','TCMSP页面声明ODbL；HERB公开访问不等同额外再分发许可。GeneCards导出版权声明要求非学术使用获得书面许可。PubChem、UniProt、GO、HGNC、STRING、KEGG各自许可/引用要求独立适用，公开接口不授予统一商业许可。原始材料保留用于内部复核，外发或商业交付前需核对相应授权。']],[128,CW-128])
add('先看这些文件','h2')
rows=[['文件 / 目录','内容'],['tables/37、32、34、36','逐味覆盖；唯一结构；完整成分-人靶点边；GeneCards交集'],['tables/50、51、54、55','完整药味-结构-靶点来源；去重组合；PPI节点回溯；6个讲解实例'],['tables/40-43、45-46','全部GO/KEGG测试、背景基因清单、显著结果；真实PPI边与节点'],['tables/52-53、results/*.graphml','Cytoscape可导入的节点、边与GraphML网络文件'],['raw/、scripts/、logs/','原始响应与用户CSV；全部脚本；请求日志、环境、最终校验和'],['results/FINAL_README_中文.txt','结果解释、复算顺序、排除规则、方法修订和许可备注']]
table(rows,[192,CW-192])
add('来源入口','h2')
add('TCMSP: https://old.tcmsp-e.com/　HERB: http://herb.ac.cn/<br/>PubChem: https://pubchem.ncbi.nlm.nih.gov/　HGNC: https://www.genenames.org/<br/>UniProt: https://www.uniprot.org/　GO: https://geneontology.org/<br/>STRING: https://string-db.org/　KEGG: https://www.kegg.jp/<br/>香港中药材标准: https://www.cmro.gov.hk/html/eng/useful_information/hkcmms/volumes.html','small')
add('本报告为Model A公开数据库初步验证客户版。全表、脚本和原始缓存随项目目录保留，支持复核及后续A/B比较。','small')
def page(c,doc):
 c.saveState();c.setFillColor(teal);c.rect(0,H-9,W,9,fill=1,stroke=0);c.setFont('CN',8);c.setFillColor(muted);c.drawString(44,H-31,'改良通便丸 · Model A · 网络药理实际结果');c.drawRightString(W-44,H-31,'客户版 | 2026-10-10');c.setStrokeColor(colors.HexColor('#D4DFDF'));c.line(44,39,W-44,39);c.drawString(44,25,'Model A 原形成分基线 | 公开数据库初步验证');c.drawRightString(W-44,25,str(doc.page));c.restoreState()
path=OUT/'改良通便丸_ModelA_网络药理初步验证_客户版.pdf'
doc=SimpleDocTemplate(str(path),pagesize=(W,H),rightMargin=44,leftMargin=44,topMargin=57,bottomMargin=55,title='改良通便丸 Model A 网络药理初步验证 客户版',author='Reproducible analysis',allowSplitting=1)
doc.build(story,onFirstPage=page,onLaterPages=page);print(path)
