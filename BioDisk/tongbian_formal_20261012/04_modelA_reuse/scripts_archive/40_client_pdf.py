"""Chinese client review draft, calculated values read from final tables."""
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
add('改良通便丸<br/>原形成分网络药理结果','title')
add('Model A · 客户原方案第一步 · 内部审阅初稿','h2')
add('数据获取与计算：2026年10月10日 | 完整20味处方 | 本机直连执行','small')
story.append(Spacer(1,14))
table([['处方覆盖','明确结构候选','人靶点','便秘交集'],['19 / 20味',f"{b['primary_unique_parent_structures']:,}",str(b['primary_human_target_genes']),str(b['disease_intersection_genes'])]],[CW/4]*4)
add('已建立从药味、原形候选结构到人靶点，再到便秘关键词基因交集的可回溯网络。完整网络包含2,410个节点、7,161条唯一边；每条成分-人靶点边均保留源数据库记录和人基因标识。')
add('本次实际结论','h2')
add('数据分析流程可执行，并已产生真实结构表、靶点交集、STRING互作、GO/KEGG计算和网络导入文件。现有结果适合提出待检验的成分-靶点假设；尚不足以确认整方通便机制或临床有效性。')
add('294个交集占331个人靶点的88.82%。以全处方可识别人靶点为背景，7,786个GO/KEGG条目统一BH校正后，仅1个宽泛GO条目显著，KEGG显著通路为0。没有为获得显著结果改变阈值或主背景。')
add('阅读主线','h2')
add('这20味药 → 实际可确认的结构与靶点 → 便秘交集与功能解释。后续页面依次给出覆盖、证据实例、真实互作及背景敏感性。')
add('本文“原形候选结构”指来源药材数据库中的明确小分子结构，不等于本丸剂实测成分、体内暴露或天然存在的逐条原始文献复核。未推导菌转化产物；未开展宏基因组、分子动力学或转录组再分析。','small')

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
add('瓜蒌仁未取得种子专属合格条目，不使用瓜蒌皮或整果代替。铁皮石斛仅取得Dendrobium officinale茎条目中的1个结构，没有合格人靶点。麦冬当前9个结构同样没有进入所用靶点证据集。')
add('源记录仍有字段问题：TCMSP黄芪、苦杏仁的物种标签冲突，相关药味归属从主集排除，采用HERB匹配药材记录；存在Liriope冲突的麦冬记录和鲜地黄记录保留在补充表。数据库药材归属不能替代本批药材鉴定与炮制核验。','small')

new('03  从数据库记录到可追溯网络')
img('02_full_formula_network.png',maxh=310)
table([['层级','实际纳入','解释'],['药味-结构','2,728条唯一边','19味有结构；20味均作为网络节点保留'],['有靶点的结构','580 / 2,059','其余1,479个结构不虚配靶点'],['结构-人基因','4,433条唯一配对','485条源库validated标签；3,948条源库SVM/RF预测']],[110,115,CW-225])
add('结构规则与证据边界','h2')
add('按来源标识符取得PubChem结构或原库结构，使用RDKit生成并核对完整InChIKey，保留立体化学差异；混合物、多片段、通配原子、结构冲突及无法解析记录不进入明确单分子网络。无OB/DL硬筛选。47个带_qt标记、缺乏独立原形来源支持的衍生结构在Model A最终核查中排除，保留排除表。')
add('主靶点边来自TCMSP已有证据，经HERB显式TCMSP靶点交叉标识及Taxonomy 9606核对人身份，再经HGNC规范化。485条validated仅是源数据库标签：未逐条核对实验亲和力、原始实验物种或本制剂效应。HERB文献挖掘和未分型关联另表保留，未混称为实测结合或本项目新预测。','small')
add('本图由NetworkX/matplotlib实际绘制；未运行Cytoscape。已提供节点表、边表和GraphML供Cytoscape导入。','small')

new('04  把靶点落回具体药味和成分')
img('06_formula_evidence_examples.png',maxh=300)
add('这6条是已存在于完整网络中的解释实例，不参与成分筛选或富集阈值选择。')
add('厚朴的magnolol、honokiol与CHRM3，以及槟榔的homoarecoline与CHRM2，为源库已有SVM/RF预测。可作为受体方向的待检验假设，不能据此认定激动/拮抗方向、有效剂量或改善肠动力。')
add('陈皮hesperidin-PTGS2、枳实naringin-TNF、生地黄γ-aminobutyric acid-SLC6A1具有源库validated标签。这里只报告来源状态，不将其解释为本项目亲和力验证或独立的便秘作用证据。')
add('如何核查任一配对','h2')
add('表55保存上述实例的完整InChIKey、来源成分ID、药材条目、成分-靶点原网页、人靶点交叉标识和原始缓存路径。表50保留全部14,022行来源证据；表51为8,745条去重的药味-结构-基因证据组合。')
add('对于任一PPI节点，可用表54回溯其全部药味与结构。共同hub反映数据库连接度，不能独立确定关键有效药味、君臣关系或配伍协同。','small')

new('05  便秘交集及真实STRING功能互作')
add('用户提供的GeneCards constipation导出含9,230条有效基因记录、9,230个唯一符号；另1条版权页脚单独排除。保留所有有效记录及相关性分数，不设阈值。HGNC映射后8,877个唯一NCBI基因，与处方331个人靶点交集294个。')
img('03_STRING_PPI.png',maxh=325)
table([['交集输入','可核验STRING节点','真实互作边','孤立节点'],['294','292','4,349','3']],[CW/4]*4)
add('STRING 12.5，Homo sapiens 9606，functional网络，score≥0.4，add_nodes=0。最大连通分量289个节点。ACP3通过HGNC的UniProt标识核实为旧名ACPP；MAPK10、VEGFA返回错误名称或无明确匹配，未接入PPI，仍保留在疾病交集及功能分析。')
top=pn.head(5);add('最高连接度：'+ '；'.join(f"{r.approved_symbol}（{r.degree}）" for r in top.itertuples())+'。这些是泛用生物学hub，图的稠密程度不代表便秘特异性或药效。STRING组合分数和各证据通道原样保留；文本挖掘与疾病数据库可共享文献，故不能当作独立验证。','small')

new('06  主背景结果：未见特异KEGG富集')
add('主问题：在所有可识别的处方人靶点中，便秘关键词交集是否具有相对偏重的功能？输入294个交集、背景331个处方靶点；每个功能体系只保留其有注释的有效背景，分母逐项留存。')
rows=[['功能体系','有效背景 / 交集','全部测试条目','BH显著条目']]
for fam,name in [('GO_BP','GO生物过程'),('GO_MF','GO分子功能'),('GO_CC','GO细胞组分'),('KEGG','KEGG')]:
 rr=bg[(bg.analysis.str.startswith('PRIMARY'))&(bg.family==fam)].iloc[0];dd=pr[pr.family==fam];rows.append([name,f'{rr.effective_universe_count} / {rr.effective_query_count}',len(dd),int((dd.p_adjust_BH_all_families<.05).sum())])
table(rows,[116,146,123,CW-385])
img('04_primary_enrichment.png',maxh=235)
add('唯一显著条目为GO:0032501“多细胞生物体过程”：交集174/294、背景180/330，富集倍数1.085，原始P=9.33×10^-7，统一BH q=0.00727。该条目范围很宽，不能据此声称识别了通便机制。')
add('KEGG共测试325条通路，显著数为0，最小统一BH q=0.9673。交集覆盖了处方靶点的大部分，剩余对照靶点少；主分析功能区分度有限。结果保持原样，不通过改变疾病分值或更换主背景制造显著性。')
add('采用单侧超几何检验；全部K≥1的条目均测试，k=0也保留；对每个分析的全部GO+KEGG条目统一BH校正，另存按功能体系BH值。GO沿is_a/part_of传播；未将regulates关系作为注释继承。','small')

new('07  全人背景：另一个问题的描述')
img('05_background_comparison.png',maxh=250)
add('敏感性分析将同一294个交集与已注释的人蛋白编码基因比较。HGNC提供19,297个有NCBI ID的蛋白编码基因；GO BP/MF/CC有效背景分别16,560/15,589/18,137，KEGG有效背景8,484。各项注释背景均完整留表。')
add('该分析的显著条目较多，但回答“这些数据库靶点与全人基因组有什么不同”，不能替代主分析的便秘条件区分，也不构成独立验证。以全部331个处方靶点作同样人背景对照，也得到大量显著条目，提示通用靶点选择和数据库注释偏倚的影响。')
ss=er[(er.analysis.str.startswith('SENSITIVITY'))&(er.family=='KEGG')].head(6)
rows=[['人背景敏感性：KEGG前6条','命中数','统一BH q']]
short={'hsa04933':'AGE-RAGE信号（糖尿病并发症）','hsa04082':'神经活性配体信号','hsa05200':'癌症通路','hsa05417':'脂质与动脉粥样硬化','hsa01100':'代谢通路','hsa05418':'流体剪切力与动脉粥样硬化'}
for rr in ss.itertuples():rows.append([rr.term_id+' '+short.get(rr.term_id,rr.term_name),rr.k_overlap,f'{rr.p_adjust_BH_all_families:.2e}'])
table(rows,[CW-150,60,90])
add('可将神经活性配体、钙信号及炎症相关候选配对作为后续验证方向，但不能只展示这些“好讲”的通路而省略癌症、代谢等同样显著的广泛信号。完整全部测试和校正结果见表40；这里的前6条仅按预设q值顺序展示。','small')

new('08  原方案执行状态与整体可行性')
table([['原方案环节','本次实际执行情况'],['TCMSP / HERB成分','完成公开药材检索、全量预定成分详情收集和范围审计；HERB 2,472/2,472条详情完成。'],['PubChem / 结构去重','实际下载CID/InChIKey性质，RDKit核对；明确结构进入主集，不为混合物编造SMILES。'],['SEA / SwissTargetPrediction','未运行新预测。SEA实际入口需验证码；Swiss条款禁止自动访问。没有绕过限制。'],['GeneCards疾病基因','采用用户手工导出CSV；保存原件和版权文字，全量纳入有效符号。HPO未进入主分析。'],['STRING PPI','实际调用官方API，导出真实组合分数及各通道；固定0.4阈值、不加节点。'],['Cytoscape网络','未运行软件；实际使用NetworkX/matplotlib，并提供可导入节点、边和GraphML。'],['GO / KEGG','真实注释、明确有效背景、全部测试及BH校正完成；主分析与人背景敏感性分别报告。']],[132,CW-132])
add('可行性判断','h2')
add('原方案第一步的数据整理、疾病交集、网络与富集部分已能产出可复算结果。能够提出有来源的具体成分-靶点候选；目前不能把整套流程描述为完成了SEA/Swiss预测，更不能由稠密网络或全人背景的显著通路推出整方疗效。')
add('集中缺口：瓜蒌仁种子记录不足，麦冬和铁皮石斛的靶点覆盖缺失；多糖、脂肪油整体及非明确单分子作用无法用当前小分子网络代表；原数据库存在物种/部位字段、衍生结构和关联证据质量问题；原始药材与本丸剂成分、含量、体内/肠腔暴露未实测；SEA/Swiss入口限制尚未解决。下一步应优先补齐制剂原形成分与关键配对的直接证据，再谈机制验证，不以追加网络分析补充这些缺口。','small')

new('09  数据版本、许可及交付索引')
add('所有网络请求在本机进程中禁用HTTP_PROXY/HTTPS_PROXY/ALL_PROXY，requests trust_env=False，NO_PROXY=*，TLS校验开启。未修改系统或Clash设置，未回退代理，未绕过登录、验证码或下载权限。','small')
table([['数据与版本','保存与用途'],['下载日期','2026-10-10 UTC；每次请求的时间、URL、状态和SHA256见logs/http_requests.jsonl。'],['GO / KEGG / STRING','GO releases/2026-07-26；KEGG pathway更新2026-10-09；STRING API 12.5。UniProt响应版本头、HGNC原文件和全人背景均留存。'],['计算环境','Python 3.12.7；pandas 2.2.2；SciPy 1.13.1；statsmodels 0.14.2；NetworkX 3.3；matplotlib 3.9.2；RDKit 2026.9.1；ReportLab 5.0.1。'],['许可备注','TCMSP页面声明ODbL；HERB公开访问不等同额外再分发许可。GeneCards导出版权声明要求非学术使用获得书面许可。PubChem、UniProt、GO、HGNC、STRING、KEGG各自许可/引用要求独立适用，公开接口不授予统一商业许可。原始材料保留用于内部复核，外发或商业交付前需核对相应授权。']],[128,CW-128])
add('先看这些文件','h2')
rows=[['文件 / 目录','内容'],['tables/37、32、34、36','逐味覆盖；唯一结构；完整成分-人靶点边；GeneCards交集'],['tables/50、51、54、55','完整药味-结构-靶点来源；去重组合；PPI节点回溯；6个讲解实例'],['tables/40-43、45-46','全部GO/KEGG测试、背景基因清单、显著结果；真实PPI边与节点'],['tables/52-53、results/*.graphml','Cytoscape可导入网络，不声称软件已运行'],['raw/、scripts/、logs/','原始响应与用户CSV；全部脚本；请求日志、环境、最终校验和'],['results/FINAL_README_中文.txt','结果解释、复算顺序、排除规则、方法修订和许可备注']]
table(rows,[192,CW-192])
add('来源入口','h2')
add('TCMSP: https://old.tcmsp-e.com/　HERB: http://herb.ac.cn/<br/>PubChem: https://pubchem.ncbi.nlm.nih.gov/　HGNC: https://www.genenames.org/<br/>UniProt: https://www.uniprot.org/　GO: https://geneontology.org/<br/>STRING: https://string-db.org/　KEGG: https://www.kegg.jp/<br/>香港中药材标准: https://www.cmro.gov.hk/html/eng/useful_information/hkcmms/volumes.html','small')
add('本稿为真实计算结果的内部审阅版本。全表和原始缓存是复核依据；未加入旧GSE结果，也未进行菌转化、宏基因组或MD扩展。','small')
def page(c,doc):
 c.saveState();c.setFillColor(teal);c.rect(0,H-9,W,9,fill=1,stroke=0);c.setFont('CN',8);c.setFillColor(muted);c.drawString(44,H-31,'改良通便丸 · Model A · 网络药理实际结果');c.drawRightString(W-44,H-31,'内部审阅初稿 | 2026-10-10');c.setStrokeColor(colors.HexColor('#D4DFDF'));c.line(44,39,W-44,39);c.drawString(44,25,'来源标签 ≠ 本项目实验证据；结果不构成疗效确认');c.drawRightString(W-44,25,str(doc.page));c.restoreState()
path=OUT/'改良通便丸_ModelA_网络药理结果_审阅初稿.pdf'
doc=SimpleDocTemplate(str(path),pagesize=(W,H),rightMargin=44,leftMargin=44,topMargin=57,bottomMargin=55,title='改良通便丸 Model A 网络药理实际结果',author='Reproducible analysis',allowSplitting=1)
doc.build(story,onFirstPage=page,onLaterPages=page);print(path)
