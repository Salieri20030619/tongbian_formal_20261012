from pathlib import Path
import json,hashlib,shutil,sys
import pandas as pd,numpy as np,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from reportlab.pdfgen import canvas
from reportlab.platypus import SimpleDocTemplate,Paragraph,Spacer,Table,TableStyle,Image,PageBreak
from reportlab.lib.styles import getSampleStyleSheet,ParagraphStyle
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.lib import colors
from reportlab.lib.enums import TA_LEFT
import pymupdf
P=Path(r'C:\Users\31937\Downloads\改良通便丸_正式分析_20261012');O=P/'10_partial_spearman_20261013/partial_spearman_v1';F=O/'figures';F.mkdir(exist_ok=True)
d=pd.read_csv(O/'tables/ALL_465_results.tsv',sep='\t');cfg=json.loads((O/'analysis_config_FROZEN.json').read_text());taxa=cfg['taxa'];short=['Bacteroides','Parabacteroides','Desulfovibrionaceae (f)','Ruminiclostridium 5','Subdoligranulum']
labels=dict(zip(taxa,short));ann=pd.read_csv(O/'inputs/metabolites_21_annotation.csv');name=dict(zip(ann.ID,ann.Name))
def plot_name(x):
 t=name[x].replace(' (peroxide free)','');return t if len(t)<=27 else t[:24]+'...'
def heat(family,models,filename):
 features=cfg[family];fig,axs=plt.subplots(1,len(models),figsize=(max(7,5*len(models)),5 if len(features)==8 else 9),squeeze=False)
 for ax,model in zip(axs[0],models):
  t=d[(d.family==family)&(d.model==model)];v=t.pivot(index='feature_id',columns='taxon',values='rho').reindex(index=features,columns=taxa);p=t.pivot(index='feature_id',columns='taxon',values='P').reindex_like(v);q=t.pivot(index='feature_id',columns='taxon',values='BH').reindex_like(v)
  im=ax.imshow(v,vmin=-.6,vmax=.6,cmap='RdBu_r',aspect='auto');ax.set_xticks(range(5),short,rotation=40,ha='right',fontsize=8);ax.set_yticks(range(len(features)),[f'{x} | {plot_name(x)}' for x in features],fontsize=8 if len(models)==1 else 7);ax.set_title({'group':'Group adjusted (primary)','unadjusted':'Unadjusted','group_age_sex_BMI':'Group + age + sex + BMI'}[model],fontsize=11)
  for i in range(len(features)):
   for j in range(5):ax.text(j,i,f'{v.iloc[i,j]:.2f}'+('**' if q.iloc[i,j]<.05 else '*' if p.iloc[i,j]<.05 else ''),ha='center',va='center',fontsize=8 if len(models)==1 else 7,color='white' if abs(v.iloc[i,j])>.43 else 'black')
 fig.tight_layout();fig.savefig(F/filename,dpi=220,bbox_inches='tight');plt.close(fig)
heat('primary40',['unadjusted','group','group_age_sex_BMI'],'Figure1_primary40_models.png');heat('primary40',['group'],'Figure1_primary40_group.png');heat('extended100',['group'],'Figure2_extended100_group.png')
fig,ax=plt.subplots(figsize=(8,3.3))
for sc,off,col in [('source_A',-.16,'#246781'),('source_B',0,'#ba6943'),('exclude_C6',.16,'#688d57')]:
 t=d[(d.family=='cholic5')&(d.model=='group')&(d.scenario==sc)].set_index('taxon').loc[taxa];ax.scatter(t.rho,np.arange(5)+off,label=sc,color=col,s=36)
ax.axvline(0,color='gray',lw=.7);ax.set_yticks(range(5),short);ax.set_xlabel('Partial Spearman rho (group adjusted)');ax.set_xlim(-.5,.05);ax.invert_yaxis();ax.legend(fontsize=8,loc='lower right');fig.tight_layout();fig.savefig(F/'Figure3_cholic_sensitivity.png',dpi=220);plt.close(fig)
pdfmetrics.registerFont(TTFont('CN',r'C:\Windows\Fonts\simsun.ttc',subfontIndex=0));styles=getSampleStyleSheet()
for key in ['Normal','BodyText','Title','Heading1','Heading2']:styles[key].fontName='CN'
styles.add(ParagraphStyle(name='BodyCN',fontName='CN',fontSize=10.3,leading=16,spaceAfter=8,wordWrap='CJK'))
styles.add(ParagraphStyle(name='SmallCN',fontName='CN',fontSize=8.1,leading=11.5,spaceAfter=6,wordWrap='CJK'))
styles['Title'].fontSize=23;styles['Title'].leading=31;styles['Heading1'].fontSize=16;styles['Heading1'].leading=22
S=[]
def p(t,style='BodyCN'):S.append(Paragraph(t.replace('−','-'),styles[style]))
def h(t):p(t,'Heading1')
def table(rows,widths,size=8.3):
 data=[[Paragraph(str(x),ParagraphStyle(name='c',fontName='CN',fontSize=size,leading=size+3,wordWrap='CJK')) for x in r] for r in rows];t=Table(data,colWidths=widths,repeatRows=1,hAlign='LEFT');t.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,0),colors.HexColor('#dfecef')),('VALIGN',(0,0),(-1,-1),'TOP'),('LINEBELOW',(0,0),(-1,0),.7,colors.HexColor('#427481')),('ROWBACKGROUNDS',(0,1),(-1,-1),[colors.white,colors.HexColor('#f5f7f7')]),('LEFTPADDING',(0,0),(-1,-1),5),('RIGHTPADDING',(0,0),(-1,-1),5),('TOPPADDING',(0,0),(-1,-1),3),('BOTTOMPADDING',(0,0),(-1,-1),3)]));S.append(t);S.append(Spacer(1,8))
def pic(n,w,maxh):
 from PIL import Image as PI
 z=PI.open(F/n);iw,ih=z.size;scale=min(w/iw,maxh/ih);S.append(Image(str(F/n),iw*scale,ih*scale));S.append(Spacer(1,6))
def fmt(x):return f'{x:.3g}'
p('改良通便丸公开多组学研究','Title');p('人群菌群—代谢物偏相关<br/>正式计算、敏感性分析与核验报告','Title');p('PMC9366892 | 60例，C/T各30 | 2026-10-13')
h('调整疾病分组后，主集合保留15对名义显著关联')
p('冻结的5个候选菌分类特征与8个无原值冲突代谢峰形成40对主检验。调整group后，15对P＜0.05，其中12对BH＜0.05；14对为负相关、1对为正相关。进一步调整年龄、性别编码和BMI后，14对保持名义显著且方向一致，11对通过BH＜0.05，40对均未发生方向翻转。')
p('主要负相关集中于牛磺酸、花生四烯酸与Desulfovibrionaceae、Ruminiclostridium 5及Subdoligranulum等候选特征。最强主关联为Desulfovibrionaceae—牛磺酸：rho=-0.514，P=3.15×10<super>-5</super>，BH=0.00126。')
table([['集合','未调整 P/BH','调整group P/BH','完整协变量 P/BH'],['主集合40对','25 / 17','15 / 12','14 / 11'],['扩展集合100对','44 / 23','19 / 5','19 / 9']],[105,120,140,135])
p('表1．计数依次为P＜0.05和BH＜0.05的检验数。主、扩展集合重叠，分别计算BH；扩展集合不构成独立验证。完整协变量为group、Age、SexCode和BMI。','SmallCN')
p('研究仍以20味全方为背景：白术、黄芪、肉苁蓉、当归、陈皮、瓜蒌仁、苦杏仁、火麻仁、枳实、槟榔、厚朴、桔梗、麦冬、玄参、生地黄、山药、铁皮石斛、黄精、党参、枸杞。本批为公开疾病队列的关联证据，供Model C分层优选。')
S.append(PageBreak());h('Methods：按样本标识连接并预先冻结检验族')
p('输入为现有去标识包的metadata、joint_selected、21峰原始峰面积、14个候选分类特征相对丰度、注释和冲突表。所有受试者矩阵先核查60个唯一sample_id与同一ID集合，再显式匹配。保留完整精度的原始峰面积和菌群数值；不归一化、不填补、不按行序拼接。Spearman使用平均秩处理并列值。')
p('joint_selected定义5个候选菌和9个代谢物。排除存在原值冲突的Cholic acid后，主集合为5×8=40对；扩展集合为相同5菌×20个无冲突峰=100对。名称Arachidonic Acid与注释Arachidonic Acid (peroxide free)显式映射至M303T39，核查跨表数值舍入误差及完全一致秩次后使用。全部集合在计算之前写入冻结配置。')
p('主模型采用ppcor::pcor.test(method="spearman")调整group；敏感性模型增加Age、SexCode和BMI。group和SexCode使用K−1哑变量，性别1/2仅保留编码。未调整结果采用cor.test(method="spearman", exact=FALSE)。均为双侧检验。')
p('每项检验报告n、调整变量列数k、设计秩、rho、统计量、自由度n−k−2、P及BH。检查原始设计和秩设计满秩、n＞k+2及联合变量矩阵可逆。主模型n=60、k=1、df=57；完整协变量n=60、k=4、df=54。排除C6的胆酸情景n=59，对应df=56或53。')
p('BH在每个预设“检验族×模型×胆酸情景”内分别计算。以P＜0.05筛选名义显著候选，保留BH；不附加效应量阈值。465行输出包括主集合120行、扩展集合300行和胆酸45行。')
p('胆酸独立使用5菌×1峰的检验族，分别计算来源A（C6=22146617.1867411）、来源B（C6=146617.187）及排除C6三种情景，每种情景保留全部三个模型。来源文件不作修改。')
h('软件与可复核性')
p('服务器osv01上R 4.5.2、ppcor 1.1；源代码包来自https://cloud.r-project.org/src/contrib/ppcor_1.1.tar.gz，安装至项目专用库。安装日志、归档包SHA256、sessionInfo、函数源码、冻结配置、输入哈希及运行PID/退出码均随附。465项结果逐项通过独立秩残差rho与n−k−2自由度P值核验；调整模型通过行序倒置检查。')
S.append(PageBreak());h('主关联在协变量敏感性分析中保持方向')
pic('Figure1_primary40_group.png',505,385)
p('图1．主集合全部40对调整group后的rho。蓝色为负相关，红色为正相关；*表示仅名义P＜0.05，**表示BH＜0.05，BH家族大小40。n=60、k=1、df=57。三个模型的完整对照另附Figure1_primary40_models.png。主模型固定为调整group，不选择更显著的模型。Desulfovibrionaceae为科级，其他4个为属级分类标签。','SmallCN')
p('从未调整的25对名义显著减少至调整group后的15对，表明部分总体关联随分组调整减弱。完整协变量分析保留14对同方向名义显著关联，支持这些候选在所检查协变量范围内具有一致方向。')
p('Subdoligranulum—4-Pyridoxic acid为主集合唯一名义显著正相关（rho=0.300，P=0.0211，BH=0.0650），其证据定位为探索性信号。其余主集合名义显著关联为负相关。')
S.append(PageBreak());h('主模型40对完整统计：不省略阴性结果')
main=d[(d.family=='primary40')&(d.model=='group')].sort_values(['taxon','feature_id'])
abbr={taxa[0]:'Bacteroides',taxa[1]:'Parabacteroides',taxa[2]:'Desulfovibrionaceae',taxa[3]:'Ruminiclostridium 5',taxa[4]:'Subdoligranulum'}
rows=[['候选菌分类特征','代谢峰ID','rho','P','BH']]+[[abbr[r.taxon],r.feature_id,f'{r.rho:.3f}',fmt(r.P),fmt(r.BH)] for r in main.itertuples()]
table(rows,[175,105,65,77,78],size=7.3)
p('表2．全部n=60，k=1，df=57；BH家族大小40。表内数值为显示舍入，完整精度及样本ID列表见TSV。代谢峰名称与注释完整保留在输入注释表及ALL_465_results.tsv。','SmallCN')
S.append(PageBreak());h('扩展集合提供补充候选，保留独立BH口径')
pic('Figure2_extended100_group.png',490,575)
p('图2．预设5菌×20无冲突代谢峰的调整group结果，n=60、k=1、df=57，BH在100对内计算。*表示仅名义P＜0.05，**表示BH＜0.05。19对名义显著、5对BH显著。与主集合重叠的关联rho/P完全相同，BH因检验族不同而变化。','SmallCN')
S.append(PageBreak());h('胆酸冲突改变部分阈值判断，原值继续隔离')
pic('Figure3_cholic_sensitivity.png',490,215)
p('图3．胆酸三个预设情景的调整group偏相关。三种情景的5个方向均为负；来源A/B使用60例，排除C6使用59例。每个情景单独在5对内BH。图仅显示rho，不表示原值冲突已解决。','SmallCN')
c=d[(d.family=='cholic5')&(d.model=='group')];rows=[['情景','菌分类特征','rho','P','BH']]+[[r.scenario,abbr[r.taxon],f'{r.rho:.3f}',fmt(r.P),fmt(r.BH)] for r in c.itertuples()];table(rows,[75,180,65,90,90],7.8)
p('表3．Parabacteroides—胆酸在来源A下P=0.0543、来源B下P=0.0394、排除C6后P=0.0503；该关联对冲突处理敏感，不能由其单独决定作用链优先级。其余完整协变量与未调整情景均在完整结果表中，不择优报告。','SmallCN')
S.append(PageBreak());h('Discussion：用于证据分层，不作转化或药效因果解释')
p('这些结果显示所选候选菌与代谢峰在疾病分组调整后的关联。负相关本身不能判定某菌消耗该底物，也不能证明其酶催化特定生成或消耗反应。Model B仍需底物、反应方向、实验证据酶及宏基因组同源候选共同支持。')
p('本研究使用原论文在同一队列预筛的代谢峰和菌特征，BH不能消除前期特征筛选带来的偏倚。候选分类相对丰度不能替代完整微生物组成，也没有校正完整组成偏差。14个嵌套分类特征未被当作14个独立菌属；主分析仅使用预设5个标签。')
p('性别编码缺少男女释义，但分类调整无需猜测其方向。调整group并非组内因果识别；未测量混杂及样本选择仍限制外推。未将PMC9366892与GSE245885进行个体层面拼接。两个Taurochenodeoxycholate峰保持独立峰ID，主分析只使用21峰表的M498T146，不合并7峰表中的M500T143。')
p('本轮复用Model A的2059个结构、331个人源靶点及294个便秘交集靶点；历史来源靶点未改称新预测。既有人RNA与代谢limma未重新拟合。整体路径仍为药材原形成分→菌群转化能力与产物→宿主靶点和功能→便秘表型，Model C仅作证据分层优选。')
h('交付与后续')
p('本批完整交付位于10_partial_spearman_20261013/partial_spearman_v1，含tables、audit、scripts、inputs、figures及本报告。计算已完成；整个多组学项目尚未完成。后续队列独立性、动物输入与大数据下载核查另列接续报告，并持续更新项目agent.md、HANDOFF.md和结构化状态。')
h('来源')
p('1. 人群原研究：DOI 10.3389/fcimb.2022.956528；PMC9366892。输入包README及SHA256SUMS保留实际来源与去标识记录。<br/>2. ppcor官方手册：https://cran.r-project.org/web/packages/ppcor/ppcor.pdf。<br/>3. 本地完整研究方案及方法进入条件流程图，2026-10-12版。<br/>4. 本批运行与核验记录：COMPUTE_COMPLETE.json、analysis_config_FROZEN.json、audit/independent_formula_and_order_checks.tsv。','SmallCN')
pdf=O/'改良通便丸_人群菌群代谢物偏相关_完整报告_20261013.pdf'
def footer(c,doc):
 c.setFont('CN',8);c.setFillColor(colors.HexColor('#567078'));c.drawString(43,25,'改良通便丸公开多组学研究 | 选定特征关联再分析 | 2026-10-13');c.drawRightString(553,25,str(doc.page))
SimpleDocTemplate(str(pdf),pagesize=(595.28,841.89),rightMargin=43,leftMargin=43,topMargin=42,bottomMargin=43).build(S,onFirstPage=footer,onLaterPages=footer)
qa=O/'audit/pdf_QA';qa.mkdir(exist_ok=True);doc=pymupdf.open(pdf)
for i,page in enumerate(doc):page.get_pixmap(matrix=pymupdf.Matrix(1.25,1.25)).save(qa/f'page_{i+1:02}.png')
(qa/'text.txt').write_text('\n'.join(p.get_text() for p in doc),encoding='utf-8');print('PDF',pdf,'PAGES',len(doc));shutil.copy2(__file__,O/'scripts'/Path(__file__).name)
