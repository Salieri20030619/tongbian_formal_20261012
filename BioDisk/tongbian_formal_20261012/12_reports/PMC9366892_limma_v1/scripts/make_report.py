from pathlib import Path
import json,hashlib,html,textwrap
import numpy as np,pandas as pd,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from reportlab.platypus import SimpleDocTemplate,Paragraph,Spacer,Table,TableStyle,PageBreak,Image
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib import colors
from reportlab.lib.pagesizes import A4
from reportlab.lib.units import mm
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.pdfgen import canvas
r=Path(__file__).resolve().parent
src=r/'server_results/06_human_selected_metabolites/limma_v1'
out=r/'report';fig=out/'figures';fig.mkdir(parents=True,exist_ok=True)
tab=pd.read_csv(src/'tables/primary20_group.tsv',sep='\t');adj=pd.read_csv(src/'tables/primary20_age_sex_BMI.tsv',sep='\t').set_index('feature_id').loc[tab.feature_id].reset_index()
ch=pd.read_csv(src/'tables/Cholic_acid_prespecified_sensitivity.tsv',sep='\t')
summary=pd.read_csv(src/'tables/analysis_summary.tsv',sep='\t')
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'savefig.facecolor':'white'})
blue='#246783';orange='#C87434';navy='#233E4F'
f,ax=plt.subplots(figsize=(9.4,11.2));y=np.arange(len(tab))[::-1]
for t,dy,col,label,marker in [(tab,.12,blue,'Group only','o'),(adj,-.12,orange,'Age / sex code / BMI adjusted','s')]:
 ax.errorbar(t.logFC,y+dy,xerr=np.vstack([t.logFC-t.CI95_low,t.CI95_high-t.logFC]),fmt=marker,color=col,markersize=4,elinewidth=1.1,capsize=2,label=label)
labels=['\n'.join(textwrap.wrap(z.Name,37))+f' [{z.feature_id}]' for z in tab.itertuples()]
ax.set_yticks(y,labels,fontsize=8.3);ax.axvline(0,color='#AAB8BF',linewidth=1);ax.set_axisbelow(True);ax.grid(axis='x',alpha=.18);ax.set_xlabel('STC - Control log2 peak area (logFC), moderated 95% CI')
ax.set_title('All 20 nonconflicting selected peaks',loc='left',fontweight='bold',pad=30,color=navy)
ax.legend(loc='lower left',bbox_to_anchor=(-.03,1.005),ncol=2,frameon=False,fontsize=9)
ax.set_ylim(-.8,19.8);f.subplots_adjust(left=.48,right=.98,top=.9,bottom=.065)
for ext in ['png','svg']:f.savefig(fig/f'Figure_01_main20_effects.{ext}',dpi=220)
plt.close(f)
orig=pd.read_csv(src/'audit/log2_peak_area_aligned.tsv',sep='\t').set_index('feature_id').loc['M407T218']
f,axes=plt.subplots(1,3,figsize=(10.4,4.0),sharey=True);rng=np.random.default_rng(20261012)
scenarios=[('Source A: 21-peak table','source_A'),('Source B: joint table C6','source_B'),('Exclude C6','exclude_C6')]
for ax,(title,scenario) in zip(axes,scenarios):
 v=orig.copy()
 if scenario=='source_B':v['C6']=np.log2(146617.187)
 if scenario=='exclude_C6':v=v.drop('C6')
 for pos,(prefix,col) in enumerate([('C',blue),('T',orange)],start=1):
  sub=v[[i.startswith(prefix) for i in v.index]];ax.boxplot(sub.to_numpy(),positions=[pos],widths=.45,showfliers=False,patch_artist=True,boxprops={'facecolor':col,'alpha':.16},medianprops={'color':col},whiskerprops={'color':col},capprops={'color':col})
  xx=np.full(len(sub),pos)+rng.uniform(-.14,.14,len(sub));ax.scatter(xx,sub,color=col,s=13,alpha=.75,zorder=3)
  if 'C6' in sub.index:
   loc=list(sub.index).index('C6');ax.scatter(xx[loc],sub['C6'],c='#B5353E',s=45,marker='D',zorder=4);ax.annotate('C6',(xx[loc],sub['C6']),xytext=(6,3),textcoords='offset points',color='#B5353E',fontsize=9)
 ax.set_xticks([1,2],[f'Control\nn={29 if scenario=="exclude_C6" else 30}','STC\nn=30']);ax.set_title(title,fontsize=10,fontweight='bold');ax.set_xlim(.5,2.5);ax.grid(axis='y',alpha=.15)
axes[0].set_ylabel('Cholic acid log2 peak area');f.tight_layout(pad=1.3)
for ext in ['png','svg']:f.savefig(fig/f'Figure_02_C6_source_sensitivity.{ext}',dpi=220)
plt.close(f)

pdfmetrics.registerFont(TTFont('CN',r'C:\Windows\Fonts\msyh.ttc'));pdfmetrics.registerFont(TTFont('CNB',r'C:\Windows\Fonts\msyhbd.ttc'))
pdfmetrics.registerFontFamily('CN',normal='CN',bold='CNB',italic='CN',boldItalic='CNB')
W,H=A4;cw=W-36*mm
base=ParagraphStyle('Body',fontName='CN',fontSize=10,leading=16,spaceAfter=7,textColor=colors.HexColor(navy),wordWrap='CJK')
small=ParagraphStyle('Small',parent=base,fontSize=8.5,leading=13,spaceAfter=5)
head=ParagraphStyle('Head',parent=base,fontName='CNB',fontSize=15,leading=22,spaceBefore=10,spaceAfter=9,keepWithNext=True)
title=ParagraphStyle('Title',parent=head,fontSize=23,leading=32,spaceAfter=12)
cell=ParagraphStyle('Cell',parent=small,fontSize=8,leading=11,spaceAfter=0)
def P(s,st=base):return Paragraph(html.escape(str(s)).replace('\n','<br/>'),st)
def grid(rows,widths,font=cell):
 t=Table([[P(x,font) for x in row] for row in rows],colWidths=widths,repeatRows=1,hAlign='LEFT')
 t.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,0),colors.HexColor('#E8F0F4')),('LINEBELOW',(0,0),(-1,0),.8,colors.HexColor('#9DB3C1')),('INNERGRID',(0,0),(-1,-1),.25,colors.HexColor('#D8E1E7')),('BOX',(0,0),(-1,-1),.4,colors.HexColor('#D8E1E7')),('VALIGN',(0,0),(-1,-1),'TOP'),('TOPPADDING',(0,0),(-1,-1),3),('BOTTOMPADDING',(0,0),(-1,-1),3)]));return t
def fnum(v):return f'{v:.4g}'
class Pages(canvas.Canvas):
 def __init__(self,*a,**kw):super().__init__(*a,**kw);self.states=[]
 def showPage(self):self.states.append(dict(self.__dict__));self._startPage()
 def save(self):
  n=len(self.states)
  for s in self.states:
   self.__dict__.update(s);self.setStrokeColor(colors.HexColor('#CAD7DE'));self.line(18*mm,16*mm,W-18*mm,16*mm);self.setFont('CN',8);self.setFillColor(colors.HexColor(navy));self.drawString(18*mm,11*mm,'PMC9366892 / limma v1 / 2026-10-12');self.drawRightString(W-18*mm,11*mm,f'{self._pageNumber} / {n}');super().showPage()
  super().save()
story=[P('改良通便丸研究',title),P('人体预选代谢物 limma 重分析',head),P('30 STC 与 30 健康对照。研究为同一队列的候选特征重分析，不是本方干预证据、独立验证或全代谢组无偏发现。'),P('20 个无冲突峰：16 个达到 P<0.05，BH<0.05 同为 16 个；其中 STC 升高 1 个、降低 15 个。调整年龄、Sex 原编码及 BMI 后，20/20 的方向与 P<0.05 状态一致。')]
story += [grid([['模型 / 结果范围','对照 / STC','P<0.05','BH<0.05'],['主模型：20 无冲突峰','30 / 30','16 / 20','16 / 20'],['协变量敏感性：20 峰','30 / 30','16 / 20','16 / 20'],['原来源全 21 峰记录*','30 / 30','17 / 21','17 / 21']],[cw*.49,cw*.19,cw*.16,cw*.16]),P('*全21峰表保留胆酸并标注C6冲突；该峰不进入主20峰汇总。主20峰单独拟合，BH按20项校正；全21峰各模型的BH按21项校正。各模型的P/BH不得混用。',small),P('输入与模型',head)]
for text in [
 '输入ZIP及11个包内文件哈希均通过；同一60个唯一sample_id按ID联结。原21峰面积范围362.741-1,154,395,322.514，无零值或缺失；直接log2，无伪计数、再归一化、插补或额外FC阈值。',
 'R 4.5.2 / limma 3.66.0；lmFit + eBayes(trend=FALSE, robust=FALSE)。主模型~group，Control为参照；logFC=STC-Control的log2峰面积差，2^logFC为几何均值比，不等于原文算术均值FC。',
 '敏感性模型~group+Age+SexCode+BMI。Sex按code1/code2分类处理，未臆测男女含义；未自动调整疾病时长、症状、Wexner、GIQLI或Bristol。60例设计秩2/5；去C6的59例设计秩仍为2/5。',
 '全部8个预设模型完成，保存166行结果、模型RDS、原值、脚本和日志。独立复算OLS效应、moderated t对应P、BH及95%区间通过，效应差最大1.87×10^-14。']:
 story.append(P(text,small))
story += [P('解释边界与待办',head),P('原始筛选发生在同一队列，当前P/BH仍受选择偏倚限制；本轮只重分析这21峰。两个Taurochenodeoxycholate峰ID保持分离；14个菌群条目含嵌套层级，不视为14个独立菌属。未运行偏相关，ppcor仍未安装；未跨RNA队列匹配个体。',small),P('来源：DOI 10.3389/fcimb.2022.956528；附件 DOI 10.3389/fcimb.2022.956528.s001。用户核验输入包：libfile_9025b7a56c008191971c58e1dd17241e v0。完整SHA-256见随附清单。',small),PageBreak()]

story += [P('主模型完整20峰结果',head),P('每峰n=30对照+30 STC；独立20峰拟合。按输入峰序展示全部结果，无P或FC筛选。正logFC为STC较高；P与BH为本轮新计算，未混用原文统计。',small)]
rows=[['峰ID / 原注释名称','logFC','P','BH']]
for z in tab.itertuples():rows.append([f'{z.feature_id}\n{z.Name}',f'{z.logFC:.4f}',fnum(z.P),fnum(z.BH)])
story += [grid(rows,[cw*.62,cw*.13,cw*.125,cw*.125]),P('完整精度、95%区间、样本数及方向见 tables/primary20_group.tsv；全部21峰原来源及协变量敏感性另表保留。',small),PageBreak()]
story += [P('20峰效应与协变量敏感性',head),Image(str(fig/'Figure_01_main20_effects.png'),width=cw,height=cw*11.2/9.4),P('图1：全部20个无冲突峰，蓝色为group主模型，橙色为年龄/SexCode/BMI敏感性；点为STC-Control的logFC，线为moderated 95%置信区间。每组n=30，排序与输入一致。该区间未作多重比较或预筛选择校正；显著性以结果表P/BH为准。',small),PageBreak()]
story += [P('胆酸 C6：预设三版本均保留',head),P('峰M407T218（Cholic acid）。C6在21峰表为22,146,617.1867411，在联合表为146,617.187。原输入均不改写；版本B仅在派生内存副本替换这一单元格，其他值与精度维持原21峰表。',small),Image(str(fig/'Figure_02_C6_source_sensitivity.png'),width=cw,height=cw*4/10.4),P('图2：log2原峰面积。红菱形为C6；各点为独立受试者。版本A/B为30对照+30 STC；去C6为29对照+30 STC。三版本预先规定，未根据结果挑选。',small)]
rows=[['胆酸版本 / 模型','C / STC','logFC','P','BH*']]
scenario_names={'source_A':'A：21峰表原值','source_B':'B：联合表C6原值','exclude_C6':'排除C6样本'}
for z in ch.itertuples():rows.append([scenario_names[z.scenario]+('\n年龄/SexCode/BMI调整' if z.design!='group' else '\ngroup主模型'),f'{z.n_Control} / {z.n_STC}',f'{z.logFC:.4f}',fnum(z.P),fnum(z.BH)])
story += [grid(rows,[cw*.38,cw*.15,cw*.15,cw*.16,cw*.16]),P('*BH取自对应全21峰拟合的21项检验家族；非仅对本表6项校正。全21峰敏感性模型亦完整保存。',small),P('六个模型均为STC较低，logFC为-2.513至-2.271，P为0.00211至0.00792；这一方向一致性不能消除源表数值冲突。胆酸不纳入主20峰计数，冲突状态继续保留。',small),P('本轮剩余阻塞：ppcor安装尚未获批准，偏相关未运行；完整未预筛代谢矩阵和完整菌群组成表仍未具备。不得由这些结果声称本方作用机制已获验证。',small)]
pdfdir=r.parent/'output/pdf';pdfdir.mkdir(exist_ok=True);pdf=pdfdir/'改良通便丸_PMC9366892_limma代谢物重分析_20261012.pdf'
SimpleDocTemplate(str(pdf),pagesize=A4,rightMargin=18*mm,leftMargin=18*mm,topMargin=18*mm,bottomMargin=22*mm,title='PMC9366892 预选代谢物limma重分析',author='Research analysis').build(story,canvasmaker=Pages)
notes={'PDF':str(pdf),'figures':[str(p) for p in fig.glob('*')],'sha256':hashlib.sha256(pdf.read_bytes()).hexdigest(),'bytes':pdf.stat().st_size}
(out/'report_generation_receipt.json').write_text(json.dumps(notes,ensure_ascii=False,indent=2),encoding='utf-8');print(json.dumps(notes,ensure_ascii=False))
