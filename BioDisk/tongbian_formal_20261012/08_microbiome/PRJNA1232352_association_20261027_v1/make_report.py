import pathlib,json,csv,hashlib,datetime,math
import sys
sys.path.append("/usr/lib/python3/dist-packages")
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.platypus import SimpleDocTemplate,Paragraph,Spacer,Table,TableStyle,Image,PageBreak
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib import colors
from xml.sax.saxutils import escape
O=pathlib.Path(__file__).parent; R=O.parent.parent
assert (O/'statistics_exit_code.txt').read_text().strip()=='0' and (O/'STATISTICS_VALIDATED.txt').exists()
def read(n):
 with (O/n).open(encoding='utf-8') as f: return list(csv.DictReader(f,delimiter='\t'))
def write(n,h,rows):
 with (O/n).open('w',encoding='utf-8',newline='') as f:
  w=csv.writer(f,delimiter='\t',lineterminator='\n'); w.writerow(h); w.writerows(rows)
meta=read('metadata.tsv'); ids=[x['sample_id'] for x in meta]; group={x['sample_id']:x['group'] for x in meta}; res=read('association_all.tsv'); summary=read('association_summary.tsv')
assert len(res)==2859 and len({(x['feature_id'],x['contrast']) for x in res})==len(res)
# Independent BH recomputation within every declared family; no pooling of classes.
for fam in {x['test_family'] for x in res}:
 rows=sorted([x for x in res if x['test_family']==fam],key=lambda x:float(x['pval'])); m=len(rows); last=1.
 for i in range(m-1,-1,-1):
  last=min(last,float(rows[i]['pval'])*m/(i+1)); assert abs(last-float(rows[i]['BH']))<1e-12
stats=[]; matrices={}
for cl,filename in [('species','species_percent.tsv'),('pathway','pathabundance_cpm.tsv')]:
 mat={x['feature']:[float(x[k]) for k in ids] for x in read(filename)}; matrices[cl]=mat
 for f,vals in mat.items():
  for g in ['H','DCTT','NCTT']:
   a=np.array([v for run,v in zip(ids,vals) if group[run]==g]); stats.append([cl,f,g,len(a),int(sum(a>0)),float(a.mean()),float(np.median(a)),float(a.std(ddof=1))])
write('feature_group_descriptives.tsv',['class','feature','group','n','positive_n','mean','median','sd'],stats)
figdir=O/'figures'; figdir.mkdir(exist_ok=True)
plt.rcParams.update({'font.size':10,'savefig.dpi':160})
# Sample-level composition preserves all samples and shows pooled top taxa without using P values.
mat=matrices['species']; top=sorted(mat,key=lambda k:-np.mean(mat[k]))[:12]; order=sorted(range(60),key=lambda i:(['H','DCTT','NCTT'].index(group[ids[i]]),meta[i]['sample_alias']))
fig,ax=plt.subplots(figsize=(12,5.4)); bottom=np.zeros(60); palette=plt.cm.tab20.colors
for j,f in enumerate(top):
 v=np.array(mat[f])[order]; ax.bar(range(60),v,bottom=bottom,width=.95,label=f.split('|')[-1].replace('s__',''),color=palette[j]); bottom+=v
other=np.array([sum(v[i] for k,v in mat.items() if k not in top) for i in order]); ax.bar(range(60),other,bottom=bottom,color='#dddddd',label='Other species')
ax.set(xlim=(-1,60),ylim=(0,101),ylabel='Relative abundance (%)',title='Species composition: all 60 samples'); ax.set_xticks([9.5,29.5,49.5],['H (n=20)','DCTT (n=20)','NCTT (n=20)']); ax.axvline(19.5,color='white'); ax.axvline(39.5,color='white'); ax.legend(loc='upper left',bbox_to_anchor=(1.01,1),fontsize=8,frameon=False); fig.tight_layout(); fig.savefig(figdir/'species_composition.png'); plt.close(fig)
plotfiles=[]
for cl in ['species','pathway']:
 for contrast in ['DCTT-H','NCTT-H','DCTT-NCTT']:
  z=sorted([x for x in res if x['class']==cl and x['contrast']==contrast],key=lambda x:(float(x['pval']),x['feature']))[:12]
  labels=[x['feature'].split('|')[-1].replace('s__','') if cl=='species' else x['feature'].split(':')[0] for x in z]
  beta=np.array([float(x['coef']) for x in z]); lo=np.array([float(x['CI95_low']) for x in z]); hi=np.array([float(x['CI95_high']) for x in z]); color=['#b33c32' if float(x['BH'])<.05 else '#28759a' if float(x['pval'])<.05 else '#888888' for x in z]
  fig,ax=plt.subplots(figsize=(11,5.5)); ax.axvline(0,color='#aaaaaa',lw=.8)
  for i in range(len(z)): ax.errorbar(beta[i],i,xerr=[[beta[i]-lo[i]],[hi[i]-beta[i]]],fmt='o',color=color[i],capsize=3)
  ax.set_yticks(range(len(z)),labels); ax.invert_yaxis(); ax.set(xlabel='Coefficient: difference in transformed log2 relative abundance (95% CI)',title=f'{cl.capitalize()}: {contrast} | top 12 by raw P')
  ax.text(.99,.01,'Red: BH < 0.05 | Blue: raw P < 0.05 only | Gray: raw P >= 0.05',transform=ax.transAxes,ha='right',fontsize=8); fig.tight_layout(); name=cl+'_'+contrast+'.png'; fig.savefig(figdir/name); plt.close(fig); plotfiles.append((cl,contrast,name))
  write('top12_'+cl+'_'+contrast+'.tsv',list(z[0]),[[x[k] for k in z[0]] for x in z])
pdfmetrics.registerFont(TTFont('CJK','/usr/share/fonts/truetype/arphic/uming.ttc'))
style=ParagraphStyle('body',fontName='CJK',fontSize=10,leading=16,spaceAfter=7,wordWrap='CJK')
title=ParagraphStyle('title',parent=style,fontSize=20,leading=27,spaceAfter=18)
head=ParagraphStyle('heading',parent=style,fontSize=14,leading=21,spaceBefore=12,spaceAfter=10)
small=ParagraphStyle('small',parent=style,fontSize=8,leading=12)
story=[]
def para(t,sty=style): story.append(Paragraph(escape(t),sty))
para('便秘项目：宏基因组验收与关联分析',title)
para('内部阶段报告｜PRJNA1232352｜版本 association_20261027_v1')
para('60 例全量结果通过验收，物种与通路关联分析已完成',head)
para('DCTT、NCTT、H 各 20 个独立样本标识。本轮完成 240 个 profiling 输出的 SHA256 复核，全部一致；逐样本退出码、输入接收凭据、读数及冻结脚本和配置通过核查。物种矩阵含 1,233 个物种，通路矩阵含 518 条通路；另保留 22,009 条物种分层通路记录。')
para('组间比较检出不同强度的物种和通路关联',head)
header=['类别 / 对比','检验数','原始 P<0.05','BH<0.05']
rows=[header]+[[x['family'].replace('species:','物种 ').replace('pathway:','通路 '),x['tested'],x['nominal'],x['BH_significant']] for x in summary]
t=Table(rows,colWidths=[210,65,105,75],repeatRows=1); t.setStyle(TableStyle([('FONTNAME',(0,0),(-1,-1),'CJK'),('FONTSIZE',(0,0),(-1,-1),9),('BACKGROUND',(0,0),(-1,0),colors.HexColor('#e7eef5')),('GRID',(0,0),(-1,-1),.4,colors.HexColor('#c9d2dc')),('TOPPADDING',(0,0),(-1,-1),8),('BOTTOMPADDING',(0,0),(-1,-1),8)])); story.append(t)
para('名义显著计数包含 BH 显著结果。物种 DCTT–H 和 DCTT–NCTT 分别有 10 和 1 个 BH 显著关联；通路 NCTT–H 有 142 个 BH 显著关联。其余对比保留原始 P<0.05 的探索性信号，不表述为校正后显著。')
para('本轮范围：全量验收、矩阵整理、MaAsLin2 关联及效应图。整个项目仍未完成；菌酶转化、Model B、Model C 与结构评价需按各自进入条件继续。',small)
story.append(PageBreak())
para('Methods',head)
para('输入为已完成的 MetaPhlAn 4.1.1 物种谱和 HUMAnN 3.9 通路表，采用同一冻结参考库版本。按 SRR 编号连接 60 份输出与原始 run manifest，核对唯一 BioSample、唯一样本别名和组别；归档全部 60 份 ENA BioSample XML。物种百分比保留原始值；HUMAnN 非分层的具名通路丰度以样本内总和归一化为 copies per million，UNMAPPED/UNINTEGRATED 单独保留；物种分层通路采用相同群落总和作为分母；覆盖度不归一化。')
para('在检验前冻结筛选规则：在全部 60 个样本中至少 6 个样本丰度大于 0，且跨样本方差大于 0。486 个物种和 467 条通路进入分析，未按 P 值或组别预筛选。模型使用 MaAsLin2 1.22.0，固定效应为 group；全部 60 个样本参与三个预设对比，残差自由度为 57。采用 normalization=NONE；LOG 变换按每个特征将零值替换为该特征最小正值的一半，再取 log2。')
para('以 H 和 NCTT 分别作为参照重拟合，提取 DCTT–H、NCTT–H、DCTT–NCTT 的目标系数、标准误、双侧 P 值和 95% 置信区间。按“特征类别×对比”分别实施 BH 校正；以原始 P<0.05 筛选名义显著，并报告 BH<0.05。每组样本量为 20，模型总样本量为 60。分析前保存参数文件；每类先用 10 个特征进行技术试运行，之后全部系数、标准误、P 值与独立 stats::lm 重算核对。')
para('Results：物种组成保留显著的样本间差异',head)
story.append(Image(str(figdir/'species_composition.png'),width=475,height=214))
para('图 1：60 个样本的物种组成。按全部样本平均丰度选择前 12 个物种，其余合并为 Other species；每根柱为一个样本，组内按样本别名排序。该图为描述性结果。',small)
for cl,contrast,name in plotfiles:
 story.append(PageBreak()); para(('物种' if cl=='species' else '通路')+'关联：'+contrast,head)
 row=next(x for x in summary if x['family']==cl+':'+contrast)
 para(f"共检验 {row['tested']} 个特征，原始 P<0.05 为 {row['nominal']} 个，其中 BH<0.05 为 {row['BH_significant']} 个。正系数表示对比左侧组的变换后相对丰度较高。")
 story.append(Image(str(figdir/name),width=475,height=238))
 para('效应图：按原始 P 值排序展示前 12 个特征，点为系数，线为 95% 置信区间。红色为 BH<0.05，蓝色仅原始 P<0.05，灰色为原始 P≥0.05。物种图保留物种名称，通路图展示稳定通路编号；完整名称、效应、P 值、BH、分组均值、中位数与检出数保存在对应统计表。',small)
 top=sorted([x for x in res if x['class']==cl and x['contrast']==contrast],key=lambda x:(float(x['pval']),x['feature']))[:3]
 for x in top:
  name=x['feature'].split('|')[-1].replace('s__','')
  para(f"{name}：系数 {float(x['coef']):.3f}；95% CI [{float(x['CI95_low']):.3f}, {float(x['CI95_high']):.3f}]；P={float(x['pval']):.3g}；BH={float(x['BH']):.3g}。",small)
story.append(PageBreak()); para('Discussion：证据边界与后续工作',head)
para('本结果确认当前队列中的组别关联，不代表改良通便丸药效或菌群变化的因果作用。公开 BioSample 属性未提供年龄、性别或 BMI，本轮为组别模型；尚未调整这些潜在混杂因素。独立样本标识来自 60 个唯一 BioSample 及样本别名，未进行跨队列个体拼接。')
para('系数表示经零值处理后的 log2 相对丰度差异，不是绝对菌量变化，也不是原始组均值的倍数差。具名通路采用排除特殊行后的样本内总和归一化；其组成性和分母选择应在后续机制解释时保持一致。多重检验按六个预设家族分别控制；仅原始 P<0.05 的结果作为探索性候选。')
para('本轮未重读全部原始 FASTQ 重算上游哈希，而是复核已冻结输入凭据、样本读数和本阶段实际输出。尚未启动组装/菌酶功能注释、反应规则、Model B/C 或结构评价；这些分支按既定方案独立验收。')
para('下一步：结合完整统计表及物种分层通路贡献筛选功能候选；核验 MEGAHIT/Prodigal/dbCAN 分支的资源与输入条件，先行技术试跑；反应转化候选必须有可追溯的底物、酶及反应依据。')
para('复核记录',head)
para('PROFILING_ACCEPTANCE.json：60 样本、240 输出哈希、阶段退出码及输入读数核查。MATRIX_ACCEPTANCE.json：矩阵全单元格保存后回读一致。PILOT_PASS.txt 与 STATISTICS_VALIDATED.txt：试运行及全量独立线性模型复核。association_all.tsv：2,859 项检验，完整保留目标系数、P 值、BH 与检验家族。')
para('方法来源：项目冻结研究方案；https://github.com/biobakery/Maaslin2；https://github.com/biobakery/humann。公开样本来源：https://www.ebi.ac.uk/ena/browser/，本阶段按 BioSample 编号保存 XML 快照。',small)
report=O/'tongbian_profiling_association_internal_20261027.pdf'
def footer(c,d):
 c.setFont('CJK',8); c.setFillColor(colors.HexColor('#667788')); c.drawString(42,25,'PRJNA1232352 | association_20261027_v1 | 内部阶段报告'); c.drawRightString(550,25,str(d.page))
SimpleDocTemplate(str(report),pagesize=(595.28,841.89),rightMargin=55,leftMargin=55,topMargin=40,bottomMargin=48).build(story,onFirstPage=footer,onLaterPages=footer)
print(json.dumps({'statistics_tests':len(res),'figures':7,'report':str(report),'BH_independent_recheck':'PASS'}))
