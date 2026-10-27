import pathlib,json,csv,datetime,sys
sys.path.append('/usr/lib/python3/dist-packages')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.platypus import SimpleDocTemplate,Paragraph,Table,TableStyle,PageBreak,Image
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib import colors
from xml.sax.saxutils import escape
N=pathlib.Path(__file__).parent;R=N.parent.parent;P=json.loads((N/'plan_frozen.json').read_text());c=json.loads((N/'FULL_SAMPLE_COMPLETE.json').read_text());q=json.loads((N/'RESOURCE_GATE_ACCEPTANCE.json').read_text());issue=json.loads((R/'00_protocol/LIMITATION_TB_COV_001.json').read_text())
with (N/'resource_samples.tsv').open() as f:obs=list(csv.DictReader(f,delimiter='\t'))
t0=datetime.datetime.fromisoformat(obs[0]['utc']);times=[(datetime.datetime.fromisoformat(x['utc'])-t0).total_seconds()/60 for x in obs]
fig,axes=plt.subplots(2,1,figsize=(9,6),sharex=True);axes[0].plot(times,[float(x['cgroup_memory_current'])/1024**3 for x in obs],label='Session memory current (includes cache)');axes[0].plot(times,[float(x['cgroup_memory_peak_including_cache'])/1024**3 for x in obs],label='Session memory peak',ls='--');axes[0].set_ylabel('GiB');axes[0].legend(fontsize=8);axes[1].plot(times,[float(x['MemAvailable_bytes'])/1024**3 for x in obs],color='#3b8b63',label='System MemAvailable');axes[1].set(xlabel='Elapsed minutes',ylabel='GiB');axes[1].legend(fontsize=8);fig.tight_layout();fig.savefig(N/'resource_trace.png',dpi=160);plt.close(fig)
pdfmetrics.registerFont(TTFont('CJK','/usr/share/fonts/truetype/arphic/uming.ttc'));body=ParagraphStyle('body',fontName='CJK',fontSize=10,leading=16,spaceAfter=9,wordWrap='CJK');head=ParagraphStyle('head',parent=body,fontSize=15,leading=22,spaceAfter=12);small=ParagraphStyle('small',parent=body,fontSize=8,leading=12);title=ParagraphStyle('title',parent=head,fontSize=19,leading=27);story=[]
def p(t,style=body):story.append(Paragraph(escape(t),style))
def table(rows,widths):
 t=Table(rows,colWidths=widths,repeatRows=1);t.setStyle(TableStyle([('FONTNAME',(0,0),(-1,-1),'CJK'),('FONTSIZE',(0,0),(-1,-1),9),('BACKGROUND',(0,0),(-1,0),colors.HexColor('#e9eff5')),('GRID',(0,0),(-1,-1),.4,colors.HexColor('#b7c4d1')),('TOPPADDING',(0,0),(-1,-1),8),('BOTTOMPADDING',(0,0),(-1,-1),8)]));story.append(t)
p('便秘项目缺陷记录与完整样本资源验收',title);p('内部阶段记录｜PRJNA1232352｜版本 fullsample_resource_gate_20261027_v1')
p('TB-COV-001 已登记并按用户要求暂缓处理',head)
p('缺陷：公开材料未提供与 60 个样本一对一匹配的年龄、性别和 BMI。样本编号与论文测序量已全部匹配，缺失的是个体协变量，并非样本身份不明。')
p('影响：当前 group-only 模型不能给出调整年龄、性别和 BMI 后的组别效应，也不能声称已排除这些混杂因素。该限制不阻断当前分组关联、功能候选整理或组装与注释。')
p('处理：暂停主动检索；保留既有结果及完整 P/BH；不用组均值填补个人记录，不合成或跨队列拼接协变量。获得可靠的逐样本个体表，或用户明确要求重新检索时再打开此项。届时在独立版本中补做调整模型，不覆盖原结果。')
p('完整样本组装、基因预测与 CAZyme 注释通过',head)
p(f"选择去宿主后读对数最多的 {c['run']}（{c['group']} 组），使用全部 {c['full_input_pairs']:,} 对读段。双端输入 SHA256 和序列数通过核查；三阶段退出码均为 0，关键输出已重新计算哈希并通过验收。")
a=c['assembly'];table([['指标','结果'],['完整输入读对',f"{c['full_input_pairs']:,}"],['contig 数（≥500 bp）',f"{a['contigs']:,}"],['组装总长度 / N50',f"{a['bases']:,} bp / {a['N50']:,} bp"],['预测蛋白数',f"{c['proteins']:,}"],['CAZyme 注释记录 / 双方法支持',f"{c['dbcan_any_method']:,} / {c['dbcan_two_methods']:,}"]],[230,235])
p('CAZyme 注释为家族层候选。双方法支持指 DIAMOND 与 dbCAN-HMM 均有支持，不等于已验证具体底物、反应方向或本方成分转化。',small)
story.append(PageBreak());p('Methods',head)
p('输入使用已验收的作者提交读段经本项目 fastp 与去宿主处理后得到的双端 FASTQ；使用完整文件，不做前缀抽样。样本选择仅依据保留读对数，与显著性结果无关。同一时间运行一个样本；MEGAHIT 和 dbCAN 使用 16 线程，Prodigal 按其原生命令单进程运行。MEGAHIT 1.2.9 内存参数为 64 GiB，k-mer 为 21,29,39,59,79,99,119,141，最短 contig 500 bp。Prodigal 2.6.3 使用 meta 模式预测基因。')
p('dbCAN 5.2.9 使用固定 CAZy-DIAMOND 和 dbCAN-HMM 参考，数据库大小及 SHA256 与技术试跑冻结记录一致。采用 DIAMOND E-value 1e-102，dbCAN-HMM E-value 1e-15、覆盖度 0.35；启用大输入处理模式，HMM 分批 5,000 条蛋白，16 线程。不运行 SignalP 或 DeepTMHMM。')
p('资源记录每 10 秒采集会话 cgroup memory.current、memory.peak、memory.events、系统 MemAvailable、内存和 I/O PSI 以及磁盘空闲量；每阶段保存 GNU time 的最大 RSS 与耗时。阶段启动前检查可用内存至少 100 GiB、磁盘空闲至少 300 GiB。FASTA 标识唯一性、基因与蛋白对应、注释记录与预测蛋白对应及输出 SHA256 均纳入验收。')
p('Results：完整样本资源实测',head)
rows=[['阶段','最大 RSS（GiB）','耗时（时:分:秒或分:秒）']]+[[x['stage'],f"{x['max_RSS_kbytes']/1024**2:.2f}",x['elapsed_wall']] for x in q['stages']];table(rows,[140,135,190])
p(f"专用 SSH 会话 cgroup 内存峰值（包含文件缓存）为 {q['cgroup_memory_peak_including_cache_GiB']:.2f} GiB；采样期间系统最低 MemAvailable 为 {q['sampled_min_MemAvailable_GiB']:.2f} GiB，最低磁盘空闲为 {q['sampled_min_disk_free_GiB']:.2f} GiB。本次输出目录占用 {q['disk_output_GiB']:.2f} GiB。会话 OOM/oom_kill 计数没有增加。")
p('GNU time 最大 RSS 与 cgroup 内存峰值采用不同口径；后者包括会话其他进程及文件缓存，不作为纯计算进程内存。MEGAHIT 的 64 GiB 参数不是系统级硬内存上限。',small)
story.append(PageBreak());p('资源轨迹与下一阶段边界',head);story.append(Image(str(N/'resource_trace.png'),width=465,height=310));p('图：上图显示专用会话 cgroup 当前与累计峰值内存，下图显示系统可用内存；横轴为启动后的分钟数。10 秒采样用于资源评估，不能保证捕获所有短时系统波动。',small)
p('本次验证完成了一个完整样本的流程和资源验收。读对数最多不保证组装复杂度最高，单样本成功不等于已证明多任务并行安全。本报告以首样本验收为范围。按用户最新授权，验收后以相同单任务配置接续其余 59 个样本，复用首样本结果；全量完成须另行核对队列完成凭据。')
p('后续保留 TB-COV-001 的解释边界。全队列组装、菌酶特异性参考核对、Model B、Model C 以及结构评价继续使用各自的进入和验收条件。')
p('记录：00_protocol/LIMITATION_TB_COV_001.json；本阶段 INPUT_ACCEPTANCE.json、ASSEMBLY_ACCEPTANCE.json、FULL_SAMPLE_COMPLETE.json、RESOURCE_GATE_ACCEPTANCE.json；逐阶段日志、退出码、资源轨迹和输出哈希。',small)
def footer(can,doc):can.setFont('CJK',8);can.drawString(42,25,'便秘项目 | TB-COV-001 与完整样本资源验收 | 内部记录');can.drawRightString(550,25,str(doc.page))
SimpleDocTemplate(str(N/'tongbian_fullsample_resource_acceptance_20261027.pdf'),pagesize=(595.28,841.89),leftMargin=55,rightMargin=55,topMargin=40,bottomMargin=48).build(story,onFirstPage=footer,onLaterPages=footer)
print('REPORT_CREATED')
