"""Export auditable Cytoscape tables and real result figures (NetworkX/matplotlib)."""
from pathlib import Path
import pandas as pd,numpy as np,json,networkx as nx,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
ROOT=Path(__file__).resolve().parents[1]
font=FontProperties(fname='C:/Windows/Fonts/simhei.ttf');plt.rcParams['font.family']=font.get_name();plt.rcParams['axes.unicode_minus']=False
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'figure.facecolor':'white','savefig.facecolor':'white'})
COL={'herb':'#147D83','compound':'#8FB7AC','gene':'#354A6A','overlap':'#CD7449','pred':'#688BAB'}
def rd(n):return pd.read_csv(ROOT/'tables'/n,sep='\t',dtype=str,keep_default_na=False)
def save(d,n):d.to_csv(ROOT/'tables'/n,sep='\t',index=False,encoding='utf-8-sig')
def figsave(name):plt.savefig(ROOT/'figures'/name,dpi=220,bbox_inches='tight');plt.close()
m=rd('30_primary_herb_compound_evidence.tsv');c=rd('32_primary_unique_parent_structures.tsv');e=rd('34_primary_compound_human_target_evidence.tsv');g=rd('35_primary_human_targets.tsv');i=rd('36_GeneCards_primary_target_intersection.tsv');f=rd('02_formula_standardization.tsv');cov=rd('37_formula_primary_coverage.tsv');I=set(i.ncbi_gene_id)
labels=dict(zip(f.herb_id,f.prescription_label))
# Explicitly separate the herb-membership source from the compound-target evidence source.
joined=m.rename(columns={'source_url':'membership_source_url','source_raw':'membership_source_raw'}).merge(e.drop(columns=['herb_id','prescription_label']).rename(columns={'source_url':'target_edge_source_url','source_raw':'target_edge_source_raw'}),on='compound_key')
joined['prescription_label']=joined.herb_id.map(labels);joined['in_GeneCards_intersection']=joined.ncbi_gene_id.isin(I)
save(joined,'50_FULL_herb_compound_human_gene_provenance.tsv')
triples=joined[['herb_id','prescription_label','compound_key','ncbi_gene_id','approved_symbol','evidence_class','in_GeneCards_intersection']].drop_duplicates();save(triples,'51_unique_herb_compound_gene_evidence.tsv')
nodes=[]
for r in f.to_dict('records'):nodes.append({'id':r['herb_id'],'label':r['prescription_label'],'node_type':'herb','relative_input':r['relative_input'],'ncbi_gene_id':'','in_GeneCards_intersection':False})
for r in c.to_dict('records'):nodes.append({'id':'C:'+r['compound_key'],'label':r['compound_names'],'node_type':'parent_structure','inchikey':r['compound_key'],'canonical_isomeric_smiles':r['canonical_isomeric_smiles'],'source_compound_ids':r['source_compound_ids'],'in_GeneCards_intersection':False})
for r in g.to_dict('records'):nodes.append({'id':'G:'+r['ncbi_gene_id'],'label':r['approved_symbol'],'node_type':'human_gene','ncbi_gene_id':r['ncbi_gene_id'],'hgnc_id':r['hgnc_id'],'in_GeneCards_intersection':r['ncbi_gene_id'] in I})
edges=[]
for (hid,key),d in m.groupby(['herb_id','compound_key']):edges.append({'source':hid,'target':'C:'+key,'interaction':'database_herb_parent_constituent','evidence_class':'database_membership','source_ids':';'.join(sorted(set(d.source_compound_id))),'source_urls':' | '.join(sorted(set(d.source_url)))})
for (key,gid),d in e.groupby(['compound_key','ncbi_gene_id']):edges.append({'source':'C:'+key,'target':'G:'+gid,'interaction':'compound_human_target','evidence_class':';'.join(sorted(set(d.evidence_class))),'source_ids':';'.join(sorted(set(d.MOL_ID))),'source_urls':' | '.join(sorted(set(d.source_url))),'new_prediction_run':False})
nd=pd.DataFrame(nodes).fillna('');ed=pd.DataFrame(edges).fillna('');save(nd,'52_Cytoscape_formula_nodes.tsv');save(ed,'53_Cytoscape_formula_edges.tsv')
G=nx.Graph()
for r in nd.to_dict('records'):G.add_node(r.pop('id'),**r)
for r in ed.to_dict('records'):s=r.pop('source');t=r.pop('target');G.add_edge(s,t,**r)
nx.write_graphml(G,ROOT/'results/ModelA_formula_network.graphml')
# All PPI hubs can be traced to specific herbs and intact structures using this table.
pn=rd('46_STRING_PPI_nodes.tsv');pe=rd('45_STRING_PPI_edges.tsv')
hub=pn[['stringId','ncbi_gene_id','approved_symbol','degree','weighted_degree']].merge(triples,on=['ncbi_gene_id','approved_symbol'])
save(hub,'54_PPI_node_herb_compound_trace.tsv')
P=nx.Graph();P.add_nodes_from(pn.stringId)
for r in pe.to_dict('records'):P.add_edge(r['source'],r['target'],weight=float(r['score']))
nx.set_node_attributes(P,dict(zip(pn.stringId,pn.approved_symbol)),'label');nx.write_graphml(P,ROOT/'results/ModelA_STRING_PPI.graphml')
# Illustrative cases selected for drug-specific interpretability, never to define analysis input.
spec=[('H11','MOL000210','CHRM3'),('H11','MOL005955','CHRM3'),('H10','MOL000555','CHRM2'),('H05','MOL007930','PTGS2'),('H09','MOL005812','TNF'),('H15','MOL000922','SLC6A1')]
examples=[]
for hid,mol,gene in spec:
 q=joined[(joined.herb_id==hid)&(joined.MOL_ID==mol)&(joined.approved_symbol==gene)]
 if q.empty and hid=='H15':q=joined[(joined.herb_id==hid)&(joined.compound_name.str.contains('aminobutyric',case=False))&(joined.approved_symbol=='SLC6A1')]
 if len(q):examples.append(q.iloc[0])
ex=pd.DataFrame(examples);save(ex,'55_reviewable_formula_examples_FULL_provenance.tsv')
# Coverage: every prescribed herb retained; zeros are database coverage counts.
fig,axs=plt.subplots(1,2,figsize=(10,7),sharey=True);y=np.arange(len(cov))
for ax,col,title,color in [(axs[0],'primary_structure_count','纳入的唯一结构',COL['herb']),(axs[1],'disease_intersection_targets','便秘交集人靶点',COL['overlap'])]:
 vals=cov[col].astype(int);ax.barh(y,vals,color=color,height=.64);ax.set_title(title,loc='left',fontweight='bold');ax.set_xlim(0,max(vals)*1.17)
 for yy,v in zip(y,vals):ax.text(v+max(vals)*.018,yy,str(v),va='center',fontsize=9)
 ax.grid(axis='x',alpha=.15);ax.set_axisbelow(True)
axs[0].set_yticks(y,labels=cov.prescription_label);axs[0].invert_yaxis();fig.suptitle('20味药逐味覆盖：未覆盖药味也保留在分母中',x=.1,ha='left',fontweight='bold');plt.tight_layout();figsave('01_herb_coverage.png')
# Full three-layer network; all input herbs, all qualified structures and human targets.
fig,ax=plt.subplots(figsize=(10,6));pos={}
groups=[list(f.herb_id),['C:'+x for x in c.compound_key],['G:'+x for x in g.ncbi_gene_id]]
for x,ids in enumerate(groups):
 for j,n in enumerate(ids):pos[n]=(x,1-j/max(1,len(ids)-1))
nx.draw_networkx_edges(G,pos,ax=ax,alpha=.035,width=.3,edge_color='#819095')
for ids,color,size in zip(groups,[COL['herb'],COL['compound'],COL['gene']],[65,2,8]):nx.draw_networkx_nodes(G,pos,nodelist=ids,node_size=size,node_color=color,ax=ax)
for n in groups[0]:ax.text(-.04,pos[n][1],labels[n],ha='right',va='center',fontsize=8)
for x,n,title in zip(range(3),[len(z) for z in groups],['处方药味','原形候选结构','人靶点']):ax.text(x,1.08,f'{title}\n{n:,}',ha='center',va='bottom',fontweight='bold')
ax.set_xlim(-.4,2.2);ax.set_ylim(-.03,1.22);ax.axis('off');ax.set_title(f'全量网络 | {G.number_of_nodes():,}节点 · {G.number_of_edges():,}条唯一边',loc='left');figsave('02_full_formula_network.png')
# True PPI, all nodes including isolates. Fixed seed; edges never invented or thresholded again.
fig,ax=plt.subplots(figsize=(10,7))
connected=[x for x,d in P.degree() if d>0];isolates=sorted(set(P)-set(connected))
# Lay out the connected component separately so isolates do not compress it.
# The same 292 nodes and all 4349 edges remain in the figure.
pos=nx.spring_layout(P.subgraph(connected),seed=20261010,weight='weight',iterations=150,k=.13)
for j,x in enumerate(isolates):pos[x]=np.array([1.6,.45-j*.3])
nx.draw_networkx_edges(P,pos,ax=ax,alpha=.065,width=.42,edge_color='#697B88')
degrees=dict(P.degree());nx.draw_networkx_nodes(P,pos,ax=ax,node_size=[10+degrees[x]*.7 for x in P],node_color=[degrees[x] for x in P],cmap='YlOrBr',linewidths=.1,edgecolors='white')
top=sorted(P,key=lambda x:(-degrees[x],x))[:12];smap=dict(zip(pn.stringId,pn.approved_symbol))
top=sorted(top,key=lambda x:np.arctan2(pos[x][1],pos[x][0]))
for j,x in enumerate(top):
 angle=-np.pi+(j+.5)*2*np.pi/len(top);lp=np.array([np.cos(angle),np.sin(angle)])*1.14
 ax.annotate(smap[x],xy=pos[x],xytext=lp,ha='center',va='center',fontsize=9,arrowprops=dict(arrowstyle='-',lw=.55,color='#546E7A'),bbox=dict(boxstyle='round,pad=.15',facecolor='white',edgecolor='none',alpha=.85))
for x in isolates:ax.text(pos[x][0]+.07,pos[x][1],smap[x],fontsize=9,va='center')
ax.text(1.6,.8,'孤立节点',ha='center',fontsize=10,color='#546E7A')
ax.set_title(f'STRING 12.5 人功能互作 | {len(P)}节点 · {P.number_of_edges():,}边',loc='left');ax.text(.0,-.02,'score ≥ 0.4；不添加节点；标注度数最高的12个节点及全部孤立节点。',transform=ax.transAxes,fontsize=9);ax.set_xlim(-1.35,2.);ax.set_ylim(-1.3,1.3);ax.set_aspect('equal');ax.axis('off');figsave('03_STRING_PPI.png')
r=pd.read_csv(ROOT/'tables/40_GO_KEGG_ALL_tested_results.tsv',sep='\t');primary=r[r.analysis.str.startswith('PRIMARY')]
q=primary.head(8).iloc[::-1].copy();fig,ax=plt.subplots(figsize=(10,4.5));colors=[COL['herb'] if x<.05 else '#A7AFB9' for x in q.p_adjust_BH_all_families]
ax.barh(np.arange(len(q)),-np.log10(q.p_adjust_BH_all_families.clip(lower=1e-300)),color=colors)
ax.set_yticks(np.arange(len(q)),labels=[a+'\n'+b for a,b in zip(q.term_id,q.term_name)],fontsize=8);ax.axvline(-np.log10(.05),color=COL['overlap'],linestyle='--',linewidth=1);ax.set_xlabel('-log10（全GO+KEGG条目BH校正值）');ax.set_title('主分析：交集靶点相对全处方靶点的功能富集',loc='left');plt.tight_layout();figsave('04_primary_enrichment.png')
sig=r.assign(sig=r.p_adjust_BH_all_families<.05).groupby(['analysis','family']).sig.sum().unstack().reindex(columns=['GO_BP','GO_MF','GO_CC','KEGG'])
sig=sig.reindex(['PRIMARY_intersection_vs_all_formula_targets','SENSITIVITY_intersection_vs_human_protein_coding','CONTROL_all_formula_targets_vs_human_protein_coding'])
fig,ax=plt.subplots(figsize=(10,4));x=np.arange(4)
for j,(name,row) in enumerate(sig.iterrows()):
 b=ax.bar(x+(j-1)*.25,row.values,width=.23,label=['主分析：交集 vs 处方','敏感性：交集 vs 人基因组','对照：全处方 vs 人基因组'][j],color=[COL['herb'],COL['overlap'],'#9FADB9'][j])
 ax.bar_label(b,fontsize=8,padding=3)
ax.set_xticks(x,labels=['GO生物过程','GO分子功能','GO细胞组分','KEGG']);ax.set_ylabel('BH < 0.05条目数');ax.legend(frameon=False,fontsize=8);ax.set_title('不同背景回答不同问题；不能把显著性差异当作验证',loc='left');plt.tight_layout();figsave('05_background_comparison.png')
# Key evidence paths, each row is an actually observed database chain.
fig,ax=plt.subplots(figsize=(10,4.2));ax.axis('off')
for j,rr in enumerate(ex.to_dict('records')):
 yy=1-j/max(1,len(ex)-1);color=COL['herb'] if 'validated' in rr['evidence_class'] else COL['pred']
 for xx,label in [(0,rr['prescription_label']),(.46,rr['compound_name']+'\n'+rr['MOL_ID']),(.91,rr['approved_symbol'])]:ax.text(xx,yy,label,ha='center',va='center',fontsize=9,bbox=dict(boxstyle='round,pad=.4',facecolor='white',edgecolor=color))
 ax.plot([.06,.33],[yy,yy],color='#9AABA9',lw=1.2);ax.plot([.62,.84],[yy,yy],color=color,lw=1.5,linestyle='-' if 'validated' in rr['evidence_class'] else '--')
ax.set_xlim(-.1,1.03);ax.set_ylim(-.15,1.15);ax.set_title('可回溯实例 | 实线：源库validated标签；虚线：源库已有预测',loc='left',fontsize=11);figsave('06_formula_evidence_examples.png')
summary={'formula_network_nodes':len(G),'formula_network_edges':G.number_of_edges(),'herb_structure_edges':len(m[['herb_id','compound_key']].drop_duplicates()),'compound_gene_edges':len(e[['compound_key','ncbi_gene_id']].drop_duplicates()),'full_source_provenance_rows':len(joined),'unique_herb_compound_gene_evidence_rows':len(triples),'example_count':len(ex),'Cytoscape_run':False,'actual_visualization_tools':['NetworkX','matplotlib'],'PPI_mapping_exclusions':sorted(I-set(pn.ncbi_gene_id))}
(ROOT/'results/network_export_summary.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2),encoding='utf-8');print(json.dumps(summary,ensure_ascii=False,indent=2))
