"""Final integrity audit, amendments, source index and Chinese delivery guide."""
from pathlib import Path
from datetime import datetime,timezone
import pandas as pd,json,hashlib,sys,platform,importlib.metadata
ROOT=Path(__file__).resolve().parents[1]
def rd(n):return pd.read_csv(ROOT/'tables'/n,sep='\t',dtype=str,keep_default_na=False)
def dump(n,o):(ROOT/n).write_text(json.dumps(o,ensure_ascii=False,indent=2),encoding='utf-8')
b=json.loads((ROOT/'results/primary_baseline_summary.json').read_text(encoding='utf-8'));ppi=json.loads((ROOT/'results/STRING_PPI_summary.json').read_text(encoding='utf-8'))
c=rd('32_primary_unique_parent_structures.tsv');m=rd('30_primary_herb_compound_evidence.tsv');e=rd('34_primary_compound_human_target_evidence.tsv');g=rd('35_primary_human_targets.tsv');i=rd('36_GeneCards_primary_target_intersection.tsv');pe=rd('45_STRING_PPI_edges.tsv');pn=rd('46_STRING_PPI_nodes.tsv');t=rd('50_FULL_herb_compound_human_gene_provenance.tsv');r=rd('40_GO_KEGG_ALL_tested_results.tsv')
checks={
 'unique_structures_match_summary':c.compound_key.nunique()==b['primary_unique_parent_structures'],
 'all_edges_have_qualified_parent':set(e.compound_key)<=set(c.compound_key),
 'all_edges_have_human_taxonomy':set(e['Taxonomy id'])=={'9606'},
 'all_edges_have_HGNC_and_GeneID':e.hgnc_id.ne('').all() and e.ncbi_gene_id.ne('').all(),
 'no_qt_memberships':not m.compound_name.str.contains(r'_qt\b',case=False,regex=True).any(),
 'no_new_predictions_claimed':set(e.new_prediction_run)=={'False'},
 'intersection_subset_of_formula':set(i.ncbi_gene_id)<=set(g.ncbi_gene_id),
 'all_compound_gene_pairs_trace_to_herbs':set(zip(e.compound_key,e.ncbi_gene_id))<=set(zip(t.compound_key,t.ncbi_gene_id)),
 'PPI_nodes_subset_of_intersection':set(pn.ncbi_gene_id)<=set(i.ncbi_gene_id),
 'PPI_no_added_nodes':set(pe.source)|set(pe.target)<=set(pn.stringId),
 'PPI_score_threshold_fixed':(pe.score.astype(float)>=.4).all(),
 'all_primary_tests_retained':len(r[r.analysis.str.startswith('PRIMARY')])==7786,
 'zero_overlap_tests_retained':(r.k_overlap.astype(int)==0).any(),
 'valid_hypergeometric_denominators':((r.k_overlap.astype(int)<=r.K_term_background.astype(int))&(r.k_overlap.astype(int)<=r.n_query.astype(int))&(r.n_query.astype(int)<=r.N_background.astype(int))).all(),
 'PPI_edges_count_matches':len(pe)==ppi['edges_score_at_least_0_4']}
checks={k:bool(v) for k,v in checks.items()};dump('results/final_integrity_checks.json',{'time_utc':datetime.now(timezone.utc).isoformat(),'checks':checks,'all_passed':all(checks.values())});assert all(checks.values()),checks
b['go_kegg_status']='completed; results/enrichment_summary.json';b['ppi_status']='completed; results/STRING_PPI_summary.json';b['final_scope_note']='47 unconfirmed _qt derivative structures excluded; final 2059 structures / 331 targets / 294 intersection supersede preliminary counts.';dump('results/primary_baseline_summary.json',b)
amend={'date_utc':datetime.now(timezone.utc).isoformat(),'supersedes_not_deletes':'results/analysis_protocol.json retains original preflight design','amendments':[
 {'item':'disease input','basis':'User manual GeneCards constipation CSV made available after initial blocked-source preflight','rule':'One fixed source search term constipation, all 9230 valid symbols no score cutoff; exact copyright footer excluded; HPO never primary.'},
 {'item':'primary enrichment background','basis':'Parent/user steering before enrichment requested formula-target background','rule':'294 intersection vs 331 formula targets; separately full human annotated protein-coding sensitivity and all-formula control. BH across all tested GO+KEGG families per analysis; family BH also saved.'},
 {'item':'material aliases','basis':'Official HKCMMS definitions with exact source pages','rule':'Add 10 Cundong ingredient records and 146 Gan Dihuang records; exclude conflicting Liriope/fresh/prepared material entries.'},
 {'item':'intact parent final audit','basis':'User pre-existing Model A only requirement; derivative tag was detected in final provenance inspection','rule':'Exclude all source membership names matching _qt\\b without independent intact membership. Removes 47 structures, 24 target-linked structures, 106 predicted edges, one target/intersection member. This was a scope correction after preliminary enrichment, not a significance-driven filter. Recomputed all analyses.'},
 {'item':'target evidence scope','basis':'Source evidence classification and explicit human ID requirement','rule':'Main target network uses TCMSP existing labelled edges only after HERB explicit TAR xref Tax9606. All 29559 HERB unclassified/ref-mining rows remain supplemental. Name-only UniProt inferred map never promoted.'},
 {'item':'STRING identity mapping','basis':'Actual mismatched preferred names returned by official API','rule':'Exact HGNC symbol preferred. ACP3 accepted through official HGNC UniProt P15309 and previous name ACPP. MAPK10/VEGFA unresolved and excluded only from PPI. No added nodes. Official stable host timeout; used official main API, no proxy.'}]
 };dump('results/final_protocol_amendments.json',amend)
req=[]
for line in (ROOT/'logs/http_requests.jsonl').read_text(encoding='utf-8').splitlines():
 try:req.append(json.loads(line))
 except json.JSONDecodeError:pass
flat=[]
for q in req:flat.append({k:q.get(k,'') for k in ['time_utc','label','url','final_url','status','bytes','sha256','raw_file','error','proxy_mode']})
pd.DataFrame(flat).to_csv(ROOT/'tables/60_download_source_dates_checksums.tsv',sep='\t',index=False,encoding='utf-8-sig')
vers={p:importlib.metadata.version(p) for p in ['pandas','numpy','scipy','statsmodels','networkx','matplotlib','requests','reportlab','pdfplumber']};vers['rdkit_workspace']='2026.9.1'
dump('logs/final_environment.json',{'python':sys.version,'platform':platform.platform(),'packages':vers,'proxy_policy':'all HTTP downloads in direct_http remove all case variants of proxy env; trust_env=False; proxies={}; NO_PROXY=*; TLS verify=True','source_response_versions':[{'label':x.get('label'),'headers':x.get('source_headers')} for x in req if x.get('source_headers',{}).get('X-UniProt-Release')]})
notes='''改良通便丸 Model A 网络药理实际结果 - 内部审阅初稿
最终口径：20味处方，19味有原形候选结构，17味有合格人靶点。
2059个唯一结构；580个结构有合格靶点；331个人靶点；4433条结构-人基因配对。
其中485条为TCMSP已有validated标签，3948条为TCMSP已有SVM/RF预测；没有本项目新预测或实验。
GeneCards：9230条有效记录/唯一原始符号，8877个唯一映射NCBI基因，294个处方交集（88.82%）。
STRING：292个可核验节点，4349条真实功能互作边，3个孤立节点；不添加节点，阈值0.4。
主GO/KEGG：7786个条目全部测试，统一BH<0.05只有GO:0032501一项；KEGG显著数0。
全人基因组背景敏感性回答不同问题，不作独立验证；全部处方靶点的人背景对照一并报告。

先读 output/pdf/改良通便丸_ModelA_网络药理结果_审阅初稿.pdf。
此稿供根审；全部数值以最终表和final_integrity_checks.json为准。
早期阶段报告中的2032/2106结构、332靶点、295交集已被最终模型审计更新，不应当作最终结果引用。

主要全表：
02_formula_standardization.tsv：精确20味、原投料、物种/部位候选及炮制约束；总量340，单位不明。
30_primary_herb_compound_evidence.tsv：全部合格药味-结构来源。药材数据库记录不等于实际制剂测定。
31_primary_membership_exclusions.tsv：结构、药材歧义、47条_qt衍生物等全部排除记录。
32_primary_unique_parent_structures.tsv：2059唯一完整InChIKey及SMILES，保留所有药味来源。
34_primary_compound_human_target_evidence.tsv：4433配对，含原始药材页靶点边来源、人Tax9606交叉标识、证据类别。
注意34中的herb_id是获取靶点边时的原药材页，不表示完整处方归属；药味归属请看50/51。
35、36、37：人靶点、GeneCards交集、逐味覆盖。
40_GO_KEGG_ALL_tested_results.tsv：所有三个分析、全部GO/KEGG条目，包括零命中及非显著结果。
41、42：每个分析和功能体系的完整基因背景及有效分母；43为显著结果便利子集，不替代全表。
44/44b、45/46：STRING原始映射与核验映射、真实加权边和节点。
50_FULL_herb_compound_human_gene_provenance.tsv：14022条完整来源证据，可从每个配对回到药味、结构和原始缓存。
51：8745条去重药味-结构-基因证据组合；54为每个PPI节点的全部药味与结构回溯；55为6条可讲实例。
52_Cytoscape_formula_nodes.tsv / 53_Cytoscape_formula_edges.tsv：导入Cytoscape，id连接source/target。
results/ModelA_formula_network.graphml 与 ModelA_STRING_PPI.graphml可直接导入；不宣称Cytoscape已运行。
60_download_source_dates_checksums.tsv：下载来源、时间和校验和；失败请求也保留。
源库已有全部候选、HERB未分型关联、名称推断人靶点、HPO备选均另表保存，不混入主结论。

复算顺序（Python 3.12.7，本机已有依赖；RDKit在vendor中，来源官方PyPI，版本/轮子SHA256留存）：
38_resolve_materials.py → 34_primary_baseline.py → 36_enrichment.py → 37_string_ppi.py → 39_exports_figures.py → 40_client_pdf.py → 41_finalize_audit.py。
37使用缓存；缓存缺失时会访问官方STRING。网络进程需明确清空HTTP_PROXY/HTTPS_PROXY/ALL_PROXY、NO_PROXY=*。
原始数据抓取与结构核查脚本01-33均保留；通常无需重新访问网站。
JSON/脚本UTF-8；表TSV为UTF-8 BOM，可用Excel文本导入，ID列按文本处理。

证据与范围：
处方一半研粉一半煎汤不视为50%化学暴露；没有按未知剂量加权，无OB/DL主硬筛选。
按药材数据库母体候选构建，未加入推导菌转化物。_qt衍生记录排除；天然母体逐条文献与丸剂含量未实测。
多糖、脂肪油混合物等无明确单结构的部分保留覆盖边界，不虚配单SMILES或靶点。
TCMSP validated是源库标签，不等同亲和力确认、原始人实验确认或本项目实验。
来源靶点身份映射到人，不证明每个源实验都在人源体系完成。
SEA需验证码，Swiss条款禁止自动访问，二者本次未运行。GeneCards使用用户CSV，无自动抓取。
瓜蒌仁未被整果/果皮代替；铁皮石斛未用泛石斛；生地黄未用熟地。
麦冬/生地黄同义条目有官方标准支撑，证据见28及material_resolution_audit.json。
PPI是功能关联网络而非物理结合专属网络，保留STRING所有通道分数；不视为疾病独立验证。
没有旧GSE/转录组、菌群、宏基因组或MD结果。

许可与外发：
TCMSP页面声明ODbL；其他数据库及标准各有独立许可/引用要求，公开可读接口不意味着无限商业再分发。
用户GeneCards原CSV包含LifeMap版权及非学术用途需书面许可声明，原件和页脚已保留。
此目录用于内部复核；未经许可核实，不应把GeneCards完整表或其他受限原始材料公开发布。
本次未发送第三方消息、未上传完整GeneCards；STRING仅使用必要交集标识符。

完整性：
logs/final_sha256_manifest.tsv保存交付目录全部文件SHA256（不含自身、临时渲染图、Python缓存、vendor运行库）。
源响应与最终成品校验和可据此复核；初始协议与方法修订分别保留，未覆盖历史设计。
'''
(ROOT/'results/FINAL_README_中文.txt').write_text(notes,encoding='utf-8-sig')
print(json.dumps(checks,ensure_ascii=False,indent=2))
