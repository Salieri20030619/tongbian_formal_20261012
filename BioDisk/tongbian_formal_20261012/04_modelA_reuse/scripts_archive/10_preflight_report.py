from pathlib import Path
import json,hashlib,pandas as pd
from datetime import datetime,timezone
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1]
cov=pd.read_csv(ROOT/'tables'/'04_tcmsp_herb_coverage.tsv',sep='\t')
s=json.loads((ROOT/'results'/'tcmsp_acquisition_summary.json').read_text(encoding='utf-8'))
h=json.loads((ROOT/'results'/'hpo_candidate_coverage.json').read_text(encoding='utf-8'))
knowledge=pd.read_csv(ROOT/'raw'/'diseases_knowledge.bin',sep='\t',header=None,names=['gene_id','gene_symbol','disease_id','disease_name','source','evidence','confidence'])
knowledge[knowledge.disease_name.str.contains('constipation',case=False,na=False)].to_csv(ROOT/'tables'/'09_DISEASES_knowledge_constipation_CANDIDATE.tsv',sep='\t',index=False,encoding='utf-8-sig')
plt.rcParams['font.sans-serif']=['Microsoft YaHei','SimHei','DejaVu Sans']
plt.rcParams['axes.unicode_minus']=False
fig,ax=plt.subplots(figsize=(10,9))
ax.barh(range(20),cov.ingredient_rows.fillna(0),color=['#218b8a' if pd.notna(v) else '#b5b5b5' for v in cov.ingredient_rows])
ax.set_yticks(range(20),cov.prescription_label);ax.invert_yaxis()
for n,val in enumerate(cov.ingredient_rows):ax.text((0 if pd.isna(val) else val)+2,n,'目录无精确匹配' if pd.isna(val) else str(int(val)),va='center',fontsize=9,color='#425466')
ax.set_xlim(0,245);ax.set_xlabel('TCMSP原始药材—成分记录数（尚未结构去重）')
ax.set_title('20味药实际获取覆盖：16味名称匹配，1,401条成分记录',loc='left',fontsize=14,pad=18)
ax.spines[['top','right']].set_visible(False)
fig.text(.13,.025,'缺失项不记作零成分；不以瓜蒌替瓜蒌仁、不以熟地黄替生地黄。\n名称匹配尚不等于物种、药用部位和炮制全部核验。',fontsize=10,color='#5a6571')
fig.subplots_adjust(left=.15,bottom=.11,right=.97,top=.92)
fig.savefig(ROOT/'figures'/'01_实际成分获取覆盖.png',dpi=180);plt.close(fig)

statuses=[
('TCMSP旧版','可用；实际下载','16/20味药精确中文名目录匹配；1401条药材—成分记录','https://old.tcmsp-e.com/tcmsp.php'),
('HERB旧版','前端可读，公开导出失败','3个公开目录下载均HTTP200但正文Sorry, unavailable path.；不是有效数据','http://herb.ac.cn/'),
('HERB2','仅首页验证','v2首页HTTP200；成分详情和导出尚未完成','http://herb.ac.cn/v2'),
('PubChem','可用；实际结构样本','TCMSP MOL000018→CID6321405→InChIKey一致；未批量获取','https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/6321405/property/InChIKey,CanonicalSMILES/JSON'),
('SEA','前端可读；预测未执行','当前前端每次最多10分子并要求验证码；挑战接口返回实际挑战，未求解、未提交；旧/search返回404','https://sea.bkslab.org/'),
('SwissTargetPrediction','条款阻止自动化；预测未执行','已下载条款，禁止自动化访问；不自动提交','https://www.swisstargetprediction.ch/termsofuse.php'),
('GeneCards','HTTP403；未使用','后续获知条款禁止自动抓取，停止；未下载疾病基因','https://www.genecards.org/'),
('CTD','验证码；未使用','返回人机验证页，停止','https://ctdbase.org/downloads/'),
('DisGeNET','仅首页可读','未访问授权数据接口，不能声称数据可用','https://disgenet.com/'),
('DISEASES知识通道','全表下载成功；便秘覆盖0','98378行，名称含constipation的行数0；不是便秘没有相关基因','https://download.jensenlab.org/human_disease_knowledge_full.tsv'),
('DISEASES实验通道','下载超时','本轮没有数据；不重刷','https://download.jensenlab.org/human_disease_experiments_full.tsv'),
('HPO','全表下载成功；候选源','1101396行；固定词精确匹配819疾病/表型—基因行，525个NCBI基因；尚未纳入正式疾病分析','https://purl.obolibrary.org/obo/hp/hpoa/phenotype_to_genes.txt'),
('NCBI Eutils','样本查询成功','MedGen constipation检索可用，尚未完整构建疾病基因集','https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=medgen&term=constipation&retmode=json'),
('ChEMBL','API状态成功','仅状态测试；尚无成分活性证据下载','https://www.ebi.ac.uk/chembl/api/data/status.json'),
('UniProt','人蛋白样本成功','P04637详情JSON可用；处方靶点映射未执行','https://rest.uniprot.org/uniprotkb/P04637.json'),
('STRING','标识符样本成功','人TP53映射预检；未构建处方PPI，样本不进入网络','https://string-db.org/api/tsv/get_string_ids?identifiers=TP53&species=9606&limit=1'),
('KEGG','人通路目录成功','尚未获取分析背景与运行富集','https://rest.kegg.jp/list/pathway/hsa'),
('g:Profiler','物种接口成功','未运行富集','https://biit.cs.ut.ee/gprofiler/api/util/organisms_list'),
('IMPPAT','仅首页可读','未纳入主分析','https://cb.imsc.res.in/imppat/'),
('SymMap','连接被拒','未下载数据','https://www.symmap.org/'),
('Open Targets','端点可达但空查询400','无有效疾病数据查询；不能计为数据获取成功','https://api.platform.opentargets.org/api/v4/graphql')]
pd.DataFrame(statuses,columns=['source','actual_status','detail','source_url']).to_csv(ROOT/'tables'/'10_reviewed_source_status.tsv',sep='\t',index=False,encoding='utf-8-sig')
text=f'''改良通便丸 Model A｜本机直连预检与真实数据获取阶段报告
状态：阶段结果，供根审阅；不是完整网络药理结果，不是疗效验证。
生成时间：{datetime.now(timezone.utc).isoformat()}

一、可行性判断
本机无代理下载路径已经验证可用，TCMSP旧版能够实际取得成分与数据库已有靶点关联，PubChem/UniProt/STRING/KEGG的样本请求成功。因此有条件继续建立原形成分证据基线。
原定全工具链目前不能照单自动完成：SwissTargetPrediction条款禁止自动访问；SEA当前页面要求验证码；GeneCards403且禁止自动抓取；CTD要求验证码。没有绕过限制，也没有声称已运行这两种预测。
当前结果不足以判断处方对便秘有效，亦不能报告人靶点、疾病交集、PPI或显著通路。上述指标为NA（尚未完成），不是0。

二、真实取得的数据和计数口径
处方20味，相对投料值合计340，单位未给定，不能写作g/day。每味一半研粉、一半煎汤，不换算为化学暴露50%。
已建立20味标准化表，基原候选和药用部位作为待验证检索范围，实际供应物种未提供，未当作已核验样本。
TCMSP旧版公开目录共502个条目，16/20味精确中文名匹配。全部取得的成分记录均保留，未按OB/DL删除。
药材—成分行数：{s['herb_ingredient_rows']}。
TCMSP唯一MOL_ID：{s['unique_tcmsp_mol_ids_NOT_structure_dedup']}（不是按InChIKey去重的原形成分数）。
药材—成分—靶点原始行数：{s['herb_compound_target_rows_UNCLASSIFIED']}。
唯一成分—靶点ID配对：{s['unique_compound_target_pairs_all_species_unverified']}。
唯一源靶点ID：{s['unique_target_ids_all_species_unverified']}（包含非人靶点，不能写成805个人靶点）。
源数据库validated标签行数：{s['database_validated_label_rows']}（尚未核对每条原始实验，不能宣称本项目实验验证）。
源SVM_score/RF_score与validated标签均保留；这些是TCMSP已有证据和预测，不是本次SEA/Swiss预测。
没有精确匹配：瓜蒌仁、麦冬、生地黄、铁皮石斛。未用瓜蒌、熟地黄或泛石斛替换。
结构小样本：MOL000018（源名(+/-)-Isoborneol）→TCMSP给出PubChem CID6321405→PubChem返回同一InChIKey DTGKSKDOIYIVQL-MRTMQBJTSA-N。源名和立体信息仍需批量一致性审核；此例只验证链路。

三、疾病源替代候选——供审阅，不冒称GeneCards
固定词：constipation、chronic constipation、functional constipation、chronic idiopathic constipation、slow transit constipation。
HPO公开phenotype_to_genes全表{h['download_rows']}行。精确名称匹配命中HP:0002019 Constipation和HP:0012450 Chronic constipation，共{h['exact_term_match_rows']}行、{h['unique_genes']}个NCBI Gene ID。
Constipation：673行、525个基因；Chronic constipation：146行、140个基因，后者基因包含在合并525中。
它们多来自遗传病/综合征的表型注释，不能等同于慢传输型/功能性便秘致病基因，也不能充作独立验证。尚未用于交集、PPI或富集。
DISEASES知识通道全表98378行，疾病名称含constipation的关联0行。实验通道下载超时。不是可以据此宣称便秘无相关基因。
NCBI MedGen检索样本可用；如需更贴近疾病的主集合，应进一步固定MedGen疾病ID及NCBI关联纳入规则，确认可下载覆盖后再决定HPO仅作补充。

四、下一步完整、可审计路线
1. 优先核验HERB2公开详情/下载，补足四味缺口；旧版公开下载失败不能推广为HERB2全部不可用。全表来源、药材部位/炮制与物种须逐项核对。
2. 批量取得1097个TCMSP源成分的结构标识，依据完整InChIKey去重，保留药材归属、源名/立体冲突和未解析结构。混合多糖/脂肪油只列边界。
3. 对源靶点逐条核验物种和正规ID，保留非人/模糊映射排除日志。将数据库validated、源预测、直接文献/ChEMBL实验分层。若原始来源不能支持实测类别，只能标数据库声称validated。
4. 根审阅疾病替代来源后，形成固定疾病集合。主富集在独立于疾病筛选的全部合格人靶点上进行；交集仅描述性。
5. STRING以9606、required_score=400、add_nodes=0获取真实网络边；GO/KEGG保留可追溯背景及完整测试家族，BH全量校正。无预测靶点宇宙时说明人类注释背景的探索性限制。
6. 导出Cytoscape节点/边表。当前Cytoscape和R未运行；本轮仅Python和matplotlib用于真实覆盖图。根审阅可信分析后再排客户中文PDF。

五、网络与可复现性
第一次沙箱内所有域名DNS失败，未据此认定网站不可用。经允许的本机联网执行复核后成功。
所有下载脚本清除大小写HTTP_PROXY/HTTPS_PROXY/ALL_PROXY，NO_PROXY=*，requests.Session.trust_env=False，proxies={{}}，TLS验证开启。不修改系统/Clash设置，不使用云web下载，不自动切代理。
logs/http_requests.jsonl保留每次时间、请求、状态、字节数、SHA256和错误；raw保留真实响应（包含错误页，不能都视为数据）。
scripts/07_reparse_tcmsp.py修正了字符串含“],”时的初始正则截断问题，改用JSONDecoder.raw_decode并重新生成完整表。当前表只将含molecule_ID和target_name的数组计作成分—靶点边，靶点—疾病数组另存且不用于便秘分析。
协议在看靶点/疾病结果前保存，analysis_protocol.json及初始哈希可查。

六、交付索引
tables/02_formula_standardization.tsv：20味药与约束。
tables/04_tcmsp_herb_coverage.tsv：实际16/20覆盖及逐味成分记录分母。
tables/05_tcmsp_ingredient_rows.tsv：全部1401条原始成分记录。
tables/06_tcmsp_target_rows_UNCLASSIFIED.tsv：11702条未分层、未限人的源关联。
tables/07_HPO_constipation_CANDIDATE_NOT_PRIMARY.tsv：819条候选疾病/表型关联。
tables/10_reviewed_source_status.tsv：所有来源实际状态与链接。
figures/01_实际成分获取覆盖.png：真实覆盖图。
results/tcmsp_acquisition_summary.json：机器可读计数与NA。
scripts/：本次实际执行脚本；logs/：请求、版本及协议哈希；raw/：原始响应。

剩余阻塞：两种网页预测未获可自动执行途径；GeneCards/CTD不可自动获取；四味精确成分覆盖待补；结构和人靶点映射未完成；疾病替代需根审阅后推进。未生成伪网络或富集图。中文客户PDF依任务要求留待根审阅后排版。
'''
(ROOT/'results'/'直连预检与数据获取阶段报告_待根审阅.txt').write_text(text,encoding='utf-8-sig')
(ROOT/'README_阶段结果.txt').write_text('本目录为预检与数据获取阶段成果，并非完整网络药理分析。请先读 results/直连预检与数据获取阶段报告_待根审阅.txt。\n当前人靶点、交集、PPI、GO/KEGG均NA，不是0。所有源靶点表标记UNCLASSIFIED，不能直接用于结论。\n',encoding='utf-8-sig')
manifest=[]
for path in sorted(ROOT.rglob('*')):
 if path.is_file() and '__pycache__' not in path.parts and path.name!='sha256_manifest.tsv':
  manifest.append({'path':path.relative_to(ROOT).as_posix(),'bytes':path.stat().st_size,'sha256':hashlib.sha256(path.read_bytes()).hexdigest()})
pd.DataFrame(manifest).to_csv(ROOT/'logs'/'sha256_manifest.tsv',sep='\t',index=False,encoding='utf-8-sig')
print(json.dumps({'files':len(manifest),'bytes':sum(x['bytes'] for x in manifest),'report':'results/直连预检与数据获取阶段报告_待根审阅.txt'},ensure_ascii=False))
