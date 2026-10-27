"""Lock formula labels and inclusion rules before examining analysis outcomes."""
from pathlib import Path
import csv,json,hashlib
from datetime import datetime,timezone
ROOT=Path(__file__).resolve().parents[1]
rows=[
('H01','生白术',30,'白术','Atractylodes macrocephala','根茎','生品，未炒制','BAI ZHU','Atractylodis Macrocephalae Rhizoma'),
('H02','生黄芪',30,'黄芪','Astragalus membranaceus; Astragalus mongholicus','根','生品，非蜜炙','HUANG QI','Astragali Radix'),
('H03','肉苁蓉',30,'肉苁蓉','Cistanche deserticola; Cistanche tubulosa','带鳞叶肉质茎','未指定炮制','ROU CONG RONG','Cistanches Herba'),
('H04','当归',15,'当归','Angelica sinensis','根','未指定炮制','DANG GUI','Angelicae Sinensis Radix'),
('H05','陈皮',10,'陈皮','Citrus reticulata','成熟果皮','干燥成熟果皮；非青皮','CHEN PI','Citri Reticulatae Pericarpium'),
('H06','瓜蒌仁',15,'瓜蒌子','Trichosanthes kirilowii; Trichosanthes rosthornii','成熟种子','非瓜蒌皮、非全瓜蒌','GUA LOU ZI','Trichosanthis Semen'),
('H07','苦杏仁',10,'苦杏仁','Prunus armeniaca; Prunus sibirica; Prunus mandshurica; Prunus armeniaca var. ansu','成熟种子','未指定燀/炒；非甜杏仁','KU XING REN','Armeniacae Semen Amarum'),
('H08','火麻仁',15,'火麻仁','Cannabis sativa','成熟果实（药名习称种仁）','未指定炮制','HUO MA REN','Cannabis Fructus'),
('H09','枳实',15,'枳实','Citrus aurantium; Citrus sinensis','幼果','非枳壳；未指定麸炒','ZHI SHI','Aurantii Fructus Immaturus'),
('H10','槟榔',10,'槟榔','Areca catechu','成熟种子','未指定焦/炒','BING LANG','Arecae Semen'),
('H11','厚朴',10,'厚朴','Magnolia officinalis; Magnolia officinalis var. biloba','干皮、根皮及枝皮','未指定姜制','HOU PO','Magnoliae Officinalis Cortex'),
('H12','桔梗',10,'桔梗','Platycodon grandiflorus','根','未指定炮制','JIE GENG','Platycodonis Radix'),
('H13','麦冬',20,'麦冬','Ophiopogon japonicus','块根','非山麦冬；未指定炮制','MAI DONG','Ophiopogonis Radix'),
('H14','玄参',20,'玄参','Scrophularia ningpoensis','根','未指定炮制','XUAN SHEN','Scrophulariae Radix'),
('H15','生地黄',20,'地黄','Rehmannia glutinosa','块根','生地黄；非熟地黄','DI HUANG','Rehmanniae Radix'),
('H16','山药',20,'山药','Dioscorea opposita; Dioscorea polystachya','根茎','未指定麸炒','SHAN YAO','Dioscoreae Rhizoma'),
('H17','铁皮石斛',15,'铁皮石斛','Dendrobium officinale','茎','不能以泛石斛或其他石斛种替代','TIE PI SHI HU','Dendrobii Officinalis Caulis'),
('H18','黄精',15,'黄精','Polygonatum kingianum; Polygonatum sibiricum; Polygonatum cyrtonema','根茎','未指定蒸/酒制','HUANG JING','Polygonati Rhizoma'),
('H19','党参',15,'党参','Codonopsis pilosula; Codonopsis pilosula var. modesta; Codonopsis tangshen','根','未指定炮制','DANG SHEN','Codonopsis Radix'),
('H20','枸杞',15,'枸杞子','Lycium barbarum','成熟果实','暂按枸杞子；原文枸杞需保留','GOU QI ZI','Lycii Fructus')
]
fields=['herb_id','prescription_label','relative_input','search_drug_name','candidate_botanical_scope_UNVERIFIED','medicinal_part_scope_UNVERIFIED','processing_constraint','pinyin_query','latin_drug_query','unit','species_status','preparation','source_status']
with (ROOT/'tables'/'02_formula_standardization.tsv').open('w',encoding='utf-8-sig',newline='') as f:
    w=csv.writer(f,delimiter='\t');w.writerow(fields)
    for r in rows:w.writerow(list(r)+['unspecified; not g/day','实际供应药材物种未提供；候选基原仅检索词，不能当作已核验样本','各取一半研粉，一半煎汤浓缩混合制丸；非50%化学暴露','用户给定处方；基原和部位为待原始数据库/标准核对的检索范围'])
protocol={
 'frozen_utc':datetime.now(timezone.utc).isoformat(),
 'model':'Model A: only intact compounds documented in the specified source medicinal material; no inferred microbial products',
 'herb_denominator':20,'relative_input_sum':sum(x[2] for x in rows),'dose_unit':'unknown; relative formulation inputs only',
 'ingredient_collection':'All records retrievable under a source database complete herb-specific query/export, before observing targets/disease results. Preserve all herb memberships and source compound IDs. No selection by effect or significance.',
 'identity':'Require documented parent chemical structure, full InChIKey preferred. Unresolved names retained in coverage, excluded from structure-dependent targets. Distinct stereochemistry kept distinct. Salts/mixtures flagged, not silently merged.',
 'herb_specificity':'Exact medicine, part, processing and allowed species checked. Broad herb/genus record excluded from strict set if exact requested material cannot be established, but retained as unresolved coverage. Do not substitute peel for seed, prepared for raw Rehmannia, or generic Dendrobium for D. officinale.',
 'absorption':'No OB/DL hard filter in primary collection; source values annotations only. Optional declared sensitivity OB>=30% and DL>=0.18 only if complete available values; missing values retained as missing.',
 'target_evidence':'Human targets only with explicit database record ID/URL and taxon. Separate experimental, source-database predicted, newly run SEA, newly run SwissTargetPrediction. Database associations never labeled as own prediction. Ambiguous species or unsupported links excluded from strict network.',
 'disease_terms':['constipation','chronic constipation','functional constipation','chronic idiopathic constipation','slow transit constipation'],
 'disease_exclusions':'Do not broaden automatically to IBS, diarrhea, general GI disease or animal-only disease annotations. Retain scores, complete source query denominator and evidence types. No effect-directed score cutoff.',
 'disease_fallback':'If GeneCards unavailable, report explicitly before using public disease/phenotype association resource. HPO constipation is a phenotype-gene fallback with syndromic bias, not a GeneCards equivalent or independent validation.',
 'ppi':'STRING Homo sapiens taxon 9606, supplied nodes only, required_score=400, add_nodes=0. Keep combined score and network type; no added connector genes.',
 'enrichment':'Primary mechanistic analysis on all eligible human drug targets independent of disease intersection. GO BP/MF/CC and KEGG each use complete tested term family and BH adjustment. Background must match assay/prediction target universe if available; otherwise explicitly label annotated human genome exploratory background and limitation. Zero-overlap terms remain in test family where universe allows.',
 'disease_intersection':'Descriptive overlap, not independent validation. Any intersection enrichment explicitly disease-conditioned and exploratory; never a proof of efficacy.',
 'mixture_boundary':'Polysaccharides, whole oils and heterogeneous extracts without a discrete identified structure remain coverage limitations. Do not invent representative SMILES or assign mixture targets.',
 'network_access':'All download processes unset HTTP_PROXY/HTTPS_PROXY/ALL_PROXY case-insensitively; requests trust_env=False; NO_PROXY=*; TLS verification; no proxy fallback, CAPTCHA evasion, login bypass, or automated retries of blocked resources.',
 'stopping_rule':'If core ingredient/target evidence cannot be acquired legitimately, report unidentifiable counts as NA rather than zero, and do not fabricate a full network or enrichment.'
}
path=ROOT/'results'/'analysis_protocol.json';path.write_text(json.dumps(protocol,ensure_ascii=False,indent=2),encoding='utf-8')
(ROOT/'logs'/'protocol_sha256.txt').write_text(hashlib.sha256(path.read_bytes()).hexdigest()+'  results/analysis_protocol.json\n',encoding='utf-8')
print(json.dumps({'herbs':len(rows),'relative_input_sum':protocol['relative_input_sum'],'protocol_sha256':hashlib.sha256(path.read_bytes()).hexdigest()}))
