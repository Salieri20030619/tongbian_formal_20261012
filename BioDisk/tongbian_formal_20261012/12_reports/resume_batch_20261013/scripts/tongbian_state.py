from pathlib import Path
import json,datetime,sys
P=Path('/mnt/BioDisk/wangmingyuan/tongbian_formal_20261012');mode=sys.argv[1];now=datetime.datetime.now(datetime.timezone.utc).isoformat()
p=P/'00_protocol/resume_state_20261013.json';s=json.loads(p.read_text()) if p.exists() else {'project_complete':False,'started_utc':now,'authorization':'User 2026-10-13: project dependencies trusted sources and SSH execution authorized; old approval pending flags superseded'}
s['updated_utc']=now
if mode=='install':
 s['ppcor']={'status':'COMPLETED_VERIFIED','version':'1.1','library':str(P/'environment/R_library_ppcor_20261013'),'exit_code':0,'test':'independent rho/P formula passed'}
 msg='ppcor 1.1从CRAN安装到项目专用库，退出0，独立rho/P测试通过；来源、归档包与sessionInfo已保留。40/100/5对检验集已冻结，尚未标记偏相关完成。'
elif mode=='partial':
 c=json.loads((P/'06_human_selected_metabolites/partial_spearman_v1/COMPUTE_COMPLETE.json').read_text());s['partial_spearman']=c
 msg='偏相关465行计算已完成并通过独立公式和顺序检查；主40对、扩展100对、胆酸5对分家族及模型/情景BH。计算完成标记已写；图文交付核验仍在进行。'
elif mode=='animal':
 s['animal_RNA']=json.loads((P/'07_animal_multiomics/GSE319612_DESeq2_v1/COMPUTE_COMPLETE.json').read_text())
 msg='动物GSE319612 RNA已完成：26个独立样本NC8/MC8/RG10，54532基因、17029预过滤、两对比各16961有效检验；MC-NC名义1247/BH9，RG-MC名义1494/BH4，98个双名义显著反向候选。跨组学配对仍未核实，动物代谢组组别仍冲突。相较人RNA参数：最小组重复过滤由5改8；预先禁用自动count替换，保留默认Cook过滤。'
else:
 s['next_inventory']={'status':'METADATA_FETCHED','manifest':str(P/'01_inventory/resume_20261013/run_manifest_summary.json'),'raw_download_started':False};msg='后续4个BioProject的ENA清单和动物GEO元数据已取回；核心60个WGS样本，FASTQ压缩总计293286068526字节。非交互直连成功，未发现代理环境或7897监听；工具及动物组别进一步审计中。'
p.write_text(json.dumps(s,ensure_ascii=False,indent=2))
with (P/'HANDOFF.md').open('a') as f:f.write('\n## '+now+' / '+mode+'\n'+msg+'\n')
print(msg)
