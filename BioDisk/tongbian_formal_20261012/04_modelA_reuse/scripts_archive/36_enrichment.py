"""Full-family hypergeometric ORA with traceable effective gene universes and BH adjustment."""
from pathlib import Path
import pandas as pd,json,re,sys
from functools import lru_cache
from scipy.stats import hypergeom
from statsmodels.stats.multitest import multipletests
from datetime import datetime,timezone
ROOT=Path(__file__).resolve().parents[1]
ontology_text=(ROOT/'raw'/'go_basic_ontology.bin').read_text(encoding='utf-8');terms={};alt={}
for block in ontology_text.split('[Term]')[1:]:
 mid=re.search(r'^id: (GO:\d+)',block,re.M);name=re.search(r'^name: (.+)',block,re.M);ns=re.search(r'^namespace: (.+)',block,re.M)
 if not mid or not name or not ns or 'is_obsolete: true' in block:continue
 ident=mid[1];family={'biological_process':'GO_BP','molecular_function':'GO_MF','cellular_component':'GO_CC'}.get(ns[1])
 if not family:continue
 parents=set(re.findall(r'^is_a: (GO:\d+)',block,re.M))|set(re.findall(r'^relationship: part_of (GO:\d+)',block,re.M))
 terms[ident]={'term_id':ident,'term_name':name[1],'family':family,'parents':parents}
 for x in re.findall(r'^alt_id: (GO:\d+)',block,re.M):alt[x]=ident
@lru_cache(None)
def ancestors(term):
 term=alt.get(term,term)
 if term not in terms:return frozenset()
 out={term}
 for p in terms[term]['parents']:out.update(ancestors(p))
 return frozenset(out)
g=pd.read_csv(ROOT/'raw'/'uniprot_human_go.bin',sep='\t',dtype=str,keep_default_na=False)
assert set(g['Organism (ID)'])=={'9606'}
term_genes={t:set() for t in terms}
go_unknown=set()
for r in g.to_dict('records'):
 geneids=set(re.findall(r'\d+',r['GeneID']));direct=re.findall(r'GO:\d+',r['Gene Ontology IDs'])
 if not geneids:continue
 for go in direct:
  aa=ancestors(go)
  if not aa:go_unknown.add(go)
  for t in aa:term_genes[t].update(geneids)
kegg_names={}
for line in (ROOT/'raw'/'kegg.bin').read_text(encoding='utf-8').splitlines():
 ident,name=line.split('\t',1);kegg_names[ident.removeprefix('path:')]=name
for line in (ROOT/'raw'/'kegg_human_gene_pathways.bin').read_text(encoding='utf-8').splitlines():
 gene,term=line.split('\t');gid=gene.removeprefix('hsa:');tid=term.removeprefix('path:')
 terms[tid]={'term_id':tid,'term_name':kegg_names.get(tid,tid),'family':'KEGG','parents':set()}
 term_genes.setdefault(tid,set()).add(gid)
targets=pd.read_csv(ROOT/'tables'/'35_primary_human_targets.tsv',sep='\t',dtype=str,keep_default_na=False)
inter=pd.read_csv(ROOT/'tables'/'36_GeneCards_primary_target_intersection.tsv',sep='\t',dtype=str,keep_default_na=False)
reg=pd.read_csv(ROOT/'raw'/'hgnc_complete_set.bin',sep='\t',dtype=str,keep_default_na=False)
reg=reg[(reg.status=='Approved')&(reg.locus_group=='protein-coding gene')&reg.entrez_id.ne('')]
symbols=dict(zip(reg.entrez_id,reg.symbol));D=set(targets.ncbi_gene_id);I=set(inter.ncbi_gene_id);H=set(reg.entrez_id)
families=['GO_BP','GO_MF','GO_CC','KEGG'];annotation_genes={f:set().union(*(term_genes[t] for t in terms if terms[t]['family']==f)) for f in families}
analyses=[('PRIMARY_intersection_vs_all_formula_targets',I,D),('SENSITIVITY_intersection_vs_human_protein_coding',I,H),('CONTROL_all_formula_targets_vs_human_protein_coding',D,H)]
rows=[];universes=[];coverage=[]
for name,query,universe in analyses:
 for fam in families:
  U=universe&annotation_genes[fam];Q=query&U;N=len(U);n=len(Q)
  for gid in sorted(universe|query):universes.append({'analysis':name,'family':fam,'ncbi_gene_id':gid,'symbol':symbols.get(gid,''),'in_supplied_universe':gid in universe,'has_family_annotation':gid in annotation_genes[fam],'in_effective_universe':gid in U,'in_supplied_query':gid in query,'in_effective_query':gid in Q})
  tested=0
  for tid,meta in terms.items():
   if meta['family']!=fam:continue
   bg=term_genes[tid]&U;K=len(bg)
   if K==0:continue
   hits=bg&Q;k=len(hits);tested+=1
   p=float(hypergeom.sf(k-1,N,K,n)) if n else 1.0
   rows.append({'analysis':name,'family':fam,'term_id':tid,'term_name':meta['term_name'],'N_background':N,'n_query':n,'K_term_background':K,'k_overlap':k,'raw_p':p,'fold_enrichment':(k/n)/(K/N) if n else 0.,'gene_ratio':f'{k}/{n}','background_ratio':f'{K}/{N}','overlap_ncbi_gene_ids':';'.join(sorted(hits,key=int)),'overlap_symbols':';'.join(symbols.get(x,x) for x in sorted(hits,key=int))})
  coverage.append({'analysis':name,'family':fam,'supplied_universe_count':len(universe),'effective_universe_count':N,'supplied_query_count':len(query),'effective_query_count':n,'tested_terms_K_at_least_1':tested,'term_filter':'none beyond nonzero universe annotations; zero query overlap retained'})
df=pd.DataFrame(rows)
df['p_adjust_BH_within_family']=df.groupby(['analysis','family']).raw_p.transform(lambda x:multipletests(x,method='fdr_bh')[1])
df['p_adjust_BH_all_families']=df.groupby('analysis').raw_p.transform(lambda x:multipletests(x,method='fdr_bh')[1])
df['significant_primary_criterion']=df.p_adjust_BH_all_families<.05
df=df.sort_values(['analysis','p_adjust_BH_all_families','raw_p','family','term_id'])
df.to_csv(ROOT/'tables'/'40_GO_KEGG_ALL_tested_results.tsv',sep='\t',index=False,encoding='utf-8-sig')
pd.DataFrame(universes).to_csv(ROOT/'tables'/'41_enrichment_gene_universes.tsv',sep='\t',index=False,encoding='utf-8-sig')
pd.DataFrame(coverage).to_csv(ROOT/'tables'/'42_enrichment_effective_backgrounds.tsv',sep='\t',index=False,encoding='utf-8-sig')
df[df.significant_primary_criterion].to_csv(ROOT/'tables'/'43_GO_KEGG_BH_significant_results.tsv',sep='\t',index=False,encoding='utf-8-sig')
summary={'time_utc':datetime.now(timezone.utc).isoformat(),'method':'upper-tail hypergeometric; GO ancestors through is_a/part_of; complete K>=1 term families including k=0; primary significance BH across all GO/KEGG terms per analysis, additional family BH retained','primary_question':'Within identifiable formula targets, which functions are relatively represented among the GeneCards keyword-overlap subset? Disease-conditioned descriptive analysis, not independent validation.','sensitivity_question':'Compare the same overlap subset with the annotated human protein-coding genome; a different question that can reflect generic drug-target/database ascertainment bias.','control_question':'All formula targets versus the same human background, independent of GeneCards filtering.','score_cutoff':'none','go_data_version':re.search(r'^data-version: (.+)',ontology_text,re.M)[1],'go_relations_propagated':['is_a','part_of'],'go_annotation_source':'UniProt reviewed Homo sapiens entries; all supplied GO IDs, no evidence-code filter possible in this export','kegg_version':(ROOT/'raw'/'kegg_database_version.bin').read_text(encoding='utf-8'),'unknown_or_obsolete_uniprot_go_ids':sorted(go_unknown),'summary':df.groupby(['analysis','family']).agg(tested_terms=('term_id','size'),significant_BH_all_families=('significant_primary_criterion','sum'),minimum_raw_p=('raw_p','min'),minimum_BH_all_families=('p_adjust_BH_all_families','min')).reset_index().to_dict('records')}
(ROOT/'results'/'enrichment_summary.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2),encoding='utf-8');print(json.dumps(summary['summary'],ensure_ascii=False,indent=2))
