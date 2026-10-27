"""STRING v12.5 functional PPI for the locked GeneCards intersection; no added nodes."""
from direct_http import *
import pandas as pd, io
import networkx as nx
# Stable v12.5 host timed out at TLS setup; use the official main API, whose
# version endpoint was downloaded as 12.5. This is an endpoint change, never a proxy.
base='https://string-db.org/api/tsv'
df=pd.read_csv(ROOT/'tables/36_GeneCards_primary_target_intersection.tsv',sep='\t',dtype=str,keep_default_na=False)
symbols=sorted(set(df.approved_symbol)); maps=[]
for i in range(0,len(symbols),150):
 label=f'string_map_{i//150:02d}';p=ROOT/'raw'/f'{label}.bin'
 if not p.exists():
  rec,res=get(base+'/get_string_ids',label,method='POST',data_body={'identifiers':'\r'.join(symbols[i:i+150]),'species':9606,'limit':1,'echo_query':1,'caller_identity':'ModelA_reproducible_research'},timeout=(15,60))
  if res is None or res.status_code!=200:raise SystemExit(f'STRING mapping stopped: {rec}')
 maps.append(pd.read_csv(p,sep='\t',dtype=str,keep_default_na=False))
mapping=pd.concat(maps,ignore_index=True)
mapping=mapping[mapping.queryItem.isin(symbols)].copy()
assert set(mapping.ncbiTaxonId)=={'9606'}
mapping['exact_approved_symbol_match']=mapping.queryItem==mapping.preferredName
mapping.to_csv(ROOT/'tables/44_STRING_identifier_mapping.tsv',sep='\t',index=False,encoding='utf-8-sig')
# Only accept exact HGNC symbol matches; alternatives remain explicitly unresolved.
accepted=mapping[mapping.exact_approved_symbol_match].drop_duplicates('queryItem')
# Resolve a renamed official gene only with independent HGNC identifier evidence.
reg=pd.read_csv(ROOT/'raw/hgnc_complete_set.bin',sep='\t',dtype=str,keep_default_na=False)
resolved=[]
for symbol in sorted(set(symbols)-set(accepted.queryItem)):
 rr=reg[reg.symbol==symbol].iloc[0]
 allowed={symbol}|set(rr.prev_symbol.split('|'))|set(rr.alias_symbol.split('|'))
 for prefix,identifier in [('uniprot',rr.uniprot_ids),('ensembl',rr.ensembl_gene_id)]:
  pp=ROOT/'raw'/f'string_{prefix}_{identifier}.bin'
  if not pp.exists():continue
  zz=pd.read_csv(pp,sep='\t',dtype=str,keep_default_na=False)
  zz=zz[zz.preferredName.isin(allowed)]
  if len(zz)==1:
   row=zz.iloc[0].to_dict();row.update(queryItem=symbol,exact_approved_symbol_match=False,mapping_resolution=f'HGNC {prefix} {identifier}; STRING preferred name is HGNC current/previous/alias symbol',resolution_source_raw=str(pp.relative_to(ROOT)))
   resolved.append(row);break
if resolved:accepted=pd.concat([accepted,pd.DataFrame(resolved)],ignore_index=True)
accepted.to_csv(ROOT/'tables/44b_STRING_accepted_verified_mapping.tsv',sep='\t',index=False,encoding='utf-8-sig')
ids=sorted(set(accepted.stringId));fingerprint=hashlib.sha256('\n'.join(ids).encode()).hexdigest()[:12]
label=f'string_final_intersection_network_0400_{fingerprint}';p=ROOT/'raw'/f'{label}.bin'
if not p.exists():
 rec,res=get(base+'/network',label,method='POST',data_body={'identifiers':'\r'.join(ids),'species':9606,'required_score':400,'add_nodes':0,'network_type':'functional','caller_identity':'ModelA_reproducible_research'},timeout=(20,90))
 if res is None or res.status_code!=200:raise SystemExit(f'STRING network stopped: {rec}')
edges=pd.read_csv(p,sep='\t',dtype={'stringId_A':str,'stringId_B':str})
assert set(edges.stringId_A)|set(edges.stringId_B)<=set(ids),'Unrequested node returned'
assert (edges.score>=.4).all()
edges['source']=edges[['stringId_A','stringId_B']].min(axis=1);edges['target']=edges[['stringId_A','stringId_B']].max(axis=1)
edges=edges[edges.source!=edges.target].drop_duplicates(['source','target'])
edges['STRING_version']='12.5';edges['network_type']='functional';edges['source_raw']=str(p.relative_to(ROOT))
edges.to_csv(ROOT/'tables/45_STRING_PPI_edges.tsv',sep='\t',index=False,encoding='utf-8-sig')
G=nx.Graph();G.add_nodes_from(ids)
for r in edges.to_dict('records'):G.add_edge(r['source'],r['target'],weight=float(r['score']),distance=1/float(r['score']))
bet=nx.betweenness_centrality(G,normalized=True,weight='distance')
nodes=accepted[['stringId','queryItem','preferredName']].merge(df,left_on='queryItem',right_on='approved_symbol',how='left')
nodes['degree']=nodes.stringId.map(dict(G.degree()));nodes['weighted_degree']=nodes.stringId.map(dict(G.degree(weight='weight')));nodes['betweenness_weighted_inverse_score']=nodes.stringId.map(bet)
nodes['isolated']=nodes.degree==0
nodes.sort_values(['degree','weighted_degree'],ascending=False).to_csv(ROOT/'tables/46_STRING_PPI_nodes.tsv',sep='\t',index=False,encoding='utf-8-sig')
unmapped=sorted(set(symbols)-set(accepted.queryItem))
summary={'STRING_version':'12.5','request_species':9606,'request_score':400,'add_nodes':0,'network_type':'functional','submitted_intersection_genes':len(symbols),'accepted_exact_symbol_mapped_genes':len(accepted),'unique_STRING_nodes':len(G),'edges_score_at_least_0_4':G.number_of_edges(),'isolated_nodes':sum(d==0 for _,d in G.degree()),'connected_components':nx.number_connected_components(G),'largest_component_nodes':max(map(len,nx.connected_components(G)),default=0),'unresolved_symbols':unmapped,'interpretation':'Functional associations; combined score is database confidence, not binding affinity. Text-mining and disease databases can overlap; PPI is not independent disease validation.','source_url':base+'/network','source_raw':str(p.relative_to(ROOT))}
summary['accepted_hgnc_verified_genes']=summary.pop('accepted_exact_symbol_mapped_genes')
summary['exact_current_symbol_matches']=int(accepted.exact_approved_symbol_match.sum())
summary['HGNC_identifier_resolved_aliases']=len(resolved)
(ROOT/'results/STRING_PPI_summary.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2),encoding='utf-8');print(json.dumps(summary,ensure_ascii=False,indent=2),flush=True)
