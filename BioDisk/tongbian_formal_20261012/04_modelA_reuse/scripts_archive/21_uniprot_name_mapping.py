"""Name-only mapping is explicitly an inferred human counterpart, never source-species proof."""
from pathlib import Path
import pandas as pd,re,json,unicodedata
ROOT=Path(__file__).resolve().parents[1]
def norm(s):
 return re.sub(r'\s+',' ',unicodedata.normalize('NFKC',str(s)).lower().replace('’',"'")).strip()
p=pd.read_csv(ROOT/'raw'/'uniprot_reviewed_human_proteome.bin',sep='\t',dtype=str).fillna('')
assert set(p['Organism (ID)'])=={'9606'}
lookup={}
for row in p.to_dict('records'):
 full=row['Protein names'];names=[full.split(' (')[0].split(' [')[0]]
 # Full parenthesized alternative names, excluding enzyme-class numbers and processing markers.
 names += [m.group(1) for m in re.finditer(r'\(([^()]+)\)',full) if not m.group(1).startswith(('EC ','Fragment','Fragments'))]
 for n in set(names):lookup.setdefault(norm(n),[]).append(row)
t=pd.read_csv(ROOT/'tables'/'06_tcmsp_target_rows_UNCLASSIFIED.tsv',sep='\t',dtype=str)[['target_ID','target_name']].drop_duplicates()
rows=[]
for r in t.to_dict('records'):
 matches={x['Entry']:x for x in lookup.get(norm(r['target_name']),[])}
 base=dict(r,mapping_type='exact_normalized_protein_name_to_reviewed_human_proteome',original_source_taxon='UNKNOWN',mapping_is_species_proof=False,match_count=len(matches))
 if not matches:rows.append(dict(base,status='unmapped'));continue
 for x in matches.values():rows.append(dict(base,status='unique_human_counterpart' if len(matches)==1 else 'ambiguous_human_counterpart',uniprot=x['Entry'],human_gene_names=x['Gene Names'],human_ncbi_gene_ids=x['GeneID'],mapped_taxon='9606',mapping_source='https://rest.uniprot.org/uniprotkb/'+x['Entry']))
out=pd.DataFrame(rows);out.to_csv(ROOT/'tables'/'19_tcmsp_uniprot_name_mapping_INFERRED.tsv',sep='\t',index=False,encoding='utf-8-sig')
print(out[['target_ID','status']].drop_duplicates().status.value_counts().to_dict())
print('proteome',len(p),'target_names',len(t))
