"""Import the user's manually exported complete GeneCards constipation table, never scrape."""
from pathlib import Path
import pandas as pd,argparse,json,hashlib,re,shutil,csv,io
from datetime import datetime,timezone
ROOT=Path(__file__).resolve().parents[1]
parser=argparse.ArgumentParser();parser.add_argument('input');parser.add_argument('--expected-total',type=int);parser.add_argument('--query',default='constipation',choices=['constipation']);args=parser.parse_args()
path=Path(args.input).resolve()
if not path.is_file():raise SystemExit('User export file does not exist')
preface=[]
if path.suffix.lower() in ['.xlsx','.xls']:df=pd.read_excel(path,dtype=str,keep_default_na=False)
else:
 text=path.read_text(encoding='utf-8-sig');lines=text.splitlines();header=next((i for i,line in enumerate(lines) if line.startswith('Symbol,') and 'Relevance Score' in line),None)
 if header is None:raise SystemExit('Cannot identify exact GeneCards CSV header')
 preface=lines[:header]
 if not any('Search results for constipation' in line for line in preface):raise SystemExit('Query title does not confirm constipation')
 df=pd.read_csv(io.StringIO('\n'.join(lines[header:])),dtype=str,keep_default_na=False)
df=df.fillna(''); norm=lambda x:re.sub(r'[^a-z]','',str(x).lower())
columns={norm(x):x for x in df.columns}
sym=next((columns[x] for x in ['symbol','genesymbol','genename'] if x in columns),None)
score=next((columns[x] for x in ['relevancescore','score'] if x in columns),None)
if sym is None or score is None:raise SystemExit('Missing unambiguous Symbol / Relevance Score columns; inspect original file rather than guessing')
all_csv_records=df.copy()
footer_mask=df[sym].str.startswith('Copyright LifeMap Sciences Inc.') & df.drop(columns=[sym]).eq('').all(axis=1)
excluded_footer=df[footer_mask].copy()
df=df[~footer_mask].copy()
if df[sym].eq('').any():raise SystemExit('Blank gene symbol requires manual row review')
if args.expected_total is not None and len(df)!=args.expected_total:raise SystemExit(f'Export row count {len(df)} does not match reported total {args.expected_total}; retain but do not promote incomplete table')
registry=pd.read_csv(ROOT/'raw'/'hgnc_complete_set.bin',sep='\t',dtype=str).fillna('');approved=registry[registry.status=='Approved']
lookup={}
for r in approved.to_dict('records'):
 for kind,values in [('approved', [r['symbol']]),('previous',r['prev_symbol'].split('|')),('alias',r['alias_symbol'].split('|'))]:
  for value in values:
   if value:lookup.setdefault(value,[]).append((kind,r))
out=[]
for ix,r in enumerate(df.to_dict('records'),1):
 symbol=r[sym].strip();matches=lookup.get(symbol,[]);exact=[x for x in matches if x[0]=='approved'];matches=exact or matches
 unique={x[1]['hgnc_id']:x for x in matches}
 rr={'original_row':ix,'input_symbol':symbol,'input_name':r.get('Name',''),'genecards_type':r.get('Type',''),'knowledge':r.get('Knowledge',''),'source_score_unfiltered':r[score],'source_query':args.query,'source':'GeneCards_user_manual_export','source_taxon':'9606','mapping_status':'unmapped','hgnc_id':'','approved_symbol':'','ncbi_gene_id':''}
 if len(unique)==1:
  kind,g=next(iter(unique.values()));rr.update(mapping_status='unique_'+kind,hgnc_id=g['hgnc_id'],approved_symbol=g['symbol'],ncbi_gene_id=g['entrez_id'],hgnc_locus_group=g['locus_group'],hgnc_locus_type=g['locus_type'])
 elif len(unique)>1:rr['mapping_status']='ambiguous'
 out.append(rr)
destination=ROOT/'raw'/'user_genecards';destination.mkdir(exist_ok=True);saved=destination/path.name
if saved.exists() and saved.read_bytes()!=path.read_bytes():raise SystemExit('An unequal previous manual export exists; use a new filename')
if path!=saved:shutil.copy2(path,saved)
all_csv_records.to_csv(ROOT/'tables'/'GeneCards_complete_original_rows.tsv',sep='\t',index=False,encoding='utf-8-sig')
excluded_footer.to_csv(ROOT/'tables'/'GeneCards_excluded_copyright_footer.tsv',sep='\t',index=False,encoding='utf-8-sig')
df.to_csv(ROOT/'tables'/'GeneCards_valid_gene_rows.tsv',sep='\t',index=False,encoding='utf-8-sig')
mapped=pd.DataFrame(out);mapped.to_csv(ROOT/'tables'/'GeneCards_human_gene_mapping.tsv',sep='\t',index=False,encoding='utf-8-sig')
scores=pd.to_numeric(df[score],errors='coerce')
metadata={'imported_utc':datetime.now(timezone.utc).isoformat(),'source_file':str(saved.relative_to(ROOT)),'source_library_id_user_supplied':'libfile_2403b75d1d4081918003c9577a90ecaf','local_origin':str(path),'sha256':hashlib.sha256(saved.read_bytes()).hexdigest(),'bytes':saved.stat().st_size,'preface_and_copyright':preface,'query':args.query,'expected_total_user_reported':args.expected_total,'actual_rows':len(df),'unique_input_symbols':df[sym].nunique(),'duplicate_symbol_rows':int(df[sym].duplicated().sum()),'missing_symbol_rows':int(df[sym].eq('').sum()),'type_distribution':df['Type'].value_counts().to_dict() if 'Type' in df else {},'numeric_score_min':float(scores.min()),'numeric_score_max':float(scores.max()),'missing_or_nonnumeric_score_rows':int(scores.isna().sum()),'score_filter':'none','mapping_counts':mapped.mapping_status.value_counts().to_dict(),'unique_mapped_ncbi_genes':mapped.loc[mapped.mapping_status.str.startswith('unique_')&mapped.ncbi_gene_id.ne(''),'ncbi_gene_id'].nunique(),'interpretation':'关键词相关集合，不是便秘因果金标准；Relevance Score不是效应大小或概率；原表未上传外部服务。'}
metadata.update(parsed_csv_records_including_footer=len(all_csv_records),excluded_copyright_footer_records=len(excluded_footer),footer_rule='only explicit Copyright LifeMap row with every other column blank')
(ROOT/'results'/'GeneCards_manual_import_metadata.json').write_text(json.dumps(metadata,ensure_ascii=False,indent=2),encoding='utf-8');print(json.dumps(metadata,ensure_ascii=False,indent=2))
