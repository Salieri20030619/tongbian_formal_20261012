import json,hashlib,zipfile,os,socket,csv,gzip,datetime
from pathlib import Path
root=Path('/mnt/BioDisk/wangmingyuan/tongbian_formal_20261012').resolve()
assert socket.gethostname()=='osv01'
assert json.loads((root/'tongbian_project_identity.json').read_text())['project_id']=='tongbian_formal_20261012'
bundle=root/'00_protocol/server_input_bundle.zip'
expected='fc525314358e57cbfe03c4a64cd2c3ffa98d6948fdc3d5e3ee8f17e7be01ed54'
assert hashlib.sha256(bundle.read_bytes()).hexdigest()==expected,'bundle SHA256 mismatch'
with zipfile.ZipFile(bundle) as z:
 for info in z.infolist():
  p=(root/info.filename).resolve()
  assert p.is_relative_to(root) and not info.is_dir()
  data=z.read(info)
  if p.exists():assert p.read_bytes()==data, 'Refuse to overwrite existing different file: '+str(p)
  else:p.parent.mkdir(parents=True,exist_ok=True);p.write_bytes(data)
manifest=json.loads((root/'01_inventory/staged_input_manifest.json').read_text())
for x in manifest:
 p=root/x['server_relative_path']
 assert p.stat().st_size==x['bytes'] and hashlib.sha256(p.read_bytes()).hexdigest()==x['sha256'],str(p)
pilot=root/'03_inputs/GSE245885_pilot'
with (pilot/'metadata/samples.tsv').open() as f:meta=list(csv.DictReader(f,delimiter='\t'))
ids=None;matrix={};samples=[]
for m in meta:
 sample=m['sample'];p=pilot/'raw'/m['supplement_file']
 with gzip.open(p,'rt') as f:
  rd=csv.DictReader(f,delimiter='\t')
  assert rd.fieldnames==['ensembl_id','gene_name',sample+'.count',sample+'.FPKM']
  these_ids=[];vals=[];symbols=[]
  for r in rd:
   v=r[sample+'.count'];assert v.isdigit(),(sample,v)
   vals.append(int(v));these_ids.append(r['ensembl_id']);symbols.append(r['gene_name'])
 assert len(set(these_ids))==len(these_ids)
 if ids is None:ids=these_ids;gene_symbols=symbols
 else:assert ids==these_ids and gene_symbols==symbols
 matrix[sample]=vals
 samples.append({'gsm':m['gsm'],'sample':sample,'condition':m['condition'],'count_column':sample+'.count','n_genes':len(vals),'count_sum':sum(vals),'input_relative_path':p.relative_to(root).as_posix(),'count_integer_nonnegative':True})
with gzip.open(pilot/'results/counts_all_genes.tsv.gz','rt') as f:
 rd=csv.DictReader(f,delimiter='\t');key=rd.fieldnames[0];n=0
 for i,r in enumerate(rd):
  assert r[key]==ids[i]
  for sample in matrix:assert int(r[sample])==matrix[sample][i]
  n+=1
 assert n==len(ids)
with (root/'02_metadata/GSE245885_samples_verified.tsv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=samples[0],delimiter='\t');w.writeheader();w.writerows(samples)
prefilter=sum(sum(matrix[s][i]>=10 for s in matrix)>=5 for i in range(len(ids)))
report={'time_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'hostname':socket.gethostname(),'project_root':str(root),'bundle_sha256':expected,'staged_files_verified':len(manifest),'verified_bytes':sum(x['bytes'] for x in manifest),'GSE245885':{'raw_tables':len(meta),'shape':[len(ids),len(meta)],'nonnegative_integer_counts':True,'identical_gene_order':True,'matches_pilot_matrix_exactly':True,'prefilter_genes':prefilter,'groups':{'Control':5,'STC':5}},'modelA_rerun':False,'formal_analysis_started':False,'operation':'input transfer integrity and numeric validation only'}
(root/'01_inventory/server_staged_input_verification.json').write_text(json.dumps(report,indent=2))
print(json.dumps(report,indent=2))
