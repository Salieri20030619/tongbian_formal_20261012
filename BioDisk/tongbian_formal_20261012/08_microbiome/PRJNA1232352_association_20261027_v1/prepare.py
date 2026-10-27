import pathlib,json,hashlib,csv,math,collections,datetime,urllib.request,xml.etree.ElementTree as ET,concurrent.futures,traceback
O=pathlib.Path(__file__).parent; R=O.parent.parent; P=O.parent/'PRJNA1232352_profiling_20261021_v1'; B=pathlib.Path('/mnt/BioData/bio/wangmingyuan/tongbian_formal_20261012/large_intermediates/PRJNA1232352_profiling_20261021_v1/full')
def now(): return datetime.datetime.now(datetime.timezone.utc).isoformat()
def save(n,d): (O/n).write_text(json.dumps(d,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
def sha(f):
 h=hashlib.sha256()
 with f.open('rb') as z:
  for c in iter(lambda:z.read(8388608),b''): h.update(c)
 return h.hexdigest()
def table(n,header,rows):
 with (O/n).open('w',encoding='utf-8',newline='') as f:
  w=csv.writer(f,delimiter='\t',lineterminator='\n'); w.writerow(header); w.writerows(rows)
def main():
 d=json.loads((P/'FULL_COMPLETE.json').read_text()); s=json.loads((P/'full_status.json').read_text()); samples=json.loads((P/'samples.json').read_text()); inp={x['run']:x for x in samples}
 assert d['status']==s['status']=='COMPLETED_VERIFIED' and (P/'full_exit_code.txt').read_text().strip()=='0'
 assert len(d['results'])==len(inp)==60 and set(inp)=={x['run'] for x in d['results']}
 frozen=json.loads((P/'launch_freeze.json').read_text())
 for key,f in [('pipeline_sha256','profile_pipeline.py'),('config_sha256','config.json'),('samples_sha256','samples.json')]: assert sha(P/f)==frozen[key],f
 records=[]; databases=set()
 for x in d['results']:
  run=x['run']; a=B/run; ir=json.loads((a/'input_receipt.json').read_text()); expected=inp[run]
  assert x==json.loads((a/'VALIDATED.json').read_text())
  assert x['reads']==x['metaphlan']['reads_processed']==ir['reads']==2*expected['retained_pairs']
  assert ir['pairs']==expected['retained_pairs'] and not ir['pilot_prefix'] and x['group']==expected['group']
  assert s['samples'][run]['status']=='VERIFIED'
  for n in ['metaphlan_exit_code.txt','humann_exit_code.txt']: assert (a/n).read_text().strip()=='0'
  databases.add(x['metaphlan']['header'][0])
  files=[(a/'metaphlan.tsv',x['profile_sha256'],None)]+[(a/'humann'/(run+'_'+k+'.tsv'),v['sha256'],v['bytes']) for k,v in x['humann'].items()]
  for f,h,size in files:
   assert f.exists() and (size is None or f.stat().st_size==size) and sha(f)==h,str(f)
   records.append([run,str(f),f.stat().st_size,h])
 assert len(databases)==1
 table('source_hashes.tsv',['run','file','bytes','sha256'],records)
 save('PROFILING_ACCEPTANCE.json',{'status':'PASS','utc':now(),'samples':60,'outputs_sha256_verified':len(records),'groups':dict(collections.Counter(x['group'] for x in samples)),'database_headers':list(databases),'input_receipts_and_read_counts_match':True,'baseline_hashes_unchanged':True,'all_stage_exit_codes':0,'raw_fastq_rehashed_this_run':False,'scope':'profiling output integrity and recorded input provenance; no rerun of upstream FASTQ validation','source_completion_sha256':sha(P/'FULL_COMPLETE.json')})
 print('PROFILING_ACCEPTANCE PASS',flush=True)
 rows=list(csv.DictReader((R/'01_inventory/resume_20261013/PRJNA1232352_runs.tsv').open(),delimiter='\t')); meta={x['run_accession']:x for x in rows}; ids=sorted(inp)
 assert len(meta)==60 and set(meta)==set(ids) and len({x['sample_accession'] for x in rows})==60 and len({x['sample_alias'] for x in rows})==60
 for run in ids: assert meta[run]['sample_alias'].rstrip('0123456789')==inp[run]['group']
 table('metadata.tsv',['sample_id','group','biosample','sample_alias'],[[run,inp[run]['group'],meta[run]['sample_accession'],meta[run]['sample_alias']] for run in ids])
 (O/'biosamples').mkdir()
 def fetch(run):
  acc=meta[run]['sample_accession']; u='https://www.ebi.ac.uk/ena/browser/api/xml/'+acc
  data=urllib.request.urlopen(u,timeout=45).read(); root=ET.fromstring(data); sample=root.find('.//SAMPLE'); assert sample is not None
  assert acc in [x.text for x in root.findall('.//EXTERNAL_ID')]+[x.text for x in root.findall('.//PRIMARY_ID')],acc
  (O/'biosamples'/ (acc+'.xml')).write_bytes(data)
  return run,[(a.findtext('TAG'),a.findtext('VALUE')) for a in root.findall('.//SAMPLE_ATTRIBUTE')]
 with concurrent.futures.ThreadPoolExecutor(max_workers=4) as ex: attrs=dict(ex.map(fetch,ids))
 table('biosample_attributes.tsv',['sample_id','attribute','value'],[[run,k,v] for run in ids for k,v in attrs[run]])
 cov={run:{k:v for k,v in vals if any(t in k.lower() for t in ['age','sex','gender','bmi','subject','individual','patient','treatment']) and 'package' not in k.lower()} for run,vals in attrs.items()}
 save('metadata_audit.json',{'samples':60,'unique_biosamples':60,'unique_aliases':60,'alias_group_matches':60,'candidate_covariates':cov,'sources':'ENA BioSample XML archived per sample; frozen run manifest','cross_cohort_join':False})
 print('METADATA PASS',json.dumps({k:v for k,v in cov.items() if v}),flush=True)
 matrices={k:{} for k in ['species_percent','pathabundance_raw','pathcoverage','pathabundance_stratified_raw']}; special=[]
 for run in ids:
  a=B/run
  with (a/'metaphlan.tsv').open() as f:
   for line in f:
    if line.startswith('#'): continue
    t=line.rstrip('\n').split('\t'); key=t[0]
    if key.split('|')[-1].startswith('s__'):
     v=float(t[2]); assert math.isfinite(v) and v>=0
     assert run not in matrices['species_percent'].setdefault(key,{})
     matrices['species_percent'][key][run]=v
  for kind in ['pathabundance','pathcoverage']:
   with (a/'humann'/(run+'_'+kind+'.tsv')).open() as f:
    for line in f:
     if line.startswith('#'): continue
     key,val=line.rstrip('\n').split('\t'); val=float(val); assert math.isfinite(val) and val>=0
     if key.split('|')[0] in ['UNMAPPED','UNINTEGRATED','UNGROUPED']:
      special.append([run,kind,key,val]); continue
     if '|' in key:
      if kind=='pathcoverage': continue
      target='pathabundance_stratified_raw'
     else: target='pathabundance_raw' if kind=='pathabundance' else 'pathcoverage'
     assert run not in matrices[target].setdefault(key,{})
     matrices[target][key][run]=val
 table('humann_special_features.tsv',['sample_id','type','feature','value'],special)
 denom={run:sum(v.get(run,0) for v in matrices['pathabundance_raw'].values()) for run in ids}; assert min(denom.values())>0
 matrices['pathabundance_cpm']={k:{run:1e6*v.get(run,0)/denom[run] for run in ids} for k,v in matrices['pathabundance_raw'].items()}
 matrices['pathabundance_stratified_cpm']={k:{run:1e6*v.get(run,0)/denom[run] for run in ids} for k,v in matrices['pathabundance_stratified_raw'].items()}
 summaries={}
 for name,mat in matrices.items():
  file=name+'.tsv'; table(file,['feature']+ids,[[k]+[format(mat[k].get(run,0),'.17g') for run in ids] for k in sorted(mat)])
  with (O/file).open() as f:
   rd=csv.reader(f,delimiter='\t'); assert next(rd)==['feature']+ids
   count=0
   for row in rd:
    assert len(row)==61; assert all(float(v)==mat[row[0]].get(run,0) for run,v in zip(ids,row[1:])); count+=1
   assert count==len(mat)
  summaries[name]={'features':len(mat),'samples':60,'sha256':sha(O/file),'bytes':(O/file).stat().st_size}
 sums={run:sum(v.get(run,0) for v in matrices['species_percent'].values()) for run in ids}; assert all(abs(v-100)<0.01 for v in sums.values())
 assert all(abs(sum(v[run] for v in matrices['pathabundance_cpm'].values())-1e6)<0.01 for run in ids)
 save('MATRIX_ACCEPTANCE.json',{'status':'PASS','utc':now(),'matrices':summaries,'species_sums_min_max':[min(sums.values()),max(sums.values())],'pathway_normalization':'CPM over named community pathways, special features excluded and retained separately; stratified rows use same community denominator','coverage_not_normalized':True,'zero_fill':'absent features in successful complete profiles only','original_features_preserved':True,'no_sample_exclusion':True,'roundtrip_all_cells_verified':True})
 save('stage1_status.json',{'status':'COMPLETED_VERIFIED','utc':now(),'scope':'profiling acceptance; sample metadata audit; species and pathway matrices','association_status':'NOT_YET_RUN'})
 print('MATRIX_ACCEPTANCE PASS',json.dumps(summaries),flush=True)
if __name__=='__main__':
 try: main(); (O/'stage1_exit_code.txt').write_text('0\n')
 except BaseException:
  (O/'stage1_exit_code.txt').write_text('1\n'); traceback.print_exc(); raise
