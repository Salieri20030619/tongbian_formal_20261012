"""Lightweight preparation audit; no installs, downloads, or model fitting."""
import os,sys,json,socket,platform,shutil,subprocess,datetime,importlib.metadata
from pathlib import Path
root=Path('/mnt/BioDisk/wangmingyuan/tongbian_formal_20261012')
assert socket.gethostname()=='osv01'
assert json.loads((root/'tongbian_project_identity.json').read_text())['project_id']=='tongbian_formal_20261012'
dest=root/'environment';dest.mkdir(exist_ok=True)
stamp=datetime.datetime.now(datetime.timezone.utc).isoformat()
mem={x.split(':')[0]:x.split(':')[1].strip() for x in Path('/proc/meminfo').read_text().splitlines() if x.startswith(('MemTotal:','MemAvailable:','MemFree:'))}
exe=['Rscript','python3','conda','metaphlan','humann','megahit','prodigal','diamond','run_dbcan','fastp','fastqc','samtools','prefetch','fasterq-dump','vina','obabel','gmx','gmx_MMPBSA','nvidia-smi']
audit={'time_utc':stamp,'hostname':socket.gethostname(),'user':os.environ.get('USER'),'project_root':str(root),'platform':platform.platform(),'cpu_logical':os.cpu_count(),'memory':mem,'project_disk':dict(shutil.disk_usage(root)._asdict()),'on_path':{k:shutil.which(k) for k in exe},'formal_computation_started':False}
(dest/'server_capacity.json').write_text(json.dumps(audit,indent=2))
print(json.dumps(audit,indent=2),flush=True)
rcode="""cat(R.version.string, '\\n')
cat('LIBPATHS:', .libPaths(), '\\n')
p <- c('DESeq2','limma','fgsea','GSVA','decoupleR','dorothea','viper','ppcor','MetaboAnalystR','Maaslin2','clusterProfiler','org.Hs.eg.db','org.Mm.eg.db','AnnotationDbi','biomaRt','jsonlite')
ip <- installed.packages()
out <- data.frame(package=p,version=NA_character_,load_status='NOT_TESTED')
for(i in seq_along(p)) {
 if(p[i] %in% rownames(ip)) {
  out$version[i]<-ip[p[i],'Version']
  if(p[i] %in% c('DESeq2','limma','fgsea','GSVA','decoupleR','dorothea','viper','ppcor','jsonlite')) {
   out$load_status[i] <- if(requireNamespace(p[i],quietly=TRUE)) 'PASS' else 'FAIL'
  }
 } else out$load_status[i]<-'MISSING'
}
write.table(out,stdout(),sep='\\t',quote=FALSE,row.names=FALSE)
"""
(dest/'check_R_packages.R').write_text(rcode)
r_candidates=[Path('/usr/bin/Rscript'),Path.home()/'miniconda3/envs/R_4.5.2/bin/Rscript',Path.home()/'miniconda3/envs/r_4.4/bin/Rscript']
records=[]
for i,p in enumerate(r_candidates):
 if not p.is_file():continue
 print('CHECK_R',str(p),flush=True)
 try:
  r=subprocess.run([str(p),'--vanilla',str(dest/'check_R_packages.R')],capture_output=True,text=True,timeout=75)
  log=dest/f'R_environment_{i}.txt';log.write_text(r.stdout+'\nSTDERR:\n'+r.stderr)
  records.append({'executable':str(p),'exit_code':r.returncode,'log':str(log)})
  print(r.stdout[-3000:],flush=True)
 except subprocess.TimeoutExpired:
  records.append({'executable':str(p),'status':'TIMEOUT_75s'})
(dest/'R_runtime_index.json').write_text(json.dumps(records,indent=2))
py_candidates=[Path(sys.executable),Path.home()/'miniconda3/bin/python']
pycode="import sys,json,importlib.metadata as m; d={'python':sys.version,'executable':sys.executable,'packages':{}}\nfor p in ['pandas','numpy','scipy','pydeseq2','rdkit','torch','dgl','networkx','requests']:\n try:d['packages'][p]=m.version(p)\n except m.PackageNotFoundError:d['packages'][p]=None\nprint(json.dumps(d))"
pr=[]
for p in py_candidates:
 if p.is_file():
  r=subprocess.run([str(p),'-c',pycode],capture_output=True,text=True,timeout=20)
  pr.append({'candidate':str(p),'exit_code':r.returncode,'output':r.stdout,'stderr':r.stderr})
(dest/'python_runtimes.json').write_text(json.dumps(pr,indent=2))
# Resource discovery is by filename only. Do not read unrelated project data.
refs=[]
for base in [Path.home()/'ref',Path.home()/'R_libs',Path.home()/'miniconda3/envs/R_4.5.2/lib/R/library']:
 if not base.is_dir():continue
 for curr,dirs,files in os.walk(base):
  dep=len(Path(curr).relative_to(base).parts)
  if dep>=4:dirs[:]=[]
  for n in files:
   if (n.endswith('.gmt') and ('c2.cp.' in n.lower() or 'go.bp.' in n.lower())) or ('dorothea' in n.lower() and n.endswith(('.rds','.RData','.tsv','.csv'))):
    p=Path(curr)/n;refs.append({'path':str(p),'bytes':p.stat().st_size})
(dest/'reference_file_candidates.json').write_text(json.dumps(refs,indent=2))
print('ENVIRONMENT_AUDIT_COMPLETE',len(records),'R runtimes; reference candidates',len(refs),flush=True)
