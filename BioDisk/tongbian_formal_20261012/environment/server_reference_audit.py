import os,subprocess,json,hashlib,datetime,socket
from pathlib import Path
root=Path('/mnt/BioDisk/wangmingyuan/tongbian_formal_20261012');assert socket.gethostname()=='osv01'
refs=root/'03_inputs/references'
R=Path.home()/'miniconda3/envs/R_4.5.2/bin/Rscript'
code="""suppressPackageStartupMessages(library(dorothea))
for(n in c('dorothea_hs','dorothea_mm')) {
 e<-new.env();data(list=n,package='dorothea',envir=e);d<-get(n,envir=e)
 d<-d[d$confidence %in% c('A','B','C'),,drop=FALSE]
 stopifnot(all(c('tf','target','mor','confidence') %in% colnames(d)))
 write.table(d,file=file.path(commandArgs(TRUE)[1],paste0(n,'_ABC_v',packageVersion('dorothea'),'.tsv')),sep='\\t',quote=FALSE,row.names=FALSE)
 cat(n,'rows',nrow(d),'TF',length(unique(d$tf)),'targets',length(unique(d$target)),'\\n')
}
"""
p=root/'environment/export_dorothea_ABC.R';p.write_text(code)
r=subprocess.run([str(R),'--vanilla',str(p),str(refs)],capture_output=True,text=True,timeout=50)
(root/'environment/dorothea_resource_export.log').write_text(r.stdout+'\n'+r.stderr)
assert r.returncode==0,r.stderr
print(r.stdout,flush=True)
found=[];skip={'miniconda3','miniconda3.bak.20260312','R','R_libs','.git','.cache','.local','.conda','node_modules','snap','__pycache__'}
for cur,dirs,files in os.walk(Path.home()):
 p=Path(cur);dep=len(p.relative_to(Path.home()).parts)
 dirs[:]=[d for d in dirs if d not in skip and not d.startswith('.')]
 if dep>=6:dirs[:]=[]
 for n in files:
  if n.endswith('.gmt') and ('c2.cp.' in n or 'c5.go.bp.' in n):
   x=p/n;found.append({'path':str(x),'bytes':x.stat().st_size})
(root/'environment/msigdb_file_candidates.json').write_text(json.dumps(found,indent=2))
print('MSIGDB_CANDIDATES',json.dumps(found,ensure_ascii=False),flush=True)
rs=[]
for p in refs.iterdir():
 if p.is_file():rs.append({'file':p.name,'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'origin':'installed dorothea 1.22.0; A-C only; species-specific' if p.name.startswith('dorothea') else 'Windows shared MSigDB file; upstream complete-release verification pending'})
(root/'01_inventory/reference_manifest.json').write_text(json.dumps(rs,indent=2))
print('REFERENCE_AUDIT_DONE')
