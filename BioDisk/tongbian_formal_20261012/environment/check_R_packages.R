cat(R.version.string, '\n')
cat('LIBPATHS:', .libPaths(), '\n')
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
write.table(out,stdout(),sep='\t',quote=FALSE,row.names=FALSE)
