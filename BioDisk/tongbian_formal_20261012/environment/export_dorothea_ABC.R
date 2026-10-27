suppressPackageStartupMessages(library(dorothea))
for(n in c('dorothea_hs','dorothea_mm')) {
 e<-new.env();data(list=n,package='dorothea',envir=e);d<-get(n,envir=e)
 d<-d[d$confidence %in% c('A','B','C'),,drop=FALSE]
 stopifnot(all(c('tf','target','mor','confidence') %in% colnames(d)))
 write.table(d,file=file.path(commandArgs(TRUE)[1],paste0(n,'_ABC_v',packageVersion('dorothea'),'.tsv')),sep='\t',quote=FALSE,row.names=FALSE)
 cat(n,'rows',nrow(d),'TF',length(unique(d$tf)),'targets',length(unique(d$target)),'\n')
}
