options(stringsAsFactors=FALSE,digits=16)
P<-'/mnt/BioDisk/wangmingyuan/tongbian_formal_20261012';A<-file.path(P,'07_animal_multiomics/GSE319612_DESeq2_v1');O<-file.path(P,'07_animal_multiomics/GSE319612_functional_integrative_v1');checks<-list()
add<-function(check,pass,detail){stopifnot(pass);checks[[length(checks)+1]]<<-data.frame(check=check,passed=pass,detail=as.character(detail))}
for(nm in c('MC_vs_NC','RG_vs_MC')){
 a<-read.delim(file.path(A,'tables',paste0(nm,'_ALL_genes.tsv')));good<-is.finite(a$pvalue)
 add(paste(nm,'Wald_P'),max(abs(a$pvalue[good]-2*pnorm(-abs(a$stat[good]))))<1e-12,sum(good))
 add(paste(nm,'BH'),max(abs(a$BH[good]-p.adjust(a$pvalue[good],'BH')))<1e-12,sum(a$BH<.05,na.rm=TRUE))
}
r<-read.delim(file.path(A,'tables/two_contrast_ALL_genes.tsv'));expected<-with(r,is.finite(MC_vs_NC_P)&is.finite(RG_vs_MC_P)&MC_vs_NC_P<.05&RG_vs_MC_P<.05&MC_vs_NC_log2FC*RG_vs_MC_log2FC<0);add('reverse98',identical(expected,r$reverse_nominal)&&sum(expected)==98,sum(expected))
for(family in c('ssGSEA','VIPER'))for(nm in c('MC_vs_NC','RG_vs_MC')){
 t<-read.delim(file.path(O,'tables',paste0(family,'_',nm,'_ALL.tsv')));add(paste(family,nm,'BH'),max(abs(t$BH-p.adjust(t$P,'BH')))<1e-12,nrow(t));add(paste(family,nm,'finite'),all(is.finite(t$P))&&all(is.finite(t$delta_score)),sum(t$P<.05))
}
for(coll in c('C2_CP','GO_BP'))for(nm in c('MC_vs_NC','RG_vs_MC')){
 t<-read.delim(file.path(O,'tables',paste0(coll,'_',nm,'_fgsea_ALL.tsv')));g<-is.finite(t$pval);add(paste(coll,nm,'BH'),max(abs(t$BH[g]-p.adjust(t$pval[g],'BH')))<1e-12,sum(g));cv<-read.delim(file.path(O,'audit',paste0(coll,'_',nm,'_coverage_ALL.tsv')));add(paste(coll,nm,'pathway_set'),setequal(t$pathway,cv$pathway[cv$eligible]),nrow(t))
}
f<-read.delim(file.path(O,'tables/human_mouse_directional_Fisher.tsv'));ref<-with(f,phyper(overlap_n-1,human_set_n,background_n-human_set_n,mouse_set_n,lower.tail=FALSE));add('Fisher_independent_hypergeometric',max(abs(ref-f$P))<1e-12,nrow(f))
b<-read.delim(file.path(O,'tables/human_mouse_common_testable_background.tsv'));add('background_unique_testable',!anyDuplicated(b$human_id)&&!anyDuplicated(b$mouse_id)&&all(is.finite(b$human_P))&&all(is.finite(b$mouse_P)),nrow(b))
write.table(do.call(rbind,checks),file.path(O,'audit/result_integrity_checks.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
jsonlite::write_json(list(status='VERIFIED',checks=length(checks),all_passed=TRUE,time_utc=format(Sys.time(),tz='UTC',usetz=TRUE)),file.path(O,'audit/VALIDATION_COMPLETE.json'),pretty=TRUE,auto_unbox=TRUE)
print(as.data.frame(do.call(rbind,checks)))
for(nm in c('MC_vs_NC','RG_vs_MC'))print(as.data.frame(read.delim(file.path(O,'tables',paste0('ssGSEA_',nm,'_ALL.tsv')))))
for(family in c('ssGSEA','VIPER')){t<-read.delim(file.path(O,'tables',paste0(family,'_opposite_direction_ALL.tsv')));cat(family,'REVERSE',sum(t$reverse_nominal),'\n')}
