root<-'/mnt/BioDisk/wangmingyuan/tongbian_formal_20261012'
out<-file.path(root,'05_human_RNA/fgsea_v1');dir.create(out,recursive=TRUE,showWarnings=FALSE)
x<-readRDS(file.path(root,'05_human_RNA/host_functional_v1/functional_inputs_frozen.rds'))
ref<-file.path(root,'03_inputs/references/msigdbr_26.1.0_20261012')
gmt<-function(p){a<-strsplit(readLines(p,warn=FALSE),'\t',fixed=TRUE);setNames(lapply(a,function(z)unique(z[-c(1,2)])),vapply(a,'[',character(1),1))}
available<-list(C2_CP=file.path(ref,'C2_CP.human.symbols.gmt'),GO_BP=file.path(ref,'GO_BP.human.symbols.gmt'))
if(!file.exists(available$C2_CP))available$C2_CP<-file.path(root,'03_inputs/references/c2.cp.v2026.1.Hs.symbols.gmt')
results<-list();resource_status<-list()
for(i in seq_along(available)){
 family<-names(available)[i];path<-available[[i]]
 if(!file.exists(path)){resource_status[[family]]<-list(status='BLOCKED_RESOURCE_MISSING',path=path);next}
 sets<-gmt(path);stopifnot(!anyDuplicated(names(sets)))
 cov<-data.frame(pathway=names(sets),family=family,total_genes=lengths(sets),matched_ranked_genes=vapply(sets,function(g)sum(g %in% names(x$ranks)),integer(1)))
 cov$tested_by_size_rule<-cov$matched_ranked_genes>=15&cov$matched_ranked_genes<=500
 write.table(cov,file.path(out,paste0(family,'_ALL_set_coverage.tsv')),sep='\t',quote=FALSE,row.names=FALSE)
 set.seed(20261012+i)
 warnings_seen<-character()
 run_fgsea<-function(ss,perms){withCallingHandlers(fgsea::fgseaMultilevel(pathways=ss,stats=x$ranks,minSize=15,maxSize=500,eps=0,sampleSize=101,nPermSimple=perms,nproc=1,scoreType='std'),warning=function(w){warnings_seen<<-c(warnings_seen,conditionMessage(w))})}
 res<-run_fgsea(sets,1000)
 saveRDS(res,file.path(out,paste0(family,'_first_pass_unmodified.rds')))
 failed<-res$pathway[!is.finite(res$pval)]
 if(length(failed)){
  set.seed(20262012+i)
  retry<-run_fgsea(sets[failed],10000)
  saveRDS(retry,file.path(out,paste0(family,'_nonfinite_P_retry_unmodified.rds')))
  good_retry<-retry[is.finite(retry$pval),]
  if(nrow(good_retry)){for(k in seq_len(nrow(good_retry)))res[match(good_retry$pathway[k],res$pathway),]<-good_retry[k,]}
 }
 res<-as.data.frame(res);res$family<-family;res$leadingEdge<-vapply(res$leadingEdge,paste,collapse=';',character(1))
 res$padj<-p.adjust(res$pval,method='BH')
 res$retried_after_nonfinite_P<-res$pathway %in% failed
 res$numerical_zero_P<-is.finite(res$pval)&res$pval==0
 res$P_display<-ifelse(!is.finite(res$pval),'NA',ifelse(res$pval==0,'numerical underflow; not exact P=0',format(res$pval,scientific=TRUE,digits=8)))
 writeLines(unique(warnings_seen),file.path(out,paste0(family,'_all_warnings.txt')))
 names(res)[names(res)=='padj']<-'BH_within_collection'
 write.table(res,file.path(out,paste0(family,'_fgsea_ALL_results.tsv')),sep='\t',quote=FALSE,row.names=FALSE,na='NA')
 results[[family]]<-res
 resource_status[[family]]<-list(status='COMPLETED',resource=path,catalog_sets=length(sets),eligible_sets=sum(cov$tested_by_size_rule),returned_sets=nrow(res),finite_P=sum(is.finite(res$pval)),NA_P=sum(!is.finite(res$pval)),retried_sets=length(failed),numerical_zero_P=sum(res$numerical_zero_P),minimum_positive_P=min(res$pval[is.finite(res$pval)&res$pval>0]),warning_count=length(unique(warnings_seen)),nominal_P_lt_005=sum(res$pval<.05,na.rm=TRUE),BH_within_lt_005=sum(res$BH_within_collection<.05,na.rm=TRUE))
 jsonlite::write_json(resource_status,file.path(out,'collection_execution_status.json'),pretty=TRUE,auto_unbox=TRUE)
 cat(family,jsonlite::toJSON(resource_status[[family]],auto_unbox=TRUE),'\n')
}
stopifnot(length(results)>0)
allres<-do.call(rbind,results);allres$BH_across_collections<-p.adjust(allres$pval,method='BH')
write.table(allres,file.path(out,'fgsea_all_collections_ALL_results.tsv'),sep='\t',quote=FALSE,row.names=FALSE,na='NA')
jsonlite::write_json(list(status=if(all(vapply(resource_status,function(z)z$status=='COMPLETED',logical(1))))'COMPLETED' else 'PARTIAL_RESOURCE_BLOCKED',ranked_symbols=length(x$ranks),all_testable_genes_used=TRUE,significance_prescreen=FALSE,collection_status=resource_status,minSize=15,maxSize=500,eps=0,seed='20261012 + collection index',P_primary=TRUE,BH_auxiliary=TRUE,interpretation='directional pathway expression evidence in an observational disease comparison'),file.path(out,'fgsea_summary.json'),pretty=TRUE,auto_unbox=TRUE)
writeLines(capture.output(sessionInfo()),file.path(out,'sessionInfo.txt'))
png(file.path(out,'Figure_fgsea_top_pathways_display_only.png'),width=1700,height=1400,res=150,type='cairo')
par(mfrow=c(length(results),1),mar=c(5,26,3,2))
for(fam in names(results)){
 rr<-results[[fam]];rr<-rr[is.finite(rr$pval),];rr<-head(rr[order(rr$pval,rr$pathway),],12);rr<-rr[order(rr$NES),]
 barplot(rr$NES,horiz=TRUE,names.arg=gsub('_',' ',rr$pathway),las=1,cex.names=.55,col=ifelse(rr$NES>0,'#BA5B3B','#2F6F8F'),xlab='NES: STC relative to cancer-surgery controls',main=paste(fam,'top 12 by P for display only'))
 abline(v=0,col='grey40')
}
dev.off()
