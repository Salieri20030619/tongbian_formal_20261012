args <- commandArgs(trailingOnly=TRUE)
o <- args[1]; setwd(o)
suppressPackageStartupMessages(library(Maaslin2))
writeLines(capture.output(sessionInfo()),'R_sessionInfo.txt')
meta <- read.delim('metadata.tsv',row.names=1,check.names=FALSE,stringsAsFactors=FALSE)
stopifnot(nrow(meta)==60, !anyDuplicated(rownames(meta)),all(table(meta$group)==20))
write_tab <- function(x,n) write.table(x,n,sep='\t',row.names=FALSE,quote=FALSE,na='NA',fileEncoding='UTF-8')
classes <- c(species='species_percent.tsv',pathway='pathabundance_cpm.tsv')
inputs <- list(); maps <- list()
for (cl in names(classes)) {
 z <- read.delim(classes[[cl]],row.names=1,check.names=FALSE)
 stopifnot(identical(colnames(z),rownames(meta)),all(is.finite(as.matrix(z))),all(as.matrix(z)>=0))
 keep <- rowSums(z>0)>=6 & apply(z,1,var)>0
 map <- data.frame(feature_id=sprintf(paste0(cl,'_%05d'),seq_len(nrow(z))),feature=rownames(z),positive_n=rowSums(z>0),retained=keep)
 write_tab(map,paste0(cl,'_feature_filter.tsv')); rownames(z)<-map$feature_id
 inputs[[cl]] <- t(as.matrix(z[keep,,drop=FALSE])); maps[[cl]]<-map
}
run_fit <- function(cl,dat,ref,tag) {
 md <- data.frame(group=factor(meta$group,levels=c(ref,setdiff(c('H','DCTT','NCTT'),ref))),row.names=rownames(meta))
 design<-model.matrix(~group,md); stopifnot(qr(design)$rank==3)
 fit<-Maaslin2(input_data=as.data.frame(dat),input_metadata=md,output=tag,fixed_effects='group',reference=paste0('group,',ref),min_abundance=0,min_prevalence=0,min_variance=0,normalization='NONE',transform='LOG',analysis_method='LM',standardize=FALSE,correction='BH',cores=1,plot_heatmap=FALSE,plot_scatter=FALSE,max_significance=0.05,save_models=TRUE)
 ans<-read.delim(file.path(tag,'all_results.tsv'),check.names=FALSE)
 stopifnot(nrow(ans)==2*ncol(dat),all(ans$metadata=='group'),all(is.finite(ans$pval)),all(ans$pval>=0 & ans$pval<=1))
 maxerr<-0
 for (f in colnames(dat)) {
  y<-dat[,f]; y[y==0]<-min(y[y>0])/2; y<-log2(y)
  sm<-coef(summary(lm(y~group,data=md)))
  ix<-which(ans$feature==f)
  for (i in ix) {
   expected<-sm[paste0('group',ans$value[i]),c('Estimate','Std. Error','Pr(>|t|)')]
   actual<-as.numeric(ans[i,c('coef','stderr','pval')]); err<-max(abs(expected-actual)); maxerr<-max(maxerr,err)
   stopifnot(err<1e-8)
  }
 }
 writeLines(paste('independent_lm_max_absolute_error',format(maxerr,digits=16)),file.path(tag,'independent_validation.txt'))
 ans
}
for (cl in names(inputs)) {
 dat<-inputs[[cl]][,seq_len(min(10,ncol(inputs[[cl]]))),drop=FALSE]
 run_fit(cl,dat,'H',paste0('pilot_',cl))
}
writeLines('PASS: two feature classes; all pilot coefficient/SE/P comparisons < 1e-8','PILOT_PASS.txt')
allresults<-list()
for (cl in names(inputs)) {
 dat<-inputs[[cl]]; map<-maps[[cl]]
 for (ref in c('H','NCTT')) {
  ans<-run_fit(cl,dat,ref,paste0('full_',cl,'_reference_',ref))
  if(ref=='NCTT') ans<-ans[ans$value=='DCTT',,drop=FALSE]
  ans$class<-cl; ans$contrast<-paste(ans$value,ref,sep='-'); ans$feature_id<-ans$feature; ans$feature<-map$feature[match(ans$feature_id,map$feature_id)]
  ans$package_qval<-ans$qval; ans$qval<-NULL
  ans$BH<-ave(ans$pval,ans$contrast,FUN=function(x)p.adjust(x,'BH'))
  ans$nominal_P_lt_0_05<-ans$pval<0.05; ans$BH_lt_0_05<-ans$BH<0.05
  ans$n_total<-60; ans$n_target<-20; ans$n_reference<-20; ans$residual_df<-57
  ans$CI95_low<-ans$coef-qt(.975,57)*ans$stderr; ans$CI95_high<-ans$coef+qt(.975,57)*ans$stderr
  ans$test_family<-paste(cl,ans$contrast,sep=':')
  allresults[[length(allresults)+1]]<-ans
 }
}
res<-do.call(rbind,allresults); write_tab(res,'association_all.tsv'); write_tab(res[res$pval<.05,],'association_nominal_P_lt_0_05.tsv'); write_tab(res[res$BH<.05,],'association_BH_lt_0_05.tsv')
summary<-do.call(rbind,lapply(split(res,res$test_family),function(z)data.frame(family=z$test_family[1],tested=nrow(z),nominal=sum(z$pval<.05),BH_significant=sum(z$BH<.05),positive_nominal=sum(z$pval<.05 & z$coef>0),negative_nominal=sum(z$pval<.05 & z$coef<0))))
write_tab(summary,'association_summary.tsv')
writeLines('PASS: all coefficients, standard errors and P values reproduced with independent lm; finite P and complete contrasts','STATISTICS_VALIDATED.txt')
print(as.data.frame(summary)); cat('STATISTICS_COMPLETE\n')
