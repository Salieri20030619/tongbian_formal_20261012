options(stringsAsFactors=FALSE,warn=1,digits=16)
root <- '/mnt/BioDisk/wangmingyuan/tongbian_formal_20261012'
.libPaths(c(file.path(root,'environment/R_library_ppcor_20261013'),.libPaths()))
library(ppcor); library(jsonlite)
out <- file.path(root,'06_human_selected_metabolites/partial_spearman_v1')
stopifnot(!file.exists(file.path(out,'COMPUTE_COMPLETE.json')))
cfg <- fromJSON(file.path(out,'analysis_config_FROZEN.json')); indir<-cfg$input
read <- function(n)read.csv(file.path(indir,n),check.names=FALSE)
write <- function(x,n)write.table(x,file.path(out,n),sep='\t',quote=FALSE,row.names=FALSE,na='NA')
meta<-read('metadata_deidentified.csv');meta<-meta[order(meta$sample_id),]
stopifnot(nrow(meta)==60,!anyDuplicated(meta$sample_id))
align <- function(x){stopifnot(!anyDuplicated(x$sample_id),setequal(x$sample_id,meta$sample_id));x<-x[match(meta$sample_id,x$sample_id),];stopifnot(identical(x$sample_id,meta$sample_id));x}
raw<-align(read('metabolites_60x21_peak_area.csv')); joint<-align(read('joint_selected_60x17.csv')); micro<-align(read('microbial_selected_taxa_60x14.csv'));ann<-read('metabolites_21_annotation.csv'); conflicts<-read('cross_sheet_conflicts.csv')
stopifnot(nrow(conflicts)==1,conflicts$feature=='Cholic acid',conflicts$sample_id=='C6')
meta$group<-factor(meta$group,levels=c('Control','STC'));meta$SexCode<-factor(meta$Sex,levels=c(1,2),labels=c('code1','code2'))
stopifnot(all(table(meta$group)==30),!anyNA(meta[,c('group','Age','SexCode','BMI')]))
taxmap<-sapply(cfg$taxa,function(t){m<-names(micro)[endsWith(names(micro),paste0('|',t))];stopifnot(length(m)==1);m})
cross<-list()
for(t in cfg$taxa){a<-micro[[taxmap[t]]];b<-joint[[t]];stopifnot(max(abs(a-b))<1e-12);cross[[length(cross)+1]]<-data.frame(feature=t,max_relative_error=max(abs(a-b)/pmax(abs(a),1e-12)),rank_identical=identical(rank(a),rank(b)))}
for(id in cfg$primary40){nm<-ann$Name[match(id,ann$ID)];if(nm %in% names(cfg$joint_name_aliases))nm<-cfg$joint_name_aliases[[nm]];a<-raw[[id]];b<-joint[[nm]];err<-max(abs(a-b)/pmax(abs(a),1e-12));stopifnot(err<1e-7,identical(rank(a),rank(b)));cross[[length(cross)+1]]<-data.frame(feature=id,max_relative_error=err,rank_identical=TRUE)}
write(do.call(rbind,cross),'audit/cross_table_numeric_agreement.tsv')
stopifnot(abs(raw$M407T218[meta$sample_id=='C6']-22146617.1867411)<1e-6,abs(joint$`Cholic acid`[meta$sample_id=='C6']-146617.187)<1e-6)
write(data.frame(taxon=cfg$taxa,full_taxonomy=unname(taxmap),rank=ifelse(startsWith(cfg$taxa,'f__'),'family','genus')),'audit/taxon_mapping.tsv')
write(meta[,c('sample_id','group','Age','SexCode','BMI')],'audit/aligned_metadata.tsv')
write(raw,'audit/aligned_original_metabolites.tsv');write(micro[,c('sample_id',taxmap),drop=FALSE],'audit/aligned_original_taxa.tsv')
write(conflicts,'audit/source_conflicts_unchanged.tsv')
tests<-list();checks<-list();designs<-list()
calc<-function(family,features,model,scenario='source_A'){
 res<-list()
 for(taxon in cfg$taxa)for(id in features){
  m<-meta;x<-micro[[taxmap[taxon]]];y<-raw[[id]]
  if(scenario=='source_B')y[m$sample_id=='C6']<-146617.187
  covs<-switch(model,unadjusted=character(),group='group',group_age_sex_BMI=c('group','Age','SexCode','BMI'))
  keep<-is.finite(x)&is.finite(y)&complete.cases(m[,covs,drop=FALSE]);if(scenario=='exclude_C6')keep<-keep&m$sample_id!='C6'
  x<-x[keep];y<-y[keep];m<-m[keep,];n<-length(x)
  form<-if(length(covs))reformulate(covs) else ~1
  d<-model.matrix(form,m);z<-d[,-1,drop=FALSE];k<-ncol(z);df<-n-k-2
  rz<-if(k)apply(z,2,rank) else matrix(numeric(0),n,0);dr<-cbind(1,rz)
  stopifnot(qr(d)$rank==ncol(d),qr(dr)$rank==ncol(dr),n>k+2,length(unique(x))>1,length(unique(y))>1)
  stopifnot(qr(cbind(dr,rank(x),rank(y)))$rank==ncol(dr)+2)
  if(k){p<-ppcor::pcor.test(x,y,z,method='spearman');rho<-p$estimate;pv<-p$p.value;stat<-p$statistic;stopifnot(p$n==n,p$gp==k)}else{p<-cor.test(x,y,method='spearman',exact=FALSE);rho<-unname(p$estimate);pv<-p$p.value;stat<-rho*sqrt(df/(1-rho^2))}
  xr<-qr.resid(qr(dr),rank(x));yr<-qr.resid(qr(dr),rank(y));ref_r<-cor(xr,yr);ref_p<-2*pt(-abs(ref_r*sqrt(df/(1-ref_r^2))),df)
  stopifnot(abs(ref_r-rho)<1e-10,abs(ref_p-pv)<1e-10)
  ord<-rev(seq_len(n));if(k){perm<-ppcor::pcor.test(x[ord],y[ord],z[ord,,drop=FALSE],method='spearman');stopifnot(abs(perm$estimate-rho)<1e-10,abs(perm$p.value-pv)<1e-10)}
  key<-paste(family,model,scenario,sep='__')
  res[[length(res)+1]]<-data.frame(family=family,model=model,scenario=scenario,taxon=taxon,feature_id=id,metabolite=ann$Name[match(id,ann$ID)],n=n,n_Control=sum(m$group=='Control'),n_STC=sum(m$group=='STC'),covariate_count=k,covariates=paste(colnames(z),collapse=';'),design_rank=qr(d)$rank,rank_design_rank=qr(dr)$rank,df=df,rho=rho,statistic=stat,P=pv,status='OK',zero_taxon_n=sum(x==0),x_unique=length(unique(x)),y_unique=length(unique(y)),sample_ids=paste(m$sample_id,collapse=';'))
  checks[[length(checks)+1]]<<-data.frame(family=family,model=model,scenario=scenario,taxon=taxon,feature_id=id,rho_difference=abs(ref_r-rho),P_difference=abs(ref_p-pv),df=df,order_invariance_pass=TRUE)
  if(is.null(designs[[key]])){designs[[key]]<<-TRUE;write(data.frame(sample_id=m$sample_id,d,check.names=FALSE),paste0('audit/design_',key,'.tsv'))}
 }
 r<-do.call(rbind,res);r$BH<-p.adjust(r$P,'BH',n=length(cfg$taxa)*length(features));r$BH_family_n<-nrow(r);r$family_id<-paste(family,model,scenario,sep='__');r$direction<-ifelse(r$rho>0,'positive','negative');r$nominal_P_lt_005<-r$P<.05;r$BH_lt_005<-r$BH<.05
 write(r,paste0('tables/',family,'__',model,'__',scenario,'.tsv'));tests[[length(tests)+1]]<<-r
}
for(family in c('primary40','extended100'))for(model in cfg$models)calc(family,cfg[[family]],model)
for(scenario in cfg$cholic_scenarios)for(model in cfg$models)calc('cholic5',cfg$cholic5,model,scenario)
all<-do.call(rbind,tests);stopifnot(nrow(all)==465,all(is.finite(all$P)),all(abs(all$rho)<=1))
write(all,'tables/ALL_465_results.tsv');write(do.call(rbind,checks),'audit/independent_formula_and_order_checks.tsv')
summary<-do.call(rbind,lapply(split(all,all$family_id),function(t)data.frame(family_id=t$family_id[1],tests=nrow(t),n=min(t$n),k=t$covariate_count[1],nominal=sum(t$P<.05),BH=sum(t$BH<.05),positive_nominal=sum(t$P<.05&t$rho>0),negative_nominal=sum(t$P<.05&t$rho<0))))
write(summary,'tables/analysis_summary.tsv')
capture.output(sessionInfo(),file=file.path(out,'audit/analysis_sessionInfo.txt'))
write(data.frame(package=c('R','ppcor','MASS','jsonlite'),version=c(as.character(getRversion()),sapply(c('ppcor','MASS','jsonlite'),function(p)as.character(packageVersion(p))))),'audit/package_versions.tsv')
write_json(list(status='COMPUTE_COMPLETE_VERIFIED',rows=nrow(all),sample_n=60,primary_family=40,extended_family=100,cholic_family=5,original_conflict_retained=TRUE,primary_model='group',independent_formula_pass=TRUE,completed_utc=format(Sys.time(),tz='UTC',usetz=TRUE)),file.path(out,'COMPUTE_COMPLETE.json'),pretty=TRUE,auto_unbox=TRUE)
print(as.data.frame(summary));cat('PARTIAL_SPEARMAN_COMPLETE_AND_VERIFIED\n')
