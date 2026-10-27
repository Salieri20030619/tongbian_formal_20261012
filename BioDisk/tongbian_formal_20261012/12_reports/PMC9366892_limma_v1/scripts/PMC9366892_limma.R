options(stringsAsFactors=FALSE, warn=1, digits=16)
args <- commandArgs(trailingOnly=TRUE)
root <- '/mnt/BioDisk/wangmingyuan/tongbian_formal_20261012'
indir <- file.path(root,'03_inputs/PMC9366892_deidentified_v1/tongbian_human_PMC9366892_deidentified_inputs_20261012')
out <- file.path(root,'06_human_selected_metabolites/limma_v1')
run <- '--run' %in% args
suppressPackageStartupMessages(library(limma))
stopifnot(as.character(packageVersion('limma'))=='3.66.0')
read <- function(n) read.csv(file.path(indir,n),check.names=FALSE,stringsAsFactors=FALSE)
meta <- read('metadata_deidentified.csv'); meta <- meta[order(factor(meta$group,levels=c('Control','STC')),as.integer(sub('^[CT]','',meta$sample_id))),]
stopifnot(nrow(meta)==60L,!anyDuplicated(meta$sample_id),setequal(meta$sample_id,c(paste0('C',1:30),paste0('T',1:30))))
meta$group <- factor(meta$group,levels=c('Control','STC'))
meta$SexCode <- factor(meta$Sex,levels=c(1,2),labels=c('code1','code2'))
stopifnot(!anyNA(meta[,c('sample_id','group','Age','SexCode','BMI')]))
raw <- read('metabolites_60x21_peak_area.csv'); ann <- read('metabolites_21_annotation.csv'); joint <- read('joint_selected_60x17.csv')
stopifnot(!anyDuplicated(raw$sample_id),!anyDuplicated(joint$sample_id),setequal(raw$sample_id,meta$sample_id),setequal(joint$sample_id,meta$sample_id),!anyDuplicated(ann$ID))
raw <- raw[match(meta$sample_id,raw$sample_id),];joint <- joint[match(meta$sample_id,joint$sample_id),]
stopifnot(identical(raw$sample_id,meta$sample_id),identical(joint$sample_id,meta$sample_id),setequal(colnames(raw)[-1],ann$ID))
Yraw <- t(as.matrix(raw[,ann$ID,drop=FALSE]));storage.mode(Yraw)<-'double';colnames(Yraw)<-meta$sample_id
stopifnot(all(is.finite(Yraw)),all(Yraw>0),nrow(Yraw)==21,ncol(Yraw)==60)
cholic <- 'M407T218';stopifnot(abs(Yraw[cholic,'C6']-22146617.1867411)<1e-6,abs(joint$`Cholic acid`[meta$sample_id=='C6']-146617.187)<1e-6)
design_audit <- list()
for(drop_c6 in c(FALSE,TRUE)) for(adjusted in c(FALSE,TRUE)) {
 m <- meta[!drop_c6 | meta$sample_id!='C6',]
 d <- if(adjusted) model.matrix(~group+Age+SexCode+BMI,m) else model.matrix(~group,m)
 stopifnot(qr(d)$rank==ncol(d),nrow(d)>ncol(d),sum(m$group=='STC')==30)
 design_audit[[length(design_audit)+1]] <- data.frame(exclude_C6=drop_c6,adjusted=adjusted,n=nrow(d),n_Control=sum(m$group=='Control'),n_STC=sum(m$group=='STC'),rank=qr(d)$rank,columns=ncol(d),residual_df=nrow(d)-qr(d)$rank)
}
cat('INPUT_CHECK_PASS: 60 unique IDs, 21 positive raw peak-area features; all four designs full rank\n')
print(do.call(rbind,design_audit));cat('PPCOR_AVAILABLE',requireNamespace('ppcor',quietly=TRUE),'PARTIAL_CORRELATION_RUN FALSE\n')
if(!run) quit(status=0)
stopifnot(!dir.exists(out));dir.create(out,recursive=TRUE)
for(d in c('tables','audit','models')) dir.create(file.path(out,d))
write <- function(x,n) write.table(x,file.path(out,n),sep='\t',quote=FALSE,row.names=FALSE,na='NA')
write(do.call(rbind,design_audit),'audit/design_rank_checks.tsv')
write(meta[,c('sample_id','group','Age','Sex','SexCode','BMI')],'audit/analysis_metadata_aligned.tsv')
write(data.frame(feature_id=rownames(Yraw),Yraw,check.names=FALSE),'audit/original_peak_area_aligned.tsv')
Y <- log2(Yraw)
write(data.frame(feature_id=rownames(Y),Y,check.names=FALSE),'audit/log2_peak_area_aligned.tsv')
write(data.frame(sample_id='C6',feature_id=cholic,source_A=Yraw[cholic,'C6'],source_B=joint$`Cholic acid`[meta$sample_id=='C6'],input_modified=FALSE,source_B_change_scope='single cell in derived in-memory copy only'),'audit/C6_source_values.tsv')
fits <- list();results <- list();summaries <- list()
run_fit <- function(id,features=rownames(Y),scenario='source_A',adjusted=FALSE) {
 yy<-Y[features,,drop=FALSE];m<-meta
 if(scenario=='source_B') yy[cholic,'C6']<-log2(joint$`Cholic acid`[meta$sample_id=='C6'])
 if(scenario=='exclude_C6') {m<-meta[meta$sample_id!='C6',];yy<-yy[,m$sample_id,drop=FALSE]}
 d<-if(adjusted) model.matrix(~group+Age+SexCode+BMI,m) else model.matrix(~group,m)
 rownames(d)<-m$sample_id;stopifnot(identical(colnames(yy),rownames(d)),qr(d)$rank==ncol(d))
 fit<-eBayes(lmFit(yy,d),trend=FALSE,robust=FALSE)
 co<-match('groupSTC',colnames(d));se<-fit$stdev.unscaled[,co]*sqrt(fit$s2.post)
 df<-rep_len(fit$df.total,nrow(yy));est<-fit$coefficients[,co];p<-fit$p.value[,co];ci<-qt(.975,df)*se
 stopifnot(all(is.finite(p)),all(is.finite(est)),all(p>0),all(p<=1))
 if(!adjusted) stopifnot(max(abs(est-(rowMeans(yy[,m$group=='STC',drop=FALSE])-rowMeans(yy[,m$group=='Control',drop=FALSE]))))<1e-10)
 ix<-match(rownames(yy),ann$ID)
 tab<-data.frame(analysis_id=id,feature_id=rownames(yy),Name=ann$Name[ix],C6_source_conflict=rownames(yy)==cholic,main_20_feature=rownames(yy)!=cholic,scenario=scenario,design=if(adjusted)'group+Age+SexCode+BMI' else 'group',n_Control=sum(m$group=='Control'),n_STC=sum(m$group=='STC'),logFC=est,geometric_mean_ratio=2^est,CI95_low=est-ci,CI95_high=est+ci,moderated_SE=se,moderated_t=fit$t[,co],df_total=df,df_residual=fit$df.residual,P=p,BH=p.adjust(p,'BH'),BH_family_n=length(p),direction=ifelse(est>0,'STC_higher',ifelse(est<0,'STC_lower','equal')),mean_log2_Control=rowMeans(yy[,m$group=='Control',drop=FALSE]),mean_log2_STC=rowMeans(yy[,m$group=='STC',drop=FALSE]),AveExpr=fit$Amean)
 write(tab,paste0('tables/',id,'.tsv'))
 write(data.frame(sample_id=rownames(d),d,check.names=FALSE),paste0('audit/design_',id,'.tsv'))
 fits[[id]]<<-list(fit=fit,design=d,features=features,sample_ids=m$sample_id,scenario=scenario,adjusted=adjusted)
 results[[id]]<<-tab
 summaries[[id]]<<-data.frame(analysis_id=id,feature_n=nrow(tab),n_Control=sum(m$group=='Control'),n_STC=sum(m$group=='STC'),nominal_P_lt_005=sum(p<.05),BH_lt_005=sum(tab$BH<.05),P_up=sum(p<.05&est>0),P_down=sum(p<.05&est<0),design_rank=qr(d)$rank,df_residual=nrow(d)-qr(d)$rank)
 invisible(tab)
}
primary_features<-setdiff(rownames(Y),cholic)
run_fit('primary20_group',features=primary_features)
run_fit('primary20_age_sex_BMI',features=primary_features,adjusted=TRUE)
for(scenario in c('source_A','source_B','exclude_C6')) for(adjusted in c(FALSE,TRUE)) run_fit(paste0('all21_',scenario,if(adjusted)'_age_sex_BMI' else '_group'),scenario=scenario,adjusted=adjusted)
allresults<-do.call(rbind,results);write(allresults,'tables/all_fits_results.tsv')
cs<-allresults[allresults$feature_id==cholic,];write(cs,'tables/Cholic_acid_prespecified_sensitivity.tsv')
write(do.call(rbind,summaries),'tables/analysis_summary.tsv')
saveRDS(fits,file.path(out,'models/all_limma_fits.rds'))
write(ann,'audit/original_21_annotation.tsv')
capture.output(sessionInfo(),file=file.path(out,'audit/R_sessionInfo.txt'))
write(data.frame(package=c('R','limma','ppcor'),version=c(as.character(getRversion()),as.character(packageVersion('limma')),if(requireNamespace('ppcor',quietly=TRUE))as.character(packageVersion('ppcor')) else 'NOT_INSTALLED')),'audit/package_versions.tsv')
capture.output(print(do.call(rbind,summaries)),print(cs[,c('analysis_id','n_Control','n_STC','logFC','P','BH')]),file=file.path(out,'audit/result_summary.txt'))
cat('ALL_EIGHT_FITS_COMPLETED\n');print(do.call(rbind,summaries));print(cs[,c('analysis_id','n_Control','n_STC','logFC','P','BH')])
