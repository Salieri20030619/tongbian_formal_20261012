# Preparation artifact. Default validates input only. --run is required to fit.
args <- commandArgs(trailingOnly=TRUE)
root <- '/mnt/BioDisk/wangmingyuan/tongbian_formal_20261012'
stopifnot(Sys.info()[['sysname']] == 'Linux', Sys.info()[['nodename']] == 'osv01')
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(jsonlite))
identity <- fromJSON(file.path(root,'tongbian_project_identity.json'))
stopifnot(identity$project_id=='tongbian_formal_20261012')
pilot <- file.path(root,'03_inputs/GSE245885_pilot')
meta <- read.delim(file.path(pilot,'metadata/samples.tsv'),check.names=FALSE,stringsAsFactors=FALSE)
stopifnot(nrow(meta)==10, sum(meta$condition=='STC')==5, sum(meta$condition=='Control')==5)
cc <- list();annotation <- NULL
for (i in seq_len(nrow(meta))) {
 s<-meta$sample[i]
 d<-read.delim(gzfile(file.path(pilot,'raw',meta$supplement_file[i])),check.names=FALSE,stringsAsFactors=FALSE)
 stopifnot(identical(colnames(d),c('ensembl_id','gene_name',paste0(s,'.count'),paste0(s,'.FPKM'))))
 stopifnot(!anyDuplicated(d$ensembl_id))
 cts<-d[[paste0(s,'.count')]]
 stopifnot(all(is.finite(cts)), all(cts>=0), all(cts==floor(cts)))
 a<-d[,c('ensembl_id','gene_name')]
 if(is.null(annotation)) annotation<-a else stopifnot(identical(annotation,a))
 cc[[s]]<-as.integer(cts)
}
counts<-do.call(cbind,cc);rownames(counts)<-annotation$ensembl_id
stopifnot(nrow(counts)==35671,identical(colnames(counts),meta$sample))
keep<-rowSums(counts>=10)>=5
stopifnot(sum(keep)==20479)
meta$condition<-factor(meta$condition,levels=c('Control','STC'));rownames(meta)<-meta$sample
stopifnot(qr(model.matrix(~condition,meta))$rank==2)
cat('INPUT_VALIDATION_PASS: 35671 x 10; eligible 20479; 5 STC and 5 cancer-surgery controls\n')
if(!('--run' %in% args)) {
 cat('CHECK_ONLY: no model fit, VST or formal output produced.\n')
 quit(status=0)
}
out<-file.path(root,'05_human_RNA','R_DESeq2_bridge_v1')
if(dir.exists(out)) stop('Refuse to overwrite existing analysis directory')
dir.create(out,recursive=TRUE)
set.seed(20261012)
writeLines(capture.output(sessionInfo()),file.path(out,'sessionInfo.txt'))
write_json(list(status='RUNNING',started_utc=format(Sys.time(),tz='UTC'),method='R DESeq2',version=as.character(packageVersion('DESeq2')),contrast=c('condition','STC','Control'),design='~condition',primary='nominal P<0.05 exploratory',BH='retained',FC_threshold=NULL,old_method='PyDESeq2 0.5.4',old_scores='mean gene-z, not ssGSEA'),file.path(out,'run_status.json'),pretty=TRUE,auto_unbox=TRUE)
dds<-DESeqDataSetFromMatrix(countData=counts[keep,,drop=FALSE],colData=meta,design=~condition)
dds<-DESeq(dds,fitType='parametric',sfType='ratio',betaPrior=FALSE,minReplicatesForReplace=7,parallel=FALSE)
rr<-as.data.frame(results(dds,contrast=c('condition','STC','Control'),cooksCutoff=TRUE,independentFiltering=FALSE,alpha=0.05))
rr$ensembl_id<-rownames(rr)
rr$log2FC_CI95_low<-rr$log2FoldChange-qnorm(.975)*rr$lfcSE
rr$log2FC_CI95_high<-rr$log2FoldChange+qnorm(.975)*rr$lfcSE
rr$test_status<-ifelse(is.na(rr$pvalue),'cooks_filtered_or_unestimable','tested')
full<-merge(annotation,rr,by='ensembl_id',all.x=TRUE,sort=FALSE)
full<-full[match(annotation$ensembl_id,full$ensembl_id),]
full$passed_prefilter<-keep
full$test_status[!keep]<-'low_count_prefilter'
write.table(full,file.path(out,'R_DESeq2_all_genes.tsv'),sep='\t',quote=FALSE,row.names=FALSE,na='NA')
write.table(data.frame(sample=colnames(dds),size_factor=sizeFactors(dds)),file.path(out,'size_factors.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
saveRDS(dds,file.path(out,'R_DESeq2_fit.rds'))
vst<-varianceStabilizingTransformation(dds,blind=TRUE)
write.table(assay(vst),file.path(out,'R_DESeq2_blind_VST_tested.tsv'),sep='\t',quote=FALSE,col.names=NA)
old<-read.delim(gzfile(file.path(pilot,'results/differential_all_genes.tsv.gz')),check.names=FALSE)
bridge<-merge(full[,c('ensembl_id','log2FoldChange','stat','pvalue','padj','test_status')],old[,c('ensembl_id','log2FoldChange','stat','pvalue','BH_FDR_all_tested','test_status')],by='ensembl_id',all=TRUE,suffixes=c('_R','_PyDESeq2'))
bridge$log2FC_difference_R_minus_Python<-bridge$log2FoldChange_R-bridge$log2FoldChange_PyDESeq2
write.table(bridge,file.path(out,'R_vs_PyDESeq2_method_bridge.tsv'),sep='\t',quote=FALSE,row.names=FALSE,na='NA')
write_json(list(status='COMPLETED',completed_utc=format(Sys.time(),tz='UTC'),method='R DESeq2',version=as.character(packageVersion('DESeq2')),eligible_genes=sum(keep),finite_primary_P=sum(is.finite(rr$pvalue)),no_FC_threshold=TRUE,VST_is_not_ssGSEA=TRUE,modelA_rerun=FALSE),file.path(out,'run_status.json'),pretty=TRUE,auto_unbox=TRUE)
