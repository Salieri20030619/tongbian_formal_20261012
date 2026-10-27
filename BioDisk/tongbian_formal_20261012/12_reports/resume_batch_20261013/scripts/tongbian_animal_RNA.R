options(stringsAsFactors=FALSE,warn=1,digits=16)
suppressPackageStartupMessages({library(DESeq2);library(jsonlite)})
root<-'/mnt/BioDisk/wangmingyuan/tongbian_formal_20261012';input<-file.path(root,'03_inputs/GSE319612_verified_v1');out<-file.path(root,'07_animal_multiomics/GSE319612_DESeq2_v1')
stopifnot(!dir.exists(out));dir.create(out,recursive=TRUE);for(n in c('tables','audit','figures','models'))dir.create(file.path(out,n))
wr<-function(x,n)write.table(x,file.path(out,n),sep='\t',quote=FALSE,row.names=FALSE,na='NA')
meta<-read.delim(file.path(input,'sample_metadata.tsv'),check.names=FALSE);raw<-read.delim(gzfile(file.path(input,'GSE319612_gene_count.txt.gz')),check.names=FALSE,quote='',comment.char='')
stopifnot(nrow(meta)==26,!anyDuplicated(meta$sample_id),!anyDuplicated(raw$gene_id),all(meta$group_verified))
meta$group<-factor(meta$group,levels=c('NC','MC','RG'));stopifnot(identical(as.integer(table(meta$group)),c(8L,8L,10L)))
Y<-as.matrix(raw[,meta$sample_id]);rownames(Y)<-raw$gene_id;stopifnot(all(is.finite(Y)),all(Y>=0),all(Y==round(Y)),all(Y<.Machine$integer.max));storage.mode(Y)<-'integer'
rownames(meta)<-meta$sample_id;stopifnot(identical(colnames(Y),rownames(meta)));D<-model.matrix(~group,meta);stopifnot(qr(D)$rank==3)
keep<-rowSums(Y>=10)>=8
cfg<-list(frozen_before_fit_utc=format(Sys.time(),tz='UTC',usetz=TRUE),design='~group',n=c(NC=8,MC=8,RG=10),prefilter='count >=10 in >=8 samples (minimum group size; extension of human minimum-group rule)',contrasts=c('MC_vs_NC','RG_vs_MC'),method='DESeq2 Wald',outlier_replacement='disabled: minReplicatesForReplace=Inf; default Cook filtering retained',independentFiltering=FALSE,BH='separately per contrast over finite P among prefiltered genes',reverse='both nominal P<0.05 and opposite log2FC; no FC threshold',crossomics_pairing=FALSE,intervention='Rothia, not Tongbian formula',VST='blind TRUE for sample QC and later functional analysis',R=as.character(getRversion()),DESeq2=as.character(packageVersion('DESeq2')))
write_json(cfg,file.path(out,'analysis_config_FROZEN.json'),pretty=TRUE,auto_unbox=TRUE)
wr(meta,'audit/metadata_verified.tsv');wr(data.frame(sample_id=meta$sample_id,D,check.names=FALSE),'audit/design_matrix.tsv');wr(data.frame(gene_id=rownames(Y),prefilter_pass=keep),'audit/prefilter_all_genes.tsv')
sha<-Sys.glob(file.path(input,'*'));hashes<-sapply(sha,function(p)strsplit(system2('sha256sum',shQuote(p),stdout=TRUE),' ')[[1]][1]);wr(data.frame(file=sha,sha256=hashes),'audit/input_sha256.tsv')
cat('PREFLIGHT_PASS',nrow(Y),'genes',sum(keep),'eligible; 26 unique sample IDs; NC8 MC8 RG10; design rank3\n')
dds<-DESeqDataSetFromMatrix(Y[keep,,drop=FALSE],meta,~group);dds<-DESeq(dds,minReplicatesForReplace=Inf,parallel=FALSE)
tabs<-list();summ<-list()
for(nm in c('MC_vs_NC','RG_vs_MC')){
 cc<-if(nm=='MC_vs_NC')c('group','MC','NC') else c('group','RG','MC')
 res<-as.data.frame(results(dds,contrast=cc,alpha=.05,independentFiltering=FALSE));ix<-match(raw$gene_id,rownames(res));tab<-data.frame(gene_id=raw$gene_id,gene_name=raw$gene_name,prefilter_pass=keep,res[ix,,drop=FALSE]);rownames(tab)<-NULL
 tab$BH<-p.adjust(tab$pvalue,'BH');tab$testable<-is.finite(tab$pvalue);tab$contrast<-nm
 stopifnot(all(abs(tab$padj[tab$testable]-tab$BH[tab$testable])<1e-12));wr(tab,paste0('tables/',nm,'_ALL_genes.tsv'));tabs[[nm]]<-tab
 r<-tab[tab$testable&is.finite(tab$stat),];r<-r[order(-r$stat,r$gene_id),];wr(r,paste0('tables/',nm,'_all_testable_signed_Wald.tsv'))
 summ[[nm]]<-data.frame(contrast=nm,total_genes=nrow(Y),prefilter=sum(keep),tested=sum(tab$testable),nominal=sum(tab$pvalue<.05,na.rm=TRUE),BH=sum(tab$BH<.05,na.rm=TRUE),up=sum(tab$pvalue<.05&tab$log2FoldChange>0,na.rm=TRUE),down=sum(tab$pvalue<.05&tab$log2FoldChange<0,na.rm=TRUE))
}
a<-tabs[[1]];b<-tabs[[2]];stopifnot(identical(a$gene_id,b$gene_id));rev<-data.frame(gene_id=a$gene_id,gene_name=a$gene_name,MC_vs_NC_log2FC=a$log2FoldChange,MC_vs_NC_P=a$pvalue,MC_vs_NC_BH=a$BH,RG_vs_MC_log2FC=b$log2FoldChange,RG_vs_MC_P=b$pvalue,RG_vs_MC_BH=b$BH);rev$reverse_nominal<-with(rev,is.finite(MC_vs_NC_P)&is.finite(RG_vs_MC_P)&MC_vs_NC_P<.05&RG_vs_MC_P<.05&MC_vs_NC_log2FC*RG_vs_MC_log2FC<0)
wr(rev,'tables/two_contrast_ALL_genes.tsv');wr(rev[rev$reverse_nominal,],'tables/opposite_direction_nominal_candidates.tsv')
vs<-varianceStabilizingTransformation(dds,blind=TRUE);V<-assay(vs);wr(data.frame(gene_id=rownames(V),V,check.names=FALSE),'tables/blind_VST_prefiltered.tsv')
pc<-prcomp(t(V));pct<-100*pc$sdev^2/sum(pc$sdev^2);wr(data.frame(sample_id=rownames(pc$x),group=meta$group,pc$x,check.names=FALSE),'tables/PCA_scores.tsv');wr(data.frame(PC=seq_along(pct),percent_variance=pct),'tables/PCA_variance.tsv')
wr(data.frame(sample_id=meta$sample_id,group=meta$group,library_size=colSums(Y),detected_genes=colSums(Y>0),size_factor=sizeFactors(dds)),'tables/sample_QC.tsv')
cols<-c(NC='#32788a',MC='#c47a4b',RG='#628c56');png(file.path(out,'figures/Figure_animal_RNA_PCA.png'),width=1800,height=1300,res=180);par(mar=c(5,5,2,2));plot(pc$x[,1],pc$x[,2],pch=19,col=cols[meta$group],xlab=sprintf('PC1 (%.1f%%)',pct[1]),ylab=sprintf('PC2 (%.1f%%)',pct[2]));text(pc$x[,1],pc$x[,2],labels=meta$sample_id,pos=3,cex=.7);legend('topright',legend=names(cols),col=cols,pch=19);dev.off()
png(file.path(out,'figures/Figure_animal_reverse.png'),width=1500,height=1300,res=180);plot(rev$MC_vs_NC_log2FC,rev$RG_vs_MC_log2FC,col=ifelse(rev$reverse_nominal,'#a74935','#c7cdd0'),pch=16,cex=.45,xlab='MC minus NC: log2FC',ylab='RG minus MC: log2FC');abline(h=0,v=0,col='gray');dev.off()
saveRDS(dds,file.path(out,'models/DESeq2_fit.rds'));saveRDS(vs,file.path(out,'models/blind_VST.rds'));wr(do.call(rbind,summ),'tables/analysis_summary.tsv');capture.output(sessionInfo(),file=file.path(out,'audit/sessionInfo.txt'))
write_json(list(status='COMPUTE_COMPLETE',contrasts=summ,reverse_nominal=sum(rev$reverse_nominal),crossomics_pairing_verified=FALSE,completed_utc=format(Sys.time(),tz='UTC',usetz=TRUE)),file.path(out,'COMPUTE_COMPLETE.json'),pretty=TRUE,auto_unbox=TRUE)
print(as.data.frame(do.call(rbind,summ)));cat('REVERSE_NOMINAL',sum(rev$reverse_nominal),'\n')
