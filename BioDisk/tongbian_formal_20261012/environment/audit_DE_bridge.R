root <- '/mnt/BioDisk/wangmingyuan/tongbian_formal_20261012'
out<-file.path(root,'05_human_RNA/R_DESeq2_bridge_v1')
stopifnot(jsonlite::fromJSON(file.path(out,'run_status.json'))$status=='COMPLETED')
d<-read.delim(file.path(out,'R_DESeq2_all_genes.tsv'),check.names=FALSE)
b<-read.delim(file.path(out,'R_vs_PyDESeq2_method_bridge.tsv'),check.names=FALSE)
valid<-is.finite(b$pvalue_R)&is.finite(b$pvalue_PyDESeq2)&is.finite(b$log2FoldChange_R)&is.finite(b$log2FoldChange_PyDESeq2)
common<-b[valid,]
nomR<-!is.na(b$pvalue_R)&b$pvalue_R<.05;nomPy<-!is.na(b$pvalue_PyDESeq2)&b$pvalue_PyDESeq2<.05
ds<-data.frame(ensembl_id=b$ensembl_id,R_status=b$test_status_R,Python_status=b$test_status_PyDESeq2,R_nominal=nomR,Python_nominal=nomPy,R_log2FC=b$log2FoldChange_R,Python_log2FC=b$log2FoldChange_PyDESeq2)
write.table(ds[is.na(b$pvalue_R)!=is.na(b$pvalue_PyDESeq2)|nomR!=nomPy,,drop=FALSE],file.path(out,'bridge_testability_and_threshold_changes.tsv'),sep='\t',quote=FALSE,row.names=FALSE,na='NA')
summary<-list(dataset='GSE245885',contrast='STC minus cancer-surgery normal-colon controls',n=c(STC=5,Control=5),method='R DESeq2 1.50.2',all_genes=nrow(d),eligible_genes=sum(d$passed_prefilter),finite_primary_P=sum(is.finite(d$pvalue)),nominal_P_lt_005=sum(d$pvalue<.05,na.rm=TRUE),nominal_up=sum(d$pvalue<.05&d$log2FoldChange>0,na.rm=TRUE),nominal_down=sum(d$pvalue<.05&d$log2FoldChange<0,na.rm=TRUE),BH_lt_005=sum(d$padj<.05,na.rm=TRUE),zero_P_values=sum(d$pvalue==0,na.rm=TRUE),common_finite_tests=nrow(common),log2FC_pearson=cor(common$log2FoldChange_R,common$log2FoldChange_PyDESeq2),log2FC_median_absolute_difference=median(abs(common$log2FC_difference_R_minus_Python)),log2FC_absolute_difference_quantiles=quantile(abs(common$log2FC_difference_R_minus_Python),c(.5,.9,.95,.99,1)),sign_agreement_fraction=mean(sign(common$log2FoldChange_R)==sign(common$log2FoldChange_PyDESeq2)),nominal_overlap=sum(nomR&nomPy),nominal_R_only=sum(nomR&!nomPy),nominal_Python_only=sum(!nomR&nomPy),R_testable_Python_not=sum(is.finite(b$pvalue_R)&!is.finite(b$pvalue_PyDESeq2)),Python_testable_R_not=sum(!is.finite(b$pvalue_R)&is.finite(b$pvalue_PyDESeq2)),no_forced_agreement=TRUE,no_FC_threshold=TRUE,interpretation='observational disease-associated tissue evidence; not formula efficacy')
jsonlite::write_json(summary,file.path(out,'bridge_audit_summary.json'),pretty=TRUE,auto_unbox=TRUE,digits=16)
cat(jsonlite::toJSON(summary,pretty=TRUE,auto_unbox=TRUE,digits=6),'\n')
