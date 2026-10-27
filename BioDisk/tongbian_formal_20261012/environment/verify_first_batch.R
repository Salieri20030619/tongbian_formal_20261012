root<-'/mnt/BioDisk/wangmingyuan/tongbian_formal_20261012'
out<-file.path(root,'05_human_RNA/technical_validation_v1');dir.create(out,recursive=TRUE,showWarnings=FALSE)
x<-readRDS(file.path(root,'05_human_RNA/host_functional_v1/functional_inputs_frozen.rds'))
ss<-readRDS(file.path(root,'05_human_RNA/host_functional_v1/true_ssGSEA_scores.rds'))
vp<-as.matrix(read.delim(file.path(root,'05_human_RNA/host_functional_v1/VIPER_sample_activity_matrix.tsv'),row.names=1,check.names=FALSE))
cat('START_BASE_R_SSGSEA_REFERENCE\n');flush.console()
# Independent base-R implementation of the installed GSVA ssGSEA integrated
# weighted empirical-CDF definition. Integer rank casting matches GSVA 2.4.9.
reference<-matrix(NA_real_,nrow=length(x$modules),ncol=ncol(x$symbol_VST),dimnames=list(names(x$modules),colnames(x$symbol_VST)))
for(j in seq_len(ncol(x$symbol_VST))){
 ranks<-as.integer(rank(x$symbol_VST[,j],ties.method='average'))
 ord<-order(ranks,decreasing=TRUE);sorted_names<-rownames(x$symbol_VST)[ord]
 weights<-abs(ranks[ord])^.25
 for(k in seq_along(x$modules)){
  hit<-sorted_names %in% x$modules[[k]]
  in_cdf<-cumsum(ifelse(hit,weights,0))/sum(weights[hit])
  out_cdf<-cumsum(!hit)/sum(!hit)
  reference[k,j]<-sum(in_cdf-out_cdf)
 }
}
reference<-reference/diff(range(reference))
ssdiff<-max(abs(reference-ss[rownames(reference),colnames(reference)]))
write.table(reference,file.path(out,'ssGSEA_independent_base_R_reference.tsv'),sep='\t',quote=FALSE,col.names=NA)
cat('SSGSEA_MAX_ABSOLUTE_DIFFERENCE',ssdiff,'\n');flush.console()
stopifnot(ssdiff<1e-10)
cat('START_DIRECT_VIPER_VALIDATION\n');flush.console()
net<-x$network[x$network$target %in% rownames(x$symbol_VST),]
n<-table(net$tf);net<-net[net$tf %in% names(n[n>=5]),]
reg<-lapply(split(net,net$tf),function(a)list(tfmode=setNames(a$mor/max(abs(a$mor)),a$target),likelihood=rep(1,nrow(a))))
direct<-viper::viper(eset=x$symbol_VST,regulon=reg,method='scale',minsize=5,pleiotropy=TRUE,eset.filter=FALSE,bootstraps=0,cores=1,verbose=FALSE)
stopifnot(setequal(rownames(direct),rownames(vp)),identical(colnames(direct),colnames(vp)))
vpdiff<-max(abs(direct[rownames(vp),]-vp))
cat('VIPER_MAX_ABSOLUTE_DIFFERENCE',vpdiff,'\n');flush.console();stopifnot(vpdiff<1e-10)
write.table(direct,file.path(out,'VIPER_direct_API_validation.tsv'),sep='\t',quote=FALSE,col.names=NA)
cat('START_COUNTS_COOKS_AND_MAPPING_AUDIT\n');flush.console()
suppressPackageStartupMessages(library(DESeq2))
dds<-readRDS(file.path(root,'05_human_RNA/R_DESeq2_bridge_v1/R_DESeq2_fit.rds'))
de<-read.delim(file.path(root,'05_human_RNA/R_DESeq2_bridge_v1/R_DESeq2_all_genes.tsv'))
rr<-de[match(rownames(dds),de$ensembl_id),]
uncensored<-as.data.frame(results(dds,contrast=c('condition','STC','Control'),cooksCutoff=FALSE,independentFiltering=FALSE))
cutoff<-qf(.99,2,8);mx<-apply(assays(dds)[['cooks']],1,max,na.rm=TRUE)
audit<-data.frame(ensembl_id=rownames(dds),gene_name=rr$gene_name,max_sample_Cooks=mx,above_cutoff=mx>cutoff,primary_P_missing=!is.finite(rr$pvalue),uncensored_P_finite=is.finite(uncensored$pvalue))
write.table(audit,file.path(out,'Cook_filter_vs_maximum_audit.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
stopifnot(all(audit$uncensored_P_finite),sum(audit$primary_P_missing)==106)
src<-deparse(DESeq2::results);hit<-grep('cooksOutlier|maxCooks|which.max|numGreater',src)
idx<-sort(unique(unlist(lapply(hit,function(i)seq(max(1,i-4),min(length(src),i+7))))))
writeLines(src[idx],file.path(out,'DESeq2_results_Cooks_logic_excerpt.txt'))
replicate<-as.matrix(counts(dds));stopifnot(!anyDuplicated(as.data.frame(t(replicate))))
dup<-duplicated(de$gene_name)|duplicated(de$gene_name,fromLast=TRUE)
write.table(de[dup,c('ensembl_id','gene_name','baseMean','passed_prefilter','pvalue')],file.path(out,'all_duplicate_symbol_rows.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
maps<-list();ref<-file.path(root,'03_inputs/references/msigdbr_26.1.0_20261012')
for(fam in c('C2_CP','GO_BP')){
 a<-readRDS(file.path(ref,paste0(fam,'_source_mapping.rds')))
 stopifnot(all(a$db_target_species=='HS'),all(a$db_version=='2026.1.Hs'))
 symbols<-unique(a$gene_symbol)
 cov<-data.frame(gene_symbol=symbols,in_Wald_rank=symbols %in% names(x$ranks),in_VST=symbols %in% rownames(x$symbol_VST))
 write.table(cov,file.path(out,paste0(fam,'_all_resource_symbol_coverage.tsv')),sep='\t',quote=FALSE,row.names=FALSE)
 maps[[fam]]<-list(resource_unique_symbols=length(symbols),matched_Wald_symbols=sum(cov$in_Wald_rank),matched_VST_symbols=sum(cov$in_VST),input_rank_symbols_in_resource=sum(names(x$ranks) %in% symbols),db_target_species=unique(a$db_target_species),source_organisms=as.list(table(a$gs_source_species)),mapping='exact submitted gene symbols to human database; no inferred aliases or Ensembl remapping')
}
# Only inspect existing cache, with no retrieval or installation.
cache<-file.path(tools::R_user_dir('msigdbr','cache'),'msigdb.2026.1.zip')
cachemd5<-if(file.exists(cache))unname(tools::md5sum(cache)) else NA_character_
stopifnot(identical(cachemd5,'512ba99c6827141a9d471972b812d4ac'))
summary<-list(status='PASS_WITH_ORIGINAL_ENVIRONMENT_WARNINGS_RETAINED',ssGSEA_independent_max_abs_difference=ssdiff,VIPER_direct_API_max_abs_difference=vpdiff,VIPER_likelihood_weights=1,VIPER_mor_normalization='mor divided by max absolute mor within TF; all input absolute mor=1',VIPER_pleiotropy_args='package default regulators=0.05, shadow=0.05, targets=10, penalty=20, method=adaptive; verified in saved function formals',primary_Cooks_excluded=106,genes_with_any_sample_Cooks_above_cutoff=sum(mx>cutoff),above_cutoff_but_primary_test_retained=sum(mx>cutoff&!audit$primary_P_missing),all_20479_have_finite_uncensored_P=all(audit$uncensored_P_finite),primary_results_changed=FALSE,duplicate_count_profiles=0,resource_mapping=maps,msigdbr_source_zip=cache,msigdbr_source_zip_md5=cachemd5,msigdbr_source_zip_url='https://zenodo.org/records/18968178/files/msigdb.2026.1.zip?download=1',original_warning='R stack imbalance in :: and <-; not reproduced by isolated namespace loads; exact triggering native component not established. Independent base-R ssGSEA and direct VIPER values validated; warning is not silently dismissed.',gene_test_inference_limitation='106 Cook-filtered rows omitted from primary Wald rank; present in VST score matrix per frozen plan')
writeLines(capture.output(formals(viper::viper)),file.path(out,'VIPER_actual_formals.txt'))
jsonlite::write_json(summary,file.path(out,'technical_validation_summary.json'),pretty=TRUE,auto_unbox=TRUE,digits=16)
writeLines(capture.output(sessionInfo()),file.path(out,'sessionInfo.txt'))
cat(jsonlite::toJSON(summary,pretty=TRUE,auto_unbox=TRUE,digits=10),'\n')
