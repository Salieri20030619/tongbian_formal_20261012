root<-'/mnt/BioDisk/wangmingyuan/tongbian_formal_20261012'
out<-file.path(root,'05_human_RNA/technical_validation_v1');dir.create(out,recursive=TRUE,showWarnings=FALSE)
cat('BEFORE_LOAD_GSVA\n');flush.console();loadNamespace('GSVA');cat('AFTER_LOAD_GSVA\n');flush.console()
cat('BEFORE_LOAD_DECOUPLER\n');flush.console();loadNamespace('decoupleR');cat('AFTER_LOAD_DECOUPLER\n');flush.console()
cat('BEFORE_LOAD_VIPER\n');flush.console();loadNamespace('viper');cat('AFTER_LOAD_VIPER\n');flush.console()
functions<-c('ssgsea','ssgseaDelayed','.fastRndWalk','.fastRndWalkNArm')
sink(file.path(out,'score_implementation_source.txt'))
cat('DECOUPLER RUN_VIPER\n');cat(deparse(decoupleR::run_viper),sep='\n')
cat('\nGSVA_SSGSEA_INTERNALS\n');print(grep('ssgsea|RndWalk',ls(asNamespace('GSVA'),all.names=TRUE),value=TRUE))
for(n in functions)if(exists(n,envir=asNamespace('GSVA'))){cat('\nFUNCTION',n,'\n');cat(deparse(get(n,envir=asNamespace('GSVA'))),sep='\n')}
cat('\nMSIGDBR_CACHE_FUNCTIONS\n')
for(n in c('check_cache','read_cached_rds','cache_dir'))if(exists(n,envir=asNamespace('msigdbr'))){cat('\nFUNCTION',n,'\n');cat(deparse(get(n,envir=asNamespace('msigdbr'))),sep='\n')}
sink()
cat(readLines(file.path(out,'score_implementation_source.txt')),sep='\n')
