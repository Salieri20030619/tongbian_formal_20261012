cat('PACKAGE msigdbr:',requireNamespace('msigdbr',quietly=TRUE),'\n')
if(requireNamespace('msigdbr',quietly=TRUE)){cat('VERSION:',as.character(packageVersion('msigdbr')),'\n');print(args(msigdbr::msigdbr))}
cat('GSVA_ssgseaParam_FORMALS\n');print(formals(GSVA::ssgseaParam))
cat('decoupleR_run_viper_FORMALS\n');print(formals(decoupleR::run_viper))
cat('fgseaMultilevel_FORMALS\n');print(formals(fgsea::fgseaMultilevel))
cat('LIMMA_VERSION',as.character(packageVersion('limma')),'\n')
