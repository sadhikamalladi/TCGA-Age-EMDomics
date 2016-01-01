# Run SAM on age data

library(samr)

baseDir <- '~/Dropbox/EMD/'

cancers <- c('BRCApos', 'BRCAneg', 'COAD', 'GBM', 'KIRC', 'KIRP', 'LGG', 'LUAD', 'LUSC', 'PRAD')

for (c in cancers) {
  path.to.exp <- paste0(baseDir,c,'/normalized_expression.RDS')
  path.to.clin <- paste0(baseDir,c,'/AgeClinical.RDS')
  
  exp <- readRDS(path.to.exp)
  clin <- readRDS(path.to.clin)
  
  # keep only common samples
  common.names <- colnames(exp)[colnames(exp)%in%names(clin)]
  clin <- clin[common.names]
  exp <- exp[,common.names]
  
  if (min(clin) != 1)
    clin <- clin-1
  
  sam <- SAM(x=exp, y=clin, resp.type='Multiclass')
  
  saveRDS(sam,paste0(baseDir,c,'/SAM.RDS'))
}