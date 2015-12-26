# Run EMDomics on normalized expression data

library(EMDomics) # downloaded from Bioconductor

baseDir <- '~/Dropbox/EMD/'

cancers <- c('BRCApos', 'BRCAneg', 'COAD', 'GBM', 'KICH', 'KIRC', 'KIRP', 'LGG', 'LUAD', 'LUSC', 'PRAD')

for (c in cancers) {
  path.to.exp <- paste0(baseDir,c,'/normalized_expression.RDS')
  path.to.clin <- paste0(baseDir,c,'/AgeClinical.RDS')
  
  exp <- readRDS(path.to.exp)
  clin <- readRDS(path.to.clin)
  
  emd <- calculate_emd(exp,clin)
  
  output.path <- paste0(baseDir,c,'/EMD.RDS')
  saveRDS(emd,output.path)
}