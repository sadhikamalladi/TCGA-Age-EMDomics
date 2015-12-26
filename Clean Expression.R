# cleans TPM (transcripts per million) data downloaded from TCGA

# PREREQUISITES ---------------------
# expression data for each of the cancers of interests ('tissues') downloaded and placed in individual
#    folders for each cancer
# Look at fileName variable to make sure you've downloaded the right kind of data
# adds one to the raw data, multiplies it by 1E6, logged 2, quantile normalization, scale by patient
# saves the resulting normalized data in the cancer folder under 'normalized_expresison.RDS'

library(preprocessCore) # for quantile normalization...must be downloaded from Bioconductor

# path to folder containing all downloaded TPM data from TCGA
baseDir <- '~/Dropbox/EMD/Original Expression/'

# path to folder that contains individual cancer folders to save normalized expression in
outputDir <- '~/Dropbox/EMD/'

tissues<-c('BRCA','COAD','GBM','KICH','KIRC','KIRP','LGG','LUAD','LUSC','PRAD')

for (tis in tissues)
{
  print(tis)
  
  # read and adjust data
  fileName<-paste0(tis,'.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes__data.data.txt')
  brca<-read.table(paste0(baseDir,fileName),sep='\t')
  print('table read in')
  
  rownames(brca)<-brca[,1]
  brca.norm<-brca[,brca[2,]=='scaled_estimate']
  brca.raw<-brca[,brca[2,]=='raw_count']
  
  # as numeric: raw
  brca.raw<-brca.raw[-2,]
  colnames(brca.raw)<-as.character(unlist(brca.raw[1,]))
  brca.raw<-brca.raw[-1,]
  genes<-rownames(brca.raw)
  exp<-apply(brca.raw,2,function(x){as.numeric(x)})
  rownames(exp)<-genes
  
  brca.raw<-exp
  sums<-apply(brca.raw,1,function(x){sum(x>20)})
  sums<-sums > dim(brca.raw)[2]/2
  
  # as numeric: norm
  brca.norm<-brca.norm[-2,]
  colnames(brca.norm)<-as.character(unlist(brca.norm[1,]))
  brca.norm<-brca.norm[-1,]
  genes<-rownames(brca.norm)
  exp<-apply(brca.norm,2,function(x){as.numeric(x)})
  rownames(exp)<-genes
  brca.norm<-exp[sums,]
  exp<-brca.norm
  
  print('data adjusted')
  
  # adjust TPM -> read count
  # +1, *1E6, log 2, quantile normalization, scale by patient
  exp<-exp*1E6
  exp<-log2(exp+1)
  genes<-rownames(exp)
  samples<-colnames(exp)
  exp<-normalize.quantiles(exp)
  rownames(exp)<-genes
  colnames(exp)<-samples
  
  # change barcode format
  colnames(exp)<-gsub('(TCGA-[0-9A-Z]+-[0-9A-Z]+)-.*','\\1',colnames(exp))
  sample.types <- as.numeric(gsub('TCGA-[A-Z0-9]+-[A-Z0-9]+-([0-9][0-9])[A-Z]-.*','\\1',samples))
  
  print('saving RDS')
  
  saveRDS(exp,paste0(outputDir,tis,'/normalized_expression.RDS'))
}
