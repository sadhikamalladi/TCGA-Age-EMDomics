# Script to extract clinical annotations for each cancer
# Fixed binning: 20 years

# PREREQUISITES -------------------------
# downloaded and normalized TPM data from TCGA for each of the cancers of interest
# normalized TPM data contained in individual folder for each cancer under file name 'normalized_expression.RDS'
# full clinical data downloaded from TCGA portal and stored in individual folder for each cancer

# directory with folders for each cancer
baseDir <- '~/Dropbox/EMD/'

# Cancers: ER+ BRCA, ER- BRCA, COAD, GBM, KIRC, KIRP, LGG, LUAD, LUSC
cancers <- c('BRCApos', 'BRCAneg', 'COAD', 'GBM', 'KIRC', 'KIRP', 'LGG', 'LUAD', 'LUSC')

# lists to keep track of the rows in the clinical files that are of interest
# age.ind is the index of the row that contains the age annotation (in days or years)
# bar.ind is the index of the row that contains the barcode for the patient (unique identifier)
# sensitive to the order of cancers list
age.ind <- as.numeric(c('14','14','3','12','9','3','7','9','7'))
bar.ind <- as.numeric(c('22','22','1','14','10','1','11','12','10'))

# maximum patient age to be included
max.age <- 80
# bin width for fixed label assignment
bin.distance <- 20

for (c in cancers) {
  
  # read in clinical data file for cancer c
  path.to.clinical <- paste0(baseDir,c,'/ClinicalFile.txt')
  clin <- read.table(path.to.clinical, sep='\t', fill=T)
  
  # index to find age/barcode indices in list
  c.ind <- grep(c,cancers)
  
  # extract age annotation
  age <- clin[age.ind[c.ind], ]
  age <- as.character(unname(unlist(age)))
  age <- as.numeric(age)
  
  # extract patient barcodes
  bar <- clin[bar.ind[c.ind], ]
  bar <- unname(unlist(bar))
  bar <- as.character(bar)
  bar <- toupper(bar)
  
  # remove missing values and unreasonable ages
  names(age) <- bar
  age <- age[!is.na(age)]
  age <- age[age<80 & age>1]
  
  # cross-reference barcodes/age with expression data availability
  path.to.exp <- paste0(baseDir,c,'/normalized_expression.RDS')
  exp <- readRDS(path.to.exp)
  age <- age[names(age) %in% colnames(exp)]
  exp <- exp[,colnames(exp) %in% names(age)]
  
  # bin patients by age
  bins <- seq(from=0, to=max.age, by=bin.distance)
  binned <- .bincode(age,bins)
  names(binned) <- names(age)
  missing <- is.na(binned)
  binned <- binned[!missing]
  age <- age[!missing]
  
  # combine age groups 0-20 and 20-40 because not many patients are placed in the first
  binned <- binned - 1
  binned[binned==0] <- 1
  
  # save annotations
  path.to.save <- paste0(baseDir,c,'/AgeClinical.RDS')
  saveRDS(binned,path.to.save)
}