#!usr/bin/Rscript

###################
#Load library
##################
require(methods)
require(tidyverse)
require(stringr)
require(VariantAnnotation) # for load vcf file
# for parallel running 
require(doParallel)
require(doMC)
registerDoMC()
require(foreach)

##################### 
#Load manifest
#####################
downloadPath <- '/extraspace/TCGA/WXS_RAW/BRCA/downloadDataList'
manifest <- read_tsv(file = file.path(downloadPath, "manifestData.info"))


#####################
#Load raw VCF
#####################
dataRootPath = "/home/cliu18/liucj/projects/1.Mutation_calling_in_non-condig_region_through_EXOME/3.calling/BRCA_reanalysis"

TrawPath = file.path(dataRootPath,"tumor/01.rawSNP")
TfilterPath = file.path(dataRootPath, "tumor/02.filterSNP")

cat("cj")