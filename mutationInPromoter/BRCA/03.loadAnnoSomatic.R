#!usr/bin/Rscript

###################
#Load library
##################
require(methods)
require(tidyverse)
require(stringr)

######################
# Load anno data
#####################
dataRootPath = "/home/cliu18/liucj/projects/1.Mutation_calling_in_non-condig_region_through_EXOME/3.calling/BRCA_reanalysis"

annoDir <- file.path(dataRootPath,'03.somatic/01.annotation')

filteroutDir <- file.path(dataRootPath, '03.somatic/02.filteroutSomatic')

forAnalysisDir <- file.path(dataRootPath, '03.somatic/03.somaticForAnalysis')

# Load to save to rds for check position
filterPositionDir <- file.path(dataRootPath, '02.filtersnp')
filterPositionFiles <- list.files(filterPositionDir, pattern = "bam.SNP.vcf.filter")

# Load mutation statistics
somaticMutationStat <- 
  read_tsv(file = file.path(dataRootPath, 'somaticStat.tsv'))

allAnnoFiles <- list.files(annoDir, pattern="avinput\\.dbsnp\\.sortByChrom\\.region")

# All tsv and rds save to this directory for next analysis
saveFilesDir <- file.path(dataRootPath, "03.somatic/03.somaticForAnalysis.saveFiles")


# Load all anno data.
loadAnnoData <- function(x, annoDir = annoDir){
    
  tmp <-
  read_tsv(file = file.path(annoDir, x)) %>%
    unite(mutation, chrom, pos, sep=":") %>%
    unite(mutation, mutation, ref) %>%
    unite(mutation, mutation, alt, sep = "/")
  
  # save analysis files
  tmp %>%
    filter(id == ".", ensr != ".") %>%
    write_tsv(path = file.path(forAnalysisDir, paste(x,"analysis", sep = ".")))
  # save filterout data.
  tmp %>%
    filter(id != "." | ensr == ".") %>%
    write_tsv(path = file.path(filteroutDir, paste(x,"filterout", sep = ".")))
    
  return(tmp)
  }

# Load all somatic mutation
# Including dbSNP position
allSomaticMutation <- 
  allAnnoFiles %>%
  lapply(loadAnnoData, annoDir = annoDir) %>%
  bind_rows()

# count analysis
somaticMutationStat <- 
  allSomaticMutation %>%
  filter(id == ".", ensr != ".") %>%
  group_by(barcode) %>%
  count() %>% rename(analysis = n ) %>% 
  right_join(somaticMutationStat, by = 'barcode') 
# count filterout
somaticMutationStat <- 
  allSomaticMutation %>%
  filter(id != "." | ensr == ".") %>%
  group_by(barcode) %>%
  count() %>% rename(filterout = n ) %>% 
  right_join(somaticMutationStat, by = 'barcode') 

# write to dataRootPath
somaticMutationStat %>%
  write_tsv(path = file.path(saveFilesDir, 'somaticMutationCountStat.tsv'))
  
# use saveRDS to save allSomaticMutation
# use readRDS to load allSomaticMutation
allSomaticMutation %>%
  saveRDS(file = file.path(saveFilesDir, "filteredSomaticMutation.rds"))
# for manual check
allSomaticMutation %>%
  write_tsv(path = file.path(saveFilesDir, "filteredSomaticMutation.tsv"))


# Manifest
downloadPath <- '/extraspace/TCGA/WXS_RAW/BRCA/downloadDataList'

typeName <- c("Primary Tumor" = "tumor", "Blood Derived Normal" = "normal", "Solid Tissue Normal" = "normal", "Metastatic" = "tumor")

manifest <- 
  read_tsv(file = file.path(downloadPath, "manifestData.info")) %>%
  mutate(type = plyr::revalue(type, typeName))
# save modified manifest file
manifest %>%
  write_tsv(path = file.path(saveFilesDir, "manifestData.info.tsv"))
manifest %>%
  saveRDS(file = file.path(saveFilesDir, "manifestData.info.rds"))

# import filter position files and save it to normal.rds and tumor.rds
loadPositionFunction <- 
  function(x, filterPositionDir = filterPositionDir){
    tmp <-
      read_tsv(file = file.path(filterPositionDir, x)) 
  }

allFilterPositions <-
  filterPositionFiles %>%
  lapply(loadPositionFunction, filterPositionDir = filterPositionDir) %>%
  bind_rows()

allFilterPositions <- 
  manifest %>%
  dplyr::select(barcode, type) %>%
  right_join(allFilterPositions, by = 'barcode')

# save normal
allFilterPositions %>%
  filter(type == "normal") %>%
  write_tsv(path = file.path(saveFilesDir, "allNormalFilterPositions.tsv"))
allFilterPositions %>%
  filter(type == "normal") %>%
  saveRDS(file = file.path(saveFilesDir, "allNormalFilterPositions.rds"))
# save tumor
allFilterPositions %>%
  filter(type == "tumor") %>%
  write_tsv(path = file.path(saveFilesDir, "allTumorFilterPositions.tsv"))
allFilterPositions %>%
  filter(type == "tumor") %>%
  saveRDS(file = file.path(saveFilesDir, "allTumorFilterPositions.rds"))

# Save the working space
save.image(file = file.path(forAnalysisDir, "z.loadAnnoSomatic.RData"))