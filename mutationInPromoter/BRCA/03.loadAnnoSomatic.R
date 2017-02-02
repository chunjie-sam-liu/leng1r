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
annoDir <- '/home/cliu18/liucj/projects/1.Mutation_calling_in_non-condig_region_through_EXOME/3.calling/BRCA_reanalysis/03.somatic/01.annotation'

filteroutDir <- '/home/cliu18/liucj/projects/1.Mutation_calling_in_non-condig_region_through_EXOME/3.calling/BRCA_reanalysis/03.somatic/02.filteroutSomatic'

forAnalysisDir <- '/home/cliu18/liucj/projects/1.Mutation_calling_in_non-condig_region_through_EXOME/3.calling/BRCA_reanalysis/03.somatic/03.somaticForAnalysis'

allAnnoFiles <- list.files(annoDir, pattern="region")


# Load all anno data.
loadAnnoData <- function(x, annoDir = annoDir){
  tmp <-
  read_tsv(file = file.path(annoDir, x)) %>%
    unite(mutation, chrom, pos, sep=":") %>%
    unite(mutation, mutation, ref) %>%
    unite(mutation, mutation, alt, sep = "/") %>%
    mutate(barcode = x)
  
  # save analysis files
  tmp %>%
    filter(id == ".", ensr != ".") %>%
    write_tsv(path = file.path(forAnalysisDir, x))
  # save filterout data.
  tmp %>%
    filter(id != "." | ensr == ".") %>%
    write_tsv(path = file.path(filteroutDir, x))
    
  return(tmp)
  }

allSomaticMutation <- 
  allAnnoFiles %>%
  lapply(loadAnnoData, annoDir = annoDir) %>%
  bind_rows()

save.image(file = file.path(forAnalysisDir, "z.loadAnnoSomatic.RData"))
