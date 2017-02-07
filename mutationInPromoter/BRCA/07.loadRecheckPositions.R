#!usr/bin/Rscript

###################
#Load library
##################
require(methods)
require(tidyverse)
require(stringr)
require(biomaRt)
######################
# Load anno data
#####################
dataRootPath = "/home/cliu18/liucj/projects/1.Mutation_calling_in_non-condig_region_through_EXOME/3.calling/BRCA_reanalysis"
recheckDir = file.path(dataRootPath, "03.somatic/03.somaticForAnalysis.recheckPositions")
tumorDir = file.path(recheckDir, "tumor")
normalDir = file.path(recheckDir, "normal")
tumorList = list.files(path = tumorDir)
normalList = list.files(path = normalDir)

maniDir = file.path(dataRootPath, "03.somatic/03.somaticForAnalysis.saveFiles")
manifest <- read_rds(path = file.path(maniDir, "manifestData.info.rds")) %>%
  dplyr::select(barcode, type, bam)

forAnalysisDir <- file.path(dataRootPath, '03.somatic/03.somaticForAnalysis.saveFiles')
# Load somatic mutation with recur >= 6
somaticMutation <- 
  read_rds(path = file.path(forAnalysisDir, 
                            'realSomaticMutation.recur5ForAnalysis.rds')) %>%
  separate(mutation, into = c("mutation", "alt"), sep = "\\/")

getPileupInfo <- function(x, path){
  read_tsv(file.path(path, x), 
           col_names = c("pileup","pos", "ref", "count", "base", "qual"),
           col_types = cols(ref = "c")) %>%
    separate(pileup, into = c("bam","chrom"), sep = ":") %>%
    mutate(bam = str_replace(bam, "/extraspace/TCGA/WXS_RAW/BRCA/regulatoryBam/",""),
           alt_pile = str_replace_all(base, pattern = "[\\,\\.\\$\\^\\-\\+\\:\\d\\@QP\"\\<\\]\\*\\!\\#\\/]","")) %>%
    mutate(bam = str_replace(bam, "tumor/", ""), alt_count = str_count(alt_pile)) %>%
    mutate(bam = str_replace(bam, "extracted\\.bam\\.mpileup","bam")) %>%
    unite(mutation, chrom, pos, sep=":") %>%
    unite(mutation, mutation, ref) %>%
    dplyr::select(-qual) %>%
    left_join(manifest, by = "bam")
}

tumorMutationBase <-
  tumorList %>%
  lapply(getPileupInfo, path = tumorDir) %>%
  bind_rows()

normalMutationBase <-
  normalList %>%
  lapply(getPileupInfo, path = normalDir) %>%
  bind_rows()

totalMutationBase <- 
  bind_rows(normalMutationBase, tumorMutationBase) 

# save the pileup info
write_tsv(totalMutationBase, path = file.path(dataRootPath, "03.somatic/03.somaticForAnalysis.saveFiles/totalMutationBase.tsv"))
write_rds(totalMutationBase, path = file.path(dataRootPath, "03.somatic/03.somaticForAnalysis.saveFiles/totalMutationBase.rds"))

totalMutationBase.count10 <- 
  totalMutationBase %>% 
  filter(count>=10)

totalMutationBase.count10 %>%
  group_by(mutation, type) %>%
  count() %>%
  View()

somaticMutation %>% 
  left_join(totalMutationBase, 
            by = c("mutation", "barcode")) %>%
  View()

# save workspace
save(list = ls(), file = file.path(recheckDir, "07.loadRecheckPositions.RData"))

