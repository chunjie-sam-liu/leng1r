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

forAnalysisDir <- file.path(dataRootPath, '03.somatic/03.somaticForAnalysis.saveFiles')

# Load rds data
# The somatic mutation is based on the paired case minus
# Considering lack of coverage on the normal sample
# And to ensure the true positive of somatic mutation
# the truePositiveSomaticMutation = somaticMutation - normalMutation
filterSomaticMutation <-
  readRDS(file = file.path(forAnalysisDir, 
                           "filteredSomaticMutation.rds")) %>%
  filter(id == ".", ensr != ".")

# Load normal mutation
normalMutation <- 
  readRDS(file = file.path(forAnalysisDir, 
                           "allNormalFilterPositions.rds")) %>%
  group_by(mutation) %>%
  filter(n() >= 2) %>%
  ungroup()
  

realSomaticMutation <- 
  filterSomaticMutation %>% 
  anti_join(normalMutation, by = "mutation")

# Save it for stat check
realSomaticMutation %>%
  write_tsv(path = file.path(forAnalysisDir, "realSomaticMutation.tsv"))
realSomaticMutation %>%
  write_rds(path = file.path(forAnalysisDir, "realSomaticMutation.rds"))

######################
# Stat for recurrent##
######################
realSomaticMutation.recurStat <-
  realSomaticMutation %>%
  group_by(mutation, feature) %>%
  count() %>%
  rename(recurrent = n) %>%
  filter(recurrent > 2)

realSomaticMutation.recurStat.barplot <- 
  realSomaticMutation.recurStat %>%
  ggplot(aes(x = as.factor(recurrent))) +
  geom_bar(aes(fill = feature)) + 
  geom_text(stat = 'count', aes(label = ..count..), 
            vjust = -0.2, hjust = 0.5) + 
  theme_minimal() + 
  xlab("Recurrent") + 
  ylab("Count") + 
  ggtitle(paste("Total # of Somatic Mutation:", sum(realSomaticMutation.recurStat$recurrent))) +
  scale_fill_brewer(palette = "Dark2")

print(realSomaticMutation.recurStat.barplot)
ggsave(filename = file.path(forAnalysisDir, "01.SomaticMutationStat.png"),
       device = "png")
realSomaticMutation.recurStat.barplot %>%
  write_rds( path = file.path(forAnalysisDir, "01.SomaticMutationStat.rds"))

########################
#Filter recurrent > 5###
########################
realSomaticMutation.recur5ForAnalysis <-
  realSomaticMutation %>%
  group_by(mutation, feature) %>%
  count() %>%
  rename(recurrent = n) %>%
  filter(recurrent > 5) %>% 
  ungroup() %>%
  inner_join(realSomaticMutation, by = c("mutation", "feature")) 

# Save it for further analysis
realSomaticMutation.recur5ForAnalysis %>%
  write_tsv(path = file.path(forAnalysisDir, "realSomaticMutation.recur5ForAnalysis.tsv"))
realSomaticMutation.recur5ForAnalysis %>%
  write_rds(path = file.path(forAnalysisDir, "realSomaticMutation.recur5ForAnalysis.rds"))

# write recur5 to tsv for revising rate and rate
saveThePositions <- function(x){
  fileName <- paste(x["chrom"], x["pos"], sep = "_")
  print(name)
}

realSomaticMutation.recur5ForAnalysis %>%
  distinct(mutation) %>%
  separate(mutation, into = c("chrom", "pos", "ref", "alt")) %>%
  write_tsv(path = file.path(forAnalysisDir, "realSomaticMutation.recur5ForAnalysis.refinePositions.tsv"), col_names = F)


save(list = ls(), file = file.path(forAnalysisDir, "04.somaticMutationAnalysis.RData"))
