#!usr/bin/Rscript

###################
#Load library
##################
require(methods)
require(tidyverse)
require(stringr)
require(broom)
require(qvalue)
######################
# Load anno data
#####################
dataRootPath = "/home/cliu18/liucj/projects/1.Mutation_calling_in_non-condig_region_through_EXOME/3.calling/BRCA_reanalysis"

# Load mutation and target genes
targetGeneDir <- 
  file.path(dataRootPath, 
            "03.somatic/04.targetGenesForMutation")

# The expression.rds was extracted by the target genes.
mutationTarget <- 
  read_rds(path = file.path(targetGeneDir, 
                            "nearestTargetProteinCodingGenes.rds"))

# get expression data
rnaseqPath <- "/extraspace/TCGA/WXS_RAW/BRCA/RNA-seq/brca.expression.tidy.44genes.rds"
rnaseq <- 
  read_rds(rnaseqPath) %>%
  filter(!is.na(type))

geneExpressionDir <- 
  file.path(dataRootPath, "03.somatic/06.geneExpression")

rnaseq <- 
  mutationTarget %>%
  tbl_df() %>%
  dplyr::select(ensid = ensg, symbol) %>%
  filter(symbol != "") %>%
  distinct(ensid, symbol) %>%
  inner_join(rnaseq, by = "ensid")

rnaseq.wilcox.test <- 
  rnaseq %>% 
  group_by(ensid, symbol) %>%
  do(tidy(wilcox.test(expression ~ type, data = .))) 

# rnaseq.wilcox.test %>%
#   ggplot(aes(p.value)) +
#   geom_histogram(binwidth = 0.5)

rnaseq.wilcox.test <-  
  rnaseq.wilcox.test %>%
  ungroup() %>%
  mutate(q.value = qvalue(p.value)$qvalues) %>%
  filter(q.value < 0.01, p.value < 0.05)

rnaseq.tumor.normal.barplot <- 
  rnaseq %>% 
  filter(symbol != "") %>%
  inner_join(rnaseq.wilcox.test, by = c("ensid", "symbol")) %>%
  ggplot(aes(type, expression)) +
  geom_boxplot() +
  theme_minimal() +
  facet_wrap( ~ ensid + symbol + p.value + q.value, 
              ncol = 4, scales = "free" )
# print(rnaseq.tumor.normal.barplot)
ggsave(filename = file.path(geneExpressionDir, 
                            "01.rnaseq.tumor.normal.barplot.png"),
       device = "png", width = 10, height = 20)
realSomaticMutation.recurStat.barplot %>%
  write_rds( path = file.path(geneExpressionDir, 
                              "01.rnaseq.tumor.normal.barplot.rds"))



# Save workspace
save(list = ls(), file = file.path(geneExpressionDir, "06.geneExpression.RData"))