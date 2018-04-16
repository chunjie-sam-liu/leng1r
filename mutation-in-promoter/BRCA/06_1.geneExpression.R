#!usr/bin/Rscript

###################
#Load library
##################
require(methods)
require(tidyverse)
require(stringr)
require(broom)
require(qvalue)
library(reshape2)
######################
# Load anno data
#####################
dataRootPath = "/home/cliu18/liucj/projects/1.Mutation_calling_in_non-condig_region_through_EXOME/3.calling/BRCA_reanalysis"
forAnalysisDir <- file.path(dataRootPath, '03.somatic/03.somaticForAnalysis.saveFiles')

geneExpressionDir <- 
  file.path(dataRootPath, "03.somatic/06.geneExpression")
# Load mutation and target genes
targetGeneDir <- 
  file.path(dataRootPath, 
            "03.somatic/04.targetGenesForMutation")

# The expression.rds was extracted by the target genes.
mutationTarget <- 
  read_rds(path = file.path(targetGeneDir, 
                            "nearestTargetProteinCodingGenes.rds")) %>% tbl_df()

candidate_mutation <- 
  read_rds(file.path(forAnalysisDir, 
                     "candidate_merged_mutation.rds")) %>% tbl_df()

candidate_mutation_target <-
  mutationTarget %>%
  semi_join(candidate_mutation, by = "mutation") %>%
  rename(ensid = ensg) %>%
  dplyr::select(ensid, symbol, entrezgene) %>%
  distinct()

  

# get expression data
rnaseqPath <- "/extraspace/TCGA/WXS_RAW/BRCA/RNA-seq/brca.expression.tidy.21genes.rds"

rnaseq <- 
  read_rds(rnaseqPath) %>%
  filter(!is.na(type)) %>%
  distinct() %>%
  left_join(candidate_mutation_target, by = 'ensid')
  

# get mutation target gene expression with samples
candidate_mutation %>% 
  select(barcode, type) %>% 
  distinct() -> 
  candidate_mutation_type

rnaseq %>%
  left_join(candidate_mutation_type, by = 'barcode') %>%
  mutate(type.y = ifelse(type.x == "normal", "NM", type.y)) %>%
  filter(!is.na(type.y)) ->
  candidate_mutation_target_expression

candidate_mutation_target_expression %>%
  filter(type.y != "WT") ->
  candidate_mutation_target_expression_only_with_normal
  
candidate_mutation_target_expression_only_with_normal %>%
  filter(!is.na(expression)) %>%
  group_by(ensid, symbol) %>%
  do(tidy(wilcox.test(expression ~ type.y, data = .))) %>%
  ungroup()%>%
  mutate(p.value = signif(p.value,3), fdr = signif(qvalue(p.value)$qvalues, 3)) ->
  candidate_mutation_target_expression_only_with_normal_model

candidate_mutation_target_expression_only_with_normal %>% 
  inner_join(candidate_mutation_target_expression_only_with_normal_model, by = c("ensid")) %>%
  rename(symbol = symbol.x) %>%
  filter(p.value < 0.05, fdr < 0.05) %>%
  ggplot(aes(x = type.y, y = expression)) +
  geom_boxplot(aes(color = type.y), outlier.shape = NA) +
  scale_color_brewer(palette = "Set1") +
  theme_bw() +
  labs(x = "Type", y = "Gene Expression (FPKM)", title = "Target Gene Expression") +
  theme(strip.background = element_rect(fill = 'white', color = 'black'),
        panel.grid = element_blank(),
        axis.line = element_line(color = 'black'),
        plot.title = element_text(hjust = 0.5)) + 
  facet_wrap( ~ paste(ensid, symbol, sep = "/") +
                paste("p-value", p.value, sep = ": "), 
              ncol = 4, scales = "free") ->
  normal_expression_barplot

ggsave(plot = normal_expression_barplot, filename = file.path(geneExpressionDir, "tumor_MT_vs_normal.png"),device = "png", width = 10, height = 10)
  
  