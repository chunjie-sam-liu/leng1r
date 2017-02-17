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
  semi_join(candidate_mutation, by = "mutation")

# save candidate mutation target genes.
candidate_mutation_target %>% 
  write_tsv(path = file.path(targetGeneDir, "candidate_mutation_target.tsv"))
candidate_mutation_target %>%
  write_rds(path = file.path(targetGeneDir, "candidate_mutation_target.rds"))


# get expression data
rnaseqPath <- "/extraspace/TCGA/WXS_RAW/BRCA/RNA-seq/brca.expression.tidy.44genes.rds"
rnaseq <- 
  read_rds(rnaseqPath) %>%
  filter(!is.na(type)) %>%
  distinct()

# get mutation target gene expression with samples
candidate_mutation_target_expression <-
  candidate_mutation %>% 
  dplyr::select(mutation, barcode, type) %>%
  left_join(candidate_mutation_target, by = "mutation") %>%
  dplyr::select(mutation, 
                barcode, 
                ensid = ensg, 
                type = type.x, 
                symbol) %>%
  left_join(rnaseq, by = c("barcode", "ensid")) %>%
  rename(type = type.x) %>%
  dplyr::select(-type.y) %>%
  filter(!is.na(expression))

# broom for wilcox.test model
expression_model <- 
  candidate_mutation_target_expression %>%
  filter(!is.na(expression)) %>%
  group_by(mutation, ensid, symbol) %>%
  do(tidy(wilcox.test(expression ~ type, data = .))) %>%
  ungroup()%>%
  mutate(p.value = round(p.value,3), fdr = round(qvalue(p.value)$qvalues, 3))

# filter p.value < 0.01
all_expression_barplot_p0.05 <-
  candidate_mutation_target_expression %>% 
  inner_join(expression_model, by = c("mutation", "ensid")) %>%
  rename(symbol = symbol.x) %>%
  filter(p.value < 0.05) %>%
  ggplot(aes(x = type, y = expression)) +
  geom_boxplot(aes(color = type), outlier.shape = NA) +
  scale_color_brewer(palette = "Set1") +
  theme_bw() +
  labs(x = "Type", y = "Gene Expression (FPKM)", title = "Target Gene Expression") +
  theme(strip.background = element_rect(fill = 'white', color = 'black'),
        panel.grid = element_blank(),
        axis.line = element_line(color = 'black'),
        plot.title = element_text(hjust = 0.5)) + 
  facet_wrap( ~ mutation + 
                paste(ensid, symbol, sep = "/") +
                paste("p-value", p.value, sep = ": ")# + paste("FDR", fdr, sep = ": ")
              , 
              ncol = 4, scales = "free")

print(all_expression_barplot_p0.05)

all_expression_barplot_p0.05 %>%
  ggsave(filename = file.path(geneExpressionDir, 
                              "02.all_expression_barplot_p0.05.png"),
         device = "png", width = 5, height = 5)

all_expression_barplot_p0.05 %>%
  write_rds( path = file.path(geneExpressionDir, 
                              "02.all_expression_barplot_p0.05.rds"))  

# I got two candidate mutation
# Add normal sample expression to compare MT WT NM expression
# c("ENSG00000103047","ENSG00000048392")
candidate_p0.01 <- 
  candidate_mutation_target_expression %>% 
  inner_join(expression_model, by = c("mutation", "ensid")) %>%
  rename(symbol = symbol.x) %>%
  filter(p.value < 0.05) %>%
  dplyr::select(-c(symbol.y:fdr))

candidate_mutation_with_normal <-
  rnaseq %>%
  filter(type == "normal") %>%
  semi_join(candidate_p0.01, by = "ensid") %>%
  mutate(symbol = ifelse(ensid == "ENSG00000048392", 
                            "RRM2B",
                            "TANGO6"),
         mutation = ifelse(ensid == "ENSG00000048392",
                           "chr8:102395112_A/C",
                           "chr16:68738315_C/T"),
         type = "NM") %>%
  dplyr::select(mutation, barcode, ensid, type, symbol, expression) %>%
  bind_rows(candidate_p0.01)
candidate_mutation_with_normal %>%
  write_tsv(file.path(geneExpressionDir, "candidate_mutation_with_normal.tsv"))
candidate_mutation_with_normal %>% 
  write_rds(file.path(geneExpressionDir, "candidate_mutation_with_normal.rds"))


candidate_mutation_expression_model<-
  candidate_mutation_with_normal %>%
  filter(!is.na(expression)) %>%
  group_by(mutation, ensid, symbol) %>%
  do(tidy(kruskal.test(expression ~ as.factor(type), data = .))) %>%
  ungroup()

candidate_mutation_2_barplot <-
  candidate_mutation_with_normal %>% 
  inner_join(candidate_mutation_expression_model, by = c("mutation", "ensid")) %>%
  rename(symbol = symbol.x) %>%
  filter(p.value < 0.05) %>%
  ggplot(aes(x = type, y = expression)) +
  geom_boxplot(aes(color = type), outlier.shape = NA) +
  scale_color_brewer(palette = "Dark2") +
  theme_bw() + 
  labs(x = "Type", y = "Gene Expression (FPKM)", title = "Target Gene Expression") +
  theme(strip.background = element_rect(fill = 'white', color = 'black'),
        panel.grid = element_blank(),
        axis.line = element_line(color = 'black'),
        plot.title = element_text(hjust = 0.5)) + 
  scale_x_discrete(limits = c("NM", "WT","MT")) +
  facet_wrap( ~ mutation + 
                paste(ensid, symbol, sep = "/"), 
              ncol = 4, scales = "free")

print(candidate_mutation_2_barplot)

candidate_mutation_2_barplot %>%
  ggsave(filename = file.path(geneExpressionDir, 
                              "03.candidate_mutation_2_barplot.png"),
         device = "png", width = 5, height = 5)

candidate_mutation_2_barplot %>%
  write_rds( path = file.path(geneExpressionDir, 
                              "03.candidate_mutation_2_barplot.rds"))  


###########################################
# over all candidate gene expression test#
##########################################
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
  do(tidy(wilcox.test(expression ~ type, data = .))) %>%
  ungroup() %>%
  mutate(q.value = round(qvalue(p.value)$qvalues,3), 
         p.value = round(p.value, 3)) %>%
  filter(q.value < 0.01, p.value < 0.05)

rnaseq.tumor.normal.barplot <- 
  rnaseq %>% 
  filter(symbol != "") %>%
  inner_join(rnaseq.wilcox.test, by = c("ensid", "symbol")) %>%
  ggplot(aes(type, expression)) +
  geom_boxplot(aes(color=type), outlier.colour = NA) +
  scale_color_brewer(palette = "Dark2") + 
  theme_minimal() +
  facet_wrap( ~ ensid + 
                symbol + 
                paste("p-value:", p.value) + 
                paste("fdr: ", q.value), 
              ncol = 4, scales = "free" )
#print(rnaseq.tumor.normal.barplot)
rnaseq.tumor.normal.barplot %>%
  ggsave(filename = file.path(geneExpressionDir, 
                            "01.rnaseq.tumor.normal.barplot.png"),
       device = "png", width = 10, height = 20)
realSomaticMutation.recurStat.barplot %>%
  write_rds( path = file.path(geneExpressionDir, 
                              "01.rnaseq.tumor.normal.barplot.rds"))



# Save workspace
save(list = ls(), file = file.path(geneExpressionDir, "06.geneExpression.RData"))
load(file.path(geneExpressionDir, "06.geneExpression.RData"))

