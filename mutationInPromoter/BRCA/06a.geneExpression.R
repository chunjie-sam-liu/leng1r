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


dataRootPath = "/home/cliu18/liucj/projects/1.Mutation_calling_in_non-condig_region_through_EXOME/3.calling/BRCA_reanalysis"

forAnalysisDir <-
  file.path(dataRootPath, '03.somatic/03.somaticForAnalysis.saveFiles')

geneExpressionDir <-
  file.path(dataRootPath, "03.somatic/06.geneExpression")

# Load mutation and target genes
targetGeneDir <-
  file.path(dataRootPath,
            "03.somatic/04.targetGenesForMutation")

# get target coding genes 
mutation_target <-
  read_rds(path = file.path(targetGeneDir,
                            "nearest_target_coding_genes_3.rds")) %>%
  tbl_df() %>%
  rename(ensid = ensg)

target_gene <- 
  mutation_target %>%
  select(ensid, symbol, entrezgene) %>%
  distinct()

rnaseq <- 
  read_tsv("/extraspace/TCGA/WXS_RAW/BRCA/RNA-seq/brca.expression.tidy.tsv") %>%
  filter(!str_detect(type, "NA")) %>%
  distinct() %>%
  inner_join(target_gene, by = "ensid")


candidate_mutation <- 
  read_rds(file.path(forAnalysisDir, 
                     "recur3_candidate_merged_mutation.rds")) %>% tbl_df()

candidate_mutation_target <-
  mutation_target %>%
  semi_join(candidate_mutation, by = "mutation") 

# save
candidate_mutation_target %>% 
  write_tsv(path = file.path(targetGeneDir, "recur3_candidate_mutation_target.tsv"))
candidate_mutation_target %>%
  write_rds(path = file.path(targetGeneDir, "recur3_candidate_mutation_target.rds"))


candidate_mutation_target_expression <-
  candidate_mutation %>% 
  dplyr::select(mutation, barcode, type) %>%
  left_join(candidate_mutation_target, by = "mutation") %>%
  dplyr::select(mutation, 
                barcode, 
                ensid, 
                type = type.x, 
                symbol) %>%
  left_join(rnaseq, by = c("barcode", "ensid")) %>%
  rename(type = type.x, symbol = symbol.x) %>%
  dplyr::select(-type.y, -symbol.y) %>%
  filter(!is.na(expression))

# broom for wilcox.test model
# expression_model <- 
#   candidate_mutation_target_expression %>%
#   filter(!is.na(expression)) %>%
#   group_by(mutation, ensid, symbol) %>%
#   do(tidy(wilcox.test(expression ~ as.factor(type), data = .))) %>%
#   ungroup()%>%
#   mutate(p.value = round(p.value,3), fdr = round(qvalue(p.value)$qvalues, 3))



candidate_mutation_target_expression %>%
  filter(!is.na(expression)) %>% 
  group_by(mutation, ensid, symbol) %>% 
  by_slice(
    safely(
      ~ wilcox.test(expression ~ as.factor(type), data = .)
    )
  ) %>%
  mutate(model = map(.out, function(l){tidy(l$result)})) %>%
  unnest(model) %>% 
  select(mutation, ensid, symbol, p.value) %>% 
  filter(!is.nan(p.value), p.value != 1) %>% 
  mutate(p.value = round(p.value, 3), fdr = round(qvalue(p.value)$qvalues, 3)) %>% 
  filter(p.value < 0.05) -> expression_model

candidate_mutation_target_expression %>% 
  inner_join(expression_model, by = c("mutation", "ensid")) %>% 
  rename(symbol = symbol.x) ->
  candidate_mutation_target_expression_model

all_expression_barplot_p0.05<-
  candidate_mutation_target_expression_model %>% 
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
                )
    
  
  print(all_expression_barplot_p0.05)
  
  all_expression_barplot_p0.05 %>%
    ggsave(filename = file.path(geneExpressionDir, 
                                "02.recur3_all_expression_barplot_p0.05.png"),
           device = "png", width = 15, height = 15)
  
  all_expression_barplot_p0.05 %>%
    write_rds( path = file.path(geneExpressionDir, 
                                "02.recur3_all_expression_barplot_p0.05.rds"))


candidate_p0.01 <- 
  candidate_mutation_target_expression_model %>% 
  dplyr::select(-c(symbol.y :fdr))

mutation_ensid <-
  candidate_p0.01 %>%  
  distinct(mutation, ensid)

candidate_mutation_with_normal <-
  rnaseq %>%
  filter(type == "normal") %>%
  semi_join(candidate_p0.01, by = "ensid") %>%
  mutate(type = "NM") %>%
  left_join(mutation_ensid, by = "ensid") %>% 
  dplyr::select(mutation, barcode, ensid, type, symbol, expression) %>%
  bind_rows(candidate_p0.01)
  
  
candidate_mutation_with_normal %>%
  write_tsv(file.path(geneExpressionDir, "recur3_candidate_mutation_with_normal.tsv"))
candidate_mutation_with_normal %>% 
  write_rds(file.path(geneExpressionDir, "recur3_candidate_mutation_with_normal.rds"))
  
  
candidate_mutation_expression_model<-
  candidate_mutation_with_normal %>%
  filter(!is.na(expression)) %>%
  group_by(mutation, ensid, symbol) %>%
  do(tidy(kruskal.test(expression ~ as.factor(type), data = .))) %>%
  ungroup()
  

label_observation_count <- function(x){
  return(c(y = 50, label = length(x)))
}


candidate_mutation_2_barplot <-
  candidate_mutation_with_normal %>% 
  inner_join(candidate_mutation_expression_model, by = c("mutation", "ensid")) %>%
  rename(symbol = symbol.x) %>%
  filter(p.value < 0.05) %>%
  ggplot(aes(x = type, y = expression)) +
  geom_boxplot(aes(color = type), outlier.shape = NA) +
  scale_color_brewer(palette = "Dark2") +
  stat_summary(fun.data = label_observation_count, geom = "text") + 
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
                              "03.reucr3_candidate_mutation_2_barplot.png"),
         device = "png", width = 10, height = 15)

candidate_mutation_2_barplot %>%
  write_rds( path = file.path(geneExpressionDir, 
                              "03.recur3_candidate_mutation_2_barplot.rds"))  



#
  
  
  
  
  
  
  