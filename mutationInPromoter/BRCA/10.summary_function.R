#!usr/bin/Rscript

###################
#Load library
##################
require(methods)
require(tidyverse)
require(stringr)
library(latex2exp)

# set data path
save_files_dir <- "/home/cliu18/liucj/projects/1.Mutation_calling_in_non-condig_region_through_EXOME/3.calling/BRCA_reanalysis/03.somatic/03.somaticForAnalysis.saveFiles"
recheck_dir <- "/home/cliu18/liucj/projects/1.Mutation_calling_in_non-condig_region_through_EXOME/3.calling/BRCA_reanalysis/03.somatic/03.somaticForAnalysis.recheckPositions"
target_gene_dir <- "/home/cliu18/liucj/projects/1.Mutation_calling_in_non-condig_region_through_EXOME/3.calling/BRCA_reanalysis/03.somatic/04.targetGenesForMutation"

# output path
out_table_dir <- "/home/cliu18/liucj/github/RstudioWithGit/mutationInPromoter/BRCA/summary/2017_01_13/tables"
out_pic_dir <- "/home/cliu18/liucj/github/RstudioWithGit/mutationInPromoter/BRCA/summary/2017_01_13/figs"
##########
# Table 1. Detailed information regarding included TCGA-BRCA including UUID codes for GDC-portal.
##########
manifest <- 
  read_rds(file.path(save_files_dir, 'manifestData.info.rds')) %>%
  rename(file_size = size, regulatory_file_size = extractedSize)

# Load hard filtered mutation number for every bam
hard_filter_dir <- "/home/cliu18/liucj/projects/1.Mutation_calling_in_non-condig_region_through_EXOME/3.calling/BRCA_reanalysis/02.filtersnp"

load_mutation_number <- function(x){
  # print(x['bam'])
  num <- as.numeric(R.utils::countLines(
    file.path(hard_filter_dir,
              paste(x['bam'], "SNP.vcf.filter", sep = "."))))
  x['hard_filter_mutation_num'] <- num
  x %>% 
    tbl_df() %>% 
    rownames_to_column() %>%
    spread(rowname, value, convert = T)
}

manifest_hard <-
  manifest %>%
  apply(1, load_mutation_number) %>%
  bind_rows() %>%
  dplyr::select(case, barcode, type, uuid,
                file_size, regulatory_file_size,
                hard_filter_mutation_num)
manifest_hard %>% 
  write_tsv(file.path(out_table_dir, 'table1_manifest_hard.tsv'))
manifest_hard %>%
  write_rds(file.path(out_table_dir, 'table1_manifest_hard.rds'))
#########
# Figure 1 Exome reads located in the regulatory region
#######
figure1_file_size_comparison_point <- 
  manifest_hard %>% 
  ggplot(aes(x = file_size / 10 ^ 10, 
             y = regulatory_file_size / 10 ^ 10)) +
  geom_point() + 
  geom_smooth(method = loess) + 
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black")
        ) + 
  labs(title = TeX("File Size Comparison ($10^{10}$ bytes)"),
       x= TeX('Raw WXS File Size'),
       y = TeX('Regulatory File Size'))
ggsave(file.path(out_pic_dir, "figure1_file_size_comparison_point.png"), plot = figure1_file_size_comparison_point, device = 'png', width = 7, height = 5)


########
#Figure 2 non-coding number vs file size ratio
#######
figure2_raw_mutation_number <-
  manifest_hard %>%
  ggplot(aes(x = regulatory_file_size/file_size, 
             y = log10(hard_filter_mutation_num))) +
  geom_point(aes(color = type)) + 
  scale_color_brewer(palette = "Set1") +
  scale_x_continuous(breaks = seq(0,.5, by = 0.02)) +
  scale_y_continuous(breaks = seq(3, 5.5, by = 0.2)) +
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black")
  ) +
  labs(title = "Mutation Number vs. Reads Ratio", 
       x = TeX('Regulatory / Raw File Size'),
       y = TeX('$log_{10}$ (# of Mutations)'))
ggsave(file.path(out_pic_dir, "figure2_raw_mutation_number.png"), plot = figure2_raw_mutation_number, device = 'png', width = 7, height = 5)



# Load mutation stata
# somaticMutationCountStat.tsv
somaticMutationCountStat <- 
  read_tsv(file = file.path(save_files_dir, 'somaticMutationCountStat.tsv')) %>%
  dplyr::select(-somatic) %>%
  rename(dbSNP = filterout, somatic = analysis)

somatic_filteration_table <-
  manifest_hard %>% 
  filter(type == 'tumor') %>%
  right_join(somaticMutationCountStat, by = 'barcode') %>%
  distinct(case, .keep_all = T)

##############
# table 2 dbsnp and feature annotation filter
#######
table2_mutation_annotation_filter <-
  manifest_hard %>%
  filter(type == "normal") %>%
  distinct(case, .keep_all = T) %>%
  right_join(somatic_filteration_table, by = 'case') %>%
  transmute(case = case, 
            tumor_barcode = barcode.y,
            tumor_uuid = uuid.y,
            tumor_file_size = file_size.y,
            tumor_regulatory_file_size = regulatory_file_size.y,
            tumor_hard_filter_mutation_num = hard_filter_mutation_num.y,
            normal_barcode = barcode.x,
            normal_uuid = uuid.x,
            normal_file_size = file_size.x,
            normal_regulatory_file_size = regulatory_file_size.x,
            normal_hard_filter_mutation_num = hard_filter_mutation_num.x,
            somatic_mutation = dbSNP + somatic,
            dbSNP = dbSNP,
            somatic_mutation_filter = somatic)

table2_mutation_annotation_filter %>%
  write_tsv(file.path(out_table_dir, 'table2_mutation_annotation_filter.tsv'))
table2_mutation_annotation_filter %>%
  write_rds(file.path(out_table_dir, 'table2_mutation_annotation_filter.rds'))
#######
# Figure 3 somatic mutation of tumor/normal regulatory file size 
######3
figure3_somatic_mutation_of_tumor_normal_file_size_point <-
  table2_mutation_annotation_filter %>%
  filter(somatic_mutation_filter < 1000) %>%
  ggplot(aes(x = tumor_regulatory_file_size / normal_regulatory_file_size, 
             y = somatic_mutation_filter)) +
  geom_point() + 
  scale_y_continuous(breaks = seq(0, 3600, by = 250)) +
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black")
  ) + 
  labs(title = TeX(""),
       x= TeX('Tumor / Normal Regulatory File Size'),
       y = TeX('# of Novel Somatic Mutation'))
ggsave(file.path(out_pic_dir, "figure3_somatic_mutation_of_tumor_normal_file_size_point.png"), plot = figure3_somatic_mutation_of_tumor_normal_file_size_point, device = 'png', width = 7, height = 5)

#######
#Figure 4 somatic mutation distributation
#######
figure4_somatic_mutation_feature_annotation_stata_barplot <- 
  read_rds(file.path(save_files_dir, '01.SomaticMutationStat.rds')) 
ggsave(file.path(out_pic_dir, "figure4_somatic_mutation_feature_annotation_stata_barplot.png"), 
       plot = figure4_somatic_mutation_feature_annotation_stata_barplot, device = 'png', width = 7, height = 5)


#######################
# Load rda 02-15-2017
#load(file.path(out_table_dir, "summary_function.rda"))
######################

my_theme <- theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black"),
        plot.title = element_text(hjust = 0.5)
  )
#######
# table 3 recurent >=5  mutation coverage 
#####


#####################################
# figure 5  candidate mutation coverage
#####################################
figure5_candidate_mutation_coverage_boxplot <-
  read_rds(file.path(recheck_dir, '01.depth_for_tumor_normal_barplot.rds')) +
  labs(title = "Candidate Somatic Mutation Coverage") +
  theme(strip.background = element_rect(fill = 'white', color = 'black'))
print(figure5_candidate_mutation_coverage_boxplot)
ggsave(file.path(out_pic_dir, "figure5_candidate_mutation_coverage_boxplot.png"),
       plot = figure5_candidate_mutation_coverage_boxplot, device = 'png', width = 9, height = 12)


###############
# table 4 mutation target gene
###############
table4_target_genes <- 
  read_rds(file.path(target_gene_dir, 'nearestTargetProteinCodingGenes.rds'))
table4_target_genes %>%
  write_tsv(file.path(out_table_dir, 'table4_target_genes.tsv'))
table4_target_genes %>%
  write_rds(file.path(out_table_dir, 'table4_target_genes.rds'))
################
# figure 6 mutation target gene expression
###############



###############
# figure 7 survival analysis of two target genes
###############


# Save sessions.
save(list = ls(), file = file.path(out_table_dir, "summary_function.rda"))
