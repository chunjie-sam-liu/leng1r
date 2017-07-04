library(magrittr)
library(ggplot2)

expr_path <- "/home/cliu18/liucj/projects/6.autophagy/02_autophagy_expr/"
rppa_path <- file.path(expr_path, "03_e_rppa")
tcga_path <- "/home/cliu18/liucj/projects/6.autophagy/TCGA"

# rppa <- readr::read_rds(path = file.path(tcga_path,"pancan_clinical_stage.rds.gz")) %>% 
  # dplyr::filter(n >= 40) %>% 
  # dplyr::select(-n)

gene_list <- readr::read_rds(file.path(expr_path, "rds_03_at_ly_comb_gene_list.rds.gz"))
gene_list_expr <- readr::read_rds(path = file.path(expr_path, ".rds_03_a_gene_list_expr.rds.gz"))
