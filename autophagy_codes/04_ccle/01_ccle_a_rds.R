library(magrittr)
tcga_path <- "/home/cliu18/liucj/projects/6.autophagy/TCGA"
ctrp_drug <- readr::read_rds(file.path(tcga_path, "ctrp_gdsc_exp_spearman.rds.gz"))
dgsc_drug <- readr::read_rds(file.path(tcga_path, "drug_gdsc_exp_spearman.rds.gz"))


