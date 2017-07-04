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

gene_list_expr_median_cluster <- readr::read_rds(path = file.path(expr_path, ".rds_01_b_rppa_median_cluster_expr.rds.gz"))
rppa_scores <- readr::read_rds(file.path(tcga_path, "pancan_rppa.rds.gz"))

fn_cluster_rppa <- function(median_cluster, rppa){
  median_cluster %>% 
    dplyr::filter(!is.na(expr)) %>% 
    dplyr::mutate(cluster = as.factor(cluster)) %>% 
    dplyr::mutate(cluster = plyr::revalue(cluster, replace = c("1" = "down", "2" = "up"))) %>% 
    dplyr::inner_join(rppa, by = "barcode") %>% 
    dplyr::distinct() -> merged_clean
}
fn_test <- function(merged_clean, cancer_types){
  print(cancer_types)
  merged_clean %>%
    dplyr::group_by(symbol, pathway) %>% 
    dplyr::do(
      broom::tidy(
        tryCatch(
          t.test(score ~ cluster, data = .),
          error = function(e){1},
          warning = function(e){1})
      )
    ) %>% 
    dplyr::select(symbol, pathway, p.value) %>%
    dplyr::group_by(symbol) %>% 
    dplyr::mutate(fdr = p.adjust(p.value, method = "fdr")) %>% 
    dplyr::mutate(bfi = p.adjust(p.value, method = "bonferroni")) %>% 
    dplyr::ungroup() -> clean_pval
  #9761, 11901
  merged_clean %>% 
    tidyr::spread(key = cluster, value = score) %>% 
    dplyr::group_by(symbol, pathway) %>% 
    dplyr::summarise(diff = mean(up, na.rm = T) - mean(down, na.rm = T)) %>% 
    dplyr::ungroup() -> clean_diff
  
  clean_pval %>% 
    dplyr::inner_join(clean_diff, by = c("symbol", "pathway")) 
}

gene_list_expr_median_cluster %>% 
  dplyr::inner_join(rppa_scores, by = "cancer_types") -> gene_rppa 
# 
# gene_rppa %>%
#   dplyr::filter(cancer_types == "GBM") %>% 
#   dplyr::mutate(merged_clean = purrr::map2(median_cluster, rppa, fn_cluster_rppa)) %>%
#   dplyr::select(cancer_types, merged_clean) %>%
#   dplyr::mutate(diff_pval = purrr::map2(merged_clean, cancer_types, fn_test))

cl <- parallel::detectCores()
cluster <- multidplyr::create_cluster(floor(cl * 5 / 6))
gene_rppa %>%
  multidplyr::partition(cluster = cluster) %>%
  multidplyr::cluster_library("magrittr") %>%
  multidplyr::cluster_assign_value("fn_cluster_rppa", fn_cluster_rppa)  %>%
  multidplyr::cluster_assign_value("fn_test", fn_test) %>% 
  dplyr::mutate(merged_clean = purrr::map2(median_cluster, rppa, fn_cluster_rppa)) %>% 
  dplyr::select(cancer_types, merged_clean) %>% 
  dplyr::mutate(diff_pval = purrr::map2(merged_clean, cancer_types, fn_test)) %>% 
  dplyr::collect() %>%
  dplyr::as_tibble() %>%
  dplyr::ungroup() %>%
  dplyr::select(-PARTITION_ID) -> gene_rppa_sig_pval
on.exit(parallel::stopCluster(cluster))

gene_rppa_sig_pval %>% 
  readr::write_rds(path = file.path(rppa_path, ".rds_03_e_rppa_gene_expr.rds.gz"), compress = "gz")

#-------------------------------------------------------------------------------------------------
pathway_replace <- c(
  "PI3KAKT"="PI3K/AKT",
  "RASMAPK"="RAS/MAPK",
  "TSCmTOR"="TSC/mTOR",
  "CellCycle"="Cell Cycle",
  "DNADamage"="DNA Damage Response"
)
gene_rppa_sig_pval %>% 
  tidyr::unnest(diff_pval) %>% 
  dplyr::filter(!is.na(p.value)) %>% 
  dplyr::mutate(pathway = plyr::revalue(pathway, pathway_replace)) %>% 
  dplyr::mutate(class = ifelse(fdr < 0.05 && diff > 0, "Activation", "None")) %>% 
  dplyr::mutate(class = ifelse(fdr < 0.05 && diff < 0, "Inhibition", class)) -> gene_rppa_sig_pval_class



save.image(file = file.path(rppa_path, ".rda_03_e_rppa_gene_expr.rda"))
load(file = file.path(rppa_path, ".rda_03_e_rppa_gene_expr.rda"))


