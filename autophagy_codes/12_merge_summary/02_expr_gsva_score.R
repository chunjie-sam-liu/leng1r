
# library -----------------------------------------------------------------
library(magrittr)
library(GSVA)
library(ggplot2)

# path --------------------------------------------------------------------
tcga_path = "/home/cliu18/liucj/projects/6.autophagy/TCGA"
expr_path <- "/home/cliu18/liucj/projects/6.autophagy/02_autophagy_expr/"
expr_path_a <- file.path(expr_path, "03_a_gene_expr")

# load data ---------------------------------------------------------------
expr <- readr::read_rds(file.path(tcga_path, "pancan33_expr.rds.gz"))
gene_list <- readr::read_rds(file.path(expr_path_a, "rds_03_a_atg_lys_gene_list.rds.gz"))

gene_sets <- list(
  LC3 = dplyr::filter(gene_list, process == "LC3"),
  ATG5_ATG12_ATG16 = dplyr::filter(gene_list, process == "ATG5-ATG12-ATG16"),
  ULK_Complex = dplyr::filter(gene_list, process == "ULK-complex"),
  ATG9_cycling_system = dplyr::filter(gene_list, process == "ATG9-cycling-system"),
  PI3K_III_complex = dplyr::filter(gene_list, process == "ATG9-PI3K-III-complex-system"),
  Negative = dplyr::filter(gene_list, process == "Negative"),
  Positive = dplyr::filter(gene_list, process == "Positive"),
  TRAPP = dplyr::filter(gene_list, process == "TRAPP"),
  Membrane_delivery = dplyr::filter(gene_list, process == "Membrane-delivery"),
  Mitophagy = dplyr::filter(gene_list, process == "Mitophagy"),
  Lys_comp = dplyr::filter(gene_list, process == "Lys_comp"),
  Lys_deg = dplyr::filter(gene_list, process == "Lys_deg"),
  atg = dplyr::filter(gene_list, type == "atg"),
  lys = dplyr::filter(gene_list, type == "lys"),
  atg_core = dplyr::filter(gene_list, pathway == "autophagesome formation-core")
) %>% purrr::map("symbol")

# gsva score --------------------------------------------------------------
fn_gsva <- function(.x, .y, gene_sets = gene_sets){
  # .x <- .te$cancer_types
  # .y <- .te$filter_expr[[1]]
  print(.x)
  .y %>% 
    tidyr::drop_na() %>% 
    dplyr::select( -entrez_id) -> .d
  
  .d_mat <- as.matrix(.d[,-1])
  rownames(.d_mat) <- .d$symbol
  
  .es_dif <- gsva(.d_mat, gene_sets, method = "gsva", mx.diff = TRUE, verbose = FALSE, parallel.sz = 1)
  
  .es_dif$es.obs %>% 
    as.data.frame() %>% 
    tibble::as_tibble() %>% 
    tibble::add_column(set = rownames(.es_dif$es.obs), .before = 1) -> .d_es
}

cluster <- multidplyr::create_cluster(nrow(expr))
expr %>% 
  multidplyr::partition(cluster = cluster) %>%
  multidplyr::cluster_library("magrittr") %>%
  multidplyr::cluster_library("GSVA") %>%
  multidplyr::cluster_assign_value("fn_gsva", fn_gsva)  %>%
  multidplyr::cluster_assign_value("gene_sets", gene_sets)  %>%
  dplyr::mutate(gsva = purrr::map2(.x = cancer_types, .y = expr, .f = fn_gsva, gene_sets = gene_sets)) %>% 
  dplyr::select(-expr) %>% 
  dplyr::collect() %>%
  dplyr::as_tibble() %>%
  dplyr::ungroup() %>%
  dplyr::select(-PARTITION_ID) -> gene_list_gsva
parallel::stopCluster(cluster)

readr::write_rds(x = gene_list_gsva, path = file.path(expr_path_a, ".rds_03_i_gsva_gene_list_gsva.rds.gz"), compress = "gz")

# tumor vs normal ---------------------------------------------------------
fn_tn <- function(.s){
  .s %>% stringr::str_sub(start = 14, end = 14)
}
fn_tn_pair <- function(.d){
  .d %>% 
    dplyr::mutate(
      sample = stringr::str_sub(
        string = barcode,
        start = 1,
        end = 12
      )
    ) %>%
    dplyr::mutate(type = purrr::map_chr(barcode, .f = fn_tn)) %>% 
    dplyr::filter(type %in% c("0", "1")) %>%
    dplyr::mutate(type = plyr::revalue(
      x = type,
      replace = c("0" = "Tumor", "1" = "Normal"),
      warn_missing = F
    )) 
  # %>%
  #   dplyr::group_by(sample) %>%
  #   dplyr::filter(n() >= 2, length(unique(type)) == 2) %>%
  #   dplyr::ungroup()
}
fn_t_test <- function(.d){
  .d %>% 
    dplyr::do(
      broom::tidy(
        tryCatch(
          expr = t.test(score ~ type, data = .),
          error = function(e){1},
          warning = function(e){1}
        ))) %>% 
    dplyr::mutate(diff = ifelse(tibble::has_name(., "estimate"), estimate, 0)) %>% 
    dplyr::mutate(pval = ifelse(tibble::has_name(., "p.value"), p.value, 1)) -> .dd
  
  if (nrow(.dd) > 1) dplyr::select(.dd, 1, diff, pval) else dplyr::select(.dd, diff, pval)
  
}
fn_gsva_tn <- function(.x, .y){
  # .x <- .te$cancer_types
  # .y <- .te$gsva[[1]]
  print(.x)
  
  .y %>% 
    tidyr::gather(key = "barcode", value = "score", -set) %>% 
    fn_tn_pair() -> .d
  
  sample_type_summary <- table(.d$type) %>% as.numeric()
  if (gtools::invalid(sample_type_summary) || any(sample_type_summary < c(10, 10))) return(NULL)
  
  .d %>% dplyr::group_by(set) %>% fn_t_test() %>% dplyr::ungroup()
}

gene_list_gsva %>% 
  dplyr::mutate(tn_diff = purrr::map2(.x = cancer_types, .y = gsva, .f = fn_gsva_tn)) %>% 
  dplyr::filter(!purrr::map_lgl(.x = tn_diff, is.null)) %>% 
  tidyr::unnest(tn_diff) -> gsva_tn_diff

gsva_tn_diff %>% 
  dplyr::mutate(pval = -log10(pval), diff = -diff) %>% 
  dplyr::filter(pval > -log10(0.05)) -> plot_ready

plot_ready %>% 
  dplyr::group_by(cancer_types) %>% 
  dplyr::summarise(m_diff = mean(diff)) %>% 
  dplyr::arrange(dplyr::desc(m_diff)) %>% 
  dplyr::pull(cancer_types) -> cancer_rank

plot_ready %>% 
  dplyr::group_by(set) %>% 
  dplyr::summarise(m_diff = mean(diff)) %>% 
  dplyr::arrange(dplyr::desc(m_diff)) %>% 
  dplyr::pull(set) -> path_rank

CPCOLS <- c("#000080", "#F7F7F7", "#CD0000", "#F8F8FF")
plot_ready %>% 
  dplyr::mutate(star = dplyr::case_when(
    dplyr::between(pval, -log10(0.05), 5) ~ "*",
    dplyr::between(pval, 5, 10) ~ "**",
    pval > 10 ~ "***",
    TRUE ~ "."
  )) %>% 
  ggplot(aes(x = cancer_types, y = set)) +
  geom_tile(aes(fill = diff)) +
  geom_text(aes(label = star), color = CPCOLS[4]) +
  scale_fill_gradient2(low = CPCOLS[1], mid = CPCOLS[2], high = CPCOLS[3], name = "Diff") +
  scale_size(name = "P-value") + 
  scale_x_discrete(limit = cancer_rank) +
  scale_y_discrete(limit = path_rank) +
  theme_bw() +
  labs(x = "", y = "Gene set", title = "ATG related gsva Norml vs. Tumor") +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5) 
  ) -> p
ggsave(filename = file.path(expr_path_a, "fig_summary_gsva_tumor_vs_normal.pdf"), plot = p)




