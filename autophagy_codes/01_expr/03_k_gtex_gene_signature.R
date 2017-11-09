
# Library -----------------------------------------------------------------
library(magrittr)
library(GSVA)
library(ggplot2)

# Path --------------------------------------------------------------------

root_path <- "/home/liucj/project/projects/6.autophagy/"
gtex_path <- file.path(root_path, "GTEx")
expr_path <- file.path(root_path, "02_autophagy_expr")
nature_path <- file.path(expr_path, "03_test_g_2015nature_paper_gene_list")
gtex_out_path <- file.path(root_path, "02_autophagy_expr/03_k_gtex_gene_signature")
tcga_path <- file.path(root_path, "TCGA")

# Load data ---------------------------------------------------------------

nature_gene_list <-
  readxl::read_xls(
    path = file.path(nature_path, 'nature14587-s1.xls'),
    skip = 3,
    sheet = 1
  ) %>%
  dplyr::select(symbol = `Gene symbol`, desc = `Description (disease relevance)`)

gtex_expr <- 
  readr::read_rds(
    path = file.path(gtex_path, "gtex_gene_tmp_annotation_phenotype_v7.rds.gz")
  )

pcc <- readr::read_tsv(file = file.path(tcga_path, "PanCanAtlasTumors_color_coded_by_organ_system_20170302.tsv"))

# Gtex gene list ----------------------------------------------------------
gtex_expr %>% 
  dplyr::mutate(
    expr = purrr::map(
      .x = expr,
      .f = function(.x) {
        .x %>% 
          dplyr::filter(symbol %in% nature_gene_list$symbol)
      }
    )
  ) -> gtex_gene_list_expr

# Calculate gsva for every sample -----------------------------------------
gene_sets <- list(atg_lys = nature_gene_list$symbol)
fn_gsva <- function(.y, gene_sets = gene_sets){
  .y %>% 
    tidyr::drop_na() %>% 
    dplyr::select( -ensembl_gene_id) -> .d
  
  .d_mat <- as.matrix(.d[,-1])
  rownames(.d_mat) <- .d$symbol
  
  .es_dif <- gsva(.d_mat, gset.idx.list = gene_sets, method = "gsva", mx.diff = TRUE, verbose = FALSE, parallel.sz = 1)
  
  .es_dif %>% 
    as.data.frame() %>% 
    tibble::as_tibble() %>% 
    tibble::add_column(set = "atg_lys", .before = 1) -> .d_es
}

gtex_expr %>% 
  dplyr::mutate(
    gsva = purrr::map(
      .x = expr,
      .f = function(.x) {
        fn_gsva(.x, gene_sets = gene_sets)
      }
    )
  ) -> gtex_expr_gsva

gtex_expr_gsva %>% 
  readr::write_rds(
    path = file.path(gtex_out_path, ".rds_01_gtex_expr_gsva.rds.gz"),
    compress = "gz"
  )


# Autophagy lysosome gene signature distribution across tissue ------------

gtex_expr_gsva %>% 
  dplyr::select(SMTS, gsva) %>% 
  dplyr::mutate(
    gsva = purrr::map(
      .x = gsva,
      .f = function(.x) {
        .x %>% 
          dplyr::select(-set) %>% 
          tidyr::gather(key = barcode, value = gsva)
      }
    )
  ) %>% 
  tidyr::unnest() -> plot_ready

# Color of tissue
gtex_expr_gsva %>% 
  dplyr::select(SMTS) %>% 
  dplyr::mutate(color = pcc$color %>% head(30)) -> tcc

plot_ready %>% 
  dplyr::group_by(SMTS) %>% 
  dplyr::summarise(m = median(gsva)) %>% 
  dplyr::arrange(m) %>% 
  dplyr::pull(SMTS) -> lev

plot_ready %>% 
  dplyr::mutate(SMTS = factor(SMTS, levels = lev)) %>% 
  ggplot(aes(x = SMTS, y = gsva)) +
  stat_boxplot(geom = 'errorbar', width = 0.3) +
  geom_boxplot(outlier.colour = NA) +
  geom_point(
    aes(color = SMTS), 
    position = position_jitter(width = 0.05), 
    alpha = 0.4,
    size = 0.8) +
  scale_color_manual(
    name = "Tissues",
    values = dplyr::slice(tcc, match(lev, SMTS)) %>% dplyr::pull(color)
  ) +
  theme(
    axis.line = element_line(color = "black"),
    axis.text.x = element_text(hjust = 1, vjust = 0.5, angle = 90),
    
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    panel.spacing.x = unit(0, "lines")
  ) +
  guides(color = F) +
  labs(
    x = "Tissue", 
    y = "Autopahgy Lysosome Signature GSVA Score", 
    title = "") -> plot_tissue_signature_distribution

ggsave(
  filename = "fig_01_tissue_signature_distribution.pdf",
  plot = plot_tissue_signature_distribution,
  device = "pdf",
  path = gtex_out_path,
  width = 8,
  height = 4
)


# Test tissue difference --------------------------------------------------

gtex_expr_gsva %>% 
  dplyr::slice(match(lev, SMTS)) %>% 
  dplyr::select(SMTS, gsva) %>% 
  dplyr::mutate(
    gsva = purrr::map(
      .x = gsva,
      .f = function(.x) {
        .x %>% 
          dplyr::select(-set) %>% 
          unlist()
      }
    )
  ) -> to_test_data

gsva_test <- function(V1, V2) {
  .slice <- c(V1, V2)
  to_test_data %>% 
    dplyr::slice(.slice) %>% 
    tibble::deframe() -> .d
  names(.d) -> .d_name
  
  t.test(.d[[.d_name[1]]], .d[[.d_name[2]]]) %>% 
    broom::tidy() %>% 
    dplyr::select(
      diff = estimate,
      rlang::UQ(.d_name[1]) := estimate1,
      rlang::UQ(.d_name[2]) := estimate2,
      pval = p.value
    )
}

embed(1:30, 2) %>% 
  as.data.frame() %>% 
  purrr::pmap(.f = gsva_test) ->
  test_result

gtools::running(1:30, fun = function(x) x, width = 2) %>% 
  t() %>% 
  as.data.frame() %>% 
  purrr::pmap(.f = gsva_test) ->
  test_result

gtools::combinations(n = 30, r = 2, v = 1:30) %>% 
  as.data.frame() %>% 
  purrr::pmap(.f = gsva_test) ->
  test_all_result

test_all_result %>% purrr::map("pval") %>% unlist() -> pvals
# From the pvalue, the signature between 

# Save immage -------------------------------------------------------------

save.image(file = file.path(gtex_out_path, ".rda_03_k_gtex.rda"))
