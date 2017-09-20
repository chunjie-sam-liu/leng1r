
# library -----------------------------------------------------------------

library(magrittr)
library(ggplot2)

# load path ---------------------------------------------------------------

tcga_path <- "/home/cliu18/liucj/projects/6.autophagy/TCGA"
afhl_path <- "/home/cliu18/liucj/projects/6.autophagy/09_afh_vs_afl"
afhl_class <- file.path(afhl_path, "01_af_h_l_classification")
expr_path <- "/home/cliu18/liucj/projects/6.autophagy/02_autophagy_expr/03_a_gene_expr"

# load data ---------------------------------------------------------------

rppa_expr <- readr::read_rds(file.path(tcga_path, "pancan32_rppa_expr.rds.gz"))
rppa_name <- readr::read_rds(file.path(tcga_path, "rppa_name_symbol.rds.gz"))
gene_list <- readr::read_rds(file.path(expr_path, "rds_03_a_atg_lys_gene_list.rds.gz"))
rppa_name %>% dplyr::semi_join(gene_list, by = "symbol") -> atg_rppa

PI3K_AKT <- "pS473|pT308|pS9|pT246|pT1462"


atg_rppa %>% 
  dplyr::inner_join(gene_list, by = "symbol") %>% 
  dplyr::select(1,2, process) %>% 
  dplyr::arrange(process) -> sym_func
knitr::kable(sym_func)

# atg involved rppa ---------------------------------------------------------
fn_select_marker <- function(.x, sym_func){
  # rppa_expr$protein_expr[[1]] -> .x
  sym_func %>% dplyr::inner_join(.x, by = c("protein", "symbol"))
}

rppa_expr %>% 
  dplyr::mutate(protein_expr = purrr::map(.x = protein_expr, .f = fn_select_marker, sym_func)) -> atg_rppa_expr


# p62 correlation with known proteins -------------------------------------

fn_corr_protein <- function(.x){
  .x %>% dplyr::select(-process) %>% dplyr::select(protein, symbol) -> .names
  setNames(.names$symbol, .names$protein) -> .names
  
  .x %>% 
    dplyr::select(-process, -symbol) %>% 
    tidyr::gather(key = barcode, value = rppa, -protein) %>% 
    tidyr::spread(key = protein, value = rppa) %>% 
    dplyr::select(-barcode) %>% 
    cor() -> .corr_matrix
  
  ggcorrplot::ggcorrplot(.corr_matrix, hc.order = TRUE, outline.col = "white")
}
atg_rppa_expr 



fn_select_gene <- function(.x, .gene){
  .x %>% 
    dplyr::filter(symbol == .gene) %>% 
    tidyr::gather(barcode, rppa, -c(protein, symbol, process)) -> .d
}

atg_rppa_expr
