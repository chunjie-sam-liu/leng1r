library(ggplot2)
`%>%` <- magrittr::`%>%`


# Path
tcga_path <- "/home/cliu18/liucj/projects/6.autophagy/TCGA"
mirna_path <- "/home/cliu18/liucj/projects/6.autophagy/07_mirna"
expr_path <- "/home/cliu18/liucj/projects/6.autophagy/02_autophagy_expr/"
expr_path_a <- file.path(expr_path, "03_a_gene_expr")

# load methylation and gene list
mirna_target <- readr::read_rds(file.path(tcga_path, "mirna_gene_target.rds.gz"))
gene_list <- readr::read_rds(file.path(expr_path_a, "rds_03_a_atg_lys_gene_list.rds.gz"))
mirna_expr <- readr::read_rds(file.path(tcga_path, "pancan33_mirna_expr.rds.gz"))

gene_list %>% 
  dplyr::filter(pathway == "autophagesome formation-core") %>% 
  dplyr::select(symbol) %>% 
  dplyr::left_join(mirna_target, by = "symbol") -> gene_list_mirna
gene_list_mirna %>% readr::write_tsv(path = file.path(mirna_path, "02_a_gene_list_mirna.tsv"))

fn_filter_gene_list <- function(.x, gene_list) {
  gene_list_mirna %>%
    dplyr::rename(name = mirna) %>% 
    tidyr::drop_na() %>% 
    dplyr::inner_join(.x, by = "name")
}

mirna_expr %>% 
  dplyr::mutate(mirna = purrr::map(.x = mirna, .f = fn_filter_gene_list)) -> gene_list_mirna_expr 

gene_list_mirna_expr %>% 
  tidyr::unnest()



save.image(file = file.path(mirna_path, ".rda_02_mirna_a_gene_list_target.rda"))
load(file = file.path(mirna_path, ".rda_02_mirna_a_gene_list_target.rda"))
