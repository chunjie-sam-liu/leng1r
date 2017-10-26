
# libarry -----------------------------------------------------------------

library(magrittr)
library(ggplot2)
# library(googlesheet)

# path --------------------------------------------------------------------

root_path <- "/home/liucj/project/projects/6.autophagy/"
tcga_path <- file.path(root_path, "TCGA")
expr_path <- file.path(root_path, "02_autophagy_expr")
nature_path <- file.path(expr_path, "03_test_g_2015nature_paper_gene_list")

nature_gene_list <- readxl::read_xls(path = file.path(nature_path, 'nature14587-s1.xls'), skip = 3, sheet = 1) %>% 
  dplyr::select(symbol = `Gene symbol`, desc = `Description (disease relevance)`)

tcga_gene_symbol <- readr::read_rds(file.path(tcga_path, 'tcga_gene_symbol.rds.gz'))

expr <- readr::read_rds(path = file.path(tcga_path, 'pancan33_expr.rds.gz'))

nature_gene_list %>% dplyr::semi_join(tcga_gene_symbol, by = "symbol") %>% dplyr::pull(symbol) -> gene_list

expr %>%
  dplyr::mutate(expr = purrr::map(
    .x = expr,
    .f = function(.x){
      .x %>% 
        dplyr::select(-entrez_id) %>% 
        dplyr::filter(symbol %in% gene_list) %>% 
        tidyr::gather(key = barcode, value = expr, -symbol) %>% 
        dplyr::mutate(sample = stringr::str_sub(barcode, start = 1, end = 12)) %>% 
        dplyr::mutate(type = stringr::str_sub(barcode, start = 14, end = 16)) %>% 
        dplyr::filter(type %in% c("01A", "11A")) %>% 
        dplyr::mutate(type = dplyr::recode(type, "01A" = "Tumor", "11A" = "Normal"))
        
    }
  )) -> gene_list_expr

gene_list_expr %>% 
  dplyr::mutate(plot = purrr::map2(
    .x = cancer_types,
    .y = expr,
    .f = function(.x, .y){
      .y %>% 
        dplyr::mutate(expr = log2(expr + 0.01)) %>% 
        dplyr::select(-c(sample, type)) %>% 
        tidyr::spread(key = barcode, value = expr) %>% 
        as.data.frame() %>% 
        tibble::remove_rownames() %>% 
        tibble::column_to_rownames(var = 'symbol') %>% 
        data.matrix() -> .mat
    }
  ))


# clustering by comparing the normal and tumor ----------------------------


# save image --------------------------------------------------------------

save.image(file = file.path(nature_path, '.rda_03_test_g.rda'))
