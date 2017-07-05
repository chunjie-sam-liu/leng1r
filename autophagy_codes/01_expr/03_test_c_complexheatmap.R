library(magrittr)
library(ggplot2)

expr_path <- "/home/cliu18/liucj/projects/6.autophagy/02_autophagy_expr/"
heatmap_path <- file.path(expr_path, "03_test_c_complexheatmap")
expr_path <- file.path(expr_path, "03_a_gene_expr")

gene_list <- readr::read_rds(file.path(expr_path, "rds_03_at_ly_comb_gene_list.rds.gz"))
gene_list_expr <- readr::read_rds(path = file.path(expr_path, ".rds_03_a_gene_list_expr.rds.gz"))

fun_barcode <- function(.b){
  stringr::str_sub(
    string = .b,
    start = 1,
    end = 12
  )
} #get short barcode from long barcode
fun_tn_type <- function(.b){
  type <- .b %>% 
    stringr::str_split(pattern = "-", simplify = T) %>% 
    .[, 4] %>% 
    stringr::str_sub(1, 2)
}

fn_extract_pair <- function(cancer_types, filter_expr){
  print(cancer_types)
  filter_expr %>% 
    dplyr::select(-entrez_id) %>% 
    tidyr::gather(key = barcode, value = expr, -symbol) %>% 
    dplyr::mutate(type = fun_tn_type(barcode)) %>% 
    dplyr::filter(type %in% c("01", "11")) %>%
    dplyr::mutate(sample = fun_barcode(barcode)) %>% 
    dplyr::mutate(type = plyr::revalue(
      x = type,
      replace = c("01" = "Tumor", "11" = "Normal"),
      warn_missing = F
    )) %>%
    dplyr::group_by(sample) %>%
    dplyr::filter(n() >= 2, length(unique(type)) == 2) %>%
    dplyr::ungroup() %>% 
    dplyr::mutate(expr = log2(expr + 1)) -> samples
  
  if (samples$sample %>% unique() %>% length() < 10) {
    return(NULL)
  }
  
  tibble::tibble(cancer_types = cancer_types, filter_expr = list(samples))
  
  }

gene_list_expr %>% 
  purrr::pmap(.f = fn_extract_pair) %>% 
  dplyr::bind_rows() -> gene_list_expr_pair
  
