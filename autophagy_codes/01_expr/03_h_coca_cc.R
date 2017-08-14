library(methods)
library(magrittr)

tcga_path = "/home/cliu18/liucj/projects/6.autophagy/TCGA"
expr_path <- "/home/cliu18/liucj/projects/6.autophagy/02_autophagy_expr/"
expr_path_a <- file.path(expr_path, "03_a_gene_expr")


load(file = file.path(expr_path_a, ".rda_03_h_coca.rda"))

# expr
expr_matrix %>% t() -> expr_matrix_t
factoextra::hcut(expr_matrix_t, k = 3, hc_func = 'hclust', hc_method = 'ward.D', hc_metric = 'pearson', stand = T) -> expr_hcut
cutree(expr_hcut, k = 3) -> expr_group

fn_encode <- function(.x){
  .d <- tibble::tibble()
  if(.x == 1) {.d <- tibble::tibble(a = 1,b = 0, c = 0)}
  if(.x == 2) {.d <- tibble::tibble(a = 0,b = 1, c = 0)}
  if(.x == 3) {.d <- tibble::tibble(a = 0,b = 0, c = 1)}
  .d
}

expr_group %>% 
  tibble::enframe(name = "sample") %>% 
  dplyr::mutate(encode = purrr::map(.x = value, .f = fn_encode)) %>% 
  tidyr::unnest() %>% 
  dplyr::select(-value) -> expr_encode


cnv_matrix %>% t() -> cnv_matrix_t
factoextra::hcut(cnv_matrix_t, k = 3, hc_func = 'hclust', hc_method = 'ward.D', hc_metric = 'pearson', stand = T) -> cnv_hcut
cutree(cnv_hcut, k = 3) -> cnv_group

cnv_group %>% 
  tibble::enframe(name = "sample") %>% 
  dplyr::mutate(encode = purrr::map(.x = value, .f = fn_encode)) %>% 
  tidyr::unnest() %>% 
  dplyr::select(-value) -> cnv_encode


methy_matrix %>% t() -> methy_matrix_t
factoextra::hcut(methy_matrix_t, k = 3, hc_func = 'hclust', hc_method = 'ward.D', hc_metric = 'pearson', stand = T) -> methy_hcut
cutree(methy_hcut, k = 3) -> methy_group

methy_group %>% 
  tibble::enframe(name = "sample") %>% 
  dplyr::mutate(encode = purrr::map(.x = value, .f = fn_encode)) %>% 
  tidyr::unnest() %>% 
  dplyr::select(-value) -> methy_encode

list(expr_encode, cnv_encode, methy_encode) %>% 
  purrr::reduce(.f = function(x, y){x %>% dplyr::inner_join(y, by = "sample")}, .init = tibble::tibble(sample = common_names[-1])) %>% 
  tidyr::gather(key = type, value = v, -sample) %>% 
  tidyr::spread(key = sample, value = v) %>% 
  dplyr::select(-type) -> cc_d


library(ConsensusClusterPlus)
ConsensusClusterPlus(cc_d %>% as.matrix(), maxK=20, reps=500,pItem=0.8,pFeature=1, title=title, clusterAlg="hc",distance="pearson",seed=1262118388.71279, plot=F) -> cc_res

cc_res %>% readr::write_rds(path = file.path(expr_path_a, ".rds_03_h_coca_cc_cc_res.rds.gz"), compress = "gz")

save.image(file = file.path(expr_path_a, ".rda_03_h_coca_cc.rda"))
load(file.path(expr_path_a, ".rda_03_h_coca_cc.rda"))






