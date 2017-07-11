library(magrittr)
tcga_path = "/home/cliu18/liucj/projects/6.autophagy/TCGA"
expr_path <- "/home/cliu18/liucj/projects/6.autophagy/02_autophagy_expr/"
expr_path_a <- file.path(expr_path, "03_a_gene_expr")
cnv_path <- "/home/cliu18/liucj/projects/6.autophagy/03_cnv"

# load cnv and gene list
cnv <- readr::read_rds(file.path(tcga_path, "pancan_cnv.rds.gz"))
gene_list <- readr::read_rds(file.path(expr_path_a, "rds_03_a_atg_lys_gene_list.rds.gz"))

filter_gene_list <- function(.x, gene_list) {
  gene_list %>%
    dplyr::select(symbol) %>%
    dplyr::left_join(.x, by = "symbol")
}

cnv %>%
  dplyr::mutate(filter_cnv = purrr::map(cnv, filter_gene_list, gene_list = gene_list)) %>%
  dplyr::select(-cnv) -> gene_list_cnv

readr::write_rds(x = gene_list_cnv, path = file.path(cnv_path, ".rds_02_cnv_a_gene_list.rds.gz"), compress = "gz")

fn_get_amplitue_threshold <- function(.x){
  ifelse(abs(.x) < 0.3, 0, .x) -> .y
  tibble::tibble(a = sum(.y > 0) / length(.y), d = sum(.y < 0) / length(.y)) 
}
fn_get_ad <- function(.d){
  .d %>% 
    unlist(use.name = F) %>% 
    fn_get_amplitue_threshold()
}
fn_get_percent <- function(cancer_types, filter_cnv){
  filter_cnv %>%
    tidyr::nest(-symbol) %>% 
    dplyr::mutate(ad = purrr::map(data, .f = fn_get_ad)) %>% 
    dplyr::select(-data) %>% 
    tidyr::unnest(ad) %>% 
    tibble::add_column(cancer_types = cancer_types, .before = 1)
}

gene_list_cnv %>% head(2) %>% 
  dplyr::mutate(rs = purrr::map2(cancer_types, filter_cnv, fn_get_percent)) %>% 
  dplyr::select(-cancer_types, -filter_cnv) %>% 
  tidyr::unnest(rs)

cl <- parallel::detectCores()
cluster <- multidplyr::create_cluster(floor(cl * 5 / 6))
gene_list_cnv %>% 
  multidplyr::partition(cluster = cluster) %>%
  multidplyr::cluster_library("magrittr") %>%
  multidplyr::cluster_assign_value("fn_get_amplitue_threshold", fn_get_amplitue_threshold)  %>%
  multidplyr::cluster_assign_value("fn_get_ad", fn_get_ad) %>% 
  multidplyr::cluster_assign_value("fn_get_percent", fn_get_percent) %>% 
  dplyr::mutate(rs = purrr::map2(cancer_types, filter_cnv, fn_get_percent)) %>% 
  dplyr::collect() %>%
  dplyr::as_tibble() %>%
  dplyr::ungroup() %>%
  dplyr::select(-PARTITION_ID) %>% 
  dplyr::select(-cancer_types, -filter_cnv) %>% 
  tidyr::unnest(rs) -> gene_list_cnv_per
on.exit(parallel::stopCluster(cluster))

library(ggplot2)
gene_list_cnv_per %>% 
  tidyr::drop_na() %>% 
  tidyr::gather(key = type, value = per, a, d) %>% 
  dplyr::mutate(type = plyr::revalue(type, replace = c("a" = "Amplification", "d" = "Deletion"))) %>% 
  dplyr::mutate(per = ifelse(per > 0.7, 0.7, per)) -> plot_ready
plot_ready %>% 
  ggplot(aes(y = symbol, x = cancer_types)) +
  geom_point(aes(size = per, color = type)) +
  scale_size_continuous(
    name = "CNV perl",
    breaks = c(0.1, 0.3, 0.5, 0.7),
    limits = c(0.1, 0.8),
    labels = c("10", "30", "50", "70")
  ) +
  ggthemes::scale_color_gdocs(
    name = "SCNA Type"
  ) +
  facet_wrap(~ type) -> p
ggsave(filename = "01_SCNV.pdf", plot = p, device = "pdf", path = cnv_path, width = 25, height = 30)


# filter_pattern based on
filter_pattern <- function(type) {
  if(type == "Amplification"){return(1)} 
  else if(type == "Deletion"){return(-1)}
  else{return(0)}
}

# get up down pattern
get_pattern <- function(.x){
  .x %>% 
    # dplyr::filter(Normal > 10 & Tumor > 10) %>%
    dplyr::mutate(pattern = purrr::map_dbl(type, filter_pattern)) %>%
    dplyr::select(cancer_types, symbol, pattern) %>% 
    dplyr::distinct(cancer_types, symbol, .keep_all = T) %>% 
    tidyr::spread(key = cancer_types, value = pattern) %>% print(n = Inf)
    # dplyr::mutate_all(.funs = dplyr::funs(ifelse(is.na(.), 0, .))) 
    dplyr::mutate_if(.predicate = is.numeric, .fun = dplyr::funs(ifelse(is.na(.), 0, .)))
  
}

# gene rank by up and down
get_gene_rank <- function(pattern){
  pattern %>% 
    dplyr::rowwise() %>%
    dplyr::do(
      symbol = .$symbol,
      rank =  unlist(.[-1], use.names = F) %>% sum(),
      up = (unlist(.[-1], use.names = F) == 1) %>% sum(),
      down = (unlist(.[-1], use.names = F) == -1) %>% sum()
    ) %>%
    dplyr::ungroup() %>%
    tidyr::unnest() %>%
    dplyr::mutate(up_p = up / 14, down_p = down / 14, none = 1 - up_p - down_p) %>% 
    dplyr::arrange(rank)
}

# cancer types rank by gene
get_cancer_types_rank <- function(pattern){
  pattern %>% 
    dplyr::summarise_if(.predicate = is.numeric, dplyr::funs(sum(abs(.)))) %>%
    tidyr::gather(key = cancer_types, value = rank) %>%
    dplyr::arrange(dplyr::desc(rank))
}

plot_ready %>% get_pattern() -> pattern

ampl <- plot_ready %>% 
  dplyr::filter(type == "Amplification")



save.image(file = file.path(cnv_path, ".rda_02_cnv_a_gene_list.rda"))
load(file = file.path(cnv_path, ".rda_02_cnv_a_gene_list.rda"))




