library(magrittr)
tcga_path <- "/home/cliu18/liucj/projects/6.autophagy/TCGA"
files.names <- list.files(path="/extraspace/TCGA/TCGA_CNV/Gene/",pattern=".txt")

cancer_types <- files.names %>% stringr::str_split(pattern = "\\.|\\_", simplify = T) %>% .[,3]
cancers_names <- tibble::tibble(names = files.names, cancer_types = cancer_types)

process_raw_data <- function(.x){
  print(.x)
  path <- "/extraspace/TCGA/TCGA_CNV/Gene/"
  d <- 
    readr::read_tsv(file = file.path(path, .x), progress = F) %>% 
    dplyr::select(- dplyr::starts_with("X")) %>% 
    dplyr::rename(symbol = `Gene Symbol`)
}

cancers_names %>%
  head(1) %>%
  dplyr::mutate(cnv = purrr::map(names, process_raw_data)) -> te


cl <- parallel::detectCores()
cluster <- multidplyr::create_cluster(floor(cl * 5 / 6))
cancers_names %>%
  multidplyr::partition(cluster = cluster) %>%
  multidplyr::cluster_library("magrittr") %>%
  multidplyr::cluster_assign_value("process_raw_data", process_raw_data) %>% 
  dplyr::mutate(cnv = purrr::map(names, process_raw_data)) %>% 
  dplyr::collect() %>%
  dplyr::as_tibble() %>%
  dplyr::ungroup() %>%
  dplyr::select(-PARTITION_ID, -names) -> pancan_cnv
on.exit(parallel::stopCluster(cluster))
pancan_cnv %>% readr::write_rds(path = file.path(tcga_path, 'pancan_cnv.rds.gz'), compress = "gz")
