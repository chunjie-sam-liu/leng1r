
# library -----------------------------------------------------------------
library(magrittr)
library(ggplot2)
library(GSVA)
library(msigdf)
library(biomaRt)

# path --------------------------------------------------------------------
tcga_path = "/home/cliu18/liucj/projects/6.autophagy/TCGA"
expr_path <- "/home/cliu18/liucj/projects/6.autophagy/02_autophagy_expr/"
expr_path_a <- file.path(expr_path, "03_a_gene_expr")

# load data ---------------------------------------------------------------
expr <- readr::read_rds(file.path(tcga_path, "pancan33_expr.rds.gz"))
gene_list_gsva <- readr::read_rds(path = file.path(expr_path_a, ".rds_03_i_gsva_gene_list_gsva.rds.gz"))

# MsigDB entrez to symbol -------------------------------------------------
# name c1-c7
msigdf.human$collection %>% unique() -> collection
name <- c("C1_POSITIONAL","C2_CURATED","C3_MOTIF","C4_COMPUTATIONAL","C5_GENE_ONTOLOGY","C6_ONCOGENIC_SIGNATURES","C7_IMMUNOLOGIC_SIGNATURES","HALLMARK")
setNames(name, collection) -> name

# entrez to symbol
GENES = useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host = "useast.ensembl.org")
GENES.ATTRIBUTES <- listAttributes(GENES)
GENES.FILTERS <- listFilters(GENES)
msigdf.human$entrez %>% unique() -> entrez
getBM(
  attributes = c("entrezgene", "hgnc_symbol", "ensembl_gene_id"), 
  filters = "entrezgene", values = entrez, mart = GENES) %>% 
  tidyr::drop_na() %>% 
  tibble::as_tibble() %>%
  dplyr::mutate(entrez = as.integer(entrezgene)) %>% 
  dplyr::select(symbol = hgnc_symbol, entrez) %>% 
  dplyr::distinct() %>% 
  dplyr::filter(symbol != "") -> entrez_symbol 
msigdf.human %>% 
  dplyr::mutate(name = plyr::revalue(collection, name)) %>% 
  dplyr::inner_join(entrez_symbol, by = "entrez") -> human
readr::write_rds(human, path = file.path(tcga_path, "MSigDB_human_collections.rds.gz"), compress = "gz")

C2_CURATED <- dplyr::filter(human, collection == "c2")
C5_GENE_ONTOLOGY <- dplyr::filter(human, collection == "c5")
C6_ONCOGENIC_SIGNATURES <- dplyr::filter(human, collection == "c6")
C7_IMMUNOLOGIC_SIGNATURES <- dplyr::filter(human, collection == "c7")
HALLMARK <- dplyr::filter(human, collection == "hallmark")

# MSigDB set gsva ---------------------------------------------------------
fn_geneset <- function(.g){
  .g %>% 
    dplyr::select(geneset, entrez) %>% 
    tidyr::nest(entrez) %>% 
    tibble::deframe() %>% purrr::map("entrez")
}
human %>% dplyr::filter(collection %in% c("c2", "c5","c6", "c7", "hallmark")) -> human_sets
human_sets %>% fn_geneset() -> gene_sets
human_sets %>% dplyr::select(collection = name, geneset) -> cg
fn_gsva <- function(.x, .y, gene_sets, cg){
  # .x <- .te$cancer_types
  # .y <- .te$expr[[1]]
  print(.x)
  
  .y %>% 
    tidyr::drop_na() %>% 
    dplyr::select(-symbol) %>% 
    dplyr::rename(entrez = entrez_id) %>% 
    dplyr::distinct(entrez, .keep_all = T) -> .d
  
  .d_mat <- as.matrix(.d[,-1])
  rownames(.d_mat) <- .d$entrez
  
  .es_dif <- gsva(.d_mat, gene_sets, method = "gsva", mx.diff = TRUE, verbose = FALSE, parallel.sz = 1)
  
  .es_dif$es.obs %>% 
    as.data.frame() %>% 
    tibble::as_tibble() %>% 
    tibble::add_column(geneset = rownames(.es_dif$es.obs), .before = 1) -> .d_es
  
  cg %>% 
    dplyr::inner_join(.d_es, by = "geneset")
}

# multidplyr
cluster <- multidplyr::create_cluster(nrow(expr))
expr %>% 
  multidplyr::partition(cluster = cluster) %>%
  multidplyr::cluster_library("magrittr") %>%
  multidplyr::cluster_library("GSVA") %>%
  multidplyr::cluster_assign_value("fn_gsva", fn_gsva)  %>%
  multidplyr::cluster_assign_value("gene_sets", gene_sets) %>% 
  multidplyr::cluster_assign_value("cg", cg) %>% 
  dplyr::mutate(gsva = purrr::map2(.x = cancer_types, .y = expr, .f = fn_gsva, gene_sets = gene_sets, cg = cg)) %>%
  dplyr::select(-expr) %>% 
  dplyr::collect() %>%
  dplyr::ungroup() %>%
  dplyr::select(-PARTITION_ID) -> MSigDB_GSVA_score
parallel::stopCluster(cluster)
readr::write_rds(x = MSigDB_GSVA_score, path = file.path(tcga, "pancan33_MSigDB_GSVA_score.rds.gz"), compress = "gz")


# MSigDB set gsva correlate with autophagy --------------------------------

gene_list_gsva <- readr::read_rds(path = file.path(expr_path_a, ".rds_03_i_gsva_gene_list_gsva.rds.gz"))
