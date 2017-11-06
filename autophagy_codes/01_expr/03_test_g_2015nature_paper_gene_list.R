
# libarry -----------------------------------------------------------------

library(magrittr)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)


# path --------------------------------------------------------------------

root_path <- "/home/liucj/project/projects/6.autophagy/"
tcga_path <- file.path(root_path, "TCGA")
expr_path <- file.path(root_path, "02_autophagy_expr")
nature_path <- file.path(expr_path, "03_test_g_2015nature_paper_gene_list")
script_path <- "/home/liucj/github/R_Leng_1/autophagy_codes/01_expr"

nature_gene_list <- readxl::read_xls(path = file.path(nature_path, 'nature14587-s1.xls'), skip = 3, sheet = 1) %>% 
  dplyr::select(symbol = `Gene symbol`, desc = `Description (disease relevance)`)

tcga_gene_symbol <- readr::read_rds(file.path(tcga_path, 'tcga_gene_symbol.rds.gz'))

expr <- readr::read_rds(path = file.path(tcga_path, 'pancan33_expr.rds.gz'))

# clustering by comparing the normal and tumor ----------------------------
nature_gene_list %>% 
  dplyr::semi_join(tcga_gene_symbol, by = "symbol") %>% 
  dplyr::pull(symbol) -> gene_list

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
  ))  %>% 
  dplyr::filter(purrr::map_lgl(expr, .f = function(.x){nrow(.x) != 0})) %>% 
  dplyr::mutate(plot = purrr::map2(
    .x = cancer_types,
    .y = expr,
    .f = function(.x, .y){
      print(.x)
        
      .y %>% 
        dplyr::arrange(type) %>% 
        dplyr::distinct(barcode, type) %>% 
        dplyr::pull(barcode) -> .bar
      
      .y %>% 
        dplyr::arrange(type) %>% 
        dplyr::distinct(barcode, type) %>% 
        dplyr::mutate(type = dplyr::recode(type, "Tumor" = 0, "Normal" = 1)) %>% 
        dplyr::pull(type) -> .sub
      if (length(unique(.sub)) == 1) return(NULL)
      .y %>% 
        # dplyr::mutate(expr = log2(expr + 0.01)) %>% 
        dplyr::select(-c(sample, type)) %>% 
        tidyr::spread(key = barcode, value = expr) %>% 
        as.data.frame() %>% 
        tibble::remove_rownames() %>% 
        tibble::column_to_rownames(var = 'symbol') %>% 
        dplyr::select(.bar) %>% 
        data.matrix() %>% 
        apply(1, scale) %>% 
        t() -> .mat
      
      cluster_col <- c("#00008B", "#E31A1C")
      names(cluster_col) <- unique(.sub)
      
      ha = HeatmapAnnotation(
        df = data.frame(cluster = .sub), 
        gap = unit(c(4,2), "mm"),
        col = list(cluster = cluster_col),
        show_legend = F,
        height = unit(0.03, "npc")
        )
      
      Heatmap(
        .mat,
        col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red"), space = "RGB"),
        name = "expression",
        
        # annotation
        top_annotation = ha, 

        # show row and columns
        show_row_names = T, 
        show_column_names = FALSE, 
        show_row_dend = F,
        show_column_dend = F, 
        
        row_names_gp = gpar(fontsize = 6),
        # clustering
        cluster_columns = F,
        clustering_distance_rows = "pearson",
        clustering_method_rows = "ward.D",
        
        show_heatmap_legend = F
      ) -> hm
    }
  )) %>% 
  dplyr::select(-expr) -> gene_heatmap

expr %>% 
  dplyr::mutate(expr = purrr::map(
    .x = expr,
    .f = function(.x) {
      .x %>% 
        dplyr::filter(symbol %in% gene_list)
    }
  )) -> gene_list_expr

# Compare autophagy lysosome signature using GSEA
gene_sets <- list(atg_lys_sig = gene_list)
gsea_path <- file.path(nature_path, "gsea")
if (!dir.exists(gsea_path)) dir.create(gsea_path)

fn_gmt <- function(.x, .path = gsea_path){
  .file_name <- file.path(.path, "atg_lys_pathway.gmt")
  .x %>% purrr::map_int(length) %>% unlist() %>% max() -> .max_length
  
  .x %>% 
    purrr::map(
      .f = function(x, .m){
          .y <- c(x, rep(NA, .m - length(x)))
          names(.y) <- 1:.m
          .y
        }, 
      .m = .max_length) %>% 
    purrr::map(.f = tibble::enframe) %>% 
    tibble::enframe() %>% 
    tidyr::unnest() %>% 
    tidyr::spread(key = name1, value = value) %>% 
    dplyr::select(1, as.character(1:.max_length)) %>% 
    tibble::add_column(des = "cj_atg_lys", .before = 2) -> .x_write
  
  .x_write %>% readr::write_delim(path = .file_name, delim = "\t", col_names = F, na = "")
}
gene_sets %>% fn_gmt(.path = gsea_path)
fn_gct_cls <- function(.x, .y, .path = gsea_path){
  # .x <- te$cancer_types
  # .y <- te$filter_expr[[1]]
  .y %>% tidyr::drop_na() -> .y
  
  tibble::tibble(
    barcode = names(.y)[-c(1,2)], 
    type = stringr::str_sub(names(.y)[-c(1,2)], start = 14, 15)) %>% 
    dplyr::filter(type %in% c("01", "11")) %>% 
    dplyr::arrange(type) %>% 
    dplyr::pull(barcode) -> .names
  
  .y %>% 
    dplyr::select(1,2, .names) %>% 
    dplyr::filter(symbol != "?") %>% 
    dplyr::distinct(symbol, .keep_all = T) %>% 
    dplyr::rename(NAME = symbol, DESCRIPTION = entrez_id) -> .y_gct
  
  .gct_filename <- file.path(.path, paste(.x,'mRNA_expression.gct', sep = "_"))
  
  # write gct files
  readr::write_lines(x = "#1.2", path = .gct_filename)
  readr::write_tsv(x = data.frame(gene = nrow(.y), sample = length(.names)), path = .gct_filename, append = T, col_names = F)
  readr::write_tsv(x = .y_gct, path = .gct_filename, append = T, col_names = T)
  
  .cls_filename <- file.path(.path, paste(.x, "mRNA_expression.cls", sep = "_"))
  readr::write_delim(x = data.frame(sample = length(.names), 2, 1), path = .cls_filename, col_names = F,  delim  = " ")
  readr::write_delim(x = data.frame("#", "Tumor", "Normal"), path = .cls_filename, col_names = F, append = T, delim = " ")
  .names %>% 
    stringr::str_sub(start = 14, 15) %>% 
    dplyr::recode("01" = "T", "11" = "N") %>% 
     as.matrix() %>% t() %>% as.data.frame() -> .label
  readr::write_delim(x = .label, path = .cls_filename, delim = " ", col_names = F, append = T)
}

# make gct cls files
cluster <- multidplyr::create_cluster(length(expr$cancer_types))
expr %>% 
  multidplyr::partition(cluster = cluster) %>%
  multidplyr::cluster_library("magrittr") %>%
  multidplyr::cluster_assign_value("fn_gct_cls", fn_gct_cls)  %>%
  multidplyr::cluster_assign_value("gsea_path", gsea_path)  %>% 
  dplyr::mutate(out = purrr::walk2(
    .x = cancer_types,
    .y = expr,
    .f = fn_gct_cls,
    .path = gsea_path
  ))
parallel::stopCluster(cluster)


# run gsea ----------------------------------------------------------------

fn_atleast2_normal <- function(cancer_types, expr){
  filter_expr <- expr
  names(filter_expr)[-c(1,2)] -> .barcode
  
  sum(.barcode %>% stringr::str_sub(14,15) %>% stringr::str_detect("11")) -> .n_normal
  
  if (.n_normal >= 10) {
    return(cancer_types)
  } else{
    return(NULL)
  }
}
cancers <- expr %>% purrr::pmap(.f = fn_atleast2_normal) %>% unlist()

fn_gsea <- function(.ds, .cls, .db, .output, .doc){
  # .ds is the gct format file
  # .cls is the phenotype format
  # .db is the gmt file format
  
  GSEA( # Input/Output Files :-------------------------------------------
        input.ds =  .ds,               # Input gene expression Affy dataset file in RES or GCT format
        input.cls = .cls,               # Input class vector (phenotype) file in CLS format
        gs.db = .db,           # Gene set database in GMT format
        output.directory = .output,            # Directory where to store output and results (default: "")
        #  Program parameters :----------------------------------------------------------------------------------------------------------------------------
        doc.string            = .doc,     # Documentation string used as a prefix to name result files (default: "GSEA.analysis")
        non.interactive.run   = F,               # Run in interactive (i.e. R GUI) or batch (R command line) mode (default: F)
        reshuffling.type      = "gene.labels", # Type of permutation reshuffling: "sample.labels" or "gene.labels" (default: "sample.labels"
        nperm                 = 1000,            # Number of random permutations (default: 1000)
        weighted.score.type   =  1,              # Enrichment correlation-based weighting: 0=no weight (KS), 1= weigthed, 2 = over-weigthed (default: 1)
        nom.p.val.threshold   = -1,              # Significance threshold for nominal p-vals for gene sets (default: -1, no thres)
        fwer.p.val.threshold  = -1,              # Significance threshold for FWER p-vals for gene sets (default: -1, no thres)
        fdr.q.val.threshold   = 0.25,            # Significance threshold for FDR q-vals for gene sets (default: 0.25)
        topgs                 = 50,              # Besides those passing test, number of top scoring gene sets used for detailed reports (default: 10)
        adjust.FDR.q.val      = F,               # Adjust the FDR q-vals (default: F)
        gs.size.threshold.min = 10,              # Minimum size (in genes) for database gene sets to be considered (default: 25)
        gs.size.threshold.max = 1000,             # Maximum size (in genes) for database gene sets to be considered (default: 500)
        reverse.sign          = F,               # Reverse direction of gene list (pos. enrichment becomes negative, etc.) (default: F)
        preproc.type          = 0,               # Preproc.normalization: 0=none, 1=col(z-score)., 2=col(rank) and row(z-score)., 3=col(rank). (def: 0)
        random.seed           = 111,             # Random number generator seed. (default: 123456)
        perm.type             = 0,               # For experts only. Permutation type: 0 = unbalanced, 1 = balanced (default: 0)
        fraction              = 1.0,             # For experts only. Subsampling fraction. Set to 1.0 (no resampling) (default: 1.0)
        replace               = F,               # For experts only, Resampling mode (replacement or not replacement) (default: F)
        save.intermediate.results = F,           # For experts only, save intermediate results (e.g. matrix of random perm. scores) (default: F)
        OLD.GSEA              = F,               # Use original (old) version of GSEA (default: F)
        use.fast.enrichment.routine = T          # Use faster routine to compute enrichment for random permutations (default: T)
  )
}
fn_run_gsea <- function(.x, .path = gsea_path, script_path = script_path){
  .gct <- file.path(.path, paste(.x, "mRNA_expression.gct", sep = "_"))
  .cls <- file.path(.path, paste(.x, "mRNA_expression.cls", sep = "_"))
  .gmt <- file.path(.path, "atg_lys_pathway.gmt")
  .output_dir <- file.path(.path, "gsea_result/")
  .doc <- paste(.x, "GSEA.analysis", sep = ".")
  if (!file.exists(.output_dir)){ dir.create(.output_dir)}
  
  source(file.path(script_path, "GSEA.1.0.r"))
  fn_gsea(.ds = .gct, .cls = .cls, .db = .gmt, .output = .output_dir, .doc = .doc)
}

tibble::tibble(cancer_types = cancers) %>% 
  head(1) %>% 
  dplyr::mutate(res = purrr::walk(.x = cancer_types, .f = fn_run_gsea, .path = gsea_path, script_path = script_path))
fn_run_gsea("THCA", .path = gsea_path, script_path = script_path)

cluster <- multidplyr::create_cluster(length(cancers))
tibble::tibble(cancer_types = cancers) %>% 
  multidplyr::partition(cluster = cluster) %>%
  multidplyr::cluster_library("magrittr") %>%
  multidplyr::cluster_assign_value("fn_run_gsea", fn_run_gsea)  %>%
  multidplyr::cluster_assign_value("fn_gsea", fn_gsea)  %>%
  multidplyr::cluster_assign_value("gsea_path", gsea_path)  %>%
  multidplyr::cluster_assign_value("script_path", script_path)  %>%
  dplyr::mutate(res = purrr::walk(.x = cancer_types, .f = fn_run_gsea, .path = gsea_path, script_path = script_path)) %>% 
  dplyr::collect()
parallel::stopCluster(cluster)
# save image --------------------------------------------------------------

save.image(file = file.path(nature_path, '.rda_03_test_g.rda'))
# load(file = file.path(nature_path, ".rda_03_test_g.rda"))


