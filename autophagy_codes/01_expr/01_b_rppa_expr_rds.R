library(magrittr)
path <- "/extraspace/TCGA/TCGA_exp_DataPortal/mRNA_exp"
files.names <- list.files("/extraspace/TCGA/TCGA_protein/", pattern="*re_20160627")
tcga_path <- "/home/cliu18/liucj/projects/6.autophagy/TCGA"

cancer_types <- files.names %>% stringr::str_split(pattern = "\\_", simplify = T) %>% .[,1]
cancers_names <- tibble::tibble(names = files.names, cancer_types = cancer_types)

process_raw_data <- function(.x){
  print(.x)
  path <- "/extraspace/TCGA/TCGA_protein/"
  d <- 
    readr::read_tsv(file = file.path(path, .x), progress = F) %>% 
    dplyr::select(- dplyr::starts_with("X"))
  
  d %>% 
    tidyr::gather(key = barcode, value = rppa, -protein) %>% 
    dplyr::mutate(type = fun_tn_type(barcode)) %>% 
    dplyr::filter(type %in% c("01", "06")) %>% 
    dplyr::mutate(barcode = fun_barcode(barcode)) %>% 
    dplyr::distinct(protein, barcode, .keep_all = T) %>% 
    tidyr::replace_na(replace = list(rppa = 0)) %>% 
    dplyr::group_by(protein) %>% 
    dplyr::mutate(z_score = (rppa - median(rppa)) / sd(rppa)) %>% 
    dplyr::ungroup() %>% 
    dplyr::select(protein, barcode, rppa, z_score) -> d_z
  
  score <- 
    tibble::tibble(pat = c("Akt_", "GSK3-alpha-beta_|GSK3_pS9", "p27_", "EGFR_", "S6_", "4E-BP1_p", "4E-BP1_p", "4E-BP1_p", "Src_", "ER-alpha"), 
                   name = c("AKT", "GSK", "p27", "pEGFR", "pS6", "EBP1.65.37", "EBP1.65.70", "EBP1.37.70", "pSrc", "ERalpha"))
  
  score %>% 
    dplyr::mutate(p_cor = purrr::map2_dbl(.x = pat, .y = name, .f = fun_cor, .d = d_z)) %>% 
    dplyr::select(-pat) -> score_cor
  #----------------------------------------------------------------------------------------------------------  
  d_z %>% 
    dplyr::select(-z_score) %>% 
    dplyr::mutate(protein = gsub("-R-C|-R-V|-M-C|-M-V|-R-E|-G-C|-R-E|-M-E","",protein)) %>% 
    dplyr::group_by(protein, barcode) %>% 
    dplyr::mutate(rppa = mean(rppa)) %>% 
    dplyr::ungroup() %>% 
    dplyr::distinct() -> d_m
  
  # AKT <- purrr::pmap(.x = "Akt_", .y = "AKT", type = "p", pathway = "PI3KAKT", .f = fun_check, .d = d_m)
  # AKT <- fun_check(pat = "Akt_", name = "AKT", type = "p", pathway = "PI3KAKT")
  # GSK <- fun_check(pat = "GSK3.alpha.beta|GSK3_pS9", name = "GSK", type = "p", pathway = "PI3KAKT")
  # pEGFR <- fun_check(pat = "EGFR_", name = "pEGFR", type = "p", pathway = "PI3KAKT")
  # pS6 <- fun_check(pat = "S6_", name = "pS6", type = "p", pathway = "PI3KAKT")
  
  #pathways -------------------------------------------------------------------------------------
  # pi3kakt
  pi3kakt_pat <- c("Akt_", "GSK3.alpha.beta|GSK3_pS9", "p27_", "PRAS40_", "Tuberin_pT", "INPP4B", "PTEN")
  pi3kakt_name <- c("AKT", "GSK", "p27", "PRAS40", "Tuberin_pT", "INPP4B", "PTEN")
  pi3kakt_type <- c("p", "p", "p", "p", "p", "n", "n")
  pi3kakt <- rep("PI3KAKT", length(pi3kakt_pat))
  
  # RASMAPK
  rasmapk_pat <- c("A.Raf_|c\\.Jun_pS73|C.Raf_|JNK_|MAPK_|MEK1_|p38_p|p90RSK_p|Shc_p|YB-1_p")
  rasmapk_name <- c("cj")
  rasmapk_type <- c("p")
  rasmapk <- c("RASMAPK")
  
  # RTK
  rtk_pat <- c("HER2_|HER3_|Ret_|Shc_p|Src_p", "EGFR_")
  rtk_name <- c("cj", "EGFR_")
  rtk_type <- c("p", "p")
  rtk <- c("RTK", "RTK")
  
  # TSCmTOR
  tscmtor_pat <- c("4E.BP1_p|mTOR_|p70S6K_|Rictor_", "S6_")
  tscmtor_name <- c("cj", "pS6")
  tscmtor_type <- c("p", "p")
  tscmtor <- c("TSCmTOR", "TSCmTOR")
  
  #Apoptosis 
  apoptosis_pat <- c("Bak|Bid|Bim|Caspase.|Bax", "Bcl.2|Bcl.xL|Bad_|cIAP")
  apoptosis_name <- c("cj", "n")
  apoptosis_type <- c("p", "n")
  apoptosis <- c("Apoptosis", "Apoptosis")
  
  # Cell cycle
  cellcycle_pat <- c("CDK1|Cyclin|p27_|PCNA")
  cellcycle_name <- c("cj")
  cellcycle_type <- c("p")
  cellcycle <- c("CellCycle")
  
  #DNA damage
  dnad_pat <- c("53BP1|ATM|BRCA2|Chk1_|Chk2_|Ku80|Mre11|PARP|Rad50|Rad51|XRCC1|p53")
  dnad_name <- c("cj")
  dnad_type <- c("p")
  dnad <- c("DNADamage")
  
  # EMT 
  emt_pat <- c("Collagen|Fibronectin|N.Cadherin", "Claudin.7|E.Cadherin")
  emt_name <- c("cj", "cjl")
  emt_type <- c("p", "n")
  emt <- c("EMT", "EMT")
  
  # ER
  er_pat <- c("ER-|PR")
  er_name <- c("cj")
  er_type <- c("p")
  er <- c("Hormone ER")
  
  # AR
  ar_pat <- c("AR|INPP4B|GATA3|Bcl.2")
  ar_name <- c("cj")
  ar_type <- c("p")
  ar <- c("Hormone AR")
  
  
  path_way <-
    tibble::tibble(
      pat = c(pi3kakt_pat, rasmapk_pat, rtk_pat, tscmtor_pat, apoptosis_pat, cellcycle_pat, dnad_pat, emt_pat, er_pat, ar_pat), 
      name = c(pi3kakt_name, rasmapk_name, rtk_name, tscmtor_name, apoptosis_name, cellcycle_name, dnad_name, emt_name, er_name, ar_name),
      type = c(pi3kakt_type, rasmapk_type, rtk_type, tscmtor_type, apoptosis_type, cellcycle_type, dnad_type, emt_type, er_type, ar_type),
      pathway = c(pi3kakt, rasmapk, rtk, tscmtor, apoptosis, cellcycle, dnad, emt, er, ar),
      .d = list(d_m)
    )
  
  path_way %>% 
    purrr::pmap(fun_check) %>% 
    dplyr::bind_rows() %>% 
    tidyr::spread(key = type, value = rppa) %>% 
    # dplyr::mutate(n = ifelse(tibble::has_name(., "n"), n, 0)) %>% 
    tidyr::replace_na(replace = list(n = 0, p = 0)) %>% 
    dplyr::group_by(barcode, pathway) %>% 
    dplyr::summarise(score = sum(p) - sum(n)) %>% 
    dplyr::ungroup()
  
}

fn_rppa_gene_symbol <- function(names, cancer_types){
  url
  
  path <- "/extraspace/TCGA/TCGA_protein/"
  d <- 
    readr::read_tsv(file = file.path(path, names), progress = F) %>% 
    dplyr::select(- dplyr::starts_with("X")) %>% 
    dplyr::select(protein)
}
cancers_names %>%
  head(1) %>% 
  purrr::pmap(fn_rppa_gene_symbol)
  