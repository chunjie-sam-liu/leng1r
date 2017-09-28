
# library -----------------------------------------------------------------

library(ggplot2)
library(magrittr)


# path --------------------------------------------------------------------

tcga_path <- "/home/cliu18/liucj/projects/6.autophagy/TCGA"
afhl_path <- "/home/cliu18/liucj/projects/6.autophagy/09_afh_vs_afl"
afhl_class <- file.path(afhl_path, "01_af_h_l_classification")
props_path <- file.path(afhl_path, "04_props")


# load data ---------------------------------------------------------------
library(dummies)
source("/home/cliu18/liucj/github/RstudioWithGit/autophagy_codes/13_afh_vs_afl/yy_code/cal.r")
clinical <- readr::read_rds(file.path(tcga_path, "pancan34_clinical.rds.gz"))
p62_sample_classification <- readr::read_rds(path = file.path(afhl_class, '.rds_01_p62_sample_classification.rds.gz')) %>% 
  dplyr::select(-barcode) %>% 
  tidyr::nest(-cancer_types)
expr <- readr::read_rds(file.path(tcga_path, "pancan33_expr.rds.gz")) 
expr %>% 
  dplyr::mutate(expr = purrr::map(
    .x = expr,
    .f = function(.x) {
      .x %>% tidyr::unite("gene", symbol, entrez_id, sep = "|") %>% as.data.frame() -> mRNAseq
      
      tibble::tibble(names = colnames(mRNAseq)[-1]) %>% 
        dplyr::filter(stringr::str_sub(names, 14, 15) == "01") %>% 
        dplyr::mutate(sample = stringr::str_sub(names, 1, 12)) %>% 
        dplyr::distinct(sample, .keep_all = T) %>% 
        tibble::deframe() -> uni_names
      if (length(uni_names) == 0) return(NULL)
      
      mRNAseq <- mRNAseq[, c("gene", names(uni_names))]
      names(mRNAseq) <- c("gene", uni_names)
      mRNAseq.pri <- data.frame(t(mRNAseq[,2:ncol(mRNAseq)]))
      colnames(mRNAseq.pri) <- as.vector(mRNAseq[,1])
      mRNAseq.pri <- rm.zero.col(mRNAseq.pri)
    }
  )) -> expr_dum

clinical %>% 
  dplyr::inner_join(expr_dum, by = "cancer_types") %>% 
  dplyr::inner_join(p62_sample_classification, by = "cancer_types") %>% 
  dplyr::mutate(
    clinical = purrr::map(
      .x = clinical,
      .f = function(.x){
        # not take pathologic stage into account.
        .x %>% 
          dplyr::select(barcode, age_at_initial_pathologic_diagnosis,gender,race) %>% 
          dplyr::filter(!race %in% c("ASIAN", NA_character_)) %>% 
          dplyr::mutate(race = as.factor(race)) %>% 
          dplyr::rename(sample = barcode) 
      }
    )
  ) %>%
  dplyr::mutate(
    merge = purrr::map2(
      .x = clinical,
      .y = data,
      .f = function(.x, .y) {
        .y %>% 
          dplyr::mutate(myclusters = ifelse(rppa >= median(rppa), 0, 1)) %>% 
          dplyr::select(sample, Z = myclusters) %>% 
          dplyr::inner_join(.x, by = "sample") 
      }
    )
  ) %>% 
  dplyr::mutate(x = purrr::pmap(
    .l = list(
      .x = merge, 
      .y = expr,
      .z = cancer_types
    )
    .f = function(.x, .y){
      data <- .x %>% as.data.frame()
      mRNAseq.pri <- .y
      
      dummy.feature <- setdiff(colnames(.x),c("Z","age_at_initial_pathologic_diagnosis"))
      data.dum <- dummy.data.frame(data, names = dummy.feature)
      dummy.list <- attr(data.dum,"dummies")
      rm.col <- c()
      for (i in 1:length(dummy.list))
      {
        rm.col <- c(rm.col, dummy.list[[i]][length(dummy.list[[i]])])
      }
      data.dum <- data.dum[,-rm.col]
      data.dum$X0 <- rep(1, nrow(data.dum))
      exclude.col <- match(c("Z","X0"), colnames(data.dum))
      colnames(data.dum) <- gsub(" ", ".", colnames(data.dum))
      form <- as.formula(paste("Z~",paste(colnames(data.dum)[-exclude.col],collapse="+"),sep=""))
      
      mRNAseq.result <- weight.test(
        data.dum, form, mRNAseq.pri, is.continuous=TRUE,
        weight="MW", mirror.plot=FALSE, cancer, data.type = "mRNAseq", 
        outdir = paste(scripts.dir, "/",cancer,"_",analysis,sep=""), perm=FALSE)
    }
  ))
