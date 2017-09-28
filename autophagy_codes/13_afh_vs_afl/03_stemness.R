# library -----------------------------------------------------------------

library(magrittr)
library(ggplot2)

# load path ---------------------------------------------------------------

tcga_path <- "/home/cliu18/liucj/projects/6.autophagy/TCGA"
afhl_path <- "/home/cliu18/liucj/projects/6.autophagy/09_afh_vs_afl"
afhl_class <- file.path(afhl_path, "01_af_h_l_classification")
stemness_path <- file.path(afhl_path, "03_stemness")

# stemness ----------------------------------------------------------------
csc <- readr::read_rds(file.path(tcga_path, "pancan33_csc_score.rds.gz"))
p62_sample_classification <- readr::read_rds(path = file.path(afhl_class, '.rds_01_p62_sample_classification.rds.gz')) %>% 
  dplyr::select(-barcode) %>% 
  tidyr::nest(-cancer_types)

csc %>% 
  dplyr::mutate(sci = purrr::map(
    .x = sci,
    .f = function(.x){
      .x %>% 
        dplyr::filter(stringr::str_sub(barcode, start = 14, 15) == "01") %>% 
        dplyr::mutate(sample = stringr::str_sub(barcode, start = 1, end = 12)) %>% 
        dplyr::distinct(sample, .keep_all = T) %>% 
        dplyr::select(-barcode)
    }
  )) %>% 
  dplyr::inner_join(p62_sample_classification, by = "cancer_types") %>% 
  dplyr::mutate(merge = purrr::map2(
    .x = sci,
    .y = data,
    .f = function(.x, .y) {
      .x %>% dplyr::inner_join(.y, by = "sample") 
    }
  )) %>% 
  dplyr::select(cancer_types, merge) -> p62_csc

p62_csc %>% 
  dplyr::mutate(pval = purrr::map(
    .x = merge,
    .f = function(.x) {
      .x %>%
        dplyr::mutate(type = dplyr::case_when(
          rppa > quantile(rppa, 0.6) ~ "L",
          rppa < quantile(rppa, 0.4) ~ "H",
          TRUE ~ "M"
        )) %>%
        dplyr::filter(type != "M") %>%
        t.test(score~type, data = .) %>%
        broom::tidy() %>% dplyr::select(estimate, p.value) -> .t_pval
      
      # t.test(score~type, data = .x) %>% 
        # broom::tidy() %>% dplyr::select(estimate, p.value) -> .t_pval
    }
  )) %>% 
  tidyr::unnest(pval) %>% 
  dplyr::mutate(cancer_types = reorder(cancer_types, estimate, sort)) -> p62_csc_pval 

CPCOLS <- c( "#CD0000","#000080")

p62_csc_pval$cancer_types %>% levels() -> levs

p62_csc_pval %>% 
  dplyr::mutate(mark = dplyr::case_when(
    dplyr::between(p.value, 0.001, 0.05) ~ "*",
    dplyr::between(p.value, 0.0001, 0.001) ~ "**",
    p.value < 0.0001 ~ "***",
    TRUE ~ ""
  )) %>% 
  dplyr::select(-merge, -p.value, -estimate) %>% 
  tibble::deframe() -> markers

p62_csc_pval %>%  
  tidyr::unnest() %>% 
  ggplot(aes(x = cancer_types, y = score, fill = type)) +
  geom_boxplot(outlier.size = NA) +
  scale_fill_manual(
    name = "Group", 
    values = CPCOLS,
    labels = c("High Autophagy", "Low Autophagy")) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    legend.position = "top",
    legend.direction = "horizontal",
    panel.grid = element_blank()
  ) +
  labs(x = "Cancer Types", y = "Stemness score")  +
  annotate(geom = "text", x = seq_along(levs), y = 2.3, label = markers[levs]) -> stemness_pattern
  
ggsave(filename = "fig_01_stemness_pattern.pdf", 
       plot = stemness_pattern, 
       width = 8,
       height = 6, 
       path = stemness_path)

# stemness marker genes ---------------------------------------------------
# cd44+, cd24-, cd133+, cd90+, aldh1a1+, FLOT2+, abcb5+, EpCAM+,OCT4, NANOG, and SOX2,FOXO1, FOXO3, FOXO4, and FOXO6
stemness_genes <- c("CD24", "CD44", "PROM1", "THY1", "ALDH1A1", "ABCB5", "EPCAM", "OCT4", "NANOG", "SOX2", "FOXO1", "FOXO3", "FOXO4", "FOXO6")

expr <- readr::read_rds(file.path(tcga_path, "pancan33_expr.rds.gz"))

expr %>% 
  dplyr::inner_join(p62_sample_classification, by = "cancer_types") %>%  
  dplyr::mutate(expr = purrr::map(
    .x = expr,
    .f = function(.x) {
      .x %>% 
        # dplyr::filter(symbol %in% stemness_genes) %>% 
        dplyr::select(-entrez_id) %>% 
        tidyr::drop_na() %>% 
        dplyr::filter(symbol != "?") %>% 
        tidyr::gather(key = barcode, value = expr, -symbol) %>% 
        dplyr::filter(stringr::str_sub(barcode, start = 14, end = 15) == "01") %>% 
        dplyr::mutate(sample = stringr::str_sub(barcode, start = 1, end = 12)) %>% 
        dplyr::select(-barcode) %>% 
        dplyr::distinct(symbol, sample, .keep_all = T) 
    }
  )) %>% 
  dplyr::mutate(merge = purrr::map2(
    .x = expr,
    .y = data,
    .f = function(.x, .y) {
      .y %>% 
        dplyr::select(sample, type) %>% 
        dplyr::inner_join(.x, by = "sample")
    }
  )) -> p62_expr_merge

p62_expr_merge %>% 
  dplyr::select(-expr, -data) %>% 
  dplyr::mutate(pval = purrr::map(
    .x = merge,
    .f = function(.x) {
      .x %>% 
        dplyr::mutate(expr = log2(expr + 0.1)) %>% 
        dplyr::group_by(symbol) %>% 
        dplyr::do(
          broom::tidy(
          tryCatch(
            t.test(expr ~ type, data = .) ,
            error = function(e){1}
          ))
        ) %>% 
        dplyr::select(symbol, fc = estimate, pval = p.value) %>% 
        dplyr::ungroup()
    }
  )) %>% 
  tidyr::unnest(pval) -> te

te %>% 
  dplyr::filter(pval < 0.001, abs(fc) > log2(10)) %>% 
  ggplot(aes(x = cancer_types, y = symbol)) +
  geom_point(aes(color = fc, size = pval)) +
  scale_color_gradient2(
    high = "red",
    mid = "white", 
    low = "blue"
  ) +
  theme(
    axis.text.x = element_text(angle = 90)
  )
save.image(file.path(stemness_path, ".rda_03_stemness.rda"))
