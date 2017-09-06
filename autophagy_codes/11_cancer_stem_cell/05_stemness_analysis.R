library(magrittr)
library(ggplot2)
tcga_dir <- "/home/cliu18/liucj/projects/6.autophagy/TCGA"
pcbc_dir <- "/home/cliu18/liucj/projects/6.autophagy/synapse/PCBC"
csc_dir <- "/home/cliu18/liucj/projects/6.autophagy/08_cancer_stem_cell"
expr_path <- "/home/cliu18/liucj/projects/6.autophagy/02_autophagy_expr/"
expr_path_a <- file.path(expr_path, "03_a_gene_expr")

scores_gsva_sci_cor <- readr::read_rds(path = file.path(csc_dir, ".rds_04_stemness_index_scores_gsva_sci_cor.rds.gz"))

# basic function -----------------------------------
fn_tn <- function(.s){
  .s %>% stringr::str_sub(start = 14, end = 14)
}
fn_tn_pair <- function(.d){
  .d %>% 
    dplyr::mutate(
      sample = stringr::str_sub(
        string = barcode,
        start = 1,
        end = 12
      )
    ) %>%
      dplyr::mutate(type = purrr::map_chr(barcode, .f = fn_tn)) %>% 
      dplyr::filter(type %in% c("0", "1")) %>%
      dplyr::mutate(type = plyr::revalue(
        x = type,
        replace = c("0" = "Tumor", "1" = "Normal"),
        warn_missing = F
      )) 
  # %>%
  #   dplyr::group_by(sample) %>%
  #   dplyr::filter(n() >= 2, length(unique(type)) == 2) %>%
  #   dplyr::ungroup()
}
fn_t_test <- function(.d){
  .d %>% 
    dplyr::do(
      broom::tidy(
        tryCatch(
          expr = t.test(score ~ type, data = .),
          error = function(e){1},
          warning = function(e){1}
        ))) %>% 
    dplyr::mutate(diff = ifelse(tibble::has_name(., "estimate"), estimate, 0)) %>% 
    dplyr::mutate(pval = ifelse(tibble::has_name(., "p.value"), p.value, 1)) -> .dd
  
  if (nrow(.dd) > 1) dplyr::select(.dd, 1, diff, pval) else dplyr::select(.dd, diff, pval)
  
}

# sci tumor vs. normal----------------------------------
sci <- scores_gsva_sci_cor %>% dplyr::select(cancer_types, sci)
fn_sci_tn <- function(.x, .y){
  print(.x)
  # .x <- .te$cancer_types
  # .y <- .te$sci[[1]]
  
  .y %>% fn_tn_pair() -> .d
  
  sample_type_summary <- table(.d$type) %>% as.numeric()
  if (gtools::invalid(sample_type_summary) || any(sample_type_summary < c(10, 10))) return(NULL)
  
  .d %>% fn_t_test()
}

sci %>%
  dplyr::mutate(tn_diff = purrr::map2(.x = cancer_types, .y = sci, .f = fn_sci_tn)) %>%
  dplyr::filter(!purrr::map_lgl(.x = tn_diff, is.null)) %>% 
  tidyr::unnest(tn_diff) -> sci_tn_diff

CPCOLS <- c("#191970", "#EE0000")

sci_tn_diff %>% 
  dplyr::filter(pval < 0.05) %>%
  tidyr::unnest(sci) %>% 
  dplyr::mutate(or = glue::glue("{cancer_types} p-value: {signif(-log10(pval), digits = 3)}")) %>% 
  dplyr::mutate(or = reorder(or, pval)) %>% 
  fn_tn_pair() %>% 
  ggplot(aes(x = type, y = score)) +
  geom_boxplot(aes(color = type)) +
  scale_color_manual(name = "Type", values = CPCOLS) +
  facet_wrap(facets =  ~or) +
  labs(y = "Stemness score", title = "Stemness score between Tumor and Normal") +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    
    strip.background = element_rect(color = "black", fill = "white"),
    
    legend.position = "bottom",
    legend.direction = "horizontal"
  ) -> stemness_tn_p

ggsave(filename = "fig_01_stemness_tn_boxplot.pdf", plot = stemness_tn_p, device = "pdf", path = csc_dir, width = 7, height = 8)
  

# gsva tumor vs. normal---------------------------------------
gsva <- scores_gsva_sci_cor %>% dplyr::select(cancer_types, gsva)
fn_gsva_tn <- function(.x, .y){
  # .x <- .te$cancer_types
  # .y <- .te$gsva[[1]]
  print(.x)
  
  .y %>% 
    tidyr::gather(key = "barcode", value = "score", -set) %>% 
    fn_tn_pair() -> .d
  
  sample_type_summary <- table(.d$type) %>% as.numeric()
  if (gtools::invalid(sample_type_summary) || any(sample_type_summary < c(10, 10))) return(NULL)
  
  .d %>% dplyr::group_by(set) %>% fn_t_test() %>% dplyr::ungroup()
}

gsva %>% 
  dplyr::mutate(tn_diff = purrr::map2(.x = cancer_types, .y = gsva, .f = fn_gsva_tn)) %>% 
  dplyr::filter(!purrr::map_lgl(.x = tn_diff, is.null)) %>% 
  tidyr::unnest(tn_diff) -> gsva_tn_diff

gsva_tn_diff %>% 
  dplyr::mutate(pval = -log10(pval), diff = -diff) %>% 
  dplyr::filter(pval > -log10(0.05)) -> plot_ready

plot_ready %>% 
  dplyr::group_by(cancer_types) %>% 
  dplyr::summarise(m_diff = mean(diff)) %>% 
  dplyr::arrange(dplyr::desc(m_diff)) %>% 
  dplyr::pull(cancer_types) -> cancer_rank

plot_ready %>% 
  dplyr::group_by(set) %>% 
  dplyr::summarise(m_diff = mean(diff)) %>% 
  dplyr::arrange(dplyr::desc(m_diff)) %>% 
  dplyr::pull(set) -> path_rank

plot_ready %>% 
  ggplot(aes(x = cancer_types, y = set, size = pval, color = diff)) +
  geom_point() +
  scale_color_gradient2(low = "blue", mid = "white", high = "red", name = "Difference") +
  scale_size(name = "P-value") + 
  scale_x_discrete(limit = cancer_rank) +
  scale_y_discrete(limit = path_rank, labels = rev(c("Pathway", "ATG", "ATG_involved", "ATG_core", "Lys_deg", "Lys", "Lys_compo"))) +
  ggthemes::theme_gdocs() +
  labs(x = "", y = "Gene set", title = "ATG related gsva Norml vs. Tumor") +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5) 
  ) -> gsva_tn_p

ggsave(filename = "fig_02_gsva_point.pdf", plot = gsva_tn_p, device = "pdf", path = csc_dir, width = 10, height = 7)

# gsea to validate gsva result --------------------------------------------
gsea <- readr::read_rds(path = file.path(expr_path_a, ".rds_02_i_gsea_cancer_es.rds.gz"))

gsea %>% dplyr::filter(q_val < 0.05) %>% dplyr::mutate(q_val = -log10(q_val)) -> plot_ready

plot_ready %>% 
  ggplot(aes(x = cancer_types, y = GS, size = q_val, color = NES)) +
  geom_point() +
  scale_color_gradient2(low = "blue", mid = "white", high = "red") +
  scale_y_discrete(limits = c("pathway", "lys_deg", "lys_com", "lys"), labels = rev(c("Pathway", "Lys_deg", "Lys_compo", "Lys"))) +
  scale_size(name = "Q-value") + 
  ggthemes::theme_gdocs() +
  labs(x = "", y = "Gene set", title = "ATG related gsea Norml vs. Tumor") +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5) 
  ) -> gsea_tn_p

ggsave(filename = "fig_03_gsea_point.pdf", plot = gsea_tn_p, device = "pdf", path = csc_dir, width = 10, height = 7)

# gsea and gsva result are not consistent.

# gsea score and stemness scores---------------------------

pathway_diff <- scores_gsva_sci_cor %>% dplyr::select(cancer_types, pathway_diff)

pathway_diff %>% 
  tidyr::unnest() %>% 
  dplyr::mutate(cor_pval = -log10(cor_pval), t_pval = -log10(t_pval)) %>% 
  dplyr::mutate(cor_pval = ifelse(cor_pval > 30, 30, cor_pval)) %>% 
  dplyr::mutate(t_pval = ifelse(t_pval > 20, 20, t_pval)) %>% 
  dplyr::mutate(pathway = plyr::revalue(pathway, replace = path_map)) %>% 
  dplyr::mutate(repel = glue::glue("{cancer_types}:{pathway}")) -> plot_ready

plot_ready %>% 
  dplyr::group_by(cancer_types) %>% 
  dplyr::summarise(m = mean(cor_coef)) %>% 
  dplyr::arrange(m) %>% 
  dplyr::pull(cancer_types) -> cancer_rank

plot_ready %>% 
  dplyr::group_by(pathway) %>% 
  dplyr::summarise(m = mean(cor_coef)) %>% 
  dplyr::arrange(m) %>% 
  dplyr::pull(pathway) -> path_rank

c("pathway" = "Pathway", "atg_core" = "ATG_core", "atg" = "ATG", "atg_form" = "ATG_involved", "lys_com" = "Lys_com", "lys" = "Lys", "lys_deg" = "Lys_deg") -> path_map

plot_ready %>%
  dplyr::filter(cor_pval > -log10(cor_pval)) %>% 
  ggplot(aes(x = cancer_types, y = pathway)) +
  geom_point(aes(color = cor_coef, size = cor_pval)) +
  scale_color_gradient2(low = "blue", mid = "white", high = "red", name = "Coef") +
  scale_size(name = "P-value") +
  scale_x_discrete(limit = cancer_rank) +
  scale_y_discrete(limit = path_rank, labels = plyr::revalue(path_rank, path_map)) +
  labs(x = "", y = "Gene set") +
  ggthemes::theme_gdocs() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  ) -> gsva_stemness_coef

ggsave(filename = "fig_04_gsea_stemness_coef.pdf", plot = gsva_stemness_coef, device = "pdf", path = csc_dir)

plot_ready %>% 
  dplyr::mutate(sig = ifelse(cor_pval > -log10(0.05) & abs(cor_coef) > 0.3, "Sig", "Not Sig")) %>% 
  ggplot(aes(x = cor_coef, y = cor_pval)) +
  geom_point(aes(color = sig)) +
  scale_color_manual(values = c("grey", "red")) +
  ggrepel::geom_text_repel(
    data = dplyr::filter(plot_ready, abs(cor_coef) > 0.3, cor_pval > -log10(0.05), pathway != "Pathway"), 
    aes(label = repel)) +
  theme_bw() 
  
plot_ready %>%
  dplyr::filter(t_pval > -log10(t_pval)) %>% 
  ggplot(aes(x = cancer_types, y = pathway)) +
  geom_point(aes(color = t_diff, size = t_pval)) +
  scale_color_gradient2(low = "blue", mid = "white", high = "red", name = "Diff") +
  scale_x_discrete(limit = cancer_rank) +
  scale_y_discrete(limit = path_rank, labels = plyr::revalue(path_rank, path_map)) +
  scale_size(name = "P-value") +
  labs(x = "", y = "Gene set") +
  ggthemes::theme_gdocs() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  ) -> gsva_stemness_diff
ggsave(filename = "fig_05_gsea_stemness_diff.pdf", plot = gsva_stemness_diff, device = "pdf", path = csc_dir)

# stemness score for the autophagy gene mrna expression----------------
fn_stemness_high_low <- function(.x, .y){
  # .x <- .te$cancer_types
  # .y <- .te$sci[[1]]
  
  .y %>% 
    fn_tn_pair() %>% 
    dplyr::group_by(type) %>% 
    dplyr::mutate(hln = ifelse(score > 0, "HS", "LS")) %>% 
    dplyr::ungroup() %>% 
    dplyr::mutate(hln = ifelse(type == "Normal", "Normal", hln))
}
fn_merge_expr_sci <- function(.x, .y){
  .x <- .te$filter_expr[[1]]
  .y <- .te$sci[[1]]
  
  .x %>% 
    dplyr::select(-entrez_id) %>% 
    tidyr::gather(barcode, expr, -symbol) %>% 
    dplyr::filter(!is.na(expr)) %>%
    dplyr::inner_join(.y, by = "barcode") %>% 
    dplyr::mutate(expr = log2(expr + 0.01)) -> .d
  
  .d %>%
    dplyr::filter(symbol == "SNCA") %>%
    dplyr::select(expr, score) %>%
    cor(method = "spearman")
  
  .d %>% 
    dplyr::filter(hln %in% c("LS", "HS")) %>% 
    dplyr::group_by(symbol) %>% 
    dplyr::do(
      broom::tidy(
        tryCatch(
          expr = t.test(expr ~ hln, data = .),
          error = function(e){1},
          warning = function(w){1}
        )
      )
    ) %>% 
    dplyr::ungroup() %>% 
    dplyr::mutate(fdr = p.adjust(p.value, method = "fdr")) %>% 
    dplyr::select(symbol, HS_expr = estimate1, LS_expr = estimate2, p.value, fdr) %>% 
    dplyr::mutate(fold = HS_expr - LS_expr) %>% 
    dplyr::mutate(sig = ifelse(fdr < 0.05 & abs(fold) > log2(1.5), "Sig", "Non-sig")) -> .dd
  .dd
}
  
sci %>% 
  dplyr::mutate(sci = purrr::map2(.x = cancer_types, .y = sci, .f = fn_stemness_high_low)) -> sci_diff

gene_list_expr <- readr::read_rds(path = file.path(expr_path_a, ".rds_03_a_gene_list_expr.rds.gz"))

gene_list_expr %>% 
  dplyr::inner_join(sci_diff, by = "cancer_types") %>% 
  dplyr::filter(cancer_types == "BRCA") -> .te
  dplyr::mutate(hl_diff = purrr::map2(.x = filter_expr, .y = sci, .f = fn_merge_expr_sci)) ->
  gene_list_sci_diff

gene_list_sci_diff %>% 
  dplyr::select(1,4) %>% 
  tidyr::unnest() %>% 
  dplyr::filter(sig == "Sig") %>% 
  dplyr::mutate(fdr = -log10(fdr)) %>% 
  dplyr::mutate(fold = dplyr::case_when(fold > 2 ~ 2, fold < -2 ~ -2, TRUE ~ fold)) -> plot_ready

plot_ready %>% 
  ggplot(aes(x = cancer_types, y = symbol, size = fdr, color = fold)) +
  geom_point() +
  scale_color_gradient2(low = "blue", mid = "white", high = "red")

# save----------------------------------------------------------------------
save.image(file = file.path(csc_dir, ".rda_05_stemness_analysis.rda"))








