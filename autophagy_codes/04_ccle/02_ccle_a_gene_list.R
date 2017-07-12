library(magrittr)
tcga_path = "/home/cliu18/liucj/projects/6.autophagy/TCGA"
expr_path <- "/home/cliu18/liucj/projects/6.autophagy/02_autophagy_expr/"
expr_path_a <- file.path(expr_path, "03_a_gene_expr")
ccle_path <- "/home/cliu18/liucj/projects/6.autophagy/05_ccle"

# load data
drug_gdsc <- readr::read_rds(file.path(tcga_path, "drug_gdsc_exp_spearman.rds.gz"))
drug_ctrp <- readr::read_rds(file.path(tcga_path, "drug_ctrp_exp_spearman.rds.gz"))
gene_list <- readr::read_rds(file.path(expr_path_a, "rds_03_a_atg_lys_gene_list.rds.gz"))

# GDSC analysis
gene_list %>% 
  dplyr::inner_join(drug_gdsc, by = "symbol") -> gdsc_gene_list

fn_filter_drug <- function(.x, .cor = 0.3, .fdr = 0.05) {
  .x %>% dplyr::filter(abs(cor_sprm) > .cor, fdr < .fdr)
}

gdsc_gene_list %>% 
  # dplyr::filter(symbol == "ATG7") %>% 
  # .$drug %>% .[[1]] -> .x
  dplyr::mutate(cor_drug = purrr::map(.x = drug, .f = fn_filter_drug)) %>% 
  tidyr::unnest(cor_drug) -> gdsc_gene_list_sig_drug

# gene rank
library(ggplot2)
gdsc_gene_list_sig_drug %>% 
  dplyr::mutate(
    fdr = ifelse(-log10(fdr) > 40, 40, -log10(fdr)),
    cor_sprm = ifelse(cor_sprm > 0.4, 0.4, cor_sprm),
    cor_sprm = ifelse(cor_sprm < -0.4, -0.4, cor_sprm)
  ) -> gdsc_plot_ready

gdsc_plot_ready %>% 
  ggplot(aes(x = symbol, y = drug_name, color = cor_sprm)) +
  geom_point(aes(size = fdr)) +
  scale_color_gradient2(
    name = "Spearman Correlation",
    high = "red",
    mid = "white",
    low = "blue"
  ) +
  scale_size_continuous(
    name = "FDR"
  ) +
  ggthemes::theme_gdocs() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.2)
  ) +
  labs(
    x = "", y = ""
  ) -> p

ggsave(filename = 'gdsc_sig_drug.pdf', plot = p, device = "pdf", path = ccle_path, width = 10, height = 12)

# CTRP analysis
gene_list %>% 
  dplyr::inner_join(drug_ctrp, by = "symbol") -> ctrp_gene_list

fn_filter_drug_ctrp <- function(.x, .cor = 0.3, .p_val = 0.05) {
  .x %>% dplyr::filter(abs(cor_sprm) > .cor, p_val < .p_val)
}
ctrp_gene_list %>% 
  dplyr::mutate(cor_drug = purrr::map(.x = drug, .f = fn_filter_drug_ctrp)) %>% 
  tidyr::unnest(cor_drug) -> ctrp_gene_list_sig_drug

ctrp_gene_list_sig_drug %>% 
  dplyr::mutate(
    p_val = ifelse(-log10(p_val) > 40, 40, -log10(p_val)),
    cor_sprm = ifelse(cor_sprm > 0.4, 0.4, cor_sprm),
    cor_sprm = ifelse(cor_sprm < -0.4, -0.4, cor_sprm)
  ) -> ctrp_plot_ready

ctrp_plot_ready %>% 
  ggplot(aes(x = symbol, y = drug_name, color = cor_sprm)) +
  geom_point(aes(size = p_val)) +
  scale_color_gradient2(
    name = "Spearman Correlation",
    high = "red",
    mid = "white",
    low = "blue"
  ) +
  scale_size_continuous(
    name = "FDR"
  ) +
  ggthemes::theme_gdocs() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.2)
  ) +
  labs(
    x = "", y = ""
  ) -> p

ggsave(filename = 'ctrp_sig_drug.pdf', plot = p, device = "pdf", path = ccle_path, width = 16, height = 22)

