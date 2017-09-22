
# library -----------------------------------------------------------------

library(magrittr)
library(ggplot2)

# load path ---------------------------------------------------------------

tcga_path <- "/home/cliu18/liucj/projects/6.autophagy/TCGA"
afhl_path <- "/home/cliu18/liucj/projects/6.autophagy/09_afh_vs_afl"
afhl_class <- file.path(afhl_path, "01_af_h_l_classification")
expr_path <- "/home/cliu18/liucj/projects/6.autophagy/02_autophagy_expr/03_a_gene_expr"

# load data ---------------------------------------------------------------

rppa_expr <- readr::read_rds(file.path(tcga_path, "pancan33_rppa_expr_l4.rds.gz")) %>% 
  dplyr::filter(cancer_types != "PANCAN19")
rppa_name <- readr::read_rds(file.path(tcga_path, "rppa_name_symbol.rds.gz"))
gene_list <- readr::read_rds(file.path(expr_path, "rds_03_a_atg_lys_gene_list.rds.gz"))
pcc <- readr::read_tsv(file = file.path(tcga_path, "PanCanAtlasTumors_color_coded_by_organ_system_20170302.tsv"))

rppa_name %>% dplyr::semi_join(gene_list, by = "symbol") -> atg_rppa

atg_rppa %>% 
  dplyr::inner_join(gene_list, by = "symbol") %>% 
  dplyr::select(1,2, process) %>% 
  dplyr::arrange(process) -> sym_func
knitr::kable(sym_func)

PI3K_AKT <- "pS473|pT308|pS9|pS21S9|pT246|pT1462" %>% stringr::str_replace_all(pattern = "\\|", replacement = "$\\|") # AKT pS473, AKT pT308, GSK3 pS9,  GSK3 pS21-pS9, PRAS40 pT246 TSC2 pT1462

mTOR <- "pS2448|pT1135|pS65|pT37T46|pT70|pT389|pS235S236|pS240S244" %>% stringr::str_replace_all(pattern = "\\|", replacement = "$\\|") # mTOR pS2448, RICTOR pT1135, 4EBP1 pS65, 4EBP1 pT37T46, 4EBP1 pT70, S6K pT389, S6 pS235S236, S6 pS240S244

glue::glue("{PI3K_AKT}$|{mTOR}|P62|BECLIN|FOXO3A|") -> PI3K_AKT_mTOR

# mTORC1 complex MTOR, RPTOR, MLST8, DEPTOR
# phosphorylated (actiavated)
p_mTOR <- c("MTOR_pS2448")

unp_mTOR <- c("MTOR", "RAPTOR")

p_ampk <- c("AMPKALPHA_pT172")

unp_ampk <- c("AMPKALPHA")

sym_names <- c("BECLIN" = "BECN1", "MTOR_pS2448" = "MTOR", "P62LCKLIGAND" = "SQSTM1",  "AMPKALPHA" = "PRKAA1", "AMPKALPHA_pT172" = "PRKAA1")
.names <- c("AMPKALPHA" = "PRKAA1", "AMPKALPHA_pT172" = "PRKAA1 pT172", "BECLIN" = "BECN1", "MTOR_pS2448" = "MTOR pS2448", "P62LCKLIGAND" = "p62", "MTOR" = "MTOR", "BCL2" = "BCL2", "RAPTOR" = "RAPTOR")

# mTOR and P62, BECN1 and BCL2---------------------------- -----------------------------
fn_one_corr <- function(.z, .d){
  .d %>% dplyr::select(.z) -> .dd
  cor.test(.dd[[1]], .dd[[2]]) %>% broom::tidy() -> .corr
  tibble::tibble(p1 = .z[1], p2 = .z[2], coef = .corr$estimate, pval = .corr$p.value)
}

fn_all_corr <- function(.d) {
  .protein <- colnames(.d)
  combn(.protein, 2) -> .combn
  as.list(data.frame(.combn)) %>% purrr::map(.f = fn_one_corr, .d)
}

fn_select_marker <- function(.x){
  
  .names <- c("AMPKALPHA" = "PRKAA1", "AMPKALPHA_pT172" = "PRKAA1 pT172", "BECLIN" = "BECN1", "MTOR_pS2448" = "MTOR pS2448", "P62LCKLIGAND" = "p62")
  
  .x %>% 
    dplyr::filter(protein %in% c(p_mTOR, unp_mTOR, p_ampk, unp_ampk, "P62LCKLIGAND", "BECLIN", "BCL2")) %>% 
    dplyr::mutate(protein = plyr::revalue(x = protein, replace = .names)) %>% 
    dplyr::select(-symbol)
  }

rppa_expr %>% 
  dplyr::mutate(protein_expr = purrr::map(.x = protein_expr, .f = fn_select_marker)) -> 
  atg_rppa_expr

fn_protein_dist <- function(.gene) {

  atg_rppa_expr %>% 
    dplyr::mutate(
      protein_expr = purrr::map(
        .x = protein_expr, 
        .f = function(.x) dplyr::filter(.x, protein == .gene) %>% 
          tidyr::gather(key = barcode, value = rppa, -protein) %>% 
          tidyr::spread(key = protein, value = rppa)
        )
      ) %>% 
    tidyr::unnest() %>% 
    dplyr::rename(p62 = rlang::UQ(.gene)) %>% 
    dplyr::group_by(cancer_types) %>% 
    dplyr::mutate(rank = rank(p62)) %>% 
    dplyr::ungroup() -> p62
  
  p62 %>% 
    dplyr::group_by(cancer_types) %>% 
    dplyr::summarise(m = mean(p62)) %>% 
    dplyr::arrange(m) %>% 
    dplyr::pull(cancer_types) -> lev
  
  lev <- c("TGCT", "LGG", "PRAD", "UCS", "GBM", "PCPG", "PAAD", "KICH", "THCA", "MESO", "CESC", "KIRP", "SARC", "COAD", "THYM", "READ", "STAD", "UCEC", "BLCA", "HNSC", "BRCA", "KIRC", "ESCA", "SKCM", "OV", "LUAD", "CHOL", "LUSC", "UVM", "ACC", "DLBC", "LIHC")
  
  p62 %>% 
    dplyr::mutate(cancer_types = factor(cancer_types, levels = lev)) %>% 
    dplyr::mutate(p62 = dplyr::case_when(p62 > 3 ~ 3, p62 < -2 ~ -2, TRUE ~ p62)) %>% 
    ggplot(aes(x = cancer_types, y = p62)) +
    stat_boxplot(geom = 'errorbar', width = 0.3) +
    geom_boxplot(outlier.colour = NA) +
    geom_point(
      aes(color = cancer_types), 
      position = position_jitter(width = 0.05), 
      alpha = 0.4,
      size = 0.8) +
    scale_color_manual(
      name = "Cancer Types",
      values = dplyr::slice(pcc, match(lev, cancer_types)) %>% dplyr::pull(color)
    ) +
    theme(
      axis.line = element_line(color = "black"),
      axis.text.x = element_text(hjust = 1, vjust = 0.5, angle = 90),
      
      panel.grid = element_blank(),
      panel.background = element_rect(fill = "white", color = NA),
      panel.spacing.x = unit(0, "lines")
    ) +
    guides(color = F) +
    labs(
      x = "", 
      y = glue::glue("{.gene} RPPA score"), 
      title = "") -> plot_p62_rppa_dist
  
  ggsave(filename = file.path(afhl_class, glue::glue("fig_01_{.gene}_rppa_distribution_across_cancer_types.pdf")), 
         plot = plot_p62_rppa_dist,
         device = "pdf", 
         width = 7, 
         height = 3)
  
  plot_p62_rppa_dist
}

.names %>% purrr::map(fn_protein_dist) -> rppa_plot_dist
gridExtra::arrangeGrob %>% purrr::lift_dl(ncol = 1) %>% purrr::invoke(rppa_plot_dist) -> rppa_plot_dist_all
ggsave(filename = file.path(afhl_class, "fig_01_rppa_plot_dist_all.pdf"), plot = rppa_plot_dist_all, device = "pdf", width = 8, height = 15)

gridExtra::arrangeGrob(grobs = rppa_plot_dist[c(4, 5, 3)]) %>% 
  ggsave(filename = "fig_01_rppa_pmtor_p62_becn1.pdf", 
         plot = .,
         device = "pdf",
         width = 7, 
         height = 6)
#

# correlation with known proteins -------------------------------------

fn_corr_protein <- function(.x, .y){
  .corr_path <- file.path(afhl_class, "01_p62_mtor_ampk")
  if (!dir.exists(.corr_path)) dir.create(.corr_path)

  .x %>% 
    tidyr::gather(key = barcode, value = rppa, -protein) %>% 
    tidyr::spread(key = protein, value = rppa) %>% 
    dplyr::select(-barcode) -> .d
  
  .d %>% cor() -> .corr_matrix
  
  .title = glue::glue("Correlation of P62, BECLIN, AMPK and mTOR in {.y}")
  .filename = glue::glue("{.y}_p62_beclin_bcl2_mtor_ampk.pdf")
  
  ggcorrplot::ggcorrplot(.corr_matrix, hc.order = TRUE, outline.col = "white") +
    labs(title = .title) -> .p
  ggsave(filename = .filename, plot = .p, device = "pdf", path = .corr_path, width = 7, height = 7)
  
  list(
    p62_vs_BECN1 = cor.test(~p62 + BECN1, data = .d) %>% broom::tidy(),
    p62_vs_MTOR = cor.test(~p62 + MTOR, data = .d) %>% broom::tidy(),
    p62_vs_MTOR_pS2448 = cor.test(~p62 + `MTOR pS2448`, data = .d) %>% broom::tidy(),
    p62_vs_BCL2 = cor.test(~p62 + BCL2, data = .d) %>% broom::tidy(),
    BECN1_vs_BCL2 = cor.test(~BECN1 + BCL2, data = .d) %>% broom::tidy(),
    MTOR_vs_MTOR_pS2448 = cor.test(~MTOR + `MTOR pS2448`, data = .d) %>% broom::tidy()
  ) %>% 
    tibble::enframe() %>% 
    tidyr::unnest() %>% 
    dplyr::select(name, coef = estimate, pval = p.value) -> .d_corr
  
  .d_corr %>% 
    dplyr::filter(abs(coef) > 0.3, pval < 0.05) %>% 
    purrr::pmap(.f = function(name, coef, pval) {
      .xy <- stringr::str_split(name, "_vs_", simplify = T) %>% as.vector()
      
      .d %>% 
        ggplot(aes(x = rlang::eval_bare(rlang::sym(.xy[1])), y = rlang::eval_bare(rlang::sym(.xy[2])))) + 
        geom_point() + 
        annotate("text", x = 0, y = 0, 
                 label = glue::glue("R = {coef}
                                    p-value = {pval}")) +
        geom_smooth(se = F, method = "lm") +
        ggthemes::theme_gdocs() +
        labs(x = .xy[1], y = .xy[2]) 
        
    })
  
  
  
  .d_corr
}

atg_rppa_expr %>% 
  dplyr::mutate(corr = purrr::map2(.y = cancer_types, .x = protein_expr, .f = fn_corr_protein)) ->
  atg_rppa_expr

atg_rppa_expr %>% 
  tidyr::unnest(corr) %>% 
  ggplot(aes(x = cancer_types, y = name, fill = coef)) +
  geom_tile() +
  scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "red"
  )
























