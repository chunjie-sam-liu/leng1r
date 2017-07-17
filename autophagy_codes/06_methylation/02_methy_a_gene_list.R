library(ggplot2)
`%>%` <- magrittr::`%>%`


# Path
tcga_path <- "/home/cliu18/liucj/projects/6.autophagy/TCGA"
methy_path <- "/home/cliu18/liucj/projects/6.autophagy/06_methylation"
expr_path <- "/home/cliu18/liucj/projects/6.autophagy/02_autophagy_expr/"
expr_path_a <- file.path(expr_path, "03_a_gene_expr")
methy_box <- file.path(methy_path, "boxplot")


# load methylation and gene list
methy <- readr::read_rds(file.path(tcga_path, "pancan33_meth.rds.gz"))
gene_list <- readr::read_rds(file.path(expr_path_a, "rds_03_a_atg_lys_gene_list.rds.gz"))

# functions
filter_gene_list <- function(.x, gene_list) {
  gene_list %>%
    dplyr::select(symbol) %>%
    dplyr::left_join(.x, by = "symbol")
}

methy %>% 
  # dplyr::slice(2:5) %>%  # tidyr::unnest()
  dplyr::mutate(filter_methy = purrr::map(methy, filter_gene_list, gene_list = gene_list)) %>% 
  dplyr::select(-methy) -> gene_list_methy
#------------------------------------------------------------
fun_barcode <- function(.b){
  stringr::str_sub(
    string = .b,
    start = 1,
    end = 12
  )
}
fun_tn_type <- function(.b){
  type <- .b %>% 
    stringr::str_split(pattern = "-", simplify = T) %>% 
    .[, 4] %>% 
    stringr::str_sub(1, 2)
}
fun_boxplot <- function(fig_name, data, path = methy_box){
  data %>% 
    ggpubr::ggboxplot(x = 'type', y = 'meth', color = 'type', pallete = 'jco' ) +
    ggpubr::stat_compare_means(
      method = "t.test", 
      label.y = 1, 
      label.x = 1.2) + 
    ggthemes::theme_gdocs() +
    scale_color_manual(values = c("#DC3912", "#3366CC")) +
    labs(y = "B-value", x = "", title = fig_name) -> p
  ggsave(filename = paste(fig_name, "pdf", sep = "."), p, device = "pdf", path = path, width = 4, height = 3)
}
fun_compare <- function(.x, .y ){
  .x %>% 
    dplyr::mutate(gene = as.character(gene )) %>% 
    tidyr::gather(key = barcode, value = meth, -symbol, -gene) %>% 
    dplyr::mutate(type = fun_tn_type(barcode)) %>% 
    dplyr::filter(type %in% c("01", "11")) %>% 
    dplyr::mutate(barcode = fun_barcode(barcode)) %>% 
    dplyr::select(symbol, gene, barcode, meth, type) %>% 
    dplyr::mutate(type = dplyr::case_when(
      type == "01" ~ "Tumor",
      type == "11" ~ "Normal"
    )) %>% 
    dplyr::filter(!is.na(gene)) -> .d
  if(nrow(.d) < 20 || length(unique(.d$type)) != 2){return(tibble::tibble())}
  # at least 10 samples
  .d %>% 
    dplyr::select(barcode, type) %>% 
    dplyr::distinct() %>%
    dplyr::group_by(type) %>% 
    dplyr::count() %>% 
    dplyr::pull(n) -> sample_num
  if(any(sample_num < 10)){return(tibble::tibble())}
  
  .d %>% 
    dplyr::group_by(symbol, gene) %>% 
    dplyr::do(
      broom::tidy(
        t.test(meth ~ type, data = .)
      )
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(fdr = p.adjust(p.value, method = "fdr")) %>% 
    dplyr::filter(fdr < 0.05) %>%
    dplyr::mutate(
      direction = dplyr::case_when(
          estimate > 0 ~ 0, # normal high
          estimate < 0 ~ 1 # tumor high
        )) %>% 
    dplyr::select(symbol, gene, direction, p_val = p.value, fdr) -> .d_out
  
  # draw every pic
  .d %>%
    dplyr::semi_join(.d_out, by = c("symbol", "gene")) %>%
    # dplyr::filter(symbol %in% c("ATP6V0D1", "ATP6V0A4")) %>%
    tidyr::nest(-symbol, -gene) %>%
    dplyr::mutate(fig_name = paste(.y, symbol, sep = "_")) %>% 
    dplyr::select(fig_name, data) %>% 
    purrr::pwalk(.f = fun_boxplot, path = methy_box)
  
  return(.d_out)
}


# gene_list_methy %>% 
#   dplyr::slice(1) %>%
#   dplyr::mutate(methy_comparison = purrr::map(.x = filter_methy, .y = cancer_types, .f = fun_compare)) %>% 
#   dplyr::select(-filter_methy) %>% 
#   tidyr::unnest() -> gene_list_methy_fdr

cl <- parallel::detectCores()
cluster <- multidplyr::create_cluster(floor(cl * 5 / 6))
gene_list_methy %>% 
  multidplyr::partition(cluster = cluster) %>%
  multidplyr::cluster_library("magrittr") %>%
  multidplyr::cluster_library("ggplot2") %>%
  multidplyr::cluster_assign_value("fun_barcode", fun_barcode)  %>%
  multidplyr::cluster_assign_value("fun_tn_type", fun_tn_type)  %>%
  multidplyr::cluster_assign_value("fun_boxplot", fun_boxplot)  %>%
  multidplyr::cluster_assign_value("fun_compare", fun_compare)  %>%
  multidplyr::cluster_assign_value("methy_box", methy_box)  %>%
  dplyr::mutate(methy_comparison = purrr::map(.x = filter_methy, .y = cancer_types, .f = fun_compare)) %>% 
  dplyr::collect() %>%
  dplyr::as_tibble() %>%
  dplyr::ungroup() %>%
  dplyr::select(-PARTITION_ID) %>% 
  dplyr::select(-filter_methy) %>% 
  tidyr::unnest() -> gene_list_methy_fdr
on.exit(parallel::stopCluster(cluster))
readr::write_rds(gene_list_methy_fdr, path = file.path(methy_path, ".rds_02_gene_list_methy_fdr.rds.gz"), compress = "gz")


gene_list_methy_fdr %>% 
  dplyr::mutate(fdr = -log10(fdr)) %>% 
  dplyr::mutate(fdr = ifelse(fdr > 50, 50, fdr)) -> plot_ready
  
plot_ready %>% 
  ggplot(aes(x = cancer_types, y = symbol)) +
  geom_point(aes(size = fdr, color = as.factor(direction))) +
  # scale_x_discrete(limit = cancer_rank$cancer_types) +
  # scale_y_discrete(limit = gene_rank$symbol) +
  scale_size_continuous(name = "FDR") +
  theme(
    panel.background = element_rect(colour = "black", fill = "white"),
    panel.grid = element_line(colour = "grey", linetype = "dashed"),
    panel.grid.major = element_line(
      colour = "grey",
      linetype = "dashed",
      size = 0.2
    ),
    
    axis.title = element_blank(),
    axis.ticks = element_line(color = "black"),
    # axis.text.y = element_text(color = gene_rank$color),
    
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14),
    legend.key = element_rect(fill = "white", colour = "black")
  ) 


#
save.image(file = file.path(methy_path, ".rda_02_methy_a_gene_list.rda"))
load(file = file.path(methy_path, ".rda_02_methy_a_gene_list.rda"))

