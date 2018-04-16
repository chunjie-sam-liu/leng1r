library(magrittr)
tcga_path = "/home/cliu18/liucj/projects/6.autophagy/TCGA"
expr_path <- "/home/cliu18/liucj/projects/6.autophagy/02_autophagy_expr/"
expr_path_a <- file.path(expr_path, "03_a_gene_expr")
cnv_path <- "/home/cliu18/liucj/projects/6.autophagy/03_cnv"

# load cnv and gene list
cnv <- readr::read_rds(file.path(tcga_path, "pancan34_cnv_threshold.rds.gz")) %>% 
  dplyr::filter(cancer_types != "COADREAD")
marker_file <- readr::read_rds(file.path(expr_path_a, "rds_03_a_atg_lys_marker.rds.gz"))
gene_list <- readr::read_rds(file.path(expr_path_a, "rds_03_a_atg_lys_gene_list.rds.gz")) %>% 
  dplyr::left_join(marker_file, by = "symbol")

gene_list %>% 
  dplyr::filter(pathway == "autophagesome formation-core") %>% 
  dplyr::mutate(color = ifelse(is.na(marker), "black", "red")) -> atg_core

filter_gene_list <- function(.x, gene_list) {
  gene_list %>%
    dplyr::select(symbol) %>%
    dplyr::left_join(.x, by = "symbol")
}

cnv %>%
  dplyr::mutate(filter_cnv = purrr::map(cnv, filter_gene_list, gene_list = gene_list)) %>%
  dplyr::select(-cnv) -> gene_list_cnv

readr::write_rds(x = gene_list_cnv, path = file.path(cnv_path, ".rds_02_cnv_a_gene_list.rds.gz"), compress = "gz")

fn_get_amplitue_threshold <- function(.x){
  tibble::tibble(
    a_total = sum(.x > 0) / length(.x), 
    d_total = sum(.x < 0) / length(.x),
    a_homo = sum(.x == 2) / length(.x),
    d_homo = sum(.x == -2) / length(.x),
    a_hete = sum(.x == 1) / length(.x),
    d_hete = sum(.x == -1) / length(.x)) 
}
fn_get_ad <- function(.d){
  .d %>% 
    unlist(use.name = F) %>% 
    fn_get_amplitue_threshold()
}
fn_get_percent <- function(cancer_types, filter_cnv){
  filter_cnv %>%
    tidyr::nest(-symbol) %>% 
    dplyr::mutate(ad = purrr::map(data, .f = fn_get_ad)) %>% 
    dplyr::select(-data) %>% 
    tidyr::unnest(ad) %>% 
    tibble::add_column(cancer_types = cancer_types, .before = 1)
}

fn_gen_combined_core_atg <- function(cancer_types, filter_cnv){
  # cancer_types <- te$cancer_types
  # filter_cnv <- te$filter_cnv[[1]]
  filter_cnv %>% 
    dplyr::semi_join(atg_core, by = "symbol") %>% 
    tidyr::drop_na() %>% 
    tidyr::gather(key = barcode, value = gistic, - symbol) %>% 
    tidyr::spread(key = symbol, value = gistic) %>% 
    dplyr::select(- barcode) -> .d
  
  n_sample <- nrow(.d)
  
  .d %>% 
    dplyr::filter_all(.vars_predicate = dplyr::any_vars(. == -2)) %>% 
    nrow() -> .del
  
  .d %>% 
    dplyr::filter_all(.vars_predicate = dplyr::any_vars(. == 2)) %>% 
    nrow() -> .amp
  
  tibble::tibble(del = .del / n_sample, amp = .amp / n_sample)
}

gene_list_cnv %>% 
  dplyr::filter(cancer_types == "TGCT") %>% 
  dplyr::mutate(rs = purrr::map2(cancer_types, filter_cnv, fn_get_percent)) %>% 
  dplyr::select(-cancer_types, -filter_cnv) %>% 
  tidyr::unnest(rs)

cl <- parallel::detectCores()
cluster <- multidplyr::create_cluster(floor(cl * 5 / 6))
gene_list_cnv %>% 
  multidplyr::partition(cluster = cluster) %>%
  multidplyr::cluster_library("magrittr") %>%
  multidplyr::cluster_assign_value("fn_get_amplitue_threshold", fn_get_amplitue_threshold)  %>%
  multidplyr::cluster_assign_value("fn_get_ad", fn_get_ad) %>% 
  multidplyr::cluster_assign_value("fn_get_percent", fn_get_percent) %>% 
  multidplyr::cluster_assign_value("cnv_path", cnv_path) %>% 
  dplyr::mutate(rs = purrr::map2(cancer_types, filter_cnv, fn_get_percent)) %>% 
  dplyr::collect() %>%
  dplyr::as_tibble() %>%
  dplyr::ungroup() %>%
  dplyr::select(-PARTITION_ID) %>% 
  dplyr::select(-cancer_types, -filter_cnv) %>% 
  tidyr::unnest(rs) -> gene_list_cnv_per
parallel::stopCluster(cluster)


library(ggplot2)
gene_list_cnv_per %>% 
  tidyr::drop_na() %>% 
  dplyr::mutate(other = 1 - a_total - d_total) -> plot_ready

plot_ready %>% 
  dplyr::semi_join(atg_core, by = "symbol") %>% 
  dplyr::select(-a_total, -d_total) -> atg_core_plot_ready

atg_core_plot_ready %>% 
  dplyr::group_by(cancer_types) %>% 
  dplyr::summarise(v = sum(1-other)) %>% 
  dplyr::arrange(dplyr::desc(v)) -> cancer_rank
  
atg_core_plot_ready %>% 
  dplyr::group_by(symbol) %>% 
  dplyr::summarise(v = sum(1 - other)) %>% 
  dplyr::left_join(atg_core, by = "symbol") %>% 
  dplyr::arrange(v) -> gene_rank

atg_core_plot_ready %>% 
  # dplyr::filter(symbol %in% c("GABARAPL3", "ULK1")) %>% 
  tidyr::gather(key = type, value = per, -c(cancer_types, symbol)) %>% 
  dplyr::mutate(
    symbol = factor(x = symbol, levels = gene_rank$symbol), 
    cancer_types = factor(x = cancer_types, levels = cancer_rank$cancer_types)) %>% 
  ggplot(aes(x = factor(1), y  = per, fill = type)) +
  geom_bar(stat = 'identity', position = "stack", color = NA) +
  # scale_y_continuous(limits = c(0,1))
  coord_polar("y") +
  facet_grid(symbol ~ cancer_types) +
  theme(axis.text=element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        
        strip.text.y = element_text(angle =0,hjust=0,color=gene_rank$color,size=11),
        strip.text.x = element_text(size=11,angle = 90,vjust = 0),
        strip.background = element_blank(),
        
        legend.title= element_blank(), 
        legend.text = element_text(size=11),
        legend.position = "bottom",
        
        panel.background = element_blank(),
        panel.spacing  = unit(0.02, "lines")) +
  scale_fill_manual(
    limits = c("a_hete", "a_homo", "d_hete", "d_homo", "other"),
    label=c("Hete Amp","Homo_Amp", "Hete Del", "Homo Del", "None"), 
    # Amp RColorBrewer name = "Spectral"
    # Del RColorBrewer name = "BrBG"
    values=c("#D53E4F", "#9E0142", "#01665E", "#003C30", "grey")) -> p

ggsave(filename = "05_core_per.pdf", plot = p, device = "pdf", path = cnv_path, width = 8, height = 10)

atg_core_plot_ready %>% 
  dplyr::select(cancer_types, symbol, a_homo, d_homo) %>% 
  tidyr::gather(key = type, value = per, -cancer_types, -symbol) %>% 
  dplyr::filter(per >= 0.05) %>% 
  dplyr::mutate(per = ifelse(per > 0.2, 0.2, per)) -> homo_core_plot_ready
homo_core_plot_ready %>% 
  dplyr::group_by(cancer_types) %>% 
  dplyr::summarise(v = sum(per)) %>% 
  dplyr::arrange(desc(v)) -> cancer_rank
homo_core_plot_ready %>% 
  dplyr::group_by(symbol) %>% 
  dplyr::summarise(v= sum(per)) %>% 
  dplyr::left_join(atg_core, by = "symbol") %>% 
  dplyr::mutate(color = ifelse(is.na(marker), "black", "red")) %>% 
  dplyr::arrange(v) -> gene_rank
homo_core_plot_ready %>% 
  ggplot(aes(y = symbol, x = cancer_types)) +
  geom_point(aes(size = per, color = type)) +
  scale_x_discrete(limits = cancer_rank$cancer_types) +
  scale_y_discrete(limits = gene_rank$symbol) +
  scale_size_continuous(
    name = "Percent",
    breaks = c(0.05, 0.1, 0.15, 0.2),
    limits = c(0.05, 0.2),
    labels = c("5", "10", "15", "20")
  ) +
  ggsci::scale_color_npg(
    name = "Type",
    limits = c("a_homo", "d_homo"),
    labels = c("Amp", "Del")
  ) +
  ggthemes::theme_gdocs() +
  theme(
    axis.text.y = element_text(color = gene_rank$color),
    axis.text.x = element_text(angle = 90, vjust = 0.3, hjust = 1)
  ) +
  labs(x = "", y = "") -> p
ggsave(filename = "05_core_homo.pdf", plot = p, device = 'pdf', path = cnv_path, width = 6, height = 5)


plot_ready %>% 
  dplyr::left_join(gene_list, by = "symbol") %>% 
  dplyr::filter(status == "a") %>% 
  dplyr::select(cancer_types, symbol, a_homo, d_homo) %>% 
  tidyr::gather(key = type, value = per, a_homo, d_homo) %>% 
  dplyr::filter(per >= 0.05) %>% 
  dplyr::mutate(per = ifelse(per > 0.4, 0.4, per)) -> atg_plot_ready

atg_plot_ready %>% 
  dplyr::group_by(cancer_types) %>% 
  dplyr::summarise(v = sum(per)) %>% 
  dplyr::arrange(desc(v)) -> cancer_rank
atg_plot_ready %>% 
  dplyr::group_by(symbol) %>% 
  dplyr::summarise(v= sum(per)) %>% 
  dplyr::left_join(atg_core, by = "symbol") %>% 
  dplyr::mutate(color = ifelse(is.na(marker), "black", "red")) %>% 
  dplyr::arrange(v) -> gene_rank

atg_plot_ready %>% 
  ggplot(aes(y = symbol, x = cancer_types)) +
  geom_point(aes(size = per, color = type)) +
  scale_x_discrete(limits = cancer_rank$cancer_types) +
  scale_y_discrete(limits = gene_rank$symbol) +
  scale_size_continuous(
    name = "Percent",
    breaks = c(0.05, 0.1, 0.25, 0.4),
    limits = c(0.05, 0.4),
    labels = c("5", "10", "25", "40")
  ) +
  ggsci::scale_color_npg(
    name = "Type",
    limits = c("a_homo", "d_homo"),
    labels = c("Amp", "Del")
  ) +
  ggthemes::theme_gdocs() +
  theme(
    axis.text.y = element_text(color = gene_rank$color),
    axis.text.x = element_text(angle = 90, vjust = 0.3, hjust = 1)
  ) +
  labs(x = "", y = "") -> p
ggsave(filename = "05_atg_homo.pdf", plot = p, device = 'pdf', path = cnv_path, width = 7, height = 8) 


plot_ready %>% 
  dplyr::left_join(gene_list, by = "symbol") %>% 
  dplyr::filter(status == "l") %>% 
  dplyr::select(cancer_types, symbol, a_homo, d_homo) %>% 
  tidyr::gather(key = type, value = per, a_homo, d_homo) %>% 
  dplyr::filter(per >= 0.05) %>% 
  dplyr::mutate(per = ifelse(per > 0.2, 0.2, per)) -> lys_plot_ready

lys_plot_ready %>% 
  dplyr::group_by(cancer_types) %>% 
  dplyr::summarise(v = sum(per)) %>% 
  dplyr::arrange(desc(v)) -> cancer_rank
lys_plot_ready %>% 
  dplyr::group_by(symbol) %>% 
  dplyr::summarise(v= sum(per)) %>% 
  dplyr::left_join(gene_list, by = "symbol") %>% 
  # dplyr::mutate(color = ifelse(is.na(marker), "black", "red")) %>% 
  dplyr::arrange(v) -> gene_rank

lys_plot_ready %>% 
  ggplot(aes(y = symbol, x = cancer_types)) +
  geom_point(aes(size = per, color = type)) +
  scale_x_discrete(limits = cancer_rank$cancer_types) +
  scale_y_discrete(limits = gene_rank$symbol) +
  scale_size_continuous(
    name = "Percent",
    breaks = c(0.05, 0.1, 0.15, 0.2),
    limits = c(0.05, 0.4),
    labels = c("5", "10", "15", "20")
  ) +
  ggsci::scale_color_npg(
    name = "Type",
    limits = c("a_homo", "d_homo"),
    labels = c("Amp", "Del")
  ) +
  ggthemes::theme_gdocs() +
  theme(
    axis.text.y = element_text(color = gene_rank$color),
    axis.text.x = element_text(angle = 90, vjust = 0.3, hjust = 1)
  ) +
  labs(x = "", y = "") -> p
ggsave(filename = "05_lys_homo.pdf", plot = p, device = 'pdf', path = cnv_path, width = 7, height = 8) 


# combined core atg

fn_gen_combined_core_atg <- function(cancer_types, filter_cnv, g_list){
  # cancer_types <- te$cancer_types
  # filter_cnv <- te$filter_cnv[[1]]
  filter_cnv %>% 
    dplyr::semi_join(g_list, by = "symbol") %>% 
    tidyr::drop_na() %>% 
    tidyr::gather(key = barcode, value = gistic, - symbol) %>% 
    tidyr::spread(key = symbol, value = gistic) %>% 
    dplyr::select(- barcode) -> .d
  
  n_sample <- nrow(.d)
  
  .d %>% 
    dplyr::filter_all(.vars_predicate = dplyr::any_vars(. == -2)) %>% 
    nrow() -> .del
  
  .d %>% 
    dplyr::filter_all(.vars_predicate = dplyr::any_vars(. == 2)) %>% 
    nrow() -> .amp
  
  tibble::tibble(del = .del / n_sample, amp = .amp / n_sample)
}

gene_list_cnv %>% 
  # dplyr::filter(cancer_types == "TGCT") %>% 
  dplyr::mutate(rs = purrr::map2(cancer_types, filter_cnv, fn_gen_combined_core_atg, g_list = atg_core)) %>% 
  dplyr::select( -filter_cnv) %>% 
  tidyr::unnest(rs) -> gene_combined_per

gene_combined_per %>% 
  dplyr::mutate(del = -del) %>% 
  tidyr::gather(key = type, value = per, -cancer_types) %>% 
  ggplot(aes(x = reorder(x = cancer_types, X = per, Fun = function(x) sum(x)), y = per, fill = type)) +
  geom_bar(stat = 'identity', position = "stack") +
  ggsci::scale_fill_npg(
    name = "Type",
    limits = c("amp", "del"),
    labels = c("Amp", "Del")
  ) +
  ggthemes::theme_gdocs() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.3)
  ) +
  labs(x = "Cancer Types", y = "Core ATG Frequency") -> p

ggsave(filename = "05_combined_atg_core.pdf", plot = p, device = 'pdf', width = 7, height = 5, path = cnv_path)


gene_list_cnv %>% 
  dplyr::mutate(rs = purrr::map2(cancer_types, filter_cnv, fn_gen_combined_core_atg, g_list = gene_list %>% dplyr::filter(status == "l"))) %>% 
  dplyr::select( -filter_cnv) %>% 
  tidyr::unnest(rs) -> lys_combined_per

lys_combined_per %>% 
  dplyr::mutate(del = -del) %>% 
  tidyr::gather(key = type, value = per, -cancer_types) %>% 
  ggplot(aes(x = reorder(x = cancer_types, X = per, Fun = function(x) sum(x)), y = per, fill = type)) +
  geom_bar(stat = 'identity', position = "stack") +
  ggsci::scale_fill_npg(
    name = "Type",
    limits = c("amp", "del"),
    labels = c("Amp", "Del")
  ) +
  ggthemes::theme_gdocs() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.3)
  ) +
  labs(x = "Cancer Types", y = "Lysosome Genes Frequency") -> p

ggsave(filename = "05_combined_lys.pdf", plot = p, device = 'pdf', width = 7, height = 5, path = cnv_path)

# mutual exclusive test
# using raphael-group R package cometExactTest
# the caculation is on in 02_cnv_b_mutual_exclusive.R
# the result is save to the rds ".rds_02_cnv_b_mutual_exclusive.rds.gz"

mutual_exclusive <- readr::read_rds(path = file.path(cnv_path, ".rds_02_cnv_b_mutual_exclusive.rds.gz"))

fn_mut_exc_draw <- function(g1, g2, p_val, fdr, filter_cnv, cancer_types){
  # g1 <- "LAMP2"
  # g2 <- "MAP1LC3A" 
  # p_val <- 3.229587^-6
  fdr <- signif(fdr, digits = 3)
  
  filter_cnv %>% 
    dplyr::filter(symbol %in% c(g1, g2)) %>% 
    tidyr::gather(key = barcode, value = gistic, -symbol) -> plot_ready

  plot_ready %>% 
    tidyr::spread(key = symbol, value = gistic) -> rank_ready
  
  rank_ready %>% 
    dplyr::filter_if(.predicate = is.integer, .vars_predicate = dplyr::any_vars(. != 0)) %>% 
    nrow() / nrow(rank_ready) -> diff
  
  fig_name <- paste(cancer_types, fdr, signif(diff, digits = 4), g1, g2, sep = "_")
  
  # 
  # sum(rank_ready %>% dplyr::pull(g1) == 2) -> g1_amp
  # sum(rank_ready %>% dplyr::pull(g2) == 2) -> g2_amp
  # 
  # if(g1_amp > g2_amp){
  #   rank_ready %>% 
  #     dplyr::filter_if(.predicate = is.integer, .vars_predicate = dplyr::any_vars(. == 2)) %>%  
  #     dplyr::arrange(dplyr::desc(get(g1, .)), dplyr::desc(get(g2, .))) %>% 
  #     dplyr::pull(barcode) -> plot_rank
  #   
  #   rank_ready %>% 
  #     dplyr::filter_if(.predicate = is.integer, .vars_predicate = dplyr::all_vars(. != 2)) %>% 
  #     dplyr::arrange(dplyr::desc(get(g2, .)), dplyr::desc(get(g1, .))) %>% 
  #     dplyr::pull(barcode) %>% 
  #     c(plot_rank, .) -> plot_rank
  # } else{
  #   rank_ready %>% 
  #     dplyr::filter_if(.predicate = is.integer, .vars_predicate = dplyr::any_vars(. == 2)) %>%  
  #     dplyr::arrange(dplyr::desc(get(g2, .)), dplyr::desc(get(g1, .))) %>% 
  #     dplyr::pull(barcode) -> plot_rank
  #   
  #   rank_ready %>% 
  #     dplyr::filter_if(.predicate = is.integer, .vars_predicate = dplyr::all_vars(. != 2)) %>% 
  #     dplyr::arrange(dplyr::desc(get(g1, .)), dplyr::desc(get(g2, .))) %>% 
  #     dplyr::pull(barcode) %>% 
  #     c(plot_rank, .) -> plot_rank
  # }
  # 
  # plot_ready %>% 
  #   ggplot(aes(x = barcode, y = factor(1), fill = as.factor(gistic))) +
  #   geom_tile() +
  #   scale_x_discrete(limits = plot_rank) +
  #   facet_grid(symbol ~ .) +
  #   scale_fill_manual(
  #     name = "Type",
  #     limits = c("2", "1", "0", "-1", "-2"),
  #     label=c("Homozygous Amplification","Heterozygous Amplification", "NO", "Heterozygous Deletion", "Homozygous Deletion"), 
  #     # Amp RColorBrewer name = "Spectral"
  #     # Del RColorBrewer name = "BrBG"
  #     values= c("#FF0000", "#FF6A6A", "#D6D6D6", "#7EC0EE", "#27408B")) +
  #   theme(
  #     axis.text=element_blank(),
  #     axis.title = element_blank(),
  #     axis.ticks = element_blank(),
  #     
  #     strip.text.y = element_text(angle =0,hjust=0,size=11),
  #     strip.text.x = element_text(size=11,angle = 90,vjust = 0),
  #     strip.background = element_blank(),
  #     
  #     legend.title= element_blank(), 
  #     legend.text = element_text(size=11),
  #     legend.position = "bottom",
  #     
  #     panel.background = element_blank(),
  #     panel.spacing  = unit(0.02, "lines")
  #   ) +
  #   labs(x = "", y = "", title = paste(g1, "and", g2, " mutual exclusive in", cancer_types, ", FDR", fdr)) -> p
  # 
  # ggsave(filename = paste(fig_name, 'pdf', sep = "."), plot = p, device = 'pdf', path = file.path(cnv_path, "mutual_exclusive"), width = 10, height = 2)
  

  tibble::tibble(fig_name = fig_name) %>% 
    tidyr::separate(col = fig_name, into = c("cancer_types", "fdr", "diff_prop", "g1", "g2"), sep = "_")
}
fn_mut_exc_filter <- function(cancer_types, comet, filter_cnv){
  # cancer_types <- te$cancer_types
  # comet <- te$comet[[1]]
  # filter_cnv <- te$filter_cnv[[1]]
  
  comet %>% 
    purrr::pmap(.f = fn_mut_exc_draw, filter_cnv = filter_cnv, cancer_types = cancer_types) %>% 
    dplyr::bind_rows()
}

mutual_exclusive %>% 
  dplyr::filter(fdr < 0.05) %>%
  tidyr::nest(-cancer_types, .key = comet) %>% 
  dplyr::inner_join(gene_list_cnv, by = "cancer_types") %>%
  purrr::pmap(.f = fn_mut_exc_filter) %>% 
  dplyr::bind_rows() -> mutual_exclusive_filter_prop
readr::write_rds(x = mutual_exclusive_filter_prop, path = file.path(cnv_path, ".rds_02_cnv_b_mutual_exclusive_prop.rds.gz"), compress = "gz")


mutual_exclusive_filter_prop %>% 
  dplyr::mutate(label = paste(g1, g2, sep = "_")) %>% 
  dplyr::mutate(fdr = as.double(fdr), diff_prop = as.double(diff_prop)) %>% 
  dplyr::mutate(fdr = -log10(fdr)) -> plot_ready

plot_ready %>% 
  ggplot(aes(y = fdr, x = diff_prop, color = cancer_types)) +
  geom_point() +
  ggsci::scale_color_npg(name = "Cancer Types") +
  ggrepel::geom_text_repel(
    data = plot_ready %>% dplyr::mutate(v = diff_prop * 8 + fdr) %>% dplyr::top_n(10, v),
    aes(label = label),
    size = 3) +
  ggthemes::theme_gdocs() +
  labs(x = "Combined Proportion", y = "FDR") -> p
ggsave(filename = '05_mutual_exclusive_proportion.pdf', plot = p, device = 'pdf', path = cnv_path, width = 8, height = 6 )


save.image(file = file.path(cnv_path, ".rda_02_cnv_a_gene_list.rda"))
load(file = file.path(cnv_path, ".rda_02_cnv_a_gene_list.rda"))


