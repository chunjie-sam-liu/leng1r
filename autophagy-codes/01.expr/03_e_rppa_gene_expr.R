library(magrittr)
library(ggplot2)

expr_path <- "/home/cliu18/liucj/projects/6.autophagy/02_autophagy_expr/"
rppa_path <- file.path(expr_path, "03_e_rppa")
tcga_path <- "/home/cliu18/liucj/projects/6.autophagy/TCGA"
expr_path <- file.path(expr_path, "03_a_gene_expr")

# rppa <- readr::read_rds(path = file.path(tcga_path,"pancan_clinical_stage.rds.gz")) %>% 
  # dplyr::filter(n >= 40) %>% 
  # dplyr::select(-n)

gene_list <- readr::read_rds(file.path(expr_path, "rds_03_a_atg_lys_gene_list.rds.gz"))
gene_list_expr <- readr::read_rds(path = file.path(expr_path, ".rds_03_a_gene_list_expr.rds.gz"))

gene_list_expr_median_cluster <- readr::read_rds(path = file.path(expr_path, ".rds_01_b_rppa_median_cluster_expr.rds.gz"))
rppa_scores <- readr::read_rds(file.path(tcga_path, "pancan_rppa_score.rds.gz"))

fn_cluster_rppa <- function(median_cluster, rppa){
  median_cluster %>% 
    dplyr::filter(!is.na(expr)) %>% 
    dplyr::mutate(cluster = as.factor(cluster)) %>% 
    dplyr::mutate(cluster = plyr::revalue(cluster, replace = c("1" = "down", "2" = "up"))) %>% 
    dplyr::inner_join(rppa, by = "barcode") %>% 
    dplyr::distinct() -> merged_clean
}
fn_test <- function(merged_clean, cancer_types){
  print(cancer_types)
  merged_clean %>%
    dplyr::group_by(symbol, pathway) %>% 
    dplyr::do(
      broom::tidy(
        tryCatch(
          t.test(score ~ cluster, data = .),
          error = function(e){1},
          warning = function(e){1})
      )
    ) %>% 
    dplyr::select(symbol, pathway, p.value) %>%
    dplyr::group_by(symbol) %>% 
    dplyr::mutate(fdr = p.adjust(p.value, method = "fdr")) %>% 
    dplyr::mutate(bfi = p.adjust(p.value, method = "bonferroni")) %>% 
    dplyr::ungroup() -> clean_pval
  #9761, 11901
  merged_clean %>% 
    tidyr::spread(key = cluster, value = score) %>% 
    dplyr::group_by(symbol, pathway) %>% 
    dplyr::summarise(diff = mean(up, na.rm = T) - mean(down, na.rm = T)) %>% 
    dplyr::ungroup() -> clean_diff
  
  clean_pval %>% 
    dplyr::inner_join(clean_diff, by = c("symbol", "pathway")) 
}

gene_list_expr_median_cluster %>% 
  dplyr::inner_join(rppa_scores, by = "cancer_types") -> gene_rppa 
# 
# gene_rppa %>%
#   dplyr::filter(cancer_types == "GBM") %>% 
#   dplyr::mutate(merged_clean = purrr::map2(median_cluster, rppa, fn_cluster_rppa)) %>%
#   dplyr::select(cancer_types, merged_clean) %>%
#   dplyr::mutate(diff_pval = purrr::map2(merged_clean, cancer_types, fn_test))

cl <- parallel::detectCores()
cluster <- multidplyr::create_cluster(floor(cl * 5 / 6))
gene_rppa %>%
  multidplyr::partition(cluster = cluster) %>%
  multidplyr::cluster_library("magrittr") %>%
  multidplyr::cluster_assign_value("fn_cluster_rppa", fn_cluster_rppa)  %>%
  multidplyr::cluster_assign_value("fn_test", fn_test) %>% 
  dplyr::mutate(merged_clean = purrr::map2(median_cluster, rppa, fn_cluster_rppa)) %>% 
  dplyr::select(cancer_types, merged_clean) %>% 
  dplyr::mutate(diff_pval = purrr::map2(merged_clean, cancer_types, fn_test)) %>% 
  dplyr::collect() %>%
  dplyr::as_tibble() %>%
  dplyr::ungroup() %>%
  dplyr::select(-PARTITION_ID) -> gene_rppa_sig_pval
on.exit(parallel::stopCluster(cluster))

gene_rppa_sig_pval %>% 
  readr::write_rds(path = file.path(rppa_path, ".rds_03_e_rppa_gene_expr.rds.gz"), compress = "gz")

#-------------------------------------------------------------------------------------------------
pathway_replace <- c(
  "PI3KAKT"="PI3K/AKT",
  "RASMAPK"="RAS/MAPK",
  "TSCmTOR"="TSC/mTOR",
  "CellCycle"="Cell Cycle",
  "DNADamage"="DNA Damage Response"
)

pathway_name <- c("Apoptosis","Cell Cycle","DNA Damage Response","EMT","Hormone AR","Hormone ER","PI3K/AKT","RAS/MAPK","RTK","TSC/mTOR")

gene_rppa_sig_pval %>% 
  tidyr::unnest(diff_pval) %>% 
  dplyr::filter(!is.na(p.value)) %>% 
  dplyr::mutate(pathway = plyr::revalue(pathway, pathway_replace)) %>% 
  dplyr::mutate(class = ifelse(fdr < 0.05 & diff > 0, "Activation", "None")) %>% 
  dplyr::mutate(class = ifelse(fdr < 0.05 & diff < 0, "Inhibition", class)) -> gene_rppa_sig_pval_class

fn_draw_pathway <- function(symbol, data){
  print(symbol)
  fig_name <- paste(symbol, "pathway.pdf", sep = "_")
  data %>%
    ggplot(aes(x = cancer_types, y = pathway)) +
    geom_tile(aes(fill = factor(class)), col = "white") +
    scale_fill_manual(
      limits=c("Activation","Inhibition","None"),
      values = c("red","blue","lightgray"),
      na.value="white",
      labels=c("Activation","Inhibition","None"),
      name= symbol) +
    scale_y_discrete(limits=pathway_name, label=pathway_name) +
    theme(panel.background=element_rect(colour=NA,fill="white",size=2),
          panel.grid=element_line(colour="white",linetype="dashed"),
          panel.grid.major=element_line(colour="white",linetype="dashed"),
          axis.title=element_blank(),
          axis.text.y=element_text(size=16,colour = "black"),
          axis.text.x=element_text(size=16,colour = "black",angle=90,vjust=0.5,hjust=1),
          axis.ticks=element_line(color="black"),
          legend.text=element_text(size=16),
          legend.title=element_text(size=14),legend.position="bottom",legend.direction="horizontal") -> p
  ggsave(filename = fig_name, plot = p, device = "pdf", path = file.path(rppa_path, "gene_tileplot"), width = 9, height = 4.5)
}
# save picture
gene_rppa_sig_pval_class %>%
  dplyr::select("symbol", "pathway", "cancer_types", "class", "p.value", "fdr") %>%
  tidyr::nest(-symbol) %>%
  purrr::pwalk(fn_draw_pathway) 
  
fn_ai_n <- function(symbol, data){
  print(symbol)
  data %>% 
    dplyr::group_by(pathway, class) %>% 
    dplyr::count() %>% 
    dplyr::ungroup() %>% 
    tidyr::spread(key = class, value = n) %>% 
    dplyr::mutate(Activation = ifelse(rep(tibble::has_name(., "Activation"), 10), Activation, 0)) %>% 
    dplyr::mutate(Inhibition = ifelse(rep(tibble::has_name(., "Inhibition"), 10), Inhibition, 0)) %>% 
    tidyr::replace_na(replace = list(Activation = 0, Inhibition = 0, None = 0)) %>% 
    dplyr::select(pathway, Activation, Inhibition, None)
}
gene_rppa_sig_pval_class %>%
  dplyr::select("symbol", "pathway", "cancer_types", "class", "p.value", "fdr") %>%
  tidyr::nest(-symbol) %>%
  dplyr::mutate(st = purrr::map2(symbol, data, .f = fn_ai_n)) %>%
  tidyr::unnest(st) -> gene_ai_n

gene_ai_n %>% 
  dplyr::mutate(a = Activation / 32, i =  Inhibition / 32 , n = None / 32 ) %>%
  dplyr::select(symbol, pathway, a, i, n) %>%
  tidyr::gather(key = type, value = percent, a, i, n) %>% 
  ggplot(aes(x = factor(1), y = percent, fill = type)) +
  geom_bar(stat = "identity", position = "stack", color = NA) +
  scale_y_continuous(limits = c(0,1)) +
  coord_polar("y") +
  facet_grid(pathway~symbol) +
  # geom_text(aes(x = factor(1), y= .5, label = val_mod, vjust = 4.5)) +
  theme(axis.text=element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        
        strip.text.y = element_text(angle =0,hjust=0,color="black",size=11),
        strip.background = element_blank(),
        strip.text.x = element_text(size=11,angle = 90,vjust = 0),
        
        legend.title= element_blank(), 
        legend.text = element_text(size=14),
        legend.position = "bottom",
        
        panel.background = element_blank(),
        panel.spacing  = unit(0.02, "lines")) +
  scale_fill_manual(label=c("Activation","Inhibition","None"), values=c("red","blue","lightgray")) -> p

# ggsave(filename = "fig_03_e_rppa_all_gene_percentage.pdf", plot = p, device = "pdf", path = rppa_path, width = 30, height = 5)



#--------------------------------------------------------------------
gene_rppa_sig_pval_class %>% 
  dplyr::filter(symbol %in% c("PARK2", "BCL2", "PINK1", "CTSW", "ULK1")) %>% 
  ggplot(aes(x=cancer_types,y=pathway))+
  geom_tile(aes(fill=factor(class)),col="white")+
  scale_fill_manual(limits=c("Activation","Inhibition","None"),values = c("red","blue","lightgray"),na.value="white",labels=c("Activation","Inhibition","None"),name="")+
  
  facet_grid(symbol~.)+
  theme(panel.background=element_rect(colour=NA,fill="white",size=2),
        panel.grid=element_line(colour="white",linetype="dashed"),
        panel.grid.major=element_line(colour="white",linetype="dashed"),
        axis.title=element_blank(),
        axis.text.y=element_text(size=9,colour = "black"),
        axis.text.x=element_text(size=10,colour = "black",angle=90,vjust=0.5,hjust=1),
        axis.ticks=element_line(color="black"),
        legend.text=element_text(size=12),strip.background = element_blank(),strip.text = element_text(size=14),
        legend.title=element_text(size=10),legend.position="bottom",legend.direction="horizontal",
        legend.key.size = unit(0.3, "cm"),
        panel.spacing = unit(0.1, "lines"))



gene_ai_n %>%
  dplyr::arrange(dplyr::desc(Activation), dplyr::desc(Inhibition)) %>%
  dplyr::filter(Activation + Inhibition> 10)

gene_ai_n %>%
  # dplyr::filter(pathway == "Apoptosis") %>%
  dplyr::mutate(a = Activation / 32 * 100, i = - Inhibition / 32 * 100) %>%
  dplyr::select(symbol, pathway, a, i) %>%
  dplyr::filter((a + abs(i)) > (5 / 32 * 100)) -> te
te %>% dplyr::arrange(a+i) %>% .$symbol -> symbol_sort

te %>%
  dplyr::left_join(gene_list, by = "symbol") %>% 
  dplyr::filter(status %in% c("a")) %>% 
  dplyr::rename(pathway = pathway.x) %>% 
  tidyr::gather(key = type, value = percent, c(a,i)) %>%
  tidyr::unite(col = pathway, pathway, type) -> te_plot
te_plot %>% 
  ggplot(aes(x = pathway, y = symbol))+
  geom_tile(aes(fill = percent), col = "white") +
  geom_text(aes(label = ceiling(percent))) +
  scale_fill_gradient2(
    high = "red",
    mid = "white",
    low = "blue"
  ) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1)
    )-> p
ggsave(filename = "rppa_seminar.pdf", plot = p, device = "pdf", path = rppa_path, height = 10, width = 12)

save.image(file = file.path(rppa_path, ".rda_03_e_rppa_gene_expr.rda"))
load(file = file.path(rppa_path, ".rda_03_e_rppa_gene_expr.rda"))


