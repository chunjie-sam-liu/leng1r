library(methods)
library(magrittr)

tcga_path = "/home/cliu18/liucj/projects/6.autophagy/TCGA"
expr_path <- "/home/cliu18/liucj/projects/6.autophagy/02_autophagy_expr/"
expr_path_a <- file.path(expr_path, "03_a_gene_expr")


load(file = file.path(expr_path_a, ".rda_03_h_coca.rda"))

# expr
expr_matrix %>% t() -> expr_matrix_t
factoextra::hcut(expr_matrix_t, k = 4, hc_func = 'hclust', hc_method = 'ward.D', hc_metric = 'pearson', stand = T) -> expr_hcut
cutree(expr_hcut, k = 4) -> expr_group

expr_group %>% 
  tibble::enframe(name = "sample", value = "group") %>% 
  dplyr::inner_join(sample_info, by = "sample") %>% 
  dplyr::distinct(sample, cancer_types, .keep_all = T) %>% 
  ggplot(aes(x = as.factor(group), y = cancer_types)) +
  geom_tile()

fn_encode <- function(.x){
  .d <- tibble::tibble()
  if(.x == 1) {.d <- tibble::tibble(a = 1,b = 0, c = 0, d= 0)}
  if(.x == 2) {.d <- tibble::tibble(a = 0,b = 1, c = 0, d = 0)}
  if(.x == 3) {.d <- tibble::tibble(a = 0,b = 0, c = 1, d = 0)}
  if(.x == 4) {.d <- tibble::tibble(a = 0,b = 0, c = 0, d = 1)}
  .d
}

expr_group %>% 
  tibble::enframe(name = "sample") %>% 
  dplyr::mutate(encode = purrr::map(.x = value, .f = fn_encode)) %>% 
  tidyr::unnest() %>% 
  dplyr::select(-value) -> expr_encode


cnv_matrix %>% t() -> cnv_matrix_t
factoextra::hcut(cnv_matrix_t, k = 4, hc_func = 'hclust', hc_method = 'ward.D', hc_metric = 'pearson', stand = T) -> cnv_hcut
cutree(cnv_hcut, k = 4) -> cnv_group

cnv_group %>% 
  tibble::enframe(name = "sample", value = "group") %>% 
  dplyr::inner_join(sample_info, by = "sample") %>% 
  dplyr::distinct(sample, cancer_types, .keep_all = T) %>% 
  ggplot(aes(x = as.factor(group), y = cancer_types)) +
  geom_tile()

cnv_group %>% 
  tibble::enframe(name = "sample") %>% 
  dplyr::mutate(encode = purrr::map(.x = value, .f = fn_encode)) %>% 
  tidyr::unnest() %>% 
  dplyr::select(-value) -> cnv_encode


methy_matrix %>% t() -> methy_matrix_t
factoextra::hcut(methy_matrix_t, k = 4, hc_func = 'hclust', hc_method = 'ward.D', hc_metric = 'pearson', stand = T) -> methy_hcut
cutree(methy_hcut, k = 4) -> methy_group

methy_group %>% 
  tibble::enframe(name = "sample", value = "group") %>% 
  dplyr::inner_join(sample_info, by = "sample") %>% 
  dplyr::distinct(sample, cancer_types, .keep_all = T) %>% 
  ggplot(aes(x = as.factor(group), y = cancer_types)) +
  geom_tile()

methy_group %>% 
  tibble::enframe(name = "sample") %>% 
  dplyr::mutate(encode = purrr::map(.x = value, .f = fn_encode)) %>% 
  tidyr::unnest() %>% 
  dplyr::select(-value) -> methy_encode

list(expr_encode, cnv_encode, methy_encode) %>% 
  purrr::reduce(.f = function(x, y){x %>% dplyr::inner_join(y, by = "sample")}, .init = tibble::tibble(sample = common_names[-1])) %>% 
  tidyr::gather(key = type, value = v, -sample) %>% 
  tidyr::spread(key = sample, value = v) %>% 
  dplyr::select(-type) -> cc_d


library(ConsensusClusterPlus)
ConsensusClusterPlus(cc_d %>% as.matrix(), maxK=20, reps=500,pItem=0.8,pFeature=1, title="ConsensusClusterplus4", clusterAlg="hc",distance="pearson",seed=1262118388.71279, plot = "pdf") -> cc_res

cc_res %>% readr::write_rds(path = file.path(expr_path_a, ".rds_03_h_coca_cc_cc_res.rds.gz"), compress = "gz")

fn_best_clust <- function(k, d = d){
  cc <- d[[k]][[3]]
  l <- length(cc)
  cm <- d[[k]][[1]]
  tibble::tibble(
    sample = names(cc),
    name = paste("V", c(1:l), sep = ""),
    group = cc) %>% 
    dplyr::arrange(group) -> rank_sample
  
  cm %>% 
    as.data.frame() %>% 
    tibble::as_tibble() %>% 
    tibble::add_column(
      sample1 = paste("V", c(1:l), sep = ""), .before = 1) %>% 
    tidyr::gather(key = sample2, value = sim, -sample1) -> plot_ready
  
  plot_ready %>% 
    ggplot(aes(x= sample1, y = sample2, fill = sim)) +
    geom_tile() +
    scale_x_discrete(limits = rank_sample$name) +
    scale_y_discrete(limits = rank_sample$name) +
    theme(
      axis.text = element_blank(),
      axis.title = element_blank(),
      axis.ticks = element_blank(),
      
      legend.position = "none"
    ) -> p
  ggsave(filename = paste("cc4_c",k, ".tif", sep = ""), plot = p, device = "tiff", width = 8, height = 8, path = "/extraspace/liucj/projects/6.autophagy/02_autophagy_expr/03_a_gene_expr")
}
cc_res -> d

cluster <- multidplyr::create_cluster(19)
tibble::tibble(k = 2:20) %>%
  multidplyr::partition(cluster = cluster) %>%
  multidplyr::cluster_library("magrittr") %>%
  multidplyr::cluster_library("ggplot2") %>%
  multidplyr::cluster_assign_value("fn_best_clust", fn_best_clust)  %>%
  multidplyr::cluster_assign_value("d", d)  %>%
  dplyr::mutate(a = purrr::walk(.x = k, .f = fn_best_clust, d = d)) 
parallel::stopCluster(cluster)

clinical <- readr::read_rds(path = file.path(tcga_path,"pancan34_clinical.rds.gz"))
sample_info <- readr::read_rds(path = file.path(tcga_path, "sample_info.rds.gz"))

d[[4]][[3]] %>% tibble::enframe(name = "sample", value = "group") -> .d3

.d3 %>%
  dplyr::inner_join(sample_info, by = "sample") %>% 
  dplyr::distinct(sample, cancer_types, .keep_all = T) -> .d3_sample

.d3_sample %>% 
  dplyr::group_by(cancer_types, group) %>% 
  dplyr::count() %>% 
  dplyr::ungroup() %>% 
  ggplot(aes(x = as.factor(group), y = cancer_types)) +
  geom_tile() +
  geom_text(aes(label = n))

clinical %>% 
  dplyr::mutate(succinct = purrr::map(.x = clinical, function(.x ){.x %>% dplyr::select(barcode, os_days, os_status)})) %>% 
  dplyr::select(-clinical) %>% 
  tidyr::unnest() %>% 
  dplyr::rename(sample = barcode) %>% 
  dplyr::inner_join(.d3_sample, by = c("cancer_types", "sample")) %>% 
  dplyr::mutate(os_status = dplyr::recode(os_status, "Dead" = 1, "Alive" = 0)) -> .d

.d_diff <- survival::survdiff(survival::Surv(os_days, os_status) ~ group, data = .d)

kmp <- 1 - pchisq(.d_diff$chisq, df = length(levels(as.factor(.d$group))) - 1)
fit_x <- survival::survfit(survival::Surv(os_days, os_status) ~ group, data = .d , na.action = na.exclude)

survminer::ggsurvplot(fit_x, data = .d, pval=T, pval.method = T,
                      xlab = "Survival in days",
                      ylab = 'Probability of survival')



save.image(file = file.path(expr_path_a, ".rda_03_h_coca_cc.rda"))
load(file.path(expr_path_a, ".rda_03_h_coca_cc4.rda"))






