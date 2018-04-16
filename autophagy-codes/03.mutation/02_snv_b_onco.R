library(magrittr)
tcga_path = "/home/cliu18/liucj/projects/6.autophagy/TCGA"
expr_path <- "/home/cliu18/liucj/projects/6.autophagy/02_autophagy_expr/"
expr_path_a <- file.path(expr_path, "03_a_gene_expr")
snv_path <- "/home/cliu18/liucj/projects/6.autophagy/04_snv"

# load mc3
mc3_file <- file.path(tcga_path, 'syn_mutation_syn7824274_mc3_public.pass.maf')
suppressPackageStartupMessages(require(maftools))
gene_list <- readr::read_rds(file.path(expr_path_a, "rds_03_a_atg_lys_gene_list.rds.gz"))

# mc3_maf <- read.maf(mc3_file)
# readr::write_rds(mc3_maf, path = file.path(tcga_path, "mutation_syn7824274_mc3_public.pass.maf.rds.gz"), compress = 'gz')
readr::read_rds(path = file.path(tcga_path, "mutation_syn7824274_mc3_public.pass.maf.rds.gz")) -> mc3_maf

subsetMaf(mc3_maf, genes = dplyr::filter(gene_list, pathway == "autophagesome formation-core") %>% dplyr::pull(symbol), mafObj = T) -> gene_list_maf

plotmafSummary(gene_list_maf)

oncodrive(maf = gene_list_maf, AACol = "HGVSp", minMut = 5, pvalMethod = 'zscore') -> gene_list_maf_sig
plotOncodrive(res = gene_list_maf_sig, fdrCutOff = 0.5, useFraction = T)
pfamDomains(maf = gene_list_maf, top = 10)
oncoplot(maf = gene_list_maf, removeNonMutated = T)

# getFields(mc3_maf)
# getGeneSummary(mc3_maf)
# pdf(file.path(snv_path, "core_atg_oncoplot.pdf"), width=10, height=10)
oncoplot(maf = mc3_maf, genes = gene_list %>% dplyr::filter(pathway == "autophagesome formation-core") %>% dplyr::pull(symbol), removeNonMutated = TRUE)
# dev.off()

pdf(file.path(snv_path, "lys_oncoplot.pdf"), width=10, height=10)
oncoplot(maf = mc3_maf, genes = gene_list %>% dplyr::filter(status == "l") %>% dplyr::pull(symbol), removeNonMutated = TRUE)
dev.off()

# lollipop 
# fn_loli <- function(.x){
#   tryCatch({
#   path <- file.path(snv_path, 'lollipop')
#   pdf(file.path(path, paste(.x, "lollipop.pdf", sep = ".")), width=12, height=3)
#   lpop = lollipopPlot(maf = mc3_maf, gene = .x,  showMutationRate = F, defaultYaxis = FALSE)
#   dev.off()},
#   error = function(e){1},
#   warning = function(e){1})
#   
# }
# gene_list %>% 
#   dplyr::filter(pathway == "autophagesome formation-core") %>% 
#   dplyr::pull(symbol) %>% 
#   purrr::walk(.f = fn_loli)




save.image(file = file.path(snv_path, ".rda_02_snv_b_onco.rda"))
load(file = file.path(snv_path, ".rda_02_snv_b_onco.rda"))


