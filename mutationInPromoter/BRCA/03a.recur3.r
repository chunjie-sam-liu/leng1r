#!usr/bin/Rscript
library(methods)
library(tidyverse)
library(biomaRt)
library(stringr)



dataRootPath = "/home/cliu18/liucj/projects/1.Mutation_calling_in_non-condig_region_through_EXOME/3.calling/BRCA_reanalysis"

forAnalysisDir <-
  file.path(dataRootPath, '03.somatic/03.somaticForAnalysis.saveFiles')

real_somatic <-
  read_rds(path = file.path(forAnalysisDir, "realSomaticMutation.rds"))


mutationDir <-
  file.path(dataRootPath, "03.somatic/04.targetGenesForMutation")

recur3 <-
  real_somatic %>%
  group_by(mutation, feature) %>%
  mutate(recurrent = n()) %>%
  filter(recurrent >= 3) %>%
  ungroup() %>%
  separate(
    mutation,
    into = c("chrom", "pos", "ref", "alt"),
    remove = F,
    convert = T
  )

write_rds(recur3, path = file.path(forAnalysisDir, "somatic_mutation_recur3_for_analysis.rds"))
write_tsv(recur3, path = file.path(forAnalysisDir, "somatic_mutation_recur3_for_analysis.tsv"))

somatic_pos <-
  recur3 %>%
  dplyr::select(mutation:alt, feature, ensr, recurrent) %>%
  distinct()

# annotation
GENES = useMart("ENSEMBL_MART_ENSEMBL",
                dataset = "hsapiens_gene_ensembl")
GENES.ATTRIBUTES <- listAttributes(GENES)
GENES.FILTERS <- listFilters(GENES)


# get target genes
getAroundGenes <- function(x, len) {
  # The function is used for apply function
  # len set big one
  # get around genes

  around.genes <-
    getBM(
      attributes = c(
        "ensembl_gene_id",
        "chromosome_name",
        "start_position",
        "end_position",
        "strand",
        "hgnc_symbol",
        "entrezgene",
        "gene_biotype"
      ),
      filters = c("chromosome_name",
                  "start",
                  "end",
                  "biotype"),
      values = list(
        chromosome_name = str_replace(x["chrom"], "chr", ""),
        start = as.numeric(x["pos"]) - as.numeric(len),
        end = as.numeric(x["pos"]) + as.numeric(len),
        biotype = "protein_coding"
      ),
      mart = GENES
    )
  # nearest gene on the plus strand
  strand.plus <-
    filter(around.genes,
           strand == 1 &
             start_position - as.numeric(x["pos"]) > 0) %>%
    slice(which.min(start_position - as.numeric(x["pos"]))) %>%
    mutate(
      distance = start_position - as.numeric(x["pos"]),
      mutation = x['mutation'],
      feature = x['feature'],
      ensr = x['ensr'],
      recurrent = x['recurrent']
    )
  # nearest gene on minus strand
  strand.minus <-
    filter(around.genes,
           strand == -1 &
             end_position - as.numeric(x["pos"]) < 0) %>%
    slice(which.max(end_position - as.numeric(x["pos"]))) %>%
    mutate(
      distance = start_position - as.numeric(x["pos"]),
      mutation = x['mutation'],
      feature = x['feature'],
      ensr = x['ensr'],
      recurrent = x['recurrent']
    )
  
  strand.plus %>%
    bind_rows(strand.minus) ->
    result
  result %>% dplyr::select(-chromosome_name)
}

# ts genes
refDir <- "/home/cliu18/liucj/reference/TSONCOG"
collectGenes <-
  read_tsv(file = file.path(refDir, "collectedOncoAndTumorsupressorGenes.tsv"))

nearest_target_coding_genes <-
  somatic_pos %>%
  apply(1, getAroundGenes, len = "100000000") %>%
  bind_rows() %>%
  rename(
    ensg = ensembl_gene_id,
    start = start_position,
    end = end_position,
    symbol = hgnc_symbol,
    biotype = gene_biotype
  ) %>%
  left_join(collectGenes, by = "symbol")

write_tsv(
  nearest_target_coding_genes,
  path = file.path(mutationDir,
                   "nearest_target_coding_genes_3.tsv")
)
write_rds(
  nearest_target_coding_genes,
  path = file.path(mutationDir,
                   "nearest_target_coding_genes_3.rds")
)


# save workspace
save(list = ls(),
     file = file.path(mutationDir,
                      "recur3.rda"))
