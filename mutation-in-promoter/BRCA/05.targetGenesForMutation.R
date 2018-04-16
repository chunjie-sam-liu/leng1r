#!usr/bin/Rscript

###################
#Load library
##################
require(methods)
require(tidyverse)
require(stringr)
require(biomaRt)
######################
# Load anno data
#####################
dataRootPath = "/home/cliu18/liucj/projects/1.Mutation_calling_in_non-condig_region_through_EXOME/3.calling/BRCA_reanalysis"

forAnalysisDir <-
  file.path(dataRootPath, '03.somatic/03.somaticForAnalysis.saveFiles')

mutationDir <-
  file.path(dataRootPath, "03.somatic/04.targetGenesForMutation")

# Load somatic mutation with recur >= 6
somaticMutation <-
  read_rds(path = file.path(forAnalysisDir,
                            'realSomaticMutation.recur5ForAnalysis.rds'))

###################
#Position analysis#
###################
somaticMutationPos <-
  somaticMutation %>%
  group_by(mutation) %>%
  mutate(recurrent = n()) %>%
  ungroup() %>%
  dplyr::select(mutation, feature, ensr, recurrent) %>%
  distinct() %>%
  separate(
    "mutation",
    into = c("chrom", "pos", "ref", "alt"),
    convert = T,
    remove = F
  )

# get nearest genes
GENES = useMart("ENSEMBL_MART_ENSEMBL",
                dataset = "hsapiens_gene_ensembl",
                host = "useast.ensembl.org")
GENES.ATTRIBUTES <- listAttributes(GENES)
GENES.FILTERS <- listFilters(GENES)

# Access from Ensembl
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
    bind_rows(strand.minus)
}

# anno gene function
refDir <- "/home/cliu18/liucj/reference/TSONCOG"
collectGenes <-
  read_tsv(file = file.path(refDir, "collectedOncoAndTumorsupressorGenes.tsv"))

nearestTargetProteinCodingGenes <-
  somaticMutationPos %>%
  apply(1, getAroundGenes, len = "100000000") %>%
  bind_rows() %>%
  dplyr::select(-chromosome_name) %>%
  rename(
    ensg = ensembl_gene_id,
    start = start_position,
    end = end_position,
    symbol = hgnc_symbol,
    biotype = gene_biotype
  ) %>%
  left_join(collectGenes, by = "symbol")

nearestTargetProteinCodingGenes %>%
  write_tsv(path = file.path(mutationDir, "nearestTargetProteinCodingGenes.tsv"))
nearestTargetProteinCodingGenes %>%
  write_rds(path = file.path(mutationDir, "nearestTargetProteinCodingGenes.rds"))


# save workspace
save(
  list = ls(),
  file = file.path(mutationDir, "05.targetGenesForMutation.RData")
)
