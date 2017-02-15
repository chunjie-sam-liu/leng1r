#!usr/bin/Rscript

###################
#Load library
##################
require(methods)
require(tidyverse)
require(stringr)

######################
# Load anno data
#####################
data_root_path = "/home/cliu18/liucj/projects/1.Mutation_calling_in_non-condig_region_through_EXOME/3.calling/BRCA_reanalysis"

binding_dir = file.path(data_root_path, '03.somatic/09.binding_region')

#Load mutation binding 
# the data was extract by the  find_motif_affect_by_nutation.sh
mutation_affect_binding <-
  read_tsv(file.path(binding_dir, 'mutation_affect_binding.tsv'),
           col_names = c("chrom", "start", "end", "TF", "mutation"))

# simplify 240 entries.
# merge same TF find the common region the one TF
common_region <- 
  mutation_affect_binding %>% 
  group_by(mutation, TF, chrom) %>% 
  summarise(start = max(start), end = min(end))

suppressMessages(require(Gviz))
suppressMessages(require(GenomicRanges))
suppressMessages(require(biomaRt))
GENES = useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host = "useast.ensembl.org")
GENES.ATTRIBUTES <- listAttributes(GENES)
GENES.FILTERS <- listFilters(GENES)

