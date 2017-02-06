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

# Load mutation and target genes
targetGeneDir <- file.path(dataRootPath, "03.somatic/04.targetGenesForMutation")
mutationTarget <- 
  read_rds(path = file.path(targetGeneDir, "nearestTargetProteinCodingGenes.rds"))



geneExpressionDir <- file.path(dataRootPath, "03.somatic/06.geneExpression")





# Save workspace
save(list = ls(), file = file.path(geneExpressionDir, "06.geneExpression.RData"))