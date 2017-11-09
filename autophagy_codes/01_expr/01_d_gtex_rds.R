
# library -----------------------------------------------------------------

library(magrittr)

# path --------------------------------------------------------------------

gtex_path <- "/home/liucj/project/projects/6.autophagy/GTEx"


# load data ---------------------------------------------------------------

gtex_gene_tmp <- readr::read_tsv(
  file = file.path(gtex_path, "GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm.gct"),
  skip = 2
)

gtex_gene_tmp %>% 
  readr::write_rds(
    path = file.path(gtex_path, ".rds_gtex_v7_gene_tmp.rds.gz"),
    compress = 'gz'
  )

sample_attributes <- readr::read_tsv(
  file = file.path(gtex_path, "GTEx_v7_Annotations_SampleAttributesDS.txt")
)

sample_attributes %>%  
  readr::write_rds(
    path = file.path(gtex_path, ".rds_gtex_v7_annotations_sample_attributes_ds.rds.gz"),
    compress = "gz"
  )

sample_head <- readxl::read_xlsx(
  path = file.path(gtex_path, "GTEx_Analysis_v7_Annotations_SampleAttributesDD.xlsx")
)

sample_head %>% 
  readr::write_rds(
    path = file.path(gtex_path, ".rds_gtex_v7_annotations_sample_attributes_dd.rds.gz"),
    compress = "gz"
  )

subject_phenotype <- readr::read_tsv(
  file = file.path(gtex_path, "GTEx_v7_Annotations_SubjectPhenotypesDS.txt")
)

subject_phenotype %>% 
  readr::write_rds(
    path = file.path(gtex_path, ".rds_gtex_v7_annotations_subject_phenotyupes_ds.rds.gz"),
    compress = "gz"
  )


# extract data ------------------------------------------------------------

sample_attributes %>% 
  dplyr::select(SAMPID, SMTS, SMTSD, SMAFRZE) ->
  sample_attributes_simple

rm(sample_attributes)
gc()

# Remain samples with attributes.
intersect(
  
  colnames(gtex_gene_tmp)[-c(1,2)],
  sample_attributes_simple %>% dplyr::pull(SAMPID)
) -> sample_barcode

sample_barcode %>% 
  stringr::str_split(pattern = "-", simplify = T) %>% 
  .[, c(1, 2)] %>% 
  as.data.frame() %>% 
  stringr::str_c(.$V1, .$V2)
  

##