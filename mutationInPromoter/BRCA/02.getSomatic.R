#!usr/bin/Rscript

###################
#Load library
##################
require(methods)
require(tidyverse)
require(stringr)
require(VariantAnnotation) # for load vcf file

##################### 
#Load manifest
#####################
downloadPath <- '/extraspace/TCGA/WXS_RAW/BRCA/downloadDataList'

dataRootPath = "/home/cliu18/liucj/projects/1.Mutation_calling_in_non-condig_region_through_EXOME/3.calling/BRCA_reanalysis"

rawsnpPath = file.path(dataRootPath,"01.rawsnp")
filtersnpPath = file.path(dataRootPath, "02.filtersnp")
somaticAnnoPath = file.path(dataRootPath, "03.somatic/01.annotation")

typeName <- c("Primary Tumor" = "tumor", "Blood Derived Normal" = "normal", "Solid Tissue Normal" = "normal", "Metastatic" = "tumor")

manifest <- 
  read_tsv(file = file.path(downloadPath, "manifestData.info")) %>%
  mutate(type = plyr::revalue(type, typeName), 
         rawsnppath = file.path(rawsnpPath,paste(bam, "SNP.vcf", sep=".")))

#################################
#Load raw VCF and filter somatic#
#################################

vcf2tibble <- function(x){

  vcf = x['rawsnppath']
  barcode = x['barcode']
  bam = x['bam']
  stype = x['type']
  
  tmp <- readVcf(vcf, genome = "GRCh38")
  
  # Convert S4Vector into data.frame
  tmpAD <- as.data.frame(geno(tmp)$AD) %>%
    mutate(mutation = rownames(.))
  colnames(tmpAD) <- c('AD',"mutation")
  
  # Get mapping quality
  tmpMQ <- as.data.frame(info(tmp)['MQ']) %>%
    mutate(mutation = rownames(.)) %>%
    tbl_df()
  
  # combine mapping and AD
  tmpCombined <-
    tmpAD %>% 
    mutate(AD = as.character(AD)) %>%
    tbl_df() %>%
    mutate(AD = str_replace_all(AD, "[c\\(|\\)]", "")) %>%
    filter(str_count(AD,",|:") == 1) %>% # remove nAD = 3
    separate(AD, c("refDepth", "altDepth"), sep = "[\\,\\:]", convert = T) %>%
    mutate(barcode = barcode, Depth = refDepth + altDepth, bam = bam) %>%
    dplyr::select(barcode, bam, mutation, Depth, refDepth, altDepth) %>%
    inner_join(tmpMQ, by = "mutation")
  
  if(stype == "tumor"){
    tmpCombined <-
      tmpCombined %>%
      filter(MQ >= 20, Depth >= 10, altDepth >= 3)
  }
  else{
    tmpCombined <-
      tmpCombined %>%
      filter(MQ >= 20, Depth >= 5, altDepth >= 3)
  }
  # write to path
  write_tsv(tmpCombined,
            file.path(filtersnpPath, paste(bam, "SNP.vcf.filter", sep = ".")
                      ))
  return(tmpCombined)
}
  
loadVcf <- function(.data){
  .data %>%
    apply(1, vcf2tibble) %>%
    bind_rows()
} 

somaticMutation <- function(mani){
  cases <- 
    unique(mani$case)
  numStat <- 
    tibble(barcode = character(), somatic = numeric())
  
  for(i in cases){
    normal <- 
      mani %>%
      filter(case == i, type == "normal") %>%
      loadVcf()
    
    tumor <-
      mani %>%
      filter(case == i, type == "tumor") %>%
      loadVcf() 
    
    somatic <- 
      tumor %>%
      anti_join(normal, by = "mutation")
    
    # Stat number
    numStat <- 
      somatic %>%
      group_by(barcode) %>%
      count() %>%
      dplyr::rename(somatic = n) %>%
      bind_rows(numStat)
      
    # for annotation
    barcode <- unique(somatic$barcode)[1]
    somatic %>% 
      separate(mutation, c("chrom", "pos", "ref", "alt"), convert = T) %>%
      mutate(end = pos) %>%
      dplyr::select(chrom, pos, end, ref, alt, 
                    Depth, refDepth, altDepth, MQ, barcode) %>% 
      write_tsv(path = file.path(somaticAnnoPath,
                                 paste(barcode,"avinput", sep = ".")),
                col_names = F)
  }
  return(numStat)
}

manifest %>% 
  somaticMutation() %>%
  write_tsv(path = file.path(dataRootPath, 'somaticStat.tsv'))
  

