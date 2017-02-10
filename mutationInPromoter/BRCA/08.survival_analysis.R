#!usr/bin/Rscript

###################
#Load library
##################
require(methods)
require(tidyverse)
require(stringr)
require(broom)
require(qvalue)
library(survival)
library(survminer)
######################
# Load anno data
#####################
dataRootPath = "/home/cliu18/liucj/projects/1.Mutation_calling_in_non-condig_region_through_EXOME/3.calling/BRCA_reanalysis"

survivalPath <- file.path(dataRootPath, "03.somatic/07.survival")
expression_path <- "/extraspace/TCGA/WXS_RAW/BRCA/RNA-seq"

# Load expression
# For the same sample, average the expression
rnaseq <- read_rds(file.path(expression_path, "brca.expression.tidy.2candidate.rds")) %>%
  mutate(barcode = str_sub(barcode, start = 1, end = 12)) %>%
  filter(type == "tumor") %>%
  distinct() %>%
  group_by(ensid, barcode, type) %>%
  summarise(expression = mean(expression)) %>%
  ungroup()

# Load clinical data
clinical <- read_rds(path = file.path(survivalPath, "BRCA_clinical_clean.rds")) %>%
  rename(time = os_days) %>%
  filter(!is.na(time), time > 0, !is.na(os_status)) %>%
  mutate(status = plyr::revalue(os_status, c(Dead = 1, Alive = 0))) 
  

# prepare for survival
for_survival <-
  rnaseq %>% 
  inner_join(clinical, by = "barcode") %>%
  dplyr::select(ensid, barcode, expression, time, status) %>%
  group_by(ensid) %>%
  mutate(status = parse_number(status),
         group = ifelse(expression < median(expression),"Low", "High")) %>%
  ungroup()
  
survival_coxph_model <-
  for_survival %>%
  group_by(ensid) %>%
  do(tidy(coxph(Surv(time, status) ~ expression, data=., na.action=na.exclude)))

survival_survdiff_model <-
  for_survival %>%
  group_by(ensid) %>%
  do(tidy(survdiff(Surv(time, status) ~ group, data= .)))

survival_survfit_model <-
  for_survival %>%
  group_by(ensid) %>%
  do(tidy(survfit(Surv(time, status) ~ group, 
                  data = ., na.action = na.exclude)))

fit_RRM2B <- survfit(Surv(time, status) ~ group, 
          data = filter(for_survival, ensid == "ENSG00000048392"),
          na.action = na.exclude)
ggsurvplot(fit_RRM2B, pval=T, ggtheme = theme_light(), xlab = "Survival in days", ylab = '', main = paste("Kaplan-Meier Curves RRM2B Expression",sep=""))





fit_TANGO6 <- survfit(Surv(time, status) ~ group, 
               data = filter(for_survival, ensid == "ENSG00000103047"),
               na.action = na.exclude)
ggsurvplot(fit_TANGO6, pval=T, ggtheme = theme_light(), xlab = "Survival in days", ylab = '', main = paste("Kaplan-Meier Curves TANGO6 Expression",sep=""))



# Save workspace
save(list = ls(), file = file.path(geneExpressionDir, "06.geneExpression.RData"))
