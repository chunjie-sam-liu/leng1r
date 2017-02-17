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
gene_expression_dir <- file.path(dataRootPath, "03.somatic/06.geneExpression")
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
for_survival %>% 
  write_tsv(path = file.path(survivalPath, "01.for_survival.tsv"))
for_survival %>% 
  write_rds(path = file.path(survivalPath, "01.for_survival.rds"))


# Survival for expression
# Get whole expression coxp
survival_coxph_model <-
  for_survival %>%
  group_by(ensid) %>%
  do(tidy(coxph(Surv(time, status) ~ expression, data=., na.action=na.exclude))) %>%
  ungroup()
survival_coxph_model %>%
  write_tsv(path = file.path(survivalPath, "survival_coxph_model.tsv"))

# survival_survdiff_model <-
#   for_survival %>%
#   group_by(ensid) %>%
#   survdiff(Surv(time, status) ~ group, data= .)

survival_survfit_model <-
  for_survival %>%
  group_by(ensid) %>%
  do(tidy(survfit(Surv(time, status) ~ group, 
                  data = ., na.action = na.exclude)))

RRM2B <- 
  for_survival %>% 
  filter(ensid == "ENSG00000048392")
fit_RRM2B <- 
  survfit(Surv(time, status) ~ group, 
          data = RRM2B,
          na.action = na.exclude)

RRM2B_survival <-
  ggsurvplot(fit_RRM2B, pval=T, ggtheme = theme_light(), 
           xlab = "Survival in days", 
           ylab = '',
           main = "Kaplan-Meier Curves RRM2B Expression (All)")

print(RRM2B_survival)
ggsave(filename = file.path(survivalPath, "01.RRM2B_survival.png"), device = "png", width = 7, height = 7)

# Chose top 1/4 and bottom 1/4
RRM2B_1_4 <- 
  RRM2B %>% 
  top_n(round(nrow(RRM2B) / 4, 0), expression)
RRM2B_1_4 <-
  RRM2B %>%
  top_n(-round(nrow(RRM2B) / 4, 0), expression) %>%
  bind_rows(RRM2B_1_4)

RRM2B_1_4 %>%
  do(tidy(coxph(Surv(time, status)~ expression, data = .)))

fit_RRM2B_1_4 <- 
  survfit(Surv(time, status) ~ group, 
          data = RRM2B_1_4,
          na.action = na.exclude)
RRM2B_1_4_survival <-
  ggsurvplot(fit_RRM2B_1_4, pval=T, 
             palette = c("#eb1a08", "#3f51b5"),
           xlab = "Overall survival (days)", 
           ylab = 'Surviving',
           legend.labs = c("High (top 1/4)", "Low (bottom 1/4)"))
print(RRM2B_1_4_survival)

RRM2B_1_4_survival %>%
  write_rds(file.path(survivalPath, '02.RRM2B_1_4_survival.rds'))
ggsave(filename = file.path(survivalPath, "02.RRM2B_1_4_survival.png"), device = "png", width = 7, height = 7)

TANGO6 <- filter(for_survival, ensid == "ENSG00000103047")
fit_TANGO6 <- 
  survfit(Surv(time, status) ~ group, 
               data = TANGO6,
               na.action = na.exclude)
TANGO6_survival <- 
  ggsurvplot(fit_TANGO6, pval=T, 
             palette = c("#eb1a08", "#3f51b5"),
             xlab = "Overall survival (days)", 
             ylab = 'Surviving',
             legend.labs = c("High", "Low"))
print(TANGO6_survival)
TANGO6_survival %>%
  write_rds(file.path(survivalPath, '03.TANGO6_survival.rds'))
ggsave(filename = file.path(survivalPath, "03.TANGO6_survival.png"), device = "png", width = 7, height = 7)
# 1/4 TANGO6 data
TANGO6_1_4 <- 
  TANGO6 %>%
  top_n(round(nrow(TANGO6) / 3, 0), expression)
TANGO6_1_4 <- 
  TANGO6 %>%
  top_n(-round(nrow(TANGO6) / 3, 0), expression) %>%
  bind_rows(TANGO6_1_4)

TANGO6_1_4 %>%
  do(tidy(coxph(Surv(time, status)~ expression, data = .)))

fit_TANGO6_1_4 <- 
  survfit(Surv(time, status) ~ group, 
          data = TANGO6_1_4,
          na.action = na.exclude)
TANGO6_1_4_survival <- 
  ggsurvplot(fit_TANGO6_1_4, pval=T, ggtheme = theme_light(), 
           xlab = "Survival in days", 
           ylab = '', 
           main = "Kaplan-Meier Curves TANGO6 Expression (1/3)")

print(TANGO6_1_4_survival)
ggsave(filename = file.path(survivalPath, "04.TANGO6_1_4_survival.png"), device = "png", width = 7, height = 7)

# survival for mutation
# candidate_mutation_with_normal.rds
candidate_mutation_with_normal <- 
  read_rds(file.path(gene_expression_dir, "candidate_mutation_with_normal.rds")) %>%
  filter(type != "NM") %>%
  mutate(barcode = str_sub(barcode, 1, 12))

mutation_for_survival <-
  candidate_mutation_with_normal %>%
  inner_join(for_survival, by = c("barcode", "ensid", "expression"))

# expression of coxp is not high
mutation_for_survival_coxph_model<-
  mutation_for_survival %>%
  group_by(mutation, ensid, symbol) %>%
  do(tidy(coxph(Surv(time, status) ~ type, .)))


RRM2B <- 
  mutation_for_survival %>% 
  filter(ensid == "ENSG00000048392")
fit_RRM2B <- 
  survfit(Surv(time, status) ~ type, 
          data = RRM2B,
          na.action = na.exclude)

RRM2B_survival <-
  ggsurvplot(fit_RRM2B, pval=T, ggtheme = theme_light(), 
             xlab = "Survival in days", 
             ylab = '',
             main = "Kaplan-Meier Curves RRM2B Expression (All)")

print(RRM2B_survival)



TANGO6 <- filter(mutation_for_survival, ensid == "ENSG00000103047")
fit_TANGO6 <- 
  survfit(Surv(time, status) ~ type, 
          data = TANGO6,
          na.action = na.exclude)
TANGO6_survival <- 
  ggsurvplot(fit_TANGO6, pval=T, 
             palette = c("#eb1a08", "#3f51b5"),
             xlab = "Overall survival (days)", 
             ylab = 'Surviving',
             legend.labs = c("High", "Low"))
print(TANGO6_survival)









# Save workspace
save(list = ls(), file = file.path(survivalPath, "08.survival_analysis.RData"))
load(file.path(survivalPath, "08.survival_analysis.RData"))
