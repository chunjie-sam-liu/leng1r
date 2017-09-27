
# library -----------------------------------------------------------------

library(magrittr)
library(ggplot2)

# load path ---------------------------------------------------------------

tcga_path <- "/home/cliu18/liucj/projects/6.autophagy/TCGA"
afhl_path <- "/home/cliu18/liucj/projects/6.autophagy/09_afh_vs_afl"
afhl_class <- file.path(afhl_path, "01_af_h_l_classification")
clinical_path <- file.path(afhl_path, "02_clinical")

# clinical data and classify ----------------------------------------------
clinical <- readr::read_rds(file.path(tcga_path, "pancan34_clinical.rds.gz"))
p62_sample_classification <- readr::read_rds(path = file.path(afhl_class, '.rds_01_p62_sample_classification.rds.gz')) %>% 
  dplyr::select(-barcode) %>% 
  tidyr::nest(-cancer_types)

clinical %>% 
  dplyr::mutate(
    clinical = purrr::map(
      .x = clinical,
      .f = function(.x){
        .x %>% 
          dplyr::select(barcode, time = os_days, status = os_status) %>% 
          dplyr::rename(sample = barcode) %>% 
          tidyr::drop_na() %>% 
          dplyr::filter(time > 0) %>% 
          dplyr::mutate(status = plyr::revalue(status, c("Dead" = 1, "Alive" = 0))) %>% 
          dplyr::mutate(status = as.integer(status))
      }
    )
  ) %>%
  dplyr::inner_join(p62_sample_classification, by = "cancer_types") %>% 
  dplyr::mutate(
    merge = purrr::map2(
      .x = clinical,
      .y = data,
      .f = function(.x, .y) {
        .x %>% dplyr::inner_join(.y, by = "sample") %>% 
          dplyr::rename(group = type)
      }
    )
  ) %>% 
  dplyr::select(cancer_types, merge) %>% 
  tidyr::unnest() -> p62_clinical

p62_clinical %>% 
  dplyr::group_by(cancer_types) %>% 
  dplyr::filter(n() > 20) %>% 
  tidyr::nest(-cancer_types) -> p62_clinical_nest

p62_clinical_nest %>% 
  dplyr::mutate(hazard_ratio = purrr::map2(
    .x = data,
    .y = cancer_types,
    .f = function(.x, .y){
      .x %>% 
        dplyr::mutate(time = time / 30) %>% 
        dplyr::mutate(time = ifelse(time > 60, 60, time)) -> .dd
      
      .dd  -> .d
      
      survival::coxph(survival::Surv(time, status) ~ rppa, data = .dd) %>% 
        broom::tidy() %>% 
        dplyr::mutate(hazard_ratio = exp(estimate)) %>% 
        dplyr::select(hazard_ratio, coxp = p.value) -> .hazard_coxp
      
      
      survival::survdiff(survival::Surv(time, status) ~ group, data = .d) -> .d_diff
      kmp <- 1 - pchisq(.d_diff$chisq, df = length(levels(as.factor(.d$group))) - 1)
      
      survival::survfit(survival::Surv(time, status) ~ group, data = .d , na.action = na.exclude) -> fit_x
      
      CPCOLS <- c("#CD0000", "#000080")
      survminer::ggsurvplot(fit_x, 
                            data = .d, 
                            xlab = "Survival in months",
                            ylab = 'Probability of survival',
                            palette = CPCOLS,
                            ggtheme = theme_bw()) -> .p
      
      .label <- glue::glue("{.y}
                           p = {signif(kmp, 2)}")
      
      tibble::tibble(
        hazard_ratio = .hazard_coxp$hazard_ratio,
        coxp = .hazard_coxp$coxp,
        kmp = kmp,
        p = list(.p$plot + annotate("text", x = 5, y = 0.25, label = .label) + theme(legend.position = ""))
      )
    }
  )) %>% 
  dplyr::select(-data) %>% 
  tidyr::unnest() -> p62_survival
p62_survival %>% dplyr::filter(coxp < 0.05, kmp > 0.05) -> p62_survival_median

# use top 40% and bottom 40%
p62_clinical_nest %>% 
  dplyr::mutate(hazard_ratio = purrr::map2(
    .x = data,
    .y = cancer_types,
    .f = function(.x, .y){
      .x %>% 
        dplyr::mutate(time = time / 30) %>% 
        dplyr::mutate(time = ifelse(time > 60, 60, time)) %>% 
        dplyr::mutate(
          group = dplyr::case_when(
            rppa > quantile(rppa, 0.6) ~ "H",
            rppa < quantile(rppa, 0.4) ~ "L",
            TRUE ~ "M"
          )
        ) -> .dd
      
      .dd %>% dplyr::filter(group != "M") -> .d
      
      survival::coxph(survival::Surv(time, status) ~ rppa, data = .dd) %>% 
        broom::tidy() %>% 
        dplyr::mutate(hazard_ratio = exp(estimate)) %>% 
        dplyr::select(hazard_ratio, coxp = p.value) -> .hazard_coxp
      
      
      survival::survdiff(survival::Surv(time, status) ~ group, data = .d) -> .d_diff
      kmp <- 1 - pchisq(.d_diff$chisq, df = length(levels(as.factor(.d$group))) - 1)
      
      survival::survfit(survival::Surv(time, status) ~ group, data = .d , na.action = na.exclude) -> fit_x
      
      CPCOLS <- c("#CD0000", "#000080")
      survminer::ggsurvplot(fit_x, 
                            data = .d, 
                            xlab = "Survival in months",
                            ylab = 'Probability of survival',
                            palette = CPCOLS,
                            ggtheme = theme_bw()) -> .p
      
      .label <- glue::glue("{.y}
                           p = {signif(kmp, 2)}")
      
      tibble::tibble(
        hazard_ratio = .hazard_coxp$hazard_ratio,
        coxp = .hazard_coxp$coxp,
        kmp = kmp,
        p = list(.p$plot + annotate("text", x = 5, y = 0.25, label = .label) + theme(legend.position = "", panel.grid = element_blank()))
      )
    }
  )) %>% 
  dplyr::select(-data) %>% 
  tidyr::unnest() %>% 
  dplyr::filter(coxp < 0.05, kmp < 0.05) -> p62_survival_quantile40

gridExtra::arrangeGrob(grobs = p62_survival_quantile40$p, nrow = 2) %>% 
  ggsave(filename = "fig_01_p62_survival.pdf",
         plot = .,
         device = "pdf", 
         width = 8,
         height = 6,
         path = clinical_path)


# p62 subtype -------------------------------------------------------------



