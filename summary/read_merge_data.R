library(haven)
library(ebal)
library(tidyverse)
library(grf)
library(car) 
library(rmarkdown) 
library(rlang)
library(formattable)
library(paletteer)
library(ggplot2)
library(ggpubr)
library(scales)
library(ggsci) 



## ==================================== ##
## read utils and analysis results ##
## ==================================== ##


# read results for pipeline

load(paste0(ROOT_DIR, "Pipeline/predictive/results_K5_stable_filtered_centered.RData"))
PP.all.studies = all.studies
load(paste0(ROOT_DIR, "Pipeline/explanatory/results_weighted_PP.RData"))
PP.all.trans = all.trans
load(paste0(ROOT_DIR, "Pipeline/explanatory/results_plain_PP.RData"))
PP.all.plain = all.plain
load(paste0(ROOT_DIR, "Pipeline/generalization/cond_KL_PIs_PP.RData"))
PP.cond.KL.PIs = cond.KL.PIs

# read results for manylabs1

load(paste0(ROOT_DIR, "Manylabs1/explanatory/results_plain_ML1.RData"))
ML.all.plain = all.plain
load(paste0(ROOT_DIR, "Manylabs1/explanatory/results_weighted_ML1.RData"))
ML.all.trans = all.trans
load(paste0(ROOT_DIR, "Manylabs1/predictive/results_stable_ML1.RData"))
ML.all.studies = stab.all
load(paste0(ROOT_DIR, "Manylabs1/generalization/cond_KL_PIs_ML1.RData"))
ML.cond.KL.PIs = cond.KL.PIs


## ==================================== ##
##   merge all results as preparation   ##
## ==================================== ##

# merge analyses of Pipeline

PP.analysis.studies = PP.all.studies %>%  
  left_join(PP.all.trans, by = c("data", "study"), suffix = c("", "._trans")) %>% # merge estimated variance for reweighting transfer
  left_join(PP.all.plain, by = c("data", "study"), suffix = c("", "._plain")) %>% # merge estimated variance for iid transfer
  left_join(PP.cond.KL.PIs %>% filter(type == 'max') %>% # merge KL-worst-case bounds
              mutate(PI_lower_KLmax = lower, PI_upper_KLmax = upper), 
            by = c("data", "study"), suffix = c("", "._KL_max"))  

PP.analysis.studies = PP.analysis.studies %>% # drop NA values
  drop_na(any_of(c('stab_x_shift', 'delta_yx')))

PP.global_qts = PP.analysis.studies %>%  
  group_by(alg) %>% 
  summarize(qt_upp_yx_global = quantile(delta_yx / delta_x, 1-alpha/2, na.rm=TRUE),
            qt_low_yx_global = quantile(delta_yx / delta_x, alpha/2, na.rm=TRUE),
            qt_upp_norm_global = quantile(cond_shift / x_shift, 1-alpha/2, na.rm=TRUE),
            qt_low_norm_global = quantile(cond_shift / x_shift, alpha/2, na.rm=TRUE), 
            qt_upp_stab_global = quantile(cond_shift / stab_x_shift, 1-alpha/2, na.rm=TRUE),
            qt_low_stab_global = quantile(cond_shift / stab_x_shift, alpha/2, na.rm=TRUE)
  ) 

PP.study_qts = PP.analysis.studies %>%  
  group_by(data, alg) %>% 
  summarize(qt_upp_yx_instudy = quantile(delta_yx / delta_x, 1-alpha/2, na.rm=TRUE),
            qt_low_yx_instudy = quantile(delta_yx / delta_x, alpha/2, na.rm=TRUE),
            qt_upp_norm_instudy = quantile(cond_shift / x_shift, 1-alpha/2, na.rm=TRUE),
            qt_low_norm_instudy = quantile(cond_shift / x_shift, alpha/2, na.rm=TRUE),
            qt_upp_stab_instudy = quantile(cond_shift / stab_x_shift, 1-alpha/2, na.rm=TRUE),
            qt_low_stab_instudy = quantile(cond_shift / stab_x_shift, alpha/2, na.rm=TRUE)
  ) 

# merge analyses of ManyLabs1

ML.all.res = ML.all.studies %>% left_join(
  ML.all.trans, by = c("data", "study", 'alg'), suffix = c("", "._trans")
) %>% left_join( 
  ML.all.plain, by = c("data", "study"), suffix = c("", "._plain")
)

ML.analysis.studies = ML.all.res %>% 
  drop_na(any_of(c('stab_x_shift', 'delta_yx'))) %>%
  filter(stab_x_shift <= 10)

ML.global_qts = ML.analysis.studies %>%  
  group_by(alg) %>% 
  summarize(qt_upp_yx_global = quantile(delta_yx / delta_x, 1-alpha/2),
            qt_low_yx_global = quantile(delta_yx / delta_x, alpha/2),
            qt_upp_norm_global = quantile(cond_shift / x_shift, 1-alpha/2),
            qt_low_norm_global = quantile(cond_shift / x_shift, alpha/2), 
            qt_upp_stab_global = quantile(cond_shift / stab_x_shift, 1-alpha/2),
            qt_low_stab_global = quantile(cond_shift / stab_x_shift, alpha/2)) 

ML.study_qts = ML.analysis.studies %>%  
  group_by(data, alg) %>% 
  summarize(qt_upp_yx_instudy = quantile(delta_yx / delta_x, 1-alpha/2),
            qt_low_yx_instudy = quantile(delta_yx / delta_x, alpha/2),
            qt_upp_norm_instudy = quantile(cond_shift / x_shift, 1-alpha/2),
            qt_low_norm_instudy = quantile(cond_shift / x_shift, alpha/2),
            qt_upp_stab_instudy = quantile(cond_shift / stab_x_shift, 1-alpha/2),
            qt_low_stab_instudy = quantile(cond_shift / stab_x_shift, alpha/2)) 



## ==================================== ##
##     construct generalization PIs     ##
## ==================================== ## 

PP.merge.res = PP.analysis.studies %>% 
  left_join(PP.global_qts, by = 'alg') %>% 
  left_join(PP.study_qts, by = c('data', 'alg')) 

PP.PI.all = PP.merge.res %>% mutate(
  PI_low_yx_instudy = theta1 + delta_x + pmin(qt_low_yx_instudy * delta_x, 
                                              qt_upp_yx_instudy * delta_x),
  PI_upp_yx_instudy = theta1 + delta_x + pmax(qt_low_yx_instudy * delta_x, 
                                              qt_upp_yx_instudy * delta_x),
  PI_low_norm_instudy = theta1 + delta_x + pmin(qt_low_norm_instudy * x_shift * cond_sensitivity,
                                                qt_upp_norm_instudy * x_shift * cond_sensitivity),
  PI_upp_norm_instudy = theta1 + delta_x + pmax(qt_low_norm_instudy * x_shift * cond_sensitivity,
                                                qt_upp_norm_instudy * x_shift * cond_sensitivity),
  PI_low_stab_instudy = theta1 + delta_x + pmin(qt_low_stab_instudy * stab_x_shift * cond_sensitivity,
                                                qt_upp_stab_instudy * stab_x_shift * cond_sensitivity),
  PI_upp_stab_instudy = theta1 + delta_x + pmax(qt_low_stab_instudy * stab_x_shift * cond_sensitivity,
                                                qt_upp_stab_instudy * stab_x_shift * cond_sensitivity), 
  PI_low_yx_global = theta1 + delta_x + pmin(as.numeric(qt_low_yx_global) * delta_x, 
                                             as.numeric(qt_upp_yx_global) * delta_x), 
  PI_upp_yx_global = theta1 + delta_x + pmax(as.numeric(qt_low_yx_global) * delta_x, 
                                             as.numeric(qt_upp_yx_global) * delta_x), 
  PI_low_norm_global = theta1 + delta_x + pmin(as.numeric(qt_low_norm_global) * x_shift * cond_sensitivity,
                                               as.numeric(qt_upp_norm_global) * x_shift * cond_sensitivity), 
  PI_upp_norm_global = theta1 + delta_x + pmax(as.numeric(qt_low_norm_global) * x_shift * cond_sensitivity,
                                               as.numeric(qt_upp_norm_global) * x_shift * cond_sensitivity),
  PI_low_stab_global = theta1 + delta_x + pmin(as.numeric(qt_low_stab_global) * stab_x_shift * cond_sensitivity,
                                               as.numeric(qt_upp_stab_global) * stab_x_shift * cond_sensitivity), 
  PI_upp_stab_global = theta1 + delta_x + pmax(as.numeric(qt_low_stab_global) * stab_x_shift * cond_sensitivity,
                                               as.numeric(qt_upp_stab_global) * stab_x_shift * cond_sensitivity), 
  PI_upper_trans = theta1+delta_x + sqrt(hat_var) * qnorm(1-alpha/2),
  PI_lower_trans = theta1+delta_x + sqrt(hat_var) * qnorm(alpha/2),
  PI_upper_plain = theta1 + sqrt(hat_var._plain) * qnorm(1-alpha/2),
  PI_lower_plain = theta1 + sqrt(hat_var._plain) * qnorm(alpha/2)
) 



ML.merge.res = ML.analysis.studies %>% 
  left_join(ML.global_qts, by = 'alg') %>% 
  left_join(ML.study_qts, by = c('data', 'alg')) 

ML.PI.all = ML.merge.res %>% mutate(
  PI_low_yx_instudy = theta1 + delta_x + pmin(qt_low_yx_instudy * delta_x, 
                                              qt_upp_yx_instudy * delta_x),
  PI_upp_yx_instudy = theta1 + delta_x + pmax(qt_low_yx_instudy * delta_x, 
                                              qt_upp_yx_instudy * delta_x),
  PI_low_norm_instudy = theta1 + delta_x + pmin(qt_low_norm_instudy * x_shift * cond_sensitivity,
                                                qt_upp_norm_instudy * x_shift * cond_sensitivity),
  PI_upp_norm_instudy = theta1 + delta_x + pmax(qt_low_norm_instudy * x_shift * cond_sensitivity,
                                                qt_upp_norm_instudy * x_shift * cond_sensitivity),
  PI_low_stab_instudy = theta1 + delta_x + pmin(qt_low_stab_instudy * stab_x_shift * cond_sensitivity,
                                                qt_upp_stab_instudy * stab_x_shift * cond_sensitivity),
  PI_upp_stab_instudy = theta1 + delta_x + pmax(qt_low_stab_instudy * stab_x_shift * cond_sensitivity,
                                                qt_upp_stab_instudy * stab_x_shift * cond_sensitivity), 
  PI_low_yx_global = theta1 + delta_x + pmin(as.numeric(qt_low_yx_global) * delta_x, 
                                             as.numeric(qt_upp_yx_global) * delta_x), 
  PI_upp_yx_global = theta1 + delta_x + pmax(as.numeric(qt_low_yx_global) * delta_x, 
                                             as.numeric(qt_upp_yx_global) * delta_x), 
  PI_low_norm_global = theta1 + delta_x + pmin(as.numeric(qt_low_norm_global) * x_shift * cond_sensitivity,
                                               as.numeric(qt_upp_norm_global) * x_shift * cond_sensitivity), 
  PI_upp_norm_global = theta1 + delta_x + pmax(as.numeric(qt_low_norm_global) * x_shift * cond_sensitivity,
                                               as.numeric(qt_upp_norm_global) * x_shift * cond_sensitivity),
  PI_low_stab_global = theta1 + delta_x + pmin(as.numeric(qt_low_stab_global) * stab_x_shift * cond_sensitivity,
                                               as.numeric(qt_upp_stab_global) * stab_x_shift * cond_sensitivity), 
  PI_upp_stab_global = theta1 + delta_x + pmax(as.numeric(qt_low_stab_global) * stab_x_shift * cond_sensitivity,
                                               as.numeric(qt_upp_stab_global) * stab_x_shift * cond_sensitivity), 
  PI_upper_trans =  theta1 + delta_x + sqrt(hat_var) * qnorm(1-alpha/2),
  PI_lower_trans =  theta1 + delta_x  + sqrt(hat_var) * qnorm(alpha/2),
  PI_upper_plain = theta1 + sqrt(hat_var._plain ) * qnorm(1-alpha/2),
  PI_lower_plain = theta1 + sqrt(hat_var._plain ) * qnorm(alpha/2)
) 
