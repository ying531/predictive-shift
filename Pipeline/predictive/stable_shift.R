#!/usr/bin/env Rscript
 
library(haven)
library(ebal)
library(tidyverse)
library(grf) 
library(rmarkdown) 
library(rlang)

head_list = c('1_BM', '2_CH', '3_BT', '4_BAC', '5_MI', '6_MC', '7_IE', 
              '8_BIH', '9_PG', '10_HS_charity')

data_list = c("1_bigot_misanthrope.csv", "2_cold_heart.csv", "3_bad_tipper.csv", 
              "4_belief_act.csv", "5_moral_inversion.csv", "6_moral_cliff.csv",
              "7_intuitive_economics.csv", "8_burn_in_hell.csv", "9_presumption_guilt.csv",
              "10_higher_standard.csv")

## ==================================== ##
## read utils and config for each study ##
## ==================================== ##

source(paste0(ROOT_DIR, "Pipeline/utils_decomp.R"))
source(paste0(ROOT_DIR, "Pipeline/utils_batch.R"))  
# PP_DATA_PATH = ...
  
## =========================== ##
## unified processing pipeline ##
## =========================== ## 

all.studies = data.frame() 
for (data_id in 1:10){
  for (dec_alg in c('ebal', 'double')){
    print(paste("Running study", data_id, "..."))
    res = stable_analysis(PP_DATA_PATH, data_id, decomp_alg = dec_alg, K=5, center=TRUE)
    res_all = res$all %>% mutate(alg = dec_alg)
    all.studies = rbind(all.studies, res_all) 
  } 
} 
  
file_name = "results_K5_stable_filtered_centered_retry.RData"
save(all.studies, file = paste0(ROOT_DIR, "Pipeline/predictive/", file_name))



 


