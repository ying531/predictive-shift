library(haven)
library(ebal)
library(tidyverse)
library(grf)
library(car) 
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
alpha = 0.05
 
## =========================================================== ##
## generalization via worst-case bounds without auxiliary data ##
## =========================================================== ## 
 
load(paste0(ROOT_DIR, "Pipeline/generalization/cond_KL_delta_PP.RData"))

# use in-study calibrated KL to compute bounds ("WorstCase" in Figure 7)

KL.max.bounds = data.frame()
for (data_id in 1:10){ 
  print(paste("Running study", data_id, "..."))
  res = batch_KL_worst_bounds_reweight(PP_DATA_PATH, data_id, delta = KL.max$kl[data_id])
  res$delta =  KL.max$kl[data_id]
  KL.max.bounds = rbind(KL.max.bounds, res)   
}
 
cond.KL.PIs = KL.max.bounds %>% mutate(type = 'max') 

save(cond.KL.PIs, file = "cond_KL_PIs_PP.RData")
 
 




 
