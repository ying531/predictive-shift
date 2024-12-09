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
 
util_path = "../utils_decomp.R"
util_path2 = "../utils_batch.R"
data_path = "/Users/ying/Desktop/Research/Dominik/multi_rep_analysis/Pipeline/clean_data/"
source(util_path)
source(util_path2) 
 

## ========================================== ##
## prediction intervals under covariate shift ##
## ========================================== ##
 
all.trans = data.frame() 
for (data_id in 1:10){   
  print(paste("Running study", data_id, "..."))
  res = compute_transfer(data_id, K=5, method = 'double', center=TRUE) 
  all.trans = rbind(all.trans, res$res)   
}

save(all.trans, file = "results_weighted.RData")
 

## ================================== ##
##    prediction intervals with iid   ##
## ================================== ## 

all.plain = data.frame() 
for (data_id in 1:10){ 
  print(paste("Running study", data_id, "..."))
  res = compute_plain(data_id) 
  all.plain = rbind(all.plain, res$res) 
}

save(all.plain, file = "results_plain_PP.RData") 
 
 
 





 
