library(haven)
library(ebal)
library(tidyverse)
library(grf)
library(car) 
library(rmarkdown) 
library(rlang)


## ==================================== ##
## read utils and config for each study ##
## ==================================== ##

source(paste0(ROOT_DIR, "ManyLabs1/utils_ml.R"))
source(paste0(ROOT_DIR, "ManyLabs1/utils_batch_ml.R"))

load(ML_DATA_PATH)
study.names = unique(pdata$original_study)
site.names = unique(pdata$site)
alpha = 0.05
 
## =========================================================== ##
## generalization via worst-case bounds without auxiliary data ##
## =========================================================== ## 
 
load(paste0(ROOT_DIR, "ManyLabs1/generalization/cond_KL_delta_PP.RData"))

# use in-study calibrated KL to compute bounds ("WorstCase" in Figure 7)

KL.max.bounds = data.frame()
for (study.id in 1:15){
  print(paste("Running study", data_id, "..."))
  res = ml_batch_KL_worst_bounds_reweight(study.id, delta = KL.bd$kl[data_id], dr_alg = 'grf', verbose = TRUE)
  res$delta = KL.bd$kl[data_id]
  KL.max.bounds = rbind(KL.max.bounds, res)   
}
 
cond.KL.PIs = KL.max.bounds %>% mutate(type = 'max') 

save(cond.KL.PIs, file = paste0(ROOT_DIR, "ManyLabs1/generalization/cond_KL_PIs_ML1.RData"))
 



 
