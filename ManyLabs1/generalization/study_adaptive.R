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

load(paste0(ROOT_DIR, "ManyLabs1/generalization/cond_KL_delta_ML1.RData"))
KL.delta.99 = matrix(0, nrow = 10, ncol = 14) 
study.orders = matrix(0, nrow = 10, ncol = 15)

for (ss in 1:10){
  set.seed(ss) 
  study.orders[ss,] = sample(1:15, 15) 
  
  for (sid in 1:14){   
    sub.delta = KL.delta.all %>% filter(data %in% study.orders[ss, 1:sid], KL != Inf) %>% select(KL)
    KL.delta.99[ss, sid] = as.numeric(quantile(sub.delta, 0.99, na.rm=TRUE)) 
  }
}

ML.KL.study.all = data.frame()
for (seed in 1:10){
  for (step.id in 1:14){
    # determine what KL bound to use 
    kl.delta = KL.delta.99[seed, step.id] 
    for (study.id in (step.id+1):15){
      # the actual data id 
      data_id = study.orders[seed, study.id]
      # generate KL bounds
      all_res = ml_batch_KL_worst_bounds_reweight(study.id = data_id, delta = kl.delta, dr_alg = 'grf', verbose = FALSE) 
      all_res$actual_data_id = data_id 
      all_res$seed = seed
      all_res$num_exist = step.id
      all_res$study_order = study.id
      
      ML.KL.study.all = rbind(ML.KL.study.all, all_res)
    }
  }
}

save(ML.KL.study.all, file=paste0(ROOT_DIR, "ManyLabs1/generalization/KL_calib_study_ML1.RData"))
  
