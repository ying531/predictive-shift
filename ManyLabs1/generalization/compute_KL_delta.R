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

# compute conditional KL bounds 

KL.all = data.frame()
for (study.id in 1:length(study.names)){
  KL.bd = ml_batch_KL_compute_weighted(study.id, dr_alg = 'grf', verbose = FALSE)     
  KL.all = rbind(KL.all, KL.bd)
} 

KL.bd = KL.all %>% group_by(data) %>% 
  summarize(kl = quantile(KL, 0.99, na.rm=TRUE)) # max is too extreme for manylabs1

KL.delta.all = KL.all
save(KL.delta.all, KL.bd, file=paste0(ROOT_DIR, "ManyLabs1/generalization/cond_KL_delta_ML1.RData"))
 





