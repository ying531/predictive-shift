library(haven)
library(ebal)
library(tidyverse)
library(grf) 
library(rmarkdown) 
library(rlang) 

## ==================================== ##
## read utils and config for each study ##
## ==================================== ##

source(paste0(ROOT_DIR, "ManyLabs1/utils_decomp_ml.R"))
source(paste0(ROOT_DIR, "ManyLabs/utils_batch_ml.R"))  
# ML_DATA_PATH = ...
load(ML_DATA_PATH)

study.names = unique(pdata$original_study)
site.names = unique(pdata$site)

## =========================== ##
## unified processing pipeline ##
## =========================== ##

stab.all = data.frame() 
for (study.id in 1:length(study.names)){
  for (alg in c("ebal", "double")){
    print(paste("running study", study.id, ", ", alg, "..."))
    res = ml_stable_analysis(study.id, decomp_alg = alg, K = 5, center=TRUE)$res
    res$alg = alg
    stab.all = rbind(stab.all, res)
  }
}

save(stab.all, file = paste0(ROOT_DIR, "Manylabs1/predictive/results_stable_ML1.RData"))

