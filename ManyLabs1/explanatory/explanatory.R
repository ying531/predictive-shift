library(haven)
library(ebal)
library(tidyverse)
library(grf)
library(car) 
library(rmarkdown) 
library(rlang)


analysis_formula = "dv ~ iv"
treatment_variable = "iv"

covariates_list = c(
  "resp_gender", "resp_age", "RACE_white", "resp_polideo", "resp_american",
  "resp_american_pid", "resp_american_ideo",  "resp_ethnicity_hisp" 
) 

## ==================================== ##
## read utils and config for each study ##
## ==================================== ##
 

source("../utils_ml.R")
source("../utils_batch_ml.R") 
load("../data/Manylabs1_data.RData") # load data 

study.names = unique(pdata$original_study)
site.names = unique(pdata$site)

options(warn=-1)
alpha = 0.05

 

## ========================================== ##
## prediction intervals under covariate shift ##
## ========================================== ##

all.trans = data.frame() 
for (study.id in 1:length(study.names)){  
  for (method in c('double', 'ebal')){
    print(paste("Running study", study.id, "..."))
    res = compute_transfer(study.id, K=5, method = method, center=TRUE) 
    all.trans = rbind(all.trans, res$res %>% mutate(alg = method))   
  } 
}

save(all.trans, file = "results_weighted_ML1.RData")


## ================================== ##
##    prediction intervals with iid   ##
## ================================== ## 
all.plain = data.frame() 
for (study.id in 1:length(study.names)){ 
  print(paste("Running study", study.id, "..."))
  res = compute_plain(study.id) 
  all.plain = rbind(all.plain, res$res) 
}

save(all.plain, file = "results_plain_ML1.RData") 
 





 
