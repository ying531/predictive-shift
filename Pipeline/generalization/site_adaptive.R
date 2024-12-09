#!/usr/bin/env Rscript
## Start of problem independent section
args <- commandArgs(trailingOnly = TRUE)
seed <- as.integer(args[1])
step.id <- as.integer(args[2])

# setwd("/Users/ying/Desktop/Stanford/Research/Dominik/github_repo/multi_rep_analysis/Pipeline/analysis")
library(haven)
library(ebal)
library(tidyverse)
library(grf)
# library(car) 
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

util_path = "../../utils_decomp.R"  
util_path2 = "../../utils_batch.R"
data_path = "../../clean_data/"
source(util_path)
source(util_path2)
decomp_alg = 'ebal'
alpha = 0.05

## ==================================== ##
## load data and existing analyses ##
## ==================================== ##
 
load("../Rdata/cond_KL_delta_PP.RData") 
 

## =================================== ##
# ======= going through studies ======= #
## =================================== ##

# cat("Running stable analysis with K = ")
# cat(K)
# cat("and center =")
# cat(if_center)

all_sites = c()
for (data_id in 1:10){
  data_all = clean_vars(read.csv(paste(data_path, data_list[data_id], sep = '')),
                        turn_numerics = turn_numeric_list[[data_id]], 
                        othername = (data_id %in% c(5,7,8))) 
  all_sites = unique(c(all_sites, data_all$datacollection)) 
}
all_sites = sort(all_sites) 

set.seed(seed)
site.res = data.frame()
site.order = sample(1:29, 29)
site.steps = seq(5, 29, by=2)
 
exist.id = all_sites[site.order[1:site.steps[step.id]]]
new.id = all_sites[site.order[(1+site.steps[step.id]):29]]
  
 
 
###===== generate KL bounds  =====###  

KL_max = max(KL.delta.all %>% filter(study1 %in% exist.id, 
                               study2 %in% exist.id) %>% select(KL), na.rm=TRUE) 
KL_qt = quantile(KL.delta.all %>% filter(study1 %in% exist.id, 
                                   study2 %in% exist.id) %>% select(KL), 0.95, na.rm=TRUE)

all_res_max = data.frame()
all_res_qt = data.frame()

for (data_id in 1:10){
  data_all = clean_vars(read.csv(paste(data_path, data_list[data_id], sep = '')),
                        turn_numerics = turn_numeric_list[[data_id]], 
                        othername = (data_id %in% c(5,7,8))) 
  # generate all pairs of sites
  new_sites = intersect(unique(data_all$datacollection), new.id)
  if (length(new_sites) > 3){
    all_grid = t(combn(1:length(new_sites),2))
    all_grid = rbind(all_grid, all_grid[,c(2,1)])    
    
    for (ii in 1:nrow(all_grid)){  
      data_org = data_all %>% filter(datacollection == new_sites[all_grid[ii,1]]) 
      data_rep = data_all %>% filter(datacollection == new_sites[all_grid[ii,2]])    
      
      ###===== pre-process =====###
      
      processed_data = process_before_analysis(data_org, data_rep, data_id)
      data_org = processed_data$data_org
      data_rep = processed_data$data_rep
      filtered_covs = processed_data$filtered_covs
      
      ###===== analysis =====###    
      
      ii.res = tryCatch({data.frame(compute_KL_bound_weighted(data_org, data_rep, delta = KL_max, alg = 'grf',
                                                              covariates = intersect(filtered_covs, KL_vars_list[[data_id]]),  
                                                              analysis_formula = formula_list[[data_id]],
                                                              treatment_variable = treat_var_list[[data_id]]))}, 
                        error = function(e){NULL}) 
      ii.res.qt = tryCatch({data.frame(compute_KL_bound_weighted(data_org, data_rep, delta = KL_qt, alg = 'grf',
                                                                 covariates = intersect(filtered_covs, KL_vars_list[[data_id]]),  
                                                                 analysis_formula = formula_list[[data_id]],
                                                                 treatment_variable = treat_var_list[[data_id]]))}, 
                           error = function(e){NULL}) 
      
      if (!is.null(ii.res)){ 
        ii.res$study = paste(all_sites[all_grid[ii,1]], "_to_", all_sites[all_grid[ii,2]], sep='')
        ii.res$study1 = new_sites[all_grid[ii,1]]
        ii.res$study2 = new_sites[all_grid[ii,2]]
        ii.res$size1 = nrow(data_org)
        ii.res$size2 = nrow(data_rep)  
        ii.res$data = data_id
        ii.res$KL = KL_max
        all_res_max = rbind(all_res_max, ii.res) 
      }
      
      if (!is.null(ii.res.qt)){ 
        ii.res.qt$study = paste(all_sites[all_grid[ii,1]], "_to_", all_sites[all_grid[ii,2]], sep='')
        ii.res.qt$study1 = new_sites[all_grid[ii,1]]
        ii.res.qt$study2 = new_sites[all_grid[ii,2]]
        ii.res.qt$size1 = nrow(data_org)
        ii.res.qt$size2 = nrow(data_rep)  
        ii.res.qt$data = data_id
        ii.res$KL = KL_qt
        all_res_qt = rbind(all_res_qt, ii.res.qt) 
      }
      
    }  
  }
  
}

all_res_max = all_res_max %>% mutate(PI_lower_KLmax = lower, PI_upper_KLmax = upper)
all_res_qt = all_res_qt %>% mutate(PI_lower_KLqt = lower, PI_upper_KLqt = upper)

write.csv(all_res_max, paste0("../output/site_adaptive/res_KLmax_", step.id, "_seed_", seed, ".csv"))
write.csv(all_res_qt, paste0("../output/site_adaptive/res_KLqt_", step.id, "_seed_", seed, ".csv"))

for (alg_name in c("double", "ebal")){
  
  qt.id = analysis.studies %>% filter(study1 %in% exist.id, 
                                      study2 %in% exist.id, alg == alg_name) %>% 
    summarize( 
      qt_upp_stab = quantile(cond_shift / stab_x_shift, 1-alpha/2),
      qt_low_stab = quantile(cond_shift / stab_x_shift, alpha/2)  
    ) 
  
  PI.remain = merge.res %>% filter(study1 %in% new.id, 
                                   study2 %in% new.id, alg==alg_name) %>% 
    mutate( 
      PI_low_stab_id = pmin( qt.id$qt_low_stab * stab_x_shift * cond_sensitivity,
                             qt.id$qt_upp_stab * stab_x_shift * cond_sensitivity),
      PI_upp_stab_id = pmax( qt.id$qt_low_stab * stab_x_shift * cond_sensitivity,
                             qt.id$qt_upp_stab * stab_x_shift * cond_sensitivity), 
      source_nums = site.steps[step.id]
    ) %>% left_join(
      all_res_max, by = c("study", "data"), suffix = c("", "._KL_max")
    ) %>% left_join(
      all_res_qt, by = c("study", "data"), suffix = c("", "._KL_qt")
    )
  
  
  
  summary.remain = rbind( 
    PI.remain %>% group_by(data) %>% 
      summarize( 
        coverage = mean((PI_low_stab_id <= delta_yx) * (delta_yx <= PI_upp_stab_id), na.rm = TRUE),
        length = mean(PI_upp_stab_id - PI_low_stab_id, na.rm = TRUE), 
        metric = 'Ours'
      ),
    PI.remain %>% group_by(data) %>% 
      summarize( 
        coverage = mean((PI_lower_KLmax <= theta2) * (theta2 <= PI_upper_KLmax), na.rm = TRUE),
        length = mean(PI_upper_KLmax - PI_lower_KLmax, na.rm = TRUE), 
        metric = 'KL (max)'
      ),
    PI.remain %>% group_by(data) %>% 
      summarize( 
        coverage = mean((PI_lower_KLqt <= theta2) * (theta2 <= PI_upper_KLqt), na.rm = TRUE),
        length = mean(PI_upper_KLqt - PI_lower_KLqt, na.rm = TRUE), 
        metric = 'KL (qt)'
      ) 
  )  %>% mutate(source_nums = site.steps[step.id], alg = alg_name, seed = seed)
  
  site.res = rbind(site.res, summary.remain)
}  


write.csv(site.res, paste0("../output/site_adaptive/sum_", step.id, "_seed_", seed, ".csv"))


