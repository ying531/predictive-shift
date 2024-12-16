#!/usr/bin/env Rscript
## Start of problem independent section
library(haven)
library(ebal)
library(tidyverse)
library(grf) 
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

## ==================================== ##
## load data and existing analyses ##
## ==================================== ##
 
load(paste0(ROOT_DIR, "Pipeline/generalization/cond_KL_delta_PP.RData"))
 
study.res = data.frame()

for (seed in 1:10){
  set.seed(seed)
  study.order = sample(1:10, 10)
  
  for (study.id in 1:9){
    exist.id = study.order[1:study.id]
    new.id = study.order[(1+study.id):10]
    
    ###===== generate KL bounds  =====###  
    
    KL_max = max(KL.delta.all %>% filter(data %in% exist.id) %>% select(KL), na.rm=TRUE)  
    
    all_res_max = data.frame() 
    
    for (data_id in new.id){
      data_all = clean_vars(read.csv(paste(PP_DATA_PATH, data_list[data_id], sep = '')),
                            turn_numerics = turn_numeric_list[[data_id]], 
                            othername = (data_id %in% c(5,7,8)))
      
      # generate all pairs of sites
      all_sites = unique(data_all$datacollection)
      all_grid = t(combn(1:length(all_sites),2))
      all_grid = rbind(all_grid, all_grid[,c(2,1)])    
      
      for (ii in 1:nrow(all_grid)){  
        data_org = data_all %>% filter(datacollection == all_sites[all_grid[ii,1]]) 
        data_rep = data_all %>% filter(datacollection == all_sites[all_grid[ii,2]])    
        
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
        if (!is.null(ii.res)){ 
          ii.res$study = paste(all_sites[all_grid[ii,1]], "_to_", all_sites[all_grid[ii,2]], sep='')
          ii.res$study1 = all_sites[all_grid[ii,1]]
          ii.res$study2 = all_sites[all_grid[ii,2]]
          ii.res$size1 = nrow(data_org)
          ii.res$size2 = nrow(data_rep)  
          ii.res$data = data_id
          ii.res$KL = KL_max 
          all_res_max = rbind(all_res_max, ii.res) 
        }
        
      }  
    }
    
    all_res_max = all_res_max %>% mutate(PI_lower_KLmax = lower, PI_upper_KLmax = upper) 
    
    
    
    for (alg_name in c("double", "ebal")){
      
      qt.id = PP.analysis.studies %>% filter(data %in% exist.id, alg == alg_name) %>% 
        summarize( 
          qt_upp_stab = quantile(cond_shift / stab_x_shift, 1-alpha/2),
          qt_low_stab = quantile(cond_shift / stab_x_shift, alpha/2)  
        ) 
      
      PI.remain = PP.merge.res %>% filter(data %in% new.id, alg==alg_name) %>% 
        mutate( 
          PI_low_stab_id = pmin( qt.id$qt_low_stab * stab_x_shift * cond_sensitivity,
                                 qt.id$qt_upp_stab * stab_x_shift * cond_sensitivity),
          PI_upp_stab_id = pmax( qt.id$qt_low_stab * stab_x_shift * cond_sensitivity,
                                 qt.id$qt_upp_stab * stab_x_shift * cond_sensitivity), 
          source_nums = study.id
        ) %>% left_join(
          all_res_max, by = c("study", "data"), suffix = c("", "._KL_max")  
        )
      
      summary.remain = PI.remain %>% group_by(data) %>% 
        summarize( 
          coverage = mean((PI_lower_KLmax <= theta2) * (theta2 <= PI_upper_KLmax), na.rm = TRUE),
          length = mean(PI_upper_KLmax - PI_lower_KLmax, na.rm = TRUE), 
          metric = 'KL'
        ) %>% mutate(source_nums = study.id, alg = alg_name, seed = seed)
      
      study.res = rbind(study.res, summary.remain)
    }  
  } 
  
}


save(study.res, file = paste0(ROOT_DIR, "Pipeline/generalization/KL_calib_hypo.RData")) 

