 

##########################################################
#####  batch analysis of stable normalized ratio  ########
########################################################## 

ml_stable_analysis <- function(study.id, decomp_alg = 'ebal', K = 5, center=FALSE){
  # generate all pairs of sites 
  all_grid = t(combn(1:length(site.names),2))
  all_grid = rbind(all_grid, all_grid[,c(2,1)])
  
  all_res = data.frame()
  for (ii in 1:nrow(all_grid)){
    set.seed(1)
    cat(ii, "..")
    data_source = pdata %>% filter(original_study == study.names[study.id],
                           site == site.names[all_grid[ii,1]])
    data_target = pdata %>% filter(original_study == study.names[study.id],
                           site == site.names[all_grid[ii,2]]) 
    
    ###===== pre-process =====###
    covariates_use = intersect(intersect(colnames(data_source), covariates_list), colnames(data_target))
    na.out = 1 * (colSums(is.na(data_source[,covariates_use])) == nrow(data_source)) + 
      1 * (colSums(is.na(data_target[,covariates_use])) == nrow(data_target[,covariates_use]))
    covariates_use = covariates_use[na.out==0]
    data_source = data_source %>% ungroup() %>% select(c("iv", "dv", covariates_use))
    data_target = data_target %>% ungroup() %>% select(c("iv", "dv", covariates_use)) 
    
    ###===== analysis =====###   
    # shortlist: c("resp_gender", "resp_age", "resp_polideo", "RACE_white")
    
    if (center){
      ii.res = tryCatch({data.frame(
        compute_stable(data_source, data_target, 
                       K=K, algorithm = 'grf', estimator = decomp_alg, 
                       covariates = covariates_use, 
                       analysis_formula = "dv ~ iv",
                       treatment_variable = 'iv',
                       center=TRUE, eps_w=0.05)
        )}, error = function(e){NULL})
    }else{
      ii.res = tryCatch({data.frame(compute_stable(data_org, data_rep, 
                                                   K=K, algorithm = 'grf', estimator = decomp_alg, 
                                                   covariates = covariates_use, 
                                                   analysis_formula = "dv ~ iv",
                                                   treatment_variable = 'iv',
                                                   center=FALSE, eps_w=0.05))}, 
                        error = function(e){NULL})
    }
    
    while (is.null(ii.res) & (length(covariates_use) > 1)){
      covariates_use = covariates_use[1:(length(covariates_use)-1)]
      
      if (center){
        ii.res = tryCatch({data.frame(
          compute_stable(data_source, data_target, 
                         K=K, algorithm = 'grf', estimator = decomp_alg, 
                         covariates = covariates_use, 
                         analysis_formula = "dv ~ iv",
                         treatment_variable = 'iv',
                         center=TRUE, eps_w=0.05)
        )}, error = function(e){NULL})
      }else{
        ii.res = tryCatch({data.frame(
          compute_stable(data_source, data_target, 
                         K=K, algorithm = 'grf', estimator = decomp_alg, 
                         covariates = covariates_use, 
                         analysis_formula = "dv ~ iv",
                         treatment_variable = 'iv',
                         center=FALSE, eps_w=0.05)
          )}, error = function(e){NULL})
      }
      
    }
     
    
    if (is.null(ii.res)){
      ii.basic = compute_basic(data_source, data_target, K=K, 
                               covariates = covariates_use,  
                               analysis_formula = "dv ~ iv",
                               treatment_variable = 'iv')
      ii.res = data.frame("theta1"=ii.basic$theta1, 
                          "theta2" =ii.basic$theta2,
                          "delta_total" = ii.basic$delta_total,
                          "delta_yx" = NA, "delta_x" =NA,   
                          "cond_sensitivity"=NA, 
                          "x_sensitivity"  =NA,
                          "cond_sens_naive" = NA, "x_sens_naive" = NA,
                          "cond_shift" =NA,      
                          "x_shift" =NA,         
                          "stab_x_shift" = ii.basic$stab_x_shift,    
                          "infl_ratio" = ii.basic$infl_ratio) 
    } 
    
     
    ii.res$study = paste( all_grid[ii,1] , "_to_", all_grid[ii,2] , sep='')
    ii.res$study1 = all_grid[ii,1]
    ii.res$study2 = all_grid[ii,2]
    ii.res$size1 = nrow(data_source)
    ii.res$size2 = nrow(data_target)   
    all_res = rbind(all_res, ii.res) 
  }
  
  all_res$data = study.id 
  return(list("res" = all_res))
}



ml_stable_analysis_invariant <- function(study.id, decomp_alg = 'ebal', K = 5, center=FALSE){
  # generate all pairs of sites 
  all_grid = t(combn(1:length(site.names),2))
  all_grid = rbind(all_grid, all_grid[,c(2,1)])
  
  all_res = data.frame()
  for (ii in 1:nrow(all_grid)){
    set.seed(1)
    cat(ii, "..")
    data_source = pdata %>% filter(original_study == study.names[study.id],
                                   site == site.names[all_grid[ii,1]])
    data_target = pdata %>% filter(original_study == study.names[study.id],
                                   site == site.names[all_grid[ii,2]]) 
    
    ###===== pre-process =====###
    covariates_use = intersect(intersect(colnames(data_source), covariates_list), colnames(data_target))
    na.out = 1 * (colSums(is.na(data_source[,covariates_use])) == nrow(data_source)) + 
      1 * (colSums(is.na(data_target[,covariates_use])) == nrow(data_target[,covariates_use]))
    covariates_use = covariates_use[na.out==0]
    data_source = data_source %>% ungroup() %>% select(c("iv", "dv", covariates_use))
    data_target = data_target %>% ungroup() %>% select(c("iv", "dv", covariates_use)) 
    
    ###===== analysis =====###   
    # shortlist: c("resp_gender", "resp_age", "resp_polideo", "RACE_white")
    
    if (center){
      ii.res = tryCatch({data.frame(
        compute_stable_invariant(data_source, data_target, 
                       K=K, algorithm = 'grf', estimator = decomp_alg, 
                       covariates = covariates_use, 
                       analysis_formula = "dv ~ iv",
                       treatment_variable = 'iv',
                       center=TRUE, eps_w=0.05)
      )}, error = function(e){NULL})
    }else{
      ii.res = tryCatch({data.frame(compute_stable(data_org, data_rep, 
                                                   K=K, algorithm = 'grf', estimator = decomp_alg, 
                                                   covariates = covariates_use, 
                                                   analysis_formula = "dv ~ iv",
                                                   treatment_variable = 'iv',
                                                   center=FALSE, eps_w=0.05))}, 
                        error = function(e){NULL})
    }
    
    while (is.null(ii.res) & (length(covariates_use) > 1)){
      covariates_use = covariates_use[1:(length(covariates_use)-1)]
      
      if (center){
        ii.res = tryCatch({data.frame(
          compute_stable(data_source, data_target, 
                         K=K, algorithm = 'grf', estimator = decomp_alg, 
                         covariates = covariates_use, 
                         analysis_formula = "dv ~ iv",
                         treatment_variable = 'iv',
                         center=TRUE, eps_w=0.05)
        )}, error = function(e){NULL})
      }else{
        ii.res = tryCatch({data.frame(
          compute_stable(data_source, data_target, 
                         K=K, algorithm = 'grf', estimator = decomp_alg, 
                         covariates = covariates_use, 
                         analysis_formula = "dv ~ iv",
                         treatment_variable = 'iv',
                         center=FALSE, eps_w=0.05)
        )}, error = function(e){NULL})
      }
      
    }
    
    
    if (is.null(ii.res)){
      ii.basic = compute_basic(data_source, data_target, K=K, 
                               covariates = covariates_use,  
                               analysis_formula = "dv ~ iv",
                               treatment_variable = 'iv')
      ii.res = data.frame("theta1"=ii.basic$theta1, 
                          "theta2" =ii.basic$theta2,
                          "delta_total" = ii.basic$delta_total,
                          "delta_yx" = NA, "delta_x" =NA,   
                          "cond_sensitivity"=NA, 
                          "x_sensitivity"  =NA,
                          "cond_sens_naive" = NA, "x_sens_naive" = NA,
                          "cond_shift" =NA,      
                          "x_shift" =NA,         
                          "stab_x_shift" = ii.basic$stab_x_shift,    
                          "infl_ratio" = ii.basic$infl_ratio) 
    } 
    
    
    ii.res$study = paste( all_grid[ii,1] , "_to_", all_grid[ii,2] , sep='')
    ii.res$study1 = all_grid[ii,1]
    ii.res$study2 = all_grid[ii,2]
    ii.res$size1 = nrow(data_source)
    ii.res$size2 = nrow(data_target)   
    all_res = rbind(all_res, ii.res) 
  }
  
  all_res$data = study.id 
  return(list("res" = all_res))
}



##########################################################
#####  batch analysis of stable normalized ratio without reweighting  ########
########################################################## 

diff_ratio_analysis <- function(data_id){
  all_grid = t(combn(1:length(site.names),2))
  all_grid = rbind(all_grid, all_grid[,c(2,1)])
  
  all_res = data.frame()
  for (ii in 1:nrow(all_grid)){
    set.seed(1)
    
    data_source = pdata %>% filter(original_study == study.names[study.id],
                                   site == site.names[all_grid[ii,1]])
    data_target = pdata %>% filter(original_study == study.names[study.id],
                                   site == site.names[all_grid[ii,2]]) 
    
    ###===== pre-process =====###
    covariates_use = intersect(intersect(colnames(data_source), covariates_list), colnames(data_target))
    na.out = 1 * (colSums(is.na(data_source[,covariates_use])) == nrow(data_source)) + 
      1 * (colSums(is.na(data_target[,covariates_use])) == nrow(data_target[,covariates_use]))
    covariates_use = covariates_use[na.out==0]
    data_source = data.frame(data_source %>% select(c("iv", "dv", covariates_use)))
    data_target = data.frame(data_target %>% select(c("iv", "dv", covariates_use)))
    
    ###===== analysis =====###    
      
    ii.res = tryCatch({data.frame(transfer_ratio(data_source, data_target, covariates_use,  
                                                 analysis_formula = "dv ~ iv",
                                                 treatment_variable = "iv"))}, 
                      error = function(e){NULL}) 
    
    if (!is.null(ii.res)){ 
      ii.res$study = paste( all_grid[ii,1] , "_to_", all_grid[ii,2] , sep='')
      ii.res$study1 = all_grid[ii,1] 
      ii.res$study2 = all_grid[ii,2] 
      ii.res$size1 = nrow(data_source)
      ii.res$size2 = nrow(data_target)  
      all_res = rbind(all_res, ii.res)
      # print(ii.res)
    }
    
  } 
  all_res$data = study.id
  
  return(all_res)
}




##################################################### 
#####  batch analysis for worst-case bounds  ########
##################################################### 



ml_batch_KL_compute <- function(study.id, dr_alg = 'grf'){
  all_grid = t(combn(1:length(site.names),2))
  all_grid = rbind(all_grid, all_grid[,c(2,1)])
  
  all_res = data.frame()
  for (ii in 1:nrow(all_grid)){
    set.seed(1)
    
    data_source = pdata %>% filter(original_study == study.names[study.id],
                                   site == site.names[all_grid[ii,1]])
    data_target = pdata %>% filter(original_study == study.names[study.id],
                                   site == site.names[all_grid[ii,2]]) 
    
    ###===== pre-process =====###
    covariates_use = intersect(intersect(colnames(data_source), covariates_list), colnames(data_target))
    na.out = 1 * (colSums(is.na(data_source[,covariates_use])) == nrow(data_source)) + 
      1 * (colSums(is.na(data_target[,covariates_use])) == nrow(data_target[,covariates_use]))
    covariates_use = covariates_use[na.out==0]
    data_source = data.frame(data_source %>% select(c("iv", "dv", covariates_use)))
    data_target = data.frame(data_target %>% select(c("iv", "dv", covariates_use)))
    
    ###===== analysis =====###    
    
    ii.res = tryCatch({compute_KL(data_source, data_target, alg= dr_alg, 
                                  covariates = covariates_use,  
                                  analysis_formula = "dv ~ iv",
                                  treatment_variable = "iv")}, 
                      error = function(e){NULL}) 
    
    if (!is.null(ii.res)){
      ii.res = data.frame("KL" = ii.res)
      ii.res$study = paste( all_grid[ii,1] , "_to_", all_grid[ii,2] , sep='')
      ii.res$study1 = all_grid[ii,1] 
      ii.res$study2 = all_grid[ii,2] 
      ii.res$size1 = nrow(data_source)
      ii.res$size2 = nrow(data_target)  
      all_res = rbind(all_res, ii.res) 
    }
    
  } 
  all_res$data = study.id
  
  return(all_res)
}




batch_KL_worst_bounds <- function(data_id, delta){
  all_grid = t(combn(1:length(site.names),2))
  all_grid = rbind(all_grid, all_grid[,c(2,1)])
  
  all_res = data.frame()
  for (ii in 1:nrow(all_grid)){
    set.seed(1)
    
    data_source = pdata %>% filter(original_study == study.names[study.id],
                                   site == site.names[all_grid[ii,1]])
    data_target = pdata %>% filter(original_study == study.names[study.id],
                                   site == site.names[all_grid[ii,2]]) 
    
    ###===== pre-process =====###
    covariates_use = intersect(intersect(colnames(data_source), covariates_list), colnames(data_target))
    na.out = 1 * (colSums(is.na(data_source[,covariates_use])) == nrow(data_source)) + 
      1 * (colSums(is.na(data_target[,covariates_use])) == nrow(data_target[,covariates_use]))
    covariates_use = covariates_use[na.out==0]
    data_source = data.frame(data_source %>% select(c("iv", "dv", covariates_use)))
    data_target = data.frame(data_target %>% select(c("iv", "dv", covariates_use)))
    
    
    ###===== analysis =====###   
    # compute_KL_bound <- function(df1, df2, delta, 
    #                              covariates, 
    #                              analysis_formula, 
    #                              treatment_variable)
    
    ii.res = tryCatch({data.frame(compute_KL_bound(data_source, data_target, delta = delta,
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
      all_res = rbind(all_res, ii.res) 
    }
    
  } 
  all_res$data = data_id  
  
  return(all_res)
}



ml_batch_KL_compute_weighted <- function(study.id, dr_alg = 'grf', verbose=FALSE){
  all_grid = t(combn(1:length(site.names),2))
  all_grid = rbind(all_grid, all_grid[,c(2,1)])
  
  all_res = data.frame()
  for (ii in 1:nrow(all_grid)){
    set.seed(1)
    
    if (verbose){
      cat(paste0("KL", ii, ".."))
    }
    
    data_source = pdata %>% filter(original_study == study.names[study.id],
                                   site == site.names[all_grid[ii,1]])
    data_target = pdata %>% filter(original_study == study.names[study.id],
                                   site == site.names[all_grid[ii,2]]) 
    
    ###===== pre-process =====###
    covariates_use = intersect(intersect(colnames(data_source), covariates_list), colnames(data_target))
    na.out = 1 * (colSums(is.na(data_source[,covariates_use])) == nrow(data_source)) + 
      1 * (colSums(is.na(data_target[,covariates_use])) == nrow(data_target[,covariates_use]))
    covariates_use = covariates_use[na.out==0]
    data_source = data.frame(data_source %>% ungroup() %>% select(c("iv", "dv", covariates_use)))
    data_target = data.frame(data_target %>% ungroup() %>% select(c("iv", "dv", covariates_use)))
    
    ###===== analysis =====###    
    
    ii.res = tryCatch({compute_KL_weighted(data_source, data_target, alg= dr_alg, 
                                           covariates = covariates_use,  
                                           analysis_formula = "dv ~ iv",
                                           treatment_variable = "iv")}, 
                      error = function(e){NULL}) 
    
    if (!is.null(ii.res)){
      ii.res = data.frame("KL" = ii.res)
      ii.res$study = paste( all_grid[ii,1] , "_to_", all_grid[ii,2] , sep='')
      ii.res$study1 = all_grid[ii,1] 
      ii.res$study2 = all_grid[ii,2] 
      ii.res$size1 = nrow(data_source)
      ii.res$size2 = nrow(data_target)  
      all_res = rbind(all_res, ii.res) 
    }
    
  } 
  all_res$data = study.id
  
  return(all_res)
}


ml_batch_KL_worst_bounds_reweight <- function(study.id, delta, dr_alg = 'grf', verbose=FALSE){
  # read data and clean some variables
  all_grid = t(combn(1:length(site.names),2))
  all_grid = rbind(all_grid, all_grid[,c(2,1)])
  
  all_res = data.frame()
  for (ii in 1:nrow(all_grid)){
    
    set.seed(1) 
    
    if (verbose){
      cat(paste0("KLPI", ii, ".."))
    }
    
    data_source = pdata %>% filter(original_study == study.names[study.id],
                                   site == site.names[all_grid[ii,1]])
    data_target = pdata %>% filter(original_study == study.names[study.id],
                                   site == site.names[all_grid[ii,2]]) 
    
    ###===== pre-process =====###
    covariates_use = intersect(intersect(colnames(data_source), covariates_list), colnames(data_target))
    na.out = 1 * (colSums(is.na(data_source[,covariates_use])) == nrow(data_source)) + 
      1 * (colSums(is.na(data_target[,covariates_use])) == nrow(data_target[,covariates_use]))
    covariates_use = covariates_use[na.out==0]
    data_source = data.frame(data_source %>% ungroup() %>% select(c("iv", "dv", covariates_use)))
    data_target = data.frame(data_target %>% ungroup() %>% select(c("iv", "dv", covariates_use)))
    
    ###===== analysis =====###    
    
    ii.res = tryCatch({data.frame(compute_KL_bound_weighted(data_source, data_target, 
                                                            delta = delta, alg = dr_alg,
                                                            covariates = covariates_use,  
                                                            analysis_formula = "dv ~ iv",
                                                            treatment_variable = 'iv'))}, 
                      error = function(e){NULL}) 
    
    if (!is.null(ii.res)){ 
      ii.res$study = paste( all_grid[ii,1] , "_to_", all_grid[ii,2] , sep='')
      ii.res$study1 = all_grid[ii,1] 
      ii.res$study2 = all_grid[ii,2] 
      ii.res$size1 = nrow(data_source)
      ii.res$size2 = nrow(data_target)  
      ii.res$KL = delta
      all_res = rbind(all_res, ii.res) 
    }
    
  } 
  all_res$data = study.id
  
  return(all_res)
}
 
ml_batch_KL_worst_bounds_reweight_subset <- function(study.id, delta, dr_alg = 'grf', verbose=FALSE, id_set=1:36){
  # read data and clean some variables
  all_grid = t(combn(id_set,2))
  all_grid = rbind(all_grid, all_grid[,c(2,1)])
  
  all_res = data.frame()
  for (ii in 1:nrow(all_grid)){
    
    set.seed(1) 
    
    if (verbose){
      cat(paste0("KLPI", ii, ".."))
    }
    
    data_source = pdata %>% filter(original_study == study.names[study.id],
                                   site == site.names[all_grid[ii,1]])
    data_target = pdata %>% filter(original_study == study.names[study.id],
                                   site == site.names[all_grid[ii,2]]) 
    
    ###===== pre-process =====###
    covariates_use = intersect(intersect(colnames(data_source), covariates_list), colnames(data_target))
    na.out = 1 * (colSums(is.na(data_source[,covariates_use])) == nrow(data_source)) + 
      1 * (colSums(is.na(data_target[,covariates_use])) == nrow(data_target[,covariates_use]))
    covariates_use = covariates_use[na.out==0]
    data_source = data.frame(data_source %>% ungroup() %>% select(c("iv", "dv", covariates_use)))
    data_target = data.frame(data_target %>% ungroup() %>% select(c("iv", "dv", covariates_use)))
    
    ###===== analysis =====###    
    
    ii.res = tryCatch({data.frame(compute_KL_bound_weighted(data_source, data_target, 
                                                            delta = delta, alg = dr_alg,
                                                            covariates = covariates_use,  
                                                            analysis_formula = "dv ~ iv",
                                                            treatment_variable = 'iv'))}, 
                      error = function(e){NULL}) 
    
    if (!is.null(ii.res)){ 
      ii.res$study = paste( all_grid[ii,1] , "_to_", all_grid[ii,2] , sep='')
      ii.res$study1 = all_grid[ii,1] 
      ii.res$study2 = all_grid[ii,2] 
      ii.res$size1 = nrow(data_source)
      ii.res$size2 = nrow(data_target)  
      ii.res$KL = delta
      all_res = rbind(all_res, ii.res) 
    }
    
  } 
  all_res$data = study.id
  
  return(all_res)
}

 
###################################################
##### Transfer estimator + CI under cov shift #####
################################################### 



compute_transfer <- function(study.id, K = 5, method = 'ebal', center=TRUE){
  # read data and clean some variables
  id_data = pdata %>% filter(original_study == study.names[study.id])
  id_sites = unique(id_data$site)
  
  # generate all pairs of sites 
  all_grid = t(combn(1:length(id_sites),2))
  all_grid = rbind(all_grid, all_grid[,c(2,1)])
  
  all_res = data.frame() 
  for (ii in 1:nrow(all_grid)){
    # cat(ii,"..")
    ###===== analysis =====###    
    set.seed(1)
    
    data_source = data.frame(pdata %>% filter(original_study == study.names[study.id],
                                              site == id_sites[all_grid[ii,1]]))
    data_target = data.frame(pdata %>% filter(original_study == study.names[study.id],
                                              site == id_sites[all_grid[ii,2]]))
 
    ii.res = tryCatch({data.frame(transfer_dr(data_source, data_target, K=K, 
                                              algorithm = 'grf', method = method,
                                              covariates = covariates_list,
                                              analysis_formula = "dv ~ iv",
                                              treatment_variable = "iv",
                                              center=center, eps_w=0.05))}, 
                      error = function(e){NULL}) 
    
    if (!is.null(ii.res)){
      ii.res$study = paste(all_grid[ii,1], "_to_", all_grid[ii,2], sep='')
      ii.res$study1 = all_grid[ii,1]
      ii.res$study2 = all_grid[ii,2]
      ii.res$size1 = nrow(data_source)
      ii.res$size2 = nrow(data_target)
      all_res = rbind(all_res, ii.res)
    }else{
      print(paste0("!==", ii, " failed"))
    }
    
  }
  
  all_res$data = study.id
  
  return(list("res" = all_res))
}


###################################################
##### Transfer CI  without adjustment #####
################################################### 

compute_plain <- function(study.id){
  
  id_data = pdata %>% filter(original_study == study.names[study.id])
  id_sites = unique(id_data$site)
  
  # generate all pairs of sites 
  all_grid = t(combn(1:length(id_sites),2))
  all_grid = rbind(all_grid, all_grid[,c(2,1)])
  
  all_res = data.frame()
  for (ii in 1:nrow(all_grid)){
    set.seed(1)
    
    data_source = data.frame(pdata %>% filter(original_study == study.names[study.id],
                                   site == id_sites[all_grid[ii,1]]))
    data_target = data.frame(pdata %>% filter(original_study == study.names[study.id],
                                   site == id_sites[all_grid[ii,2]]))
    
    ###===== pre-process =====###
    covariates_use = intersect(intersect(colnames(data_source), covariates_list), colnames(data_target))
    na.out = 1 * (colSums(is.na(data_source[,covariates_use])) == nrow(data_source)) + 
      1 * (colSums(is.na(data_target[,covariates_use])) == nrow(data_target[,covariates_use]))
    covariates_use = covariates_use[na.out==0]
    data_source = data_source[,c("iv", "dv", covariates_use)]
    data_target = data_target[,c("iv", "dv", covariates_use)]
    
    
    ###===== analysis =====###    
    
    ii.res = tryCatch({data.frame(transfer_plain(data_source, data_target,  
                                              analysis_formula = "dv ~ iv",
                                              treatment_variable = "iv"))}, 
                      error = function(e){NULL}) 
    
    if (!is.null(ii.res)){
      ii.res$study = paste( all_grid[ii,1] , "_to_", all_grid[ii,2] , sep='')
      ii.res$study1 = all_grid[ii,1]
      ii.res$study2 = all_grid[ii,2]
      ii.res$size1 = nrow(data_source)
      ii.res$size2 = nrow(data_target)   
      all_res = rbind(all_res, ii.res)
    }else{
      print(ii)
    }
    
  }
  
  all_res$data = study.id
  
  return(list("res" = all_res))
}
 