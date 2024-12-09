
turn_numeric_list = list(
  # 1_BM
  c('pltclideo', 'gender', 'age', 'pincome', 'age2', 'parented', 'parented2'),
  # 2_CH
  c('pltclideo', 'gender', 'age', 'pincome', 'age2', 'parented', 'parented2'),
  # 3_BT
  c('pltclideo', 'gender', 'age', 'pincome', 'age2', 'parented', 'parented2'),
  # 4_BAC
  c('pltclideo', 'gender', 'age', 'pincome', 'age2', 'parented', 'parented2'),
  # 5_MI
  c('poltclid', 'gender', 'age', 'pincome', 'age2', 'parented', 'parented2'),
  # 6_MC
  c('pltclideo', 'gender', 'age', 'pincome', 'age2', 'parented', 'parented2'),
  # 7_IE 
  c('poltclid', 'gender', 'age', 'pincome', 'age2', 'parented', 'parented2'),
  # 8_BIH
  c('poltclid', 'gender', 'age', 'pincome', 'age2', 'parented', 'parented2'),
  # 9_PG
  c('pltclideo', 'gender', 'age', 'pincome', 'age2', 'parented', 'parented2'),
  # 10_HS_charity
  c('pltclideo', 'gender', 'age', 'pincome', 'age2', 'parented', 'parented2') 
)

sel_var_list = list(
  # 1_BM
  c('condition', 'bigot_personjudge', 'pltclideo', 'gender', 'age', 'pincome', 'pincome2', 'age2', 'parented', 'parented2'),
  # 2_CH
  c('condition', 'tdiff', 'pltclideo', 'gender', 'age', 'pincome', 'pincome2','age2', 'parented', 'parented2'),
  # 3_BT
  c('condition', 'tipper_personjudg', 'pltclideo', 'gender', 'age', 'pincome', 'pincome2','age2', 'parented', 'parented2'),
  # 4_BAC
  c('condition', 'beliefact_mrlblmw_rec', 'beliefact_trustw_rec', 'beliefact_hypocrite', 'pltclideo', 'gender', 'age', 'pincome', 'pincome2','age2', 'parented', 'parented2'),
  # 5_MI
  c('moralgood', 'condition', 'poltclid', 'gender', 'age', 'pincome','pincome2', 'age2', 'parented', 'parented2'),
  # 6_MC
  c('mc_dishonesty', 'mc_ps_dishonesty', 'pltclideo', 'gender', 'age', 'pincome', 'pincome2','age2', 'parented', 'parented2'),
  # 7_IE 
  c('yz', 'condition', 'poltclid', 'gender', "age", 'pincome','pincome2', 'age2', 'parented', 'parented2'),
  # 8_BIH
  c('tdiff', 'poltclid', 'gender', 'age', 'pincome','pincome2', 'age2', 'parented', 'parented2'),
  # 9_PG
  c('companyevaluation', 'condition', 'pltclideo', 'gender', 'age', 'pincome', 'pincome2', 'age2', 'parented', 'parented2'),
  # 10_HS_charity
  c('standard_condition', 'standard_evalu_7items', 'pltclideo', 'gender', 'age', 'pincome', 'pincome2','age2', 'parented', 'parented2') 
)

  

covariates_list = list(
  # 1_BM
  c("pltclideo", "gender", "age", "parented", "age2",  'pincome'),
  # 2_CH
  c("pltclideo", "gender", "age", "parented", "age2", "parented2"),
  # 3_BT
  c("pltclideo", "gender", "age", "parented", "age2", "parented2"),
  # 4_BAC
  c("pltclideo", "gender", "age", "parented", "age2", "parented2"),
  # 5_MI
  c("poltclid", "gender", "age", "parented", "age2", "parented2"),
  # 6_MC
  c("pltclideo", "gender", "age", "parented", "pincome", "age2", "parented2" ),
  # 7_IE 
  c("poltclid", "gender", "age", "pincome", "parented", "age2", "parented2" ),
  # 8_BIH
  c("poltclid", "gender", "age", "parented", "pincome", "age2", "parented2" ),
  # 9_PG
  c("pltclideo", "gender", "age", "parented", "pincome", "age2", "parented2" ),
  # 10_HS_charity
  c("pltclideo", "gender", "age", "parented", "age2", "parented2" ) 
)


KL_vars_list = list(
  # 1_BM
  c("pltclideo", "gender", "age", "parented", 'pincome'),
  # 2_CH
  c("pltclideo", "gender", "age", "parented"),
  # 3_BT
  c("pltclideo", "gender", "age", "parented"),
  # 4_BAC
  c("pltclideo", "gender", "age", "parented"),
  # 5_MI
  c("poltclid", "gender", "age", "parented"),
  # 6_MC
  c("pltclideo", "gender", "age", "parented", "pincome"),
  # 7_IE 
  c("poltclid", "gender", "age", "pincome", "parented"),
  # 8_BIH
  c("poltclid", "gender", "age", "parented", "pincome"),
  # 9_PG
  c("pltclideo", "gender", "age", "parented", "pincome"),
  # 10_HS_charity
  c("pltclideo", "gender", "age", "parented") 
)

formula_list = list(
  # 1_BM
  "bigot_personjudge ~ condition",
  # 2_CH
  "tdiff ~ 1",
  # 3_BT
  "tipper_personjudg ~ condition",
  # 4_BAC
  "beliefact_mrlblmw_rec ~ condition_13",
  # 5_MI
  "moralgood ~ condition",
  # 6_MC
  "diff ~ 1",
  # 7_IE 
  "yz ~ condition",
  # 8_BIH
  "tdiff ~ 1",
  # 9_PG
  "companyevaluation ~ condition",
  # 10_HS_charity
  "standard_evalu_7items ~ condition" 
)

response_var_list = list(
  # 1_BM
  "bigot_personjudge",
  # 2_CH
  "tdiff",
  # 3_BT
  "tipper_personjudg",
  # 4_BAC
  "beliefact_mrlblmw_rec",
  # 5_MI
  "moralgood",
  # 6_MC
  "diff",
  # 7_IE 
  "yz",
  # 8_BIH
  "tdiff",
  # 9_PG
  "companyevaluation",
  # 10_HS_charity
  "standard_evalu_7items" 
)

treat_var_list = list(
  # 1_BM
  "condition",
  # 2_CH
  "(Intercept)",
  # 3_BT
  "condition",
  # 4_BAC
  "condition_13",
  # 5_MI
  "condition",
  # 6_MC
  "(Intercept)",
  # 7_IE 
  "condition",
  # 8_BIH
  "(Intercept)",
  # 9_PG
  "condition",
  # 10_HS_charity
  "condition" 
) 

## ============================================= ##
# unified preprocess of site data before analysis #
## ============================================= ##
# 
process_before_analysis <- function(df_org, df_rep, data_id){

  if (data_id %in% c(2)){
    df_org = df_org %>% mutate(tdiff = cold_moral - cold_traits)
    df_rep = df_rep %>% mutate(tdiff = cold_moral - cold_traits)
  }
  if (data_id %in% c(5)){
    df_org = df_org %>% filter(mi_condition == 1 | mi_condition == 3)
    df_rep = df_rep %>% filter(mi_condition == 1 | mi_condition == 3)
    df_org = df_org %>% mutate(condition = 1 * (mi_condition == 1))
    df_rep = df_rep %>% mutate(condition = 1 * (mi_condition == 1))
  }
  if (data_id %in% c(9)){
    df_org = df_org %>% filter(pg_condition == 1 | pg_condition == 4) %>% mutate(condition = 1 * (pg_condition == 1))
    df_rep = df_rep %>% filter(pg_condition == 1 | pg_condition == 4) %>% mutate(condition = 1 * (pg_condition == 1))
  }
  if (data_id %in% c(7)){
    # generate peusdo outcomes for correlation analysis
    df_org$yz = 0
    df_org$yz[df_org$condition==1] = scale(df_org$htxfair[df_org$condition==1], center=T, scale=T ) * scale(df_org$htxgood[df_org$condition==1], center=T, scale=T )
    df_org$yz[df_org$condition==2] = scale(df_org$htxfair[df_org$condition==2], center=T, scale=T ) * scale(df_org$htxgood[df_org$condition==2], center=T, scale=T )

    df_rep$yz = 0
    df_rep$yz[df_rep$condition==1] = scale(df_rep$htxfair[df_rep$condition==1], center=T, scale=T ) * scale(df_rep$htxgood[df_rep$condition==1], center=T, scale=T )
    df_rep$yz[df_rep$condition==2] = scale(df_rep$htxfair[df_rep$condition==2], center=T, scale=T ) * scale(df_rep$htxgood[df_rep$condition==2], center=T, scale=T )
  }
  if (data_id %in% c(8)){
    # generate peusdo outcomes for correlation analysis
    df_org$tdiff = df_org$BIH_executives - df_org$BIH_vandals
    df_rep$tdiff = df_rep$BIH_executives - df_rep$BIH_vandals
  }

  sel_vars = sel_var_list[[data_id]]

  for (jj in 1:length(sel_vars)){
    if ((sum(is.na(df_org[,sel_vars[jj]])) == nrow(df_org)) || (sum(is.na(df_rep[,sel_vars[jj]])) == nrow(df_rep))){
      next
    }
    df_org[,sel_vars[jj]][is.na(df_org[,sel_vars[jj]])] = median(df_org[,sel_vars[jj]], na.rm=TRUE)
    df_rep[,sel_vars[jj]][is.na(df_rep[,sel_vars[jj]])] = median(df_rep[,sel_vars[jj]], na.rm=TRUE)
  }

  if (data_id %in% c(4)){
    df_org = df_org %>% filter(condition != 2) %>% mutate(condition_13 = 1 * (condition == 1))
    df_rep = df_rep %>% filter(condition != 2) %>% mutate(condition_13 = 1 * (condition == 1))
  }
  if (data_id %in% c(6)){
    df_org$diff = df_org$mc_dishonesty - df_org$mc_ps_dishonesty
    df_rep$diff = df_rep$mc_dishonesty - df_rep$mc_ps_dishonesty
  }
  if (data_id %in% c(10)){
    df_org = df_org %>% mutate(condition = 1 * (standard_condition == 6))
    df_rep = df_rep %>% mutate(condition = 1 * (standard_condition == 6))
  }


  ###===== remove outliers, fill na's =====###
  covariates = covariates_list[[data_id]]
  remove_id = c()
  flag_org = rep(0, nrow(df_org))
  flag_rep = rep(0, nrow(df_rep))
  for (j in 1:length(covariates)){
    # remove variables without any value
    if ((sum(is.na(df_org[,covariates[j]])) == nrow(df_org)) || (sum(is.na(df_rep[,covariates[j]])) == nrow(df_rep))){
      remove_id = c(remove_id, j)
      next
    }else{ # impossible to balance variables (out of support)
      if ((min(df_org[,covariates[j]], na.rm=TRUE) > mean(df_rep[,covariates[j]], na.rm=TRUE)) || (max(df_org[,covariates[j]], na.rm=TRUE) < mean(df_rep[,covariates[j]], na.rm=TRUE))){
        remove_id = c(remove_id, j)
        next
      }

    }
    # remove outliers
    if (covariates[j] %in% c("pincome")){
      flag_org = flag_org + 1 * (df_org[,covariates[j]] >= 3 * sd(df_org[,covariates[j]]) + quantile(df_org[,covariates[j]], max(0.99, 1-1/nrow(df_org)))) +
        1 * (df_org[,covariates[j]] <= -3 * sd(df_org[,covariates[j]]) + quantile(df_org[,covariates[j]], min(0.01, 1/nrow(df_org))))
      flag_rep = flag_rep + 1 * (df_rep[,covariates[j]] >= 3 * sd(df_rep[,covariates[j]]) + quantile(df_rep[,covariates[j]], max(0.99, 1-1/nrow(df_rep)))) +
        1 * (df_rep[,covariates[j]] <= -3 * sd(df_rep[,covariates[j]]) + quantile(df_rep[,covariates[j]], max(0.01, 1/nrow(df_rep))))
    }
    }
  df_org = df_org[flag_org == 0, ]
  df_rep = df_rep[flag_rep == 0, ]

  if (length(remove_id) > 0){
    filtered_covs = covariates[-remove_id]
  }else{
    filtered_covs = covariates
  }

  return(list("data_org" = df_org, "data_rep" = df_rep, "filtered_covs" = filtered_covs))
}

 
## =========================================== ##
##  batch analysis of stable normalized ratio  ## 
## =========================================== ## 

stable_analysis <- function(data_path, data_id, decomp_alg = 'ebal', K = 5, center=FALSE){
  # read data and clean some variables
  data_all = clean_vars(read.csv(paste(data_path, data_list[data_id], sep = '')),
                        turn_numerics = turn_numeric_list[[data_id]], 
                        othername = (data_id %in% c(5,7,8)))
  
  # generate all pairs of sites
  all_sites = unique(data_all$datacollection)
  all_grid = t(combn(1:length(all_sites),2))
  all_grid = rbind(all_grid, all_grid[,c(2,1)])
  
  all_res = data.frame()
  for (ii in 1:nrow(all_grid)){
    set.seed(1)
    data_org = data_all %>% filter(datacollection == all_sites[all_grid[ii,1]]) 
    data_rep = data_all %>% filter(datacollection == all_sites[all_grid[ii,2]])    
    
    processed_data = process_before_analysis(data_org, data_rep, data_id)
    data_org = processed_data$data_org
    data_rep = processed_data$data_rep
    filtered_covs = processed_data$filtered_covs 
    
    ###===== analysis =====###   
    
    if (center){
      ii.res = tryCatch({data.frame(compute_stable(data_org, data_rep, 
                                                   K=K, algorithm = 'grf', estimator = decomp_alg, 
                                                   covariates = filtered_covs, 
                                                   analysis_formula = formula_list[[data_id]],
                                                   treatment_variable = treat_var_list[[data_id]],
                                                   center=TRUE, eps_w=0.05))}, 
                        error = function(e){NULL})
    }else{
      ii.res = tryCatch({data.frame(compute_stable(data_org, data_rep, 
                                                   K=K, algorithm = 'grf', estimator = decomp_alg, 
                                                   covariates = filtered_covs, 
                                                   analysis_formula = formula_list[[data_id]],
                                                   treatment_variable = treat_var_list[[data_id]],
                                                   center=FALSE, eps_w=0.05))}, 
                        error = function(e){NULL})
    }
    
    if (is.null(ii.res)){
      if (center){
        ii.res = tryCatch({data.frame(compute_stable(data_org, data_rep, 
                                                     K=K, algorithm = 'grf', estimator = decomp_alg, 
                                                     covariates = filtered_covs[!(filtered_covs%in% c('pincome', 'pincome2', 'parented2', 'age2'))], 
                                                     analysis_formula = formula_list[[data_id]],
                                                     treatment_variable = treat_var_list[[data_id]],
                                                     center=TRUE, eps_w=0.05))}, 
                          error = function(e){NULL})
      }else{
        ii.res = tryCatch({data.frame(compute_stable(data_org, data_rep, 
                                                     K=K, algorithm = 'grf', estimator = decomp_alg, 
                                                     covariates = filtered_covs[!(filtered_covs%in% c('pincome', 'pincome2', 'parented2', 'age2'))], 
                                                     analysis_formula = formula_list[[data_id]],
                                                     treatment_variable = treat_var_list[[data_id]],
                                                     center=FALSE, eps_w=0.05))}, 
                          error = function(e){NULL})
      }
    }else{
      used_covs_list = list(filtered_covs)
    }
    
    if (is.null(ii.res)){
      ii.basic = compute_basic(data_org, data_rep, K=K, 
                               covariates = filtered_covs,  
                               analysis_formula = formula_list[[data_id]],
                               treatment_variable = treat_var_list[[data_id]])
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
      used_covs_list = NA
    }else{
      used_covs_list = list(filtered_covs[!(filtered_covs%in% c('pincome', 'pincome2'))])
    }
    
     
    ii.res$study = paste(all_sites[all_grid[ii,1]], "_to_", all_sites[all_grid[ii,2]], sep='')
    ii.res$study1 = all_sites[all_grid[ii,1]]
    ii.res$study2 = all_sites[all_grid[ii,2]]
    ii.res$size1 = nrow(data_org)
    ii.res$size2 = nrow(data_rep)
    ii.res$covs = used_covs_list
    cov_ratio = compute_covariate(data_org, data_rep, covariates_list[[data_id]]) 
    ii.res$cov_total = sqrt(mean(cov_ratio^2))
    all_res = rbind(all_res, ii.res) 
  }
  
  all_res$data = data_id
  sub_res = all_res %>% filter(abs(cond_shift) < 100, abs(x_shift) < 100, abs(cond_shift / x_shift) < 10)
  
  return(list("all" = all_res, "sub" = sub_res))
}

  
 

## ============================================ ##
## === batch submission for  
# compute KL divergence E_P[Q/P * log(Q/P)] 
# where P is the !weighted! source distribution
## ============================================ ##

batch_KL_compute_reweight <- function(data_path, data_id, dr_alg = 'grf'){
  # read data and clean some variables
  data_all = clean_vars(read.csv(paste(data_path, data_list[data_id], sep = '')),
                        turn_numerics = turn_numeric_list[[data_id]], 
                        othername = (data_id %in% c(5,7,8)))
  
  # generate all pairs of sites
  all_sites = unique(data_all$datacollection)
  all_grid = t(combn(1:length(all_sites),2))
  all_grid = rbind(all_grid, all_grid[,c(2,1)]) 
  
  all_res = data.frame()
  for (ii in 1:nrow(all_grid)){
    
    # print(paste("running", all_sites[all_grid[ii,1]], "and", all_sites[all_grid[ii,2]]))
    
    data_org = data_all %>% filter(datacollection == all_sites[all_grid[ii,1]]) 
    data_rep = data_all %>% filter(datacollection == all_sites[all_grid[ii,2]])    
    
    ###===== pre-process =====###
    
    processed_data = process_before_analysis(data_org, data_rep, data_id)
    data_org = processed_data$data_org
    data_rep = processed_data$data_rep
    filtered_covs = processed_data$filtered_covs
    
    ###===== analysis =====###    
    
    ii.res = tryCatch({compute_KL_weighted(data_org, data_rep, alg= dr_alg, 
                                  covariates = intersect(filtered_covs, KL_vars_list[[data_id]]),  
                                  analysis_formula = formula_list[[data_id]],
                                  treatment_variable = treat_var_list[[data_id]])}, 
                      error = function(e){NULL}) 
    
    if (!is.null(ii.res)){
      ii.res = data.frame("KL" = ii.res)
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


## ============================================ ##
## === batch submission for  
# compute worst-case intervals under constraint 
# E_P[Q/P * log(Q/P)]  <= delta
# where P is the !weighted! source distribution
## ============================================ ##

batch_KL_worst_bounds_reweight <- function(data_path, data_id, delta, dr_alg = 'grf'){
  # read data and clean some variables
  data_all = clean_vars(read.csv(paste(data_path, data_list[data_id], sep = '')),
                        turn_numerics = turn_numeric_list[[data_id]], 
                        othername = (data_id %in% c(5,7,8)))
  
  # generate all pairs of sites
  all_sites = unique(data_all$datacollection)
  all_grid = t(combn(1:length(all_sites),2))
  all_grid = rbind(all_grid, all_grid[,c(2,1)]) 
  
  all_res = data.frame()
  for (ii in 1:nrow(all_grid)){
    
    cat(paste0(ii, "..."))
    
    data_org = data_all %>% filter(datacollection == all_sites[all_grid[ii,1]]) 
    data_rep = data_all %>% filter(datacollection == all_sites[all_grid[ii,2]])    
    
    ###===== pre-process =====###
    processed_data = process_before_analysis(data_org, data_rep, data_id)
    data_org = processed_data$data_org
    data_rep = processed_data$data_rep
    filtered_covs = processed_data$filtered_covs
    
    ###===== analysis =====###    
    
    ii.res = tryCatch({data.frame(compute_KL_bound_weighted(data_org, data_rep, delta = delta, alg = dr_alg,
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


 
## ============================================ ##
## === batch submission for  
# transfer estimator + PI uncer covariate shift
# method = 'double': doubly robust estimator 
# method = 'ebal': entropy balancing estimator
# K is cross-fitting folder number
## ============================================ ##  

compute_transfer <- function(data_id, K = 5, method = 'double', center=TRUE){
  # read data and clean some variables
  data_all = clean_vars(read.csv(paste(data_path, data_list[data_id], sep = '')),
                        turn_numerics = turn_numeric_list[[data_id]], 
                        othername = (data_id %in% c(5,7,8)))
  
  # generate all pairs of sites
  all_sites = unique(data_all$datacollection)
  all_grid = t(combn(1:length(all_sites),2))
  all_grid = rbind(all_grid, all_grid[,c(2,1)])
  
  all_res = data.frame()
  for (ii in 1:nrow(all_grid)){
    data_org = data_all %>% filter(datacollection == all_sites[all_grid[ii,1]]) 
    data_rep = data_all %>% filter(datacollection == all_sites[all_grid[ii,2]])   
    
    ###===== pre-process =====###
    processed_data = process_before_analysis(data_org, data_rep, data_id) 
    
    ###===== analysis =====###    
 
    ii.res = tryCatch({data.frame(transfer_dr(processed_data$data_org, 
                                              processed_data$data_rep, K=K, 
                                              algorithm = 'grf', method = method,
                                              covariates = processed_data$filtered_covs,
                                              analysis_formula = formula_list[[data_id]],
                                              treatment_variable = treat_var_list[[data_id]],
                                              center=center, eps_w=0.05))}, 
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
  
  return(list("res" = all_res))
}



## ============================================ ##
## === batch submission for  
# transfer estimator + PI uncer iid assumption 
## ============================================ ##  

compute_plain <- function(data_id){
  # read data and clean some variables
  data_all = clean_vars(read.csv(paste(data_path, data_list[data_id], sep = '')),
                        turn_numerics = turn_numeric_list[[data_id]], 
                        othername = (data_id %in% c(5,7,8)))
  
  # generate all pairs of sites
  all_sites = unique(data_all$datacollection)
  all_grid = t(combn(1:length(all_sites),2))
  all_grid = rbind(all_grid, all_grid[,c(2,1)])
  
  all_res = data.frame()
  for (ii in 1:nrow(all_grid)){
    
    data_org = data_all %>% filter(datacollection == all_sites[all_grid[ii,1]]) 
    data_rep = data_all %>% filter(datacollection == all_sites[all_grid[ii,2]])    
    
    ###===== pre-process =====###
    processed_data = process_before_analysis(data_org, data_rep, data_id)
    data_org = processed_data$data_org
    data_rep = processed_data$data_rep
    filtered_covs = processed_data$filtered_covs
    
    ###===== analysis =====###    
    
    ii.res = tryCatch({data.frame(transfer_plain(data_org, data_rep,  
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
  
  return(list("res" = all_res))
}
 
