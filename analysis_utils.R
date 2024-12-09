# created Oct 5, 2024
# utils for plots for both pipeline and manylabs1



shift_trans = function(d = 0) {
  scales::trans_new("shift", transform = function(x) x - d, inverse = function(x) x + d)
}
 

clean_vars <- function(data, turn_numerics, othername=FALSE){
  data_all = data %>% mutate(age = clean_age(yearbirth)) %>% filter(age < 100)
  if (othername){
    data_all$pincome = suppressWarnings(as.numeric(gsub(">", "", gsub("\\+", "", gsub(",", "", gsub("\\$", "", data_all$familyinc)))))) / 10^4 
    # data_all$pincome[is.na(data_all$pincome)] = median(data_all$pincome, na.rm=TRUE)
  }else{
    data_all$pincome = suppressWarnings(as.numeric(gsub(">", "", gsub("\\+", "", gsub(",", "", gsub("\\$", "", data_all$faminc)))))) / 10^4
    # data_all$pincome[is.na(data_all$pincome)] = median(data_all$pincome, na.rm=TRUE)
  }
  # clean up extreme values of pincome 
  data_all$pincome = pmin(data_all$pincome, 1000)
  
  data_all = data_all %>% mutate(gender = gender - 1, pincome2 = pincome^2, parented2 = parented^2, age2 = age^2)
  
  data_all[turn_numerics] = sapply(data_all[turn_numerics], as.numeric)
  
  return(data_all)
}

clean_age <- function(yearbirth){ 
  age = yearbirth
  age[yearbirth>=1000 & yearbirth <= 2014 & !is.na(yearbirth)] = 2014 - yearbirth[yearbirth>=1000 & yearbirth <= 2014 & !is.na(yearbirth)]
  age[age >= 80 & age <= 99 & !is.na(age)] = 114 - age[age >= 80 & age <= 99 & !is.na(age)]
  return(age)
} 

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


 

head_list = c('1_BM', '2_CH', '3_BT', '4_BAC', '5_MI', '6_MC', '7_IE', 
              '8_BIH', '9_PG', '10_HS_charity')

data_list = c("1_bigot_misanthrope.csv", "2_cold_heart.csv", "3_bad_tipper.csv", 
              "4_belief_act.csv", "5_moral_inversion.csv", "6_moral_cliff.csv",
              "7_intuitive_economics.csv", "8_burn_in_hell.csv", "9_presumption_guilt.csv",
              "10_higher_standard.csv") 



generate_ratio_plot <- function(ratio.data, alg_name = 'ebal'){
  to.plot = rbind(ratio.data %>% 
                    mutate(discrepancy = abs(delta_x), name = 'covariate', 
                           target = abs(delta_yx), target_name = 'conditional', 
                           type = "Non-standardized") %>% 
                    select(discrepancy, name, target, target_name, type, alg, study1, study2, size1, size2), 
                  ratio.data %>% 
                    mutate(discrepancy = abs(delta_yx), name = 'conditional',
                           target = abs(delta_x), target_name = 'covariate', 
                           type = "Non-standardized") %>% 
                    select(discrepancy, name, target, target_name, type, alg, study1, study2, size1, size2),
                  ratio.data %>% 
                    mutate(discrepancy = abs(x_shift / sqrt(1/size1 + 1/size2)), name = 'covariate', 
                           target = abs(cond_shift / sqrt(1/size1 + 1/size2)), target_name = 'conditional', 
                           type = "Standardized") %>% 
                    select(discrepancy, name, target, target_name, type, alg, study1, study2, size1, size2), 
                  ratio.data %>% 
                    mutate(discrepancy = abs(cond_shift / sqrt(1/size1 + 1/size2)), name = 'conditional',
                           target = abs(x_shift / sqrt(1/size1 + 1/size2)), target_name = 'covariate', 
                           type = "Standardized") %>% 
                    select(discrepancy, name, target, target_name, type, alg, study1, study2, size1, size2),
                  ratio.data %>% 
                    mutate(discrepancy = stab_x_shift / sqrt(1/size1 + 1/size2), name = 'covariate', 
                           target = abs(cond_shift / sqrt(1/size1 + 1/size2)), target_name = 'conditional', 
                           type = "Standardized + stabilized") %>% 
                    select(discrepancy, name, target, target_name, type, alg, study1, study2, size1, size2), 
                  ratio.data %>% 
                    mutate(discrepancy = abs(cond_shift / sqrt(1/size1 + 1/size2)), name = 'conditional',
                           target = stab_x_shift / sqrt(1/size1 + 1/size2), target_name = 'covariate', 
                           type = "Standardized + stabilized") %>% 
                    select(discrepancy, name, target, target_name, type, alg, study1, study2, size1, size2)) 
  to.plot$Method = factor(to.plot$alg)
  to.plot$Method = plyr::revalue(to.plot$Method, c("double" = "Doubly robust", "ebal" = "Entropy balancing"))
  
  plt.ratio  = to.plot %>% filter(alg== alg_name) %>% 
    ggplot(aes(x = name, y = discrepancy, group = name)) + theme_bw() + 
    geom_violin(aes(linetype = name)) + 
    geom_point(aes(x = name, y = discrepancy), 
               size = 0.3, col = cbPalette[3]) +  
    geom_segment(aes(x = name, y = discrepancy, 
                     xend = target_name, yend = target ), 
                 linewidth = 0.2, alpha = 0.1, col = cbPalette[3]) +  
    facet_wrap(vars(type), scales = 'free') +   
    theme(legend.position="none") + 
    scale_linetype_manual(values = c("dashed", "solid")) + 
    ylab("Measure of distribution shift") + xlab("")
  
  return(list("df" = to.plot, "plot" = plt.ratio))
}




