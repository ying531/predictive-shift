library(haven)
library(ebal)
library(tidyverse)
library(grf)
library(car) 
library(rmarkdown) 
library(rlang)

# prepare for plotting Figure 7


# summarize analysis results for Pipeline

PP.summary.const = data.frame()
for (alg_name in c("ebal", "double")){
  for (study_id in 1:10){ 
    # apply bounds to new generalization tasks
    PI.remain = PP.merge.res %>% filter(data==study_id, alg==alg_name) %>% 
      mutate( 
        PI_low_stab_id = pmin( - 1 * stab_x_shift * cond_sensitivity,
                               1 * stab_x_shift * cond_sensitivity),
        PI_upp_stab_id = pmax(- 1 * stab_x_shift * cond_sensitivity,
                              1 * stab_x_shift * cond_sensitivity), 
        source_study = study_id
      )
    
    # our method (constant calibration) and KL-bounds 
    summary.remain = rbind( 
      PI.remain %>%  
        summarize( 
          coverage = mean((PI_low_stab_id <= delta_yx) * (delta_yx <= PI_upp_stab_id), na.rm = TRUE),
          length = mean(PI_upp_stab_id - PI_low_stab_id, na.rm = TRUE), 
          metric = 'Ours'
        ),
      PI.remain %>%  
        summarize( 
          coverage = mean((PI_lower_KLmax <= theta2) * (theta2 <= PI_upper_KLmax), na.rm = TRUE),
          length = mean(PI_upper_KLmax - PI_lower_KLmax, na.rm = TRUE), 
          metric = 'KL'
        ) 
    )  %>% mutate(source_study = study_id, alg = alg_name)
    
    # transfer assuming IID
    plain.res = PP.PI.all %>% filter(data== study_id, alg == alg_name) %>% 
      mutate(len = PI_upper_plain - PI_lower_plain, 
             cover = (theta2 <= PI_upper_plain) * (PI_lower_plain <= theta2)) %>% 
      group_by(alg, data) %>% summarize(length = mean(len), coverage = mean(cover)) %>% 
      group_by(alg) %>% summarize(length = mean(length), coverage = mean(coverage)) %>%
      mutate(metric = "i.i.d.", source_study = study_id)
    
    plain.res = plain.res[,c("coverage", "length", "metric", "source_study", "alg")]
    
    # Our method with oracle in-study quantiles
    oracle.res = PP.PI.all %>% filter(data == study_id, alg == alg_name) %>% 
      mutate(len = PI_upp_stab_instudy - PI_low_stab_instudy, 
             cover = (theta2 <= PI_upp_stab_instudy) * (PI_low_stab_instudy <= theta2)) %>% 
      group_by(alg, data) %>% 
      summarize(length = mean(len), coverage = mean(cover)) %>% 
      group_by(alg) %>% 
      summarize(length = mean(length), coverage = mean(coverage)) %>%
      mutate(metric = "Oracle", source_study = study_id)
    
    oracle.res = oracle.res[,c("coverage", "length", "metric", "source_study", "alg")]
    
    PP.summary.const = rbind(PP.summary.const, summary.remain, plain.res, oracle.res)
  }
} 


PP.summary.const$metric = factor(PP.summary.const$metric, 
                                 levels = c("i.i.d.", "KL", "Ours", "Oracle")) 
PP.summary.const$metric = plyr::revalue(PP.summary.const$metric, c("KL" = "WorstCase", "i.i.d." = "IID", 'Ours' = "Ours_Const"))
PP.summary.const$source_name = factor(sapply(PP.summary.const$source_study, function(x) paste("Study", x)),
                                      levels = sapply(1:10, function(x) paste("Study", x)))
PP.summary.const$source_study = factor(PP.summary.const$source_study, 
                                       levels = seq(1,10))
PP.max.const = PP.summary.const %>% filter(alg %in% c('ebal', 'grf')) %>% 
  group_by(source_study) %>% summarize(max.length = max(length))
PP.summary.const = PP.summary.const %>% left_join(PP.max.const, by = 'source_study')


# summarize analysis results for ManyLabs1

ML.summary.const = data.frame()
for (alg_name in c("ebal", "double")){
  for (study_id in 1:15){  
    PI.remain = ML.merge.res %>% filter(data==study_id, alg==alg_name) %>% 
      mutate( 
        PI_low_stab_id = pmin( theta1 + delta_x - 1 * stab_x_shift * cond_sensitivity,
                               theta1 + delta_x + 1 * stab_x_shift * cond_sensitivity),
        PI_upp_stab_id = pmax(theta1 + delta_x - 1 * stab_x_shift * cond_sensitivity,
                              theta1 + delta_x +1 * stab_x_shift * cond_sensitivity), 
        source_study = study_id
      ) 
    summary.remain = rbind( 
      PI.remain %>% #group_by(data) %>% 
        summarize( 
          coverage = mean((PI_low_stab_id <= theta2) * (theta2 <= PI_upp_stab_id), na.rm = TRUE),
          length = mean(PI_upp_stab_id - PI_low_stab_id, na.rm = TRUE), 
          metric = 'Ours'
        )
    )  %>% mutate(source_study = study_id, alg = alg_name)
    
    plain.res = ML.PI.all %>% filter(data== study_id, alg == alg_name) %>% 
      mutate(len = PI_upper_plain - PI_lower_plain, 
             cover = (theta2 <= PI_upper_plain) * (PI_lower_plain <= theta2)) %>% 
      group_by(alg, data) %>% 
      summarize(length = mean(len), coverage = mean(cover)) %>% 
      group_by(alg) %>% 
      summarize(length = mean(length), coverage = mean(coverage)) %>%
      mutate(metric = "i.i.d.", source_study = study_id)
    
    oracle.res = ML.PI.all %>% filter(data == study_id, alg == alg_name) %>% 
      mutate(len = PI_upp_stab_instudy - PI_low_stab_instudy, 
             cover = (theta2 <= PI_upp_stab_instudy) * (PI_low_stab_instudy <= theta2)) %>% 
      group_by(alg, data) %>% 
      summarize(length = mean(len), coverage = mean(cover)) %>% 
      group_by(alg) %>% 
      summarize(length = mean(length), coverage = mean(coverage)) %>%
      mutate(metric = "Oracle", source_study = study_id)
    
    # calibrate KL with in-study bounds
    KL.res = ML.cond.KL.PIs %>% filter(data == study_id) %>%
      mutate(len = upper - lower, 
             cover = (theta2 <= upper) * (theta2 >= lower), alg = 'grf') %>% 
      group_by(alg) %>% 
      summarize(length = mean(len), coverage = mean(cover)) %>% 
      mutate(metric = "KL", source_study = study_id)
    
    ML.summary.const = rbind(ML.summary.const, summary.remain, oracle.res, plain.res, KL.res)  
    
  }
} 


ML.summary.const$metric = factor(ML.summary.const$metric, 
                                 levels = c("i.i.d.", "KL","Ours",  "Oracle"))
ML.summary.const$metric = plyr::revalue(ML.summary.const$metric, c("KL" = "WorstCase", "i.i.d." = "IID", 'Ours' = "Ours_Const"))
ML.summary.const$source_name = factor(sapply(ML.summary.const$source_study, function(x) paste("Study", x)),
                                      levels = sapply(1:15, function(x) paste("Study", x)))
ML.summary.const$source_study = factor(ML.summary.const$source_study, 
                                       levels = seq(1,15))

ML.summary.const.eb = ML.summary.const %>% filter(alg %in% c('ebal', 'grf'))
ML.summary.const.eb.max = ML.summary.const.eb %>% group_by(source_study) %>% 
  summarize(max.length = max(length))
ML.summary.const.eb = ML.summary.const.eb %>% left_join(ML.summary.const.eb.max, by = 'source_study')

 
