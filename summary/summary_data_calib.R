alpha = 0.05
load(paste0(ROOT_DIR, "Pipeline/generalization/KL_calib_hypo.RData"))
PP.KL.study.res = study.res

PP.study.res = data.frame()

for (irep in 1:10){ 
  set.seed(irep)
  study.order = sample(1:10, 10)
  
  for (jj in 1:9){
    exist.id = study.order[1:jj]
    new.id = study.order[(1+jj):10]
    
    qt.id = PP.analysis.studies %>% filter(data %in% exist.id, alg == 'ebal') %>% 
      summarize( 
        qt_upp_stab = quantile(cond_shift / stab_x_shift, 1-alpha/2),
        qt_low_stab = quantile(cond_shift / stab_x_shift, alpha/2)  
      )   
    
    plain.res = PP.PI.all %>% filter(data %in% new.id) %>% 
      mutate(len = PI_upper_plain - PI_lower_plain, 
             cover = (theta2 <= PI_upper_plain) * (PI_lower_plain <= theta2)) %>% 
      group_by(alg, data) %>% summarize(length = mean(len), coverage = mean(cover)) %>% 
      group_by(alg) %>% summarize(length = mean(length), coverage = mean(coverage)) %>%
      mutate(metric = "i.i.d.", source_nums = jj, seed = irep)
    
    plain.res = plain.res[,c("coverage", "length", "metric", "source_nums", "alg", "seed")]
    
    oracle.res = PP.PI.all %>% filter(data %in% new.id, alg == 'ebal') %>% 
      mutate(len = PI_upp_stab_instudy - PI_low_stab_instudy,
             cover = (theta2 <= PI_upp_stab_instudy) * (PI_low_stab_instudy <= theta2)) %>% 
      group_by(alg, data) %>% summarize(length = mean(len), coverage = mean(cover)) %>% 
      group_by(alg) %>% summarize(length = mean(length), coverage = mean(coverage)) %>% 
      mutate(metric = "Oracle", source_nums = jj, seed = irep)
    oracle.res = oracle.res[,c("coverage", "length", "metric", "source_nums", "alg", "seed")]
    
    our.res = data.frame()
    
    our.res = rbind(our.res, 
                    PP.merge.res %>% filter(data %in% new.id, alg=='ebal') %>% 
                      mutate( 
                        PI_low_stab_id = pmin( qt.id$qt_low_stab * stab_x_shift * cond_sensitivity,
                                               qt.id$qt_upp_stab * stab_x_shift * cond_sensitivity),
                        PI_upp_stab_id = pmax( qt.id$qt_low_stab * stab_x_shift * cond_sensitivity,
                                               qt.id$qt_upp_stab * stab_x_shift * cond_sensitivity), 
                        source_nums = jj
                      )  %>% group_by(data) %>% 
                      summarize( 
                        coverage = mean((PI_low_stab_id <= delta_yx) * (delta_yx <= PI_upp_stab_id), na.rm = TRUE),
                        length = mean(PI_upp_stab_id - PI_low_stab_id, na.rm = TRUE), 
                        metric = 'Ours'
                      ) %>% mutate(source_nums = jj, alg = 'ebal', seed = irep)) 
    
    our.res = our.res[,c("coverage", "length", "metric", "source_nums", "alg", "seed")]
    
    kl.res = PP.KL.study.res %>% filter(source_nums == jj, seed == irep, alg == 'ebal') %>% 
      group_by(metric, source_nums, alg, seed) %>% 
      summarize(coverage = mean(coverage), length = mean(length)) 
    
    kl.res = kl.res[,c("coverage", "length", "metric", "source_nums", "alg", "seed")]
    
    PP.study.res = rbind(PP.study.res, our.res, kl.res, plain.res, oracle.res) 
    
  }
}
PP.study.res$metric = factor(PP.study.res$metric, 
                             levels = c("i.i.d.", "KL", "Ours", "Oracle"))
PP.study.res$metric = plyr::revalue(PP.study.res$metric, c("KL" = "WorstCase", "i.i.d." = "IID"))



# study.res$source_nums = as.factor(study.res$source_nums)
PP.study.sum = PP.study.res %>% 
  group_by(metric, source_nums, alg) %>% 
  summarize(coverage = mean(coverage), 
            length = mean(length))

PP.study.max = PP.study.sum %>% group_by(source_nums, alg) %>% 
  summarize(max.length = max(length))
PP.study.sum = PP.study.sum %>% left_join(PP.study.max, by = c('source_nums', 'alg'), suffix = c('', '_max'))





# ML1

study.orders = matrix(0, nrow = 10, ncol = 15)

for (ss in 1:10){
  set.seed(ss) 
  study.orders[ss,] = sample(1:15, 15)  
}

load(paste0(ROOT_DIR, "ManyLabs1/generalization/KL_calib_hypo_ML1.RData"))

  
ML.study.res = data.frame()

for (irep in 1:10){  
  study.order = study.orders[irep,]
  
  for (alg_name in c("ebal", "double")){
    for (jj in 1:14){
      exist.id = study.order[1:jj]
      new.id = study.order[(1+jj):15]
      
      KL.res = ML.KL.study.all %>% filter(seed == irep, num_exist == jj) %>% 
        mutate(len = upper - lower, 
               cover = (theta2 <= upper) * (lower <= theta2))  %>% 
        group_by(data) %>% summarize(length = mean(len), coverage = mean(cover)) %>% 
        summarize(length = mean(length), coverage = mean(coverage)) %>%
        mutate(metric = 'KL', source_nums = jj, seed = irep, alg = 'ebal')
      
      plain.res = ML.PI.all %>% filter(data %in% new.id) %>% 
        mutate(len = PI_upper_plain - PI_lower_plain, 
               cover = (theta2 <= PI_upper_plain) * (PI_lower_plain <= theta2)) %>% 
        group_by(alg, data) %>% summarize(length = mean(len), coverage = mean(cover)) %>% 
        summarize(length = mean(length), coverage = mean(coverage)) %>%
        mutate(metric = "i.i.d.", source_nums = jj, seed = irep)
      
      plain.res = plain.res[,c("coverage", "length", "metric", "source_nums", "alg", "seed")]
      
      oracle.res = ML.PI.all %>% filter(data %in% new.id ) %>% 
        mutate(len = PI_upp_stab_instudy - PI_low_stab_instudy, 
               cover = (theta2 <= PI_upp_stab_instudy) * (PI_low_stab_instudy <= theta2)) %>% 
        group_by(alg, data) %>% 
        summarize(length = mean(len), coverage = mean(cover)) %>% 
        group_by(alg) %>% 
        summarize(length = mean(length), coverage = mean(coverage)) %>%
        mutate(metric = "Oracle", source_nums = jj, seed = irep)
      
      oracle.res = oracle.res[,c("coverage", "length", "metric", "source_nums", "alg", "seed")]
      
      ML.study.res = rbind(ML.study.res, KL.res, plain.res, oracle.res)  
      
      for (alg_name in c('ebal', 'double')){
        qt.exist = ML.merge.res %>% filter(data %in% exist.id, alg == alg_name) %>% 
          group_by(alg) %>% 
          summarize(qt_upp = quantile(cond_shift / stab_x_shift, 1-alpha/2),
                    qt_low = quantile(cond_shift / stab_x_shift, alpha/2))
        PI.remain = ML.merge.res %>% filter(data %in% new.id, alg==alg_name) %>% 
          mutate( 
            PI_low_stab = pmin( theta1 + delta_x + qt.exist$qt_upp * stab_x_shift * cond_sensitivity,
                                theta1 + delta_x + qt.exist$qt_low * cond_sensitivity),
            PI_upp_stab = pmax(theta1 + delta_x + qt.exist$qt_upp * stab_x_shift * cond_sensitivity,
                               theta1 + delta_x + qt.exist$qt_low * stab_x_shift * cond_sensitivity) 
          )
        our.res = PI.remain %>% 
          mutate(len = PI_upp_stab - PI_low_stab, 
                 cover = (theta2 <= PI_upp_stab) * (PI_low_stab <= theta2)) %>% 
          group_by(data) %>% 
          summarize(length = mean(len), coverage = mean(cover)) %>%   
          summarize(length = mean(length), coverage = mean(coverage)) %>%
          mutate(metric = "Ours", source_nums = jj, seed = irep, alg = alg_name)
        ML.study.res = rbind(ML.study.res, our.res)
      }  
    } 
  }
}
ML.study.res$metric = factor(ML.study.res$metric, 
                             levels = c("i.i.d.", "KL", "Ours", "Oracle"))
ML.study.res$metric = plyr::revalue(ML.study.res$metric, c("KL" = "WorstCase", "i.i.d." = "IID")) 
ML.study.sum = ML.study.res %>% 
  group_by(metric, source_nums, alg) %>% 
  summarize(coverage = mean(coverage), 
            length = mean(length))
ML.study.max = ML.study.sum %>% group_by(source_nums, alg) %>% 
  summarize(max.length = max(length))
ML.study.sum = ML.study.sum %>% left_join(ML.study.max, by = c('source_nums', 'alg'), suffix = c('', '_max'))

