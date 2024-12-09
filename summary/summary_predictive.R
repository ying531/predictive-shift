library(haven)
library(ebal)
library(tidyverse)
library(grf)
library(car) 
library(rmarkdown) 
library(rlang)

# prepare for plotting Figure 4


# across location, Pipeline

us.id = c(1,2,3,4,5,7,9,10,11,12,15,17,19,22,24,25,26) 
euro.id = c(6,8,13,16,18,27) 

data.id = 6
PP.res.within.us = generate_ratio_plot(PP.PI.all %>% filter(data == data.id, study1 %in% us.id, study2 %in% us.id), 
                                       alg_name = 'ebal')  
PP.df.within.us = PP.res.within.us$df

PP.res.between = generate_ratio_plot(PP.PI.all %>% 
                                       filter(data == data.id, ((study1 %in% euro.id) & (study2 %in% us.id)) | ((study2 %in% euro.id) & (study1 %in% us.id))), 
                                     alg_name = 'ebal')  
PP.df.between = PP.res.between$df 

PP.df.both = rbind(PP.df.within.us %>% mutate("compare" = "Within_US"),
                   PP.df.between %>% mutate("compare" = "Between_US_Euro"))
PP.df.both$type_short = plyr::revalue(factor(PP.df.both$type), 
                                      c("Non-standardized" = "Non-stand.",
                                        "Standardized" = "Stand.", 
                                        "Standardized + stabilized" = "Stand. + stab."))
# across hypothesis, Pipeline

PP.res.6 = generate_ratio_plot(PP.PI.all %>% filter(data == 6), alg_name = 'ebal')  
PP.df.6 = PP.res.6$df

PP.res.7 = generate_ratio_plot(PP.PI.all %>% filter(data == 5), alg_name = 'ebal')  
PP.df.7 = PP.res.7$df 

PP.df.67 = rbind(PP.df.6 %>% mutate(study = "Hypothesis 6"),
                 PP.df.7 %>% mutate(study = "Hypothesis 5"))
PP.df.67$type_short = plyr::revalue(factor(PP.df.67$type), 
                                    c("Non-standardized" = "Non-stand.",
                                      "Standardized" = "Stand.", 
                                      "Standardized + stabilized" = "Stand. + stab."))


# across location, Manylabs1

ML.study.names = unique(pdata$original_study)
ML.site.names = unique(pdata$site)
ML.site.ifus = c(1,0,0,1,1,0,1,1,0,0,0,1,1,0, rep(1,8), 0,0, 1,1,1,0,1,0,rep(1,6))
us.id = (1:36)[ML.site.ifus==1]
nonus.id = (1:36)[ML.site.ifus==0] 

data.id = 1
ML.res.within.us = generate_ratio_plot(ML.PI.all %>% 
                                         filter(data == data.id, study1 %in% us.id, study2 %in% us.id), 
                                       alg_name = 'ebal')  
ML.df.within.us = ML.res.within.us$df

ML.res.between = generate_ratio_plot(ML.PI.all %>% 
                                       filter(data == data.id, 
                                              ((study1 %in% nonus.id) & (study2 %in% us.id)) | ((study2 %in% nonus.id) & (study1 %in% us.id))), 
                                     alg_name = 'ebal')  
ML.df.between = ML.res.between$df 

ML.df.both = rbind(ML.df.within.us %>% mutate("compare" = "Within_US"),
                   ML.df.between %>% mutate("compare" = "Between_US_NonUS"))
ML.df.both$type_short = plyr::revalue(factor(ML.df.both$type), 
                                      c("Non-standardized" = "Non-stand.",
                                        "Standardized" = "Stand.", 
                                        "Standardized + stabilized" = "Stand. + stab."))

# across hypothesis, Manylabs1

ML.res.6 = generate_ratio_plot(ML.PI.all %>% filter(data == 3), alg_name = 'ebal')  
ML.df.6 = ML.res.6$df 
ML.res.7 = generate_ratio_plot(ML.PI.all %>% filter(data == 4), alg_name = 'ebal')  
ML.df.7 = ML.res.7$df 

ML.df.67 = rbind(ML.df.6 %>% mutate(study = "Hypothesis 3"),
                 ML.df.7 %>% mutate(study = "Hypothesis 4"))
ML.df.67$type_short = plyr::revalue(factor(ML.df.67$type), 
                                    c("Non-standardized" = "Non-stand.",
                                      "Standardized" = "Stand.", 
                                      "Standardized + stabilized" = "Stand. + stab."))

# remove extreme values for easier visualization (doesn't change takeaway message)

ML.both.remove = ML.df.both %>% filter(
  ( (compare == 'Within_US') & (type == 'Standardized + stabilized') &
      (name == 'conditional') & (discrepancy > 3)) | 
    ( (compare == 'Between_US_NonUS') & (type == 'Standardized + stabilized') &
        (name == 'covariate') & (discrepancy > 20))
) 

ML.67.remove = ML.df.67 %>% filter(
  !(( (study == 'Hypothesis 3') & (type == 'Standardized + stabilized') &
        (name == 'covariate') & (discrepancy > 15)) | 
      ( (study == 'Hypothesis 3') & (type == 'Standardized + stabilized') &
          (name == 'conditional') & (target > 15)) |
      ( (study == 'Hypothesis 4') & (type == 'Standardized + stabilized') &
          (name == 'covariate') & (discrepancy > 20)) | 
      ( (study == 'Hypothesis 4') & (type == 'Standardized + stabilized') &
          (name == 'conditional') & (target > 20)))
) 




## ==========================================
### empirical quantiles in panel (c)

PP.qt.all = data.frame()
for (alpha in seq(0.02, 0.2, by=0.01)){
  PP.qt.all = rbind(PP.qt.all, PP.analysis.studies %>%  
                      group_by(data, alg) %>% 
                      summarize(upp = quantile(cond_shift / stab_x_shift, 1-alpha/2),
                                low = quantile(cond_shift / stab_x_shift,  alpha/2)) %>% mutate("alpha" = alpha))
}

ML.qt.all = data.frame()
for (alpha in seq(0.02, 0.2, by=0.01)){
  ML.qt.all = rbind(ML.qt.all, ML.analysis.studies %>%  
                      group_by(data, alg) %>% 
                      summarize(upp = quantile(cond_shift / stab_x_shift, 1-alpha/2),
                                low = quantile(cond_shift / stab_x_shift,  alpha/2)) %>% mutate("alpha" = alpha))
}

PP.qt.all$upp.ref = qnorm(1-PP.qt.all$alpha/2)
PP.qt.all$data = factor(PP.qt.all$data)  
ML.qt.all$upp.ref = qnorm(1-ML.qt.all$alpha/2)
ML.qt.all$data = factor(ML.qt.all$data) 

PPML.qt.all = rbind(PP.qt.all %>% mutate("project" = "Pipeline"), 
                    ML.qt.all %>% mutate("project" = "ManyLabs1"))
PPML.qt.all$project = factor(PPML.qt.all$project, levels=c("Pipeline", "ManyLabs1"))


 
