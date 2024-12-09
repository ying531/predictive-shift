library(haven)
library(ebal)
library(tidyverse)
library(grf)
library(car) 
library(rmarkdown) 
library(rlang)

# prepare for plotting Figure 3


PP.df.base = rbind(
  PP.PI.all %>% filter(alg == 'double') %>% 
    mutate(lower = PI_lower_trans, upper = PI_upper_trans, method = 'CovShift (DR)') %>% 
    select(data, study, study1, study2, theta1, theta2, lower, upper, method),
  PP.PI.all %>% filter(alg == 'ebal') %>% 
    mutate(lower = PI_lower_trans, upper = PI_upper_trans, method = 'CovShift (EB)') %>% 
    select(data, study, study1, study2, theta1, theta2, lower, upper, method),
  PP.PI.all %>% filter(alg == 'double') %>% 
    mutate(lower = PI_lower_plain, upper = PI_upper_plain, method = 'IID') %>% 
    select(data, study, study1, study2, theta1, theta2, lower, upper, method) 
) %>% group_by(data, method) %>% 
  summarize(coverage = mean((theta2 <= upper) * (theta2 >= lower), na.rm = TRUE),
            length = mean(upper-lower), na.rm = TRUE)

PP.df.base$data = factor(PP.df.base$data, level = seq(1, 10))
PP.df.base$method = factor(PP.df.base$method, level = c("IID", "CovShift (DR)", "CovShift (EB)"))


ML.df.base = rbind(
  ML.PI.all %>% filter(alg == 'double') %>%
    mutate(lower = PI_lower_trans, upper = PI_upper_trans, method = 'CovShift (DR)') %>%
    select(data, study, study1, study2, theta1, theta2, lower, upper, method),
  ML.PI.all %>% filter(alg == 'ebal') %>%
    mutate(lower = PI_lower_trans, upper = PI_upper_trans, method = 'CovShift (EB)') %>%
    select(data, study, study1, study2, theta1, theta2, lower, upper, method),
  ML.PI.all %>% filter(alg == 'double') %>%
    mutate(lower = PI_lower_plain, upper = PI_upper_plain, method = 'IID') %>%
    select(data, study, study1, study2, theta1, theta2, lower, upper, method)
) %>% group_by(data, method) %>%
  summarize(coverage = mean((theta2 <= upper) * (theta2 >= lower), na.rm = TRUE),
            length = mean(upper-lower, na.rm = TRUE))

ML.df.base$method = factor(ML.df.base$method, level = c("IID", "CovShift (DR)", "CovShift (EB)"))
ML.df.base$data = factor(ML.df.base$data, level = seq(1, 15))
ML.df.base.m = ML.df.base %>% group_by(data) %>% summarize(max.length = max(length))
ML.df.base = ML.df.base %>% left_join(ML.df.base.m, by = 'data')
 





 
