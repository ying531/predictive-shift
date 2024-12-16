# created Oct 5, 2024
# final plotting code combining pipeline and manylabs1 results
 
library(haven)
library(ebal)
library(tidyverse)
library(grf)
library(car) 
library(rmarkdown) 
library(rlang)
library(formattable)
library(paletteer)
library(ggplot2)
library(ggpubr)
library(scales)
library(ggsci) 



cbPalette <- c(   "#B0B0B0", paletteer_c("viridis::viridis", n=3)[2], "#895C96") 
const.palette = c("#B0B0B0", paletteer_c("viridis::viridis", n=3)[1], 
                  "#C52A20", paletteer_c("viridis::viridis", n=3)[3])  
study.palette = const.palette
 
alpha = 0.05

ROOT_DIR = ""
 
## ===================================== ##
## read utils and merge analysis results ##
## ===================================== ##

source(paste0(ROOT_DIR, "analysis_utils.R"))
source(paste0(ROOT_DIR, "summary/read_merge_data.R"))

ML_DATA_PATH = "./Manylabs1_data.RData"
PP_DATA_PATH = "./clean_data/"
load(ML_DATA_PATH)

## ========================================= ##
##        Code for creating Figure 3         ##
## ========================================= ##

##=== Figure 3: coverage of IID and CovShift

source(paste0(ROOT_DIR, "summary/summary_explanatory.R"))

# plot for Pipeline 

PP.base.plt = PP.df.base %>%
  ggplot(aes(x = data, y = coverage, group = method)) + theme_bw() + 
  geom_bar(aes(fill = method), width=0.7,  
           col = 'black', linewidth = 0.15,
           stat='identity', position=position_dodge(0.72, preserve = "single")) + 
  geom_hline(yintercept = 0.95, col = 'red', linetype = 'dashed') + 
  theme(legend.position = 'top', legend.title = element_text(size = 14)) +
  xlab("Hypothesis index (Pipeline)") + ylab("Empirical coverage") + 
  scale_color_manual(values = cbPalette) + 
  scale_fill_manual(values = cbPalette) +
  labs(color='Method', fill = "Method")  

# plot for ManyLabs1 
 
ML.base.plt = ML.df.base %>%  
  ggplot(aes(x = data, y = coverage, group = method)) + theme_bw() + 
  geom_bar(aes(fill = method), width=0.7,  
           col = 'black', linewidth = 0.2,
           stat='identity', position=position_dodge(0.72, preserve = "single")) + 
  geom_hline(yintercept = 0.95, col = 'red', linetype = 'dashed') + 
  theme(legend.position = 'top', legend.title = element_text(size = 14)) +
  xlab("Hypothesis index (ManyLabs1)") + ylab("Empirical coverage") + 
  scale_color_manual(values = cbPalette) + 
  scale_fill_manual(values = cbPalette) +
  labs(color='Method', fill = "Method") 

    
##=== Figure 3: estimate for reweight

PP.df.rw.pi = PP.PI.all %>% filter(data == 5, study2 == 5) %>% 
  select(theta1, theta2, study, study1, study2, alg, delta_x, delta_yx, delta_total) 

PP.plt.rw = PP.df.rw.pi %>% filter(alg=='double') %>% 
  ggplot() + theme_bw() + 
  geom_hline(yintercept = PP.df.rw.pi$theta2[1], col = 'red', linetype = 'dashed') + 
  geom_segment(aes(x = 1, y = theta1, 
                   xend = 2, yend = theta1+delta_x), linewidth = 0.3, col = cbPalette[2]) + 
  geom_segment(aes(x = 1, y = theta1, 
                   xend = 0, yend = theta1 + delta_x), 
               data = PP.df.rw.pi %>% filter(alg=='ebal'),
               linewidth = 0.3, col = cbPalette[3]) + 
  geom_point(aes(x = 1, y = theta1), col = cbPalette[1]) + 
  geom_point(aes(x = 2, y = theta1+delta_x), col = cbPalette[2]) + 
  geom_point(aes(x = 0, y = theta1 + delta_x), 
             data = PP.df.rw.pi %>% filter(alg=='ebal'),
             col = cbPalette[3]) + 
  scale_x_continuous(breaks=c(0, 1,2), labels = c("Reweight (EB)", "Source", "Reweight (DR)"),
                     expand = expansion(mult = c(0.15, 0.15))) + 
  xlab("Method") + ylab("Value of estimate")  

ML.df.rw.pi = ML.PI.all %>% filter(data == 4, study2 == 4) %>% 
  select(theta1, theta2, study, study1, study2, hat_theta_x, alg, delta_x, delta_yx, delta_total) 

ML.plt.rw = ML.df.rw.pi %>% filter(alg=='double') %>% 
  ggplot() + theme_bw() + 
  geom_hline(yintercept = ML.df.rw.pi$theta2[1]/1000, col = 'red', linetype = 'dashed') + 
  geom_segment(aes(x = 1, y = theta1/1000, 
                   xend = 2, yend = hat_theta_x/1000), linewidth = 0.3, col = cbPalette[2]) + 
  geom_segment(aes(x = 1, y = theta1/1000, 
                   xend = 0, yend = theta1/1000 + delta_x/1000), 
               data = ML.df.rw.pi %>% filter(alg=='ebal'),
               linewidth = 0.3, col = cbPalette[3]) + 
  geom_point(aes(x = 1, y = theta1/1000), col = cbPalette[1]) + 
  geom_point(aes(x = 2, y = hat_theta_x/1000), col = cbPalette[2]) + 
  geom_point(aes(x = 0, y = theta1/1000 + delta_x/1000), 
             data = ML.df.rw.pi %>% filter(alg=='ebal'),
             col = cbPalette[3]) + 
  scale_x_continuous(breaks=c(0, 1,2), labels = c("CovShift (EB)", "Source", "CovShift (DR)"),
                     expand = expansion(mult = c(0.15, 0.15))) + 
  xlab("Method") + ylab("Value of estimate")  
 
# put everything together 

PPML.plt.1 = ggarrange(
  NULL, PP.base.plt, NULL, PP.plt.rw, NULL, ML.base.plt, NULL, ML.plt.rw, 
  labels = c("", "(P, a)", "", "(P, b)", "", "(M, a)", "", "(M, b)"), 
  widths = c(0.15,1.2, 0.1, 1, 0.15,1.2, 0.1, 1), 
  vjust = 1.5, nrow=2, ncol = 4, hjust=0.6,
  common.legend = TRUE, legend = "top"
) 

ggsave(paste0(ROOT_DIR, "plots/explanatory_PPML.pdf"), PPML.plt.1, width = 8.5, height=6, units='in')




## ============================================== ##
##  Code for creating Figure 4  (bounding role)   ##
## ============================================== ##

source(paste0(ROOT_DIR, "summary/summary_predictive.R"))

###=== Create four cross-context plots in panels (a,b)

ML.plt.measure.context.ebal = 
  ML.df.both %>% filter(type == 'Standardized + stabilized') %>% 
  filter(
    !( (compare == 'Between_US_NonUS') & 
         ( paste0(study1,"_",study2) %in% paste0(ML.both.remove$study1, "_", ML.both.remove$study2 ))) 
  ) %>%
  filter(alg== 'ebal') %>% 
  ggplot(aes(x = name, y = discrepancy, group = interaction(name, compare))) + theme_bw() + 
  geom_violin(aes(linetype = name)) + 
  geom_point(aes(x = name, y = discrepancy), 
             size = 0.3, col = cbPalette[3]) +  
  geom_segment(aes(x = name, y = discrepancy, 
                   xend = target_name, yend = target ), 
               linewidth = 0.2, alpha = 0.1, col = cbPalette[3]) +  
  facet_wrap(vars(compare)) +  
  theme(legend.position="none") + 
  scale_linetype_manual(values = c("dashed", "solid")) + 
  ylab("Measure of distr. shift") + xlab("Type of distribution shift (ManyLabs1)")  

ML.plt.measure.site.ebal = 
  ML.67.remove %>% filter(type == 'Standardized + stabilized') %>% 
  filter(!(  (type == 'Standardized') & (name == 'covariate') & (discrepancy > 500) ), 
         !(  (type == 'Standardized') & (name == 'conditional') & (target > 500))
  ) %>% 
  filter(alg== 'ebal') %>% 
  ggplot(aes(x = name, y = discrepancy, group = interaction(name, study))) + theme_bw() + 
  geom_violin(aes(linetype = name)) + 
  geom_point(aes(x = name, y = discrepancy), 
             size = 0.3, col = cbPalette[3]) +  
  geom_segment(aes(x = name, y = discrepancy, 
                   xend = target_name, yend = target ), 
               linewidth = 0.2, alpha = 0.01, col = cbPalette[3]) +  
  facet_wrap(vars(study)) +  
  theme(legend.position="none") + 
  scale_linetype_manual(values = c("dashed", "solid")) + 
  ylab("Measure of distr. shift") + xlab("Type of distribution shift (ManyLabs1)") 


PP.plt.measure.context.ebal = 
  PP.df.both %>% filter(type == 'Standardized + stabilized') %>% 
  filter(!( (compare == 'Within_US') & (type == 'Standardized + stabilized') & (name == 'conditional') & (discrepancy > 3) ),
         !( (compare == 'Within_US') & (type == 'Non-standardized') & (name == 'conditional') & (discrepancy > 0.5) ),
         !( (compare == 'Within_US') & (type == 'Standardized + stabilized') & (name == 'covariate') & (target > 3) ),
         !( (compare == 'Within_US') & (type == 'Non-standardized') & (name == 'covariate') & (target > 0.5) ),
         !(  (type == 'Standardized') & (name == 'covariate') & (discrepancy > 500) ), 
         !(  (type == 'Standardized') & (name == 'conditional') & (target > 500))
  ) %>% 
  filter(alg== 'ebal') %>% 
  ggplot(aes(x = name, y = discrepancy, group = interaction(name, compare))) + theme_bw() + 
  geom_violin(aes(linetype = name)) + 
  geom_point(aes(x = name, y = discrepancy), 
             size = 0.3, col = cbPalette[3]) +  
  geom_segment(aes(x = name, y = discrepancy, 
                   xend = target_name, yend = target ), 
               linewidth = 0.2, alpha = 0.1, col = cbPalette[3]) +  
  facet_wrap(vars(compare) ) +  
  theme(legend.position="none") + 
  scale_linetype_manual(values = c("dashed", "solid")) + 
  ylab("Measure of distr. shift") + xlab("Type of distribution shift (Pipeline)")  


PP.plt.measure.site.ebal = 
  PP.df.67 %>% filter(type == 'Standardized + stabilized') %>% 
  filter(!(  (type == 'Standardized') & (name == 'covariate') & (discrepancy > 500) ), 
         !(  (type == 'Standardized') & (name == 'conditional') & (target > 500))
  ) %>% 
  filter(alg== 'ebal') %>% 
  ggplot(aes(x = name, y = discrepancy, group = interaction(name, study))) + theme_bw() + 
  geom_violin(aes(linetype = name)) + 
  geom_point(aes(x = name, y = discrepancy), 
             size = 0.3, col = cbPalette[3]) +  
  geom_segment(aes(x = name, y = discrepancy, 
                   xend = target_name, yend = target ), 
               linewidth = 0.2, alpha = 0.1, col = cbPalette[3]) +  
  facet_wrap(vars(study) ) + 
  # facet_grid(rows = vars(type_short), cols = vars(study), scales = 'free_y') +   
  theme(legend.position="none") + 
  scale_linetype_manual(values = c("dashed", "solid")) + 
  ylab("Measure of distr. shift") + xlab("Type of distribution shift (Pipeline)")  


context.plt.PPML = 
  ggarrange(
    NULL, PP.plt.measure.context.ebal, NULL, PP.plt.measure.site.ebal, 
    NULL, ML.plt.measure.context.ebal, NULL, ML.plt.measure.site.ebal, 
    labels = c("", "(P, a)", "", "(P, b)" ,"", "(M, a)", "", "(M, b)"  ), 
    widths = rep(c( 0.08, 1, 0.1, 1), 2), heights = rep(1, 8),
    vjust = 1.6, nrow=2, ncol = 4, hjust= 0.2 
  )   

 
###=== Create empirical quantiles in panel (c) 

PPML.qt.plt = PPML.qt.all %>% filter(alg=='ebal') %>% 
  ggplot(aes(x = alpha, y = upp, group = data, linetype = "Empirical upper quantile within hypothesis", 
             col="Empirical upper quantile within hypothesis")) + theme_bw() + 
  geom_line(aes(x=alpha, y=low, group=data, linetype = "Empirical lower quantile within hypothesis", 
                col="Empirical lower quantile within hypothesis"),  alpha = 0.4, size=0.4) + 
  geom_line(   alpha = 0.4, size=0.4) + 
  geom_line(aes(x=alpha, y=upp.ref, linetype = "Upper/lower normal quantile", 
                col="Upper/lower normal quantile"), size=0.4)+   
  geom_line(aes(x=alpha, y=upp.ref/2, linetype = "0.5 * Upper/lower normal quantile", 
                col="0.5 * Upper/lower normal quantile"),   size=0.3) + 
  geom_line(aes(x=alpha, y=upp.ref*0.75, linetype = "0.75 * Upper/lower normal quantile", 
                col="0.75 * Upper/lower normal quantile"),   size=0.3)+
  geom_line(aes(x=alpha, y=upp.ref/4, linetype = "0.25 * Upper/lower normal quantile", 
                col="0.25 * Upper/lower normal quantile"),  size=0.3)+
  geom_line(aes(x=alpha, y=-upp.ref, linetype = "Upper/lower normal quantile", 
                col="Upper/lower normal quantile"), size=0.4)+   
  geom_line(aes(x=alpha, y=-upp.ref/2,  linetype = "0.5 * Upper/lower normal quantile", 
                col="0.5 * Upper/lower normal quantile"), size=0.3) + 
  geom_line(aes(x=alpha, y=-upp.ref*0.75, linetype = "0.75 * Upper/lower normal quantile", 
                col="0.75 * Upper/lower normal quantile"), size=0.3)+
  geom_line(aes(x=alpha, y=-upp.ref/4, linetype = "0.25 * Upper/lower normal quantile", 
                col="0.25 * Upper/lower normal quantile"), size=0.3)+
  ylim(c(-3, 3))  + facet_wrap(vars(project)) +  
  ylab("Quantiles") + xlab("Confidence level (alpha)") + 
  scale_color_manual( 
    values = c(
      "Empirical upper quantile within hypothesis" = "#0B132B",
      "Empirical lower quantile within hypothesis" = "#6A4E55",
      "Upper/lower normal quantile"= "red",
      "0.5 * Upper/lower normal quantile" = "red",
      "0.75 * Upper/lower normal quantile" = "red",
      "0.25 * Upper/lower normal quantile" = "red" 
    )
  ) +
  scale_linetype_manual( 
    values = c(
      "Empirical upper quantile within hypothesis" = "solid",
      "Empirical lower quantile within hypothesis" = "solid",
      "Upper/lower normal quantile"= "solid",
      "0.5 * Upper/lower normal quantile" = "dotdash",
      "0.75 * Upper/lower normal quantile" = "dashed",
      "0.25 * Upper/lower normal quantile" = "dotted" 
    )) +  
  guides(color = guide_legend(title = "Curve Type", nrow=6), 
         linetype = guide_legend(title = "Curve Type", nrow=6)) + 
  theme(legend.position = 'right', strip.text = element_text(size = 10))

 
# put everything together 

all.ratio.plt.PPML =
  ggarrange(
    context.plt.PPML, NULL,
    ggarrange(NULL, PPML.qt.plt, NULL, 
              labels = c("", "(c)", ""),
              widths = c(0.03, 0.8, 0.02), nrow=1),  
    heights = c(1,0.02, 0.55),
    vjust = 2, nrow=3, ncol = 1, hjust= 0.2 
  )  

ggsave(paste0(ROOT_DIR, "plots/combine_ratio_plot_PPML.pdf"), all.ratio.plt.PPML, width=9, height=7.5, units='in')
 

## ============================================== ##
##    Code for Figure 7 (constant calibration)    ##
## ============================================== ##

source(paste0(ROOT_DIR, "summary/summary_const_calib.R"))

ML.plt.const.cover = ML.summary.const.eb %>% 
  ggplot(aes(x = source_study, y = coverage, group = metric)) + theme_bw() + 
  geom_bar(aes(fill = metric), width=0.7,  
           col = 'black', linewidth = 0.2,
           stat='identity', position=position_dodge(0.725, preserve = "single")) + 
  geom_hline(yintercept=0.95, col = 'red', linetype = 'dashed') + 
  scale_fill_manual(values =  const.palette) + 
  xlab("Hypothesis index (ManyLabs1)") + ylab("Empirical coverage") + 
  guides(fill=guide_legend(title="Measure")) + 
  scale_y_continuous(trans = shift_trans(0.2)) + theme(plot.margin = margin(0.15,0.1,0.15,0.1, "in")) 

ML.plt.const.length = ML.summary.const.eb %>% 
  ggplot(aes(x = source_study, y = length /max.length, group = metric)) + theme_bw() + 
  geom_bar(aes(fill = metric), width=0.7,  
           col = 'black', linewidth = 0.2,
           stat='identity', position=position_dodge(0.725, preserve = "single")) +  
  scale_fill_manual(values = const.palette) + 
  xlab("Hypothesis index (ManyLabs1)") + ylab("Relative PI lengths") +  
  guides(fill=guide_legend(title="Measure"))+ theme(plot.margin = margin(0.15,0.1,0.15,0.1, "in")) 
 
PP.plt.const.cover = PP.summary.const %>% filter(alg == 'ebal', metric != 'KL (qt)') %>% 
  ggplot(aes(x = source_study, y = coverage, group = metric)) + theme_bw() + 
  geom_bar(aes(fill = metric), width=0.85,  
           col = 'black', linewidth = 0.1,
           stat='identity', position=position_dodge(0.875, preserve = "single")) + 
  geom_hline(yintercept=0.95, col = 'red', linetype = 'dashed') + 
  scale_fill_manual(values = const.palette) + 
  xlab("Hypothesis index (Pipeline)") + ylab("Empirical coverage") + 
  guides(fill='none')+ theme(plot.margin = margin(0.15,0.1,0.15,0.1, "in")) 

PP.plt.const.length = PP.summary.const %>% filter(alg == 'ebal', metric != 'KL (qt)') %>% 
  ggplot(aes(x = source_study, y = length/max.length, group = metric)) + theme_bw() + 
  geom_bar(aes(fill = metric), width=0.85,  
           col = 'black', linewidth = 0.1,
           stat='identity', position=position_dodge(0.875, preserve = "single")) +  
  scale_fill_manual(values = const.palette) + 
  xlab("Hypothesis index (Pipeline)") + ylab("Relative PI lengths") +  
  guides(fill='none') + theme(plot.margin = margin(0.15,0.1,0.15,0.1, "in")) 


# putting everything together 

PPML.plt.const.PI = ggarrange(
  NULL, PP.plt.const.cover, ML.plt.const.cover,  
  NULL,  PP.plt.const.length, ML.plt.const.length, 
  labels = c("", "(a)", "",  "", "(b)", ""), 
  widths = c(0.05,1, 1.5,   0.05,1, 1.5), 
  heights = c(1,1,1,1,1,1),
  vjust = 1.5, nrow=2, ncol = 3, hjust= 0,
  common.legend = TRUE, legend = "top"
)
ggsave(paste0(ROOT_DIR, "plots/const_eb_KL_PPML.pdf"), PPML.plt.const.PI, width = 8.5, height=6, units='in')



## =================================================== ##
##    Code for Figure 7 (data-adaptive calibration)    ##
## =================================================== ##

source(paste0(ROOT_DIR, "summary/summary_data_calib.R"))

ML.study.cover.plt = ML.study.res %>% filter(alg != 'double', source_nums > 2 ) %>%
  ggplot(aes(x = source_nums, y = coverage , group = interaction(metric, irep))) + theme_bw() + 
  geom_bar(aes(x = source_nums, y = coverage , fill = metric),
           data = ML.study.sum %>% filter(alg == 'ebal', source_nums > 2), width=0.7,
           col = 'black', linewidth = 0.2, 
           stat='identity', position=position_dodge(0.725, preserve = "single")) +
  geom_hline(yintercept = 1-alpha, col = 'red', linetype = 'dashed') +
  scale_alpha_manual(values = c(rep(0.3, 2), 0.3,0.3 )) +  
  scale_color_manual(values = study.palette) +
  scale_fill_manual(values = study.palette) + 
  xlab("# Existing hypotheses (ManyLabs1)") + ylab("Empirical coverage")  + 
  scale_y_continuous(trans = shift_trans(0.5))+ 
  guides(fill='none') + 
  theme(plot.margin = margin(0.1,0.1,0.1,0.1, "in")) 


ML.study.len.plt = ML.study.sum %>% filter(alg == 'ebal', source_nums > 2) %>%
  ggplot(aes(x = source_nums, y = length/max.length, group = metric)) + theme_bw() +
  geom_bar(aes(fill = metric), width=0.7,
           col = 'black', linewidth = 0.2,
           stat='identity', position=position_dodge(0.725, preserve = "single")) + 
  scale_color_manual(values = study.palette) +
  scale_fill_manual(values = study.palette) + 
  xlab("# Existing hypotheses (ManyLabs1)") + ylab("Relative PI lengths")+ 
  guides(fill='none') + 
  theme(plot.margin = margin(0.1,0.1,0.1,0.1, "in")) 


PP.study.cover.plt = PP.study.res %>% filter(alg == 'ebal', source_nums > 2 ) %>%
  ggplot(aes(x = source_nums, y = coverage , group = interaction(metric, irep))) + theme_bw() + 
  geom_bar(aes(x = source_nums, y = coverage , fill = metric),
           data = PP.study.sum %>% filter(alg == 'ebal', source_nums > 2), width=0.7,
           col = 'black', linewidth = 0.2, 
           stat='identity', position=position_dodge(0.725, preserve = "single")) +
  geom_hline(yintercept = 1-alpha, col = 'red', linetype = 'dashed') +
  scale_alpha_manual(values = c(rep(0.3, 2), 0.3,0.3 )) +  
  scale_color_manual(values = study.palette) +
  scale_fill_manual(values = study.palette) + 
  xlab("# Existing hypotheses (Pipeline)") + ylab("Empirical coverage")  + 
  scale_y_continuous(trans = shift_trans(0.5))+ 
  theme(plot.margin = margin(0.1,0.1,0.1,0.1, "in")) 

PP.study.len.plt = PP.study.sum %>% filter(alg == 'ebal', source_nums > 2) %>%
  ggplot(aes(x = source_nums, y = length/max.length, group = metric)) + theme_bw() +
  geom_bar(aes(fill = metric), width=0.7,
           col = 'black', linewidth = 0.2,
           stat='identity', position=position_dodge(0.725, preserve = "single")) + 
  scale_color_manual(values = study.palette) +
  scale_fill_manual(values = study.palette) + 
  xlab("# Existing hypotheses (Pipeline)") + ylab("Relative PI lengths") + 
  theme(plot.margin = margin(0.1,0.1,0.1,0.1, "in")) 



xy.grid = data.frame(expand.grid(1:5, 1:10)) %>% mutate(order = Var1)
colnames(xy.grid) = c("Study", "Site", "order") 

study.order.plt = xy.grid %>% ggplot(aes(Study, Site, fill=order)) + theme_void() + 
  geom_tile(color = "white",
            lwd = 0.5,
            linetype = 1) + coord_fixed() + xlab("Hypothesis ID") + ylab("Site ID") + 
  theme(legend.position = "none",
        axis.title.x = element_text(size=10),
        axis.title.y = element_text(angle = 90, size=10)) + 
  scale_fill_gradient(low="black", high= "#D8D8D8") + 
  guides(fill="none")

study.PI.plt.PPML =
  ggarrange(
    NULL, study.order.plt, 
    ggarrange(NULL, PP.study.cover.plt,  ML.study.cover.plt ,
              NULL, PP.study.len.plt,   ML.study.len.plt,
              labels = c("", "(a)", "", 
                         "", "(b)", ""),
              widths = rep(c( 0.1, 1, 1.5),2),
              vjust = 1.5, nrow=2, ncol = 3, hjust= 0.6,
              common.legend = TRUE, legend = "top"
    ),
    widths = c(0.03, 0.3, 2.2),  labels = c("", "",""),
    vjust = 13, nrow=1, ncol=3, hjust= 0.2 
  ) 

ggsave(paste0(ROOT_DIR, "plots/study_PI_eb_KL_PPML.pdf"), study.PI.plt.PPML, width=9, height=6, unit = 'in')



 
## ========================================= ##
##   Code for Figure 2 (preview of results)
## ========================================= ## 

df.main.plt = rbind(
  PP.summary.const %>% filter(alg == 'ebal', source_study == 7) %>% 
    select(coverage, metric) %>% 
    mutate(method = metric, Metric = 'PI Coverage', y = coverage) %>% select(y, method, Metric),
  PP.summary.const %>% filter(alg == 'ebal', source_study == 7) %>% 
    select(length, metric) %>% 
    mutate(method = metric, Metric = 'PI Length', y = length/7.92) %>% select(y, method, Metric))

df.main.plt$method = plyr::revalue(df.main.plt$method, c("Ours_Const" = "Ours_Const", 
                                                  "Oracle" = "Ours_Oracle", 
                                                  "WorstCase" = "WorstCase", "IID"="IID"))

plt.main.cal = df.main.plt %>% ggplot(aes(x = Metric, y = y, group = method)) + theme_bw() + 
  geom_bar(aes(fill = method), width=0.5,  
           col = 'black', linewidth = 0.2,
           stat='identity', position=position_dodge(0.51, preserve = "single")) +  
  geom_hline(yintercept = 0.95, col = 'red', linetype = 'dashed') + 
  scale_fill_manual(values = const.palette) + 
  xlab("") + ylab("Coverage / Relative PI lengths") + 
  theme(legend.title = element_text(size=10),
        legend.text = element_text(size=7),
        legend.spacing.x = unit(0.1, "cm"), 
        legend.justification = 'right',
        legend.key.width=unit(0.5,"cm"),
        legend.key.height=unit(0.5,"cm")) +
  labs(color='Method', fill = "Method") 
   

base.plt = PP.df.base %>% mutate(dataname = paste0("H", data)) %>% 
  filter(!(data %in% c(1,2,8,9,10))) %>% 
  ggplot(aes(x = dataname, y = coverage+0.005, group = method)) + theme_bw() + 
  geom_bar(aes(fill = method), width=0.65,  
           col = 'black', linewidth = 0.2,
           stat='identity', position=position_dodge(0.66, preserve = "single")) + 
  geom_hline(yintercept = 0.95, col = 'red', linetype = 'dashed') + 
  theme(legend.position = 'top', legend.title = element_text(size = 14)) +
  xlab("") + ylab("Coverage across all site pairs") + 
  scale_color_manual(values = cbPalette) + 
  scale_fill_manual(values = cbPalette) +
  theme(legend.title = element_text(size=10),
        legend.text = element_text(size=7), 
        legend.key.width=unit(0.5,"cm"),
        legend.key.height=unit(0.5,"cm")) +
  labs(color='Method', fill = "Method")  + ylim(c(0,1))

plt.main = ggarrange(
  NULL, base.plt, NULL, plt.main.cal, 
  labels = c("", "(a)", "", "(b)"), widths = c(0.05,1, 0.1, 1), 
  vjust = 1.2, nrow=1, hjust= 0.5, 
  legend = "top"
) 
 
ggsave(paste0(ROOT_DIR, "plots/main_plot.pdf"), plt.main, width = 8, height=3, units='in')

 

