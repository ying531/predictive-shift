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


alpha = 0.05 
ROOT_DIR = "/Users/ying/Desktop/Research/Dominik/predictive-shift/"

## ========================================= ##
#      pre-process datasets 
## ========================================= ##

# Pipeline data: run 
#   1. ROOT_DIR/Pipeline/preprocess/process.Rmd
# Manylabs1 data: run 
#   1. ROOT_DIR/ManyLabs1/pre-process/ML1_data_process.R and then 
#   2. ROOT_DIR/ManyLabs1/pre-process/ML1_data_process_2.R

# set directories that save cleaned data files 

# PP_DATA_PATH = ...
# ML_DATA_PATH = ...

## ========================================= ##
#      basic analysis for explanatory role 
## ========================================= ##

# this script saves 
# 1. "results_plain_PP.RData" for iid transfer 
# 2. "results_weighted_PP.RData" for covariate shift transfer 
# at directory ROOT_DIR/Pipeline/exploratory/
source(paste0(ROOT_DIR, "Pipeline/explanatory/explanatory.R"))

# this script saves 
# 1. "results_plain_ML1.RData" for iid transfer 
# 2. "results_weighted_ML1.RData" for covariate shift transfer 
# at directory ROOT_DIR/ManyLabs1/exploratory/
source(paste0(ROOT_DIR, "ManyLabs1/explanatory/explanatory.R"))

# this script summarize the above two projects to create data frames ready to plot Figure 3
source(paste0(ROOT_DIR, "summary/summary_explanatory.R"))

## ================================================ ##
#     basic analysis for bounding/predictive role 
## ================================================ ##

# this script saves 
#    "results_K5_stable_filtered_centered.RData" for shift measures based on our proposals 
# at directory ROOT_DIR/Pipeline/predictive/
source(paste0(ROOT_DIR, "Pipeline/predictive/stable_shift.R"))

# this script saves 
#    "results_stable_ML1.RData" for shift measures based on our proposals 
# at directory ROOT_DIR/ManyLabs1/predictive/
source(paste0(ROOT_DIR, "ManyLabs/predictive/stable_shift.R"))

# this script uses the above two results to create data frames ready to plot Figure 4
source(paste0(ROOT_DIR, "summary/summary_predictive.R"))


## ======================================================= ##
#     analysis for generalization (constant calibration)
## ======================================================= ##

# this script saves 
#    "cond_KL_delta_PP.RData" for estimated conditional KL measures across sites
# at directory ROOT_DIR/Pipeline/generalization/
source(paste0(ROOT_DIR, "Pipeline/generalization/compute_KL_delta.R"))

# this script saves 
#    "cond_KL_PIs_PP.RData" for KL-worst-case PIs with KL bounds calibrated in-study ("WorstCase" in Figure 7)
# at directory ROOT_DIR/Pipeline/generalization/
source(paste0(ROOT_DIR, "Pipeline/generalization/KL_no_aux_data.R"))


# this script saves 
#    "cond_KL_delta_ML1.RData" for estimated conditional KL measures across sites
# at directory ROOT_DIR/ManyLabs1/generalization/
source(paste0(ROOT_DIR, "ManyLabs1/generalization/compute_KL_delta.R"))

# this script saves 
#    "cond_KL_PIs_PP.RData" for KL-worst-case PIs with KL bounds calibrated in-study ("WorstCase" in Figure 7)
# at directory ROOT_DIR/Pipeline/generalization/
source(paste0(ROOT_DIR, "ManyLabs1/generalization/KL_no_aux_data.R"))


# this script uses the computed distribution shift measures and KL bounds to 
# create data frames ready to plot Figure 7
source(paste0(ROOT_DIR, "summary/summary_const_calib.R"))


## =========================================================== ##
#     analysis for generalization (data-adaptive calibration)
## =========================================================== ##

# this script saves 
#    "KL_calib_study.RData" for KL-worst-case PIs with KL bounds calibrated with data ("WorstCase" in Figure 8)
# at directory ROOT_DIR/Pipeline/generalization/
source(paste0(ROOT_DIR, "Pipeline/generalization/study_adaptive.R"))


# this script saves 
#    "KL_calib_study_ML1.RData" for KL-worst-case PIs with KL bounds calibrated with data ("WorstCase" in Figure 8)
# at directory ROOT_DIR/ManyLabs1/generalization/
source(paste0(ROOT_DIR, "ManyLabs1/generalization/study_adaptive.R"))


# this scripts combines the above two results to produce data frames ready for plotting Figure 8
source(paste0(ROOT_DIR, "summary/summary_data_calib.R"))



## ============================== ##
#     generate plots  
## ============================== ##

# this script uses all above results to produce plots in the main text
# saving files in ROOT_DIR/plots/.
source(paste0(ROOT_DIR, "plots_main.R"))

