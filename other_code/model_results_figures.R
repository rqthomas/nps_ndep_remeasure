# Model results figure
# Author: Mary Lofton
# Date: 03APR25

# Purpose: two panel figure, first panel is delta growth vs delta Ndep,
# second panel is histograms of marginal effect of delta Ndep on growth
# per year per individual

# load packages
library(tidyverse)
library(lubridate)
library(arrow)
library(ggpubr)
library(ggthemes)
library(stringr)

# delta growth vs delta Ndep
source("./other_code/delta_growth_vs_delta_Ndep.R")
p1 <- delta_growth_vs_delta_Ndep(data = "./data/McDonnell_etal_InPrep_TreeData_2024_10_11.csv", 
                           model_output_folder = "./experiments/new_delta_Ndep_only_saveTreeEffect",
                           experiment_name = "new_delta_Ndep_only_saveTreeEffect")
p1
experiment_name = "new_delta_Ndep_only_saveTreeEffect"
ggsave(plot = p1, filename = paste0("./visualizations/delta_growth_vs_delta_Ndep-",experiment_name,".tif"),
       device = "tiff", height = 5, width = 8, units = "in",bg = "white")

p2 <- delta_growth_vs_delta_Ndep(data = "./data/McDonnell_etal_InPrep_TreeData_2024_10_11.csv", 
                                 model_output_folder = "./experiments/delta_env_saveTreeEffect",
                                 experiment_name = "delta_env_saveTreeEffect")
p2
experiment_name = "delta_env_saveTreeEffect"
ggsave(plot = p2, filename = paste0("./visualizations/delta_growth_vs_delta_Ndep-",experiment_name,".tif"),
       device = "tiff", height = 4, width = 6, units = "in",bg = "white")

source("./other_code/marginal_effect_Ndep.R")
p3 <- marginal_effect_Ndep(data = "./data/McDonnell_etal_InPrep_TreeData_2024_10_11.csv")
p3
ggsave(plot = p3, filename = "./visualizations/marginal_effect_Ndep.tif",
       device = "tiff", height = 8, width = 10, units = "in",bg = "white")

source("./other_code/marginal_effect_Ndep_species.R")
p4 <- marginal_effect_Ndep_species(data = "./data/McDonnell_etal_InPrep_TreeData_2024_10_11.csv")
p4
ggsave(plot = p4, filename = "./visualizations/marginal_effect_Ndep_species.tif",
       device = "tiff", height = 8, width = 10, units = "in",bg = "white")

source("./other_code/marginal_effect_Ndep_longterm_shortterm.R")
p5 <- marginal_effect_Ndep_longterm_shortterm(data = "./data/McDonnell_etal_InPrep_TreeData_2024_10_11.csv")
p5
ggsave(plot = p5, filename = "./visualizations/marginal_effect_Ndep_short-term_long-term.tif",
       device = "tiff", height = 8, width = 10, units = "in",bg = "white")

source("./other_code/growth_vs_size.R")
p6 <- growth_vs_size(data = "./data/McDonnell_etal_InPrep_TreeData_2024_10_11.csv", 
                     model_output_folder = "./experiments/short-term_long-term")
p6
ggsave(plot = p6, filename = "./visualizations/growth_vs_size_short-term.tif",
       device = "tiff", height = 10, width = 10, units = "in",bg = "white")

source("./other_code/growth_vs_size_longtermDepN.R")
p7 <- growth_vs_size_longtermDepN(data = "./data/McDonnell_etal_InPrep_TreeData_2024_10_11.csv", 
                     model_output_folder = "./experiments/short-term_long-term")
p7
ggsave(plot = p7, filename = "./visualizations/growth_vs_size_long-term.tif",
       device = "tiff", height = 10, width = 10, units = "in",bg = "white")

source("./other_code/short-term_vs_long-term.R")
p8 <- shortterm_vs_longterm(data = "./data/McDonnell_etal_InPrep_TreeData_2024_10_11.csv", 
                                  model_output_folder = "./experiments/short-term_long-term",
                            zoom = FALSE)
p8
ggsave(plot = p8, filename = "./visualizations/short-term_vs_long-term.tif",
       device = "tiff", height = 5, width = 7, units = "in",bg = "white")

p9 <- shortterm_vs_longterm(data = "./data/McDonnell_etal_InPrep_TreeData_2024_10_11.csv", 
                            model_output_folder = "./experiments/short-term_long-term",
                            zoom = TRUE)
p9
ggsave(plot = p9, filename = "./visualizations/short-term_vs_long-term_zoom.tif",
       device = "tiff", height = 5, width = 7, units = "in",bg = "white")
