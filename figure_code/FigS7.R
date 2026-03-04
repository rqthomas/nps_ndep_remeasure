# Marginal effect of N deposition on growth
# Author: Mary Lofton
# Date: 03APR25
# Last update: 22DEC25

# Purpose: figure with marginal effect of N deposition (decrease of 1 kg
# per hectare per year) on growth (kg C per year per individual)

# load packages
library(tidyverse)
library(lubridate)
library(arrow)
library(ggpubr)
library(ggthemes)
library(stringr)
library(scales)
library(metR)
library(grid)
library(alphahull)

# needed for creating factorial combinations of S and N dep and for orthogonalization
mem.maxVSize(vsize = Inf)

# needed to convert ashapes into standard polygons so can get grid points
source("./other_code/hull2poly.R")

# read in raw data
data = "./data/McDonnell_etal_InPrep_TreeData_2024_10_11.csv"

# list files in model output folder
out <- list.files("./experiments/ortho_log_t_interaction_adj_priors",pattern = "mcmc.parquet",
                   full.names = TRUE)

# read in and combine model output files to get posterior distributions of parameters
for(i in 1:length(out)){
  
  spp_id = str_split(out[i], pattern = "-")[[1]][2]
  model_id = str_split(out[i], pattern = "/")[[1]][3]
  temp <- read_parquet(file = out[i]) %>%
    mutate(spp_id = spp_id,
           model_id = model_id)
  
  if(i == 1){
    final <- temp
  } else {
    final <- bind_rows(final, temp)
  }
  
}

# a little data wrangling to get mean values of parameters and correct species names
final1 <- final %>%
  mutate(spp_id = ifelse(spp_id == "yellow","yellow poplar",spp_id)) %>%
  select(model_id, spp_id, global_tree_effect, p2, p3, p5, p6, p7, p8, p9, p10, p11, p12) %>%
  group_by(model_id, spp_id) %>%
  summarise(across(global_tree_effect:p12, \(x) mean(x, na.rm = TRUE)))

# read in pre-processed data (used in models)
df <- read_csv("./data/processed_data.csv")

# read in key of common and scientific names so can use both throughout the script
spp_df <- read_csv(data) %>%
  select(common_name, species) %>%
  distinct(.)

# get mean values needed for model predictions (e.g., mean tree size by ecoregion)
df1 <- left_join(df, spp_df, by = "common_name") %>%
  group_by(species, common_name, ecoregion_ortho) %>%
  summarize(log_mean_size = log(mean(AG_carbon_m1, na.rm = TRUE)),
            log_mean_ba_gt = log(mean(subp_BA_GT_m1, na.rm = TRUE) + 1),
            mean_Dep_Nhistoric = mean(Dep_Nhistoric, na.rm = TRUE),
            mean_Dep_Shistoric = mean(Dep_Shistoric, na.rm = TRUE),
            c = mean(c, na.rm = TRUE),
            mean_size = mean(AG_carbon_m1, na.rm = TRUE),
            n_trees_ecoregion = n_distinct(tree_ID),
            mean_growth = mean(AG_carbon_pYear, na.rm = TRUE),
            mean_perc_growth = round(mean_growth / mean_size * 100, 2)) %>%
  mutate(common_name = ifelse(common_name == "yellow-poplar","yellow poplar",common_name)) %>%
  ungroup()

##### Get model predictions at various levels of recent and antecedent N dep ##################

# vector of ecoregions for for-loop
ecoregions <- unique(df1$ecoregion_ortho)

# placeholders to hold model prediction output
final_pred_df <- NULL
N_ranges <- NULL

# vector of alpha values for creating ahulls for each ecoregion
# these were found through trial and error by running the following:
# N_ahull <- ahull(x = Nrange_dat$Dep_Nhistoric_SO, 
#                  y = Nrange_dat$Dep_Ndiff_SO,
#                  alpha = 3)
# plot(N_ahull, wpoints = TRUE, col = "blue", do.shape = TRUE)
alpha_vals <- c(7, 3, 12, 2, 5, 10, 3)

# for-loop; first loop through ecoregions, then species
tiff(file = file.path("./visualizations/final_figures/FigureS7.png"),
     width = 8, height = 12, units = "in", res = 300)
par(mfrow = c(4,2))#,mgp = c(2.5,1,0), mar = c(3,4,2,0)+0.1)

for(i in 1:length(ecoregions)){
  
  # where are we
  print(ecoregions[i])
  
  # get wrapped title
  wrapped_title = paste(strwrap(ecoregions[i], width = 30), collapse = "\n")
  
  # data from which we will retrieve ranges of N deposition for figure
  Nrange_dat <- df %>%
    mutate(common_name = ifelse(common_name == "yellow-poplar","yellow poplar",common_name)) %>%
    filter(ecoregion_ortho == ecoregions[i]) %>%
    select(plot_ID, ecoregion_ortho, Dep_Nhistoric_SO, Dep_Ndiff_SO) %>%
    distinct(.)
  
  N_ahull <- ahull(x = Nrange_dat$Dep_Nhistoric_SO,
                   y = Nrange_dat$Dep_Ndiff_SO,
                   alpha = alpha_vals[i])
  plot(N_ahull, wpoints = TRUE, col = c("gray","blue","black"), do.shape = TRUE,
       main = wrapped_title, 
       xlab = expression(paste("residual antecedent N deposition (kg N ", ha^-1," ",yr^-1,")")), 
       ylab = expression(paste("residual recent N dep. change  (kg N ", ha^-1," ",yr^-1,")")), 
       lwd = c(1,2,1))
  if(i == 1){
    legend("topleft", # Position (e.g., "topright", "bottomleft")
           legend = c("alpha shape", "alpha hull","observations"), # Text for the legend
           col = c("blue", "gray","black"),      # Colors corresponding to the plot elements
           lty = c(1, 1, NA),              # Line types (optional)
           pch = c(NA, NA, 1),  
           lwd = c(2, 1, NA),
           cex = 0.8,
           bty = "n" 
    )
  } else {
  legend("topright", # Position (e.g., "topright", "bottomleft")
         legend = c("alpha shape", "alpha hull","observations"), # Text for the legend
         col = c("blue", "gray","black"),      # Colors corresponding to the plot elements
         lty = c(1, 1, NA),              # Line types (optional)
         pch = c(NA, NA, 1),  
         lwd = c(2, 1, NA),
         cex = 0.8,
         bty = "n" 
  )
  }
  
}

dev.off()

