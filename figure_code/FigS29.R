# Marginal effect of N deposition on growth
# Author: Mary Lofton
# Date: 03APR25
# Last updated: 22DEC25

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
library(ggh4x)
library(ggpattern)

# needed for creating factorial combinations of S and N dep and for orthogonalization
mem.maxVSize(vsize = Inf)

# needed to convert ashapes into standard polygons so can get grid points
source("./other_code/hull2poly.R")

# read in raw data
data = "./data/McDonnell_etal_InPrep_TreeData_2024_10_11.csv"

# list files in model output folder
out <- list.files("./experiments/ortho_log_t_interaction_01JAN26",pattern = "mcmc.parquet",
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
  filter(!(spp_id == "sugar maple" & .iteration <= 1500)) %>%
  filter(!(spp_id == "ponderosa pine" & .iteration <= 1000)) %>%
  select(model_id, spp_id, global_tree_effect, p2, p3, p5, p6, p7, p8, p9, p10, p11, p12)

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
            mean_growth = round(mean(AG_carbon_pYear, na.rm = TRUE), 2),
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
for(i in 1:length(ecoregions)){
  
  # where are we
  print(ecoregions[i])
  
  # data from which we will retrieve ranges of N deposition for figure
  Nrange_dat <- df %>%
    mutate(common_name = ifelse(common_name == "yellow-poplar","yellow poplar",common_name)) %>%
    filter(ecoregion_ortho == ecoregions[i]) %>%
    select(plot_ID, ecoregion_ortho, Dep_Nhistoric_SO, Dep_Ndiff_SO) %>%
    distinct(.)
  
  # get alpha convex hull
  eco_ashape <- ashape(x = Nrange_dat$Dep_Nhistoric_SO, 
                             y = Nrange_dat$Dep_Ndiff_SO,
                             alpha = alpha_vals[i])
  
  # convert to regular polygon
  eco_poly <- hull2poly(eco_ashape)
  
  # get grid polygons
  grid_polygons <- st_make_grid(eco_poly, 
                                cellsize = c(0.1, 0.1), # Define the size of grid cells
                                what = "polygons") |> 
    st_sf() # Convert to an sf object
  
  # filter grid cells that intersect with the ashape polygon
  grid_within_ahull <- st_intersection(grid_polygons, eco_poly)
  
  # convert the resulting grid polygons to centroids
  grid_points <- st_centroid(grid_within_ahull)
  
  # convert to dataframe
  N_grid <- data.frame(st_coordinates(grid_points)) %>%
    rename(Dep_Nhistoric_SO = X,
           Dep_Ndiff_SO = Y) %>%
    add_row(Dep_Nhistoric_SO = c(0,0,-1),
            Dep_Ndiff_SO = c(0,-1,0)) 
  
  # map from SO to DO
  rtilde_n = N_grid$Dep_Ndiff_SO
  atilde_n = N_grid$Dep_Nhistoric_SO
  
  n_n = length(rtilde_n)
  
  P_n = diag(n_n) - (rtilde_n %*% solve(t(rtilde_n) %*% rtilde_n, t(rtilde_n)))
  
  astar_n = as.numeric(P_n %*% atilde_n)
  
  N_grid$Dep_Nhistoric_DO_scaled = scale(astar_n)
  N_grid$Dep_Ndiff_DO_scaled = scale(rtilde_n)
  
  st_N_DO_values <- N_grid %>%
    filter(Dep_Ndiff_SO %in% c(0,-1) & Dep_Nhistoric_SO == 0)
  lt_N_DO_values <- N_grid %>%
    filter(Dep_Ndiff_SO == 0 & Dep_Nhistoric_SO %in% c(0,-1))
  
  # get names of species in that ecoregion
  common_names = df %>%
    filter(ecoregion_ortho == ecoregions[i]) %>% 
    select(common_name) %>%
    distinct(.) %>%
    mutate(common_name = ifelse(common_name == "yellow-poplar","yellow poplar",common_name)) %>%
    pull(common_name) 
  
  # loop through species in each ecoregion
  for(j in 1:length(common_names)){
    
    # where are we
    print(common_names[j])
    
    # get mean values from pre-processed data (e.g., mean tree size per ecoregion) for that species
    model_data <- df1 %>% filter(common_name == common_names[j] & ecoregion_ortho == ecoregions[i])
    
    # get values of model parameters for that species
    model_output <- final1 %>% filter(spp_id == common_names[j],
                                      model_id == "ortho_log_t_interaction_01JAN26")
    params <- model_output[sample(nrow(model_output), 1000), ]
    
    # back transformation of model parameters to SO space
    # p5 is coefficient on recent N dep; p12 is coefficient on recent N dep quadratic term
    # this is not needed if we work in DO space
    # p5_so <- params$p5 - model_data$c * params$p9
    # p12_so <- params$p12 - model_data$c * params$p11
  
    # get model predictions!
    
    # short-term
    pred_log_baseline_st <- ((params$global_tree_effect + st_N_DO_values$Dep_Ndiff_DO_scaled[1]*params$p5 + 0*params$p6 + 0*params$p7 + 0*params$p8 + st_N_DO_values$Dep_Nhistoric_DO_scaled[1]*params$p9 + model_data$mean_Dep_Shistoric*params$p10 + st_N_DO_values$Dep_Nhistoric_DO_scaled[1]*st_N_DO_values$Dep_Ndiff_DO_scaled[1]*params$p11 + st_N_DO_values$Dep_Ndiff_DO_scaled[1]*st_N_DO_values$Dep_Ndiff_DO_scaled[1]*params$p12) 
                        + params$p2*model_data$log_mean_size) - params$p3*log(1)
    m2_baseline_st = exp(pred_log_baseline_st + model_data$log_mean_size)
    pred_baseline_st = m2_baseline_st - model_data$mean_size
    
    pred_baseline <- round(mean(pred_baseline_st),2)
    
    pred_log_increase_st <- ((params$global_tree_effect + st_N_DO_values$Dep_Ndiff_DO_scaled[2]*params$p5 + 0*params$p6 + 0*params$p7 + 0*params$p8 + st_N_DO_values$Dep_Nhistoric_DO_scaled[2]*params$p9 + model_data$mean_Dep_Shistoric*params$p10 + st_N_DO_values$Dep_Nhistoric_DO_scaled[2]*st_N_DO_values$Dep_Ndiff_DO_scaled[2]*params$p11 + st_N_DO_values$Dep_Ndiff_DO_scaled[2]*st_N_DO_values$Dep_Ndiff_DO_scaled[2]*params$p12) 
                             + params$p2*model_data$log_mean_size) - params$p3*log(1)
    m2_increase_st = exp(pred_log_increase_st + model_data$log_mean_size)
    pred_increase_st = m2_increase_st - model_data$mean_size
    
    growth_delta_st <- pred_increase_st - pred_baseline_st
    
    # long-term
    pred_log_baseline_lt <- ((params$global_tree_effect + lt_N_DO_values$Dep_Ndiff_DO_scaled[1]*params$p5 + 0*params$p6 + 0*params$p7 + 0*params$p8 + lt_N_DO_values$Dep_Nhistoric_DO_scaled[1]*params$p9 + model_data$mean_Dep_Shistoric*params$p10 + lt_N_DO_values$Dep_Nhistoric_DO_scaled[1]*lt_N_DO_values$Dep_Ndiff_DO_scaled[1]*params$p11 + lt_N_DO_values$Dep_Ndiff_DO_scaled[1]*lt_N_DO_values$Dep_Ndiff_DO_scaled[1]*params$p12) 
                             + params$p2*model_data$log_mean_size) - params$p3*log(1)
    m2_baseline_lt = exp(pred_log_baseline_lt + model_data$log_mean_size)
    pred_baseline_lt = m2_baseline_lt - model_data$mean_size
    
    pred_log_increase_lt <- ((params$global_tree_effect + lt_N_DO_values$Dep_Ndiff_DO_scaled[2]*params$p5 + 0*params$p6 + 0*params$p7 + 0*params$p8 + lt_N_DO_values$Dep_Nhistoric_DO_scaled[2]*params$p9 + model_data$mean_Dep_Shistoric*params$p10 + lt_N_DO_values$Dep_Nhistoric_DO_scaled[2]*lt_N_DO_values$Dep_Ndiff_DO_scaled[2]*params$p11 + lt_N_DO_values$Dep_Ndiff_DO_scaled[2]*lt_N_DO_values$Dep_Ndiff_DO_scaled[2]*params$p12) 
                             + params$p2*model_data$log_mean_size) - params$p3*log(1)
    m2_increase_lt = exp(pred_log_increase_lt + model_data$log_mean_size)
    pred_increase_lt = m2_increase_lt - model_data$mean_size
    
    growth_delta_lt <- pred_increase_lt - pred_baseline_lt
    
    pred <- c(growth_delta_st, growth_delta_lt)
    
    # assemble final df for this species-ecoregion combination
    pred_df <- data.frame(species = model_data$species[1],
                          model_id = params$model_id[1],
                          ecoregion = model_data$ecoregion_ortho[1],
                          pred = pred,
                          pred_baseline = pred_baseline,
                          change_type = rep(c("short-term","antecedent"),each = 1000)
                          )
  
    if(i == 1 & j == 1){
      final_pred_df <- pred_df
    } else {
      final_pred_df <- bind_rows(final_pred_df, pred_df)
    }
  
  }
  
}

# get predictions for all ecoregions with < 100 trees
figs20_eco <- df1 %>%
  filter(n_trees_ecoregion < 100) %>%
  rename(ecoregion = ecoregion_ortho) %>%
  select(species, ecoregion, n_trees_ecoregion)

letters_df <- final_pred_df %>%
  select(species, ecoregion) %>%
  distinct(.) %>%
  right_join(., figs20_eco) %>%
  arrange(species, ecoregion) %>%
  mutate(letters = paste0(letters[seq_len(nrow(figs20_eco))],"."))

facet_label_df <- final_pred_df %>%
  select(species, ecoregion, pred_baseline) %>%
  distinct(.) %>%
  right_join(., letters_df, by = c("species","ecoregion")) %>%
  mutate(pred_baseline = round(pred_baseline, 2)) %>%
  separate_wider_delim(species, delim = " ", names = c("genus","spp"), cols_remove = FALSE) %>%
  separate_wider_delim(ecoregion, delim = "&", names = c("er1","er2"), too_few = "align_start",cols_remove = FALSE) %>%
  mutate(er1 = str_replace_all(er1, " ", "~"),
         er2 = str_replace_all(er2, " ", "~"),
         er1 = str_remove(er1, "~$"),
         er2 = str_remove(er2, "^~"),
         facet_labels = paste0("atop(italic(",letters,"~",genus,"~",spp,"), atop(atop(",er1, ",",er2,"),expected~growth:~",pred_baseline,"~kg~C~y^-1~ind^-1))")) %>%
  select(species, ecoregion, facet_labels)

plot_data <- right_join(final_pred_df, facet_label_df, by = c("species", "ecoregion") ) %>%
  arrange(species) %>%
  mutate(perc_change = pred / pred_baseline * 100) 

my_col <- c(RColorBrewer::brewer.pal(3, "Blues")[3],"orange")

spp_col <- data.frame(panel_fill_cols = colorblind_pal()(8)[c(2,6)],
                      species = unique(plot_data$species))
show_col(spp_col$panel_fill_cols)

strip_cols <- right_join(facet_label_df, spp_col, by = c("species")) %>%
  arrange(species, ecoregion) %>%
  pull(panel_fill_cols)
better_strip_cols <- alpha(strip_cols, 0.3)

x_lims <- plot_data %>%
  select(species, perc_change) %>%
  group_by(species) %>%
  summarize(min = round(min(perc_change, na.rm = TRUE),2),
            max = round(max(perc_change, na.rm = TRUE),2)) %>%
  right_join(facet_label_df, by = "species")

facet_lims_lst <- NULL
for(i in 1:nrow(x_lims)){
  facet_lims_lst[[i]] <- scale_x_continuous(limits = c(x_lims$min[i],x_lims$max[i]))
}

p <- ggplot(data = plot_data)+
  ggpattern::geom_density_pattern(aes(x = perc_change, group = change_type, pattern = change_type, fill = change_type),
                                  color = "black", pattern_color = "white")+
  ggh4x::facet_wrap2(~ facet_labels, scales = "free", labeller = label_parsed,
                     # Use strip_themed() to apply aesthetics
                     strip = ggh4x::strip_themed(
                       # Map 'fill' to the 'cyl_factor' variable
                       background_x = lapply(better_strip_cols, element_rect, color = "black"),
                     ),
                     nrow = 7, ncol = 4
  ) +
  # Manually set colors for the factor levels for clarity
  theme_bw()+
  scale_fill_manual(values = c("white","gray"))+
  geom_vline(xintercept = 0)+
  labs(color = "", fill = "", pattern = "", x = expression(paste("% change in growth")))+
  ggtitle(expression(paste("Marginal effect of N deposition decrease (per kg N ", ha^-1," ",y^-1,")")))+
  theme(strip.text = element_text(size = 12))+
  ggh4x::facetted_pos_scales(x = facet_lims_lst)+
  scale_pattern_manual(values = c("none","circle"))

ggsave(p,filename = "./visualizations/final_figures/FigureS28.png",
       device = "png", height = 3, width = 10, units = "in")

