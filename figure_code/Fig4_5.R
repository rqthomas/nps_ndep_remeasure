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

#### Data wrangling

# get predictions for most common ecoregion
letters_df <- data.frame(species = sort(unique(final_pred_df$species)),
                         letters = paste0(letters[seq_len(length(unique(final_pred_df$species)))],"."))

facet_label_df <- final_pred_df %>%
  select(species, ecoregion, pred_baseline) %>%
  distinct(.) %>%
  left_join(., letters_df, by = "species") %>%
  mutate(pred_baseline = round(pred_baseline, 2)) %>%
  separate_wider_delim(species, delim = " ", names = c("genus","spp"), cols_remove = FALSE) %>%
  separate_wider_delim(ecoregion, delim = "&", names = c("er1","er2"), too_few = "align_start",cols_remove = FALSE) %>%
  mutate(er1 = str_replace_all(er1, " ", "~"),
         er2 = str_replace_all(er2, " ", "~"),
         er1 = str_remove(er1, "~$"),
         er2 = str_remove(er2, "^~"),
         facet_labels = paste0("atop(italic(",letters,"~",genus,"~",spp,"), atop(atop(",er1, ",",er2,"),expected~growth:~",pred_baseline,"~kg~C~y^-1~ind^-1))")) %>%
  select(species, ecoregion, facet_labels)

df2 <- df1 %>%
  group_by(species, common_name) %>%
  summarize(n_trees_total = sum(n_trees_ecoregion))

spp <- unique(df1$species)

for(i in 1:length(spp)){
  
  n_trees_total <- df2 %>%
    filter(species == spp[i]) %>%
    pull(n_trees_total)
  
  most_n_eco <- df1 %>%
    filter(species == spp[i]) %>%
    filter(n_trees_ecoregion == max(n_trees_ecoregion)) %>%
    pull(ecoregion_ortho)
  

    species_pred_eco <- final_pred_df %>%
      mutate(perc_change = pred / pred_baseline * 100) %>%
      filter(species == spp[i] & ecoregion == most_n_eco) %>%
      pull(perc_change)
  
  final_pred_species <- data.frame(species = spp[i],
                                   change_type = rep(c("short-term","antecedent"),each = 1000),
                                   pred = species_pred_eco, 
                                   ecoregion = most_n_eco)
  
  if(i == 1){
    final_pred <- final_pred_species
  } else {
    final_pred <- bind_rows(final_pred, final_pred_species)
  }
  
  
  
}

plot_data <- left_join(final_pred, facet_label_df, by = c("species", "ecoregion") ) %>%
  arrange(species)

#### Figure 4

p <- ggplot(data = plot_data)+
  ggpattern::geom_density_pattern(aes(x = pred, group = change_type, pattern = change_type, fill = change_type),
                                  color = "black", pattern_color = "white")+
  facet_wrap(~facet_labels, scales = "free_y", labeller = label_parsed)+
  theme_bw()+
  scale_fill_manual(values = c("white","gray"))+
  geom_vline(xintercept = 0)+
  labs(color = "", fill = "",pattern = "", x = expression(paste("% change in growth")))+
  ggtitle(expression(paste("Marginal effect of N deposition decrease (per kg N ", ha^-1," ",y^-1,")")))+
  scale_pattern_manual(values = c("none","circle"))+
  theme(strip.text = element_text(size = 12),
        strip.background = element_rect(fill = "white"),
        legend.position = c(0.85, 0.2),
        legend.key.size = unit(1, "cm"),
        legend.text = element_text(size = 16))
ggsave(p,filename = "./visualizations/final_figures/Figure4.png",
       device = "png", height = 6.5, width = 9.5, units = "in")

fig4_ecos <- final_pred %>%
  select(species, ecoregion) %>%
  distinct(.)

mean_pred_responses <- final_pred_df %>%
  group_by(species, ecoregion, change_type, pred_baseline) %>%
  summarize(mean_pred = mean(pred, na.rm = TRUE)) %>%
  mutate(perc_change = mean_pred / pred_baseline * 100) %>%
  arrange(change_type, perc_change) %>%
  right_join(., fig4_ecos)

##### Figure 5

plot_data2 <- plot_data %>%
  group_by(species, change_type, ecoregion) %>%
  summarize(pred = mean(pred, na.rm = TRUE)) %>%
  pivot_wider(names_from = change_type, values_from = pred) %>%
  add_column(spp_abbrevs = c("Acsa","Bepa","Litu","Piru","Pipo","Pode","Potr","Prse"))

plot_data5 <- plot_data %>%
  arrange(species, ecoregion, change_type) %>%
  add_column(draw = rep(c(1:1000),times = 16)) %>%
  group_by(species, change_type, ecoregion, draw) %>%
  pivot_wider(names_from = change_type, values_from = pred) %>%
  ungroup() %>%
  select(-draw) %>%
  mutate(species = factor(species, levels = c("Liriodendron tulipifera","Populus deltoides","Prunus serotina","Acer saccharum","Betula papyrifera","Picea rubens","Populus tremuloides","Pinus ponderosa")))

my.cols.og <- hue_pal()(3)
my.cols.og
my.cols <- c("#F8766D","#faa49e","#fdd1ce","#008729","#00BA38","#00fa4b","#7affa2","#619CFF")

spp_eco <- plot_data2 %>%
  select(ecoregion, species, spp_abbrevs) %>%
  arrange(ecoregion) %>%
  distinct(.)

plot_data3 <- plot_data2 %>%
  mutate(antecedent_lab = ifelse(spp_abbrevs %in% c("Potr","Bepa"),antecedent + 0.4,
                             ifelse(spp_abbrevs == "Acsa",antecedent + 0.4, antecedent))) %>%
  mutate(`short-term_lab` = ifelse(spp_abbrevs %in% c("Potr","Bepa"),`short-term` + 0.4,
                             ifelse(spp_abbrevs == "Acsa",`short-term` - 0.4, `short-term`))) %>%
  mutate(species = factor(species, levels = c("Liriodendron tulipifera","Populus deltoides","Prunus serotina","Acer saccharum","Betula papyrifera","Picea rubens","Populus tremuloides","Pinus ponderosa")))


seg_dat <- plot_data3 %>%
  filter(spp_abbrevs %in% c("Acsa","Bepa","Potr"))

fig5 <- ggplot()+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  geom_point(data = plot_data5, aes(x = antecedent, y = `short-term`, group = species, color = species), size = 0.5, alpha = 0.5)+
  geom_segment(data = seg_dat, aes(x = antecedent_lab, y = `short-term_lab`, xend = antecedent, yend = `short-term`))+
  geom_label(data = plot_data3, aes(x = antecedent_lab, y = `short-term_lab`, group = species, label = spp_abbrevs, fill = species))+
  #xlim(c(-0.4, 0.4))+
  #ylim(c(-0.4, 0.4))+
  theme_classic()+
  ylab(expression(atop(paste("% growth change w/ 1 kg ",ha^-1," ",yr^-1),"short-term N dep. decrease")))+
  xlab(expression(atop(paste("% growth change w/ 1 kg ",ha^-1," ",yr^-1),"antecedent N dep. decrease")))+
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 16))+
  annotate(
    "text",
    x = -Inf, y = Inf,
    hjust = -0.5, vjust = 1,
    label = "H1",
    size = 6,         # Adjust text size
    color = "black",    # Adjust text color
    fontface = 2
  )+
  annotate(
    "text",
    x = -Inf, y = -Inf,
    hjust = -0.5, vjust = -1,
    label = "H3",
    size = 6,         # Adjust text size
    color = "black",    # Adjust text color
    fontface = 2
  )+
  annotate(
    "text",
    x = Inf, y = -Inf,
    hjust = 1, vjust = -1,
    label = "H4",
    size = 6,         # Adjust text size
    color = "black",    # Adjust text color
    fontface = 2
  )+
  annotate(
    "text",
    x = Inf, y = Inf,
    hjust = 1, vjust = 1,
    label = "H2",
    size = 6,         # Adjust text size
    color = "black",    # Adjust text color
    fontface = 2
  )+
  scale_y_continuous(breaks = breaks_width(1))+
  scale_fill_manual(name = "",
                      labels = function(x) str_wrap(x, width = 30),
                    values = my.cols)+
  scale_color_manual(name = "",
                      labels = function(x) str_wrap(x, width = 30),
                     values = my.cols)+
  guides(color = guide_legend(override.aes = list(shape = c(16,16,16,16,16,16,16,16),
                                                  size = 3)),
         fill = "none")

fig5

ggsave(fig5,filename = "./visualizations/final_figures/Figure5.png",
       device = "png", height = 7, width = 9, units = "in")

ggsave(fig5,filename = "./visualizations/final_figures/Figure5_graphicalAbstract.png",
       device = "png", height = 4.5, width = 8, units = "in")

#### END OF CODE FOR FINAL FIGURE INCLUDED IN MANUSCRIPT


##### code for quadrant figure for CLAD presentation

plot_data2 <- plot_data %>%
  group_by(species, change_type, ecoregion) %>%
  summarize(pred = mean(pred, na.rm = TRUE)) %>%
  pivot_wider(names_from = change_type, values_from = pred) %>%
  add_column(spp_abbrevs = c("Acsa","Bepa","Litu","Piru","Pipo","Pode","Potr","Prse"))

plot_data3 <- final_pred_df %>%
  select(species, ecoregion, pred_baseline) %>%
  distinct(.) %>%
  mutate(pred_baseline = round(pred_baseline, 2))

plot_data4 <- left_join(plot_data2, plot_data3, by = c("species","ecoregion")) %>%
  mutate(antecedent = antecedent / pred_baseline * 100,
         `short-term` = `short-term` / pred_baseline * 100)


quad00 <- ggplot(data = plot_data4, aes(x = antecedent, y = `short-term`, group = species))+
  #geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  geom_point(color = "white")+
  #geom_label(aes(label = spp_abbrevs))+
  #xlim(c(-0.4, 0.4))+
  #ylim(c(-0.4, 0.4))+
  theme_classic()+
  ylab(expression(paste("% growth change w/ 1 kg ",ha^-1," ",yr^-1," short-term N dep. decrease")))+
  xlab(expression(paste("% growth change w/ 1 kg ",ha^-1," ",yr^-1," antecedent N dep. decrease")))+
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        axis.line.y = element_line(color = "white"),
        axis.text.y = element_text(color = "white"),
        axis.title.y = element_text(color = "white"),
        axis.ticks.y = element_line(color = "white"))

quad00

ggsave(quad00,filename = "./visualizations/quad_figure00.png",
       device = "png", height = 8, width = 8, units = "in")

quad0 <- ggplot(data = plot_data4, aes(x = antecedent, y = `short-term`, group = species))+
  #geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  geom_point()+
  geom_label(aes(label = spp_abbrevs, fill = ecoregion))+
  #xlim(c(-0.4, 0.4))+
  #ylim(c(-0.4, 0.4))+
  theme_classic()+
  ylab(expression(paste("% growth change w/ 1 kg ",ha^-1," ",yr^-1," short-term N dep. decrease")))+
  xlab(expression(paste("% growth change w/ 1 kg ",ha^-1," ",yr^-1," antecedent N dep. decrease")))+
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        axis.line.y = element_line(color = "white"),
        axis.text.y = element_text(color = "white"),
        axis.title.y = element_text(color = "white"),
        axis.ticks.y = element_line(color = "white"),
        legend.position = "none")

quad0

ggsave(quad0,filename = "./visualizations/quad_figure0.png",
       device = "png", height = 8, width = 8, units = "in")

quad <- ggplot(data = plot_data4, aes(x = antecedent, y = `short-term`, group = species))+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  geom_point()+
  geom_label(aes(label = spp_abbrevs, fill = ecoregion))+
  #xlim(c(-0.4, 0.4))+
  #ylim(c(-0.4, 0.4))+
  theme_classic()+
  ylab(expression(paste("% growth change w/ 1 kg ",ha^-1," ",yr^-1," short-term N dep. decrease")))+
  xlab(expression(paste("% growth change w/ 1 kg ",ha^-1," ",yr^-1," antecedent N dep. decrease")))+
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        legend.position = "none")

quad

ggsave(quad,filename = "./visualizations/quad_figure.png",
       device = "png", height = 8, width = 8, units = "in")

quad_leg_plot <- ggplot(data = plot_data4, aes(x = antecedent, y = `short-term`, group = species))+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  geom_point()+
  geom_label(aes(label = spp_abbrevs, fill = ecoregion))+
  #xlim(c(-0.4, 0.4))+
  #ylim(c(-0.4, 0.4))+
  theme_classic()+
  ylab(expression(paste("% growth change w/ 1 kg ",ha^-1," ",yr^-1," short-term N dep. decrease")))+
  xlab(expression(paste("% growth change w/ 1 kg ",ha^-1," ",yr^-1," antecedent N dep. decrease")))+
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 18))+
  scale_fill_discrete(name = "Ecoregion",
                      labels = function(x) str_wrap(x, width = 30))

quad_leg_plot

# Extract the legend. Returns a gtable
quad_leg <- get_legend(quad_leg_plot)

# Convert to a ggplot and print
quad_leg2 <- as_ggplot(quad_leg)
quad_leg2


ggsave(quad_leg2,filename = "./visualizations/quad_figure_legend.png",
       device = "png", height = 2, width = 4, units = "in")