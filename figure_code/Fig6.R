# Marginal effect of N deposition on growth
# Author: Mary Lofton
# Date: 03APR25
# Date last updated: 22DEC25

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
eco_polys <- NULL

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
  eco_polys[[i]] <- hull2poly(eco_ashape)
  
  # get grid polygons
  grid_polygons <- st_make_grid(eco_polys[[i]], 
                                cellsize = c(0.1, 0.1), # Define the size of grid cells
                                what = "polygons") |> 
    st_sf() # Convert to an sf object
  
  # filter grid cells that intersect with the ashape polygon
  grid_within_ahull <- st_intersection(grid_polygons, eco_polys[[i]])
  
  # convert the resulting grid polygons to centroids
  grid_points <- st_centroid(grid_within_ahull)
  
  # convert to dataframe
  N_grid <- data.frame(st_coordinates(grid_points)) %>%
    rename(Dep_Nhistoric_SO = X,
           Dep_Ndiff_SO = Y) %>%
    add_row(Dep_Nhistoric_SO = c(0),
            Dep_Ndiff_SO = c(0)) 
  
  # map from SO to DO
  rtilde_n = N_grid$Dep_Ndiff_SO
  atilde_n = N_grid$Dep_Nhistoric_SO
  
  n_n = length(rtilde_n)
  
  P_n = diag(n_n) - (rtilde_n %*% solve(t(rtilde_n) %*% rtilde_n, t(rtilde_n)))
  
  astar_n = as.numeric(P_n %*% atilde_n)
  
  N_grid$Dep_Nhistoric_DO_scaled = scale(astar_n)
  N_grid$Dep_Ndiff_DO_scaled = scale(rtilde_n)
  
  baseline_N_DO_values <- N_grid %>%
    filter(Dep_Ndiff_SO == 0 & Dep_Nhistoric_SO == 0)
  
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
    
    # get mean values of model parameters for that species
    params <- final1 %>% filter(spp_id == common_names[j],
                                model_id == "ortho_log_t_interaction_01JAN26")
    
    # back transformation of model parameters to SO space
    # p5 is coefficient on recent N dep; p12 is coefficient on recent N dep quadratic term
    # this is not needed if we work in DO space
    # p5_so <- params$p5 - model_data$c * params$p9
    # p12_so <- params$p12 - model_data$c * params$p11
    
    # get baseline model predictions!
    pred_log_baseline <- ((params$global_tree_effect + baseline_N_DO_values$Dep_Ndiff_DO_scaled*params$p5 + 0*params$p6 + 0*params$p7 + 0*params$p8 + baseline_N_DO_values$Dep_Nhistoric_DO_scaled*params$p9 + model_data$mean_Dep_Shistoric*params$p10 + baseline_N_DO_values$Dep_Nhistoric_DO_scaled*baseline_N_DO_values$Dep_Ndiff_DO_scaled*params$p11 + baseline_N_DO_values$Dep_Ndiff_DO_scaled*baseline_N_DO_values$Dep_Ndiff_DO_scaled*params$p12) 
                 + params$p2*model_data$log_mean_size) - params$p3*log(1)
    
    # exponentiate from log space and calculate percent increase/decrease in growth
    m2_baseline = exp(pred_log_baseline + model_data$log_mean_size)
    pred_baseline = m2_baseline - model_data$mean_size
    pred_baseline <- pred_baseline[1,1]
  
    # get other model predictions!
    pred_log <- ((params$global_tree_effect + N_grid$Dep_Ndiff_DO_scaled*params$p5 + 0*params$p6 + 0*params$p7 + 0*params$p8 + N_grid$Dep_Nhistoric_DO_scaled*params$p9 + model_data$mean_Dep_Shistoric*params$p10 + N_grid$Dep_Nhistoric_DO_scaled*N_grid$Dep_Ndiff_DO_scaled*params$p11 + N_grid$Dep_Ndiff_DO_scaled*N_grid$Dep_Ndiff_DO_scaled*params$p12) 
                        + params$p2*model_data$log_mean_size) - params$p3*log(1)
  
    # exponentiate from log space and calculate percent increase/decrease in growth
    m2 = exp(pred_log + model_data$log_mean_size)
    pred = m2 - model_data$mean_size
    
    pred_perc = (pred - pred_baseline) / pred_baseline * 100

    # assemble final df for this species-ecoregion combination
    pred_df <- data.frame(species = model_data$species[1],
                          model_id = params$model_id[1],
                          ecoregion = model_data$ecoregion_ortho[1],
                          pred = pred_perc,
                          pred_log = pred_log,
                          pred_baseline = pred_baseline,
                          N_ante_SO = N_grid$Dep_Nhistoric_SO,
                          N_rec_SO = N_grid$Dep_Ndiff_SO)
  
    if(i == 1 & j == 1){
      final_pred_df <- pred_df
    } else {
      final_pred_df <- bind_rows(final_pred_df, pred_df)
    }
  
  }
  
}

###### create figure

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

most_n_ecos <- df1 %>%
  group_by(species) %>%
  filter(n_trees_ecoregion == max(n_trees_ecoregion)) %>%
  select(species, common_name, ecoregion_ortho) %>%
  rename(ecoregion = ecoregion_ortho)

plot_data <- left_join(final_pred_df, facet_label_df, by = c("species", "ecoregion") ) %>%
  arrange(species) %>%
  right_join(., most_n_ecos, by = c("species", "ecoregion")) %>%
  mutate(bins = cut(
    pred,
    breaks = seq(-60, 40, by = 10),
    labels = c("(-59.9)-(-50)","(-49.9)-(-40)", "(-39.9)-(-30)", "(-29.9)-(-20)", "(-19.9)-(-10)",
               "(-9.9)-0","0.1-10","10.1-20","20.1-30","30.1-40"),
    right = TRUE # (0, 29] means up to and including 29
  )) %>%
  mutate(syndrome = ifelse(species %in% c("Acer saccharum","Betula papyrifera"),1,
                           ifelse(species %in% c("Liriodendron tulipifera","Prunus serotina"),2,
                                  ifelse(species == "Pinus ponderosa",3,
                                         ifelse(species %in% c("Populus deltoides","Populus tremuloides"),4,5)))))

check <- plot_data %>%
  filter(is.na(bins))

# Define number of bins
n_bins <- length(unique(plot_data$bins))

# Get hex codes for n bins
color_hex_codes <- colorRampPalette(c("#D5B60A","lightyellow","darkblue"))(n_bins) # Option "D" is default viridis
names(color_hex_codes) <- c("(-59.9)-(-50)","(-49.9)-(-40)", "(-39.9)-(-30)", "(-29.9)-(-20)", "(-19.9)-(-10)",
                            "(-9.9)-0","0.1-10","10.1-20","20.1-30","30.1-40")
show_col(color_hex_codes)
color_hex_codes

color_axis_lim_data <- left_join(most_n_ecos, plot_data, by = c("species", "ecoregion"))
color_axis_lims <- range(color_axis_lim_data$pred)

plots <- NULL

species_list <- unique(plot_data$species)

l=1

for(j in 1:length(species_list)){
  
  most_n_eco <- df1 %>%
    filter(species == species_list[j]) %>%
    filter(n_trees_ecoregion == max(n_trees_ecoregion)) %>%
    pull(ecoregion_ortho)
  
  mean_growth_eco <- df1 %>%
    filter(species == species_list[j]) %>%
    filter(n_trees_ecoregion == max(n_trees_ecoregion)) %>%
    pull(mean_perc_growth) %>%
    round(., 1)
  
  plot_dat <- plot_data %>%
    filter(species == species_list[j] & ecoregion == most_n_eco) %>%
    arrange(pred)
  
    plots[[l]] <- ggplot(data = plot_dat)+
      geom_point(aes(x = N_ante_SO, y = N_rec_SO, fill = bins, color = bins), shape = 21)+
      annotate("text",
               x = Inf, y = Inf,
               label = paste0("RS ",plot_dat$syndrome),
               hjust = 1.2, vjust = 2,
               size = 4, color = "black")+
      scale_color_manual(values = color_hex_codes)+
      scale_fill_manual(values = color_hex_codes)+
      facet_wrap(~facet_labels, scales = "free", labeller = label_parsed)+
      theme_bw()+
      labs(fill = "% change \nin growth", color = "% change \nin growth")+
      xlab(NULL)+
      ylab(NULL)+
      geom_hline(yintercept = 0, linetype = 2)+
      geom_vline(xintercept = 0, linetype = 2)+
      theme(strip.text = element_text(size = 15),
            panel.grid = element_blank(),
            strip.background = element_rect(fill = "white"),
            legend.position = "none")+
      guides(color = guide_legend(override.aes = list(size = 3)))
    l <- l + 1
  
}

# panel i guide to quadrants
plot_dat0 <- plot_dat[1,] %>%
  mutate(facet_labels = paste0("italic(i.~Plot~quadrant~guide)"))
plots[[9]] <- ggplot(data = plot_dat0)+
  geom_point(aes(x = N_ante_SO, y = N_rec_SO), fill = "white", color = "white", shape = 21)+
  facet_wrap(~facet_labels, scales = "free", labeller = label_parsed)+
  theme_bw()+
  labs(fill = "% change \nin growth", color = "% change \nin growth")+
  xlab(NULL)+
  ylab(NULL)+
  geom_hline(yintercept = 0, linetype = 2)+
  geom_vline(xintercept = 0, linetype = 2)+
  theme(strip.text = element_text(size = 15),
        panel.grid = element_blank(),
        strip.background = element_rect(fill = "white"),
        legend.position = "none")+
  guides(color = guide_legend(override.aes = list(size = 3)))+
  annotate(
    "text",
    x = -Inf, y = Inf,
    hjust = -0.05, vjust = 1.5,
    label = "antecedent N dep. higher \nthan expected & short-term \nchange in N dep. lower \nthan expected given S dep.",
    size = 2.5,         # Adjust text size
    color = "black",    # Adjust text color
    fontface = 2
  )+
  annotate(
    "text",
    x = -Inf, y = -Inf,
    hjust = -0.05, vjust = -1,
    label = "antecedent and short-term \nchange in N dep. both lower \nthan expected given S dep.",
    size = 2.5,         # Adjust text size
    color = "black",    # Adjust text color
    fontface = 2
  )+
  annotate(
    "text",
    x = Inf, y = -Inf,
    hjust = 1.05, vjust = -0.5,
    label = "antecedent N dep. higher \nthan expected & short-term \nchange in N dep. lower \nthan expected given S dep.",
    size = 2.5,         # Adjust text size
    color = "black",    # Adjust text color
    fontface = 2
  )+
  annotate(
    "text",
    x = Inf, y = Inf,
    hjust = 1.05, vjust = 1.7,
    label = "antecedent and short-term \nchange in N dep. both higher \nthan expected given S dep.",
    size = 2.5,         # Adjust text size
    color = "black",    # Adjust text color
    fontface = 2
  )+
  scale_y_continuous(limits = c(-20, 20), breaks = c(-20, -10, 0, 10, 20), labels = c("-","",0,"","+"))+
  scale_x_continuous(limits = c(-20, 20), breaks = c(-20, -10, 0, 10, 20), labels = c("-","",0,"","+"))
plots[[9]]

leg_plot <-  ggplot(data = plot_data)+
  geom_point(aes(x = N_ante_SO, y = N_rec_SO, fill = bins, color = bins), shape = 21)+
  scale_color_manual(values = color_hex_codes)+
  scale_fill_manual(values = color_hex_codes)+
  facet_wrap(~facet_labels, scales = "free", labeller = label_parsed)+
  theme_bw()+
  xlab(c)+
  ylab(expression(paste("recent N dep. change  (kg N ", ha^-1," ",yr^-1,")")))+
  labs(fill = "% change \nin growth", color = "% change \nin growth")+
  xlab(NULL)+
  ylab(NULL)+
  geom_hline(yintercept = 0, linetype = 2)+
  geom_vline(xintercept = 0, linetype = 2)+
  theme(strip.text = element_text(size = 15),
        panel.grid = element_blank(),
        strip.background = element_rect(fill = "white"),
        legend.position = "bottom",
        legend.key.size = unit(1.5, "cm"),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 16))+
  guides(color = guide_legend(override.aes = list(size = 5)))
leg_plot

# Extract the legend. Returns a gtable
fig5_leg <- cowplot::get_plot_component(leg_plot, 'guide-box-bottom', return_all = TRUE)

# Convert to a ggplot and print
fig5_leg2 <- as_ggplot(fig5_leg)
fig5_leg2


money_plot <- ggarrange(ggarrange(plotlist = plots, ncol = 3, nrow = 3))#, labels = LETTERS[1:length(plots)])

labeled_plot <- annotate_figure(
  money_plot,
  left = textGrob(expression(paste("residual short-term change in N deposition (kg N ", ha^-1," ",yr^-1,")")), rot = 90, gp = gpar(cex = 1.3)),
  bottom = textGrob(expression(paste("residual antecedent N deposition (kg N ", ha^-1," ",yr^-1,")")), gp = gpar(cex = 1.3))
)

final_plot <- ggarrange(labeled_plot,
                        fig5_leg2, ncol = 1, nrow = 2, heights = c(7,1))

ggsave(final_plot,filename = "./visualizations/final_figures/Figure6.png",
       device = "png", bg = "white", height = 10, width = 11, units = "in")
