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

###### different attempt at figure

letters_df <- final_pred_df %>%
  select(species, ecoregion) %>%
  distinct(.) %>%
  arrange(species, ecoregion) %>%
  mutate(letters = paste0(c(letters[seq_len(26)],"aa","bb","cc"),"."))

facet_label_df <- final_pred_df %>%
  select(species, ecoregion, pred_baseline) %>%
  distinct(.) %>%
  left_join(., letters_df, by = c("species","ecoregion")) %>%
  mutate(pred_baseline = round(pred_baseline, 2)) %>%
  separate_wider_delim(species, delim = " ", names = c("genus","spp"), cols_remove = FALSE) %>%
  separate_wider_delim(ecoregion, delim = "&", names = c("er1","er2"), too_few = "align_start",cols_remove = FALSE) %>%
  mutate(er1 = str_replace_all(er1, " ", "~"),
         er2 = str_replace_all(er2, " ", "~"),
         er1 = str_remove(er1, "~$"),
         er2 = str_remove(er2, "^~"),
         facet_labels = paste0("atop(italic(",letters,"~",genus,"~",spp,"), atop(atop(",er1, ",",er2,"),expected~growth:~",pred_baseline,"~kg~C~y^-1~ind^-1))")) %>%
  select(species, ecoregion, facet_labels)

plot_data <- left_join(final_pred_df, facet_label_df, by = c("species", "ecoregion") ) %>%
  arrange(species) %>%
  mutate(bins = cut(
    pred,
    breaks = seq(-60, 70, by = 10),
    labels = c("(-59.9)-(-50)","(-49.9)-(-40)", "(-39.9)-(-30)", "(-29.9)-(-20)", "(-19.9)-(-10)",
               "(-9.9)-0","0.1-10","10.1-20","20.1-30","30.1-40","40.1-50","50.1-60","60.1-70"),
    right = TRUE # (0, 29] means up to and including 29
  ))

check <- plot_data %>%
  filter(is.na(bins))

# Define number of bins
n_bins <- length(unique(plot_data$bins))

# Get hex codes for n bins
color_hex_codes <- colorRampPalette(c("#D5B60A","lightyellow","darkblue"))(n_bins) # Option "D" is default viridis
names(color_hex_codes) <- c("(-59.9)-(-50)","(-49.9)-(-40)", "(-39.9)-(-30)", "(-29.9)-(-20)", "(-19.9)-(-10)",
                            "(-9.9)-0","0.1-10","10.1-20","20.1-30","30.1-40","40.1-50","50.1-60","60.1-70")
show_col(color_hex_codes)
color_hex_codes

plots <- NULL

species_list <- unique(plot_data$species)

l=1

for(j in 1:length(species_list)){
  
  plot_dat <- plot_data %>%
    filter(species == species_list[j])
  
  ecoreg <- unique(plot_dat$ecoregion)
  
  for(k in 1:length(ecoreg)){
    
    plot_dat2 <- plot_dat %>%
      filter(ecoregion == ecoreg[k]) %>%
      arrange(pred)
    
    plots[[l]] <- ggplot(data = plot_dat2)+
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
            legend.position = "none")+
      guides(color = guide_legend(override.aes = list(size = 3)))
    l <- l + 1
  }
  
}

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
        strip.background = element_rect(fill = "white"))+
  guides(color = guide_legend(override.aes = list(size = 3)))
leg_plot

# Extract the legend. Returns a gtable
figs19_leg <- get_legend(leg_plot)

# Convert to a ggplot and print
figs19_leg2 <- as_ggplot(figs19_leg)
figs19_leg2

money_plot <- ggarrange(ggarrange(plotlist = plots, ncol = 6, nrow = 5),
                        figs19_leg2, ncol = 2, nrow = 1, widths = c(7,1))#, labels = LETTERS[1:length(plots)])

final_plot <- annotate_figure(
  money_plot,
  left = textGrob(expression(paste("residual recent N deposition (kg N ", ha^-1," ",yr^-1,")")), rot = 90, gp = gpar(cex = 1.3)),
  bottom = textGrob(expression(paste("residual antecedent N deposition (kg N ", ha^-1," ",yr^-1,")")), gp = gpar(cex = 1.3))
)

ggsave(final_plot,filename = "./visualizations/final_figures/FigureS27.png",
       device = "png", bg = "white", height = 12, width = 22, units = "in")



##################################


# get weighted average of predictions
df2 <- df1 %>%
  group_by(species, common_name) %>%
  summarize(n_trees_total = sum(n_trees_ecoregion))

spp <- unique(df1$species)

for(i in 1:length(spp)){
  
  n_trees_total <- df2 %>%
    filter(species == spp[i]) %>%
    pull(n_trees_total)
  
  n_trees_ecoregion <- df1 %>%
    filter(species == spp[i]) %>%
    select(ecoregion_ortho, n_trees_ecoregion)
  
  tree_ecos <- unique(n_trees_ecoregion$ecoregion_ortho)
  
  for(j in 1:length(tree_ecos)){ 
    
    species_pred_eco <- final_pred_df %>%
      filter(species == spp[i] & ecoregion == tree_ecos[j]) %>%
      pull(pred)
    
    n_trees_this_eco <- n_trees_ecoregion %>%
      filter(ecoregion_ortho == tree_ecos[j]) %>%
      pull(n_trees_ecoregion)
    
    adj_pred <- species_pred_eco*n_trees_this_eco / n_trees_total

    if(j == 1){
      pred <- adj_pred
    } else {
      pred <- pred + adj_pred
    }
    
  }
  
  spp_N_ranges <- final_pred_df %>%
    filter(species == spp[i] & ecoregion == tree_ecos[j])
  
  final_pred_species <- data.frame(species = spp[i],
                                   pred = pred,
                                   N_ante_SO = spp_N_ranges$N_ante_SO,
                                   N_rec_SO = spp_N_ranges$N_rec_SO)
  
  if(i == 1){
    final_pred <- final_pred_species
  } else {
    final_pred <- bind_rows(final_pred, final_pred_species)
  }
  
  
  
}

plot_labels <- spp_df %>%
  dplyr::filter(!common_name %in% c("Douglas-fir","western hemlock")) %>%
  arrange(species) %>%
  mutate(labels = paste0(letters[seq_len(length(unique(species)))],". ")) 

plot_data <- left_join(final_pred, plot_labels, by = "species") %>%
  mutate(final_labels = paste0(labels, species))

mean_levels <- df %>%
  select(common_name, Dep_Ndiff, Dep_Nhistoric, plot_ID) %>%
  distinct(.) %>%
  group_by(common_name) %>%
  summarize(Dep_Ndiff_mean = mean(Dep_Ndiff, na.rm = TRUE),
            Dep_Nhistoric_mean = mean(Dep_Nhistoric, na.rm = TRUE)) %>%
  left_join(spp_df, by = "common_name")
  
  plots <- NULL
  
  species_list <- unique(final_pred$species)
  
  for(j in 1:length(species_list)){
    
    plot_dat <- plot_data %>%
      filter(species == species_list[j]) 
    
    spp_mean_levels <- mean_levels %>%
      filter(species == species_list[j])
    
    plots[[j]] <- ggplot()+
      geom_tile(data = plot_dat, aes(x = N_ante, y = N_rec, fill = pred))+
      geom_contour(data = plot_dat, aes(x = N_ante, y = N_rec, z = pred), color = "lightgray", alpha = 0.3)+
      metR::geom_text_contour(data = plot_dat, aes(x = N_ante, y = N_rec, z = pred),
                              label.placer = label_placer_fraction(frac = 0.5),
                              size = 3,
                              color = "lightgray")+
      scale_fill_viridis_c()+
      #geom_hline(yintercept = spp_mean_levels$Dep_Ndiff_mean, color = "white",linetype = 2)+
      #geom_vline(xintercept = spp_mean_levels$Dep_Nhistoric_mean, color = "white",linetype = 2)+
      ggtitle(plot_dat$final_labels[1])+
      theme_bw()+
      xlab(expression(paste("antecedent N deposition (kg N ", ha^-1," ",yr^-1,")")))+
      ylab(expression(paste("recent N dep. change  (kg N ", ha^-1," ",yr^-1,")")))+
      labs(fill = "growth (%)")
    
  }
  
  money_plot <- ggarrange(plotlist = plots, ncol = 3, nrow = 3)#, labels = LETTERS[1:length(plots)])

  
  ggsave(money_plot,filename = "./visualizations/final_figures/Figure6_test2.png",
         device = "png", height = 9, width = 12, units = "in", bg = "white")
  


