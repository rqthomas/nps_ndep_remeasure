# Marginal effect of N deposition on growth
# Author: Mary Lofton
# Date: 03APR25

# Purpose: figure with marginal effect of N deposition (decrease of 1 kg
# per hectare per year) on growth (kg C per year per individual)

library(scales)

data = "./data/McDonnell_etal_InPrep_TreeData_2024_10_11.csv"

growth_change_baseline_DO <- function(data){
  
# list files in each model output folder
out <- list.files("./experiments/space_vs_time_ortho",pattern = "mcmc.parquet",
                   full.names = TRUE)

# read in and combine files
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

final1 <- final %>%
  mutate(spp_id = ifelse(spp_id == "yellow","yellow poplar",spp_id)) %>%
  select(model_id, spp_id, global_tree_effect, p2, p3, p5, p6, p7, p8, p9, p10) %>%
  group_by(model_id, spp_id) %>%
  summarise(across(global_tree_effect:p10, \(x) mean(x, na.rm = TRUE)))

df <- read_csv("./data/processed_data.csv")

spp_df <- read_csv(data) %>%
  select(common_name, species) %>%
  distinct(.)

df1 <- left_join(df, spp_df, by = "common_name") %>%
  group_by(species, common_name, ecoregion_ortho) %>%
  summarize(mean_size = mean(AG_carbon_m1, na.rm = TRUE),
            mean_Dep_Nhistoric = mean(Dep_Nhistoric, na.rm = TRUE),
            mean_Dep_Shistoric = mean(Dep_Shistoric, na.rm = TRUE),
            c = mean(c, na.rm = TRUE)) %>%
  mutate(common_name = ifelse(common_name == "yellow-poplar","yellow poplar",common_name)) %>%
  ungroup()

# get model predictions for recent N dep changes

common_names <- unique(df1$common_name)

final_pred_df <- NULL

N_ranges <- NULL

for(i in 1:length(common_names)){
  
  model_data <- df1 %>% filter(common_name == common_names[i])
  params <- final1 %>% filter(spp_id == common_names[i],
                              model_id == "space_vs_time_ortho")

  for(j in 1:length(model_data$ecoregion_ortho)){
    
    curr_dat <- df %>%
      mutate(common_name = ifelse(common_name == "yellow-poplar","yellow poplar",common_name)) %>%
      filter(common_name == common_names[i] & ecoregion_ortho == model_data$ecoregion_ortho[j])
    
    if(i == 2 & j == 3) next
    
    min_Nhistoric <- min(curr_dat$Dep_Nhistoric, na.rm = TRUE)
    max_Nhistoric <- max(curr_dat$Dep_Nhistoric, na.rm = TRUE)
    min_Ndiff <- min(curr_dat$Dep_Ndiff, na.rm = TRUE)
    max_Ndiff <- max(curr_dat$Dep_Ndiff, na.rm = TRUE)
    
    ante_Ndep_values = seq(min_Nhistoric, max_Nhistoric, by = 0.1)
    rec_Ndep_values = seq(min_Ndiff, max_Ndiff, by = 0.1)
    
    N_ranges <- data.frame(ante_Ndep = rep(ante_Ndep_values, times = length(rec_Ndep_values)),
                           rec_Ndep = rep(rec_Ndep_values, each = length(ante_Ndep_values)),
                           species = model_data$species[1],
                           spp_id = model_data$common_name[1]) %>%
      mutate(spp_id = ifelse(spp_id == "yellow-poplar","yellow poplar",spp_id)) %>%
      filter(ante_Ndep + rec_Ndep >= 0)
  
  p5_og <- params$p5 - model_data$c[j] * params$p9
  
  pred <- ((params$global_tree_effect + N_ranges$rec_Ndep*p5_og + 0*params$p6 + 0*params$p7 + 0*params$p8 + N_ranges$ante_Ndep*params$p9 + model_data$mean_Dep_Shistoric[j]*params$p10) 
           * model_data$mean_size[j] ^ params$p2) * 1 #exp(-model_data$mean_ba_gt*params$p3) 
  
  pred_perc = (pred / model_data$mean_size[j]) * 100
  
  pred_df <- data.frame(species = model_data$species[1],
                     model_id = params$model_id[1],
                     ecoregion = model_data$ecoregion_ortho[j],
                     pred = pred_perc,
                     N_ante = N_ranges$ante_Ndep,
                     N_rec = N_ranges$rec_Ndep)
  
  if(i == 1 & j == 1){
    final_pred_df <- pred_df
  } else {
    final_pred_df <- bind_rows(final_pred_df, pred_df)
  }
  
  }
  
}

ecoregions <- unique(df$ecoregion_ortho)
spp <- unique(df1$species)

for(i in 1:length(ecoregions)){
  
  plot_pred <- final_pred_df %>%
    filter(ecoregion == ecoregions[i]) 
  
  plots <- NULL
  
  ecoreg_species <- unique(plot_pred$species)
  
  for(j in 1:length(ecoreg_species)){
    
    plot_dat <- plot_pred %>%
      filter(species == ecoreg_species[j]) 
    
    plots[[j]] <- ggplot(data = plot_dat)+
      geom_tile(aes(x = N_ante, y = N_rec, fill = pred))+
      scale_fill_viridis_c()+
      ggtitle(ecoreg_species[j])+
      theme_bw()+
      xlab(expression(paste("antecedent N deposition (kg N ", ha^-1," ",yr^-1,")")))+
      ylab(expression(paste("N dep. change  (kg N ", ha^-1," ",yr^-1,")")))+
      labs(fill = "growth (%)")
    
  }
  
  money_plot <- ggarrange(plotlist = plots, ncol = 2, nrow = ceiling(length(plots) / 2))#, labels = LETTERS[1:length(plots)])
  money_plot <- annotate_figure(money_plot, top = text_grob(ecoregions[i], face = "bold", size = 16))
  
  
  ggsave(money_plot,filename = paste0("./visualizations/growth_change_baseline_DO/ecoregions/",ecoregions[i],".png"),
         device = "png", scale = 2)
  
  
}

for(i in 1:length(spp)){
  
  plot_pred <- final_pred_df %>%
    filter(species == spp[i]) 
  
  plots <- NULL
  
  spp_ecoreg <- unique(plot_pred$ecoregion)
  
  for(j in 1:length(spp_ecoreg)){
    
    plot_dat <- plot_pred %>%
      filter(ecoregion == spp_ecoreg[j]) 
    
    plots[[j]] <- ggplot(data = plot_dat)+
      geom_tile(aes(x = N_ante, y = N_rec, fill = pred))+
      scale_fill_viridis_c()+
      ggtitle(spp_ecoreg[j])+
      theme_bw()+
      xlab(expression(paste("antecedent N deposition (kg N ", ha^-1," ",yr^-1,")")))+
      ylab(expression(paste("N dep. change  (kg N ", ha^-1," ",yr^-1,")")))+
      labs(fill = "growth (%)")
    
  }
  
  money_plot <- ggarrange(plotlist = plots, ncol = 2, nrow = ceiling(length(plots) / 2))#, labels = LETTERS[1:length(plots)])
  money_plot <- annotate_figure(money_plot, top = text_grob(spp[i], face = "bold", size = 16))
  
  
  ggsave(money_plot,filename = paste0("./visualizations/growth_change_baseline_DO/species/",spp[i],".png"),
         device = "png", scale = 2)
  
  
}

}

