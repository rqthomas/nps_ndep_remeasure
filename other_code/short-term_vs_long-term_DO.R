# Marginal effect of N deposition on growth
# Author: Mary Lofton
# Date: 03APR25

# Purpose: figure with marginal effect of N deposition (decrease of 1 kg
# per hectare per year) on growth (kg C per year per individual)

library(scales)

data = "./data/McDonnell_etal_InPrep_TreeData_2024_10_11.csv"

marginal_effect_Ndep_antecedent_recent <- function(data){
  
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
  select(model_id, spp_id, global_tree_effect, p2, p3, p5, p6, p7, p8, p9, p10) 

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

final_pred_shortterm <- NULL

for(i in 1:length(common_names)){
  model_data <- df1 %>% filter(common_name == common_names[i])
  model_output <- final1 %>% filter(spp_id == common_names[i],
                              model_id == "space_vs_time_ortho")
  params <- model_output[sample(nrow(model_output), 1000), ]
  
  for(j in 1:length(model_data$ecoregion_ortho)){
  
  p5_og <- params$p5 - model_data$c[j] * params$p9
  
  pred_baseline <- ((params$global_tree_effect + 0*p5_og + 0*params$p6 + 0*params$p7 + 0*params$p8 + model_data$mean_Dep_Nhistoric[j]*params$p9 + model_data$mean_Dep_Shistoric[j]*params$p10) 
           * model_data$mean_size[j] ^ params$p2) * 1 #exp(-model_data$mean_ba_gt*params$p3) 
  
  pred_decrease <- ((params$global_tree_effect + 1*p5_og + 0*params$p6 + 0*params$p7 + 0*params$p8 + model_data$mean_Dep_Nhistoric[j]*params$p9 + model_data$mean_Dep_Shistoric[j]*params$p10) 
                    * model_data$mean_size[j] ^ params$p2) * 1 #exp(-model_data$mean_ba_gt*params$p3) 
  
  growth_delta <- pred_decrease - pred_baseline
  
  pred <- data.frame(species = model_data$species[1],
                     model_id = model_output$model_id[1],
                     ecoregion = model_data$ecoregion_ortho[j],
                     pred = growth_delta)
  
  if(i == 1 & j == 1){
    final_pred_shortterm <- pred
  } else {
    final_pred_shortterm <- bind_rows(final_pred_shortterm, pred)
  }
  
  }
  
}


# get model predictions for antecedent Ndep changes 

final_pred_longterm <- NULL

for(i in 1:length(common_names)){
  model_data <- df1 %>% filter(common_name == common_names[i])
  model_output <- final1 %>% filter(spp_id == common_names[i],
                                    model_id == "space_vs_time_ortho")
  params <- model_output[sample(nrow(model_output), 1000), ]
  
  for(j in 1:length(model_data$ecoregion_ortho)){
    
    p5_og <- params$p5 - model_data$c[j] * params$p9
    
    pred_baseline <- ((params$global_tree_effect + 0*p5_og + 0*params$p6 + 0*params$p7 + 0*params$p8 + model_data$mean_Dep_Nhistoric[j]*params$p9 + model_data$mean_Dep_Shistoric[j]*params$p10) 
                      * model_data$mean_size[j] ^ params$p2) * 1 #exp(-model_data$mean_ba_gt*params$p3) 
    
    pred_decrease <- ((params$global_tree_effect + 0*p5_og + 0*params$p6 + 0*params$p7 + 0*params$p8 + (model_data$mean_Dep_Nhistoric[j]+1)*params$p9 + model_data$mean_Dep_Shistoric[j]*params$p10) 
                      * model_data$mean_size[j] ^ params$p2) * 1 #exp(-model_data$mean_ba_gt*params$p3) 
    
    growth_delta <- pred_decrease - pred_baseline
    
    pred <- data.frame(species = model_data$species[1],
                       model_id = model_output$model_id[1],
                       ecoregion = model_data$ecoregion_ortho[j],
                       pred = growth_delta)
    
    if(i == 1 & j == 1){
      final_pred_longterm <- pred
    } else {
      final_pred_longterm <- bind_rows(final_pred_longterm, pred)
    }
    
  }
  
}

final_pred_longterm <- final_pred_longterm %>%
  rename(longterm = pred)
final_pred_shortterm <- final_pred_shortterm %>%
  rename(shortterm = pred)

final_pred <- bind_cols(final_pred_longterm, final_pred_shortterm$shortterm) 
colnames(final_pred)[5] <- "shortterm"

ecoregions <- unique(df$ecoregion_ortho)

for(i in 1:length(ecoregions)){
  
  plot_pred <- final_pred %>%
    filter(ecoregion == ecoregions[i]) 
  
  cb_palette <- colorblind_pal()(8)
  
  p <- ggplot(data = plot_pred)+
    geom_hline(yintercept = 0)+
    geom_vline(xintercept = 0)+
    geom_point(aes(x = longterm, y = shortterm, col = species), size = 0.1, alpha = 0.1)+
    stat_ellipse(aes(x = longterm, y = shortterm, col = species), level = 0.95, alpha = 0.5)+
    theme_bw()+
    xlim(-2,2)+
    ylim(-2,2)+
    scale_color_manual(values = c("Acer saccharum" = cb_palette[1],
                                  "Betula papyrifera" = cb_palette[2],
                                  "Liriodendron tulipifera" = cb_palette[3],
                                  "Picea rubens" = cb_palette[4],
                                  "Pinus ponderosa" = cb_palette[5],
                                  "Populus deltoides" = cb_palette[6],
                                  "Populus tremuloides" = cb_palette[7],
                                  "Prunus serotina" = cb_palette[8]))+ # Adjust width as needed
    ylab("Growth change given 1 kg ha^-1 yr^-1 N 'short-term' increase")+
    xlab("Growth change given 1 kg ha^-1 yr^-1 N 'long-term' increase")+
    labs(color = "Species", title = ecoregions[i])
  
  ggsave(p,filename = paste0("./visualizations/short-term_vs_long-term_DO/ecoregions/",ecoregions[i],".png"),
         device = "png")
  
  
}

spp <- unique(df1$species)

for(i in 1:length(spp)){
  
  plot_pred <- final_pred %>%
    filter(species == spp[i])
  
  cb_palette <- hue_pal()(10)
  
  p <- ggplot(data = plot_pred)+
    geom_hline(yintercept = 0)+
    geom_vline(xintercept = 0)+
    geom_point(aes(x = longterm, y = shortterm, col = ecoregion), size = 0.1, alpha = 0.1)+
    stat_ellipse(aes(x = longterm, y = shortterm, col = ecoregion), level = 0.95, alpha = 0.5)+
    theme_bw()+
    xlim(-2,2)+
    ylim(-2,2)+
    scale_color_manual(values = c("EASTERN TEMPERATE FORESTS" = cb_palette[1],
                                  "NORTHWESTERN FORESTED MOUNTAINS & MARINE WEST COAST FOREST" = cb_palette[2],
                                  "MEDITERRANEAN CALIFORNIA" = cb_palette[3],
                                  "NORTH AMERICAN DESERTS" = cb_palette[4],
                                  "GREAT PLAINS" = cb_palette[5],
                                  "NORTHERN FORESTS" = cb_palette[6],
                                  "TEMPERATE SIERRAS & SOUTHERN SEMIARID HIGHLANDS" = cb_palette[7]), labels = function(x) str_wrap(x, width = 15))+ # Adjust width as needed
    ylab("Growth change given 1 kg ha^-1 yr^-1 N 'short-term' increase")+
    xlab("Growth change given 1 kg ha^-1 yr^-1 N 'long-term' increase")+
    labs(color = "Ecoregion", title = spp[i])
  
  ggsave(p,filename = paste0("./visualizations/short-term_vs_long-term_DO/species/",spp[i],".png"),
         device = "png")
  
  
}

}

