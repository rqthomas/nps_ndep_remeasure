# Marginal effect of N deposition on growth
# Author: Mary Lofton
# Date: 03APR25

# Purpose: figure with marginal effect of N deposition (decrease of 1 kg
# per hectare per year) on growth (kg C per year per individual)

source("./other_code/get_model_inputs.R")

marginal_effect_Ndep <- function(data){
  
# list files in each model output folder
out1 <- list.files("./experiments/new_delta_Ndep_only_saveTreeEffect",pattern = "mcmc.parquet",
                   full.names = TRUE)
out2 <- list.files("./experiments/delta_env_saveTreeEffect",pattern = "mcmc.parquet",
                   full.names = TRUE)
out3 <- list.files("./experiments/N_species_saveTreeEffect",pattern = "mcmc.parquet",
                   full.names = TRUE)
out <- c(out1,out2, out3)


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
  select(model_id, spp_id, global_tree_effect, p2, p3, p4, p5, p6, p7, p8) 

df <- get_model_inputs(data)

df1 <- df %>%
  group_by(species, common_name) %>%
  summarize(mean_size = mean(AG_carbon_m1, na.rm = TRUE),
            mean_Dep_Sdelta = mean(Dep_Sdelta, na.rm = TRUE),
            mean_MAT_delta = mean(MAT_delta, na.rm = TRUE),
            mean_MAP_delta_dm = mean(MAP_delta_dm, na.rm = TRUE),
            mean_ba_gt = mean(subp_BA_GT_m1, na.rm = TRUE)) %>%
  mutate(common_name = ifelse(common_name == "yellow-poplar","yellow poplar",common_name)) %>%
  ungroup()

# get model predictions for delta_env

common_names <- unique(df1$common_name)

final_pred_env <- NULL

for(i in 1:length(common_names)){
  model_data <- df1 %>% filter(common_name == common_names[i])
  model_output <- final1 %>% filter(spp_id == common_names[i],
                              model_id == "delta_env_saveTreeEffect")
  params <- model_output[sample(nrow(model_output), 1000), ]
  
  pred_baseline <- ((params$global_tree_effect + 0*params$p5 + 0*params$p6 + 0*params$p7 + 0*params$p8) 
           * model_data$mean_size ^ params$p2) * 1 #exp(-model_data$mean_ba_gt*params$p3) 
  
  pred_decrease <- ((params$global_tree_effect + -1*params$p5 + 0*params$p6 + 0*params$p7 + 0*params$p8) 
                    * model_data$mean_size ^ params$p2) * 1 #exp(-model_data$mean_ba_gt*params$p3) 
  
  growth_delta <- pred_decrease - pred_baseline
  
  pred <- data.frame(species = model_data$species[1],
                     model_id = model_output$model_id[1],
                     pred = growth_delta)
  
  if(i == 1){
    final_pred_env <- pred
  } else {
    final_pred_env <- bind_rows(final_pred_env, pred)
  }
}


# get model predictions for delta Ndep only 

final_pred_Ndep_only <- NULL

for(i in 1:length(common_names)){
  model_data <- df1 %>% filter(common_name == common_names[i])
  model_output <- final1 %>% filter(spp_id == common_names[i],
                                    model_id == "new_delta_Ndep_only_saveTreeEffect")
  params <- model_output[sample(nrow(model_output), 1000), ]
  
  pred_baseline <- ((params$global_tree_effect + 0*params$p5) 
                    * model_data$mean_size ^ params$p2) * 1 #exp(-model_data$mean_ba_gt*params$p3) 
  
  pred_decrease <- ((params$global_tree_effect + -1*params$p5) 
                    * model_data$mean_size ^ params$p2) * 1 #exp(-model_data$mean_ba_gt*params$p3) 
  
  growth_delta <- pred_decrease - pred_baseline

  pred <- data.frame(species = model_data$species[1],
                     model_id = model_output$model_id[1],
                     pred = growth_delta)
  
  if(i == 1){
    final_pred_Ndep_only <- pred
  } else {
    final_pred_Ndep_only <- bind_rows(final_pred_Ndep_only, pred)
  }
}

final_pred <- bind_rows(final_pred_env, final_pred_Ndep_only) %>%
  mutate(model_label = ifelse(model_id == "delta_env_saveTreeEffect","N dep. + other env. variables","N dep. only"))

plot_pred <- final_pred %>%
  filter(model_label == "N dep. + other env. variables")

my_col <- RColorBrewer::brewer.pal(3, "Blues")[3]

p <- ggplot(data = plot_pred)+
  geom_density(aes(x = pred, group = model_label, 
                   color = model_label, fill = model_label),
               alpha = 0.5)+
  facet_wrap(facets = vars(species), scales = "free")+
  theme_bw()+
  scale_color_manual(values = my_col)+
  scale_fill_manual(values = my_col)+
  geom_vline(xintercept = 0)+
  labs(color = "", fill = "", x = expression(paste("change in growth (kg C ", y^-1," ",ind^-1,")")))+
  ggtitle(expression(paste("Marginal effect of N deposition decrease (per kg N ", ha^-1," ",y^-1,")")))

return(p)

}

