# Marginal effect of N deposition on growth
# Author: Mary Lofton
# Date: 03APR25

# Purpose: figure with marginal effect of N deposition (decrease of 1 kg
# per hectare per year) on growth (kg C per year per individual)

source("./other_code/get_model_inputs.R")

marginal_effect_Ndep_longterm_shortterm <- function(data){
  
# list files in each model output folder
out <- list.files("./experiments/short-term_long-term",pattern = "mcmc.parquet",
                   full.names = TRUE)

# read in and combine files
for(i in 1:length(out)){
  
  spp_id = str_split(out[i], pattern = "-")[[1]][6]
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

df <- get_model_inputs(data)

df1 <- df %>%
  group_by(species, common_name) %>%
  summarize(mean_size = mean(AG_carbon_m1, na.rm = TRUE),
            mean_Dep_Nhistoric = mean(Dep_Nhistoric, na.rm = TRUE)) %>%
  mutate(common_name = ifelse(common_name == "yellow-poplar","yellow poplar",common_name)) %>%
  ungroup()

# get model predictions for longterm N dep changes

common_names <- unique(df1$common_name)

final_pred_longterm <- NULL

for(i in 1:length(common_names)){
  model_data <- df1 %>% filter(common_name == common_names[i])
  model_output <- final1 %>% filter(spp_id == common_names[i],
                              model_id == "short-term_long-term")
  params <- model_output[sample(nrow(model_output), 1000), ]
  
  pred_baseline <- ((params$global_tree_effect + 0*params$p5 + 0*params$p6 + 0*params$p7 + 0*params$p8 + 0*params$p9 + model_data$mean_Dep_Nhistoric*params$p10) 
           * model_data$mean_size ^ params$p2) * 1 #exp(-model_data$mean_ba_gt*params$p3) 
  
  pred_decrease <- ((params$global_tree_effect + 0*params$p5 + 0*params$p6 + 0*params$p7 + 0*params$p8 + -1*params$p9 + model_data$mean_Dep_Nhistoric*params$p10) 
                    * model_data$mean_size ^ params$p2) * 1 #exp(-model_data$mean_ba_gt*params$p3) 
  
  growth_delta <- pred_decrease - pred_baseline
  
  pred <- data.frame(species = model_data$species[1],
                     model_id = model_output$model_id[1],
                     pred = growth_delta)
  
  if(i == 1){
    final_pred_longterm <- pred
  } else {
    final_pred_longterm <- bind_rows(final_pred_longterm, pred)
  }
}


# get model predictions for short-term Ndep changes 

final_pred_shortterm <- NULL

for(i in 1:length(common_names)){
  model_data <- df1 %>% filter(common_name == common_names[i])
  model_output <- final1 %>% filter(spp_id == common_names[i],
                                    model_id == "short-term_long-term")
  params <- model_output[sample(nrow(model_output), 1000), ]
  
  pred_baseline <- ((params$global_tree_effect + 0*params$p5 + 0*params$p6 + 0*params$p7 + 0*params$p8 + 0*params$p9 + model_data$mean_Dep_Nhistoric*params$p10) 
                    * model_data$mean_size ^ params$p2) * 1 #exp(-model_data$mean_ba_gt*params$p3) 
  
  pred_decrease <- ((params$global_tree_effect + -1*params$p5 + 0*params$p6 + 0*params$p7 + 0*params$p8 + 0*params$p9 + model_data$mean_Dep_Nhistoric*params$p10) 
                    * model_data$mean_size ^ params$p2) * 1 #exp(-model_data$mean_ba_gt*params$p3) 
  
  growth_delta <- pred_decrease - pred_baseline

  pred <- data.frame(species = model_data$species[1],
                     model_id = model_output$model_id[1],
                     pred = growth_delta)
  
  if(i == 1){
    final_pred_shortterm <- pred
  } else {
    final_pred_shortterm <- bind_rows(final_pred_shortterm, pred)
  }
}

final_pred_longterm$change_type <- "long-term"
final_pred_shortterm$change_type <- "short-term"


final_pred <- bind_rows(final_pred_longterm, final_pred_shortterm) 

plot_pred <- final_pred 

my_col <- RColorBrewer::brewer.pal(3, "Blues")[2:3]

p <- ggplot(data = plot_pred)+
  geom_density(aes(x = pred, group = change_type, 
                   color = change_type, fill = change_type),
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

