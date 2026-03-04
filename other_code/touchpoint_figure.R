# Tree growth vs antecedent N deposition
# Author: Mary Lofton
# Date: 30OCT25

# Purpose: "touchpoint" figure showing link between Horn et al and current approach
# x axis is antecedent N deposition, y axis is tree growth, show curves with no 
# short-term change, and plus/minus 1 kg/ha/yr N dep change

# use sugar maple in eastern temperate forest as an example

# load packages
library(tidyverse)
library(lubridate)
library(arrow)
library(ggpubr)
library(ggthemes)
library(stringr)

data = "./data/McDonnell_etal_InPrep_TreeData_2024_10_11.csv"


  # list files in each model output folder
  out <- list.files("./experiments/ortho_log_t_interaction_adj_priors",pattern = "mcmc.parquet",
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
    select(model_id, spp_id, global_tree_effect, p2, p3, p5, p6, p7, p8, p9, p10, p11, p12) 
  
  df <- read_csv("./data/processed_data.csv")
  
  spp_df <- read_csv(data) %>%
    select(common_name, species) %>%
    distinct(.)
  
  df1 <- left_join(df, spp_df, by = "common_name") %>%
    group_by(species, common_name, ecoregion_ortho) %>%
    summarize(log_mean_size = log(mean(AG_carbon_m1, na.rm = TRUE)),
              log_mean_ba_gt = log(mean(subp_BA_GT_m1, na.rm = TRUE) + 1),
              min_Dep_Nhistoric = min(Dep_Nhistoric, na.rm = TRUE),
              max_Dep_Nhistoric = max(Dep_Nhistoric, na.rm = TRUE),
              mean_Dep_Shistoric = mean(Dep_Shistoric, na.rm = TRUE),
              c = mean(c, na.rm = TRUE),
              mean_size = mean(AG_carbon_m1, na.rm = TRUE),
              mean_growth = mean(AG_carbon_pYear, na.rm = TRUE),
              n_trees_ecoregion = n_distinct(tree_ID)) %>%
    mutate(common_name = ifelse(common_name == "yellow-poplar","yellow poplar",common_name)) %>%
    ungroup()
  
  # get model predictions for recent N dep changes
  
  common_names <- unique(df1$common_name)
  
  final_pred <- NULL
  
    model_data <- df1 %>% filter(common_name == "sugar maple")
    model_output <- final1 %>% filter(spp_id == "sugar maple",
                                      model_id == "ortho_log_t_interaction_adj_priors")
    params <- model_output %>%
      select(global_tree_effect:p12) %>%
      summarize(across(global_tree_effect:p12, ~ mean(.x, na.rm = TRUE)))
    
      p5_og <- params$p5 - model_data$c[1] * params$p9
      
      ante_N_dep_values <- seq(model_data$min_Dep_Nhistoric[1], model_data$max_Dep_Nhistoric[1], 0.1)
      
      pred_baseline_log <- ((params$global_tree_effect + 0*p5_og + 0*params$p6 + 0*params$p7 + 0*params$p8 + ante_N_dep_values*params$p9 + model_data$mean_Dep_Shistoric[1]*params$p10 + ante_N_dep_values*0*params$p11 + 0*0*params$p12) 
                            + params$p2*model_data$log_mean_size[1]) - params$p3*log(1)
      m2_baseline = exp(pred_baseline_log + model_data$log_mean_size[1])
      pred_baseline = m2_baseline - model_data$mean_size[1]
      
      pred_increase_log <- ((params$global_tree_effect + 1*p5_og + 0*params$p6 + 0*params$p7 + 0*params$p8 + ante_N_dep_values*params$p9 + model_data$mean_Dep_Shistoric[1]*params$p10 + ante_N_dep_values*1*params$p11 + 1*1*params$p12) 
                            + params$p2*model_data$log_mean_size[1]) - params$p3*log(1)
      m2_increase = exp(pred_increase_log + model_data$log_mean_size[1])
      pred_increase = m2_increase - model_data$mean_size[1]
      
      pred_decrease_log <- ((params$global_tree_effect + -1*p5_og + 0*params$p6 + 0*params$p7 + 0*params$p8 + ante_N_dep_values*params$p9 + model_data$mean_Dep_Shistoric[1]*params$p10 + ante_N_dep_values*-1*params$p11 + -1*-1*params$p12) 
                            + params$p2*model_data$log_mean_size[1]) - params$p3*log(1)
      m2_decrease = exp(pred_decrease_log + model_data$log_mean_size[1])
      pred_decrease = m2_decrease - model_data$mean_size[1]
      
      pred <- data.frame(species = model_data$species[1],
                         model_id = model_output$model_id[1],
                         ecoregion = model_data$ecoregion_ortho[1],
                         scenario = rep(c("baseline","plus_one","minus_one"), each = length(ante_N_dep_values)),
                         pred = c(pred_baseline, pred_increase, pred_decrease),
                         Dep_Nhistoric = rep(ante_N_dep_values, times = 3))
      
        final_pred <- pred %>%
          mutate(scenario = factor(scenario, levels = c("minus_one", "baseline", "plus_one")))
  
  my_col <- RColorBrewer::brewer.pal(3, "Blues")
  
  mean_growth <- df1[1,]
  
  p <- ggplot()+
    geom_hline(data = mean_growth, aes(yintercept = mean_growth), linetype = 2, linewidth = 1)+
    geom_line(data = final_pred, aes(x = Dep_Nhistoric, y = pred, group = scenario, color = scenario), linewidth = 1)+
    theme_classic()+
    scale_color_manual(values = my_col, name = "Short-term change in N dep.", 
                       labels = c(expression(paste("-1 kg N ", ha^-1," ",y^-1)), 
                                  "no change",
                                  expression(paste("+1 kg N ", ha^-1," ",y^-1))))+
    xlab(expression(paste("antecedent N deposition (kg N ", ha^-1," ",y^-1,")")))+
    ylab(expression(paste("tree growth (kg C ", ind^-1," ",y^-1,")")))+
    guides(color = guide_legend(reverse = TRUE))+
    annotate("text", x = 19, y = mean_growth$mean_growth + 0.3, label = "mean observed \ngrowth rate") +
    ggtitle(expression(paste(italic("Acer saccharum"), " in EASTERN TEMPERATE FORESTS")))
    
  
  ggsave(p,filename = "./visualizations/touchpoint.png",
         device = "png")
  
  
  
  
  


