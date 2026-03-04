# Model predictions vs observations
# Author: Mary Lofton
# Date: 03MAR25
# Last updated: 22DEC25

# Purpose: calculate model predictions using JAGS model output and compare to observations

#source("./other_code/get_model_inputs.R")

# load packages
library(tidyverse)
library(lubridate)
library(arrow)
library(ggpubr)
library(ggthemes)
library(ggpmisc)
library(ggh4x)

data = "./data/McDonnell_etal_InPrep_TreeData_2024_10_11.csv" 
model_output_folder = "./experiments/ortho_log_t_interaction_01JAN26"

  ## READ IN AND WRANGLE MODEL OUTPUT
  # list files
  out <- list.files(model_output_folder,pattern = "treeEffect.parquet",
                    full.names = TRUE)
  
  if(!length(grep("short-term_long-term",model_output_folder)) == 0){
    for(i in 1:length(out)){
      
      spp_name = str_split(out[i], pattern = "-")[[1]][6]
      model_name = str_split(out[i], pattern = "/")[[1]][3]
      iters <- sample(c(1:2000),100)
      
      temp <- read_parquet(file = out[i], as_data_frame = FALSE,col_select = c("n","tree_effect",".iteration")) |>
        filter(.iteration %in% iters) |>
        select(-.iteration) |>
        collect() |>
        mutate(model_id = model_name,
               spp_id = spp_name)
      
      if(i == 1){
        final <- temp
      } else {
        final <- bind_rows(final, temp)
      }
      
    }
  } else {
  
  for(i in 1:length(out)){
    
    model_name = str_split(out[i], pattern = "-")[[1]][2]
    iters <- sample(c(1:2000),100)
    
    message(model_name)
    
    temp <- read_parquet(file = out[i], as_data_frame = FALSE,col_select = c("n","tree_effect",".iteration")) |>
      filter(.iteration %in% iters) |>
      select(-.iteration) |>
      collect() |>
      mutate(spp_id = model_name)
    
    if(i == 1){
      final <- temp
    } else {
      final <- bind_rows(final, temp)
    }
    
  }
  }
  
  final <- final %>%
    mutate(spp_id = ifelse(spp_id == "yellow","yellow poplar",spp_id))
  
  tree_effects <- final %>%
    group_by(spp_id, n) %>%
    summarize(mean_tree_effect = mean(tree_effect, na.rm = TRUE)) %>%
    rename(tree_index = n,
           common_name = spp_id)
  
  # list files
  out1 <- list.files(model_output_folder,pattern = "mcmc.parquet",
                    full.names = TRUE)
  
  if(!length(grep("short-term_long-term",model_output_folder)) == 0){
    
    for(i in 1:length(out1)){
    
    model_name = str_split(out1[i], pattern = "-")[[1]][6]
    
    temp <- read_parquet(file = out1[i]) %>%
      mutate(model_id = model_name)
    
    if(i == 1){
      final1 <- temp
    } else {
      final1 <- bind_rows(final1, temp)
    }
    }
    
  } else {
  
  for(i in 1:length(out1)){
    
    model_name = str_split(out1[i], pattern = "-")[[1]][2]

    temp <- read_parquet(file = out1[i]) %>%
      mutate(model_id = model_name)
    
    if(i == 1){
      final1 <- temp
    } else {
      final1 <- bind_rows(final1, temp)
    }
    
  }
  }
  
  final1 <- final1 %>%
    mutate(model_id = ifelse(model_id == "yellow","yellow poplar",model_id))
  
  param_sum <- final1 %>%
    select(-c(.chain,.iteration,.draw)) %>%
    pivot_longer(-model_id,names_to = "param", values_to = "param_value") %>%
    group_by(model_id, param) %>%
    summarize(mean_param_value = mean(param_value, na.rm = TRUE)) %>%
    pivot_wider(names_from = "param", values_from = "mean_param_value") %>%
    ungroup() %>%
    rename(common_name = model_id) 
  
  final_param_sum <- param_sum 
  
  
  ## READ IN AND WRANGLE DATA
  df <- read_csv("./data/processed_data.csv") %>%
    filter(!AG_carbon_pYear < -2000)
  
  trees_index <- df |> 
    group_by(common_name) |> 
    distinct(tree_ID) |> 
    mutate(tree_index = 1:n())
  
  df1 <- df %>%
    left_join(trees_index, by = join_by(common_name, tree_ID)) |>
    mutate(common_name = ifelse(common_name == "yellow-poplar","yellow poplar",common_name)) |>
    left_join(tree_effects, by = join_by(common_name, tree_index)) |>
    left_join(final_param_sum, by = join_by(common_name)) 
  
  df2 <- df1 %>%
    mutate(dt = as.numeric(date_m2 - date_m1)/365,
           log_AG_carbon_pYear = (log(AG_carbon_m2) - log(AG_carbon_m1)) / dt,
           log_AG_carbon_m1 = log(AG_carbon_m1),
           log_subp_BA_GT_m1 = log(subp_BA_GT_m1 + 1),
           pred_log = ((mean_tree_effect + Dep_Ndiff_ortho*p5 + Dep_Sdiff*p6 + MAT_delta*p7 + MAP_delta_dm*p8 + Dep_Nhistoric_ortho*p9 + Dep_Shistoric*p10 + Dep_Nhistoric_ortho*Dep_Ndiff_ortho*p11 + Dep_Ndiff_ortho*Dep_Ndiff_ortho*p12) + p2*log_AG_carbon_m1) - p3*log_subp_BA_GT_m1,
           m2_not_log = exp(pred_log + log_AG_carbon_m1),
           pred = m2_not_log - AG_carbon_m1)
  
  rsq <- function(pred, obs){
    1 - (sum((obs - pred)^2, na.rm = TRUE) / sum((obs - mean(obs, na.rm = TRUE))^2, na.rm = TRUE))
  }
  
  mod_assess0 <- df2 %>%
    group_by(common_name) %>%
    summarize(r2 = rsq(pred = pred, obs = AG_carbon_pYear),
              rmse = sqrt(mean((AG_carbon_pYear - pred)^2, na.rm = TRUE)),
              mae = mean(abs(pred - AG_carbon_pYear), na.rm = TRUE))
  
  spp_df <- read_csv(data) %>%
    dplyr::filter(!common_name %in% c("Douglas-fir","western hemlock")) %>%
    select(species, common_name) %>%
    distinct(.) %>%
    mutate(common_name = ifelse(common_name == "yellow-poplar","yellow poplar",common_name)) 
  
  mod_assess <- left_join(mod_assess0, spp_df, by = "common_name")
  df3 <- left_join(df2, spp_df, by = "common_name")
  
  write.csv(mod_assess, "./experiments/ortho_log_t_interaction_adj_priors/model_performance.csv", row.names = FALSE)
    
  r2 <- ggplot(mod_assess, aes(x=reorder(species, r2), y=r2, color=as.factor(species))) + 
    geom_point() +
    geom_segment(aes(x=species,xend=species,y=0,yend=r2)) +
    ylab(expression(paste(R^2))) +
    xlab("") +
    coord_flip() +
    scale_color_colorblind()+
    theme_bw()+
    theme(legend.position = "none")
  
  rmse <- ggplot(mod_assess, aes(x=reorder(species, -rmse), y=rmse, color=as.factor(species))) + 
    geom_point() +
    geom_segment(aes(x=species,xend=species,y=0,yend=rmse)) +
    ylab(expression(paste("RMSE (kg C ", y^-1," ",ind^-1,")"))) +
    xlab("") +
    coord_flip() +
    scale_color_colorblind()+
    theme_bw()+
    theme(legend.position = "none")
  
  mae <- ggplot(mod_assess, aes(x=reorder(species, -mae), y=mae, color=as.factor(species))) + 
    geom_point() +
    geom_segment(aes(x=species,xend=species,y=0,yend=mae)) +
    ylab(expression(paste("MAE log(kg C ", y^-1," ",ind^-1,")"))) +
    xlab("") +
    coord_flip() +
    scale_color_colorblind()+
    theme_bw()+
    theme(legend.position = "none")
  
  p1 <- ggarrange(r2, rmse, ncol = 2)
  
  grouped_models <- df3 %>%
    select(species, AG_carbon_pYear, pred) %>%
    group_by(species) %>%
    nest() %>%
    mutate(model = map(data, ~lm(pred ~ AG_carbon_pYear, data = .)))
  
  model_results <- grouped_models %>%
    mutate(
      # Extracts coefficients, SE, p-values, and CIs
      tidy_results = map(model, broom::tidy, conf.int = TRUE)
    ) %>%
    unnest(tidy_results) %>%
    select(species, term, estimate, conf.low, conf.high, p.value) # Select relevant columns
  
  # take-home: confidence intervals are tiny and not going to be visible on figure
  
  # remove outliers following request from Chris
  # Calculate IQR limits for observed
  Q1 <- quantile(df3$AG_carbon_pYear, 0.25)
  Q3 <- quantile(df3$AG_carbon_pYear, 0.75)
  IQR <- Q3 - Q1
  lower_limit <- Q1 - 1.5 * IQR
  upper_limit <- Q3 + 1.5 * IQR
  
  df4 <- df3 %>%
    filter(AG_carbon_pYear > lower_limit & AG_carbon_pYear < upper_limit)
  
  custom_scales <- list(
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    scale_y_continuous(limits = c(-60, 40)), # Custom limits for Potr
    NULL
  )
  
  p2 <- ggplot(data = df4, aes(x = AG_carbon_pYear, y = pred))+
    geom_point(alpha = 0.2)+
    geom_abline(slope = 1)+
    theme_bw()+
    theme(legend.position = "none")+
    facet_wrap(facets = vars(species), scales = "free")+
    geom_smooth(method = "lm", formula = y ~ x, se = TRUE)+
    stat_poly_eq(use_label(c("eq")), formula = y~x,label.y = 0.95)+
    stat_poly_eq(use_label(c("R2")), formula = y~x, label.y = 0.85)+
    xlab(expression(paste("observed tree growth (kg C ", y^-1," ",ind^-1,")")))+
    ylab(expression(paste("predicted tree growth (kg C ", y^-1," ",ind^-1,")")))+
    facetted_pos_scales(y = custom_scales)

  p3 <- ggplot(data = df3)+
    geom_point(aes(x = log_AG_carbon_pYear, y = pred_log))+
    geom_abline(slope = 1)+
    theme_bw()+
    theme(legend.position = "none")+
    #ggtitle("Displaying repeated measures")+
    facet_wrap(facets = vars(species), scales = "free")+
    xlab(expression(paste("observed tree growth log(kg C ", y^-1," ",ind^-1,")")))+
    ylab(expression(paste("predicted tree growth log(kg C ", y^-1," ",ind^-1,")")))
  
  ggsave(plot = p1, filename = "./visualizations/final_figures/FigureS24.tif",
         device = "tiff", height = 5, width = 6.6, units = "in",bg = "white")
  ggsave(plot = p2, filename = "./visualizations/final_figures/FigureS25.tif",
         device = "tiff", height = 8, width = 8, units = "in",bg = "white")
  