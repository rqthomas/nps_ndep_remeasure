# Model predictions vs observations
# Author: Mary Lofton
# Date: 03MAR25

# Purpose: calculate model predictions using JAGS model output and compare to observations

source("./other_code/get_model_inputs.R")

growth_vs_size_longtermDepN <- function(data = "./data/McDonnell_etal_InPrep_TreeData_2024_10_11.csv", 
                        model_output_folder = "./experiments/short-term_long-term"){
  
  ## READ IN AND WRANGLE MODEL OUTPUT
  # list files
  out <- list.files(model_output_folder,pattern = "mcmc.parquet",
                    full.names = TRUE)
  
  for(i in 1:length(out)){
    
    model_name = str_split(out[i], pattern = "-")[[1]][6]
    temp <- read_parquet(file = out[i]) %>%
      mutate(model_id = model_name)
    
    if(i == 1){
      final <- temp
    } else {
      final <- bind_rows(final, temp)
    }
    
  }
  
  final <- final %>%
    mutate(model_id = ifelse(model_id == "yellow","yellow poplar",model_id))
  
  param_sum <- final %>%
    select(-c(.chain,.iteration,.draw)) %>%
    pivot_longer(-model_id,names_to = "param", values_to = "param_value") %>%
    group_by(model_id, param) %>%
    summarize(mean_param_value = mean(param_value, na.rm = TRUE)) %>%
    pivot_wider(names_from = "param", values_from = "mean_param_value") %>%
    ungroup() %>%
    rename(common_name = model_id) 
  
  final_param_sum <- param_sum %>% 
    select(common_name, global_tree_effect:p9)
  
  ## READ IN AND WRANGLE DATA
  df <- get_model_inputs(data)
  
  df1 <- df %>%
    group_by(common_name, species) %>%
    summarize(min_size = min(AG_carbon_m1, na.rm = TRUE),
              max_size = max(AG_carbon_m1, na.rm = TRUE),
              mean_ba_gt = mean(subp_BA_GT_m1, na.rm = TRUE),
              ndep_historic = mean(Dep_Nhistoric, na.rm = TRUE)) %>%
    mutate(common_name = ifelse(common_name == "yellow-poplar","yellow poplar",common_name))
  
  growth_v_size <- function(global_tree_effect,
                            ndep_longterm_delta,
                            ndep_historic,
                            p9,
                            x,
                            p2,
                            p3,
                            p10,
                            ba_gt){
    y =  ((global_tree_effect + ndep_longterm_delta*p9 + ndep_historic*p10) * x ^ p2) * 1 #exp(-ba_gt*p3) 
  }
  
  species_list <- sort(unique(df1$species))
  
  plot_list <- NULL
  
  for(i in 1:length(species_list)){
    
    current_species_data <- df1 %>%
      filter(species == species_list[i])
    
    current_species_params <- final_param_sum %>%
      filter(common_name == current_species_data$common_name[1])
    
    min_size <- current_species_data$min_size
    max_size <- current_species_data$max_size
    
    base <-
      ggplot() +
      xlim(min_size, max_size) +
      geom_function(fun = growth_v_size, args = list(global_tree_effect = current_species_params$global_tree_effect,
                                                     ba_gt = current_species_data$mean_ba_gt,
                                                     ndep_longterm_delta = 0,
                                                     ndep_historic = current_species_data$ndep_historic,
                                                     p2 = current_species_params$p2,
                                                     p3 = current_species_params$p3,
                                                     p9 = current_species_params$p9,
                                                     p10 = current_species_params$p10),
                    aes(col = "baseline N dep.")) +
      geom_function(fun = growth_v_size, args = list(global_tree_effect = current_species_params$global_tree_effect,
                                                     ba_gt = current_species_data$mean_ba_gt,
                                                     ndep_longterm_delta = -1,
                                                     ndep_historic = current_species_data$ndep_historic,
                                                     p2 = current_species_params$p2,
                                                     p3 = current_species_params$p3,
                                                     p9 = current_species_params$p9,
                                                     p10 = current_species_params$p10),
                    aes(col = "marginal N dep. decrease"))+
      geom_function(fun = growth_v_size, args = list(global_tree_effect = current_species_params$global_tree_effect,
                                                     ba_gt = current_species_data$mean_ba_gt,
                                                     ndep_longterm_delta = 1,
                                                     ndep_historic = current_species_data$ndep_historic,
                                                     p2 = current_species_params$p2,
                                                     p3 = current_species_params$p3,
                                                     p9 = current_species_params$p9,
                                                     p10 = current_species_params$p10),
                    aes(col = "marginal N dep. increase"))+
      scale_color_manual(values = c("baseline N dep." = "black",
                                    "marginal N dep. decrease" = "blue",
                                    "marginal N dep. increase" = "red"))+
      xlab(expression(paste("tree size (kg C ", ind^-1,")")))+
      ylab(expression(paste("growth (kg C ", y^-1," ",ind^-1,")")))+
      ggtitle(current_species_data$species)+
      theme_bw()+
      theme(legend.position = "bottom")+
      labs(color = "")
    
    if(i == 1){
      # Extract the legend. Returns a gtable
      legend <- get_legend(base)
      
      # Convert to a ggplot and print
      leg <- as_ggplot(legend)
    }
    
    base <- base + theme(legend.position = "none")
    
    plot_list[[i]] <- base
  }
  
  final_plot <- ggpubr::ggarrange(ggarrange(plotlist=plot_list,ncol = 2, nrow = 4),
                                  leg,
                                  nrow = 2,
                                  heights = c(1, 0.1)) 
  
  return(final_plot)
  
  message("all done!")
    
  }
