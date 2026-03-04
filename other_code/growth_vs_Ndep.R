# Model predictions vs observations
# Author: Mary Lofton
# Date: 03MAR25

# Purpose: calculate model predictions using JAGS model output and compare to observations

growth_vs_Ndep <- function(data = "./data/McDonnell_etal_InPrep_TreeData_2024_10_11.csv", 
                        model_output_folder = "./experiments/delta_Ndep"){
  
  ## READ IN AND WRANGLE MODEL OUTPUT
  # list files
  out <- list.files(model_output_folder,pattern = ".parquet",
                    full.names = TRUE)
  
  for(i in 1:length(out)){
    
    model_name = str_split(out[i], pattern = "-")[[1]][2]
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
    select(common_name, global_tree_effect:p5)
  
  ## READ IN AND WRANGLE DATA
  df <- read_csv(data, show_col_types = FALSE) %>%
    filter(!common_name %in% c("Douglas-fir","western hemlock"))
  
  live_tree_ids <- df |> 
    summarise(count = n(),
              sum = sum(live_m2), .by = tree_ID) |> 
    dplyr::filter(count >= sum) |> 
    pull(tree_ID)
  
  focal_data <- df |> 
    dplyr::filter(tree_ID %in% live_tree_ids) |> 
    group_by(common_name, tree_ID) |> 
    tidyr::fill(subp_BA_GT_m1, .direction = "down") |> 
    dplyr::filter(!is.na(subp_BA_GT_m1) & !is.na(Dep_N)) |>
    ungroup() 
  
  baseline_Ndep <- focal_data %>%
    select(common_name, tree_ID, interval_no, Dep_N15, Dep_N, date_m2, date_m1) %>%
    filter(interval_no == 1) %>%
    mutate(dt = as.numeric(date_m2 - date_m1)/365) %>%
    mutate(total_years = 15 + dt) %>%
    mutate(Dep_Nbaseline = (Dep_N15*total_years - Dep_N*dt) / 15) %>%
    select(tree_ID, Dep_Nbaseline)
  
  focal_data2 <- left_join(focal_data, baseline_Ndep, by = "tree_ID") %>%
    mutate(Dep_Ndelta = Dep_N - Dep_Nbaseline) %>%
    filter(!is.na(Dep_Nbaseline) & !is.na(Dep_Ndelta))
  
  df1 <- focal_data2 %>%
    group_by(common_name) %>%
    summarize(min_baseline_Ndep = min(Dep_Nbaseline, na.rm = TRUE),
              max_baseline_Ndep = max(Dep_Nbaseline, na.rm = TRUE),
              mean_delta_Ndep = mean(Dep_Ndelta, na.rm = TRUE),
              mean_size = mean(AG_carbon_m1, na.rm = TRUE)) %>%
    mutate(common_name = ifelse(common_name == "yellow-poplar","yellow poplar",common_name))
  
  growth_v_Ndep <- function(global_tree_effect,
                            x,
                            p4,
                            ndep_delta,
                            p5,
                            tree_agb_obs,
                            p2){
    y =  ((global_tree_effect + x*p4 + ndep_delta*p5) * tree_agb_obs ^ p2) 
  }
  
  species_list <- unique(df1$common_name)
  
  plot_list <- NULL
  
  for(i in 1:length(species_list)){
    
    current_species_data <- df1 %>%
      filter(common_name == species_list[i])
    
    current_species_params <- final_param_sum %>%
      filter(common_name == species_list[i])
    
    min_Ndep <- current_species_data$min_baseline_Ndep
    max_Ndep <- current_species_data$max_baseline_Ndep
    
    base <-
      ggplot() +
      xlim(min_Ndep, max_Ndep) +
      geom_function(fun = growth_v_Ndep, args = list(global_tree_effect = current_species_params$global_tree_effect,
                                                     tree_agb_obs = current_species_data$mean_size,
                                                     ndep_delta = 0,
                                                     p2 = current_species_params$p2,
                                                     p4 = current_species_params$p4,
                                                     p5 = current_species_params$p5),
                    aes(col = "mean baseline N dep.")) +
      geom_function(fun = growth_v_Ndep, args = list(global_tree_effect = current_species_params$global_tree_effect,
                                                     tree_agb_obs = current_species_data$mean_size,
                                                     ndep_delta = current_species_data$mean_delta_Ndep,
                                                     p2 = current_species_params$p2,
                                                     p4 = current_species_params$p4,
                                                     p5 = current_species_params$p5),
                    aes(col = "mean baseline N dep. +/- mean change"))+
      scale_color_manual(values = c("mean baseline N dep." = "black",
                                    "mean baseline N dep. +/- mean change" = "red"))+
      xlab("baseline N deposition")+
      ylab("AG_carbon_pYear")+
      ggtitle(paste0(species_list[i],": mean delta Ndep = ",round(current_species_data$mean_delta_Ndep,2)))+
      theme_bw()+
      theme(legend.position = "bottom")
    base
    
    plot_list[[i]] <- base
  }
  
  final_plot <- ggpubr::ggarrange(plotlist=plot_list,ncol = 2, nrow = 4) +
    bgcolor("white")

  ggsave(plot = final_plot, filename = "./visualizations/delta_Ndep_growth_vs_Ndep.tif",
         device = "tiff", height = 12, width = 10, units = "in")
  
  message("all done!")
    
  }
