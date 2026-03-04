# Model predictions vs observations
# Author: Mary Lofton
# Date: 03MAR25

# Purpose: calculate model predictions using JAGS model output and compare to observations

delta_growth_vs_delta_Ndep <- function(data = "./data/McDonnell_etal_InPrep_TreeData_2024_10_11.csv", 
                        model_output_folder = "./experiments/delta_Ndep",
                        experiment_name = "new_delta_Ndep_only_saveTreeEffect"){
  
  ## READ IN AND WRANGLE MODEL OUTPUT
  # list files
  out <- list.files(model_output_folder,pattern = "mcmc.parquet",
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
    select(-procErr)
  
  ## READ IN AND WRANGLE DATA
  og_df <- read_csv(data, show_col_types = FALSE) %>%
    dplyr::filter(!common_name %in% c("Douglas-fir","western hemlock")) 
  #filter(common_name %in% c("eastern cottonwood"))
  
  focal_df <- og_df %>%
    select(species, common_name, plot_ID, tree_ID, interval_no, Dep_N, Dep_Noxi, Dep_Nred, 
           Dep_S, MAT, MAP, date_m2, date_m1, AG_carbon_pYear, AG_carbon_m1, 
           AG_carbon_m2, subp_BA_GT_m1, live_m2) 
  
  baseline_vars <- focal_df %>%
    select(plot_ID, date_m1, date_m2, Dep_N, Dep_Noxi, Dep_Nred, Dep_S, MAT, MAP) %>%
    distinct(.) %>%
    group_by(plot_ID) %>%
    summarize(Dep_Nbaseline = mean(Dep_N, na.rm = TRUE),
              Dep_Noxibaseline = mean(Dep_Noxi, na.rm = TRUE),
              Dep_Nredbaseline = mean(Dep_Nred, na.rm = TRUE),
              Dep_Sbaseline = mean(Dep_S, na.rm = TRUE),
              MAT_baseline = mean(MAT, na.rm = TRUE),
              MAP_baseline = mean(MAP, na.rm = TRUE)) %>%
    ungroup()
  
  focal_df2 <- left_join(focal_df, baseline_vars, by = "plot_ID") %>%
    group_by(plot_ID) %>%
    mutate(Dep_Ndelta = Dep_N - Dep_Nbaseline,
           Dep_Noxidelta = Dep_Noxi - Dep_Noxibaseline,
           Dep_Nreddelta = Dep_Nred - Dep_Nredbaseline,
           Dep_Sdelta = Dep_S - Dep_Sbaseline,
           MAP_delta_dm = (MAP - MAP_baseline) * 0.01,
           MAT_delta = MAT - MAT_baseline) %>%
    ungroup() 
  
  live_tree_ids <- focal_df2 |> 
    summarise(count = n(),
              sum = sum(live_m2), .by = tree_ID) |> 
    dplyr::filter(count >= sum) |> 
    pull(tree_ID)
  
  focal_df3 <- focal_df2 |> 
    dplyr::filter(tree_ID %in% live_tree_ids) |> 
    group_by(tree_ID) |> 
    tidyr::fill(subp_BA_GT_m1, .direction = "down") |> 
    ungroup() 
  
  df <- focal_df3 %>%
    dplyr::filter(complete.cases(.)) %>%
    group_by(tree_ID) %>%
    mutate(num_intervals = max(interval_no, na.rm = TRUE)) %>%
    ungroup() %>%
    dplyr::filter(num_intervals >= 2)
  
  df1 <- df %>%
    group_by(species, common_name) %>%
    summarize(mean_size = mean(AG_carbon_m1, na.rm = TRUE),
              mean_Dep_Sdelta = mean(Dep_Sdelta, na.rm = TRUE),
              mean_MAT_delta = mean(MAT_delta, na.rm = TRUE),
              mean_MAP_delta_dm = mean(MAP_delta_dm, na.rm = TRUE),
              mean_ba_gt = mean(subp_BA_GT_m1, na.rm = TRUE)) %>%
    mutate(common_name = ifelse(common_name == "yellow-poplar","yellow poplar",common_name)) %>%
    ungroup()
  
  species_list <- unique(df1$common_name)
  
  Ndep_range <- df %>%
    summarize(max_delta_Ndep = max(Dep_Ndelta, na.rm = TRUE),
              min_delta_Ndep = min(Dep_Ndelta, na.rm = TRUE),
              mean_delta_Ndep = mean(Dep_Ndelta, na.rm = TRUE),
              median_delta_Ndep = median(Dep_Ndelta, na.rm = TRUE))
  
  max_delta_Ndep = Ndep_range$max_delta_Ndep
  min_delta_Ndep = Ndep_range$min_delta_Ndep
  
  plot_data <- left_join(df1, final_param_sum, by = join_by(common_name))
  
  my.cols <- colorblind_pal()(8)
  names(my.cols) <- unique(df1$species)
    
  if(experiment_name == "new_delta_Ndep_only_saveTreeEffect"){
    plot_params <- plot_data %>%
      select(-common_name, -species, -mean_Dep_Sdelta, -mean_MAT_delta, -mean_MAP_delta_dm) %>%
      rename(tree_agb_obs = mean_size,
             ba_gt = mean_ba_gt) %>%
      mutate(color = my.cols) 
    
    
  p1 <- ggplot(data = plot_params) +
    purrr::pmap(
      plot_params,
      function(global_tree_effect,p5,tree_agb_obs,p2,ba_gt,p3,color) geom_function(data = plot_params,fun = function(x){
        growth_baseline =  ((global_tree_effect + 0*p5) * tree_agb_obs ^ p2) * exp(-ba_gt*p3)
        growth_delta =  ((global_tree_effect + x*p5) * tree_agb_obs ^ p2) * 1 #exp(-ba_gt*p3) 
        y = growth_delta - growth_baseline
        return(y)
      }, col = color)
    ) +
    xlim(min_delta_Ndep, 0) +
    theme_bw() +
    xlab(expression(paste("N deposition change (kg N ", ha^-1," ",y^-1,")"))) +
    ylab(expression(paste("growth (kg C ", y^-1," ",ind^-1,")")))
  }
  if(experiment_name == "delta_env_saveTreeEffect"){
    
    plot_params <- plot_data %>%
      select(-common_name, -species) %>%
      rename(tree_agb_obs = mean_size,
             ba_gt = mean_ba_gt,
             Dep_Sdelta = mean_Dep_Sdelta,
             MAT_delta = mean_MAT_delta,
             MAP_delta_dm = mean_MAP_delta_dm) %>%
      mutate(color = my.cols) 
    
    p1 <- ggplot(data = plot_params) +
      purrr::pmap(
        plot_params,
        function(global_tree_effect,p5,Dep_Sdelta,p6,MAT_delta,p7,MAP_delta_dm,p8,tree_agb_obs,p2,ba_gt, p3, color) geom_function(data = plot_params,fun = function(x){
          growth_baseline =  ((global_tree_effect + 0*p5 + 0*p6 + 0*p7 + 0*p8) * tree_agb_obs ^ p2) * exp(-ba_gt*p3)
          growth_delta =  ((global_tree_effect + x*p5 + 0*p6 + 0*p7 + 0*p8) * tree_agb_obs ^ p2) * 1 #exp(-ba_gt*p3)
          y = growth_delta - growth_baseline
          return(y)
        }, col = color)
      ) +
      xlim(min_delta_Ndep, 0) +
      theme_bw() +
      xlab(expression(paste("N deposition deviation (kg N ", ha^-1," ",y^-1,")"))) +
      ylab(expression(paste("change in growth (kg C ", y^-1," ",ind^-1,")")))
  }
  p1
  
  plot_params_legend <- plot_data 
  
  leg_p1_plot <- ggplot(data = plot_params_legend, aes(x = p2, y = p5,
                                                  group = species, color = species))+
    geom_line()+
    theme_bw()+
    scale_color_manual(values = my.cols)+
    labs(color = "Species")
  
  # Extract the legend. Returns a gtable
  leg_p1 <- get_legend(leg_p1_plot)
  
  # Convert to a ggplot and print
  leg <- as_ggplot(leg_p1)
  leg
  
  p <- ggarrange(p1,leg,
                  ncol = 2,nrow = 1,
                  widths = c(1,0.5)
  ) 
  
  return(p)
  
  message("all done!")
    
  }
