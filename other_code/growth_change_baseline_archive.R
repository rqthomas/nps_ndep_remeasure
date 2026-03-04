# Tree growth given long-term change and historic baseline
# Author: Mary Lofton
# Date: 02JUN25

# Purpose: plot tree growth given different levels of long-term change in N deposition
# and historic baseline values

source("./other_code/get_model_inputs.R")

growth_change_baseline <- function(data = "./data/McDonnell_etal_InPrep_TreeData_2024_10_11.csv",
                                     model_output_folder = "./experiments/short-term_long-term",
                                   temporal_scale = "longterm"){
  
  # get data
  df <- get_model_inputs(data)
  
  focal_data <- df %>%
    select(common_name, species, Dep_Nhistoric,
           Dep_N_LTchange, Dep_Ndelta) %>%
    mutate(spp_id = ifelse(common_name == "yellow-poplar","yellow poplar",common_name))
  
  mean_data <- df %>%
    select(common_name, AG_carbon_m1, subp_BA_GT_m1) %>%
    group_by(common_name) %>%
    summarize_all(mean, na.rm = TRUE) %>%
    rename(spp_id = common_name) %>%
    mutate(spp_id = ifelse(spp_id == "yellow-poplar","yellow poplar",spp_id)) 
  
  N_ranges <- NULL
  
  spp_list <- unique(focal_data$species)
  
  for(i in 1:length(spp_list)){
    
    curr_dat <- focal_data %>%
      filter(species == spp_list[i])
    
    min_Nbaseline <- min(curr_dat$Dep_Nhistoric, na.rm = TRUE)
    max_Nbaseline <- max(curr_dat$Dep_Nhistoric, na.rm = TRUE)
    min_Ndelta <- min(curr_dat$Dep_N_LTchange, na.rm = TRUE)
    max_Ndelta <- max(curr_dat$Dep_N_LTchange, na.rm = TRUE)

    baseline_Ndep_values = seq(min_Nbaseline, max_Nbaseline, by = 0.1)
    delta_Ndep_values = seq(min_Ndelta, max_Ndelta, by = 0.1)
    
    temp <- data.frame(baseline_Ndep = rep(baseline_Ndep_values, times = length(delta_Ndep_values)),
                           delta_Ndep = rep(delta_Ndep_values, each = length(baseline_Ndep_values)),
                           species = spp_list[i],
                       spp_id = curr_dat$common_name[1]) %>%
      mutate(spp_id = ifelse(spp_id == "yellow-poplar","yellow poplar",spp_id)) %>%
      filter(baseline_Ndep + delta_Ndep >= 0)
    
    if(i == 1){
      N_ranges <- temp
    } else {
      N_ranges <- bind_rows(N_ranges, temp)
    }
    
  }
  
  # get model output
  out <- list.files(model_output_folder,pattern = "mcmc.parquet",
                    full.names = TRUE)
  
  for(i in 1:length(out)){
    
    spp_name = str_split(out[i], pattern = "-")[[1]][6]
    model_name = str_split(out[i], pattern = "/")[[1]][3]
    temp <- read_parquet(file = out[i]) %>%
      mutate(spp_id = spp_name,
             model_id = model_name)
    
    if(i == 1){
      final <- temp
    } else {
      final <- bind_rows(final, temp)
    }
    
  }
  
  summary_final <- final %>%
    mutate(spp_id = ifelse(spp_id == "yellow","yellow poplar",spp_id)) %>%
    select(-.draw,-.chain,-.iteration, -model_id) %>%
    group_by(spp_id) %>%
    summarize_all(mean, na.rm = TRUE)
  
  plot_data <- left_join(N_ranges, summary_final, by = "spp_id") %>%
    left_join(mean_data, by = "spp_id") %>%
    distinct() %>%
    group_by(spp_id) %>%
    mutate(pred = ((global_tree_effect + 0*p5 + 0*p6 + 0*p7 + 0*p8 + delta_Ndep*p9 + baseline_Ndep*p10) * AG_carbon_m1 ^ p2 * 1 ) ) %>%
    mutate(pred_perc = (pred / AG_carbon_m1) * 100)
  
  species_list <- unique(plot_data$species)
  
  plots <- NULL
  
  for(i in 1:length(species_list)){
    
    plot_dat <- plot_data %>%
      filter(species == species_list[i]) 
    
    plots[[i]] <- ggplot(data = plot_dat)+
      geom_tile(aes(x = baseline_Ndep, y = delta_Ndep, fill = pred_perc))+
      scale_fill_viridis_c()+
      ggtitle(species_list[i])+
      theme_bw()+
      xlab(expression(paste("historic N deposition (kg N ", ha^-1," ",yr^-1,")")))+
      ylab(expression(paste("N deposition change  (kg N ", ha^-1," ",yr^-1,")")))+
      labs(fill = expression(paste(atop(paste("Percent change\nin aboveground\nbiomass ")),paste("(kg C ",ind^-1," ",yr^-1,")"))))
    
  }
  
  money_plot <- ggarrange(plots[[1]], plots[[2]], plots[[3]], plots[[4]], plots[[5]], plots[[6]], plots[[7]], plots[[8]], 
                          nrow = 3, ncol = 3,
                          labels = c("A","B","C","D","E","F","G","H")
  )
  return(money_plot)
  
  
  
}

ggsave(plot = money_plot, filename = "./visualizations/deltaGrowth_deltaN_baselineN_tile.tif",
       device = "tiff", height = 10, width = 18, units = "in",bg = "white")
