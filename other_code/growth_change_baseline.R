# Tree growth given long-term change and historic baseline
# Author: Mary Lofton
# Date: 02OCT25

# Purpose: plot tree growth given different levels of long-term change in N deposition
# and historic baseline values

growth_change_baseline <- function(data = "./data/McDonnell_etal_InPrep_TreeData_2024_10_11.csv",
                                     model_output_folder = "./experiments/space_vs_time_ortho"){
  
  # get data
  df <- read_csv("./data/processed_data.csv")
  
  # get species
  spp_df <- read_csv(data) %>%
    dplyr::filter(!common_name %in% c("Douglas-fir","western hemlock")) %>%
    select(species, common_name) %>%
    distinct(.)
  
  focal_data <- left_join(df, spp_df, by = "common_name") %>%
    select(common_name, species, Dep_Nhistoric,
           Dep_Ndiff) %>%
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
    
    min_Nhistoric <- min(curr_dat$Dep_Nhistoric, na.rm = TRUE)
    max_Nhistoric <- max(curr_dat$Dep_Nhistoric, na.rm = TRUE)
    min_Ndiff <- min(curr_dat$Dep_Ndiff, na.rm = TRUE)
    max_Ndiff <- max(curr_dat$Dep_Ndiff, na.rm = TRUE)

    ante_Ndep_values = seq(min_Nhistoric, max_Nhistoric, by = 0.1)
    rec_Ndep_values = seq(min_Ndiff, max_Ndiff, by = 0.1)
    
    temp <- data.frame(ante_Ndep = rep(ante_Ndep_values, times = length(rec_Ndep_values)),
                           rec_Ndep = rep(rec_Ndep_values, each = length(ante_Ndep_values)),
                           species = spp_list[i],
                       spp_id = curr_dat$common_name[1]) %>%
      mutate(spp_id = ifelse(spp_id == "yellow-poplar","yellow poplar",spp_id)) %>%
      filter(ante_Ndep + rec_Ndep >= 0)
    
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
    
    spp_name = str_split(out[i], pattern = "-")[[1]][2]
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
