# Growth change given short-term vs. long-term N deposition decreases
# Author: Mary Lofton
# Date: 02JUN25

# Purpose: visualize changes in tree growth across parameter values at both short-
# and long-term

source("./other_code/get_model_inputs.R")

shortterm_vs_longterm <- function(data = "./data/McDonnell_etal_InPrep_TreeData_2024_10_11.csv",
                                    model_output_folder = "./experiments/short-term_long-term",
                                  zoom = FALSE){
  
  # get data
  df <- get_model_inputs(data)
  
  mean_data <- df %>%
    select(common_name, species, AG_carbon_m1, subp_BA_GT_m1, Dep_Nhistoric) %>%
    group_by(common_name, species) %>%
    summarize_all(mean, na.rm = TRUE) %>%
    mutate(spp_id = ifelse(common_name == "yellow-poplar","yellow poplar",common_name))
  
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
  
  final <- final %>%
    mutate(spp_id = ifelse(spp_id == "yellow","yellow poplar",spp_id))
  
  plot_data <- left_join(final, mean_data, by = "spp_id") %>%
    group_by(spp_id) %>%
    sample_n(2000) %>%
    mutate(st_decrease = (((global_tree_effect + -1*p5 + 0*p6 + 0*p7 + 0*p8 + 0*p9 + Dep_Nhistoric*p10) * AG_carbon_m1 ^ p2 * 1 ) - ((global_tree_effect + 0*p5 + 0*p6 + 0*p7 + 0*p8 + 0*p9 + Dep_Nhistoric*p10) * AG_carbon_m1 ^ p2 * 1 )),
           lt_decrease = (((global_tree_effect + 0*p5 + 0*p6 + 0*p7 + 0*p8 + -1*p9 + Dep_Nhistoric*p10) * AG_carbon_m1 ^ p2 * 1 ) - ((global_tree_effect + 0*p5 + 0*p6 + 0*p7 + 0*p8 + 0*p9 + Dep_Nhistoric*p10) * AG_carbon_m1 ^ p2 * 1 )))
  
  p <- ggplot(data = plot_data)+
    geom_hline(yintercept = 0)+
    geom_vline(xintercept = 0)+
    geom_point(aes(x = lt_decrease, y = st_decrease, col = species), size = 0.1, alpha = 0.1)+
    stat_ellipse(aes(x = lt_decrease, y = st_decrease, col = species), level = 0.95, alpha = 0.5)+
    theme_bw()+
    scale_color_colorblind()+
    ylab("Growth change given 1 kg ha^-1 yr^-1 N 'short-term' decrease")+
    xlab("Growth change given 1 kg ha^-1 yr^-1 N 'long-term' decrease")+
    labs(color = "Species")
  
  if(zoom == TRUE){
    plot_data <- plot_data %>%
      filter(!common_name == "eastern cottonwood")
    
    my.cols <- colorblind_pal()
    my.real.cols <- my.cols(8)[c(1:5,7:8)]
    
    p <- ggplot(data = plot_data)+
      geom_hline(yintercept = 0)+
      geom_vline(xintercept = 0)+
      geom_point(aes(x = lt_decrease, y = st_decrease, col = species), size = 0.1, alpha = 0.1)+
      stat_ellipse(aes(x = lt_decrease, y = st_decrease, col = species), level = 0.95, alpha = 0.5)+
      theme_bw()+
      scale_color_manual(values = my.real.cols)+
      ylab("Growth change given 1 kg ha^-1 yr^-1 N 'short-term' decrease")+
      xlab("Growth change given 1 kg ha^-1 yr^-1 N 'long-term' decrease")+
      labs(color = "Species") +
      ylim(c(-1, 1))+
      xlim(c(-1, 1))
    
  }
  
  return(p)
    
  }
  

  