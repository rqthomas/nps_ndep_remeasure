# Marginal effect of N deposition on growth
# Author: Mary Lofton
# Date: 03APR25
# Last updated: 22DEC25

# Purpose: figure with marginal effect of N deposition (decrease of 1 kg
# per hectare per year) on growth (kg C per year per individual)

library(tidyverse)
library(lubridate)
library(scales)

data = "./data/McDonnell_etal_InPrep_TreeData_2024_10_11.csv"

# list files in each model output folder
out <- list.files("./experiments/ortho_log_t_interaction_01JAN26",pattern = "mcmc.parquet",
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
  mutate(spp_id = ifelse(spp_id == "yellow","yellow-poplar",spp_id)) %>%
  select(model_id, spp_id, nu, .chain, .iteration, .draw) %>%
  rename(common_name = spp_id)

spp_df <- read_csv(data) %>%
  select(common_name, species) %>%
  distinct(.)

plot_dat <- left_join(final1, spp_df, by = "common_name")

p1 <- ggplot(plot_dat, aes(x = .iteration, y = nu, group = as.factor(.chain), color = as.factor(.chain)))+
  geom_line()+
  facet_wrap(facets = vars(species), scales = "free")+
  theme_bw()+
  labs(color = NULL)+
  xlab("iteration")+
  ylab(expression(paste(nu)))
  
  ggsave(p1,filename = "./visualizations/final_figures/FigureS26.png",
         device = "png", bg = "white", height = 7, width = 7.5, units = "in")
  

