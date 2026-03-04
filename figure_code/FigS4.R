# Boxplot of growth rates by species to show outliers
# Author: Mary Lofton
# Date: 04NOV25
# Last updated: 22DEC25

library(tidyverse)
library(lubridate)
library(scales)

spp_df <- read_csv("./data/McDonnell_etal_InPrep_TreeData_2024_10_11.csv", show_col_types = FALSE) %>%
  dplyr::filter(!common_name %in% c("Douglas-fir","western hemlock")) %>%
  select(common_name, species) %>%
  distinct(.)

df <- read_csv("./data/processed_data.csv") %>%
  left_join(., spp_df, by = "common_name")

means <- df %>%
  group_by(species) %>%
  summarize(mean_growth = mean(AG_carbon_pYear, na.rm = TRUE),
            sd_growth = sd(AG_carbon_pYear, na.rm = TRUE))

plot_data <- left_join(df, means, by = "species") %>%
  mutate(outlier = ifelse((AG_carbon_pYear > (mean_growth + 3*sd_growth) | AG_carbon_pYear < (mean_growth - 3*sd_growth)),"yes","no"))

p <- ggplot(data = plot_data, aes(x = species, y = AG_carbon_pYear))+
  geom_boxplot(outlier.shape = 1)+
  theme_classic()+
  xlab("")+
  ylab(expression(paste("observed tree growth (kg C ", y^-1," ",ind^-1,")")))+
  scale_x_discrete(labels = label_wrap(width = 12))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

ggsave(p, filename = "./visualizations/final_figures/FigS4.png", dev = "png",
       height = 4, width = 4, units = "in")
