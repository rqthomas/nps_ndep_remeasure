# McDonnell data figure
# Author: Mary Lofton
# Date: 20OCT25
# Last updated: 22DEC25

# Purpose: data figure to illustrate:

#'A. map of mean N deposition during measurement period (including 15 yr prior to 1st measure)
#'B. timeseries of short-term differences in N deposition by plot and ecoregion

# load packages
library(tidyverse)
library(lubridate)
library(sf)
library(viridis)
library(ggthemes)
library(ggpubr)
library(rnaturalearth)

og_df <- read_csv("./data/McDonnell_etal_InPrep_TreeData_2024_10_11.csv", show_col_types = FALSE) %>%
  dplyr::filter(!common_name %in% c("Douglas-fir","western hemlock")) 
#filter(common_name %in% c("eastern cottonwood"))

df <- read_csv("./data/processed_data.csv")

spp_df <- og_df %>%
  select(species, common_name) %>%
  distinct(.)

# timeseries of short-term differences in N deposition by plot and ecoregion

dat_plot <- df %>%
  select(plot_ID, date_m1, Dep_Ndiff, ecoregion_ortho) %>%
  distinct(.)

fig2_supp <- ggplot(dat_plot, aes(x = date_m1, y = Dep_Ndiff, group = plot_ID)) + 
  geom_line(color = "gray") +
  geom_point(color = "black") +
  theme_classic() +
  ylab(expression(paste("short-term N deposition change (kg N ", ha^-1," ",y^-1,")"))) +
  xlab("measurement period start date") + 
  geom_hline(yintercept = 0, color = "blue")+
  facet_wrap(facets = vars(ecoregion_ortho), labeller = labeller(ecoregion_ortho = label_wrap_gen(width = 30)),
             nrow = 4, scales = "free_y")

fig2_supp


ggsave(plot = fig2_supp, filename = "./visualizations/final_figures/FigureS5.tif",
       device = "tiff", height = 7.5, width = 5.5, units = "in", bg = "white")
