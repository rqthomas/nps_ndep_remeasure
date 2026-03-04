# McDonnell data figure
# Author: Mary Lofton
# Date: 17OCT25
# Last updated: 22DEC25

# Purpose: data figure to illustrate:

#'A. plots by ecoregion
#'B. number of remeasurements plot
#'C. number of trees measured per species over time

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

# number of remeasurements by species

spp_df <- og_df %>%
  select(common_name, species) %>%
  distinct(.)

plot_dat <- left_join(df, spp_df, by = "common_name") %>%
  select(species, num_intervals) %>%
  group_by(species, num_intervals) %>%
  count(num_intervals) %>% 
  ungroup() %>%
  group_by(species) %>%
  mutate(perc = n / sum(n) * 100) 

fig1_supp2 <- ggplot(plot_dat, aes(x = num_intervals, y = perc)) + 
  geom_bar(stat = "identity") +
  theme_classic() +
  ylab("Percent of total trees of a species") +
  xlab("Number of measurement intervals") +
  facet_wrap(facets = vars(species))

fig1_supp2

ggsave(plot = fig1_supp2, filename = "./visualizations/final_figures/FigureS2.tif",
       device = "tiff", bg = "white")
