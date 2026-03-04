# McDonnell data figure
# Author: Mary Lofton
# Date: 17OCT25
# Last updated: 22DEC25

# Purpose: supplemental data figure to show the number of trees of each species per plot


# load packages
library(tidyverse)
library(lubridate)
library(sf)
library(viridis)
library(ggthemes)
library(ggpubr)

og_df <- read_csv("./data/McDonnell_etal_InPrep_TreeData_2024_10_11.csv", show_col_types = FALSE) %>%
  dplyr::filter(!common_name %in% c("Douglas-fir","western hemlock")) 
#filter(common_name %in% c("eastern cottonwood"))

df <- read_csv("./data/processed_data.csv")

spp_df <- og_df %>%
  select(species, common_name) %>%
  distinct(.)

# number of trees of each species per plot map

usa_coordinates <- map_data("state")
map_data <- left_join(df, spp_df, by = "common_name") %>%
  select(plot_ID, lat, lon, tree_ID, species) %>%
  distinct(.) %>%
  group_by(plot_ID, lat, lon) %>%
  count(species)

fig1_supp <- ggplot() +
  geom_map(
    data = usa_coordinates, map = usa_coordinates,
    aes(long, lat, map_id = region),
    color = "black", fill = "white")+
  geom_point(
    data = map_data,
    aes(lon, lat, color = n),
    shape = 16, alpha = 1, size = 1
  ) +
  scale_color_viridis(option = "H", trans = "log10")+
  xlab("")+
  ylab("")+
  theme_classic()+
  labs(color = "Number of trees per plot")+
  theme(legend.position = "bottom",
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.line.x = element_blank(),
        axis.line.y = element_blank())+
  facet_wrap(facets = vars(species), nrow = 4, ncol = 2)

fig1_supp


ggsave(plot = fig1_supp, filename = "./visualizations/final_figures/FigureS1.tif",
       device = "tiff", height = 8, width = 6, units = "in",bg = "white")
