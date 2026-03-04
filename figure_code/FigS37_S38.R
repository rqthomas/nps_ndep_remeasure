# McDonnell data figure
# Author: Mary Lofton
# Date: 18FEB26
# Last updated: 18FEB26

# Purpose: data figure to illustrate:

#'A. map of residual antecedent N deposition averaged across study period
#'B. map of residual short-term change in N deposition averaged across study period

# load packages
library(tidyverse)
library(lubridate)
library(sf)
library(viridis)
library(ggthemes)
library(ggpubr)
library(rnaturalearth)
library(grid)

og_df <- read_csv("./data/McDonnell_etal_InPrep_TreeData_2024_10_11.csv", show_col_types = FALSE) %>%
  dplyr::filter(!common_name %in% c("Douglas-fir","western hemlock")) 
#filter(common_name %in% c("eastern cottonwood"))

df <- read_csv("./data/processed_data.csv")

spp_df <- og_df %>%
  select(species, common_name) %>%
  distinct(.)

usa_coordinates <- map_data("state")

map_data_ante <- left_join(df, spp_df, by = "common_name") %>%
  select(plot_ID, lat, lon, Dep_Nhistoric_SO, date_m1, date_m2) %>%
  distinct(.) %>%
  group_by(plot_ID, lat, lon) %>%
  summarize(Dep_Nhistoric_SO_mean = mean(Dep_Nhistoric_SO, na.rm = TRUE)) %>%
  select(plot_ID, lat, lon, Dep_Nhistoric_SO_mean) %>%
  ungroup()

map_data_st <- left_join(df, spp_df, by = "common_name") %>%
  select(plot_ID, lat, lon, Dep_Ndiff_SO, date_m1, date_m2) %>%
  distinct(.) %>%
  group_by(plot_ID, lat, lon) %>%
  summarize(Dep_Ndiff_SO_mean = mean(Dep_Ndiff_SO, na.rm = TRUE)) %>%
  select(plot_ID, lat, lon, Dep_Ndiff_SO_mean) %>%
  ungroup()

fig1 <- ggplot() +
  geom_map(
    data = usa_coordinates, map = usa_coordinates,
    aes(long, lat, map_id = region),
    color = "black", fill = "white")+
  geom_point(
    data = map_data_ante,
    aes(lon, lat, color = Dep_Nhistoric_SO_mean),
    shape = 16, alpha = 1, size = 1
  ) +
  scale_color_viridis(option = "H", limits = c(-10,25))+
  xlab("")+
  ylab("")+
  theme_classic()+
  labs(color = expression(paste("mean residual antecedent N deposition (kg N ", ha^-1," ",y^-1,")")))+
  guides(color = guide_colorbar(title.position = "bottom", title.hjust = 0.0))+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        legend.background = element_rect(fill = "transparent"))+
  theme(legend.position = "inside", legend.position.inside =  c(0.3, 0.1),
        legend.direction = "horizontal", plot.margin = unit(c(0.1,0,0.1,0), "pt"))+
  ylim(c(min(usa_coordinates$lat),max(usa_coordinates$lat)))

fig1

fig2 <- ggplot() +
  geom_map(
    data = usa_coordinates, map = usa_coordinates,
    aes(long, lat, map_id = region),
    color = "black", fill = "white")+
  geom_point(
    data = map_data_st,
    aes(lon, lat, color = Dep_Ndiff_SO_mean),
    shape = 16, alpha = 1, size = 1
  ) +
  scale_color_viridis(option = "H", limits = c(-10,25))+
  xlab("")+
  ylab("")+
  theme_classic()+
  labs(color = expression(paste("mean residual short-term change in N deposition (kg N ", ha^-1," ",y^-1,")")))+
  guides(color = guide_colorbar(title.position = "bottom", title.hjust = 0.0))+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        legend.background = element_rect(fill = "transparent"))+
  theme(legend.position = "inside", legend.position.inside =  c(0.35, 0.1),
        legend.direction = "horizontal", plot.margin = unit(c(0.1,0,0.1,0), "pt"))+
  ylim(c(min(usa_coordinates$lat),max(usa_coordinates$lat)))

fig2


# Assemble figure

fig <- ggarrange(fig1, fig2,
                nrow = 2, ncol = 1,
                labels = c("(a)","(b)")
)
fig
ggsave(plot = fig, filename = "./visualizations/final_figures/FigureSX.tif",
       device = "tiff", height = 8, width = 7.5, units = "in", bg = "white")

map_data_ante_hist0 <- left_join(df, spp_df, by = "common_name") %>%
  select(plot_ID, lat, lon, Dep_Nhistoric_SO, date_m1, date_m2) %>%
  filter(year(date_m1) %in% c(1997:2000))

range(map_data_ante_hist0$date_m2, na.rm = TRUE)

map_data_ante_hist <- map_data_ante_hist0 %>%
  distinct(.) %>%
  group_by(plot_ID, lat, lon) %>%
  summarize(Dep_Nhistoric_SO_mean = mean(Dep_Nhistoric_SO, na.rm = TRUE)) %>%
  select(plot_ID, lat, lon, Dep_Nhistoric_SO_mean) %>%
  ungroup()

map_data_st_hist0 <- left_join(df, spp_df, by = "common_name") %>%
  select(plot_ID, lat, lon, Dep_Ndiff_SO, date_m1, date_m2) %>%
  filter(year(date_m1) %in% c(1997:2000)) 

range(map_data_st_hist0$date_m2, na.rm = TRUE)

map_data_st_hist <- map_data_st_hist0 %>%
  distinct(.) %>%
  group_by(plot_ID, lat, lon) %>%
  summarize(Dep_Ndiff_SO_mean = mean(Dep_Ndiff_SO, na.rm = TRUE)) %>%
  select(plot_ID, lat, lon, Dep_Ndiff_SO_mean) %>%
  ungroup()

fig3 <- ggplot() +
  geom_map(
    data = usa_coordinates, map = usa_coordinates,
    aes(long, lat, map_id = region),
    color = "black", fill = "white")+
  geom_point(
    data = map_data_ante_hist,
    aes(lon, lat, color = Dep_Nhistoric_SO_mean),
    shape = 16, alpha = 1, size = 1
  ) +
  scale_color_viridis(option = "H", limits = c(-10,25))+
  xlab("")+
  ylab("")+
  theme_classic()+
  labs(color = expression(paste("residual antecedent N deposition (2001-2013; kg N ", ha^-1," ",y^-1,")")))+
  guides(color = guide_colorbar(title.position = "bottom", title.hjust = 0.0))+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        legend.background = element_rect(fill = "transparent"))+
  theme(legend.position = "inside", legend.position.inside =  c(0.3, 0.1),
        legend.direction = "horizontal", plot.margin = unit(c(0.1,0,0.1,0), "pt"))+
  ylim(c(min(usa_coordinates$lat),max(usa_coordinates$lat)))

fig3

fig4 <- ggplot() +
  geom_map(
    data = usa_coordinates, map = usa_coordinates,
    aes(long, lat, map_id = region),
    color = "black", fill = "white")+
  geom_point(
    data = map_data_st_hist,
    aes(lon, lat, color = Dep_Ndiff_SO_mean),
    shape = 16, alpha = 1, size = 1
  ) +
  scale_color_viridis(option = "H", limits = c(-10,25))+
  xlab("")+
  ylab("")+
  theme_classic()+
  labs(color = expression(paste("residual short-term change in N deposition (2001-2013; kg N ", ha^-1," ",y^-1,")")))+
  guides(color = guide_colorbar(title.position = "bottom", title.hjust = 0.0))+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        legend.background = element_rect(fill = "transparent"))+
  theme(legend.position = "inside", legend.position.inside =  c(0.35, 0.1),
        legend.direction = "horizontal", plot.margin = unit(c(0.1,0,0.1,0), "pt"))+
  ylim(c(min(usa_coordinates$lat),max(usa_coordinates$lat)))

fig4

map_data_ante_rec <- left_join(df, spp_df, by = "common_name") %>%
  select(plot_ID, lat, lon, Dep_Nhistoric_SO, date_m1, date_m2) %>%
  filter(year(date_m2) %in% c(2020:2022)) %>%
  distinct(.) %>%
  group_by(plot_ID, lat, lon) %>%
  summarize(Dep_Nhistoric_SO_mean = mean(Dep_Nhistoric_SO, na.rm = TRUE)) %>%
  select(plot_ID, lat, lon, Dep_Nhistoric_SO_mean) %>%
  ungroup()

map_data_st_rec <- left_join(df, spp_df, by = "common_name") %>%
  select(plot_ID, lat, lon, Dep_Ndiff_SO, date_m1, date_m2) %>%
  filter(year(date_m2) %in% c(2020:2022)) %>%
  distinct(.) %>%
  group_by(plot_ID, lat, lon) %>%
  summarize(Dep_Ndiff_SO_mean = mean(Dep_Ndiff_SO, na.rm = TRUE)) %>%
  select(plot_ID, lat, lon, Dep_Ndiff_SO_mean) %>%
  ungroup()

fig5 <- ggplot() +
  geom_map(
    data = usa_coordinates, map = usa_coordinates,
    aes(long, lat, map_id = region),
    color = "black", fill = "white")+
  geom_point(
    data = map_data_ante_rec,
    aes(lon, lat, color = Dep_Nhistoric_SO_mean),
    shape = 16, alpha = 1, size = 1
  ) +
  scale_color_viridis(option = "H", limits = c(-10,25))+
  xlab("")+
  ylab("")+
  theme_classic()+
  labs(color = expression(paste("residual antecedent N deposition (2020-2022; kg N ", ha^-1," ",y^-1,")")))+
  guides(color = guide_colorbar(title.position = "bottom", title.hjust = 0.0))+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        legend.background = element_rect(fill = "transparent"))+
  theme(legend.position = "inside", legend.position.inside =  c(0.3, 0.1),
        legend.direction = "horizontal", plot.margin = unit(c(0.1,0,0.1,0), "pt"))+
  ylim(c(min(usa_coordinates$lat),max(usa_coordinates$lat)))

fig5

fig6 <- ggplot() +
  geom_map(
    data = usa_coordinates, map = usa_coordinates,
    aes(long, lat, map_id = region),
    color = "black", fill = "white")+
  geom_point(
    data = map_data_st_rec,
    aes(lon, lat, color = Dep_Ndiff_SO_mean),
    shape = 16, alpha = 1, size = 1
  ) +
  scale_color_viridis(option = "H", limits = c(-10,25))+
  xlab("")+
  ylab("")+
  theme_classic()+
  labs(color = expression(paste("residual short-term change in N deposition (2020-2022; kg N ", ha^-1," ",y^-1,")")))+
  guides(color = guide_colorbar(title.position = "bottom", title.hjust = 0.0))+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        legend.background = element_rect(fill = "transparent"))+
  theme(legend.position = "inside", legend.position.inside =  c(0.35, 0.1),
        legend.direction = "horizontal", plot.margin = unit(c(0.1,0,0.1,0), "pt"))+
  ylim(c(min(usa_coordinates$lat),max(usa_coordinates$lat)))

fig6


# Assemble figure

fig <- ggarrange(fig3, fig4, fig5, fig6,
                 nrow = 2, ncol = 2,
                 labels = c("(a)","(b)","(c)","(d)")
)
fig
ggsave(plot = fig, filename = "./visualizations/final_figures/FigureSY.tif",
       device = "tiff", height = 8, width = 15, units = "in", bg = "white")

