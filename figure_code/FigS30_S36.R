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
library(grid)
library(scales)
library(plotly)

source("./other_code/Fig6_supp_fig_function.R")

og_df <- read_csv("./data/McDonnell_etal_InPrep_TreeData_2024_10_11.csv", show_col_types = FALSE) %>%
  dplyr::filter(!common_name %in% c("Douglas-fir","western hemlock")) 
#filter(common_name %in% c("eastern cottonwood"))

spp_df <- og_df %>%
  select(species, common_name) %>%
  distinct(.)

# need Dep_Nhistoric_SO and c for map figures and SI version of Fig. 6
df <- read_csv("./data/processed_data.csv") %>%
  left_join(., spp_df, by = "common_name")

# Read in Class I data
# Set the path to your downloaded and unzipped shapefile
shapefile_path <- "./data/npsClassI_20071119/npsClassI.shp"

# Read the shapefile into R and limit to CONUS
nps_class1 <- st_read(shapefile_path)
crs_c1 <- st_crs(nps_class1)

us_boundaries <- ne_states(country = "United States of America", returnclass = "sf") %>%
  filter(!(iso_3166_2 %in% c("US-AK","US-HI")))
crs_us <- st_crs(us_boundaries)

us_sf_transformed <- st_transform(us_boundaries, crs = crs_c1)
conus_c1_sf <- st_intersection(nps_class1, us_sf_transformed) 

# Which C1 area do the most trees fall into for each species?

spp <- unique(spp_df$species)
final_df <- NULL

for(i in 1:length(spp)){

  plots_location <- df %>%
    select(plot_ID, lat, lon, species, tree_ID) %>%
    filter(species == spp[i]) %>%
    distinct(.) 

  plots_sf <- st_as_sf(
    x = plots_location,
    coords = c("lon", "lat"), # Specify the longitude column first
    crs = crs_c1
  )

  c1_fia_intersect <- conus_c1_sf |> 
    st_filter(y = plots_sf, .predicate = st_intersects)

  c1_areas <- unique(c1_fia_intersect$ID)
  
  if(length(c1_areas) == 0) next
  
  for(j in 1:length(c1_areas)){
    
    curr_c1 <- c1_fia_intersect %>%
      filter(ID == c1_areas[j])
    
    fia_in_curr_c1 <- st_filter(plots_sf, curr_c1)
    
    temp_df <- data.frame(species = spp[i],
                          c1_area = c1_areas[j],
                          num_trees = length(unique(fia_in_curr_c1$tree_ID)),
                          num_plots = length(unique(fia_in_curr_c1$plot_ID)))
    
    if(i == 1 & j == 1){
      final_df <- temp_df
    } else {
      final_df <- bind_rows(final_df, temp_df)
    }
    
  }

}

final_df <- final_df %>%
  group_by(species) %>%
  arrange(species, desc(num_trees)) %>%
  ungroup()

c1_most_trees <- final_df %>%
  group_by(species) %>%
  filter(num_trees == max(num_trees, na.rm = TRUE))

# Now plot

# first get Fig. 6 panels
Fig6_supp_plots <- Fig6_supp_fig_function()
fig6_panels_to_use <- Fig6_supp_plots[c(1:5,7,8)]

usa_coordinates <- map_data("state")

c1_to_plot <- conus_c1_sf %>%
  filter(ID %in% c1_most_trees$c1_area)

spp_in_c1 <- unique(c1_most_trees$species)

for(i in 1:length(spp_in_c1)){
  
  print(spp_in_c1[i])
  
  # some data wrangling
  focal_c1 <- c1_most_trees %>%
    filter(species == spp_in_c1[i])
  
  focal_c1_poly <- c1_to_plot %>%
    filter(ID == focal_c1$c1_area)
  
  bbox_c1 <- st_bbox(focal_c1_poly)
  
  if(i %in% c(1,2,6)){
    bbox_c1[1] = bbox_c1[1]-0.1
    bbox_c1[3] = bbox_c1[3]+0.1
  }
  
  fill_lab <- paste0(focal_c1_poly$Name, " (", focal_c1_poly$State, ", USA)")
  
  focal_c1_poly$fill_label <- fill_lab
  
  # then make locater map
  p0 <- ggplot() +
    geom_map(
      data = usa_coordinates, map = usa_coordinates,
      aes(long, lat, map_id = region),
      color = "black", fill = "white")+
    geom_sf(data = focal_c1_poly, aes(fill = fill_label), color = "black", linewidth = 0.1)+
    theme_bw()+
    ylab("")+
    xlab("")+
    scale_fill_manual(name = "", values = c("orange"))+
    ylim(c(bbox_c1[2] - 2.5,bbox_c1[4] + 2.5))+
    xlim(c(bbox_c1[1] - 5,bbox_c1[3] + 5))+
    theme(legend.position = "bottom",
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          legend.text = element_text(face = "bold"),
          panel.background = element_rect(fill = "lightgray"),
          legend.box.spacing = unit(0, "pt"), # Reduce space between legend and plot panel
          legend.margin = margin(0, 0, 0, 0, "pt"),
          title = element_text(face = "bold"),
          plot.margin = margin(0,0.5,0,0, "cm"))+
    ggtitle("b. focal Class I area",
            subtitle = paste0("(",focal_c1$num_trees," trees across ",focal_c1$num_plots," plots)"))
  p0
  
  # then plot N deposition
  
  # get lims from Fig. 6 panel
  panel <- ggplot_build(fig6_panels_to_use[[i]])
  
  # Extract x-axis limits (xmin, xmax)
  ante_limits <- panel$layout$panel_params[[1]]$x.range
  
  # Extract y-axis limits (ymin, ymax)
  st_limits <- panel$layout$panel_params[[1]]$y.range
  
  mean_Ndep0 <- df %>%
    select(plot_ID, lat, lon, Dep_Nhistoric_SO, Dep_Ndiff_SO, date_m2, species) %>%
    filter(species == spp_in_c1[i]) %>%
    distinct(.) 
    
  mean_Ndep <- mean_Ndep0 %>%
    group_by(plot_ID, lat, lon, species) %>%
    summarize(Dep_Nhistoric_SO = mean(Dep_Nhistoric_SO, na.rm = TRUE),
              Dep_Ndiff_SO = mean(Dep_Ndiff_SO, na.rm = TRUE)) %>%
    ungroup()
  
  date_range_mean_Ndep <- c(min(mean_Ndep0$date_m2, na.rm = TRUE), max(mean_Ndep0$date_m2, na.rm = TRUE))
  
  mean_Ndep_sf <- st_as_sf(
    x = mean_Ndep,
    coords = c("lon", "lat"), # Specify the longitude column first
    crs = crs_c1
  )
  
  mean_fia_in_focal_c1 <- st_filter(mean_Ndep_sf, focal_c1_poly)

  p1 <- ggplot() +
    geom_map(
      data = usa_coordinates, map = usa_coordinates,
      aes(long, lat, map_id = region),
      color = "black", fill = "white")+
    geom_sf(data = focal_c1_poly, color = "black", linewidth = 0.1)+
    geom_sf(data = mean_fia_in_focal_c1, aes(color = Dep_Nhistoric_SO))+
    ylim(c(bbox_c1[2],bbox_c1[4]))+
    xlim(c(bbox_c1[1],bbox_c1[3]))+
    ggtitle(paste0("c. Mean measurements: ",year(date_range_mean_Ndep[1]),"-",year(date_range_mean_Ndep[2])))+
    ylab("")+
    xlab("")+
    scale_color_viridis(option = "A", name = expression(atop("residual \nantecedent \nN deposition ",(kg~N~ha^-1~yr^-1))), limits = ante_limits)+
    theme(panel.border = element_rect(fill = NA, color = "darkgreen", linewidth = 2, linetype = "solid"),
          axis.title.x = element_text(hjust = 1),
          title = element_text(face = "bold"),
          plot.margin = margin(0.5,0,0,0.5, "cm"),
          plot.background = element_rect(fill = "transparent", colour = NA),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = "lightcyan"))+
    labs(x = bquote(atop(bold("observed res. ante. N dep. range:"),~~~~~~~~~~~~~~~~~~~~~~~~.(round(min(mean_fia_in_focal_c1$Dep_Nhistoric_SO),2)) ~ "-" ~ .(round(max(mean_fia_in_focal_c1$Dep_Nhistoric_SO),2)) ~ "kg N"~ha^-1~yr^1)))
  p1
  
  p2 <- ggplot() +
    geom_map(
      data = usa_coordinates, map = usa_coordinates,
      aes(long, lat, map_id = region),
      color = "black", fill = "white")+
    geom_sf(data = focal_c1_poly, color = "black", linewidth = 0.1)+
    geom_sf(data = mean_fia_in_focal_c1, aes(color = Dep_Nhistoric_SO))+
    ylim(c(bbox_c1[2],bbox_c1[4]))+
    xlim(c(bbox_c1[1],bbox_c1[3]))+
    ggtitle(paste0("d. Mean measurements: ",year(date_range_mean_Ndep[1]),"-",year(date_range_mean_Ndep[2])))+
    ylab("")+
    xlab("")+
    scale_color_viridis(option = "A", name = expression(atop("residual \nshort-term \nchange in \nN deposition ",(kg~N~ha^-1~yr^-1))), limits = st_limits)+
    theme(panel.border = element_rect(fill = NA, color = "darkgreen", linewidth = 2, linetype = "solid"),
          axis.title.x = element_text(hjust = 1),
          title = element_text(face = "bold"),
          plot.margin = margin(0.5,0,0,0.5, "cm"),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = "lightcyan"))+
    labs(x = bquote(atop(bold("obs. res. ST change in N dep. range:"),~~~~~~~~~~~~~~~~~~~~~~~~~~.(round(min(mean_fia_in_focal_c1$Dep_Ndiff_SO),2)) ~ "-" ~ .(round(max(mean_fia_in_focal_c1$Dep_Ndiff_SO),2)) ~ "kg N"~ha^-1~yr^1)))
  p2
  
  last_Ndep <- df %>%
    select(plot_ID, lat, lon, Dep_Nhistoric_SO, Dep_Ndiff_SO, date_m2, species) %>%
    filter(species == spp_in_c1[i]) %>%
    distinct(.) %>%
    group_by(plot_ID) %>%
    filter(year(date_m2) >= 2012) %>%
    ungroup()
  
  date_range_last_Ndep <- range(last_Ndep$date_m2, na.rm = TRUE)
  
  last_Ndep_sf <- st_as_sf(
    x = last_Ndep,
    coords = c("lon", "lat"), # Specify the longitude column first
    crs = crs_c1
  )
  
  last_fia_in_focal_c1 <- st_filter(last_Ndep_sf, focal_c1_poly)
  
  p3 <- ggplot() +
    geom_map(
      data = usa_coordinates, map = usa_coordinates,
      aes(long, lat, map_id = region),
      color = "black", fill = "white")+
    geom_sf(data = focal_c1_poly, color = "black", linewidth = 0.1)+
    geom_sf(data = last_fia_in_focal_c1, aes(color = Dep_Nhistoric_SO))+
    ylim(c(bbox_c1[2],bbox_c1[4]))+
    xlim(c(bbox_c1[1],bbox_c1[3]))+
    ggtitle(paste0("e. Recent measurements: ",year(date_range_last_Ndep[1]),"-",year(date_range_last_Ndep[2])))+
    ylab("")+
    xlab("")+
    scale_color_viridis(option = "A", name = expression(atop("residual \nantecedent \nN deposition ",(kg~N~ha^-1~yr^-1))), limits = ante_limits)+
    theme(panel.border = element_rect(fill = NA, color = "lightgreen", linewidth = 2, linetype = "solid"),
          axis.title.x = element_text(hjust = 1),
          title = element_text(face = "bold"),
          plot.margin = margin(0.5,0,0,0.5, "cm"),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = "lightcyan"))+
    labs(x = bquote(atop(bold("observed res. ante. N dep. range:"),~~~~~~~~~~~~~~~~~~~~~~~~.(round(min(last_fia_in_focal_c1$Dep_Nhistoric_SO),2)) ~ "-" ~ .(round(max(last_fia_in_focal_c1$Dep_Nhistoric_SO),2)) ~ "kg N"~ha^-1~yr^1)))
  p3
  
  p4 <- ggplot() +
    geom_map(
      data = usa_coordinates, map = usa_coordinates,
      aes(long, lat, map_id = region),
      color = "black", fill = "white")+
    geom_sf(data = focal_c1_poly, color = "black", linewidth = 0.1)+
    geom_sf(data = last_fia_in_focal_c1, aes(color = Dep_Nhistoric_SO))+
    ylim(c(bbox_c1[2],bbox_c1[4]))+
    xlim(c(bbox_c1[1],bbox_c1[3]))+
    ggtitle(paste0("f. Recent measurements: ",year(date_range_last_Ndep[1]),"-",year(date_range_last_Ndep[2])))+
    ylab("")+
    scale_color_viridis(option = "A", name = expression(atop("residual \nshort-term \nchange in \nN deposition ",(kg~N~ha^-1~yr^-1))), limits = st_limits)+
    theme(panel.border = element_rect(fill = NA, color = "lightgreen", linewidth = 2, linetype = "solid"),
          axis.title.x = element_text(hjust = 1),
          title = element_text(face = "bold"),
          plot.margin = margin(0.5,0,0,0.5, "cm"),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = "lightcyan"))+
    labs(x = bquote(atop(bold("obs. res. ST change in N dep. range:"),~~~~~~~~~~~~~~~~~~~~~~~~~~.(round(min(last_fia_in_focal_c1$Dep_Ndiff_SO),2)) ~ "-" ~ .(round(max(last_fia_in_focal_c1$Dep_Ndiff_SO),2)) ~ "kg N"~ha^-1~yr^1)))
  p4
  
  # finally adjust Fig. 6 panel
  p00 <- fig6_panels_to_use[[i]]+
    annotate("rect", 
             xmin = round(min(mean_fia_in_focal_c1$Dep_Nhistoric_SO),2), 
             xmax = round(max(mean_fia_in_focal_c1$Dep_Nhistoric_SO),2), 
             ymin = round(min(mean_fia_in_focal_c1$Dep_Ndiff_SO),2), 
             ymax = round(max(mean_fia_in_focal_c1$Dep_Ndiff_SO),2),
             color = "darkgreen", fill = NA, linewidth = 1)+
    annotate("rect", 
             xmin = round(min(last_fia_in_focal_c1$Dep_Nhistoric_SO),2), 
             xmax = round(max(last_fia_in_focal_c1$Dep_Nhistoric_SO),2), 
             ymin = round(min(last_fia_in_focal_c1$Dep_Ndiff_SO),2), 
             ymax = round(max(last_fia_in_focal_c1$Dep_Ndiff_SO),2),
             color = "lightgreen", fill = NA, linewidth = 1)+
    theme(plot.margin = margin(0.5,0,0,0.5, "cm"))

  p <- ggarrange(ggarrange(p00, p0, ncol = 2, nrow = 1, widths = c(2,1)), 
                 ggarrange(p1, p2, p3, p4, ncol = 2, nrow = 2, align = "hv"),
                 ncol = 1, nrow = 2, heights = c(1.3,2))
  p
 
  ggsave(p, filename = paste0("./visualizations/final_figures/classI_all/",spp_in_c1[i],".png"),
        height = 10, width = 12, units = "in", bg = "white")
 
  
}
  
###############################################################################

plots_location <- df %>%
  select(plot_ID, lat, lon, Dep_Nhistoric_SO, Dep_Ndiff_SO, date_m2, species, tree_ID) %>%
  distinct(.) %>%
  group_by(plot_ID) %>%
  filter(date_m2 == max(date_m2, na.rm = TRUE)) %>%
  ungroup()
plots_sf <- st_as_sf(
  x = plots_location,
  coords = c("lon", "lat"), # Specify the longitude column first
  crs = crs_c1
)

c1_fia_intersect <- conus_c1_sf |> 
  st_filter(y = plots_sf, .predicate = st_intersects)

#'A. map of mean N deposition during measurement period (including 15 yr prior to 1st measure)

usa_coordinates <- map_data("state")

c1_areas <- unique(c1_fia_intersect$ID)
curr_c1 <- c1_fia_intersect %>%
  filter(ID == "shen")

fia_in_curr_c1 <- st_filter(plots_sf, curr_c1)

bbox_c1 <- st_bbox(curr_c1)

p <- ggplot() +
  geom_map(
    data = usa_coordinates, map = usa_coordinates,
    aes(long, lat, map_id = region),
    color = "white", fill = "white")+
  geom_sf(data = curr_c1, color = "black", linewidth = 0.1)+
  geom_sf(data = fia_in_curr_c1, aes(color = Dep_Nhistoric_SO))+
  ylim(c(bbox_c1[2],bbox_c1[4]))+
  xlim(c(bbox_c1[1],bbox_c1[3]))+
  ggtitle(curr_c1$Name)
p

  geom_point(
    data = map_data,
    aes(lon, lat, color = Dep_Nmean),
    shape = 16, alpha = 1, size = 1
  ) +
  scale_color_viridis(option = "H")+
  xlab("")+
  ylab("")+
  theme_classic()+
  labs(color = expression(paste("mean N deposition (kg N ", ha^-1," ",y^-1,")")))+
  guides(color = guide_colorbar(title.position = "bottom", title.hjust = 0.0))+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        legend.background = element_rect(fill = "transparent"))+
  theme(legend.position = "bottom",
        legend.direction = "horizontal", plot.margin = unit(c(0.1,0,0.1,0), "pt"))+
  ylim(c(min(usa_coordinates$lat),max(usa_coordinates$lat)))


fig3_a


# Assemble figure

p2 <- ggarrange(fig3_a, fig3_b,
                nrow = 2, ncol = 1,
                labels = c("(a)","(b)")
)
p2
ggsave(plot = p2, filename = "./visualizations/final_figures/Figure3.tif",
       device = "tiff", height = 8, width = 7.5, units = "in", bg = "white")



