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
library(scales)

og_df <- read_csv("./data/McDonnell_etal_InPrep_TreeData_2024_10_11.csv", show_col_types = FALSE) %>%
  dplyr::filter(!common_name %in% c("Douglas-fir","western hemlock")) 
#filter(common_name %in% c("eastern cottonwood"))

df <- read_csv("./data/processed_data.csv")

#'A. plots by ecoregion

# Set the path to your downloaded and unzipped shapefile
shapefile_path <- "./data/na_cec_eco_L1/NA_CEC_Eco_Level1.shp"

# Read the shapefile into R
na_ecoregions <- st_read(shapefile_path)

# Filter for the United States ecoregions.
# The 'COUNTRY' or 'NA_L1_CODE' column can be used for filtering.
# We will combine US ecoregions by using the 'NA_L1_CODE' column, which uniquely identifies each Level I ecoregion across North America.
us_ecoregions <- na_ecoregions %>%
  filter(grepl("United States", NA_L1NAME) ) # A quick filter, adjust as needed

us_boundaries <- ne_states(country = "United States of America", returnclass = "sf") %>%
  filter(!iso_3166_2 == "US-AK")
crs_us <- st_crs(us_boundaries)
crs_na <- st_crs(na_ecoregions)
us_sf_transformed <- st_transform(us_boundaries, crs = crs_na)
us_ecoregions_sf <- st_intersection(na_ecoregions, us_sf_transformed) %>%
  filter(!NA_L1NAME == "WATER")

plots <- read_csv("./data/plot_ecoregions.csv") %>%
  filter(!level1_ecoregion %in% c("WATER",NA) & plot_ID %in% df$plot_ID)

plotting_points <- st_as_sf(plots,
                            coords = c("lon", "lat"),
                            crs = 4269) #NAD 83 as specified in FIA docs

# Create the plot
fig2_a <- ggplot() +
  # Add the ecoregions layer, using the ecoregion name for the fill color
  geom_sf(data = us_ecoregions_sf, aes(fill = NA_L1NAME), color = "white", linewidth = 0.1) +
  geom_sf(data = plotting_points, aes(fill = level1_ecoregion), color = "black", shape = 21) +
  # Set the title and remove axis labels
  labs(
    title = "Measurement plots by USA Level I Ecoregion",
    fill = "Ecoregion"
  ) +
  
  # Remove the gridlines and axis text for a cleaner map
  theme_minimal() +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank()
  )

# B. number of remeasurements by species

dat_b <- df %>%
  select(plot_ID, num_intervals) %>%
  distinct(.) %>%
  count(num_intervals) %>% 
  mutate(perc = n / sum(n) * 100) 

fig2_b <- ggplot(dat_b, aes(x = num_intervals, y = perc)) + 
  geom_bar(stat = "identity") +
  theme_classic() +
  ylab("Percent of total plots in dataset") +
  xlab("Number of measurement intervals")

# C. stacked timeseries of number of measurements per year by species

spp_df <- og_df %>%
  select(common_name, species) %>%
  distinct(.)

dat_c <- left_join(df, spp_df, by = "common_name") %>%
  mutate(year = year(date_m2)) %>%
  group_by(species, year) %>%
  summarize(num_measures = length(unique(tree_ID)))

fig2_c <- ggplot(data = dat_c)+
  geom_area(aes(x = year, y = num_measures, group = species,
                fill = species), color = "black")+
  theme_classic()+
  xlab("Study year")+
  ylab("Number of trees measured")+
  labs(fill = "Species")+
  scale_fill_colorblind()

# Assemble figure

p1 <- ggarrange(fig2_a, 
                ggarrange(fig2_b, fig2_c, ncol = 2, labels = c("(b)","(c)"),
                          hjust = c(0.03, 0.03), widths = c(1, 1.5)), 
                nrow = 2,
                labels = "(a)",
                hjust = c(0.03)
)
p1
ggsave(plot = p1, filename = "./visualizations/final_figures/Figure2.tif",
       device = "tiff", height = 7, width = 10, units = "in", bg = "white")

#### END OF CODE FOR FINAL FIGURE INCLUDED IN MANUSCRIPT


#### additional plots for CLAD presentation

# Generate the default discrete color palette for 3 colors
default_colors <- hue_pal()(3)

# Extract the first two colors
first_two_colors <- default_colors[1:2]

# Print the colors
print(first_two_colors)

# Create the plot
clad_map0 <- ggplot() +
  # Add the ecoregions layer, using the ecoregion name for the fill color
  geom_sf(data = us_ecoregions_sf, aes(fill = NA_L1NAME), color = "black", linewidth = 0.3) +
  #geom_sf(data = plotting_points, aes(fill = level1_ecoregion), color = "black", shape = 21) +
  # Remove the gridlines and axis text for a cleaner map
  theme_minimal() +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank()
  ) +
  theme(legend.position = "none")+
  scale_fill_manual(values = c(first_two_colors[1],"white","white","white",
                               "white","white",
                               "white","white","white","white"))
clad_map0
ggsave(plot = clad_map0, filename = "./visualizations/clad_map0.tif",
       device = "tiff", height = 4, width = 6, units = "in", bg = "white")

# Create the plot
clad_map <- ggplot() +
  # Add the ecoregions layer, using the ecoregion name for the fill color
  geom_sf(data = us_ecoregions_sf, aes(fill = NA_L1NAME), color = "black", linewidth = 0.3) +
  #geom_sf(data = plotting_points, aes(fill = level1_ecoregion), color = "black", shape = 21) +
  # Remove the gridlines and axis text for a cleaner map
  theme_minimal() +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank()
  ) +
  theme(legend.position = "none")+
  scale_fill_manual(values = c(first_two_colors[1],"white","white","white",
                                "white",first_two_colors[2],
                                "white","white","white","white"))
clad_map
ggsave(plot = clad_map, filename = "./visualizations/clad_map.tif",
       device = "tiff", height = 4, width = 6, units = "in", bg = "white")

# Create the plot
clad_map_nf <- ggplot() +
  # Add the ecoregions layer, using the ecoregion name for the fill color
  geom_sf(data = us_ecoregions_sf, aes(fill = NA_L1NAME), color = "black", linewidth = 0.3) +
  #geom_sf(data = plotting_points, aes(fill = level1_ecoregion), color = "black", shape = 21) +
  # Remove the gridlines and axis text for a cleaner map
  theme_minimal() +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank()
  ) +
  theme(legend.position = "none")+
  scale_fill_manual(values = c("white","white","white","white",
                               "white","gray",
                               "white","white","white","white"))
clad_map_nf
ggsave(plot = clad_map_nf, filename = "./visualizations/clad_map_nf.tif",
       device = "tiff", height = 4, width = 6, units = "in", bg = "white")
