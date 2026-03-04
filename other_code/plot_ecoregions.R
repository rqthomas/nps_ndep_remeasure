# Define ecoregions of tree locations
# Author: Mary Lofton
# Date: 12SEP25

library(sf)
library(rnaturalearth)

plot_ecoregions <- function(){
# verify this worked with plots

# Set the path to your downloaded and unzipped shapefile
shapefile_path <- "na_cec_eco_L1/NA_CEC_Eco_Level1.shp"

# Read the shapefile into R
na_ecoregions <- st_read(shapefile_path)

# Filter for the United States ecoregions.
# The 'COUNTRY' or 'NA_L1_CODE' column can be used for filtering.
# We will combine US ecoregions by using the 'NA_L1_CODE' column, which uniquely identifies each Level I ecoregion across North America.
us_ecoregions <- na_ecoregions %>%
  filter(grepl("United States", NA_L1NAME) ) # A quick filter, adjust as needed

# You can examine the ecoregion names to decide how to filter
# unique(na_ecoregions$NA_L1NAME) 

# For a more robust filtering, you might consider merging with state boundaries.
# Here's an example to isolate US ecoregion parts more directly:
us_boundaries <- ne_states(country = "United States of America", returnclass = "sf") %>%
  filter(!iso_3166_2 == "US-AK")
crs_us <- st_crs(us_boundaries)
crs_na <- st_crs(na_ecoregions)
us_sf_transformed <- st_transform(us_boundaries, crs = crs_na)
us_ecoregions_sf <- st_intersection(na_ecoregions, us_sf_transformed)

plots <- read_csv("./data/plot_ecoregions.csv") %>%
  filter(!level1_ecoregion %in% c("WATER",NA))

plotting_points <- st_as_sf(plots,
                         coords = c("lon", "lat"),
                         crs = 4269) #NAD 83 as specified in FIA docs

# Create the plot
p <- ggplot() +
  # Add the ecoregions layer, using the ecoregion name for the fill color
  geom_sf(data = us_ecoregions_sf, aes(fill = NA_L1NAME), color = "white", linewidth = 0.1) +
  geom_sf(data = plotting_points, aes(fill = level1_ecoregion), color = "black", shape = 21) +
  # Set the title and remove axis labels
  labs(
    title = "USA Level I Ecoregions",
    fill = "Ecoregion"
  ) +
  
  # Remove the gridlines and axis text for a cleaner map
  theme_minimal() +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank()
  )

return(p)

}
