# List of R packages to install for nps_ndep_remeasure code repository
# Author: Mary Lofton
# Date: 22DEC25

# Purpose: install all packages needed to run the NPS N dep remeasure models
# and analysis

# install packages
install.packages(c("tidyverse", "lubridate", "rjags", "tidybayes", "bayesplot",
                   "furrr", "sf", "viridis", "ggthemes", "ggpubr", "rnaturalearth",
                   "scales", "arrow", "stringr", "metR", "grid", "alphahull",
                   "ggh4x","ggpattern","ggpmisc","zen4R"))

# load packages
library(tidyverse)
library(rjags)
library(tidybayes)
library(bayesplot)
library(furrr)
library(lubridate)
library(sf)
library(viridis)
library(ggthemes)
library(ggpubr)
library(rnaturalearth)
library(scales)
library(arrow)
library(stringr)
library(metR)
library(grid)
library(alphahull)
library(ggh4x)
library(ggpattern)
library(ggpmisc)
library(zen4R)
