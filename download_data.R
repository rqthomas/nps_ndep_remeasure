# Download raw data from Zenodo for use in the nps_ndep_remeasure code repository
# Author: Mary Lofton
# Date: 25FEB26

# Purpose: Download data from McDonnell (2026): https://doi.org/10.5281/zenodo.18701698

# Note: The total file size for download is large (473.8 MB) and may take a few minutes
# depending on your internet speed. You can increase the 'timeout' argument in the 
# 'download_zenodo' function (provided in seconds) if the download times out before it
# is complete.

library(zen4R)

doi <- "10.5281/zenodo.18701698"
download_zenodo(doi = doi, path = "./data/", timeout = 600) 
