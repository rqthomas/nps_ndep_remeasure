# Code repository for Lofton et al. 2026, "Contrasting effects of short- versus long-term responses to nitrogen deposition in U.S. tree species"
Mary E. Lofton<sup>1</sup>, Michael R. Schwob<sup>2</sup>, Michael D. Bell<sup>3</sup>, Christopher M. Clark<sup>4</sup>, Emmi Felker-Quinn<sup>3</sup>, Todd McDonnell<sup>5</sup>, R. Quinn Thomas<sup>1</sup> 

*<sup>1</sup>Center for Ecosystem Forecasting, Virginia Tech, Blacksburg, VA, USA*<br>
*<sup>2</sup>Department of Statistics, Virginia Tech, Blacksburg, VA, USA*<br>
*<sup>3</sup>US National Park Service Air Resources Division, Denver, CO, USA*<br> 
*<sup>4</sup>US Environmental Protection Agency Office of Research and Development, Washington DC, USA*<br>
*<sup>5</sup>E&S Environmental, Corvallis, OR, USA* 

### Accessing data associated with this repository:
**Due to GitHub file size limitations, data associated with this manuscript are published in the Zenodo repository:**

*McDonnell, T. (2026). Tree aboveground biomass, basal area, atmospheric deposition, soil, and climate data for sugar maple (Acer saccharum), paper birch (Betula papyrifera), tuliptree (Liriodendron tulipifera), red spruce (Picea rubens), ponderosa pine (Pinus ponderosa), cottonwood (Populus deltoides), quaking aspen (Populus tremuloides), and black cherry (Prunus serotina) measured in United States Department of Agriculture Forest Inventory Analysis plots from 1997-2022. [Data set]. Zenodo.*
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18701698.svg)](https://doi.org/10.5281/zenodo.18701698)

Run the `download_data.R` script to access and download the data onto your local computer for use with other repository files.

### Guide to files and folders:

**data**: contains all data used for the project<sup>a</sup>

**experiments**: contains results from modeling experiments, including JAGS model .txt files, model output as parquet files, and default model output visualizations

**figure_code**: code used to generate figures used either in the main manuscript or the supplement

**install.R**: R code to install and load packages (libraries) that are needed to run the models and generate figures

**modeling_code**: workflow code used to run different modeling experiments

**other_code**: code used, e.g., for data exploration or figures for presentations

**visualizations**: all visualizations other than the default model output visualizations housed in the `./experiments/` folder

<sup>a</sup>*The raw data is not provided in this code repository because it is too big to upload to GitHub; see the McDonnell (2026) Zenodo publication (https://doi.org/10.5281/zenodo.18701697) to access the raw data. If you wish to work with the processed data (e.g., orthogonalized N deposition values) used as model inputs for this project but do not wish to run any models or generate any figures, you can run: 1) the `download_data.R` script to download the raw data onto your local computer; then 2) the `modeling_code/generate_processed_data.R` in this repository to generate processed data from the raw data.*

### Guide for reviewers:

#### Recommended workflow for review

1. Run `install.R` to be sure you have the necessary libraries loaded.
2. Run `download_data.R` to download the raw data from Zenodo onto your computer. If you are working in an RProject or have set your working directory to be the directory containing all files associated with this repository, the data will automatically be downloaded into the `./data/` subdirectory. Otherwise, you will need to move the data into this subdirectory after download. Be sure to adjust the download directory as appropriate for your computer.
3. Because the Bayesian hierarchical models are compute-intensive (we used a high-performance cluster to complete model runs for the manuscript), we recommend to run the example model for 100 individuals for one tree species, located in `modeling_code/run_model_example.R` to review the modeling workflow. This script will load USFS FIA data from our data publication, generate a `.csv` file of processed data which will be needed to generate the final figures for the manuscript, and do an example run of the Bayesian hierarchical model for 100 individuals of *Prunus serotina*.
4. Complete all analyses and generate final figures in the manuscript using the scripts in the `figure_code` folder. The scripts in this folder will use the `processed_data.csv` that you generated in Step 2 as well as the actual model output used in the manuscript from the `ortho_log_t_interaction_01JAN26` folder. Specifically, to generate the manuscript figures (as opposed to supplemental figures), run:
    a. `Fig1.R`
    b. `Fig2.R`
    c. `Fig3.R`
    d. `Fig4.R`
    e. `Fig5.R`
    f. `Fig6.R`
    g. `Fig7.R`
    
#### Funding
Funding for this project was provided by the U.S. National Park Service (task agreement P23AC01502).

#### Disclaimer text
The views expressed here are those of the authors and do not necessarily represent the views or policies of the U.S. Environmental Protection Agency or any other federal agency.