# Explore data to inform model structure
# Author: Mary Lofton
# Date: 20JAN25

# Purpose: determine how many trees there are per plot and whether species
# are experiencing increases and decreases in N dep

# load packages
library(tidyverse)
library(lubridate)
library(sf)
library(viridis)
library(ggthemes)
library(ggpubr)

# load data
df <- read_csv("./data/McDonnell_etal_InPrep_TreeData_2024_10_11.csv", show_col_types = FALSE)
colnames(df)

# how many trees per plot for each species?
unique(df$common_name)
length(unique(df$plot_ID))
df1 <- df %>%
  group_by(common_name, plot_ID) %>%
  summarize(num_per_plot = n_distinct(tree_ID)) 

num_per_plot_fig <- df1 %>%
  ggplot(aes(x = num_per_plot, group = common_name)) +
  geom_histogram() +
  facet_wrap(facets = vars(common_name), scales = "free_y") + 
  theme_bw()
num_per_plot_fig

median_per_plot <- df1 %>%
  ungroup() %>%
  group_by(common_name) %>%
  summarize(med_per_plot = median(num_per_plot, na.rm = TRUE))

num_plots_species <- df1 %>%
  ungroup() %>%
  group_by(common_name) %>%
  summarize(total_plots = sum(n_distinct(plot_ID)))

num_plots_multiple <- df1 %>%
  ungroup() %>%
  filter(!num_per_plot == 1) %>%
  group_by(common_name) %>%
  summarize(num_plots_multiple = sum(n_distinct(plot_ID)))

prop_plots_multiple <- left_join(num_plots_species, num_plots_multiple,
                                 by = "common_name") %>%
  mutate(prop_plots_multiple = num_plots_multiple / total_plots)

# how many trees per species experience both increases and decreases in Ndep?
colnames(df)

df2 <- df %>%
  select(common_name, tree_ID, interval_no, Dep_N) %>%
  group_by(tree_ID) %>%
  mutate(num_intervals = max(interval_no, na.rm = TRUE))

num_intervals_fig <- df2 %>%
  ggplot(aes(x = num_intervals, group = common_name)) +
  geom_histogram() +
  facet_wrap(facets = vars(common_name), scales = "free_y") + 
  theme_bw()
num_intervals_fig

inc_dec <- df2 %>%
  filter(num_intervals > 1) %>%
  select(-num_intervals) %>%
  pivot_wider(names_from = interval_no, values_from = Dep_N) %>%
  rename(interval_1 = `1`,
         interval_2 = `2`,
         interval_3 = `3`,
         interval_4 = `4`) %>%
  mutate(inc_1_2 = interval_2 - interval_1,
         inc_2_3 = interval_3 - interval_2,
         inc_3_4 = interval_4 - interval_3) %>%
  select(-interval_1, -interval_2, -interval_3, -interval_4) %>%
  pivot_longer(inc_1_2:inc_3_4, names_to = "interval_span",
               values_to = "delta_ndep") %>%
  filter(!is.na(delta_ndep)) %>%
  ungroup()

min <- inc_dec %>% filter(delta_ndep == min(delta_ndep))
max <- inc_dec %>% filter(delta_ndep == max(delta_ndep))

inc_dec_fig <- inc_dec %>%
  ggplot(aes(x = delta_ndep, group = common_name)) +
  geom_density(color = "black", fill = "white") +
  geom_vline(xintercept = 0, color = "blue") +
  xlim(c(-5,5)) +
  facet_wrap(facets = vars(common_name), scales = "free_y") + 
  theme_bw()
inc_dec_fig

baseline_Ndep <- df %>%
  select(common_name, tree_ID, interval_no, Dep_N15, Dep_N, date_m2, date_m1) %>%
  filter(interval_no == 1) %>%
  mutate(dt = as.numeric(date_m2 - date_m1)/365) %>%
  mutate(total_years = 15 + dt) %>%
  mutate(Dep_Nbaseline = (Dep_N15*total_years - Dep_N*dt) / 15) %>%
  select(tree_ID, Dep_Nbaseline)
hist(baseline_Ndep$Dep_Nbaseline)

delta_Ndep <- df %>%
  select(common_name, tree_ID, interval_no, Dep_N) %>%
  left_join(baseline_Ndep, by = "tree_ID") %>%
  mutate(Dep_Ndelta = Dep_N - Dep_Nbaseline) %>%
  filter(!common_name %in% c("Douglas-fir","western hemlock")) 

plot_delta_Ndep <- ggplot(data = delta_Ndep, aes(x = Dep_Ndelta,
                                                 group = interval_no, fill = as.factor(interval_no)))+
  geom_density()+
  facet_wrap(facets = vars(common_name), scales = "free_y")+
  xlim(c(-15,15))+
  geom_vline(xintercept = 0)+
  theme_bw()
plot_delta_Ndep  

check_intervals <- df %>%
  select(interval_no, date_m1, date_m2) %>%
  pivot_longer(date_m1:date_m2, names_to = "m1_m2", values_to = "datetime")

plot_interval_times <- ggplot(data = check_intervals, 
                              aes(x = datetime, group = m1_m2,
                                  fill = m1_m2))+
  geom_histogram()+
  facet_wrap(facets = vars(interval_no))+
  theme_bw()
plot_interval_times

int4_trees <- df %>%
  group_by(common_name, tree_ID) %>%
  summarise(count = n()) %>%
  filter(count == 4) %>%
  group_by(common_name) %>%
  summarise(num_int4_trees = n())

# additional maps to address collaborator questions 17APR25

og_df <- read_csv("./data/McDonnell_etal_InPrep_TreeData_2024_10_11.csv", show_col_types = FALSE) %>%
  dplyr::filter(!common_name %in% c("Douglas-fir","western hemlock")) 
#filter(common_name %in% c("eastern cottonwood"))

focal_df <- og_df %>%
  select(species, common_name, plot_ID, tree_ID, interval_no, Dep_N, Dep_Noxi, Dep_Nred, 
         Dep_S, MAT, MAP, date_m2, date_m1, AG_carbon_pYear, AG_carbon_m1, 
         AG_carbon_m2, subp_BA_GT_m1, live_m2, lat, lon) 

baseline_vars <- focal_df %>%
  select(plot_ID, date_m1, date_m2, Dep_N, Dep_Noxi, Dep_Nred, Dep_S, MAT, MAP) %>%
  distinct(.) %>%
  group_by(plot_ID) %>%
  summarize(Dep_Nbaseline = mean(Dep_N, na.rm = TRUE),
            Dep_Noxibaseline = mean(Dep_Noxi, na.rm = TRUE),
            Dep_Nredbaseline = mean(Dep_Nred, na.rm = TRUE),
            Dep_Sbaseline = mean(Dep_S, na.rm = TRUE),
            MAT_baseline = mean(MAT, na.rm = TRUE),
            MAP_baseline = mean(MAP, na.rm = TRUE)) %>%
  ungroup()

focal_df2 <- left_join(focal_df, baseline_vars, by = "plot_ID") %>%
  group_by(plot_ID) %>%
  mutate(Dep_Ndelta = Dep_N - Dep_Nbaseline,
         Dep_Noxidelta = Dep_Noxi - Dep_Noxibaseline,
         Dep_Nreddelta = Dep_Nred - Dep_Nredbaseline,
         Dep_Sdelta = Dep_S - Dep_Sbaseline,
         MAP_delta_dm = (MAP - MAP_baseline) * 0.01,
         MAT_delta = MAT - MAT_baseline) %>%
  ungroup() 

live_tree_ids <- focal_df2 |> 
  summarise(count = n(),
            sum = sum(live_m2), .by = tree_ID) |> 
  dplyr::filter(count >= sum) |> 
  pull(tree_ID)

focal_df3 <- focal_df2 |> 
  dplyr::filter(tree_ID %in% live_tree_ids) |> 
  group_by(tree_ID) |> 
  tidyr::fill(subp_BA_GT_m1, .direction = "down") |> 
  ungroup() 

df <- focal_df3 %>%
  dplyr::filter(complete.cases(.)) %>%
  group_by(tree_ID) %>%
  mutate(num_intervals = max(interval_no, na.rm = TRUE)) %>%
  ungroup() 

species_guide <- df %>%
  select(species, common_name) %>%
  distinct(.)

usa_coordinates <- map_data("state")
map_data <- df %>%
  select(plot_ID, lat, lon, Dep_Nbaseline, Dep_Sbaseline,
         MAT_baseline, MAP_baseline, Dep_Noxibaseline, Dep_Nredbaseline,
         num_intervals) %>%
  mutate(MAP_baseline_dm = 0.01*MAP_baseline) %>%
  distinct(.)

my.cols <- colorblind_pal()(8)
names(my.cols) <- unique(as.factor(df$num_intervals))

all_int_map <- ggplot() +
  geom_map(
    data = usa_coordinates, map = usa_coordinates,
    aes(long, lat, map_id = region),
    color = "black", fill = "white")+
  geom_point(
    data = map_data,
    aes(lon, lat, color = as.factor(num_intervals)),
    shape = 16, alpha = 1, size = 1
  ) +
  scale_color_manual(values = my.cols)+
  xlab("")+
  ylab("")+
  theme_classic()+
  labs(color = "number of measurement intervals (re-measurements)")+
  theme(legend.position = "bottom",
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.line.x = element_blank(),
        axis.line.y = element_blank())
all_int_map

int_2_df <- map_data %>%
  filter(num_intervals >= 2)

int_2_map <- ggplot() +
  geom_map(
    data = usa_coordinates, map = usa_coordinates,
    aes(long, lat, map_id = region),
    color = "black", fill = "white")+
  geom_point(
    data = int_2_df,
    aes(lon, lat, color = as.factor(num_intervals)),
    shape = 16, alpha = 1, size = 1
  ) +
  scale_color_manual(values = my.cols)+
  xlab("")+
  ylab("")+
  theme_classic()+
  labs(color = "number of measurement intervals (re-measurements)")+
  theme(legend.position = "bottom",
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.line.x = element_blank(),
        axis.line.y = element_blank())
int_2_map

for(i in 1:length(species_guide$species)){
  map_dat <- df %>%
    filter(num_intervals >=2,
           species == species_guide$species[i]) %>%
    select(plot_ID, lat, lon, num_intervals, species) %>%
    distinct(.)
  
  species_int_2_map <- ggplot() +
    geom_map(
      data = usa_coordinates, map = usa_coordinates,
      aes(long, lat, map_id = region),
      color = "black", fill = "white")+
    geom_point(
      data = map_dat,
      aes(lon, lat, color = as.factor(num_intervals)),
      shape = 16, alpha = 1, size = 1
    ) +
    scale_color_manual(values = my.cols)+
    xlab("")+
    ylab("")+
    ggtitle(paste(species_guide$species[i],"(n sites = ",nrow(map_dat),")"))+
    theme_classic()+
    labs(color = "number of measurement intervals (re-measurements)")+
    theme(legend.position = "bottom",
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          axis.line.x = element_blank(),
          axis.line.y = element_blank())
  species_int_2_map
  ggsave(species_int_2_map, filename = paste0("./visualizations/site_map_",species_guide$species[i],".png"),
         device = "png")
}

# percent of sites vs. num remeasurements

quick_df <- read_csv("./data/processed_data.csv") %>%
  select(plot_ID, num_intervals) %>%
  distinct(.) %>%
  count(num_intervals) %>% 
  mutate(perc = n / sum(n) * 100) 

ggplot(quick_df, aes(x = num_intervals, y = perc)) + 
  geom_bar(stat = "identity") +
  theme_classic() +
  ylab("percent of total plots in dataset") +
  xlab("number of measurement intervals)")
  
