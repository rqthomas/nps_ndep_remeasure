# model_output_explore
# Author: Mary Lofton
# Date: 26FEB25

# load packages
library(tidyverse)
library(lubridate)
library(arrow)
library(ggpubr)
library(ggthemes)


# list files ----
out <- list.files("./experiments/delta_Ndep_only",pattern = ".parquet",
                  full.names = TRUE)

for(i in 1:length(out)){
  
  spp_name = str_split(out[i], pattern = "-")[[1]][2]
  model_name = str_split(out[i], pattern = "/")[[1]][3]
  temp <- read_parquet(file = out[i]) %>%
    mutate(spp_id = spp_name,
           model_id = model_name)
  
  if(i == 1){
    final <- temp
  } else {
    final <- bind_rows(final, temp)
  }
  
}

final <- final %>%
  mutate(spp_id = ifelse(spp_id == "yellow","yellow poplar",spp_id))

ggplot(data = final)+
  geom_density(aes(x = p5, group = model_id, color = model_id, fill = model_id),
               alpha = 0.5)+
  theme_classic()+
  facet_wrap(facets = vars(spp_id), scales = "free")+
  ggtitle("")

Ndep_param_sum <- final %>%
  select(model_id, p4, p5) %>%
  pivot_longer(p4:p5,names_to = "param", values_to = "param_value") %>%
  group_by(model_id, param) %>%
  summarize(mean_param_value = mean(param_value, na.rm = TRUE)) %>%
  pivot_wider(names_from = "param", values_from = "mean_param_value") %>%
  ggplot(aes(x = p4, y = p5, group = model_id, color = model_id))+
  geom_point(size = 2)+
  theme_bw()+
  geom_vline(xintercept = 0)+
  geom_hline(yintercept = 0)
Ndep_param_sum

# compare obs vs pred
source("./other_code/pred_vs_obs.R")
compare_pred_vs_obs <- pred_vs_obs(data = "./data/McDonnell_etal_InPrep_TreeData_2024_10_11.csv", 
                                    model_output_folder = "./experiments/new_delta_Ndep_only")
compare_pred_vs_obs$plot2 
compare_pred_vs_obs$plot1
rsq <- function(pred, obs){
  1 - (sum((obs - pred)^2, na.rm = TRUE) / sum((obs - mean(obs, na.rm = TRUE))^2, na.rm = TRUE))
}

r2 <- rsq(pred = compare_pred_vs_obs$df$pred, obs = compare_pred_vs_obs$df$AG_carbon_pYear)

# growth vs size
growth_vs_size(data = "./data/McDonnell_etal_InPrep_TreeData_2024_10_11.csv", 
            model_output_folder = "./experiments/delta_Ndep")

# growth vs Ndep
growth_vs_Ndep(data = "./data/McDonnell_etal_InPrep_TreeData_2024_10_11.csv", 
               model_output_folder = "./experiments/delta_Ndep")

# delta growth vs delta Ndep
delta_growth_vs_delta_Ndep(data = "./data/McDonnell_etal_InPrep_TreeData_2024_10_11.csv", 
               model_output_folder = "./experiments/delta_Ndep",
               experiment_name = "delta_Ndep_only")

# checking param distributions with Ndep only vs all env variables vs
# separating N into different species

# list files
out1 <- list.files("./experiments/new_delta_Ndep_only_saveTreeEffect",pattern = "mcmc.parquet",
                  full.names = TRUE)
out2 <- list.files("./experiments/delta_env_saveTreeEffect",pattern = "mcmc.parquet",
                   full.names = TRUE)
out3 <- list.files("./experiments/N_species_saveTreeEffect",pattern = "mcmc.parquet",
                   full.names = TRUE)
out <- c(out1,out2, out3)

for(i in 1:length(out)){
  
  spp_id = str_split(out[i], pattern = "-")[[1]][2]
  model_id = str_split(out[i], pattern = "/")[[1]][3]
  temp <- read_parquet(file = out[i]) %>%
    mutate(spp_id = spp_id,
           model_id = model_id)
  
  if(i == 1){
    final <- temp
  } else {
    final <- bind_rows(final, temp)
  }
  
}

final1 <- final %>%
  mutate(spp_id = ifelse(spp_id == "yellow","yellow poplar",spp_id)) %>%
  select(model_id, spp_id, p2, p3, p4, p5, p6, p7, p8) %>%
  pivot_longer(p2:p8, names_to = "parameter_name", values_to = "parameter_value") %>%
  mutate(parameter_name = ifelse(model_id %in% c("new_delta_Ndep_only_saveTreeEffect") & parameter_name == "p5","delta_Ndep_only",
                                 ifelse(model_id == "N_species_saveTreeEffect" & parameter_name == "p4","delta_Ndep_oxi",
                                        ifelse(model_id == "N_species_saveTreeEffect" & parameter_name == "p5","delta_Ndep_red",
                                               ifelse(model_id == "delta_env_saveTreeEffect" & parameter_name == "p5","delta_Ndep_env",parameter_name))))) %>%
  filter(parameter_name %in% c("delta_Ndep_only","delta_Ndep_oxi","delta_Ndep_red", "delta_Ndep_env"))

ggplot(data = final1)+
  geom_density(aes(x = parameter_value, group = parameter_name, 
                   color = parameter_name, fill = parameter_name),
               alpha = 0.5)+
  facet_wrap(facets = vars(spp_id), scales = "free")+
  theme_bw()

final2 <- final1 %>%
  filter(parameter_name %in% c("delta_Ndep_only","delta_Ndep_env"))

ggplot(data = final2)+
  geom_density(aes(x = parameter_value, group = parameter_name, 
                   color = parameter_name, fill = parameter_name),
               alpha = 0.5)+
  facet_wrap(facets = vars(spp_id), scales = "free")+
  theme_bw()

final3 <- final1 %>%
  filter(parameter_name %in% c("delta_Ndep_oxi","delta_Ndep_red","delta_Ndep_env"))

ggplot(data = final3)+
  geom_density(aes(x = parameter_value, group = parameter_name, 
                   color = parameter_name, fill = parameter_name),
               alpha = 0.5)+
  facet_wrap(facets = vars(spp_id), scales = "free")+
  theme_bw()

ggplot(data = final)+
  geom_density(aes(x = p4, group = model_id, color = model_id, fill = model_id),
               alpha = 0.5)+
  theme_classic()

ggplot(data = final)+
  geom_density(aes(x = p5, group = model_id, color = model_id, fill = model_id),
               alpha = 0.5)+
  theme_classic()+
  ggtitle("Delta Ndep only")+
  xlim(c(-0.12, 0.12))

# visualizing coefs for new terms added in 24APR25 ----

# for unpacking the plot effect experiment

# list files
out <- list.files("./experiments/unpacking_plot_effect",pattern = "mcmc.parquet",
                  full.names = TRUE)

for(i in 1:length(out)){
  
  spp_name = str_split(out[i], pattern = "-")[[1]][2]
  model_name = str_split(out[i], pattern = "/")[[1]][3]
  temp <- read_parquet(file = out[i]) %>%
    mutate(spp_id = spp_name,
           model_id = model_name)
  
  if(i == 1){
    final <- temp
  } else {
    final <- bind_rows(final, temp)
  }
  
}

final <- final %>%
  mutate(spp_id = ifelse(spp_id == "yellow","yellow poplar",spp_id))

ggplot(data = final)+
  geom_density(aes(x = p9, group = model_id, color = model_id, fill = model_id),
               alpha = 0.5)+
  theme_classic()+
  facet_wrap(facets = vars(spp_id), scales = "free_y")+
  xlab("ozone coefficient")+
  geom_vline(xintercept = 0)+
  ggtitle("")

multi_params <- final %>%
  select(-c(.chain,.iteration,.draw)) %>%
  pivot_longer(p2:global_tree_effect, names_to = "param_name", 
               values_to = "param_value")

unpack_data1 <- multi_params %>%
  filter(param_name %in% c("p4","p10")) %>%
  mutate(param_name = ifelse(param_name == "p10","average","deviation"))

ggplot(data = unpack_data1)+
  geom_density(aes(x = param_value, group = param_name, color = param_name, fill = param_name),
               alpha = 0.5)+
  theme_classic()+
  facet_wrap(facets = vars(spp_id), scales = "free_y")+
  xlab("Noxi interval average and deviation coefficients")+
  geom_vline(xintercept = 0)+
  ggtitle("unpacking plot effect")

unpack_data2 <- multi_params %>%
  filter(param_name %in% c("p5","p11")) %>%
  mutate(param_name = ifelse(param_name == "p11","average","deviation"))

ggplot(data = unpack_data2)+
  geom_density(aes(x = param_value, group = param_name, color = param_name, fill = param_name),
               alpha = 0.5)+
  theme_classic()+
  facet_wrap(facets = vars(spp_id), scales = "free_y")+
  xlab("Nred interval average and deviation coefficients")+
  geom_vline(xintercept = 0)+
  ggtitle("unpacking plot effect")

# for historic deviation experiment

# list files
out <- list.files("./experiments/historic_deviation",pattern = "mcmc.parquet",
                  full.names = TRUE)

for(i in 1:length(out)){
  
  spp_name = str_split(out[i], pattern = "-")[[1]][2]
  model_name = str_split(out[i], pattern = "/")[[1]][3]
  temp <- read_parquet(file = out[i]) %>%
    mutate(spp_id = spp_name,
           model_id = model_name)
  
  if(i == 1){
    final <- temp
  } else {
    final <- bind_rows(final, temp)
  }
  
}

final <- final %>%
  mutate(spp_id = ifelse(spp_id == "yellow","yellow poplar",spp_id))

multi_params <- final %>%
  select(-c(.chain,.iteration,.draw)) %>%
  pivot_longer(p2:global_tree_effect, names_to = "param_name", 
               values_to = "param_value")

unpack_data1 <- multi_params %>%
  filter(param_name %in% c("p4","p6")) %>%
  mutate(param_name = ifelse(param_name == "p4","baseline","difference"))

ggplot(data = unpack_data1)+
  geom_density(aes(x = param_value, group = param_name, color = param_name, fill = param_name),
               alpha = 0.5)+
  theme_classic()+
  facet_wrap(facets = vars(spp_id), scales = "free_y")+
  xlab("Noxi baseline and interval difference coefficients")+
  geom_vline(xintercept = 0)+
  ggtitle("historic deviation")

ggplot(data = final)+
  geom_point(aes(x = p4, y = p6))+
  facet_wrap(facets = vars(spp_id), scales = "free")+
  theme_bw()

unpack_data2 <- multi_params %>%
  filter(param_name %in% c("p5","p7")) %>%
  mutate(param_name = ifelse(param_name == "p5","baseline","difference"))

ggplot(data = unpack_data2)+
  geom_density(aes(x = param_value, group = param_name, color = param_name, fill = param_name),
               alpha = 0.5)+
  theme_classic()+
  facet_wrap(facets = vars(spp_id), scales = "free_y")+
  xlab("Nred baseline and interval difference coefficients")+
  geom_vline(xintercept = 0)+
  ggtitle("historic deviation")

ggplot(data = final)+
  geom_point(aes(x = p5, y = p7))+
  facet_wrap(facets = vars(spp_id), scales = "free")+
  theme_bw()

ggplot(data = final)+
  geom_point(aes(x = p6, y = p7))+
  facet_wrap(facets = vars(spp_id), scales = "free")+
  theme_bw()

# what do differences actually look like?

ggplot(data = df)+
  geom_density(aes(x = Dep_Noxidiff))+
  facet_wrap(facets = vars(common_name), scales = "free_y")+
  geom_vline(xintercept = 0, color = "red")+
  theme_bw()

ggplot(data = df)+
  geom_density(aes(x = Dep_Nreddiff))+
  facet_wrap(facets = vars(common_name), scales = "free_y")+
  geom_vline(xintercept = 0, color = "red")+
  theme_bw()

# visualizing coefs for baseline proportion experiment 01MAY25 ----

# list files
out <- list.files("./experiments/baseline_proportion",pattern = "mcmc.parquet",
                  full.names = TRUE)

for(i in 1:length(out)){
  
  spp_name = str_split(out[i], pattern = "-")[[1]][2]
  model_name = str_split(out[i], pattern = "/")[[1]][3]
  temp <- read_parquet(file = out[i]) %>%
    mutate(spp_id = spp_name,
           model_id = model_name)
  
  if(i == 1){
    final <- temp
  } else {
    final <- bind_rows(final, temp)
  }
  
}

final <- final %>%
  mutate(spp_id = ifelse(spp_id == "yellow","yellow poplar",spp_id))

ggplot(data = final)+
  geom_density(aes(x = p5, group = model_id, color = model_id, fill = model_id),
               alpha = 0.5)+
  theme_classic()+
  facet_wrap(facets = vars(spp_id), scales = "free")+
  geom_vline(xintercept = 0)+
  ggtitle("")

ggplot(data = final)+
  geom_density(aes(x = p9, group = model_id, color = model_id, fill = model_id),
               alpha = 0.5)+
  theme_classic()+
  facet_wrap(facets = vars(spp_id), scales = "free")+
  geom_vline(xintercept = 0)+
  ggtitle("")

long_short_term_df <- final %>%
  select(spp_id, p5, p9) %>%
  pivot_longer(p5:p9, names_to = "param_name", values_to = "param_value") %>%
  group_by(spp_id, param_name) %>%
  summarize(mean = mean(param_value, na.rm = TRUE),
            q97.5 = quantile(param_value, probs = c(0.975)),
            q2.5 = quantile(param_value, probs = c(0.025)))
p5 <- long_short_term_df %>%
  filter(param_name == "p5") %>%
  rename(p5_mean = mean,
         `p5_q97.5` = `q97.5`,
         `p5_q2.5` = `q2.5`) %>%
  select(-param_name)
p9 <- long_short_term_df %>%
  filter(param_name == "p9") %>%
  rename(p9_mean = mean,
         `p9_q97.5` = `q97.5`,
         `p9_q2.5` = `q2.5`) %>%
  select(-param_name)
long_short_term_df2 <- full_join(p5, p9, by = c("spp_id"))

ggplot(data = long_short_term_df2)+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  geom_point(aes(x = p9_mean, y = p5_mean, col = spp_id))+
  geom_segment(aes(x = p9_q2.5, y = p5_mean, xend = p9_q97.5, yend = p5_mean,
                   col = spp_id))+
  geom_segment(aes(y = p5_q2.5, x = p9_mean, yend = p5_q97.5, xend = p9_mean,
                   col = spp_id))+
  xlim(c(-2.5,2.5))+
  ylim(c(-0.2,0.2))+
  theme_bw()+
  ylab("Coefficient on N dep deviation")+
  xlab("Coefficient on proportion of baseline N dep")

# $$ figure: delta growth at different values of baseline N and delta N

og_df <- read_csv("./data/McDonnell_etal_InPrep_TreeData_2024_10_11.csv", show_col_types = FALSE) %>%
  dplyr::filter(!common_name %in% c("Douglas-fir","western hemlock")) 
#filter(common_name %in% c("eastern cottonwood"))

focal_df <- og_df %>%
  select(common_name, plot_ID, tree_ID, interval_no, Dep_N, Dep_N15, Dep_Noxi15, Dep_Nred15, Dep_Noxi, Dep_Nred, 
         Dep_S, Dep_S15, MAT, MAP, date_m2, date_m1, AG_carbon_pYear, AG_carbon_m1, 
         AG_carbon_m2, subp_BA_GT_m1, live_m2, Ozone_avg) 

baseline_vars <- focal_df %>%
  select(plot_ID, date_m1, date_m2, Dep_N, Dep_Noxi, Dep_Nred, Dep_S, MAT, MAP, Ozone_avg) %>%
  distinct(.) %>%
  group_by(plot_ID) %>%
  summarize(Dep_Nbaseline = mean(Dep_N, na.rm = TRUE),
            Dep_Noxibaseline = mean(Dep_Noxi, na.rm = TRUE),
            Dep_Nredbaseline = mean(Dep_Nred, na.rm = TRUE),
            Dep_Sbaseline = mean(Dep_S, na.rm = TRUE),
            MAT_baseline = mean(MAT, na.rm = TRUE),
            MAP_baseline = mean(MAP, na.rm = TRUE),
            Ozone_avg_baseline = mean(Ozone_avg, na.rm = TRUE)) %>%
  ungroup()

baseline_vars2 <- focal_df %>%
  select(plot_ID, date_m1, date_m2, Dep_N15, Dep_N, Dep_Noxi,
         Dep_Noxi15, Dep_Nred, Dep_Nred15, Dep_S, Dep_S15) %>%
  distinct(.) %>%
  group_by(plot_ID) %>%
  filter(date_m1 == min(date_m1, na.rm = TRUE)) %>%
  arrange(plot_ID) %>%
  mutate(dt = as.numeric(date_m2 - date_m1)/365) %>%
  mutate(total_years = 15 + dt) %>%
  mutate(Dep_Nhistoric = (Dep_N15*total_years - Dep_N*dt) / 15,
         Dep_Noxihistoric = (Dep_Noxi15*total_years - Dep_Noxi*dt) / 15,
         Dep_Nredhistoric = (Dep_Nred15*total_years - Dep_Nred*dt) / 15,
         Dep_Shistoric = (Dep_S15*total_years - Dep_S*dt) / 15) %>%
  mutate(Dep_NpropBaseline = Dep_N / Dep_Nhistoric,
         Dep_NoxipropBaseline = Dep_Noxi / Dep_Noxihistoric,
         Dep_NredpropBaseline = Dep_Nred / Dep_Nredhistoric) %>%
  select(plot_ID, Dep_Nhistoric, Dep_Noxihistoric, Dep_Nredhistoric, Dep_Shistoric,
         Dep_NpropBaseline, Dep_NoxipropBaseline, Dep_NredpropBaseline)

focal_df2 <- left_join(focal_df, baseline_vars, by = "plot_ID") %>%
  left_join(baseline_vars2, by = "plot_ID") %>%
  group_by(plot_ID) %>%
  mutate(Dep_Ndelta = Dep_N - Dep_Nbaseline,
         Dep_Noxidelta = Dep_Noxi - Dep_Noxibaseline,
         Dep_Nreddelta = Dep_Nred - Dep_Nredbaseline,
         Dep_Sdelta = Dep_S - Dep_Sbaseline,
         MAP_delta_dm = (MAP - MAP_baseline) * 0.01,
         MAT_delta = MAT - MAT_baseline,
         Ozone_avg_delta = Ozone_avg - Ozone_avg_baseline,
         Dep_Ndiff = Dep_N - Dep_Nhistoric,
         Dep_Noxidiff = Dep_Noxi - Dep_Noxihistoric,
         Dep_Nreddiff = Dep_Nred - Dep_Nredhistoric,
         Dep_Sdiff = Dep_S - Dep_Shistoric,
         MAP_baseline_dm = MAP_baseline * 0.01) %>%
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
  ungroup() %>%
  dplyr::filter(num_intervals >= 2) %>%
  mutate(diff = Dep_N - Dep_Nhistoric) %>%
  mutate(common_name = ifelse(common_name == "yellow-poplar","yellow poplar",common_name))

range(df$Dep_Nhistoric)
# [1]  2.94734 33.29947
range(df$diff)
# [1] -13.44028  23.39295

mean_data <- df %>%
  select(common_name, AG_carbon_m1, subp_BA_GT_m1) %>%
  group_by(common_name) %>%
  summarize_all(mean, na.rm = TRUE)

mean_params <- final %>%
  select(spp_id, p2:p9, global_tree_effect) %>%
  group_by(spp_id) %>%
  summarize_all(mean, na.rm = TRUE)

baseline_Ndep_values = seq(2.0, 35.0, by = 0.1)
delta_Ndep_values = seq(-14, 24, by = 0.1)

pred_df <- data.frame(baseline_Ndep = rep(rep(baseline_Ndep_values, times = length(delta_Ndep_values)), times = 8),
                      delta_Ndep = rep(delta_Ndep_values, each = length(baseline_Ndep_values), times = 8),
                      spp_id = rep(unique(df$common_name), each = length(baseline_Ndep_values)*length(delta_Ndep_values))) %>%
  mutate(prop_baseline = (delta_Ndep + baseline_Ndep) / baseline_Ndep) %>%
  filter(prop_baseline >= 0.4 & prop_baseline <= 3) # this is to match the observed range
range(pred_df$prop_baseline)

species <- unique(pred_df$spp_id)

pred <- NULL

for(i in 1:length(species)){
  
  dat <- mean_data %>%
    filter(common_name == species[i])
  params <- mean_params %>%
    filter(spp_id == species[i])
  
  prop_baseline <- pred_df %>%
    filter(spp_id == species[i])
  message(species[i])
  message(length(prop_baseline$prop_baseline))
  
  pred_no_change <- ((params$global_tree_effect + 1*params$p9) * dat$AG_carbon_m1 ^ params$p2) * exp(-dat$subp_BA_GT_m1*params$p3)
  pred_change <- ((params$global_tree_effect + prop_baseline$prop_baseline*params$p9) * dat$AG_carbon_m1 ^ params$p2) * exp(-dat$subp_BA_GT_m1*params$p3)
  
  prop_baseline$pred <- (pred_change - pred_no_change)
  
  if(i == 1){
    pred <- prop_baseline
  } else {
    pred <- bind_rows(pred, prop_baseline)
  }
}
unique(pred$spp_id)
unique(pred_df$spp_id)
final_pred <- left_join(pred_df, pred)

plots <- NULL

for(i in 1:length(species)){
  
  plot_dat <- final_pred %>%
    filter(spp_id == species[i])
  
  plots[[i]] <- ggplot(data = plot_dat)+
              geom_contour_filled(aes(x = baseline_Ndep, y = delta_Ndep, z = pred))+
              ggtitle(species[i])+
    theme_bw()
}

money_plot <- ggarrange(plots[[1]], plots[[2]], plots[[3]], plots[[4]], plots[[5]], plots[[6]], plots[[7]], plots[[8]], 
                      nrow = 3, ncol = 3,
                      labels = c("A","B","C","D","E","F","G","H")
)
money_plot
ggsave(plot = money_plot, filename = "./visualizations/deltaGrowth_deltaN_baselineN.tif",
       device = "tiff", height = 10, width = 14, units = "in",bg = "white")

# same thing but for deviation

range(df$Dep_Nbaseline)
# [1]  2.781566 32.111455
range(df$Dep_Ndelta)
# [1] -10.32194  11.39654

mean_Ndep_values = seq(2.0, 33.0, by = 0.1)
dev_Ndep_values = seq(-11, 12, by = 0.1)

pred_df2 <- data.frame(mean_Ndep = rep(rep(mean_Ndep_values, times = length(dev_Ndep_values)), times = 8),
                      dev_Ndep = rep(dev_Ndep_values, each = length(mean_Ndep_values), times = 8),
                      spp_id = rep(unique(df$common_name), each = length(mean_Ndep_values)*length(dev_Ndep_values))) 

species2 <- unique(pred_df2$spp_id)

pred2 <- NULL

for(i in 1:length(species2)){
  
  dat <- mean_data %>%
    filter(common_name == species2[i])
  params <- mean_params %>%
    filter(spp_id == species2[i])
  
  N_dev <- pred_df2 %>%
    filter(spp_id == species2[i])
  message(species[i])
  message(length(N_dev$dev_Ndep))
  
  pred_no_change <- ((params$global_tree_effect + 0*params$p5) * dat$AG_carbon_m1 ^ params$p2) * exp(-dat$subp_BA_GT_m1*params$p3)
  pred_change <- ((params$global_tree_effect + N_dev$dev_Ndep*params$p5) * dat$AG_carbon_m1 ^ params$p2) * exp(-dat$subp_BA_GT_m1*params$p3)
  
  N_dev$pred <- (pred_change - pred_no_change)
  
  if(i == 1){
    pred2 <- N_dev
  } else {
    pred2 <- bind_rows(pred2, N_dev)
  }
}
unique(pred2$spp_id)
unique(pred_df2$spp_id)
final_pred2 <- left_join(pred_df2, pred2)

plots2 <- NULL

for(i in 1:length(species2)){
  
  plot_dat <- final_pred2 %>%
    filter(spp_id == species2[i])
  
  plots2[[i]] <- ggplot(data = plot_dat)+
    geom_contour_filled(aes(x = mean_Ndep, y = dev_Ndep, z = pred))+
    ggtitle(species[i])+
    theme_bw()
}

money_plot2 <- ggarrange(plots2[[1]], plots2[[2]], plots2[[3]], plots2[[4]], plots2[[5]], plots2[[6]], plots2[[7]], plots2[[8]], 
                        nrow = 3, ncol = 3,
                        labels = c("A","B","C","D","E","F","G","H")
)
money_plot2
ggsave(plot = money_plot2, filename = "./visualizations/deltaGrowth_devN_meanN.tif",
       device = "tiff", height = 10, width = 14, units = "in",bg = "white")

# visualizing coefs for baseline proportion experiment 14MAY25 ----

# list files
out <- list.files("./experiments/short-term_long-term",pattern = "mcmc.parquet",
                  full.names = TRUE)

for(i in 1:length(out)){
  
  spp_name = str_split(out[i], pattern = "-")[[1]][6]
  model_name = str_split(out[i], pattern = "/")[[1]][3]
  temp <- read_parquet(file = out[i]) %>%
    mutate(spp_id = spp_name,
           model_id = model_name)
  
  if(i == 1){
    final <- temp
  } else {
    final <- bind_rows(final, temp)
  }
  
}

final <- final %>%
  mutate(spp_id = ifelse(spp_id == "yellow","yellow poplar",spp_id))

ggplot(data = final)+
  geom_density(aes(x = p5, group = model_id, color = model_id, fill = model_id),
               alpha = 0.5)+
  theme_classic()+
  facet_wrap(facets = vars(spp_id), scales = "free")+
  geom_vline(xintercept = 0)+
  ggtitle("")

ggplot(data = final)+
  geom_density(aes(x = p9, group = model_id, color = model_id, fill = model_id),
               alpha = 0.5)+
  theme_classic()+
  facet_wrap(facets = vars(spp_id), scales = "free")+
  geom_vline(xintercept = 0)+
  ggtitle("")

ggplot(data = final)+
  geom_density(aes(x = p10, group = model_id, color = model_id, fill = model_id),
               alpha = 0.5)+
  theme_classic()+
  facet_wrap(facets = vars(spp_id), scales = "free")+
  geom_vline(xintercept = 0)+
  ggtitle("")

long_short_term_df <- final %>%
  select(spp_id, p5, p9) %>%
  pivot_longer(p5:p9, names_to = "param_name", values_to = "param_value") %>%
  group_by(spp_id, param_name) %>%
  summarize(mean = mean(param_value, na.rm = TRUE),
            q97.5 = quantile(param_value, probs = c(0.975)),
            q2.5 = quantile(param_value, probs = c(0.025)))
p5 <- long_short_term_df %>%
  filter(param_name == "p5") %>%
  rename(p5_mean = mean,
         `p5_q97.5` = `q97.5`,
         `p5_q2.5` = `q2.5`) %>%
  select(-param_name)
p9 <- long_short_term_df %>%
  filter(param_name == "p9") %>%
  rename(p9_mean = mean,
         `p9_q97.5` = `q97.5`,
         `p9_q2.5` = `q2.5`) %>%
  select(-param_name)
long_short_term_df2 <- full_join(p5, p9, by = c("spp_id")) %>%
  mutate(spp_id = ifelse(spp_id == "yellow","yellow poplar",spp_id))

og_df <- read_csv("./data/McDonnell_etal_InPrep_TreeData_2024_10_11.csv", show_col_types = FALSE) %>%
  dplyr::filter(!common_name %in% c("Douglas-fir","western hemlock")) 
#filter(common_name %in% c("eastern cottonwood"))

focal_df <- og_df %>%
  select(common_name, plot_ID, tree_ID, interval_no, Dep_N, Dep_N15, Dep_Noxi15, Dep_Nred15, Dep_Noxi, Dep_Nred, 
         Dep_S, Dep_S15, MAT, MAP, date_m2, date_m1, AG_carbon_pYear, AG_carbon_m1, 
         AG_carbon_m2, subp_BA_GT_m1, live_m2, Ozone_avg) 

baseline_vars <- focal_df %>%
  select(plot_ID, date_m1, date_m2, Dep_N, Dep_Noxi, Dep_Nred, Dep_S, MAT, MAP, Ozone_avg) %>%
  distinct(.) %>%
  group_by(plot_ID) %>%
  summarize(Dep_Nbaseline = mean(Dep_N, na.rm = TRUE),
            Dep_Noxibaseline = mean(Dep_Noxi, na.rm = TRUE),
            Dep_Nredbaseline = mean(Dep_Nred, na.rm = TRUE),
            Dep_Sbaseline = mean(Dep_S, na.rm = TRUE),
            MAT_baseline = mean(MAT, na.rm = TRUE),
            MAP_baseline = mean(MAP, na.rm = TRUE),
            Ozone_avg_baseline = mean(Ozone_avg, na.rm = TRUE)) %>%
  ungroup()

baseline_vars2 <- focal_df %>%
  select(plot_ID, date_m1, date_m2, Dep_N15, Dep_N, Dep_Noxi,
         Dep_Noxi15, Dep_Nred, Dep_Nred15, Dep_S, Dep_S15) %>%
  distinct(.) %>%
  group_by(plot_ID) %>%
  filter(date_m1 == min(date_m1, na.rm = TRUE)) %>%
  arrange(plot_ID) %>%
  mutate(dt = as.numeric(date_m2 - date_m1)/365) %>%
  mutate(total_years = 15 + dt) %>%
  mutate(Dep_Nhistoric = (Dep_N15*total_years - Dep_N*dt) / 15,
         Dep_Noxihistoric = (Dep_Noxi15*total_years - Dep_Noxi*dt) / 15,
         Dep_Nredhistoric = (Dep_Nred15*total_years - Dep_Nred*dt) / 15,
         Dep_Shistoric = (Dep_S15*total_years - Dep_S*dt) / 15) %>%
  mutate(Dep_NpropBaseline = Dep_N / Dep_Nhistoric,
         Dep_NoxipropBaseline = Dep_Noxi / Dep_Noxihistoric,
         Dep_NredpropBaseline = Dep_Nred / Dep_Nredhistoric) %>%
  select(plot_ID, Dep_Nhistoric, Dep_Noxihistoric, Dep_Nredhistoric, Dep_Shistoric,
         Dep_NpropBaseline, Dep_NoxipropBaseline, Dep_NredpropBaseline)

focal_df2 <- left_join(focal_df, baseline_vars, by = "plot_ID") %>%
  left_join(baseline_vars2, by = "plot_ID") %>%
  group_by(plot_ID) %>%
  mutate(Dep_Ndelta = Dep_N - Dep_Nbaseline,
         Dep_Noxidelta = Dep_Noxi - Dep_Noxibaseline,
         Dep_Nreddelta = Dep_Nred - Dep_Nredbaseline,
         Dep_Sdelta = Dep_S - Dep_Sbaseline,
         MAP_delta_dm = (MAP - MAP_baseline) * 0.01,
         MAT_delta = MAT - MAT_baseline,
         Ozone_avg_delta = Ozone_avg - Ozone_avg_baseline,
         Dep_Ndiff = Dep_N - Dep_Nhistoric,
         Dep_Noxidiff = Dep_Noxi - Dep_Noxihistoric,
         Dep_Nreddiff = Dep_Nred - Dep_Nredhistoric,
         Dep_Sdiff = Dep_S - Dep_Shistoric,
         MAP_baseline_dm = MAP_baseline * 0.01,
         Dep_N_LTchange = Dep_Nbaseline - Dep_Nhistoric) %>%
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
  ungroup() %>%
  dplyr::filter(num_intervals >= 2)

mean_data <- df %>%
  select(common_name, AG_carbon_m1, subp_BA_GT_m1, Dep_Nhistoric) %>%
  group_by(common_name) %>%
  summarize_all(mean, na.rm = TRUE) %>%
  mutate(common_name = ifelse(common_name == "yellow-poplar","yellow poplar",common_name))

mean_params <- final %>%
  select(spp_id, p2:p10, global_tree_effect) %>%
  group_by(spp_id) %>%
  summarize_all(mean, na.rm = TRUE) %>%
  mutate(spp_id = ifelse(spp_id == "yellow","yellow poplar",spp_id))

species <- unique(long_short_term_df2$spp_id)

pred_df <- matrix(nrow = 8, ncol = 7, data = NA)

for(i in 1:length(species)){
  
  dat <- mean_data %>%
    filter(common_name == species[i])
  params <- mean_params %>%
    filter(spp_id == species[i])
  
  delta_params <- long_short_term_df2 %>%
    filter(spp_id == species[i])
  message(species[i])
  
  for(j in 1:6){
    
    if(j %in% c(1:3)){
      pred_no_change <- ((params$global_tree_effect + 0*delta_params[1,j+1] + dat$Dep_Nhistoric*params$p10) * dat$AG_carbon_m1 ^ params$p2) * exp(-dat$subp_BA_GT_m1*params$p3)
      pred_change <- ((params$global_tree_effect + -1*delta_params[1,j+1] + dat$Dep_Nhistoric*params$p10) * dat$AG_carbon_m1 ^ params$p2) * exp(-dat$subp_BA_GT_m1*params$p3)
      pred <- unlist((pred_change - pred_no_change))
      pred_df[i,j+1] <- pred
    }
    else {
      pred_no_change <- ((params$global_tree_effect + 0*delta_params[1,j+1] + dat$Dep_Nhistoric*params$p10) * dat$AG_carbon_m1 ^ params$p2) * exp(-dat$subp_BA_GT_m1*params$p3)
      pred_change <- ((params$global_tree_effect + -1*delta_params[1,j+1] + dat$Dep_Nhistoric*params$p10) * dat$AG_carbon_m1 ^ params$p2) * exp(-dat$subp_BA_GT_m1*params$p3)
      pred <- unlist((pred_change - pred_no_change))
      pred_df[i,j+1] <- pred
    }
  }

}

pred_df[,1] <- species
colnames(pred_df) <- colnames(long_short_term_df2)
pred_df <- data.frame(pred_df) %>%
  mutate_at(vars(p5_mean:`p9_q2.5`), as.numeric)

ggplot(data = pred_df)+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  geom_point(aes(x = p9_mean, y = p5_mean, col = spp_id))+
  geom_segment(aes(x = p9_q2.5, y = p5_mean, xend = p9_q97.5, yend = p5_mean,
                   col = spp_id))+
  geom_segment(aes(y = p5_q2.5, x = p9_mean, yend = p5_q97.5, xend = p9_mean,
                   col = spp_id))+
  xlim(c(-4,4))+
  ylim(c(-4,4))+
  theme_bw()+
  ylab("Growth change given 1 kg ha^-1 yr^-1 N 'short-term' decrease")+
  xlab("Growth change given 1 kg ha^-1 yr^-1 N 'long-term' decrease")

ggplot(data = pred_df)+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  geom_point(aes(x = p9_mean, y = p5_mean, col = spp_id))+
  geom_segment(aes(x = p9_q2.5, y = p5_mean, xend = p9_q97.5, yend = p5_mean,
                   col = spp_id))+
  geom_segment(aes(y = p5_q2.5, x = p9_mean, yend = p5_q97.5, xend = p9_mean,
                   col = spp_id))+
  xlim(c(-0.6,0.6))+
  ylim(c(-0.6,0.6))+
  theme_bw()+
  ylab("Growth change given 1 kg ha^-1 yr^-1 N 'short-term' decrease")+
  xlab("Growth change given 1 kg ha^-1 yr^-1 N 'long-term' decrease")

# $$ figure: delta growth at different values of baseline N and long term delta N


range(df$Dep_Nhistoric)
# [1]  2.94734 33.29947
range(df$Dep_N_LTchange)
# [1] -11.58009  15.08363

ggplot(data = df)+
  geom_density(aes(x = Dep_N_LTchange), alpha = 0.5)+
  theme_classic()+
  facet_wrap(facets = vars(common_name), scales = "free")+
  geom_vline(xintercept = 0)+
  ggtitle("")
  
baseline_Ndep_values = seq(2.0, 35.0, by = 0.1)
delta_Ndep_values = seq(-12, 16, by = 0.1)

pred_df <- data.frame(baseline_Ndep = rep(rep(baseline_Ndep_values, times = length(delta_Ndep_values)), times = 8),
                      delta_Ndep = rep(delta_Ndep_values, each = length(baseline_Ndep_values), times = 8),
                      spp_id = rep(unique(df$common_name), each = length(baseline_Ndep_values)*length(delta_Ndep_values))) %>%
  mutate(spp_id = ifelse(spp_id == "yellow-poplar","yellow poplar",spp_id)) %>%
  filter(baseline_Ndep + delta_Ndep >= 0)

species <- unique(pred_df$spp_id)

pred <- NULL

for(i in 1:length(species)){
  
  dat <- mean_data %>%
    filter(common_name == species[i])
  params <- mean_params %>%
    filter(spp_id == species[i])
  
  prop_baseline <- pred_df %>%
    filter(spp_id == species[i])
  message(species[i])

  pred_no_change <- ((params$global_tree_effect + 0*params$p9 + prop_baseline$baseline_Ndep*params$p10) * dat$AG_carbon_m1 ^ params$p2) * exp(-dat$subp_BA_GT_m1*params$p3)
  pred_change <- ((params$global_tree_effect + prop_baseline$delta_Ndep*params$p9 + prop_baseline$baseline_Ndep*params$p10) * dat$AG_carbon_m1 ^ params$p2) * exp(-dat$subp_BA_GT_m1*params$p3)
  
  #prop_baseline$pred <- (pred_change - pred_no_change)
  prop_baseline$pred <- pred_change
  
  if(i == 1){
    pred <- prop_baseline
  } else {
    pred <- bind_rows(pred, prop_baseline)
  }
}
unique(pred$spp_id)
unique(pred_df$spp_id)
final_pred <- left_join(pred_df, pred)

plots <- NULL

for(i in 1:length(species)){
  
  plot_dat <- final_pred %>%
    filter(spp_id == species[i])
  
  plots[[i]] <- ggplot(data = plot_dat)+
    geom_tile(aes(x = baseline_Ndep, y = delta_Ndep, fill = pred))+
    ggtitle(species[i])+
    scale_fill_viridis_c()+
    theme_bw()
}

money_plot <- ggarrange(plots[[1]], plots[[2]], plots[[3]], plots[[4]], plots[[5]], plots[[6]], plots[[7]], plots[[8]], 
                        nrow = 3, ncol = 3,
                        labels = c("A","B","C","D","E","F","G","H")
)
money_plot
ggsave(plot = money_plot, filename = "./visualizations/deltaGrowth_deltaN_baselineN_tile.tif",
       device = "tiff", height = 10, width = 14, units = "in",bg = "white")

# same thing but for deviation

range(df$Dep_Nhistoric)
# [1]  2.94734 33.29947
range(df$Dep_Ndelta)
# [1] -10.32194  11.39654

mean_Ndep_values = seq(2.0, 34.0, by = 0.1)
dev_Ndep_values = seq(-11, 12, by = 0.1)

pred_df2 <- data.frame(mean_Ndep = rep(rep(mean_Ndep_values, times = length(dev_Ndep_values)), times = 8),
                       dev_Ndep = rep(dev_Ndep_values, each = length(mean_Ndep_values), times = 8),
                       spp_id = rep(unique(df$common_name), each = length(mean_Ndep_values)*length(dev_Ndep_values))) %>%
  mutate(spp_id = ifelse(spp_id == "yellow-poplar","yellow poplar",spp_id)) %>%
  filter(mean_Ndep + dev_Ndep >= 0)

species2 <- unique(pred_df2$spp_id)

pred2 <- NULL

for(i in 1:length(species2)){
  
  dat <- mean_data %>%
    filter(common_name == species2[i])
  params <- mean_params %>%
    filter(spp_id == species2[i])
  
  N_dev <- pred_df2 %>%
    filter(spp_id == species2[i])
  message(species[i])

  pred_no_change <- ((params$global_tree_effect + 0*params$p5 + N_dev$mean_Ndep*params$p10) * dat$AG_carbon_m1 ^ params$p2) * exp(-dat$subp_BA_GT_m1*params$p3)
  pred_change <- ((params$global_tree_effect + N_dev$dev_Ndep*params$p5 + N_dev$mean_Ndep*params$p10) * dat$AG_carbon_m1 ^ params$p2) * exp(-dat$subp_BA_GT_m1*params$p3)
  
  #N_dev$pred <- (pred_change - pred_no_change)
  N_dev$pred <- pred_change
  
  if(i == 1){
    pred2 <- N_dev
  } else {
    pred2 <- bind_rows(pred2, N_dev)
  }
}
unique(pred2$spp_id)
unique(pred_df2$spp_id)
final_pred2 <- left_join(pred_df2, pred2)

plots2 <- NULL

for(i in 1:length(species2)){
  
  plot_dat <- final_pred2 %>%
    filter(spp_id == species2[i])
  
  plots2[[i]] <- ggplot(data = plot_dat)+
    geom_tile(aes(x = mean_Ndep, y = dev_Ndep, fill = pred))+
    ggtitle(species[i])+
    scale_fill_viridis_c()+
    theme_bw()
}

money_plot2 <- ggarrange(plots2[[1]], plots2[[2]], plots2[[3]], plots2[[4]], plots2[[5]], plots2[[6]], plots2[[7]], plots2[[8]], 
                         nrow = 3, ncol = 3,
                         labels = c("A","B","C","D","E","F","G","H")
)
money_plot2
ggsave(plot = money_plot2, filename = "./visualizations/deltaGrowth_devN_meanN.tif",
       device = "tiff", height = 10, width = 14, units = "in",bg = "white")

