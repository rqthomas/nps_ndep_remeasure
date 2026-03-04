# Script to generate processed data from raw (McDonnell et al) data
# Author: Mary Lofton
# Date: 18FEB26

# Purpose: workflow script to run ss model for one species

library(tidyverse)
library(rjags)
library(tidybayes)
library(bayesplot)
library(furrr)

mem.maxVSize(vsize = Inf)

og_df <- read_csv("./data/McDonnell_etal_InPrep_TreeData_2024_10_11.csv", show_col_types = FALSE) %>%
  dplyr::filter(!(common_name %in% c("Douglas-fir","western hemlock"))) %>%
  dplyr::filter(!(AG_carbon_pYear < -2000))
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
  #dplyr::filter(date_m1 == min(date_m1, na.rm = TRUE)) %>%
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
  select(plot_ID, date_m1, Dep_Nhistoric, Dep_Noxihistoric, Dep_Nredhistoric, Dep_Shistoric,
         Dep_NpropBaseline, Dep_NoxipropBaseline, Dep_NredpropBaseline)

focal_df2 <- left_join(focal_df, baseline_vars, by = "plot_ID") %>%
  left_join(baseline_vars2, by = c("plot_ID","date_m1")) %>%
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

# read in ecoregion info
ecoreg <- read_csv("./data/plot_ecoregions.csv")

focal_df3 <- left_join(focal_df2, ecoreg, by = c("plot_ID"))

reg <- unique(ecoreg$level1_ecoregion)[!is.na(unique(ecoreg$level1_ecoregion))]

final <- NULL

#### orthogonalization code
for(i in 1:length(reg)){
  
  print(reg[i])
  
  if(i == 3) next
  if(i == 9) next
  
  current_dat <- focal_df3 %>%
    dplyr::filter(level1_ecoregion == reg[i]) %>%
    select(level1_ecoregion, plot_ID, Dep_Shistoric, Dep_Sdiff, Dep_Nhistoric, Dep_Ndiff) %>%
    distinct(.) %>%
    dplyr::filter(complete.cases(.)) %>%
    mutate(ecoregion_ortho = reg[i])
  
  if(i == 2){
    current_dat <- focal_df3 %>%
      dplyr::filter(level1_ecoregion == reg[i] | level1_ecoregion == reg[3]) %>%
      select(level1_ecoregion, plot_ID, Dep_Shistoric, Dep_Sdiff, Dep_Nhistoric, Dep_Ndiff) %>%
      distinct(.) %>%
      dplyr::filter(complete.cases(.)) %>%
      mutate(ecoregion_ortho = paste(reg[i],reg[3],sep = " & "))
  }
  if(i == 8){
    current_dat <- focal_df3 %>%
      dplyr::filter(level1_ecoregion == reg[i] | level1_ecoregion == reg[9]) %>%
      select(level1_ecoregion, plot_ID, Dep_Shistoric, Dep_Sdiff, Dep_Nhistoric, Dep_Ndiff) %>%
      distinct(.) %>%
      dplyr::filter(complete.cases(.)) %>%
      mutate(ecoregion_ortho = paste(reg[i],reg[9],sep = " & "))
  }
  
  r_s <- current_dat %>%
    pull(Dep_Sdiff)
  
  a_s <- current_dat %>%
    pull(Dep_Shistoric)
  
  r_n <- current_dat %>%
    pull(Dep_Ndiff)
  
  a_n <- current_dat %>%
    pull(Dep_Nhistoric)
  
  n_s = length(r_s)
  n_n = length(r_n)
  
  S = cbind(1, r_s, a_s) # design matrix needs to include intercept [1, S]
  P_s = diag(n_s) - (S %*% solve(t(S) %*% S, t(S))) # projection onto span{1,S}
  rtilde_n = as.numeric(P_s %*% r_n)
  atilde_n = as.numeric(P_s %*% a_n)
  
  P_n = diag(n_n) - (rtilde_n %*% solve(t(rtilde_n) %*% rtilde_n, t(rtilde_n)))
  
  astar_n = as.numeric(P_n %*% atilde_n)
  
  c0 = (t(rtilde_n) %*% atilde_n) / (t(rtilde_n) %*% rtilde_n)
  c = c0[1,]
  
  current_dat$Dep_Nhistoric_ortho = scale(astar_n)
  current_dat$Dep_Ndiff_ortho = scale(rtilde_n)
  current_dat$c = c
  current_dat$Dep_Nhistoric_ortho_mean = mean(astar_n)
  current_dat$Dep_Nhistoric_ortho_sd = sd(astar_n)
  current_dat$Dep_Ndiff_ortho_mean = mean(rtilde_n)
  current_dat$Dep_Ndiff_ortho_sd = sd(rtilde_n)
  current_dat$Dep_Nhistoric_SO = atilde_n
  current_dat$Dep_Ndiff_SO = rtilde_n
  
  if(i == 1){
    final <- current_dat
  } else {
    final <- bind_rows(final, current_dat)
  }
  
}

focal_df4 <- left_join(focal_df3, final, by = c("plot_ID","level1_ecoregion","Dep_Shistoric","Dep_Sdiff",
                                                "Dep_Nhistoric","Dep_Ndiff"))

####

live_tree_ids <- focal_df4 |>
  summarise(count = n(),
            sum = sum(live_m2), .by = tree_ID) |>
  dplyr::filter(count >= sum) |>
  pull(tree_ID)

focal_df5 <- focal_df4 |>
  dplyr::filter(tree_ID %in% live_tree_ids) |>
  group_by(tree_ID) |>
  tidyr::fill(subp_BA_GT_m1, .direction = "down") |>
  ungroup()

df <- focal_df5 %>%
  dplyr::filter(complete.cases(.)) %>%
  group_by(tree_ID) %>%
  mutate(num_intervals = max(interval_no, na.rm = TRUE)) %>%
  ungroup() %>%
  dplyr::filter(num_intervals >= 1 & !(Dep_Nhistoric > 60))

# these are the data used for input into Bayesian hierarchical models;
# the file size is ~420 MB
write.csv(df, "./data/processed_data.csv", row.names = FALSE)