# Job script to run model
# Author: Mary Lofton
# Date: 13JAN24
# Last updated: 22DEC25

# Purpose: example workflow script to run model for 100 individuals of one species

library(tidyverse)
library(rjags)
library(tidybayes)
library(bayesplot)

mem.maxVSize(vsize = Inf)

# eventually users will pull this from our data publication, but just pulling from the data folder for now
og_df <- read_csv("./data/McDonnell_etal_InPrep_TreeData_2024_10_11.csv", show_col_types = FALSE) %>%
  dplyr::filter(!(common_name %in% c("Douglas-fir","western hemlock"))) %>%
  dplyr::filter(!(AG_carbon_pYear < -2000))

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

# write processed data to file
write.csv(df, "./data/processed_data.csv", row.names = FALSE)

# optionally, start from here and read in file
# df <- read_csv("./data/processed_data_example.csv") 

total_species <- length(unique(df$common_name))

sim <- "ortho_log_t_interaction_example"
print(sim)

# model code

# set counter k = 1 because we are not looping, we are just running one species
k = 1

tree_species <- unique(df$common_name)[k]

print(tree_species)

focal_data <- df %>%
  dplyr::filter(common_name == tree_species) 

trees_index <- focal_data |> 
  distinct(tree_ID) |> 
  mutate(tree_index = 1:n())

plots <- focal_data |> 
  distinct(plot_ID) |> 
  mutate(plot_index = 1:n()) 

df1 <- focal_data |> 
  left_join(plots, by = join_by(plot_ID)) |> 
  left_join(trees_index, by = join_by(tree_ID)) |>
  mutate(dt = as.numeric(date_m2 - date_m1)/365) |>
  mutate(log_AG_carbon_pYear = (log(AG_carbon_m2) - log(AG_carbon_m1)) / dt,
         log_AG_carbon_m1 = log(AG_carbon_m1),
         log_subp_BA_GT_m1 = log(subp_BA_GT_m1 + 1)) |>
  slice(c(1:100)) # limiting to 100 trees here for example!

log_mean_annual_avg_growth <- c(df1$log_AG_carbon_pYear)
log_start_measures <- c(df1$log_AG_carbon_m1)
plot_index <- c(df1$plot_index)
tree_index <- c(df1$tree_index)
log_ba_gt <- c(df1$log_subp_BA_GT_m1)
sdep_historic <- c(df1$Dep_Shistoric)
sdep_diff <- c(df1$Dep_Sdiff)
mat_delta <- c(df1$MAT_delta)
map_delta_dm <- c(df1$MAP_delta_dm)
ndep_historic <- c(df1$Dep_Nhistoric_ortho)
ndep_diff <- c(df1$Dep_Ndiff_ortho)

n_measures <- nrow(df1)

ssData   <- list(log_tree_agb_obs = log_start_measures,
                 log_tree_growth_obs = log_mean_annual_avg_growth,
                 ntrees = length(unique(df1$tree_ID)), 
                 n_plots = length(unique(df1$plot_index)),
                 n_measures = n_measures,
                 log_ba_gt = log_ba_gt,
                 sdep_historic = sdep_historic,
                 sdep_diff = sdep_diff,
                 mat_delta = mat_delta,
                 map_delta_dm = map_delta_dm,
                 ndep_historic = ndep_historic,
                 ndep_diff = ndep_diff,
                 plot_index = plot_index,
                 tree_index = tree_index)

cat( "model {

  global_tree_effect ~ dnorm(0,0.001)
  tau_global ~ dgamma(0.01,0.01)
  tau_plot ~ dgamma(0.01,0.01)
  procErr ~ dgamma(0.01,0.01) # consider using a gamma here and for other tau parameters
  nu ~ dexp(1/30) T(2,) # truncated exponential
  p2 ~ dnorm(0,1e-12)
  p3 ~ dnorm(0,1e-12)
  p5 ~ dnorm(0,1e-12)
  p6 ~ dnorm(0,1e-12)
  p7 ~ dnorm(0,1e-12)
  p8 ~ dnorm(0,1e-12)
  p9 ~ dnorm(0,1e-12)
  p10 ~ dnorm(0,1e-12)
  p11 ~ dnorm(0,1e-12)
  p12 ~ dnorm(0,1e-12)
  
  for(p in 1:n_plots){
    plot_effect[p] ~ dnorm(global_tree_effect, tau_global)
  }

  for(n in 1:ntrees){
    tree_effect[n] ~ dnorm(plot_effect[plot_index[n]], tau_plot)
  }
  
  for(t in 1:n_measures){

     log_tree_growth_mean[t] <-  ((tree_effect[tree_index[t]] + ndep_diff[t]*p5 + sdep_diff[t]*p6 + mat_delta[t]*p7 + map_delta_dm[t]*p8 + ndep_historic[t]*p9 + sdep_historic[t]*p10 + ndep_historic[t]*ndep_diff[t]*p11 + ndep_diff[t]*ndep_diff[t]*p12) #if on log scale would interpret these in terms of multiplicative change
     + p2*log_tree_agb_obs[t]) 
         - p3*log_ba_gt[t] # how much do we want to penalize the log growth based on basal area

    lambda[t] ~ dgamma(nu/2, nu/2)

    log_tree_growth_obs[t] ~ dnorm(log_tree_growth_mean[t], procErr*lambda[t]) # multiplied b/c JAGS uses precision
  
  }

  }  ", fill=T,file=paste0("./experiments/",sim,"/",sim,"-growthModel-",tree_species,".txt"))

#black cherry, eastern cottonwood, sugar maple, yellow-poplar, quaking aspen
#ponderosa pine, paper birch, red spruce
init_values <- data.frame(species = c("black cherry","eastern cottonwood","sugar maple",
                                      "yellow-poplar","quaking aspen","ponderosa pine",
                                      "paper birch","red spruce"),
                          p5 = c(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
                          p5_delta = c(0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5),
                          p6 = c(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
                          p6_delta = c(0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5),
                          p7 = c(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
                          p7_delta = c(0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5),
                          p8 = c(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
                          p8_delta = c(0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5),
                          p9 = c(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
                          p9_delta = c(0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5),
                          p10 = c(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
                          p10_delta = c(0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5),
                          p11 = c(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
                          p11_delta = c(0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5),
                          p12 = c(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
                          p12_delta = c(0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5))

inits <- list(list("global_tree_effect" = 1,
                   "tau_global" = 1,
                   "tau_plot" = 1,
                   "procErr" = 0.001,
                   "nu" = 3,
                   "p2" = 0.6,
                   "p3" = 1,
                   "p5" = init_values[k,"p5"] + init_values[k,"p5_delta"],
                   "p6" = init_values[k,"p6"] + init_values[k,"p6_delta"],
                   "p7" = init_values[k,"p7"] + init_values[k,"p7_delta"],
                   "p8" = init_values[k,"p8"] + init_values[k,"p8_delta"],
                   "p9" = init_values[k,"p9"] + init_values[k,"p9_delta"],
                   "p10" = init_values[k,"p10"] + init_values[k,"p10_delta"],
                   "p11" = init_values[k,"p11"] + init_values[k,"p11_delta"],
                   "p12" = init_values[k,"p12"] + init_values[k,"p12_delta"]),
              list("global_tree_effect" = 1,
                   "tau_global" = 1,
                   "tau_plot" = 1,
                   "procErr" = 0.001,
                   "nu" = 10,
                   "p2" = 0.4,
                   "p3" = 1,
                   "p5" = init_values[k,"p5"],
                   "p6" = init_values[k,"p6"],
                   "p7" = init_values[k,"p7"],
                   "p8" = init_values[k,"p8"],
                   "p9" = init_values[k,"p9"],
                   "p10" = init_values[k,"p10"],
                   "p11" = init_values[k,"p11"],
                   "p12" = init_values[k,"p12"]),
              list("global_tree_effect" = 1,
                   "tau_global" = 1,
                   "tau_plot" = 1,
                   "procErr" = 0.001,
                   "nu" = 100,
                   "p2" = 0.4,
                   "p3" = 1,
                   "p5" = init_values[k,"p5"] - init_values[k,"p5_delta"],
                   "p6" = init_values[k,"p6"] - init_values[k,"p6_delta"],
                   "p7" = init_values[k,"p7"] - init_values[k,"p7_delta"],
                   "p8" = init_values[k,"p8"] - init_values[k,"p8_delta"],
                   "p9" = init_values[k,"p9"] - init_values[k,"p9_delta"],
                   "p10" = init_values[k,"p10"] - init_values[k,"p10_delta"],
                   "p11" = init_values[k,"p11"] - init_values[k,"p11_delta"],
                   "p12" = init_values[k,"p12"] - init_values[k,"p12_delta"]))

ssFit    <- jags.model(data=ssData, file=paste0("./experiments/",sim,"/",sim,"-growthModel-",tree_species,".txt"), n.chains = 3, inits = inits, n.adapt = 10000)
parNames <- c("tree_effect", 
              "p2", 
              "p3", 
              "p5",
              "p6",
              "p7",
              "p8",
              "p9",
              "p10",
              "p11",
              "p12",
              "procErr", 
              "global_tree_effect",
              "plot_effect",
              "nu")

if(tree_species == "sugar maple"){
  
  ssFit <- coda.samples(ssFit, variable.names = parNames, n.iter=20000, thin = 5)
  start.iter <- seq(from = 10005, to = 35000, by = 5)[2001]
  ssFit <- window(ssFit, start = start.iter)
  
} else {
  
  ssFit <- coda.samples(ssFit, variable.names = parNames, n.iter=30000, thin = 5)
  start.iter <- seq(from = 10005, to = 45000, by = 5)[4001]
  ssFit <- window(ssFit, start = start.iter)
  
}

# get ESS
ess_vars <- ssFit[, varnames(ssFit) %in% parNames]
ess <- effectiveSize(ess_vars)
write.csv(ess, paste0("./experiments/",sim,"/",sim, "-",tree_species,"-ess.csv"))

mcmc <- spread_draws(ssFit,
                     p2, 
                     p3, 
                     p5,
                     p6,
                     p7,
                     p8,
                     p9,
                     p10,
                     p11,
                     p12,
                     procErr, 
                     global_tree_effect,
                     nu) 
mcmc_treeEffect <- spread_draws(ssFit,
                                `tree_effect`[n])


arrow::write_parquet(mcmc, sink = paste0("./experiments/",sim,"/",sim, "-",tree_species, "-mcmc.parquet"))
arrow::write_parquet(mcmc_treeEffect, sink = paste0("./experiments/",sim,"/",sim, "-",tree_species, "-mcmc-treeEffect.parquet"))

p <- mcmc |> 
  rename(chain = .chain, iteration = .iteration, draw = .draw) |> 
  select(-draw) |> 
  pivot_longer(-c("chain", "iteration"), names_to = "par", values_to = "value") |> 
  ggplot(aes(x = iteration, y = value, color = factor(chain))) + 
  geom_line() + 
  facet_wrap(~par, scales = "free")

ggsave(filename = paste0("./experiments/",sim,"/",sim, "-",tree_species,"-mcmc.pdf"), plot = p, device = "pdf")

mcmc_burned <- mcmc |> 
  dplyr::filter(.iteration > 1000)

p <- mcmc_burned |> 
  rename(chain = .chain,
         iteration = .iteration) |> 
  select(chain:p12) |> 
  select(-.draw, -iteration) |> 
  mcmc_pairs()


ggsave(filename = paste0("./experiments/",sim,"/",sim, "-", tree_species,"-pairs.pdf"), plot = p, device = "pdf", height = 12, width = 12)

  


