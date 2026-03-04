library(tidyverse)
library(rjags)
library(tidybayes)
library(bayesplot)

df <- read_csv("./data/McDonnell_etal_InPrep_TreeData_2024_10_11.csv", show_col_types = FALSE) %>%
  filter(common_name %in% c("eastern cottonwood"))

total_species <- length(unique(df$common_name))

sim <- "random_effects"

k = 1

run_model <- function(k, df, sim){
  
  tree_species <- unique(df$common_name)[k]
  
  print(tree_species)
  
  live_tree_ids <- df |> 
    dplyr::filter(common_name == tree_species) |> 
    summarise(count = n(),
              sum = sum(live_m2), .by = tree_ID) |> 
    dplyr::filter(count >= sum) |> 
    pull(tree_ID)
  
  focal_data <- df |> 
    dplyr::filter(tree_ID %in% live_tree_ids) |> 
    group_by(tree_ID) |> 
    tidyr::fill(subp_BA_GT_m1, .direction = "down") |> 
    dplyr::filter(!is.na(subp_BA_GT_m1) & !is.na(Dep_N)) |>
    ungroup() 
  
  baseline_vars <- focal_data %>%
    select(common_name, tree_ID, interval_no, Dep_N15, Dep_N, Dep_S15, Dep_S, MAT, MAP, date_m2, date_m1) %>%
    group_by(tree_ID) %>%
    summarize(Dep_Nbaseline = mean(Dep_N, na.rm = TRUE),
              Dep_Sbaseline = mean(Dep_S, na.rm = TRUE),
              MAT_baseline = mean(MAT, na.rm = TRUE),
              MAP_baseline = mean(MAP, na.rm = TRUE)) %>%
    ungroup() %>%
    select(tree_ID, Dep_Nbaseline, Dep_Sbaseline, MAT_baseline, MAP_baseline)
  
  focal_data2 <- left_join(focal_data, baseline_vars, by = "tree_ID") %>%
    mutate(Dep_Ndelta = Dep_N - Dep_Nbaseline,
           Dep_Sdelta = Dep_S - Dep_Sbaseline,
           MAP_delta_dm = (MAP - MAP_baseline) * 0.01,
           MAT_delta = MAT - MAT_baseline) %>%
    filter(!is.na(Dep_Ndelta) & !is.na(Dep_Sdelta) & !is.na(MAP_delta_dm) & !is.na(MAT_delta))
  
  trees_index <- focal_data2 |> 
    dplyr::filter(common_name == tree_species) |> 
    distinct(tree_ID) |> 
    dplyr::filter(tree_ID %in% live_tree_ids) |> 
    mutate(tree_index = 1:n())
  
  plots <- focal_data2 |> 
    dplyr::filter(common_name == tree_species) |> 
    dplyr::filter(tree_ID %in% live_tree_ids) |> 
    distinct(plot_ID) |> 
    mutate(plot_index = 1:n()) #USE THIS TO GET PLOT_INDEX BELOW
  
  df1 <- focal_data2 |> 
    mutate(dt = as.numeric(date_m2 - date_m1)/365) |> 
    select(AG_carbon_pYear, AG_carbon_m1, AG_carbon_m2, tree_ID, plot_ID, dt, Dep_Ndelta, subp_BA_GT_m1, MAP_delta_dm, Dep_Sdelta, MAT_delta) |>
    dplyr::filter(tree_ID %in% live_tree_ids) |>
    left_join(plots, by = join_by(plot_ID)) |> 
    left_join(trees_index, by = join_by(tree_ID)) 
  
  mean_annual_avg_growth <- df1$AG_carbon_pYear
  start_measures <- df1$AG_carbon_m1
  plot_index <- df1$plot_index
  tree_index <- df1$tree_index
  #dt <- df1$dt
  ba_gt <- df1$subp_BA_GT_m1
  ndep_delta <- df1$Dep_Ndelta
  sdep_delta <- df1$Dep_Sdelta
  mat_delta <- df1$MAT_delta
  map_delta_dm <- df1$MAP_delta_dm
  n_measures <- nrow(df1)
  
  # # code to check for plots and trees
  # nplots <- length(unique(plot_index))
  # ntrees_p <- df1 %>%
  #   filter(plot_index == 3)
  # length(unique(ntrees_p$tree_ID))
  
  # use plot = 3 and plot = 13
  # use tree index 39, 40, 41 and 3, 4

  ssData   <- list(tree_agb_obs = start_measures,
                   tree_growth_obs = mean_annual_avg_growth,
                   ntrees = length(unique(df1$tree_ID)), 
                   n_plots = length(unique(df1$plot_index)),
                   n_measures = n_measures,
                   ba_gt = ba_gt,
                   #dt = dt,
                   ndep_delta = ndep_delta,
                   sdep_delta = sdep_delta,
                   mat_delta = mat_delta,
                   map_delta_dm = map_delta_dm,
                   plot_index = plot_index,
                   tree_index = tree_index)
                   #max_ndep = max(c(ndep), na.rm = TRUE),
                   #min_ndep = min(c(ndep), na.rm = TRUE))
                   #temp = temp,
                   #max_temp = max(c(temp), na.rm = TRUE),
                   #min_temp = min(c(temp), na.rm = TRUE),
                   #precip = precip,
                   #max_precip = max(c(precip), na.rm = TRUE),
                   #min_precip = min(c(precip), na.rm = TRUE),
                   #sdep = sdep,
                   #max_sdep = max(c(sdep), na.rm = TRUE),
                   #min_sdep = min(c(sdep), na.rm = TRUE),
                   #obsErr = 1)
  
  cat( "model {

  global_tree_effect ~ dunif(-10,10)
  tau_global ~ dunif(0.0001,10)
  tau_plot ~ dunif(0.0001,10)
  procErr ~ dunif(0.0001,10)
  p2 ~ dunif(0,2) 
  p3 ~ dunif(0,100) #Sample in log space
  p5 ~ dnorm(0,0.001)
  p6 ~ dnorm(0,0.001)
  p7 ~ dnorm(0,0.001)
  p8 ~ dnorm(0,0.001)
  #lp4 ~ dunif(log(min_ndep), log(max_ndep*2))
  #p5inv2 ~ dexp(0.001)  
  #p6 ~ dunif(min_temp*0.5, max_temp*2)
  #p7 ~ dunif(0.3, 20)  
  #p8 ~ dunif(min_precip*0.5, max_precip*2)
  #p9 ~ dunif(0.3, 20)   
  #p10 ~ dunif(min_sdep*0.5, max_sdep*2)
  #p10 ~ dunif(min_sdep, max_sdep*2)
  #p11 ~ dunif(0.3, 10)    
  
  for(p in 1:n_plots){
    plot_effect[p] ~ dnorm(global_tree_effect, tau_global)
  }

  for(n in 1:ntrees){
    # tree_effect[n] ~ dnorm(global_tree_effect, tau_global)
    tree_effect[n] ~ dnorm(plot_effect[plot_index[n]], tau_plot)
  }
  
  for(t in 1:n_measures){

    tree_growth_mean[t] <-  ((tree_effect[tree_index[t]] + ndep_delta[t]*p5 + sdep_delta[t]*p6 + mat_delta[t]*p7 + map_delta_dm[t]*p8) 
    * tree_agb_obs[t] ^ p2) 
        * exp(-ba_gt[t]*p3)  

    tree_growth_obs[t] ~ dnorm(tree_growth_mean[t], procErr)

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
                            p8_delta = c(0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5))
  
  inits <- list(list("global_tree_effect" = 1,
                     "tau_global" = 1,
                     "tau_plot" = 1,
                     "procErr" = 0.001,
                     "p2" = 0.6,
                     "p3" = 1,
                     "p5" = init_values[k,"p5"] + init_values[k,"p5_delta"],
                     "p6" = init_values[k,"p6"] + init_values[k,"p6_delta"],
                     "p7" = init_values[k,"p7"] + init_values[k,"p7_delta"],
                     "p8" = init_values[k,"p8"] + init_values[k,"p8_delta"]),
                list("global_tree_effect" = 1,
                     "tau_global" = 1,
                     "tau_plot" = 1,
                     "procErr" = 0.001,
                     "p2" = 0.4,
                     "p3" = 1,
                     "p5" = init_values[k,"p5"],
                     "p6" = init_values[k,"p6"],
                     "p7" = init_values[k,"p7"],
                     "p8" = init_values[k,"p8"]),
                list("global_tree_effect" = 1,
                     "tau_global" = 1,
                     "tau_plot" = 1,
                     "procErr" = 0.001,
                     "p2" = 0.4,
                     "p3" = 1,
                     "p5" = init_values[k,"p5"] - init_values[k,"p5_delta"],
                     "p6" = init_values[k,"p6"] - init_values[k,"p6_delta"],
                     "p7" = init_values[k,"p7"] - init_values[k,"p7_delta"],
                     "p8" = init_values[k,"p8"] - init_values[k,"p8_delta"]))
  
  ssFit    <- jags.model(data=ssData, file=paste0("./experiments/",sim,"/",sim,"-growthModel-",tree_species,".txt"), n.chains = 3, inits = inits, n.adapt = 5000)
  parNames <- c("tree_effect", 
                "p2", 
                "p3", 
                "p5",
                "p6",
                "p7",
                "p8",
                "procErr", 
                "global_tree_effect",
                "plot_effect")
  
  ssFit <- coda.samples(ssFit, variable.names = parNames, n.iter=15000, thin = 5)
  mcmc <- spread_draws(ssFit,
                       p2, 
                       p3, 
                       p5,
                       p6,
                       p7,
                       p8,
                       procErr, 
                       `tree_effect[3]`,
                       `tree_effect[4]`,
                       `tree_effect[39]`,
                       `tree_effect[40]`,
                       `tree_effect[41]`,
                       global_tree_effect,
                       `plot_effect[3]`,
                       `plot_effect[13]`) 
  
  
  arrow::write_parquet(mcmc, sink = paste0("./experiments/",sim,"/",sim, "-",tree_species, "-mcmc.parquet"))
  
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
    select(chain:p8) |> 
    select(-.draw, -iteration) |> 
    mcmc_pairs()
  
  
  ggsave(filename = paste0("./experiments/",sim,"/",sim, "-", tree_species,"-pairs.pdf"), plot = p, device = "pdf", height = 12, width = 12)
  

}
