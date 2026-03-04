library(tidyverse)
library(rjags)
library(tidybayes)

df <- read_csv("../McDonnell_etal_InPrep_TreeData_2024_10_11.csv", show_col_types = FALSE)

k <- 1
sim <- "testS2"

tree_species <- unique(df$common_name)[k]

print(tree_species)

live_tree_ids <- df |> 
  filter(common_name == tree_species) |> 
  summarise(count = n(),
            sum = sum(live_m2), .by = tree_ID) |> 
  filter(count >= sum) |> 
  pull(tree_ID)

focal_data <- df |> 
  filter(tree_ID %in% live_tree_ids) |> 
  group_by(tree_ID) |> 
  tidyr::fill(subp_BA_GT_m1, .direction = "down") |> 
  ungroup()


trees_index <- focal_data |> 
  filter(common_name == tree_species) |> 
  distinct(tree_ID) |> 
  filter(tree_ID %in% live_tree_ids) |> 
  mutate(tree_index = 1:n())

plots <- focal_data |> 
  filter(common_name == tree_species) |> 
  filter(tree_ID %in% live_tree_ids) |> 
  distinct(plot_ID) |> 
  mutate(plot_index = 1:n())

tree_measures <- focal_data |> 
  filter(common_name == tree_species) |> 
  filter(tree_ID %in% live_tree_ids) |> 
  group_by(tree_ID) |> 
  summarise(count = n())

tree_info <- focal_data |> 
  filter(common_name == tree_species) |> 
  filter(tree_ID %in% live_tree_ids) |> 
  distinct(tree_ID, plot_ID) |> 
  left_join(plots, by = join_by(plot_ID)) |> 
  left_join(tree_measures, by = join_by(tree_ID)) |> 
  left_join(trees_index, by = join_by(tree_ID))

max_measures <- max(tree_info$count)


df1 <- focal_data |> 
  mutate(dt = as.numeric(date_m2 - date_m1)/365) |> 
  select(AG_carbon_m1, AG_carbon_m2, tree_ID, dt, Dep_N, subp_BA_GT_m1, MAT, MAP, Dep_S)


tree_matrix <- matrix(NA, nrow(tree_info), max_measures)
n_measures <- rep(NA, nrow(tree_info))
dt <- matrix(NA, nrow(tree_info), max_measures)
ba_gt <- matrix(NA, nrow(tree_info), max_measures)
ndep <- matrix(NA, nrow(tree_info), max_measures)
sdep <- matrix(NA, nrow(tree_info), max_measures)
temp <- matrix(NA, nrow(tree_info), max_measures)
precip <- matrix(NA, nrow(tree_info), max_measures)

for(i in 1:nrow(tree_info)){
  measures <- tree_info |> 
    filter(tree_ID == tree_info$tree_ID[i]) |> 
    pull(count)
  
  n_measures[i] <- measures
  
  curr_tree <- df1 |> 
    filter(tree_ID == tree_info$tree_ID[i])
  
  if(measures > 1){
    for(j in 2:(measures)){
      tree_matrix[i, j-1] <-  curr_tree$AG_carbon_m1[j]
      tree_matrix[i, j] <-  curr_tree$AG_carbon_m2[j]
      dt[i,j] <- curr_tree$dt[j]
      ba_gt[i, j] <- curr_tree$subp_BA_GT_m1[j] 
      ndep[i, j] <- curr_tree$Dep_N[j] 
      sdep[i, j] <- curr_tree$Dep_S[j]
      temp[i, j] <- curr_tree$MAT[j]
      precip[i, j] <- curr_tree$MAP[j]
    }
  }else{
    tree_matrix[i, 1] <-  curr_tree$AG_carbon_m1[1]
    tree_matrix[i, 2] <-  curr_tree$AG_carbon_m2[1]
    dt[i,2] <- curr_tree$dt[1]
    ba_gt[i, 2] <- curr_tree$subp_BA_GT_m1[1] 
    ndep[i, 2] <- curr_tree$Dep_N[1] 
    sdep[i, 2] <- curr_tree$Dep_S[1]
    temp[i, 2] <- curr_tree$MAT[1]
    precip[i, 2] <- curr_tree$MAP[1]
  }
}

precip <- precip / 100

global_tree_effect <- uniform(-10,10)
tau_global <- uniform(0.0001,10)
procErr <- uniform(0.0001,10)
p2 <- uniform(0,2) 
p3 <- uniform(0,100) #Sample in log space
p4 <- uniform(min_ndep*0.5, max_ndep*2)
p5 <- uniform(0.3, 20)  
#lp4 ~ dunif(log(min_ndep), log(max_ndep*2))
#p5inv2 ~ dexp(0.001)  
p6 <- uniform(min_temp*0.5, max_temp*2)
p7 <- uniform(0.3, 20)  
p8 <- uniform(min_precip*0.5, max_precip*2)
p9 <- uniform(0.3, 20)   
#p10 ~ dunif(min_sdep*0.5, max_sdep*2)
p10 <- uniform(min_sdep, max_sdep*2)
p11 <- uniform(0.3, 10)   


tree_agb_latent <- normal(, obsErr)

for(n in 1:ntrees){
  tree_agb_latent[n,1] ~ dnorm(tree_agb_obs[n, 1], obsErr)
}

for(n in 1:ntrees){
  
  tree_effect[n] ~ dnorm(global_tree_effect, tau_global)
  
  
  for(t in 2:n_measures[n]){
    tree_agb_mean[n, t] <-  tree_agb_latent[n, t-1] +  dt[n,t] * (tree_effect[n] * tree_agb_latent[n, t-1] ^ p2) 
    * exp(-ba_gt[n,t]*p3)  
    * exp(-0.5 * (log(ndep[n,t] / p4)/ p5) ^ 2)  
    #* exp(-p5inv2 * log(ndep[n,t] - lp4))  
    * exp(-0.5 * (log(temp[n,t] / p6)/ p7) ^ 2)
    #* exp(-0.5 * (log(sdep[n,t] / p10)/ p11) ^ 2)
    * 1/ (1 + ((sdep[n,t]-min_sdep)/p10)^p11)
    * exp(-0.5 * (log(precip[n,t] / p8)/ p9) ^ 2)
    
    tree_agb_latent[n,t] ~ dnorm(tree_agb_mean[n, t], procErr/dt[n,t])
    tree_agb_obs[n,t] ~ dnorm(tree_agb_latent[n,t], obsErr)
  }
}