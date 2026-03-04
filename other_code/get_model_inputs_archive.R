# Wrangle raw data to obtain model input df
# Author: Mary Lofton
# Date: 03APR25

get_model_inputs <- function(data){
  ## READ IN AND WRANGLE DATA
  og_df <- read_csv("./data/McDonnell_etal_InPrep_TreeData_2024_10_11.csv", show_col_types = FALSE) %>%
    dplyr::filter(!common_name %in% c("Douglas-fir","western hemlock")) 
  #filter(common_name %in% c("eastern cottonwood"))
  
  focal_df <- og_df %>%
    select(common_name, plot_ID, tree_ID, interval_no, Dep_N, Dep_N15, Dep_Noxi15, Dep_Nred15, Dep_Noxi, Dep_Nred, 
           Dep_S, Dep_S15, MAT, MAP, date_m2, date_m1, AG_carbon_pYear, AG_carbon_m1, 
           AG_carbon_m2, subp_BA_GT_m1, live_m2, Ozone_avg, species) 
  
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
  
  return(df)
}