# Make table of ESS for all species and parameters
# Author: Mary Lofton
# Date: 14NOV25

# list ess files in model output folder
out <- list.files("./experiments/ortho_log_t_interaction_adj_priors",pattern = "ess.csv",
                  full.names = TRUE)

# get species names
data = "./data/McDonnell_etal_InPrep_TreeData_2024_10_11.csv"

spp_df <- read_csv(data) %>%
  select(common_name, species) %>%
  distinct(.)

# read and collate files
for(i in 1:length(out)){
  
  spp_id = str_split(out[i], pattern = "-")[[1]][2]
  
  if(spp_id == "yellow"){spp_id = "yellow-poplar"}
  
  sci_name = spp_df %>%
    filter(common_name == spp_id) %>%
    pull(species)
  
  temp <- read_csv(out[i]) %>%
    rename(Parameter = `...1`,
           ESS = x) %>%
    mutate(ESS = round(ESS)) %>%
    pivot_wider(names_from = Parameter, values_from = ESS) %>%
    add_column(species = sci_name) %>%
    rename(p1 = global_tree_effect,
           p4 = p5,
           p5 = p6,
           p6 = p7,
           p7 = p8,
           p8 = p9,
           p9 = p10,
           p10 = p11,
           p11 = p12) %>%
    select(species, p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, procErr, nu)
  
  if(i == 1){
    final <- temp
  } else {
    final <- bind_rows(final, temp)
  }
  
}

write.csv(final, "./experiments/ortho_log_t_interaction_01JAN26/ess_all.csv", row.names = FALSE)
