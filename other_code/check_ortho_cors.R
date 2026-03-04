# Check correlations of orthogonalized variables
# Author: Mary Lofton
# Date: 12SEP25

# Purpose: check whether "historic" and "short-term" N deposition variables are
# still orthogonal if grouped by ecoregion

# read in processed data
processed_dat <- read_csv("./data/processed_data.csv")

# join dataframes
dat <- processed_dat %>%
  select(ecoregion_ortho, plot_ID, Dep_Nhistoric_ortho, Dep_Ndiff_ortho, Dep_Shistoric, Dep_Sdiff) %>%
  distinct(.)

# check correlation by interval only - should be 0!

# check correlation by ecoregion and interval
regions <- unique(dat$ecoregion_ortho)[!is.na(unique(dat$ecoregion_ortho))]

final_df <- data.frame(ecoregion = regions,
                       n = NA,
                       cor_astar_n_rtilde_n = NA,
                       cor_astar_n_r_s = NA,
                       cor_astar_n_a_s = NA,
                       cor_rtilde_n_r_s = NA,
                       cor_rtilde_n_a_s = NA)

for(i in 1:length(regions)){
  
  current_dat <- dat %>%
    filter(ecoregion_ortho == regions[i]) 
  
  current_cor_astar_n_rtilde_n <- cor(current_dat$Dep_Nhistoric_ortho, current_dat$Dep_Ndiff_ortho)
  current_cor_astar_n_r_s <- cor(current_dat$Dep_Nhistoric_ortho, current_dat$Dep_Sdiff)
  current_cor_astar_n_a_s <- cor(current_dat$Dep_Nhistoric_ortho, current_dat$Dep_Shistoric)
  current_cor_rtilde_n_r_s <- cor(current_dat$Dep_Ndiff_ortho, current_dat$Dep_Sdiff)
  current_cor_rtilde_n_a_s <- cor(current_dat$Dep_Ndiff_ortho, current_dat$Dep_Shistoric)
  
  final_df$n[i] <- length(current_dat$Dep_Ndiff_ortho)
  final_df$cor_astar_n_rtilde_n[i] <- current_cor_astar_n_rtilde_n
  final_df$cor_astar_n_r_s[i] <- current_cor_astar_n_r_s
  final_df$cor_astar_n_a_s[i] <- current_cor_astar_n_a_s
  final_df$cor_rtilde_n_r_s[i] <- current_cor_rtilde_n_r_s
  final_df$cor_rtilde_n_a_s[i] <- current_cor_rtilde_n_a_s
  
  # if(length(current_dat$Dep_Ndeviation) > 0){
  #   p <- ggplot(data = current_dat, aes(x = Dep_Ndeviation, y = Dep_N_wgt_avg_to_date_ortho))+
  #     geom_point()+
  #     ggtitle(paste(regions[i],intervals[j]))
  #   plot_name <- paste0("./visualizations/corr_plots/",regions[i],"_",intervals[j],".png")
  #   ggsave(p, filename = plot_name, device = "png")
  # }
  
}

write.csv(final_df, "./data/ortho_correlations_by_ecoregion.csv",row.names = FALSE)

slide_df <- final_df %>%
  filter(!n == 0)
# visualize correlations
og_df <- read_csv("./data/McDonnell_etal_InPrep_TreeData_2024_10_11.csv", show_col_types = FALSE) %>%
  dplyr::filter(!common_name %in% c("Douglas-fir","western hemlock")) 
#filter(common_name %in% c("eastern cottonwood"))
interval_hist <- og_df %>%
  filter(interval_no == 1)
hist(interval_hist$date_m1, breaks = "years", main = "Histogram of interval start dates",
     xlab = "Start date of first measurement interval for tree")

plot_ids <- unique(og_df$plot_ID)
plot1 <- og_df %>%
  filter(plot_ID == plot_ids[4])
