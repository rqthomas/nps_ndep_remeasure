# McDonnell data figure
# Author: Mary Lofton
# Date: 20OCT25
# Last updated: 22DEC25

# Purpose: data figure to illustrate:

#'A. map of mean N deposition during measurement period (including 15 yr prior to 1st measure)
#'B. timeseries of short-term differences in N deposition by plot and ecoregion

# load packages
library(tidyverse)
library(lubridate)
library(sf)
library(viridis)
library(ggthemes)
library(ggpubr)
library(rnaturalearth)
library(grid)

og_df <- read_csv("./data/McDonnell_etal_InPrep_TreeData_2024_10_11.csv", show_col_types = FALSE) %>%
  dplyr::filter(!common_name %in% c("Douglas-fir","western hemlock")) 
#filter(common_name %in% c("eastern cottonwood"))

df <- read_csv("./data/processed_data.csv")

ann_ndep <- read_csv("./data/Fig3_annual_plot_Ndep.csv")

spp_df <- og_df %>%
  select(species, common_name) %>%
  distinct(.)

#'A. map of mean N deposition during measurement period (including 15 yr prior to 1st measure)

usa_coordinates <- map_data("state")
map_data_ante <- left_join(df, spp_df, by = "common_name") %>%
  select(plot_ID, lat, lon, Dep_N15, Dep_N, date_m1, date_m2) %>%
  distinct(.) %>%
  group_by(plot_ID, lat, lon) %>%
  dplyr::filter(date_m1 == min(date_m1, na.rm = TRUE)) %>%
  mutate(dt = as.numeric(date_m2 - date_m1)/365) %>%
  mutate(total_years = 15 + dt) %>%
  mutate(Dep_Nante = (Dep_N15*total_years - Dep_N*dt)) %>%
  select(plot_ID, lat, lon, Dep_Nante) %>%
  distinct(.)

map_data <- left_join(df, spp_df, by = "common_name") %>%
  select(plot_ID, lat, lon, Dep_N, date_m1, date_m2) %>%
  distinct(.) %>%
  mutate(dt = as.numeric(date_m2 - date_m1)/365,
         Dep_Ncum = Dep_N*dt) %>%
  group_by(plot_ID, lat, lon) %>%
  summarize(Dep_Ncumall = sum(Dep_Ncum),
            dt_all = sum(dt)) %>%
  left_join(map_data_ante, by = c("plot_ID","lat","lon")) %>%
  mutate(Dep_Nmean = (Dep_Nante + Dep_Ncumall) / (dt_all + 15)) %>%
  select(plot_ID, lat, lon, Dep_Nmean)

plot <- df %>%
  filter(level1_ecoregion == "EASTERN TEMPERATE FORESTS" & num_intervals == 4) %>%
  select(plot_ID, Dep_Nhistoric, Dep_N, date_m1, date_m2, interval_no, lat, lon) %>%
  distinct() %>%
  filter(Dep_Nhistoric > 20 & interval_no == 1 & date_m1 == "2000-05-10") %>%
  pull(plot_ID)

plot_data <- df %>%
  select(plot_ID, Dep_Nhistoric, Dep_N, date_m1, date_m2, interval_no, lat, lon) %>%
  filter(plot_ID == plot) %>%
  distinct()

segment_df <- plot_data[1,] %>%
  mutate(lat_start = lat - 5, lon_start = lon + 5, lon = lon + 0.35, lat = lat - 0.35)
  
fig3_a <- ggplot() +
  geom_map(
    data = usa_coordinates, map = usa_coordinates,
    aes(long, lat, map_id = region),
    color = "black", fill = "white")+
  geom_point(
    data = map_data,
    aes(lon, lat, color = Dep_Nmean),
    shape = 16, alpha = 1, size = 1
  ) +
  geom_point(
    data = plot_data,
    aes(lon, lat),
    shape = 0, size = 2, stroke = 2
  ) +
  geom_segment(
    data = segment_df,
    aes(x = lon_start, y = lat_start, xend = lon, yend = lat),
    arrow = arrow(length = unit(0.2,"cm"))
  ) +
  annotate("text", x = segment_df$lon_start, y = segment_df$lat_start - 1.2, label = "plot location \nused in (b)") +
  scale_color_viridis(option = "H")+
  xlab("")+
  ylab("")+
  theme_classic()+
  labs(color = expression(paste("mean N deposition (kg N ", ha^-1," ",y^-1,")")))+
  guides(color = guide_colorbar(title.position = "bottom", title.hjust = 0.0))+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        legend.background = element_rect(fill = "transparent"))+
  theme(legend.position = "inside", legend.position.inside =  c(0.2, 0.15),
        legend.direction = "horizontal", plot.margin = unit(c(0.1,0,0.1,0), "pt"))+
  ylim(c(min(usa_coordinates$lat),max(usa_coordinates$lat)))


fig3_a

# B. conceptual figure of how antecedent and recent N dep are calculated for example plot

plot_data2 <- plot_data %>%
  mutate(datetime_ante_finish = date_m1) %>%
  pivot_longer(date_m1:date_m2, names_to = "measurement_time", values_to = "datetime") %>%
  pivot_longer(Dep_Nhistoric:Dep_N, names_to = "dep_type", values_to = "Dep_N") %>%
  mutate(interval_no_by_type = paste(dep_type, interval_no, sep = "_"),
         datetime_ante_start = datetime %m-% years(15),
         datetime = ifelse((measurement_time == "date_m1" & dep_type == "Dep_Nhistoric"), datetime_ante_start,
                           ifelse((measurement_time == "date_m2" & dep_type == "Dep_Nhistoric"), datetime_ante_finish, datetime)),
         datetime = as.Date(datetime, origin = "1970-01-01"))

interval_starts <- plot_data2 %>%
  filter(measurement_time == "date_m1" & dep_type == "Dep_N") %>%
  select(datetime, measurement_time) %>%
  mutate(datetime = datetime + 35) %>%
  add_column(interval_no = c(1:4),
             legend_var = "a_interval")
interval_ends <- plot_data2 %>%
  filter(measurement_time == "date_m2" & dep_type == "Dep_N") %>%
  select(datetime, measurement_time) %>%
  mutate(datetime = datetime - 35) %>%
  add_column(interval_no = c(1:4),
             legend_var = "a_interval") 
intervals <- bind_rows(interval_starts, interval_ends)

diffs <- plot_data %>%
  select(date_m1, Dep_Nhistoric, Dep_N, interval_no) %>%
  add_column(linewidth_var = "Dep_Ndiff",
             direction = "down") %>%
  mutate(date_m1 = date_m1 + 35)

# Custom drawing function to reverse the vertical path
GeomArrow <- ggproto(
  NULL, GeomSegment,
  required_aes = c("x", "y", "xend", "yend", "direction"),
  draw_layer = function (self, data, params, layout, coord) 
  {
    swap <- data$direction == "down"
    v1 <- data$x
    v2 <- data$xend
    data$x[swap] <- v2[swap]
    data$xend[swap] <- v1[swap]
    GeomSegment$draw_layer(data, params, layout, coord)
  }
)

draw_key_vpath_reversed <- function(data, params, size) {
  # Original vpath goes from y=0.1 to y=0.9
  # This version goes from y=0.9 to y=0.1
  segmentsGrob(
    x0 = 0.5, y0 = 0.9, 
    x1 = 0.5, y1 = 0.1, # Reversed y coordinates
    gp = grid::gpar(
      col  = alpha(data$colour %||% data$fill %||% "black", data$alpha),
      fill = alpha(params$arrow.fill %||% data$colour %||% data$fill %||% 
                     "black", data$alpha),
      lwd = (data$size %||% 0.5) * .pt,
      lty = data$linetype %||% 1, lineend = "butt"
    ),
    arrow = params$arrow
  )
}

my.cols <- c("#DEEBF7", "#9ECAE1", "#4292C6", "#084594")

fig3_b <- ggplot()+
  geom_point(data = ann_ndep, aes(x = datetime, y = ndep, shape = "annual N dep."))+
  geom_vline(data = intervals, aes(xintercept = datetime, linetype = legend_var, alpha = legend_var), color = "black")+
  geom_vline(data = intervals, aes(xintercept = datetime, linetype = legend_var, alpha = legend_var, color = as.factor(interval_no)), show.legend = FALSE)+
  geom_line(data = plot_data2, aes(x = datetime, y = Dep_N, group = interval_no_by_type, alpha = dep_type, linetype = dep_type), linewidth = 1, color = "black")+
  geom_line(data = plot_data2, aes(x = datetime, y = Dep_N, group = interval_no_by_type, color = as.factor(interval_no), alpha = dep_type, linetype = dep_type), linewidth = 1, show.legend = FALSE)+
  geom_point(data = plot_data2, aes(x = datetime, y = Dep_N, group = interval_no_by_type, color = as.factor(interval_no)))+
  #geom_segment(data = diffs, aes(x = date_m1, y = Dep_Nhistoric, xend = date_m1, yend = Dep_N, linewidth = linewidth_var), arrow = arrow(length = unit(0.2,"cm")), color = "gray", key_glyph = draw_key_vpath_reversed)+
  #geom_segment(data = diffs, aes(x = date_m1, y = Dep_Nhistoric, xend = date_m1, yend = Dep_N, color = as.factor(interval_no), linewidth = linewidth_var), arrow = arrow(length = unit(0.2,"cm")), show.legend = FALSE, key_glyph = draw_key_vpath_reversed)+
  stat_identity(
    geom = GeomArrow,
    data = diffs,
    aes(x = date_m1, y = Dep_Nhistoric, xend = date_m1, yend = Dep_N, linewidth = linewidth_var, direction = direction),
    color = "black",
    arrow = arrow(length = unit(0.3, "cm"), type = "closed"),
    key_glyph = draw_key_vpath_reversed
  )+
  stat_identity(
    geom = GeomArrow,
    data = diffs,
    aes(x = date_m1, y = Dep_Nhistoric, xend = date_m1, yend = Dep_N, linewidth = linewidth_var, direction = direction, color = as.factor(interval_no)),
    arrow = arrow(length = unit(0.3, "cm"), type = "closed"),
    key_glyph = draw_key_vpath_reversed, show.legend = FALSE
  )+
  theme_classic()+
  scale_color_manual(values = my.cols, name = "Interval")+
  scale_alpha_manual(values = c(a_interval = 1, Dep_N = 0.5, Dep_Nhistoric = 1), name = NULL, labels = c("measurement interval \nstart/end date","current interval N dep.","antecedent N dep."))+
  scale_linetype_manual(values = c(a_interval = 2, Dep_N = 2, Dep_Nhistoric = 1), name = NULL, labels = c("measurement interval \nstart/end date","current interval N dep.","antecedent N dep."))+
  scale_linewidth_manual(values = c("Dep_Ndiff" = 0.5), name = NULL, labels = c("short-term change in N dep."))+
  scale_shape_manual(values = c("annual N dep." = 16), name = NULL)+
  xlab("Study year")+
  ylab(expression(paste("N deposition (kg N ", ha^-1," ",y^-1,")")))+
  theme(legend.spacing.y = unit(0.0, "cm"),
        legend.key.width = unit(3, "line"),
        legend.key.height = unit(1.5, "line"))+
  guides(
    color = guide_legend(order = 1,
                         override.aes = list(
                           size = c(5, 5, 5, 5),
                           shape = c(15, 15, 15, 15)
                         )),
    linetype = guide_legend(override.aes = list(
      linetype = c(2, 2, 1), # Example values
      alpha = c(1, 1, 1) # Example values
    ), reverse = TRUE),
    alpha = guide_legend(override.aes = list(
      linetype = c(1, 2, 2), # Example values
      alpha = c(1, 1, 1) # Example values
    ), reverse = TRUE),
    linewidth = guide_legend(order = 4)
  )

fig3_b

ggsave(plot = fig3_b, filename = "./visualizations/final_figures/Figure3b.tif",
       device = "tiff", height = 4, width = 7, units = "in", bg = "white")

# Assemble figure

p2 <- ggarrange(fig3_a, fig3_b,
                nrow = 2, ncol = 1,
                labels = c("(a)","(b)")
)
p2
ggsave(plot = p2, filename = "./visualizations/final_figures/Figure3.tif",
       device = "tiff", height = 8, width = 7.5, units = "in", bg = "white")

#### END OF CODE FOR FINAL FIGURE INCLUDED IN MANUSCRIPT


# Additional figure versions for CLAD presentation

fig3_b_1 <- ggplot()+
  geom_vline(data = interval_starts, aes(xintercept = datetime, linetype = measurement_time, alpha = measurement_time), color = "gray")+
  geom_line(data = subset(plot_data2, dep_type == "Dep_Nhistoric" & interval_no == 1), aes(x = datetime, y = Dep_N, group = interval_no_by_type, color = as.factor(interval_no), alpha = dep_type, linetype = dep_type), linewidth = 1)+
  geom_point(data = subset(plot_data2, dep_type == "Dep_Nhistoric" & interval_no == 1), aes(x = datetime, y = Dep_N, group = interval_no_by_type, color = as.factor(interval_no)))+
  #geom_segment(data = diffs, aes(x = date_m1, y = Dep_Nhistoric, xend = date_m1, yend = Dep_N, color = as.factor(interval_no), linewidth = linewidth_var), arrow = arrow(length = unit(0.2,"cm")))+
  theme_classic()+
  scale_color_discrete(name = "Measurement interval")+
  scale_alpha_manual(values = c(date_m1 = 1, Dep_Nhistoric = 1), name = NULL, labels = c("measurement \ninterval start date","antecedent N dep."))+
  scale_linetype_manual(values = c(date_m1 = 2, Dep_Nhistoric = 1), name = NULL, labels = c("measurement \ninterval start date","antecedent N dep."))+
  scale_linewidth_manual(values = c("Dep_Ndiff" = 0.5), name = NULL, labels = c("short-term change in N dep."))+
  xlab("")+
  ylab(expression(paste("N deposition (kg N ", ha^-1," ",y^-1,")")))+
  theme(legend.spacing.y = unit(0.0, "cm"),
        legend.key.width = unit(3, "line"),
        legend.key.height = unit(1.5, "line"))+
  guides(
    color = guide_legend(order = 1),
    linetype = guide_legend(override.aes = list(
      linetype = c(2, 1), # Example values
      alpha = c(1, 1) # Example values
    ), reverse = TRUE),
    alpha = guide_legend(override.aes = list(
      linetype = c(1, 2), # Example values
      alpha = c(1, 1) # Example values
    ), reverse = TRUE),
    linewidth = guide_legend(order = 4)
  )+
  xlim(c(min(plot_data2$datetime_ante_start), max(plot_data2$datetime)))+
  ylim(c(min(plot_data2$Dep_N), max(plot_data2$Dep_N)))

fig3_b_1

ggsave(plot = fig3_b_1, filename = "./visualizations/final_figures/Figure3b_1.tif",
       device = "tiff", height = 4, width = 6.5, units = "in", bg = "white")

fig3_b_2 <- ggplot()+
  geom_vline(data = interval_starts, aes(xintercept = datetime, linetype = measurement_time, alpha = measurement_time), color = "gray")+
  geom_line(data = subset(plot_data2, dep_type == "Dep_Nhistoric"), aes(x = datetime, y = Dep_N, group = interval_no_by_type, color = as.factor(interval_no), alpha = dep_type, linetype = dep_type), linewidth = 1)+
  geom_point(data = subset(plot_data2, dep_type == "Dep_Nhistoric"), aes(x = datetime, y = Dep_N, group = interval_no_by_type, color = as.factor(interval_no)))+
  #geom_segment(data = diffs, aes(x = date_m1, y = Dep_Nhistoric, xend = date_m1, yend = Dep_N, color = as.factor(interval_no), linewidth = linewidth_var), arrow = arrow(length = unit(0.2,"cm")))+
  theme_classic()+
  scale_color_discrete(name = "Measurement interval")+
  scale_alpha_manual(values = c(date_m1 = 1, Dep_Nhistoric = 1), name = NULL, labels = c("measurement \ninterval start date","antecedent N dep."))+
  scale_linetype_manual(values = c(date_m1 = 2, Dep_Nhistoric = 1), name = NULL, labels = c("measurement \ninterval start date","antecedent N dep."))+
  scale_linewidth_manual(values = c("Dep_Ndiff" = 0.5), name = NULL, labels = c("short-term change in N dep."))+
  xlab("Study year")+
  ylab(expression(paste("N deposition (kg N ", ha^-1," ",y^-1,")")))+
  theme(legend.spacing.y = unit(0.0, "cm"),
        legend.key.width = unit(3, "line"),
        legend.key.height = unit(1.5, "line"))+
  guides(
    color = guide_legend(order = 1),
    linetype = guide_legend(override.aes = list(
      linetype = c(2, 1), # Example values
      alpha = c(1, 1) # Example values
    ), reverse = TRUE),
    alpha = guide_legend(override.aes = list(
      linetype = c(1, 2), # Example values
      alpha = c(1, 1) # Example values
    ), reverse = TRUE),
    linewidth = guide_legend(order = 4)
  )+
  xlim(c(min(plot_data2$datetime_ante_start), max(plot_data2$datetime)))+
  ylim(c(min(plot_data2$Dep_N), max(plot_data2$Dep_N)))

fig3_b_2

ggsave(plot = fig3_b_2, filename = "./visualizations/final_figures/Figure3b_2.tif",
       device = "tiff", height = 4, width = 6.5, units = "in", bg = "white")
