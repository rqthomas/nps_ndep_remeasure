# model conceptual figure
# Author: Mary Lofton
# Date: 04APR25

# Purpose: develop a conceptual figure to explain how the random effects work
# in the N dep models

# load packages
library(tidyverse)
library(lubridate)
library(arrow)

# list files
out <- list.files("./experiments/random_effects",pattern = ".parquet",
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

plot_a <- mean(final$`plot_effect[3]`, na.rm = TRUE)
plot_b <- mean(final$`plot_effect[13]`, na.rm = TRUE)
rect_color = "black"

p1 <- ggplot(data = final)+
  geom_density(aes(x = global_tree_effect))+
  theme_classic()+
  xlab("Global species effect")+
  geom_vline(xintercept = mean(final$global_tree_effect, na.rm = TRUE), col = "black")+
  geom_vline(xintercept = plot_a, col = "red", linetype = "dashed")+
  geom_vline(xintercept = plot_b, col = "blue", linetype = "dashed")+
  theme(panel.background=element_rect(colour=rect_color, linewidth = 2),
        axis.line.x.bottom=element_line(color=rect_color),
        axis.line.y.left=element_line(color=rect_color))
ggsave(plot = p1, filename = "./visualizations/global_tree_effect.png",
       device = "png", height = 3, width = 3, units = "in")

tree1 <- mean(final$`tree_effect[3]`, na.rm = TRUE)
tree2 <- mean(final$`tree_effect[4]`, na.rm = TRUE)

rect_color = "red"

p2 <- ggplot(data = final)+
  geom_density(aes(x = `plot_effect[3]`))+
  theme_classic()+
  xlab("Plot A effect")+
  geom_vline(xintercept = plot_a, col = "red")+
  geom_vline(xintercept = tree1, col = "pink", linetype = "dashed")+
  geom_vline(xintercept = tree2, col = "darkred", linetype = "dashed")+
  theme(panel.background=element_rect(colour=rect_color, linewidth = 2),
        axis.line.x.bottom=element_line(color=rect_color),
        axis.line.y.left=element_line(color=rect_color))
ggsave(plot = p2, filename = "./visualizations/plot_a_effect.png",
       device = "png", height = 3, width = 3, units = "in")

tree3 <- mean(final$`tree_effect[39]`, na.rm = TRUE)
tree4 <- mean(final$`tree_effect[40]`, na.rm = TRUE)
tree5 <- mean(final$`tree_effect[41]`, na.rm = TRUE)

rect_color = "blue"

p3 <- ggplot(data = final)+
  geom_density(aes(x = `plot_effect[13]`))+
  theme_classic()+
  xlab("Plot B effect")+
  geom_vline(xintercept = plot_b, col = "blue")+
  geom_vline(xintercept = tree3, col = "skyblue", linetype = "dashed")+
  geom_vline(xintercept = tree4, col = "cornflowerblue", linetype = "dashed")+
  geom_vline(xintercept = tree5, col = "darkblue", linetype = "dashed")+
  theme(panel.background=element_rect(colour=rect_color, linewidth = 2),
        axis.line.x.bottom=element_line(color=rect_color),
        axis.line.y.left=element_line(color=rect_color))
ggsave(plot = p3, filename = "./visualizations/plot_b_effect.png",
       device = "png", height = 3, width = 3, units = "in")

rect_color = "pink"

p4 <- ggplot(data = final)+
  geom_density(aes(x = `tree_effect[3]`))+
  theme_classic()+
  xlab("Tree 1 effect")+
  geom_vline(xintercept = tree1, col = "pink")+
  theme(panel.background=element_rect(colour=rect_color, linewidth = 2),
        axis.line.x.bottom=element_line(color=rect_color),
        axis.line.y.left=element_line(color=rect_color))
ggsave(plot = p4, filename = "./visualizations/tree1_effect.png",
       device = "png", height = 3, width = 3, units = "in")

rect_color = "darkred"

p5 <- ggplot(data = final)+
  geom_density(aes(x = `tree_effect[4]`))+
  theme_classic()+
  xlab("Tree 2 effect")+
  geom_vline(xintercept = tree2, col = "darkred")+
  theme(panel.background=element_rect(colour=rect_color, linewidth = 2),
        axis.line.x.bottom=element_line(color=rect_color),
        axis.line.y.left=element_line(color=rect_color))
ggsave(plot = p5, filename = "./visualizations/tree2_effect.png",
       device = "png", height = 3, width = 3, units = "in")

rect_color = "skyblue"

p6 <- ggplot(data = final)+
  geom_density(aes(x = `tree_effect[39]`))+
  theme_classic()+
  xlab("Tree 3 effect")+
  geom_vline(xintercept = tree3, col = "skyblue")+
  theme(panel.background=element_rect(colour=rect_color, linewidth = 2),
        axis.line.x.bottom=element_line(color=rect_color),
        axis.line.y.left=element_line(color=rect_color))
ggsave(plot = p6, filename = "./visualizations/tree3_effect.png",
       device = "png", height = 3, width = 3, units = "in")

rect_color = "cornflowerblue"

p7 <- ggplot(data = final)+
  geom_density(aes(x = `tree_effect[40]`))+
  theme_classic()+
  xlab("Tree 4 effect")+
  geom_vline(xintercept = tree4, col = "cornflowerblue")+
  theme(panel.background=element_rect(colour=rect_color, linewidth = 2),
        axis.line.x.bottom=element_line(color=rect_color),
        axis.line.y.left=element_line(color=rect_color))
ggsave(plot = p7, filename = "./visualizations/tree4_effect.png",
       device = "png", height = 3, width = 3, units = "in")

rect_color = "darkblue"

p8 <- ggplot(data = final)+
  geom_density(aes(x = `tree_effect[41]`))+
  theme_classic()+
  xlab("Tree 5 effect")+
  geom_vline(xintercept = tree5, col = "darkblue")+
  theme(panel.background=element_rect(colour=rect_color, linewidth = 2),
        axis.line.x.bottom=element_line(color=rect_color),
        axis.line.y.left=element_line(color=rect_color))
ggsave(plot = p8, filename = "./visualizations/tree5_effect.png",
       device = "png", height = 3, width = 3, units = "in")
