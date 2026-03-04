# Fake histogram plot

my_col <- c("white","gray")

# all fert
h1 <- data.frame(ante = rnorm(1000, mean = -0.5, sd = 0.1),
                 st = rnorm(1000, mean = -0.4, sd = 0.1)) %>%
  pivot_longer(ante:st, names_to = "dep_type", values_to = "growth_change")

p_h1 <- ggplot(data = h1)+
  geom_density_pattern(aes(x = growth_change, group = dep_type, fill = dep_type, pattern = dep_type), color = "black")+
  xlim(c(-1,1))+
  geom_vline(xintercept = 0)+
  scale_color_manual(values = my_col, labels = c("antecedent N dep.","short-term change in N dep."), name = "")+
  scale_fill_manual(values = my_col, labels = c("antecedent N dep.","short-term change in N dep."), name = "")+
  theme_bw()+
  xlab(expression(paste("change in growth (kg C ", y^-1," ",ind^-1,")")))+
  scale_pattern_manual(values = c("none","circle"), labels = c("antecedent N dep.","short-term change in N dep."), name = "")

ggsave(p_h1,filename = "./visualizations/other_visualizations/H1.png",
       device = "png", height = 2, width = 5, units = "in")

# all stress
h2 <- data.frame(ante = rnorm(1000, mean = 0.5, sd = 0.1),
                 st = rnorm(1000, mean = 0.4, sd = 0.1)) %>%
  pivot_longer(ante:st, names_to = "dep_type", values_to = "growth_change")

p_h2 <- ggplot(data = h2)+
  geom_density_pattern(aes(x = growth_change, group = dep_type, fill = dep_type, pattern = dep_type), color = "black")+
  xlim(c(-1,1))+
  geom_vline(xintercept = 0)+
  scale_color_manual(values = my_col, labels = c("antecedent N dep.","short-term change in N dep."), name = "")+
  scale_fill_manual(values = my_col, labels = c("antecedent N dep.","short-term change in N dep."), name = "")+
  theme_bw()+
  xlab(expression(paste("change in growth (kg C ", y^-1," ",ind^-1,")")))+
  scale_pattern_manual(values = c("none","circle"), labels = c("antecedent N dep.","short-term change in N dep."), name = "")

ggsave(p_h2,filename = "./visualizations/other_visualizations/H2.png",
       device = "png", height = 2, width = 5, units = "in")

# long-term stress, short term fertilization
h3 <- data.frame(ante = rnorm(1000, mean = 0.5, sd = 0.1),
                 st = rnorm(1000, mean = -0.5, sd = 0.1)) %>%
  pivot_longer(ante:st, names_to = "dep_type", values_to = "growth_change")

p_h3 <- ggplot(data = h3)+
  geom_density_pattern(aes(x = growth_change, group = dep_type, fill = dep_type, pattern = dep_type), color = "black")+
  xlim(c(-1,1))+
  geom_vline(xintercept = 0)+
  scale_color_manual(values = my_col, labels = c("antecedent N dep.","short-term change in N dep."), name = "")+
  scale_fill_manual(values = my_col, labels = c("antecedent N dep.","short-term change in N dep."), name = "")+
  theme_bw()+
  xlab(expression(paste("change in growth (kg C ", y^-1," ",ind^-1,")")))+
  scale_pattern_manual(values = c("none","circle"), labels = c("antecedent N dep.","short-term change in N dep."), name = "")

ggsave(p_h3,filename = "./visualizations/other_visualizations/H3.png",
       device = "png", height = 2, width = 5, units = "in")

# long-term fertilization, long-term stress
h4 <- data.frame(ante = rnorm(1000, mean = -0.5, sd = 0.1),
                 st = rnorm(1000, mean = 0.5, sd = 0.1)) %>%
  pivot_longer(ante:st, names_to = "dep_type", values_to = "growth_change")

p_h4 <- ggplot(data = h4)+
  geom_density_pattern(aes(x = growth_change, group = dep_type, fill = dep_type, pattern = dep_type), color = "black")+
  xlim(c(-1,1))+
  geom_vline(xintercept = 0)+
  scale_color_manual(values = my_col, labels = c("antecedent N dep.","short-term change in N dep."), name = "")+
  scale_fill_manual(values = my_col, labels = c("antecedent N dep.","short-term change in N dep."), name = "")+
  theme_bw()+
  xlab(expression(paste("change in growth (kg C ", y^-1," ",ind^-1,")")))+
  scale_pattern_manual(values = c("none","circle"), labels = c("antecedent N dep.","short-term change in N dep."), name = "")

ggsave(p_h4,filename = "./visualizations/other_visualizations/H4.png",
       device = "png", height = 2, width = 5, units = "in")
