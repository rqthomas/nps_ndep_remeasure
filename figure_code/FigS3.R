# Correlation between antecedent and short-term change in N dep
# Author: Mary Lofton
# Date: 04NOV25
# Last updated: 22DEC25

# Load packages
library(tidyverse)
library(lubridate)

# Read in data
dat <- read_csv("./data/processed_data.csv") 

# Data wrangling
mod_dat <- dat %>%
  select(Dep_Nhistoric, Dep_Ndiff) %>%
  distinct(.)
mod <- lm(mod_dat$Dep_Ndiff ~ mod_dat$Dep_Nhistoric)
sum <- summary(mod)
r2 <- round(sum$r.squared,2)
fitted <- data.frame(Dep_Nhistoric = mod_dat$Dep_Nhistoric, fitted = mod$fitted.values)
pcor <- round(cor(mod_dat$Dep_Ndiff, mod_dat$Dep_Nhistoric),2)

mod_dat2 <- dat %>%
  select(Dep_Nhistoric, Dep_Shistoric) %>%
  distinct(.)
mod2 <- lm(mod_dat2$Dep_Nhistoric ~ mod_dat2$Dep_Shistoric)
sum2 <- summary(mod2)
r2_2 <- round(sum2$r.squared,2)
fitted2 <- data.frame(Dep_Shistoric = mod_dat2$Dep_Shistoric, fitted = mod2$fitted.values)
pcor2 <- round(cor(mod_dat2$Dep_Nhistoric, mod_dat2$Dep_Shistoric),2)

mod_dat3 <- dat %>%
  select(Dep_Ndiff, Dep_Sdiff) %>%
  distinct(.)
mod3 <- lm(mod_dat3$Dep_Ndiff ~ mod_dat3$Dep_Sdiff)
sum3 <- summary(mod3)
r2_3 <- round(sum3$r.squared,2)
fitted3 <- data.frame(Dep_Sdiff = mod_dat3$Dep_Sdiff, fitted = mod3$fitted.values)
pcor3 <- round(cor(mod_dat3$Dep_Ndiff, mod_dat3$Dep_Sdiff),2)


# Plotting
p <- ggplot()+
  geom_point(data = mod_dat, aes(x = Dep_Nhistoric, y = Dep_Ndiff),alpha = 0.3)+
  #geom_line(data = fitted, aes(x = Dep_Nhistoric, y = fitted), color = "blue")+
  annotate("text", label = "r = -0.59", x = 30, y = 20, color = "black")+
  #annotate("text", label = expression(paste("p-value of slope term < 2.2 x ",10^-16)), x = 40, y = 20, color = "black")+
  theme_classic()+
  xlab(expression(paste("antecedent N deposition (kg N ", ha^-1," ",y^-1,")")))+
  ylab(expression(paste("short-term N deposition change (kg N ", ha^-1," ",y^-1,")")))
p

p2 <- ggplot()+
  geom_point(data = mod_dat2, aes(x = Dep_Shistoric, y = Dep_Nhistoric),alpha = 0.3)+
  #geom_line(data = fitted2, aes(x = Dep_Shistoric, y = fitted), color = "blue")+
  annotate("text", label = "r = 0.84", x = 35, y = 30, color = "black")+
  #annotate("text", label = expression(paste("p-value of slope term < 2.2 x ",10^-16)), x = 30, y = 55, color = "black")+
  theme_classic()+
  xlab(expression(paste("antecedent S deposition (kg S ", ha^-1," ",y^-1,")")))+
  ylab(expression(paste("antecedent N deposition (kg N ", ha^-1," ",y^-1,")")))
p2

p3 <- ggplot()+
  geom_point(data = mod_dat3, aes(x = Dep_Sdiff, y = Dep_Ndiff),alpha = 0.3)+
  #geom_line(data = fitted3, aes(x = Dep_Sdiff, y = fitted), color = "blue")+
  annotate("text", label = "r = 0.76", x = -18, y = 20, color = "black")+
  #annotate("text", label = expression(paste("p-value of slope term < 2.2 x ",10^-16)), x = -13, y = 15, color = "black")+
  theme_classic()+
  xlab(expression(paste("short-term S deposition change (kg S ", ha^-1," ",y^-1,")")))+
  ylab(expression(paste("short-term N deposition change (kg N ", ha^-1," ",y^-1,")")))
p3

p_final <- ggarrange(p, p2, p3, ncol = 2, nrow = 2, labels = c("a","b","c"))

ggsave(p_final, filename = "./visualizations/final_figures/FigureS3.png", dev = "png",
       height = 8, width = 8, units = "in", bg = "white")
