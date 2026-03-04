# Correlation between antecedent and short-term change in N dep
# Author: Mary Lofton
# Date: 04NOV25
# Last updated: 22DEC25

# Load packages
library(tidyverse)
library(lubridate)
library(ggpubr)

# Read in data
dat <- read_csv("./data/processed_data.csv") 

# Data wrangling
mod_dat <- dat %>%
  select(Dep_Nhistoric_SO, Dep_Nred) %>%
  distinct(.)
pcor <- round(cor(mod_dat$Dep_Nhistoric_SO, mod_dat$Dep_Nred),2)

mod_dat2 <- dat %>%
  select(Dep_Ndiff_SO, Dep_Nred) %>%
  distinct(.)
pcor2 <- round(cor(mod_dat2$Dep_Ndiff_SO, mod_dat2$Dep_Nred),2)

mod_dat3 <- dat %>%
  select(Dep_Nhistoric_SO, Dep_Noxi) %>%
  distinct(.)
pcor3 <- round(cor(mod_dat3$Dep_Nhistoric_SO, mod_dat3$Dep_Noxi),2)

mod_dat4 <- dat %>%
  select(Dep_Ndiff_SO, Dep_Noxi) %>%
  distinct(.)
pcor4 <- round(cor(mod_dat4$Dep_Ndiff_SO, mod_dat4$Dep_Noxi),2)

# Plotting
p <- ggplot()+
  geom_point(data = mod_dat, aes(x = Dep_Nred, y = Dep_Nhistoric_SO),alpha = 0.3)+
  #geom_line(data = fitted, aes(x = Dep_Nhistoric, y = fitted), color = "blue")+
  annotate("text", label = "r = 0.64", x = 5, y = 20, color = "black")+
  #annotate("text", label = expression(paste("p-value of slope term < 2.2 x ",10^-16)), x = 40, y = 20, color = "black")+
  theme_classic()+
  ylab(expression(paste("residual antecedent N deposition (kg N ", ha^-1," ",y^-1,")")))+
  xlab(expression(paste("reduced N deposition (N",H[3],"/N",H[4]^{"+"},"; kg N ", ha^-1," ",y^-1,")")))+
  ggtitle("a.")+
  theme(title = element_text(face = "bold"))
p

p2 <- ggplot()+
  geom_point(data = mod_dat2, aes(x = Dep_Nred, y = Dep_Ndiff_SO),alpha = 0.3)+
  #geom_line(data = fitted2, aes(x = Dep_Shistoric, y = fitted), color = "blue")+
  annotate("text", label = "r = 0.4", x = 5, y = 20, color = "black")+
  #annotate("text", label = expression(paste("p-value of slope term < 2.2 x ",10^-16)), x = 30, y = 55, color = "black")+
  theme_classic()+
  ylab(expression(paste("residual short-term change in N deposition (kg N ", ha^-1," ",y^-1,")")))+
  xlab(expression(paste("reduced N deposition (N",H[3],"/N",H[4]^{"+"},"; kg N ", ha^-1," ",y^-1,")")))+
  ggtitle("b.")+
  theme(title = element_text(face = "bold"))
p2

p3 <- ggplot()+
  geom_point(data = mod_dat3, aes(x = Dep_Noxi, y = Dep_Nhistoric_SO),alpha = 0.3)+
  #geom_line(data = fitted, aes(x = Dep_Nhistoric, y = fitted), color = "blue")+
  annotate("text", label = "r = 0.23", x = 10, y = 20, color = "black")+
  #annotate("text", label = expression(paste("p-value of slope term < 2.2 x ",10^-16)), x = 40, y = 20, color = "black")+
  theme_classic()+
  ylab(expression(paste("residual antecedent N deposition (kg N ", ha^-1," ",y^-1,")")))+
  xlab(expression(paste("oxidized N deposition (N",O[x],"; kg N ", ha^-1," ",y^-1,")")))+
  ggtitle("c.")+
  theme(title = element_text(face = "bold"))
p3

p4 <- ggplot()+
  geom_point(data = mod_dat4, aes(x = Dep_Noxi, y = Dep_Ndiff_SO),alpha = 0.3)+
  #geom_line(data = fitted2, aes(x = Dep_Shistoric, y = fitted), color = "blue")+
  annotate("text", label = "r = 0.01", x = 10, y = 20, color = "black")+
  #annotate("text", label = expression(paste("p-value of slope term < 2.2 x ",10^-16)), x = 30, y = 55, color = "black")+
  theme_classic()+
  ylab(expression(paste("residual short-term change in N deposition (kg N ", ha^-1," ",y^-1,")")))+
  xlab(expression(paste("oxidized N deposition (N",O[x],"; kg N ", ha^-1," ",y^-1,")")))+
  ggtitle("d.")+
  theme(title = element_text(face = "bold"))
p4

p_final <- ggarrange(p, p2, p3, p4, ncol = 2, nrow = 2)

ggsave(p_final, filename = "./visualizations/final_figures/FigureSX.png", dev = "png",
       height = 10, width = 12, units = "in", bg = "white")
