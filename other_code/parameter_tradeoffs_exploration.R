# parameter tradeoffs exploration
# Author: Mary Lofton
# Date: 07APR25

# Purpose: explore tradeoffs in estimated model parameters
# 1. nitrogen and sulfur
# 2. reduced N and precipitation

# load packages 
library(tidyverse)
library(lubridate)
library(ggpubr)

# source functions
source("./other_code/get_model_inputs.R")

# load data
dat <- get_model_inputs(data = "./data/McDonnell_etal_InPrep_TreeData_2024_10_11.csv")
colnames(dat)

plot_data <- dat %>%
  mutate(MAP_baseline_dm = 0.01*MAP_baseline) %>%
  select(plot_ID, Dep_Ndelta, Dep_Sdelta, Dep_Noxidelta, Dep_Nreddelta, 
         MAP_delta_dm, MAP_baseline_dm, Dep_Nredbaseline) %>%
  distinct()

mod <- lm(Dep_Ndelta ~ Dep_Sdelta, data = plot_data)
intercept = mod$coefficients[1]
slope = mod$coefficients[2]

fig1 <- ggplot(data = plot_data, aes(x = Dep_Sdelta, y = Dep_Ndelta))+
  geom_point(color = "gray")+
  geom_hline(yintercept = 0, linetype = "dashed")+
  geom_abline(slope = slope, intercept = intercept, color = "blue", linewidth = 1)+
  theme_bw()+
  ylab(expression(paste("N deposition deviation (kg N ", ha^-1," ",y^-1,")")))+
  xlab(expression(paste("S deposition deviation (kg S ", ha^-1," ",y^-1,")")))
fig1
ggsave(plot = fig1, filename = "./visualizations/delta_Ndep_vs_delta_Sdep.png",
       device = "png", height = 5, width = 5.5, units = "in",bg = "white")

mod <- lm(Dep_Noxidelta ~ Dep_Sdelta, data = plot_data)
intercept = mod$coefficients[1]
slope = mod$coefficients[2]

fig2 <- ggplot(data = plot_data, aes(x = Dep_Sdelta, y = Dep_Noxidelta))+
  geom_point(color = "gray")+
  geom_hline(yintercept = 0, linetype = "dashed")+
  geom_abline(slope = slope, intercept = intercept, color = "blue", linewidth = 1)+
  theme_bw()+
  ylab(expression(paste("oxidized N deposition deviation (kg N ", ha^-1," ",y^-1,")")))+
  xlab(expression(paste("S deposition deviation (kg S ", ha^-1," ",y^-1,")")))
fig2
ggsave(plot = fig2, filename = "./visualizations/delta_Noxidep_vs_delta_Sdep.png",
       device = "png", height = 5, width = 5.5, units = "in",bg = "white")

mod <- lm(Dep_Nreddelta ~ Dep_Sdelta, data = plot_data)
intercept = mod$coefficients[1]
slope = mod$coefficients[2]

fig3 <- ggplot(data = plot_data, aes(x = Dep_Sdelta, y = Dep_Nreddelta))+
  geom_point(color = "gray")+
  geom_hline(yintercept = 0, linetype = "dashed")+
  geom_abline(slope = slope, intercept = intercept, color = "blue", linewidth = 1)+
  theme_bw()+
  ylab(expression(paste("reduced N deposition deviation (kg N ", ha^-1," ",y^-1,")")))+
  xlab(expression(paste("S deposition deviation (kg S ", ha^-1," ",y^-1,")")))
fig3
ggsave(plot = fig3, filename = "./visualizations/delta_Nreddep_vs_delta_Sdep.png",
       device = "png", height = 5, width = 5.5, units = "in",bg = "white")

# list files in each model output folder
out1 <- list.files("./experiments/new_delta_Ndep_only_saveTreeEffect",pattern = "mcmc.parquet",
                   full.names = TRUE)
out2 <- list.files("./experiments/delta_env_saveTreeEffect",pattern = "mcmc.parquet",
                   full.names = TRUE)
out3 <- list.files("./experiments/N_species_saveTreeEffect",pattern = "mcmc.parquet",
                   full.names = TRUE)
out <- c(out1,out2, out3)


# read in and combine files
for(i in 1:length(out)){
  
  spp_id = str_split(out[i], pattern = "-")[[1]][2]
  model_id = str_split(out[i], pattern = "/")[[1]][3]
  temp <- read_parquet(file = out[i]) %>%
    mutate(spp_id = spp_id,
           model_id = model_id)
  
  if(i == 1){
    final <- temp
  } else {
    final <- bind_rows(final, temp)
  }
  
}

final1 <- final %>%
  mutate(spp_id = ifelse(spp_id == "yellow","yellow poplar",spp_id)) %>%
  select(model_id, spp_id, global_tree_effect, p2, p3, p4, p5, p6, p7, p8) 

delta_env <- final1 %>%
  filter(model_id == "delta_env_saveTreeEffect" & spp_id == "red spruce")
a <- ggplot(data = delta_env, aes(x = p5))+
  geom_density()+
  theme_classic()+
  xlab("Coefficient on N deposition deviation")
a

b <- ggplot(data = delta_env, aes(x = p6))+
  geom_density()+
  theme_classic()+
  xlab("Coefficient on S deposition deviation")
b

c <- ggplot(data = delta_env, aes(y = p5, x = p6))+
  geom_point()+
  theme_bw()+
  xlab("Coefficient on S deposition deviation")+
  ylab("Coefficient on N deposition deviation")
c

fig4 <- ggarrange(a, b, c,
                nrow = 1, ncol = 3,
                labels = c("A","B","C")
)
fig4
ggsave(plot = fig4, filename = "./visualizations/param_tradeoff_S_N.png",
       device = "png", height = 4, width = 12, units = "in")

mod <- lm(Dep_Noxidelta ~ MAP_delta_dm, data = plot_data)
intercept = mod$coefficients[1]
slope = mod$coefficients[2]

fig5 <- ggplot(data = plot_data, aes(x = MAP_delta_dm, y = Dep_Noxidelta))+
  geom_point(color = "gray")+
  geom_hline(yintercept = 0, linetype = "dashed")+
  geom_abline(slope = slope, intercept = intercept, color = "blue", linewidth = 1)+
  theme_bw()+
  ylab(expression(paste("oxidized N deposition deviation (kg N ", ha^-1," ",y^-1,")")))+
  xlab(expression(paste("precipitation deviation (dm ", y^-1,")")))
fig5
ggsave(plot = fig5, filename = "./visualizations/delta_Nreddep_vs_delta_precip.png",
       device = "png", height = 5, width = 5.5, units = "in",bg = "white")

mod <- lm(Dep_Noxibaseline ~ MAP_baseline_dm, data = plot_data)
intercept = mod$coefficients[1]
slope = mod$coefficients[2]

fig6 <- ggplot(data = plot_data, aes(x = MAP_baseline_dm, y = Dep_Noxibaseline))+
  geom_point(color = "gray")+
  geom_abline(slope = slope, intercept = intercept, color = "blue", linewidth = 1)+
  theme_bw()+
  ylab(expression(paste("mean reduced N deposition (kg N ", ha^-1," ",y^-1,")")))+
  xlab(expression(paste("mean precipitation (dm ", y^-1,")")))
fig6
ggsave(plot = fig6, filename = "./visualizations/baseline_Nreddep_vs_baseline_precip.png",
       device = "png", height = 5, width = 5.5, units = "in",bg = "white")

N_species <- final1 %>%
  filter(model_id == "N_species_saveTreeEffect" & spp_id == "red spruce")
a <- ggplot(data = N_species, aes(x = p5))+
  geom_density()+
  theme_classic()+
  xlab("Coefficient on reduced N deposition deviation")
a

b <- ggplot(data = N_species, aes(x = p8))+
  geom_density()+
  theme_classic()+
  xlab("Coefficient on precipitation deviation")
b

c <- ggplot(data = N_species, aes(y = p5, x = p8))+
  geom_point()+
  theme_bw()+
  xlab("Coefficient on precipitation deviation")+
  ylab("Coefficient on reduced N deposition deviation")
c

fig7 <- ggarrange(a, b, c,
                  nrow = 1, ncol = 3,
                  labels = c("A","B","C")
)
fig7
ggsave(plot = fig7, filename = "./visualizations/param_tradeoff_Nred_precip.png",
       device = "png", height = 5, width = 15, units = "in")

