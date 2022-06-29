#Sensitivity

library(mizer)
library(viridisLite)
library(viridis)
library(ggplot2)
library(ggrepel)
library(tidyverse)
library(mizerExperimental)
library(mizerHowTo)

sim_loop4 <- readRDS("sim_loop4.rds")

#grid for +-5% of all 
library(data.table)
kappa_5 <- c(0.256,0.2688,0.243)
lambda_5 <- c(2.05,2.153,1.948)
r_pp_5 <- c(10, 10.5,9.5)

grid <- expand.grid(kappa = kappa_5, lambda= lambda_5, r_pp =r_pp_5)
grid$name <- as.factor(1:27)

#getting output from experiments
params <- sim_loop4@params #  replace with your param file
max_t <- 50 # How long you want sims to run
output <- # create output object with grid
  grid %>%
  group_by(name)%>% # group by model scenario (i.e. name)
  nest() %>% # create nested object
  mutate(param_list = map(data, ~setResource(params, kappa = .$kappa,lambda = .$lambda))) %>% 
  mutate(param_list2 = map(param_list, ~projectToSteady(., effort=0, t_max= max_t, progress_bar= F))) %>% # run model projection
  mutate(mod_output = map(param_list2, ~project(., effort=0, t_max= max_t, progress_bar= F))) %>% # run model projection
  mutate(biomass = map(mod_output, ~reshape::melt(getBiomass(.))))
output %>%
  select(name, biomass) %>% #choose only the name and biomass dataframes
  unnest(biomass) %>%
  fwrite("df_Biomass_sens.csv") # save a biomass dataframe with all the combined scenarios

df_biomass <- read.csv("df_Biomass_sens.csv")
df_biomass$name <- as.factor(df_biomass$name)
grid$name <- as.factor(grid$name)
df_plot <- left_join(df_biomass, grid)


df_plot_1 <- df_plot %>%
  filter(time == max_t) %>% # only use the final year
  filter(name == 1 | name == 2 | name == 3 | name== 4 | name == 7 | name == 10 | name == 19) %>% 
  ggplot(aes(x=value, y =sp, colour = name)) +
  geom_point() +
  # facet_wrap(~name) + # separate plot for each scenario
  scale_x_log10()
df_plot_1

#get total biomass (not species specific)
df_plot_2 <- df_plot %>% 
  filter(time== max_t) %>% 
  group_by(name) %>% 
  summarise(totalbiomass= sum(value)) %>% 
  ungroup() %>% 
  dplyr::mutate(name_type = case_when(name == "1" ~ "Base",
                                       name == "2" ~ "kappa",
                                       name == "3" ~ "kappa",
                                       name == "4" ~ "lambda",
                                       name == "7" ~ "lambda",
                                       name == "10" ~ "r_pp",
                                       name == "19" ~ "r_pp",
                                       name == "5" ~ "KL",
                                       name == "6" ~ "KL",
                                       name == "8" ~ "KL",
                                       name == "9" ~ "KL",
                                       name == "11" ~ "KR",
                                       name == "12" ~ "KR",
                                       name == "13" ~ "KL",
                                      name == "14" ~ "KLR",
                                      name == "15" ~ "KLR",
                                      name == "16" ~ "LR",
                                      name == "17" ~ "KLR",
                                      name == "18" ~ "KLR",
                                      name == "20" ~ "KR",
                                      name == "21" ~ "KR",
                                      name == "22" ~ "LR",
                                      name == "23" ~ "KLR",
                                      name == "24" ~ "KLR",
                                      name == "25" ~ "LR",
                                      name == "26" ~ "KLR",
                                      name == "27" ~ "KLR",
                                      ))

#no interactions
plot_3 <- df_plot_2 %>% filter(name == 1 | name == 2 | name == 3 | name== 4 | name == 7 | name == 10 | name == 19) %>% 
ggplot(aes(x= name, y= totalbiomass, colour= name_type))+
  geom_point()
plot_3
  # facet_wrap(~name) + # separate plot for each scenario


#make table for biomass outputs







