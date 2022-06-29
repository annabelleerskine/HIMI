library(tidyverse)
library(dplyr)

lme_resource <- readRDS("lme_resource_ts.rds")

lme_resource_antarctica <- filter(lme_resource, lme== "61")


#historical 
lme_antarctica_hist <- filter(lme_resource_antarctica, scenario== "historical")
average_hist_L <- lme_antarctica_hist %>% 
  mutate(average_lamda = ave(lambda))

average_hist_K <- lme_antarctica_hist %>% 
  mutate(average_kappa = ave(kappa))


#RCP8.5
lme_antarctica_RCP85 <- filter(lme_resource_antarctica, scenario== "rcp85")

average_RCP85_L <- lme_antarctica_RCP85 %>% 
  mutate(average_lamda = ave(lambda))

average_RCP85_K <- lme_antarctica_RCP85 %>% 
  mutate(average_kappa = ave(kappa))

#RCP all scenarios average

Lme_antarctica_future <- filter(lme_resource_antarctica, !scenario== "historical")
Average_future_L <- Lme_antarctica_future %>% 
  mutate(average_lamda = ave(lambda))
Average_future_K <- Lme_antarctica_future %>% 
  mutate(average_kappa = ave(kappa))
