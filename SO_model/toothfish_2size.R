library(ggplot2)
library(ggrepel)
library(tidyverse)
library(igraph)

df_diet <- read.csv("df_diet_d.csv")

df_diet_Base <- filter(df_diet, name== "Base")

df_diet_LI <- filter(df_diet, name== "LI")
df_diet_mawsoni_LI <- filter(df_diet_LI, predator== "D.mawsoni")


df_diet_wend <- filter(df_diet, w== 94800)

