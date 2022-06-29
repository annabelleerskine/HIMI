##########################################################################################################
### Provide information from HIMI RSTS and LL to calibrate/validate size-spectrum models for TREV      ###
### Calculate total biomass and total effort through time for key HIMI species                         ###
### Calculate spatial overlap of key species                                                           ###
### Author: N.Hill Modified Feb 2021                                                                   ###
##########################################################################################################


## 1) Set up  ----
library(tidyverse)
library(cooccur)

path<- "C:\\Users\\hillna\\OneDrive - University of Tasmania\\UTAS_work\\Projects\\Toothfish FRDC\\"


# load HIMI RSTS biomass data from Joel/ SO Eco Change project
rsts<-read.csv(paste0(path, "Size_spectra\\himi_weights_wide_2003-2016.csv"))
table(is.na(rsts))# 3

##load haul info from TREV
haul<-read.csv(paste0(path, "Analysis\\Fisheries_Data\\Raw_Data\\Raw_Data Haul.csv"), as.is=TRUE)
haul$StartDate<-as.Date(haul$StartDate, format= "%Y-%m-%d")


## 2) Biomass of key species per year -----
##Note: Effort= Area swept by trawl in km2
## biomass of each species in kg

#extract relevant rows + columns, create new ones
haul2<- haul %>% 
        filter(Purpose =="RSTS") %>%
        select( BatchID, StartDate, AreaID, Season, Effort) %>%
        mutate(Month= format(StartDate, "%b"), Year= format(StartDate, "%Y"))

#merge haul and biomass info and check
rsts_bio <- right_join(haul2, rsts)
dim(rsts)
dim(rsts_bio)
table(rsts_bio$Year, rsts_bio$Month)
        
# choose species of interest and aggregate biomass by year
# all strata
rsts_trends<-rsts_bio %>%
        select( Year, Effort,
              Dissostichus.eleginoides, Champsocephalus.gunnari, Channichthys.rhinoceratus, Lepidonotothen.squamifrons, 
              Macrourus.sp, Bathyraja.eatonii, Bathyraja.irrasa, Gobionotothen.acuta, Bathyraja.murrayi) %>%
              
              group_by(Year) %>%
              summarise_all(sum, na.rm=TRUE)

plot(rsts_trends$Year, rsts_trends$Dissostichus.eleginoides/rsts_trends$Effort, type="l")


## 3) Calculate overlap/co-occurrence of each species of interest in hauls ----

# METHOD 1: proportion of hauls across entire dataset where pairs of species co-occur
# convert biomass to PA
PA<- function(x) ifelse(x > 0, 1, 0)

m1_occur<-rsts_bio %>%
          select( Dissostichus.eleginoides, Champsocephalus.gunnari, Channichthys.rhinoceratus, Lepidonotothen.squamifrons, 
                Macrourus.sp, Bathyraja.eatonii, Bathyraja.irrasa, Gobionotothen.acuta, Bathyraja.murrayi) %>%
          mutate_all( PA)

#calculate number of co-occuring sites and divide by total number of sites
rsts_m1_cooccur<-round( create.N.matrix(t(na.omit(m1_occur)))/ dim(t(na.omit(m1_occur)))[2], 3)
row.names(rsts_m1_cooccur) <-colnames(rsts_m1_cooccur) <-names(rsts_m1_cooccur)



# METHOD 2: Following Kempf, A., Dingsør, G.E., Huse, G., Vinther, M., Floeter, J. & Temming,
# A. (2010) The importance of predator – prey overlap: predicting North
# Sea cod recovery with a multispecies assessment model. ICES Journal of Marine Science, 67, 1989–1997.


# unclear if the proportions are of all possible species or only species of interest. 
# is cpue the raw data???
# Need to check with Julia and Romain



## 4) Calculate trends for  HIMI longline data (only have toothfish at this stage) ----

##NOTE: effort is number of hooks and biomass in in T

ll<-read.csv(paste0(path, "Analysis\\Fisheries_Data\\Processed_Data\\HIMI_ll_Haul_QC.csv"))

ll_trends<-ll %>%
  filter(Purpose =="COMM") %>%
  filter(Year<2021) %>%
  select( Year, Effort, CatchT) %>%
  group_by(Year) %>%
  summarise_all(sum)

plot(ll_trends$Year, ll_trends$CatchT/ll_trends$Effort, type="l")


## 4) export data ----
save(rsts_trends, rsts_m1_cooccur, ll_trends, file=paste0(path, "Size_spectra\\size_spec_inputs.RData"))


