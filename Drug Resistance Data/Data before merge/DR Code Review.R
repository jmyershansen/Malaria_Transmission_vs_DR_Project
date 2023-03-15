library(tidyverse)
library(dplyr)
library(lme4)
library(lmerTest)
library(ggplot2)
library(lattice)


# PFCRT & PFMDR1 POINT MUTATIONS

WWARN <- read_csv("WWARN_ACT_Partner_Drug_Mol_Surveyor_Data(6).csv")

WWARN <- WWARN %>% distinct(lon,lat)

LATLONG <- as.data.frame(cbind(WWARN$lat,WWARN$lon))

colnames(LATLONG)<-c("latitude","longitude")

write_csv(LATLONG,file="LatitudeLongitudeWWARN.csv")

EXTR <- read.csv("Extracted-points-data.csv")

WWARN_full <- read_csv("WWARN_ACT_Partner_Drug_Mol_Surveyor_Data(6).csv")

WWARN_full <- cbind(WWARN_full$`study start`,WWARN_full)

colnames(WWARN_full)[1] <- "year"

colnames(WWARN_full)[6] <- "long"

WWARN_full <- WWARN_full %>% arrange(year) %>% group_by(lat,long,year)

WWARN_full$long <- round(WWARN_full$long, digits = 2)
EXTR$long <- round (EXTR$long, digits =2)

WWARN_full$lat <- round(WWARN_full$lat, digits = 2)
EXTR$lat <- round(EXTR$lat, digits = 2)

Complete <- right_join(EXTR,WWARN_full,by = c("lat","long","year"))

Complete <- Complete %>% mutate(PR = as.numeric(value)*100)

Complete_pfcrt <- Complete %>% filter(`marker group`=="pfcrt 76T")

Complete_pfmdr1 <- Complete %>% filter(`marker group`=="pfmdr1 86Y")

#write_csv(Complete, file = "pfcrt_pfmdr1_pm.csv")
#write_csv(Complete_pfcrt, file = "pfcrt_pm.csv")
#write_csv(Complete_pfmdr1, file = "pfmdr1_pm.csv")

######
Old_data_pfcrt <- read_csv("WWARN_pfcrt.csv")

Country_region <- Old_data_pfcrt[c(3,5)] %>% distinct(country, Region) %>% arrange(country)

Complete_pfcrt_region <- right_join(Country_region, Complete_pfcrt, by = "country")

#write_csv(Complete_pfcrt_region, file = "pfcrt_pm_region.csv")
write_csv(Complete_pfcrt_region, file = "pfcrt_SNP.csv")
######

Complete_pfmdr1_region <- right_join(Country_region, Complete_pfmdr1, by = "country")
#write_csv(Complete_pfmdr1_region, file = "pfmdr1_pm_region.csv")
write_csv(Complete_pfmdr1_region, file = "pfmdr1_SNP.csv")
######


# PFDHPS & PFDHFR POINT MUTATIONS

DHPS.DHFR <- read_csv("WWARN_pfdhfr_pfdhps_data.csv")

DHPS.DHFR <- cbind(DHPS.DHFR$`study start year`,DHPS.DHFR)

colnames(DHPS.DHFR)[1] <- "year"

colnames(DHPS.DHFR)[7] <- "long"

DHPS.DHFR <- DHPS.DHFR %>% arrange(year) %>% group_by(lat,long,year)

DHPS.DHFR$long <- round(DHPS.DHFR$long, digits = 2)
DHPS.DHFR$lat <- round(DHPS.DHFR$lat, digits = 2)

Complete3 <- right_join(EXTR,DHPS.DHFR,by = c("lat","long","year"))

Complete3 <- Complete3 %>% mutate(PR = as.numeric(value)*100)

Complete_pfdhps <- Complete3 %>% filter(`mutation`=="dhps 437G")
Complete_pfdhfr <- Complete3 %>% filter(`mutation`=="dhfr 108N")

Complete_pfdhps_region <- right_join(Country_region, Complete_pfdhps, by = "country")
#write_csv(Complete_pfdhps_region, file = "pfdhps_pm_region.csv")
write_csv(Complete_pfdhps_region, file = "pfdhps_SNP.csv")

Complete_pfdhfr_region <- right_join(Country_region, Complete_pfdhfr, by = "country")
#write_csv(Complete_pfdhfr_region, file = "pfdhfr_pm_region.csv")
write_csv(Complete_pfdhfr_region, file = "pfdhfr_SNP.csv")


### Filtering out Data from 2000-2005
pcrt_SNP <- read_csv("pfcrt_SNP.csv")
pcrt_SNP_year1a <- pcrt_SNP %>% filter (year%in%(2000:2005))
pcrt_SNP_year1b <- pcrt_SNP %>% filter (year%in%(2006:2021))
write_csv(pcrt_SNP_year1a, file = "pfcrt_SNP.1a.csv")
write_csv(pcrt_SNP_year1b, file = "pfcrt_SNP.1b.csv")

Numstudies0 <- length(unique(pcrt_SNP$`study Id`))
Numstudies1 <- length(unique(pcrt_SNP_year1a$`study Id`))
Numstudies2 <- length(unique(pcrt_SNP_year1b$`study Id`))


pfmdr1_SNP <- read_csv("pfmdr1_SNP.csv")
pfmdr1_SNP_year1a <- pmdr1_SNP %>% filter (year%in%(2000:2005))
pfmdr1_SNP_year1b <- pmdr1_SNP %>% filter (year%in%(2006:2021))
write_csv(pfmdr1_SNP_year1a, file = "pfmdr1_SNP.1a.csv")
write_csv(pfmdr1_SNP_year1b, file = "pfmdr1_SNP.1b.csv")

Numstudies0 <- length(unique(pfmdr1_SNP$`study Id`))
Numstudies1 <- length(unique(pfmdr1_SNP_year1a$`study Id`))
Numstudies2 <- length(unique(pfmdr1_SNP_year1b$`study Id`))

#######
pfdhps_SNP.1<- read_csv("pfdhps_SNP.csv")
write_csv(pfdhps_SNP.1, file = "pfdhps_SNP.1.csv")
Numstudies0 <- length(unique(pfdhps_SNP.1$`study id`))

pfdhfr_SNP.1<- read_csv("pfdhfr_SNP.csv")
write_csv(pfdhfr_SNP.1, file = "pfdhfr_SNP.1.csv")
Numstudies0 <- length(unique(pfdhps_SNP.1$`study id`))


pfcrt_SNP1a<- read.csv("pfmdr1_SNP.1a.csv")

#######
pfcrt_p<- read.csv("pfcrt_SNP.1a.csv")
pfcrt_p <- pfcrt_p %>% mutate(Prev_DR=(present/tested)*100)
write_csv(pfcrt_p, file = "pfcrt_SNP.1a1.csv")

pfcrt_p<- read.csv("pfcrt_SNP.1b.csv")
pfcrt_p <- pfcrt_p %>% mutate(Prev_DR=(present/tested)*100)
write_csv(pfcrt_p, file = "pfcrt_SNP.1b1.csv")

pfmdr1_p<- read.csv("pfmdr1_SNP.1a.csv")
pfmdr1_p <- pfmdr1_p %>% mutate(Prev_DR=(present/tested)*100)
write_csv(pfmdr1_p, file = "pfmdr1_SNP.1a1.csv")

pfmdr1_prop<- read.csv("pfmdr1_SNP.1b.csv")
pfmdr1_prop <- pfmdr1_prop %>% mutate(Prev_DR=(present/tested)*100)
write_csv(pfmdr1_prop, file = "pfmdr1_SNP.1b1.csv")

pfdhps_p<- read.csv("pfdhps_SNP.1.csv")
pfdhps_p<- pfdhps_p %>% mutate(Prev_DR=(present/tested)*100)
write_csv(pfdhps_p, file = "pfdhps_SNP.1I.csv")


pfdhfr_p<- read.csv("pfdhfr_SNP.1.csv")
pfdhfr_p<- pfdhfr_p %>% mutate(Prev_DR=(present/tested)*100)
write_csv(pfdhfr_p, file = "pfdhfr_SNP.1I.csv")


######
######
#Models for all sample sizes

###pfcrt-76T

pfcrt_SNP1a<- read.csv("pfcrt_SNP.1a1.csv")
# assign factors
pfcrt_SNP1a$study.Id <- as.factor(pfcrt_SNP1a$study.Id)
pfcrt_SNP1a$country <- as.factor(pfcrt_SNP1a$country)

Meta1.pfcrt_SNP1a <- lmer(Prev_DR ~ year + PR + (1|country/study.Id), data=pfcrt_SNP1a)
summary(Meta1.pfcrt_SNP1a)

Meta2.pfcrt_SNP1a <- lmer(Prev_DR ~ year + PR + (1|study.Id), data=pfcrt_SNP1a)
summary(Meta2.pfcrt_SNP1a)

###pfmdr1-86Y

pfmdr1_SNP1a<- read.csv("pfmdr1_SNP.1a1.csv")
# assign factors
pfmdr1_SNP1a$study.Id <- as.factor(pfmdr1_SNP1a$study.Id)
pfmdr1_SNP1a$country <- as.factor(pfmdr1_SNP1a$country)

Meta1.pfmdr1_SNP1a <- lmer(Prev_DR ~ year + PR + (1|country/study.Id), data=pfmdr1_SNP1a)
summary(Meta1.pfmdr1_SNP1a)

Meta2.pfmdr1_SNP1a <- lmer(Prev_DR ~ year + PR + (1|study.Id), data=pfmdr1_SNP1a)
summary(Meta2.pfmdr1_SNP1a)

###pfdhps-437G

pfdhps_SNP1<- read.csv("pfdhps_SNP.1I.csv")
# assign factors
pfdhps_SNP1$study.id <- as.factor(pfdhps_SNP1$study.id)
pfdhps_SNP1$country <- as.factor(pfdhps_SNP1$country)

Meta1.pfdhps_SNP1 <- lmer(Prev_DR ~ year + PR + (1|country/study.id), data=pfdhps_SNP1)
summary(Meta1.pfdhps_SNP1)
anova(Meta1.pfdhps_SNP1, ddf= "Kenward-Roger")

Meta2.pfdhps_SNP1 <- lmer(Prev_DR ~ year + PR + (1|study.id), data=pfdhps_SNP1)
summary(Meta2.pfdhps_SNP1)


###pfdhfr-108N

pfdhfr_SNP1<- read.csv("pfdhfr_SNP.1I.csv")
# assign factors
pfdhfr_SNP1$study.id <- as.factor(pfdhfr_SNP1$study.id)
pfdhfr_SNP1$country <- as.factor(pfdhfr_SNP1$country)

Meta1.pfdhfr_SNP1 <- lmer(Prev_DR ~ year + PR + (1|country/study.id), data=pfdhfr_SNP1)
summary(Meta1.pfdhfr_SNP1)
anova(Meta1.pfdhfr_SNP1, ddf= "Kenward-Roger")

Meta2.pfdhfr_SNP1 <- lmer(Prev_DR ~ year + PR + (1|study.id), data=pfdhfr_SNP1)
summary(Meta2.pfdhfr_SNP1)


####
####
#Models for sample sizes >50

pfcrt_76T_new_data_1<- read_csv("pfcrt_SNP.1a1.csv")
pfcrt_76T_new_data_2 <- pfcrt_76T_new_data_1 %>% filter(tested>50)
write_csv(pfcrt_76T_new_data_2, file = "pfcrt_SNP.1a3.csv")

pfcrt_SNP.1a3<- read_csv("pfcrt_SNP.1a3.csv")
Meta1.pfcrt_SNP1a3 <- lmer(Prev_DR ~ year + PR + (1|country/study.Id), data=pfcrt_SNP.1a3)
summary(Meta1.pfcrt_SNP1a3)
Meta2.pfcrt_SNP1a3 <- lmer(Prev_DR ~ year + PR + (1|study.Id), data=pfcrt_SNP.1a3)
summary(Meta2.pfcrt_SNP1a3)


pfmdr1_86Y_new_data_1<- read_csv("pfmdr1_SNP.1a1.csv")
pfmdr1_86Y_new_data_2 <- pfmdr1_86Y_new_data_1 %>% filter(tested>50)
write_csv(pfmdr1_86Y_new_data_2, file = "pfmdr1_SNP.1a3.csv")

pfmdr1_SNP.1a3<- read_csv("pfmdr1_SNP.1a3.csv")
Meta1.pfmdr1_SNP1a3 <- lmer(Prev_DR ~ year + PR + (1|country/study.Id), data=pfmdr1_SNP.1a3)
summary(Meta1.pfmdr1_SNP1a3)
Meta2.pfmdr1_SNP1a3 <- lmer(Prev_DR ~ year + PR + (1|study.Id), data=pfmdr1_SNP.1a3)
summary(Meta2.pfmdr1_SNP1a3)


pfdhps_437G_new_data_1<- read_csv("pfdhps_SNP.1I.csv")
pfdhps_437G_new_data_2 <- pfdhps_437G_new_data_1 %>% filter(tested>50)
write_csv(pfdhps_437G_new_data_2, file = "pfdhps_SNP.1II.csv")

pfdhps_SNP.1II<- read_csv("pfdhps_SNP.1II.csv")
Meta1.pfdhps_SNP1II <- lmer(Prev_DR ~ year + PR + (1|country/study.id), data=pfdhps_SNP.1II)
summary(Meta1.pfdhps_SNP1II)
Meta2.pfdhps_SNP1II <- lmer(Prev_DR ~ year + PR + (1|study.id), data=pfdhps_SNP.1II)
summary(Meta2.pfdhps_SNP1II)


pfdhfr_108N_new_data_1<- read_csv("pfdhfr_SNP.1I.csv")
pfdhfr_108N_new_data_2 <- pfdhfr_108N_new_data_1 %>% filter(tested>50)
write_csv(pfdhfr_108N_new_data_2, file = "pfdhfr_SNP.1II.csv")

pfdhfr_SNP.1II<- read_csv("pfdhfr_SNP.1II.csv")
Meta1.pfdhfr_SNP1II <- lmer(Prev_DR ~ year + PR + (1|country/study.id), data=pfdhfr_SNP.1II)
summary(Meta1.pfdhfr_SNP1II)
Meta2.pfdhfr_SNP1II <- lmer(Prev_DR ~ year + PR + (1|study.id), data=pfdhfr_SNP.1II)
summary(Meta2.pfdhfr_SNP1II)

