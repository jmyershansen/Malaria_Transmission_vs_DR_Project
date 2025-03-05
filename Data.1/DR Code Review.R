install.packages ("tidyverse")
install.packages ("dplyr")
install.packages ("lme4")
install.packages ("lmerTest")
install.packages ("ggplot2")
install.packages ("betareg")
install.packages("imputeTS")

library(imputeTS)
library(tidyverse)
library(dplyr)
library(lme4)
library(lmerTest)
library(ggplot2)
library(lattice)
library(betareg)


#setwd("/Users/jmyershansengmail.com/Documents/NEW COMPUTER FILES_RECENT/DATA/Purdue University/Purdue_Fall 2022/BIO 69500 (Data Science for Biologists)/Project/BIO 69500_Project/CQ Markers/New folder/Drug Resistance Data_vital/Data.1")

pfcrt_ID <- read_csv("pfcrt_SNP.csv")
pfcrt_ID <- pfcrt_ID %>% distinct(`study Id`,`author`,`publication year`,`publication URL`)
pfcrt_ID  <- as.data.frame(cbind(pfcrt_ID$`study Id`,pfcrt_ID$`author`,pfcrt_ID$`publication year`,pfcrt_ID$`publication URL`))
colnames(pfcrt_ID )<-c("Id","author","year","URL")
NStudies <- length(unique(pfcrt_ID$`Id`))
write_csv(pfcrt_ID,file="pfcrt_year.csv")
###OR
pfcrt_ID <- read_csv("pfcrt_SNP.csv") %>% as_tibble()
pfcrt_ID <- pfcrt_ID %>% select(`study Id`, country, author, `publication year`, `publication URL`)
pfcrt_ID <- pfcrt_ID %>% distinct()
#pfcrt_ID <- cbind(pfcrt_ID$`study Id`,pfcrt_ID$`author`,pfcrt_ID$`publication year`)
NStudies <- length(unique(pfcrt_ID$`study Id`))
write_csv(pfcrt_ID,file="author_pfcrt.csv")
#!duplicated(pfcrt_ID$`study Id`)

pfmdr1_ID <- read_csv("pfmdr1_SNP.csv") %>% as_tibble()
pfmdr1_ID <- pfmdr1_ID %>% select(`study Id`, country, author, `publication year`, `publication URL`)
pfmdr1_ID <- pfmdr1_ID %>% distinct()
NStudies <- length(unique(pfmdr1_ID$`study Id`))
write_csv(pfmdr1_ID,file="author_pfmdr1.csv")

pfdhps_ID <- read_csv("pfdhps_SNP.csv") %>% as_tibble()
pfdhps_ID <- pfdhps_ID %>% select(`study id`, country, authors, `publication year`, `publication Url`)
pfdhps_ID <- pfdhps_ID %>% distinct()
NStudies <- length(unique(pfdhps_ID$`study id`))
write_csv(pfdhps_ID,file="author_pfdhps.csv")

pfdhfr_ID <- read_csv("pfdhfr_SNP.csv") %>% as_tibble()
pfdhfr_ID <- pfdhfr_ID %>% select(`study id`, country, authors, `publication year`, `publication Url`)
pfdhfr_ID <- pfdhfr_ID %>% distinct()
NStudies <- length(unique(pfdhfr_ID$`study id`))
write_csv(pfdhfr_ID,file="author_pfdhfr.csv")
#######
#######
#######
#######

# PFCRT & PFMDR1 POINT MUTATIONS

WWARN <- read_csv("WWARN_ACT_Partner_Drug_Mol_Surveyor_Data(6).csv")

WWARN <- WWARN %>% distinct(lon,lat)

LATLONG <- as.data.frame(cbind(WWARN$lat,WWARN$lon))

colnames(LATLONG)<-c("latitude","longitude")

write_csv(LATLONG,file="LatitudeLongitudeWWARN.csv")

EXTR <- read.csv("Extracted-points-data.csv")

############## 
WWARN1 <- read_csv("WWARN_pfdhfr_pfdhps_data.csv")

WWARN1 <- WWARN1 %>% distinct(lon,lat)

LATLONG1 <- as.data.frame(cbind(WWARN1$lat,WWARN1$lon))

colnames(LATLONG1)<-c("latitude","longitude")

write_csv(LATLONG1,file="LatitudeLongitudeWWARN.1.csv")

EXTR1 <- read.csv("Extracted-points-data.1.csv")

#########
WWARN2 <- read_csv("K13_surveyor_data.csv")

WWARN2 <- WWARN2 %>% distinct(lon,lat)

LATLONG2 <- as.data.frame(cbind(WWARN2$lat,WWARN2$lon))

colnames(LATLONG2)<-c("latitude","longitude")

write_csv(LATLONG2,file="LatitudeLongitudeWWARN.2.csv")

EXTR2 <- read.csv("Extracted-points-data.2.csv")

#########
WWARN_full <- read_csv("WWARN_ACT_Partner_Drug_Mol_Surveyor_Data(6).csv")

WWARN_full <- cbind(WWARN_full$`study end`,WWARN_full)

colnames(WWARN_full)[1] <- "year"

colnames(WWARN_full)[6] <- "long"

WWARN_full <- WWARN_full %>% arrange(year) %>% group_by(lat,long,year)

WWARN_full$long <- round(WWARN_full$long, digits = 2)
EXTR$long <- round (EXTR$long, digits =2)

WWARN_full$lat <- round(WWARN_full$lat, digits = 2)
EXTR$lat <- round(EXTR$lat, digits = 2)

Complete <- right_join(EXTR,WWARN_full,by = c("lat","long","year"))

Complete <- Complete %>% mutate(PR = as.numeric(value)*100)
View(Complete)
#Numstudies6 <- length(unique(Complete$`study Id`))
#Numstudies6i <- length(unique(Complete_pfcrt$`study Id`))
#Numstudies6ii <- length(unique(Complete_pfmdr1$`study Id`))


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
#####
#####

pfcrt_SNP <- read_csv ("pfcrt_SNP.csv")
pfcrt_SNP <- pfcrt_SNP %>% mutate(Prev_DR=(present/tested)*100)
write_csv(pfcrt_SNP, file = "pfcrt_SNP.csv")

DUdata <- read_csv ("FeveredChildrenDrugTreatProp_CQ.0.csv")
pfcrt_complete <- right_join(DUdata,pfcrt_SNP, by = c("country","year"))
pfcrt_complete <- pfcrt_complete %>% mutate(CQ_use= as.numeric(CQ)*100)
write_csv(pfcrt_complete, file = "pfcrt_SNP_complete.csv")
pfcrt_SNP.1<- read_csv("pfcrt_SNP_complete.csv")
#write_csv(pfcrt_SNP.1, file = "pfcrt_SNP.1.csv")
Numstudies0 <- length(unique(pfcrt_SNP.1$`CQ_use`))


pfmdr1_SNP <- read_csv ("pfmdr1_SNP.csv")
pfmdr1_SNP <- pfmdr1_SNP %>% mutate(Prev_DR=(present/tested)*100)
write_csv(pfmdr1_SNP, file = "pfmdr1_SNP.csv")

DUdata <- read_csv ("FeveredChildrenDrugTreatProp_CQ.0.csv")
pfmdr1_complete <- right_join(DUdata,pfmdr1_SNP, by = c("country","year"))
pfmdr1_complete <- pfmdr1_complete %>% mutate(CQ_use= as.numeric(CQ)*100)
write_csv(pfmdr1_complete, file = "pfmdr1_SNP_complete.csv")
pfmdr1_SNP.1<- read_csv("pfmdr1_SNP_complete.csv")
#write_csv(pfmdr1_SNP.1, file = "pfmdr1_SNP.1.csv")
Numstudies0 <- length(unique(pfmdr1_SNP.1$`CQ_use`))
######
#####
#####


# PFDHPS & PFDHFR POINT MUTATIONS

DHPS.DHFR <- read_csv("WWARN_pfdhfr_pfdhps_data.csv")

DHPS.DHFR <- cbind(DHPS.DHFR$`study end year`,DHPS.DHFR)

colnames(DHPS.DHFR)[1] <- "year"

colnames(DHPS.DHFR)[7] <- "long"

DHPS.DHFR <- DHPS.DHFR %>% arrange(year) %>% group_by(lat,long,year)

DHPS.DHFR$long <- round(DHPS.DHFR$long, digits = 2)
EXTR1$long <- round (EXTR1$long, digits =2)

DHPS.DHFR$lat <- round(DHPS.DHFR$lat, digits = 2)
EXTR1$lat <- round(EXTR1$lat, digits = 2)

Complete3 <- right_join(EXTR1,DHPS.DHFR,by = c("lat","long","year"))

Complete3 <- Complete3 %>% mutate(PR = as.numeric(value)*100)
#Numstudies7 <- length(unique(Complete3$`study id`))

Complete_pfdhps <- Complete3 %>% filter(`mutation`=="dhps 437G")
Complete_pfdhfr <- Complete3 %>% filter(`mutation`=="dhfr 108N")
#Numstudies8 <- length(unique(Complete_pfdhfr$`study id`))
#Numstudies9 <- length(unique(Complete_pfdhps$`study id`))

Complete_pfdhps_region <- right_join(Country_region, Complete_pfdhps, by = "country")
#write_csv(Complete_pfdhps_region, file = "pfdhps_pm_region.csv")
write_csv(Complete_pfdhps_region, file = "pfdhps_SNP.csv")

pfdhps_SNP <- read_csv ("pfdhps_SNP.csv")
pfdhps_SNP <- pfdhps_SNP %>% mutate(Prev_DR=(present/tested)*100)
write_csv(pfdhps_SNP, file = "pfdhps_SNP.csv")
Numstudies0 <- length(unique(pfdhps_SNP$`study id`))
Numstudies0 <- length(unique(pfdhps_SNP$`PR`))
Numstudies0 <- length(unique(pfdhps_SNP$`Prev_DR`))

DUdata <- read_csv ("FeveredChildrenDrugTreatProp_SP.0.csv")
pfdhps_complete <- right_join(DUdata,pfdhps_SNP, by = c("country","year"))
pfdhps_complete <- pfdhps_complete %>% mutate(SP_use= as.numeric(SP)*100)
write_csv(pfdhps_complete, file = "pfdhps_SNP_complete.csv")
pfdhps_SNP.1<- read_csv("pfdhps_SNP_complete.csv")
write_csv(pfdhps_SNP.1, file = "pfdhps_SNP.1.csv")
Numstudies0 <- length(unique(pfdhps_SNP.1$`SP_use`))


Complete_pfdhfr_region <- right_join(Country_region, Complete_pfdhfr, by = "country")
#write_csv(Complete_pfdhfr_region, file = "pfdhfr_pm_region.csv")
write_csv(Complete_pfdhfr_region, file = "pfdhfr_SNP.csv")

pfdhfr_SNP <- read_csv ("pfdhfr_SNP.csv")
pfdhfr_SNP <- pfdhfr_SNP %>% mutate(Prev_DR=(present/tested)*100)
write_csv(pfdhfr_SNP, file = "pfdhfr_SNP.csv")
Numstudies0 <- length(unique(pfdhfr_SNP$`study id`))
Numstudies0 <- length(unique(pfdhfr_SNP$`PR`))
Numstudies0 <- length(unique(pfdhfr_SNP$`Prev_DR`))

DUdata <- read_csv ("FeveredChildrenDrugTreatProp_SP.0.csv")
pfdhfr_complete <- right_join(DUdata,pfdhfr_SNP, by = c("country","year"))
pfdhfr_complete <- pfdhfr_complete %>% mutate(SP_use= as.numeric(SP)*100)
write_csv(pfdhfr_complete, file = "pfdhfr_SNP_complete.csv")
pfdhfr_SNP.1<- read_csv("pfdhfr_SNP_complete.csv")
write_csv(pfdhfr_SNP.1, file = "pfdhfr_SNP.1.csv")
Numstudies0 <- length(unique(pfdhfr_SNP.1$`SP_use`))
#####
#####

# PFK13 POINT MUTATIONS

K13 <- read_csv("K13_surveyor_data.csv")

colnames(K13)[7] <- "long"

K13 <- K13 %>% arrange(year) %>% group_by(lat,long,year)

K13$long <- round(K13$long, digits = 2)
EXTR2$long <- round (EXTR2$long, digits =2)

K13$lat <- round(K13$lat, digits = 2)
EXTR2$lat <- round(EXTR2$lat, digits = 2)

Complete4 <- right_join(EXTR2,K13,by = c("lat","long","year"))

Complete4 <- Complete4 %>% mutate(PR = as.numeric(value)*100)

Complete_pfK13 <- Complete4 %>% filter(`mutation`=="C580Y")
#Numstudies8 <- length(unique(Complete_pfK13$`sid`))

Complete_pfK13_region <- right_join(Country_region, Complete_pfK13, by = "country")
write_csv(Complete_pfK13_region, file = "pfK13_SNP.csv")

pfK13_SNP <- read_csv ("pfK13_SNP.csv")
pfK13_SNP <- pfK13_SNP %>% mutate(Prev_DR=(present/tested)*100)
write_csv(pfK13_SNP, file = "pfK13_SNP.csv")
Numstudies0 <- length(unique(pfK13_SNP$`sid`))
Numstudies1 <- length(unique(pfK13_SNP$`PR`))
Numstudies2 <- length(unique(pfK13_SNP$`Prev_DR`))


DUdata <- read_csv ("ACT_use.0.csv")
pfK13_complete <- right_join(DUdata,pfK13_SNP, by = c("country","year"))
pfK13_complete <- pfK13_complete %>% mutate(ACT_use= as.numeric(ACT)*100)
write_csv(pfK13_complete, file = "pfK13_SNP_complete.csv")
pfK13_SNP.1<- read_csv("pfK13_SNP_complete.csv")
write_csv(pfK13_SNP.1, file = "pfK13_SNP.1.csv")
Numstudies0 <- length(unique(pfK13_SNP.1$`ACT_use`))


#DUdata <- read_csv ("FeveredChildrenDrugTreatProp_ACT.csv")
#pfK13_complete <- right_join(DUdata,pfK13_SNP, by = c("country","year"))
#pfK13_complete <- pfK13_complete %>% mutate(ACT_use= as.numeric(ACT)*100)
#write_csv(pfK13_complete, file = "pfK13_SNP_complete.csv")
#pfK13_SNP.1<- read_csv("pfK13_SNP_complete.csv")
#write_csv(pfK13_SNP.1, file = "pfK13_SNP.1.csv")
#Numstudies0 <- length((pfK13_SNP.1$`ACT_use`))

#DUdata1 <- read_csv ("GLOBAL_DATAFLOW_1970-2022.csv")
#colnames(DUdata1)[1] <- "country"
#colnames(DUdata1)[4] <- "year"
#colnames(DUdata1)[5] <- "ACT_use"
#DUdata1 <- DUdata1 %>% filter(`Indicator`=="Malaria, first line treatment - percentage of febrile children (under age 5) receiving ACT (first line antimalarial drug), among those receiving any antimalarial drugs")
#DUdata1 <- DUdata1 %>% filter(`Sex`=="Total")
#View(DUdata1)
#Country_drug <- DUdata1[c(1,4,5)] %>% distinct(country, year, ACT_use) %>% arrange(country)
#Country_drug <- Country_drug %>% mutate(ACT_use= as.numeric(ACT_use)/100)
#write_csv(Country_drug , file = "Country_drug.csv")
#ACT_use <- right_join(Country_drug,DUdata, by = c("country","year"))
#ACT_use1 <- right_join(DUdata,Country_drug, by = c("country","year"))
#write_csv(ACT_use , file = "ACT_use.csv")
#write_csv(ACT_use1 , file = "ACT_use1.csv")
#pfK13_complete <- right_join(Country_drug,pfK13_SNP, by = c("country","year"))
#pfK13_complete <- right_join(DUdata,pfK13_SNP, by = c("country","year"))
#View(Country_drug)
#write_csv(pfK13_complete, file = "pfK13_SNP_complete.csv")
#pfK13_SNP.1<- read_csv("pfK13_SNP_complete.csv")
#write_csv(pfK13_SNP.1, file = "pfK13_SNP.1.csv")
#Numstudies0 <- length(unique(pfK13_SNP.1$`ACT_use`))


### Filtering out Data from 2000-2005
pcrt_SNP <- read_csv("pfcrt_SNP_complete.csv")
pcrt_SNP_year1a <- pcrt_SNP %>% filter (year%in%(2000:2005))
pcrt_SNP_year1b <- pcrt_SNP %>% filter (year%in%(2006:2021))
write_csv(pcrt_SNP_year1a, file = "pfcrt_SNP.1a.csv")
write_csv(pcrt_SNP_year1b, file = "pfcrt_SNP.1b.csv")
Numstudies0 <- length(unique(pcrt_SNP$`study Id`))
Numstudies1 <- length(unique(pcrt_SNP_year1a$`study Id`))
Numstudies2 <- length(unique(pcrt_SNP_year1b$`study Id`))
Numstudies3 <- length(unique(pcrt_SNP_year1a$`PR`))
Numstudies4 <- length(unique(pcrt_SNP_year1a$`Prev_DR`))
Numstudies5 <- length(unique(pcrt_SNP_year1a$`CQ_use`))
Numstudies3 <- length(unique(pcrt_SNP_year1b$`PR`))
Numstudies4 <- length(unique(pcrt_SNP_year1b$`Prev_DR`))
Numstudies5 <- length(unique(pcrt_SNP_year1b$`CQ_use`))


pfmdr1_SNP <- read_csv("pfmdr1_SNP_complete.csv")
pfmdr1_SNP_year1a <- pfmdr1_SNP %>% filter (year%in%(2000:2005))
pfmdr1_SNP_year1b <- pfmdr1_SNP %>% filter (year%in%(2006:2021))
write_csv(pfmdr1_SNP_year1a, file = "pfmdr1_SNP.1a.csv")
write_csv(pfmdr1_SNP_year1b, file = "pfmdr1_SNP.1b.csv")
view(pfmdr1_SNP_year1a)
Numstudies0 <- length(unique(pfmdr1_SNP$`study Id`))
Numstudies1 <- length(unique(pfmdr1_SNP_year1a$`study Id`))
Numstudies2 <- length(unique(pfmdr1_SNP_year1b$`study Id`))
Numstudies3 <- length(unique(pfmdr1_SNP_year1a$`PR`))
Numstudies4 <- length(unique(pfmdr1_SNP_year1a$`Prev_DR`))
Numstudies5 <- length(unique(pfmdr1_SNP_year1a$`CQ_use`))
Numstudies3 <- length(unique(pfmdr1_SNP_year1b$`PR`))
Numstudies4 <- length(unique(pfmdr1_SNP_year1b$`Prev_DR`))
Numstudies5 <- length(unique(pfmdr1_SNP_year1b$`CQ_use`))


#######
#pfdhps_SNP.1<- read_csv("pfdhps_SNP.csv")
#write_csv(pfdhps_SNP.1, file = "pfdhps_SNP.1.csv")
#Numstudies0 <- length(unique(pfdhps_SNP.1$`study id`))

#pfdhfr_SNP.1<- read_csv("pfdhfr_SNP.csv")
#write_csv(pfdhfr_SNP.1, file = "pfdhfr_SNP.1.csv")
#Numstudies0 <- length(unique(pfdhfr_SNP.1$`study id`))


#######
#pfcrt_p<- read.csv("pfcrt_SNP.1a.csv")
#pfcrt_p <- pfcrt_p %>% mutate(Prev_DR=(present/tested)*100)
#write_csv(pfcrt_p, file = "pfcrt_SNP.1a1.csv")

#pfcrt_p<- read.csv("pfcrt_SNP.1b.csv")
#pfcrt_p <- pfcrt_p %>% mutate(Prev_DR=(present/tested)*100)
#write_csv(pfcrt_p, file = "pfcrt_SNP.1b1.csv")

#pfmdr1_p<- read.csv("pfmdr1_SNP.1a.csv")
#pfmdr1_p <- pfmdr1_p %>% mutate(Prev_DR=(present/tested)*100)
#write_csv(pfmdr1_p, file = "pfmdr1_SNP.1a1.csv")

#pfmdr1_prop<- read.csv("pfmdr1_SNP.1b.csv")
#pfmdr1_prop <- pfmdr1_prop %>% mutate(Prev_DR=(present/tested)*100)
#write_csv(pfmdr1_prop, file = "pfmdr1_SNP.1b1.csv")

#pfdhps_p<- read.csv("pfdhps_SNP.1.csv")
#pfdhps_p<- pfdhps_p %>% mutate(Prev_DR=(present/tested)*100)
#write_csv(pfdhps_p, file = "pfdhps_SNP.1I.csv")


#pfdhfr_p<- read.csv("pfdhfr_SNP.1.csv")
#pfdhfr_p<- pfdhfr_p %>% mutate(Prev_DR=(present/tested)*100)
#write_csv(pfdhfr_p, file = "pfdhfr_SNP.1I.csv")
#summary(pfdhfr_p)




#####
#####
#####
###Plots
#1 Plots for pfcrt_76T

#pfcrt_SNP1a<- read.csv("pfcrt_SNP.1a1.csv")
pfcrt_SNP1a<- read.csv("pfcrt_SNP.1a.csv")
view (pfcrt_SNP1a)

### Pfcrt_76T Prevalence vs Year_All regions
ggplot(data=pfcrt_SNP1a%>% filter(Region %in% c("West Africa","Central Africa", "East Africa", "Southern Africa", "SEA", "South America")))+ aes(x = year, y = Prev_DR) + geom_point(aes(color=as.factor(Region), size=tested))+ geom_smooth(method = "lm") + xlab("Year") + ylab("Pfcrt 76T Prevalence") 
ggsave("pfcrt_76T_p1.png")

###Pfcrt_76T Prevalence vs PR_All regions
ggplot(data=pfcrt_SNP1a %>% filter(Region %in% c("West Africa","Central Africa", "East Africa", "Southern Africa", "SEA", "South America")))+ aes(x = PR, y = Prev_DR) + geom_point(aes(color=as.factor(Region), size=tested))+ geom_smooth(method = "lm") + xlab("Parasite Rate") + ylab("Pfcrt 76T Prevalence")
#ggplot(data=pfcrt_SNP1a %>% filter(Region %in% c("West Africa","Central Africa", "East Africa", "Southern Africa", "SEA", "South America")))+ aes(x = PR, y = Prev_DR) + geom_point(aes(color=as.factor(Region), size=tested))+ geom_smooth(method = "lm") + xlab("Parasite Rate") + ylab("Pfcrt 76T Prevalence") + scale_color_brewer(palette = "YlOrBr")
#ggplot(data=pfcrt_SNP1a %>% filter(Region %in% c("West Africa","Central Africa", "East Africa", "Southern Africa", "SEA", "South America")))+ aes(x = PR, y = Prev_DR) + geom_point(aes(color=as.factor(Region), size=tested))+ geom_smooth(method = "lm") + xlab("Parasite Rate") + ylab("Pfcrt 76T Prevalence") + scale_color_brewer(name = "as.factor(Region)", palette ="YlOrBr")
ggsave("pfcrt_76T_p2.png")

### Plot of the relationship between Region and Year varied between pfcrt76T Prevalence.
#pfcrt_76T_I<-pfcrt_SNP1a%>%mutate(Prev_DRCat = cut(Prev_DR,breaks=10*(0:10),include.lowest=T))
#ggplot(data = pfcrt_76T_I, aes (y = Region, x = year,colour = Prev_DRCat)) + geom_point(aes(size = 3)) + scale_color_viridis_d(name="Pfcrt 76T Prevalence") + theme_classic() + xlab ("Year") + ylab ("Region")
#ggplot(data = pfcrt_76T_I, aes (y = Region, x = year,colour = Prev_DRCat)) + geom_point() + scale_color_viridis_d(name="Pfcrt 76T Prevalence") + theme_classic() + xlab ("Year") + ylab ("Region")

ggplot(data = pfcrt_SNP1a %>% group_by(Region,year) %>% mutate(mean=mean(Prev_DR)), aes (y = Region, x = year,colour = mean)) + geom_point(aes(size=tested)) + scale_color_viridis_c(name="Pfcrt 76T Prevalence") + theme_classic() + xlab ("Year") + ylab ("Region")
#ggplot(data = pfcrt_SNP1a, aes (y = Region, x = year,colour = Prev_DR)) + geom_point() + scale_color_viridis_c(name="Pfcrt 76T Prevalence") + theme_classic() + xlab ("Year") + ylab ("Region")
#ggplot(data = pfcrt_76T_I, aes (y = Region, x = year,colour = Prev_DRCat)) + geom_point(aes(size = 0.1)) + scale_color_brewer(name = "Pfcrt 76T Prevalence", palette ="YlOrRd") + theme_classic() + xlab ("Year") + ylab ("Region")
ggsave("pfcrt 76T_Timep1.png")

# Plot of the relationship between Region and Year varied between Parasite rate.
pfcrt_76T_II<-pfcrt_SNP1a%>%mutate(PRCat = cut(PR,breaks=10*(0:10),include.lowest=T))
ggplot(data = pfcrt_76T_II, aes (y = Region, x = year,colour = PRCat)) + geom_point(size = 3) + scale_color_viridis_d(name="PR") + theme_classic() + xlab ("Year") + ylab ("Region")
ggplot(data = pfcrt_76T_II, aes (y = Region, x = year,colour = PRCat)) + geom_point(size = 3) + scale_color_viridis_d(name="PR") + theme_classic() + xlab ("Year") + ylab ("Region")

ggplot(data = pfcrt_SNP1a %>% group_by(Region,year) %>% mutate(mean=mean(PR)), aes (y = Region, x = year,colour = mean)) + geom_point(aes(size=tested)) + scale_color_viridis_c(name="PR") + theme_classic() + xlab ("Year") + ylab ("Region")
ggsave("pfcrt 76T_Timep2.png")


####
####
#pfmdr1_SNP1a<- read.csv("pfmdr1_SNP.1a1.csv")
pfmdr1_SNP1a<- read.csv("pfmdr1_SNP.1a.csv")

#### Pfmdr1_86Y Prevalence vs Year_All regions
ggplot(data=pfmdr1_SNP1a %>% filter(Region %in% c("West Africa","Central Africa", "East Africa", "Southern Africa", "SEA", "South America")))+ aes(x = year, y = Prev_DR) + geom_point(aes(color=as.factor(Region), size=tested))+ geom_smooth(method = "lm") + xlab("Year") + ylab("Pfmdr1 86Y Prevalence") 
ggsave("pfmdr1_86Y_p1.png")

###Pfmdr1_86Y Prevalence vs PR_All regions
ggplot(data=pfmdr1_SNP1a %>% filter(Region %in% c("West Africa","Central Africa", "East Africa", "Southern Africa", "SEA", "South America")))+ aes(x = PR, y = Prev_DR) + geom_point(aes(color=as.factor(Region), size=tested))+ geom_smooth(method = "lm") + xlab("Parasite Rate") + ylab("Pfmdr1 86Y Prevalence")
ggsave("pfmdr1_86Y_p2.png")

### Plot of the relationship between Region and Year varied between pfmdr1 86Y Prevalence.
pfmdr1_86T_I<-pfmdr1_SNP1a%>%mutate(Prev_DRCat = cut(Prev_DR,breaks=10*(0:10),include.lowest=T))
ggplot(data = pfmdr1_86T_I, aes (y = Region, x = year,colour = Prev_DRCat)) + geom_point(aes(size = 3)) + scale_color_viridis_d(name="Pfmdr1 86Y Prevalence") + theme_classic() + xlab ("Year") + ylab ("Region")

ggplot(data = pfmdr1_SNP1a %>% group_by(Region,year) %>% mutate(mean=mean(Prev_DR)), aes (y = Region, x = year,colour = mean)) + geom_point(aes(size=tested)) + scale_color_viridis_c(name="Pfmdr1 86Y Prevalence") + theme_classic() + xlab ("Year") + ylab ("Region")
ggsave("pfmdr1 86Y_Timep1.png")

# Plot of the relationship between Region and Year varied between Parasite rate.
pfmdr1_86T_II<-pfmdr1_SNP1a%>%mutate(PRCat = cut(PR,breaks=10*(0:10),include.lowest=T))
ggplot(data = pfmdr1_86T_II, aes (y = Region, x = year,colour = PRCat)) + geom_point(size = 3) + scale_color_viridis_d(name="PR") + theme_classic() + xlab ("Year") + ylab ("Region")

ggplot(data = pfmdr1_SNP1a %>% group_by(Region,year) %>% mutate(mean=mean(PR)), aes (y = Region, x = year,colour = mean)) + geom_point(aes(size=tested)) + scale_color_viridis_c(name="PR") + theme_classic() + xlab ("Year") + ylab ("Region")
ggsave("pfmdr1 86Y_Timep2.png")

####
####

#pfdhps_SNP1I<- read.csv("pfdhps_SNP.1I.csv")
pfdhps_SNP1I<- read.csv("pfdhps_SNP.csv")

#### Pfdhps_437G Prevalence vs Year_All regions
ggplot(data=pfdhps_SNP1I %>% filter(year %in% (2000:2021), Region %in% c("West Africa","Central Africa", "East Africa", "Southern Africa", "SEA", "South America")))+ aes(x = year, y = Prev_DR) + geom_point(aes(color=as.factor(Region), size=tested))+ geom_smooth(method = "lm") + xlab("Year") + ylab("Pfdhps 437G Prevalence") 
ggsave("pfdhps_437G_p1.png")

###Pfdhps_437G Prevalence vs PR_All regions
ggplot(data=pfdhps_SNP1I %>% filter(Region %in% c("West Africa","Central Africa", "East Africa", "Southern Africa", "SEA", "South America")))+ aes(x = PR, y = Prev_DR) + geom_point(aes(color=as.factor(Region), size=tested))+ geom_smooth(method = "lm") + xlab("Parasite Rate") + ylab("Pfdhps 437G Prevalence")
ggsave("pfdhps_437G_p2.png")

### Plot of the relationship between Region and Year varied between pfdhps_437G Prevalence.
pfdhps_437G_I<-pfdhps_SNP1I %>%mutate(Prev_DRCat = cut(Prev_DR,breaks=10*(0:10),include.lowest=T))
ggplot(data = pfdhps_437G_I, aes (y = Region, x = year,colour = Prev_DRCat)) + geom_point(aes(size = 3)) + scale_color_viridis_d(name="Pfdhps 437G Prevalence") + theme_classic() + xlab ("Year") + ylab ("Region")

ggplot(data = pfdhps_SNP1I %>% filter(year %in% (2000:2021)) %>%  group_by(Region,year) %>% mutate(mean=mean(Prev_DR)), aes (y = Region, x = year,colour = mean)) + geom_point(aes(size=tested)) + scale_color_viridis_c(name="Pfdhps 437G Prevalence") + theme_classic() + xlab ("Year") + ylab ("Region")
ggsave("pfdhps 437G_Timep1.png")

# Plot of the relationship between Region and Year varied between Parasite rate.
pfdhps_437G_II<-pfdhps_SNP1I%>%mutate(PRCat = cut(PR,breaks=10*(0:10),include.lowest=T))
ggplot(data = pfdhps_437G_II, aes (y = Region, x = year,colour = PRCat)) + geom_point(size = 3) + scale_color_viridis_d(name="PR") + theme_classic() + xlab ("Year") + ylab ("Region")

ggplot(data = pfdhps_SNP1I %>%  filter(year %in% (2000:2021)) %>%  group_by(Region,year) %>% mutate(mean=mean(PR)), aes (y = Region, x = year,colour = mean)) + geom_point(aes(size=tested)) + scale_color_viridis_c(name="PR") + theme_classic() + xlab ("Year") + ylab ("Region")
ggsave("pfdhps 437G_Timep2.png")

####
####

#pfdhfr_SNP1I<- read.csv("pfdhfr_SNP.1I.csv")
pfdhfr_SNP1I<- read.csv("pfdhfr_SNP.csv")

#### Pfdhfr_108N Prevalence vs Year_All regions
ggplot(data=pfdhfr_SNP1I %>% filter(year %in% (2000:2021), Region %in% c("West Africa","Central Africa", "East Africa", "Southern Africa", "SEA", "South America")))+ aes(x = year, y = Prev_DR) + geom_point(aes(color=as.factor(Region), size=tested))+ geom_smooth(method = "lm") + xlab("Year") + ylab("Pfdhfr 108N Prevalence") 
ggsave("pfdhfr_108N_p1.png")

###Pfdhfr_108N Prevalence vs PR_All regions
ggplot(data=pfdhfr_SNP1I %>% filter(Region %in% c("West Africa","Central Africa", "East Africa", "Southern Africa", "SEA", "South America")))+ aes(x = PR, y = Prev_DR) + geom_point(aes(color=as.factor(Region), size=tested))+ geom_smooth(method = "lm") + xlab("Parasite Rate") + ylab("Pfdhfr 108N Prevalence")
ggsave("pfdhfr_108N_p2.png")

### Plot of the relationship between Region and Year varied between pfdhfr_108N Prevalence.
pfdhfr_108N_I<-pfdhfr_SNP1I %>%mutate(Prev_DRCat = cut(Prev_DR,breaks=10*(0:10),include.lowest=T))
ggplot(data = pfdhfr_108N_I, aes (y = Region, x = year,colour = Prev_DRCat)) + geom_point(aes(size = 3)) + scale_color_viridis_d(name="Pfdhfr 108N Prevalence") + theme_classic() + xlab ("Year") + ylab ("Region")

ggplot(data = pfdhfr_SNP1I %>% filter(year %in% (2000:2021)) %>% group_by(Region,year) %>% mutate(mean=mean(Prev_DR)), aes (y = Region, x = year,colour = mean)) + geom_point(aes(size=tested)) + scale_color_viridis_c(name="Pfdhfr 108N Prevalence") + theme_classic() + xlab ("Year") + ylab ("Region")

ggsave("pfdhfr 108N_Timep1.png")

# Plot of the relationship between Region and Year varied between Parasite rate.
pfdhfr_108N_II<-pfdhfr_SNP1I %>%mutate(PRCat = cut(PR,breaks=10*(0:10),include.lowest=T))
ggplot(data = pfdhfr_108N_II, aes (y = Region, x = year,colour = PRCat)) + geom_point(size = 3) + scale_color_viridis_d(name="PR") + theme_classic() + xlab ("Year") + ylab ("Region")

ggplot(data = pfdhfr_SNP1I %>% filter(year %in% (2000:2021)) %>% group_by(Region,year) %>% mutate(mean=mean(PR)), aes (y = Region, x = year,colour = mean)) + geom_point(aes(size=tested)) + scale_color_viridis_c(name="PR") + theme_classic() + xlab ("Year") + ylab ("Region")
ggsave("pfdhfr 108N_Timep2.png")


######
######
#Models for all sample sizes

###pfcrt-76T

pfcrt_SNP1a<- read.csv("pfcrt_SNP.1a.csv")
summary(pfcrt_SNP1a)
view((pfcrt_SNP1a))
#pfcrt_SNP1a$CQ_use[is.na(pfcrt_SNP1a$CQ_use)] <- 0

pfcrt_SNP1<- read.csv("pfcrt_SNP_complete.csv")
summary(pfcrt_SNP1)
view((pfcrt_SNP1))

# assign factors
pfcrt_SNP1$study.Id <- as.factor(pfcrt_SNP1$study.Id)
pfcrt_SNP1$country <- as.factor(pfcrt_SNP1$country)
#a$x[is.na(a$x)] <- 0
#pfcrt_SNP1$CQ_use <- as.factor(pfcrt_SNP1$CQ_use)

Model <- glm(Prev_DR ~ year + PR + CQ_use, data = pfcrt_SNP1, family = gaussian)
summary(Model)

Meta0.pfcrt_SNP1 <- lmer(Prev_DR ~ year + PR + CQ_use + (1|Region/study.Id), data=pfcrt_SNP1)
summary(Meta0.pfcrt_SNP1)
plot(Meta0.pfcrt_SNP1)

pfcrt_SNP1 <- pfcrt_SNP1 %>% mutate(PR = PR / 100, Prev_DR = Prev_DR / 100, CQ_use = CQ_use / 100)
Meta0.pfcrt_SNP1 <- glmer(Prev_DR ~  year + PR + CQ_use + (1|Region/study.Id), data=pfcrt_SNP1, family = binomial(link = "logit"))
summary(Meta0.pfcrt_SNP1)
plot(Meta0.pfcrt_SNP1)

Meta1.pfcrt_SNP1 <- lmer(Prev_DR ~  year + PR + CQ_use + (1|study.Id), data=pfcrt_SNP1)
summary(Meta1.pfcrt_SNP1)
plot(Meta1.pfcrt_SNP1)

Meta1.pfcrt_SNP1 <- glmer(Prev_DR ~ year + PR + CQ_use + (1|study.Id), data=pfcrt_SNP1, family = binomial(link = "logit"))
summary(Meta1.pfcrt_SNP1)
plot(Meta1.pfcrt_SNP1)


?logit
library(boot)
pfcrt_SNP1_sub <-pfcrt_SNP1 %>% filter(PR<100,PR>0,Prev_DR<100,Prev_DR>0, CQ_use<100,CQ_use>0)
Meta2.pfcrt_SNP1 <- lmer(logit(Prev_DR/100) ~ year + logit(PR/100) + logit(CQ_use/100) + (1|country/study.Id), weights=tested, data=pfcrt_SNP1_sub)
summary(Meta2.pfcrt_SNP1)
plot(Meta2.pfcrt_SNP1)

pfcrt_SNP1_sub <-pfcrt_SNP1 %>% filter(PR<100,PR>0,Prev_DR<100,Prev_DR>0, CQ_use<100,CQ_use>0)
Meta3.pfcrt_SNP1 <- lmer(logit(Prev_DR/100) ~ year + logit(PR/100) + logit(CQ_use/100) + (1|study.Id), weights=tested, data=pfcrt_SNP1_sub)
summary(Meta3.pfcrt_SNP1)
plot(Meta3.pfcrt_SNP1a)

###pfmdr1-86Y

pfmdr1_SNP1<- read.csv("pfmdr1_SNP_complete.csv")
summary(pfmdr1_SNP1)
# assign factors
pfmdr1_SNP1$study.Id <- as.factor(pfmdr1_SNP1$study.Id)
pfmdr1_SNP1$country <- as.factor(pfmdr1_SNP1$country)

Meta0.pfmdr1_SNP1 <- lmer(Prev_DR ~ year + PR + CQ_use + (1|Region/study.Id), data=pfmdr1_SNP1)
summary(Meta0.pfmdr1_SNP1)
plot(Meta0.pfmdr1_SNP1)

Meta1.pfmdr1_SNP1 <- lmer(Prev_DR ~ year + PR + CQ_use + (1|study.Id), data=pfmdr1_SNP1)
summary(Meta1.pfmdr1_SNP1)
plot(Meta1.pfmdr1_SNP1)

library(boot)
pfmdr1_SNP1_sub <-pfmdr1_SNP1 %>% filter(PR<100,PR>0,Prev_DR<100,Prev_DR>0, CQ_use<100,CQ_use>0)
Meta2.pfmdr1_SNP1 <- lmer(logit(Prev_DR/100) ~ year + logit(PR/100) + logit(CQ_use/100) + (1|country/study.Id), weights=tested, data=pfmdr1_SNP1_sub)
summary(Meta2.pfmdr1_SNP1)
plot(Meta2.pfmdr1_SNP1)

pfmdr1_SNP1_sub <-pfmdr1_SNP1 %>% filter(PR<100,PR>0,Prev_DR<100,Prev_DR>0,CQ_use<100,CQ_use>0)
Meta3.pfmdr1_SNP1 <- lmer(logit(Prev_DR/100) ~ year + logit(PR/100) + logit(CQ_use/100) + (1|study.Id), weights=tested, data=pfmdr1_SNP1_sub)
summary(Meta3.pfmdr1_SNP1)
plot(Meta3.pfmdr1_SNP1)


###pfdhps-437G

pfdhps_SNP1<- read.csv("pfdhps_SNP.1.csv")
# assign factors
pfdhps_SNP1$study.id <- as.factor(pfdhps_SNP1$study.id)
pfdhps_SNP1$country <- as.factor(pfdhps_SNP1$country)

Meta0.pfdhps_SNP1 <- lmer(Prev_DR ~ year + PR + SP_use + (1|Region/study.id), data=pfdhps_SNP1)
summary(Meta0.pfdhps_SNP1)
plot(Meta0.pfdhps_SNP1)
anova(Meta0.pfdhps_SNP1, ddf= "Kenward-Roger")

Meta1.pfdhps_SNP1 <- lmer(Prev_DR ~ year + PR + SP_use + (1|study.id), data=pfdhps_SNP1)
summary(Meta1.pfdhps_SNP1)
plot(Meta1.pfdhps_SNP1)


library(boot)
pfdhps_SNP1_sub <-pfdhps_SNP1 %>% filter(PR<100,PR>0,Prev_DR<100,Prev_DR>0, SP_use<100,SP_use>0)
Meta2.pfdhps_SNP1 <- lmer(logit(Prev_DR/100) ~ year + logit(PR/100) + logit(SP_use/100) + (1|country/study.id), weights=tested, data=pfdhps_SNP1_sub)
summary(Meta2.pfdhps_SNP1)
plot(Meta2.pfdhps_SNP1)

pfdhps_SNP1_sub <-pfdhps_SNP1 %>% filter(PR<100,PR>0,Prev_DR<100,Prev_DR>0, SP_use<100,SP_use>0)
Meta3.pfdhps_SNP1 <- lmer(logit(Prev_DR/100) ~ year + logit(PR/100) + logit(SP_use/100) + (1|study.id), weights=tested, data=pfdhps_SNP1_sub)
summary(Meta3.pfdhps_SNP1)
plot(Meta3.pfdhps_SNP1)

###pfdhfr-108N

pfdhfr_SNP1<- read.csv("pfdhfr_SNP.1.csv")
# assign factors
pfdhfr_SNP1$study.id <- as.factor(pfdhfr_SNP1$study.id)
pfdhfr_SNP1$country <- as.factor(pfdhfr_SNP1$country)

Meta0.pfdhfr_SNP1 <- lmer(Prev_DR ~ year + PR + SP_use + (1|Region/study.id), data=pfdhfr_SNP1)
summary(Meta0.pfdhfr_SNP1)
plot(Meta0.pfdhfr_SNP1)
anova(Meta0.pfdhfr_SNP1, ddf= "Kenward-Roger")

Meta1.pfdhfr_SNP1 <- lmer(Prev_DR ~ year + PR + SP_use + (1|study.id), data=pfdhfr_SNP1)
summary(Meta1.pfdhfr_SNP1)
plot(Meta1.pfdhfr_SNP1)

library(boot)
pfdhfr_SNP1_sub <-pfdhfr_SNP1 %>% filter(PR<100,PR>0,Prev_DR<100,Prev_DR>0, SP_use<100,SP_use>0)
Meta2.pfdhfr_SNP1 <- lmer(logit(Prev_DR/100) ~ year + logit(PR/100) + logit(SP_use/100) + (1|country/study.id), weights=tested, data=pfdhfr_SNP1_sub)
summary(Meta2.pfdhfr_SNP1)
plot(Meta2.pfdhfr_SNP1)

pfdhfr_SNP1_sub <-pfdhfr_SNP1 %>% filter(PR<100,PR>0,Prev_DR<100,Prev_DR>0, SP_use<100,SP_use>0)
Meta3.pfdhfr_SNP1 <- lmer(logit(Prev_DR/100) ~ year + logit(PR/100) + logit(SP_use/100) + (1|study.id), weights=tested, data=pfdhfr_SNP1_sub)
summary(Meta3.pfdhfr_SNP1)
plot(Meta3.pfdhfr_SNP1)


###pfK13-580Y

pfK13_SNP1<- read.csv("pfK13_SNP.1.csv")
# assign factors
pfK13_SNP1$sid <- as.factor(pfK13_SNP1$sid)
pfK13_SNP1$country <- as.factor(pfK13_SNP1$country)

Meta0.pfK13_SNP1 <- lmer(Prev_DR ~ year + PR + ACT_use + (1|Region/sid), data=pfK13_SNP1)
summary(Meta0.pfK13_SNP1)
plot(Meta0.pfK13_SNP1)
anova(Meta0.pfK13_SNP1, ddf= "Kenward-Roger")

Meta1.pfK13_SNP1 <- lmer(Prev_DR ~ year + PR + ACT_use + (1|sid), data=pfK13_SNP1)
summary(Meta1.pfK13_SNP1)
plot(Meta1.pfK13_SNP1)

library(boot)
pfK13_SNP1_sub <-pfK13_SNP1 %>% filter(PR<100,PR>0,Prev_DR<100,Prev_DR>0, ACT_use<100,ACT_use>0)
Meta2.pfK13_SNP1 <- lmer(logit(Prev_DR/100) ~ year + logit(PR/100) + logit(ACT_use/100) + (1|country/sid), weights=tested, data=pfK13_SNP1_sub)
summary(Meta2.pfK13_SNP1)
plot(Meta2.pfK13_SNP1)

pfK13_SNP1_sub <-pfK13_SNP1 %>% filter(PR<100,PR>0,Prev_DR<100,Prev_DR>0, ACT_use<100,ACT_use>0)
Meta3.pfK13_SNP1 <- lmer(logit(Prev_DR/100) ~ year + logit(PR/100) + logit(ACT_use/100) + (1|sid), weights=tested, data=pfK13_SNP1_sub)
summary(Meta3.pfK13_SNP1)
plot(Meta3.pfK13_SNP1)



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


library(boot)
pfcrt_SNP1a3_sub <-pfcrt_SNP.1a3 %>% filter(PR<100,PR>0,Prev_DR<100,Prev_DR>0)
Meta3.pfcrt_SNP1a3 <- lmer(logit(Prev_DR/100) ~ year + logit(PR/100) + (1|study.Id), weights=tested, data=pfcrt_SNP1a3_sub)
summary(Meta3.pfcrt_SNP1a3)
plot(Meta3.pfcrt_SNP1a3)


pfmdr1_86Y_new_data_1<- read_csv("pfmdr1_SNP.1a1.csv")
pfmdr1_86Y_new_data_2 <- pfmdr1_86Y_new_data_1 %>% filter(tested>50)
write_csv(pfmdr1_86Y_new_data_2, file = "pfmdr1_SNP.1a3.csv")

pfmdr1_SNP.1a3<- read_csv("pfmdr1_SNP.1a3.csv")
Meta1.pfmdr1_SNP1a3 <- lmer(Prev_DR ~ year + PR + (1|country/study.Id), data=pfmdr1_SNP.1a3)
summary(Meta1.pfmdr1_SNP1a3)
Meta2.pfmdr1_SNP1a3 <- lmer(Prev_DR ~ year + PR + (1|study.Id), data=pfmdr1_SNP.1a3)
summary(Meta2.pfmdr1_SNP1a3)

library(boot)
pfmdr1_SNP1a3_sub <-pfmdr1_SNP.1a3 %>% filter(PR<100,PR>0,Prev_DR<100,Prev_DR>0)
Meta3.pfmdr1_SNP1a3 <- lmer(logit(Prev_DR/100) ~ year + logit(PR/100) + (1|study.Id), weights=tested, data=pfmdr1_SNP1a3_sub)
summary(Meta3.pfmdr1_SNP1a3)
plot(Meta3.pfmdr1_SNP1a3)


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

library(boot)
pfdhfr_SNP1_sub <-pfdhfr_SNP1 %>% filter(PR<100,PR>0,Prev_DR<100,Prev_DR>0)
Meta3.pfdhfr_SNP1 <- lmer(logit(Prev_DR/100) ~ year + logit(PR/100) + (1|study.id), weights=tested, data=pfdhfr_SNP1_sub)
summary(Meta3.pfdhfr_SNP1)
plot(Meta3.pfdhfr_SNP1)
