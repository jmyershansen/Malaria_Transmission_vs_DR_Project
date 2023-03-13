# DATE: 03/10/2023

library(tidyverse)
library(dplyr)
library (lme4)
library(lme4)
library(lattice)
library(ggplot2)
library(gridExtra)

setwd("D:/DATA/Purdue University/Purdue_Fall 2022/BIO 69500 (Data Science for Biologists)/Project/BIO 69500_Project/CQ Markers")


#1 Plots for pfcrt_76T

pfcrt_76T<- read_csv("pfcrt_SNP.3a1.csv")

### Pfcrt_76T Prevalence vs Year_All regions
ggplot(data=pfcrt_76T %>% filter(Region %in% c("West Africa","Central Africa", "East Africa", "Southern Africa", "SEA", "South America")))+ aes(x = year, y = Prev_DR) + geom_point(aes(color=as.factor(Region)))+ geom_smooth(method = "lm") + xlab("Year") + ylab("Pfcrt 76T Prevalence") 
ggsave("pfcrt_76T_p1.png")

# Pfcrt_76T Prevalence vs Year_Facetted for Africa
ggplot(data=pfcrt_76T %>% filter(Region %in% c("West Africa","Central Africa", "East Africa", "Southern Africa"))) + aes(x = year, y = Prev_DR) + geom_point() + facet_wrap(~Region) + geom_smooth(method = "lm")+ xlab("Year") + ylab("Pfcrt 76T Prevalence")
ggsave("pfcrt_76T_p2.png")


###Pfcrt_76T Prevalence vs PR_All regions
ggplot(data=pfcrt_76T %>% filter(Region %in% c("West Africa","Central Africa", "East Africa", "Southern Africa", "SEA", "South America")))+ aes(x = PR, y = Prev_DR) + geom_point(aes(color=as.factor(Region)))+ geom_smooth(method = "lm") + xlab("Parasite Rate") + ylab("Pfcrt 76T Prevalence")
ggsave("pfcrt_76T_p3.png")

# Pfcrt_76T Prevalence vs PR_Facetted for Africa
ggplot(data=pfcrt_76T  %>% filter(Region %in% c("West Africa","Central Africa", "East Africa", "Southern Africa"))) + aes(x = PR, y = Prev_DR) + geom_point() + facet_wrap(~Region) +  geom_smooth(method = "lm")+ xlab("Parasite Rate") + ylab("Pfcrt 76T Prevalence")
ggsave("pfcrt_76T_p4.png")


### Plot of the relationship between Region and Year varied between pfcrt76T Prevalence.
pfcrt_76T_I<-pfcrt_76T%>%mutate(Prev_DRCat = cut(Prev_DR,breaks=10*(0:10),include.lowest=T))
ggplot(data = pfcrt_76T_I, aes (y = Region, x = year,colour = Prev_DRCat)) + geom_point(size = 3) + scale_color_viridis_d(name="Pfcrt 76T Prevalence") + theme_classic() + xlab ("Year") + ylab ("Region")
ggsave("pfcrt 76T_Timep1.png")

# Plot of the relationship between Region and Year varied between Parasite rate.
pfcrt_76T_II<-pfcrt_76T%>%mutate(PRCat = cut(PR,breaks=10*(0:10),include.lowest=T))
ggplot(data = pfcrt_76T_II, aes (y = Region, x = year,colour = PRCat)) + geom_point(size = 3) + scale_color_viridis_d(name="PR") + theme_classic() + xlab ("Year") + ylab ("Region")
#ggplot(data = pfcrt_pm3_region, aes (y = Region, x = Year,colour = PRCat)) + geom_point(size = 3) + scale_color_viridis_d(name="PR") + theme_classic() + xlim(2010,2020) + xlab ("Year") + ylab ("Region")
ggsave("pfcrt 76T_Timep2.png")

###Caterpillar plot of random effects
randoms1 <- ranef(Meta1.pfcrt_SNP3a)
dotplot(randoms1,scales=list(cex=0.7))
#OR
qqmath(randoms)
ggsave("pfcrt 76T_Randomp.png")


#2 Plots for pfmdr1_86Y

pfmdr1_86Y<- read_csv("pfmdr1_SNP.3a1.csv")

### Pfmdr1_86Y Prevalence vs Year_All regions
ggplot(data=pfmdr1_86Y %>% filter(Region %in% c("West Africa","Central Africa", "East Africa", "Southern Africa", "SEA", "South America")))+ aes(x = year, y = Prev_DR) + geom_point(aes(color=as.factor(Region)))+ geom_smooth(method = "lm") + xlab("Year") + ylab("Pfmdr1 86Y Prevalence") 
ggsave("pfmdr1_86Y_p1.png")

# Pfmdr1_86Y Prevalence vs Year_Facetted for Africa
ggplot(data=pfmdr1_86Y %>% filter(Region %in% c("West Africa","Central Africa", "East Africa", "Southern Africa"))) + aes(x = year, y = Prev_DR) + geom_point() + facet_wrap(~Region) + geom_smooth(method = "lm")+ xlab("Year") + ylab("Pfmdr1 86Y Prevalence")
ggsave("pfmdr1_86Y_p2.png")


###Pfmdr1_86Y Prevalence vs PR_All regions
ggplot(data=pfmdr1_86Y %>% filter(Region %in% c("West Africa","Central Africa", "East Africa", "Southern Africa", "SEA", "South America")))+ aes(x = PR, y = Prev_DR) + geom_point(aes(color=as.factor(Region)))+ geom_smooth(method = "lm") + xlab("Parasite Rate") + ylab("Pfmdr1 86Y Prevalence")
ggsave("pfmdr1_86Y_p3.png")


# Pfmdr1_86Y Prevalence vs PR_Facetted for Africa
ggplot(data=pfmdr1_86Y   %>% filter(Region %in% c("West Africa","Central Africa", "East Africa", "Southern Africa"))) + aes(x = PR, y = Prev_DR) + geom_point() + facet_wrap(~Region) +  geom_smooth(method = "lm")+ xlab("Parasite Rate") + ylab("Pfmdr1 86Y  Prevalence")
ggsave("pfmdr1_86Y_p4.png")


### Plot of the relationship between Region and Year varied between Pfmdr186Y Prevalence.
pfmdr1_86Y_I<-pfmdr1_86Y%>%mutate(Prev_DRCat = cut(Prev_DR,breaks=10*(0:10),include.lowest=T))
ggplot(data = pfmdr1_86Y_I, aes (y = Region, x = year,colour = Prev_DRCat)) + geom_point(size = 3) + scale_color_viridis_d(name="Pfmdr1 86Y Prevalence") + theme_classic() + xlab ("Year") + ylab ("Region")
ggsave("pfmdr1 86Y_Timep1.png")

# Plot of the relationship between Region and Year varied between Parasite rate.
pfmdr1_86Y_II<-pfmdr1_86Y%>%mutate(PRCat = cut(PR,breaks=10*(0:10),include.lowest=T))
ggplot(data = pfmdr1_86Y_II, aes (y = Region, x = year,colour = PRCat)) + geom_point(size = 3) + scale_color_viridis_d(name="PR") + theme_classic() + xlab ("Year") + ylab ("Region")
#ggplot(data = pfcrt_pm3_region, aes (y = Region, x = Year,colour = PRCat)) + geom_point(size = 3) + scale_color_viridis_d(name="PR") + theme_classic() + xlim(2010,2020) + xlab ("Year") + ylab ("Region")
ggsave("pfmdr1 86Y_Timep2.png")

###Caterpillar plot of random effects
randoms2 <- ranef(Meta1.pfmdr1_SNP3a)
dotplot(randoms2,scales=list(cex=0.7))
#OR
qqmath(randoms)
ggsave("pfmdr1 86Y_Randomp.png")


#3 Plots for pfdhfr_108N

pfdhfr_108N<- read_csv("pfdhfr_SNP.3I.csv")

### Pfdhfr_108N Prevalence vs Year_All regions
ggplot(data=pfdhfr_108N %>% filter(Region %in% c("West Africa","Central Africa", "East Africa", "Southern Africa", "SEA", "South America")))+ aes(x = year, y = Prev_DR) + geom_point(aes(color=as.factor(Region)))+ geom_smooth(method = "lm") + xlab("Year") + ylab("Pfdhfr 108N Prevalence") 
ggsave("pfdhfr_108N_p1.png")

# Pfdhfr_108N Prevalence vs Year_Facetted for Africa
ggplot(data=pfdhfr_108N %>% filter(Region %in% c("West Africa","Central Africa", "East Africa", "Southern Africa"))) + aes(x = year, y = Prev_DR) + geom_point() + facet_wrap(~Region) + geom_smooth(method = "lm")+ xlab("Year") + ylab("Pfdhfr 108N Prevalence")
ggsave("pfdhfr_108N_p2.png")


###Pfdhfr_108N Prevalence vs PR_All regions
ggplot(data=pfdhfr_108N %>% filter(Region %in% c("West Africa","Central Africa", "East Africa", "Southern Africa", "SEA", "South America")))+ aes(x = PR, y = Prev_DR) + geom_point(aes(color=as.factor(Region)))+ geom_smooth(method = "lm") + xlab("Parasite Rate") + ylab("Pfdhfr 108N Prevalence")
ggsave("pfdhfr_108N_p3.png")


# Pfdhfr_108N Prevalence vs PR_Facetted for Africa
ggplot(data=pfdhfr_108N %>% filter(Region %in% c("West Africa","Central Africa", "East Africa", "Southern Africa"))) + aes(x = PR, y = Prev_DR) + geom_point() + facet_wrap(~Region) +  geom_smooth(method = "lm")+ xlab("Parasite Rate") + ylab("Pfdhfr 108N Prevalence")
ggsave("pfdhfr_108N_p4.png")


### Plot of the relationship between Region and Year varied between Pfdhfr108N  Prevalence.
pfdhfr_108N_I<-pfdhfr_108N%>%mutate(Prev_DRCat = cut(Prev_DR,breaks=10*(0:10),include.lowest=T))
ggplot(data = pfdhfr_108N_I, aes (y = Region, x = year,colour = Prev_DRCat)) + geom_point(size = 3) + scale_color_viridis_d(name="Pfdhfr 108N Prevalence") + theme_classic() + xlab ("Year") + ylab ("Region")
ggsave("pfdhfr 108N_Timep1.png")

# Plot of the relationship between Region and Year varied between Parasite rate.
pfdhfr_108N_II<-pfdhfr_108N%>%mutate(PRCat = cut(PR,breaks=10*(0:10),include.lowest=T))
ggplot(data = pfdhfr_108N_II, aes (y = Region, x = year,colour = PRCat)) + geom_point(size = 3) + scale_color_viridis_d(name="PR") + theme_classic() + xlab ("Year") + ylab ("Region")
#ggplot(data = pfcrt_pm3_region, aes (y = Region, x = Year,colour = PRCat)) + geom_point(size = 3) + scale_color_viridis_d(name="PR") + theme_classic() + xlim(2010,2020) + xlab ("Year") + ylab ("Region")
ggsave("pfdhfr 108N_Timep2.png")

###Caterpillar plot of random effects
randoms3 <- ranef(Meta1.pfdhfr_SNP3)
dotplot(randoms3,scales=list(cex=0.7))
#OR
qqmath(randoms)
ggsave("pfdhfr 108N_Randomp.png")


#4 Plots for pfdhps_437G

pfdhps_437G<- read_csv("pfdhps_SNP.3I.csv")

### Pfdhps_437G Prevalence vs Year_All regions
ggplot(data=pfdhps_437G %>% filter(Region %in% c("West Africa","Central Africa", "East Africa", "Southern Africa", "SEA", "South America")))+ aes(x = year, y = Prev_DR) + geom_point(aes(color=as.factor(Region)))+ geom_smooth(method = "lm") + xlab("Year") + ylab("Pfdhps 437G Prevalence") 
ggsave("pfdhps_437G_p1.png")

# Pfdhps_437G Prevalence vs Year_Facetted for Africa
ggplot(data=pfdhps_437G %>% filter(Region %in% c("West Africa","Central Africa", "East Africa", "Southern Africa"))) + aes(x = year, y = Prev_DR) + geom_point() + facet_wrap(~Region) + geom_smooth(method = "lm")+ xlab("Year") + ylab("Pfdhps 437G Prevalence")
ggsave("pfdhps_437G_p2.png")


###Pfdhps_437G Prevalence vs PR_All regions
ggplot(data=pfdhps_437G %>% filter(Region %in% c("West Africa","Central Africa", "East Africa", "Southern Africa", "SEA", "South America")))+ aes(x = PR, y = Prev_DR) + geom_point(aes(color=as.factor(Region)))+ geom_smooth(method = "lm") + xlab("Parasite Rate") + ylab("Pfdhps 437G  Prevalence")
ggsave("pfdhps_437G_p3.png")


# Pfdhps_437G Prevalence vs PR_Facetted for Africa
ggplot(data=pfdhps_437G %>% filter(Region %in% c("West Africa","Central Africa", "East Africa", "Southern Africa"))) + aes(x = PR, y = Prev_DR) + geom_point() + facet_wrap(~Region) +  geom_smooth(method = "lm")+ xlab("Parasite Rate") + ylab("Pfdhps 437G Prevalence")
ggsave("pfdhps_437G_p4.png")


### Plot of the relationship between Region and Year varied between Pfdhps437G Prevalence.
pfdhps_437G_I<-pfdhps_437G%>%mutate(Prev_DRCat = cut(Prev_DR,breaks=10*(0:10),include.lowest=T))
ggplot(data = pfdhps_437G_I, aes (y = Region, x = year,colour = Prev_DRCat)) + geom_point(size = 3) + scale_color_viridis_d(name="Pfdhps 437G Prevalence") + theme_classic() + xlab ("Year") + ylab ("Region")
ggsave("pfdhps 437G_Timep1.png")

# Plot of the relationship between Region and Year varied between Parasite rate.
pfdhps_437G_II<-pfdhps_437G%>%mutate(PRCat = cut(PR,breaks=10*(0:10),include.lowest=T))
ggplot(data = pfdhps_437G_II, aes (y = Region, x = year,colour = PRCat)) + geom_point(size = 3) + scale_color_viridis_d(name="PR") + theme_classic() + xlab ("Year") + ylab ("Region")
#ggplot(data = pfcrt_pm3_region, aes (y = Region, x = Year,colour = PRCat)) + geom_point(size = 3) + scale_color_viridis_d(name="PR") + theme_classic() + xlim(2010,2020) + xlab ("Year") + ylab ("Region")
ggsave("pfdhps 437G_Timep2.png")

###Caterpillar plot of random effects
randoms4 <- ranef(Meta1.pfdhps_SNP3)
dotplot(randoms4,scales=list(cex=0.7))
#OR
qqmath(randoms)
ggsave("pfdhps 437G_Randomp.png")


######
######
#Models

###pfcrt-76T

pfcrt_SNP3a<- read.csv("pfcrt_SNP.3a1.csv")
# assign factors
pfcrt_SNP3a$study.Id <- as.factor(pfcrt_SNP3a$study.Id)
pfcrt_SNP3a$country <- as.factor(pfcrt_SNP3a$country)

Meta1.pfcrt_SNP3a <- lmer(Prev_DR ~ year + PR + (1|country/study.Id), data=pfcrt_SNP3a)
summary(Meta1.pfcrt_SNP3a)


###pfmdr1-86Y

pfmdr1_SNP3a<- read.csv("pfmdr1_SNP.3a1.csv")
# assign factors
pfmdr1_SNP3a$study.Id <- as.factor(pfmdr1_SNP3a$study.Id)
pfmdr1_SNP3a$country <- as.factor(pfmdr1_SNP3a$country)

Meta1.pfmdr1_SNP3a <- lmer(Prev_DR ~ year + PR + (1|country/study.Id), data=pfmdr1_SNP3a)
summary(Meta1.pfmdr1_SNP3a)


###pfdhps-437G

pfdhps_SNP3<- read.csv("pfdhps_SNP.3I.csv")
# assign factors
pfdhps_SNP3$study.Id <- as.factor(pfdhps_SNP3$study.Id)
pfdhps_SNP3$country <- as.factor(pfdhps_SNP3$country)

Meta1.pfdhps_SNP3 <- lmer(Prev_DR ~ year + PR + (1|country/study.Id), data=pfdhps_SNP3)
summary(Meta1.pfdhps_SNP3)
anova(Meta1.pfdhps_SNP3, ddf= "Kenward-Roger")


###pfdhfr-108N

pfdhfr_SNP3<- read.csv("pfdhfr_SNP.3I.csv")
# assign factors
pfdhfr_SNP3$study.Id <- as.factor(pfdhfr_SNP3$study.Id)
pfdhfr_SNP3$country <- as.factor(pfdhfr_SNP3$country)

Meta1.pfdhfr_SNP3 <- lmer(Prev_DR ~ year + PR + (1|country/study.Id), data=pfdhfr_SNP3)
summary(Meta1.pfdhfr_SNP3)
anova(Meta1.pfdhfr_SNP3, ddf= "Kenward-Roger")










