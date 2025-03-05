
library(tidyverse)
library(dplyr)
library(lme4)
library(lmerTest)
library(ggplot2)

#Summary:
#Multilevel mixed regression models coupled with time lag framework was used to access the effect of key malaria indicators: Parasite rate, PR (transmission intensity) and Drug use (Proportion of fevered children treated with CQ, SP and ACTs) on the reported prevalence of key drug (CQ, SP and ACTs) resistant markers: pfcrt-K76T, pfmdr1-N86Y, pfdhps-A437G, pfdhfr-S108N and pfk13-C580Y. independently.


###PFCRT K76T
#Prepare Database

cq<- read_csv("CQ_imputed.csv")
#Set up lag years
cq<-cq%>%group_by(country)%>%mutate(drug_1yr_ago = lag(mean_cq,n=1),drug_2yr_ago = lag(mean_cq,n=2), drug_3yr_ago = lag(mean_cq,n=3),drug_4yr_ago = lag(mean_cq,n=4),drug_5yr_ago = lag(mean_cq,n=5))

# join original observation
crt<-read_csv ("pfcrt_SNP.csv")%>%left_join(cq)
crt_na_rem <- crt %>% filter(!is.na(mean_cq))
crt_na_rem$Region[is.na(crt_na_rem$Region)] <- "West Africa"
crt_na_rem<- crt_na_rem%>% mutate(PR = PR / 100, Prev_DR = Prev_DR / 100)
write_csv(crt_na_rem, file = "pfcrt_database.csv")

###Model
Mod2_pfcrt <- lmer(Prev_DR ~ PR +mean_cq+PR *mean_cq +drug_1yr_ago+drug_2yr_ago+drug_3yr_ago+(1|country), weights=log(tested),data=crt_na_rem)
summary(Mod2_pfcrt)
plot(Mod2_pfcrt)
qqnorm(residuals(Mod2_pfcrt))

#Interpretation:
#The effect of CQ use on the prevalence of pfcrt-K76T mutation over 3 years lag period, suggested that CQ use negatively influenced pfcrt-K76T mutation in 3years, but this was not found to be significant. On the other hand, Parasite rate, and interaction between Parasite rate and CQ use were significantly (p<0.001 and p<0.05) associated with pfcrt-K76T mutatation, but with positive and negative estimates, respectively. 

##Plots
crt_na_rem<- crt_na_rem%>% mutate(PR = PR * 100, Prev_DR = Prev_DR * 100)
#### 1. Plot of Pfcrt_76T Prevalence vs Mean_cq_All regions
ggplot(crt_na_rem, aes(cut(mean_cq,10),Prev_DR,col=Region))+geom_boxplot()
ggsave("pfcrt_76T_p1.png") 

### 2. Plot of Pfcrt_76T Prevalence vs PR_All regions
ggplot(crt_na_rem, aes(cut(PR,10),Prev_DR,col=Region))+geom_boxplot()
ggsave("pfcrt_76T_p2.png")

### 3. Plot of the relationship between Region and Year varied between pfcrt 76T Prevalence.
ggplot(data = crt_na_rem %>% group_by(Region,year) %>% mutate(mean=mean(Prev_DR)), aes (y = Region, x = year,colour = mean)) + geom_point(aes(size=log(tested))) + scale_color_viridis_c(name="Pfcrt 76T Prevalence") + theme_classic() + xlab ("Year") + ylab ("Region")
ggsave("pfcrt 76T_Timep1.png")

### 4. Plot of the relationship between Region and Year varied between Parasite rate.
ggplot(data = crt_na_rem %>% group_by(Region,year) %>% mutate(mean=mean(PR)), aes (y = Region, x = year,colour = mean)) + geom_point(aes(size=log(tested))) + scale_color_viridis_c(name="PR") + theme_classic() + xlab ("Year") + ylab ("Region")
ggsave("pfcrt 76T_Timep2.png")

### 5. Plot of the relationship between Region and Year varied between Mean_cq.
ggplot(data = crt_na_rem %>% group_by(Region,year) %>% mutate(mean=mean(mean_cq)), aes (y = Region, x = year,colour = mean)) + geom_point(aes(size=log(tested))) + scale_color_viridis_c(name="Mean_cq") + theme_classic() + xlab ("Year") + ylab ("Region")
ggsave("pfcrt 76T_Timep3.png")


###PFMDR1 N86Y
#Prepare Database 

# join original observation
mdr<-read_csv ("pfmdr1_SNP.csv")%>%left_join(cq)
mdr_na_rem <- mdr %>% filter(!is.na(mean_cq))
mdr_na_rem$Region[is.na(mdr_na_rem$Region)] <- "West Africa"
mdr_na_rem<- mdr_na_rem%>% mutate(PR = PR / 100, Prev_DR = Prev_DR / 100)
write_csv(mdr_na_rem, file = "pfmdr1_database.csv")

###Model
Mod2_pfmdr1 <- lmer(Prev_DR ~ PR +mean_cq+PR *mean_cq+drug_1yr_ago+drug_2yr_ago+drug_3yr_ago+(1|country), weights=log(tested),data=mdr_na_rem)
summary(Mod2_pfmdr1)
plot(Mod2_pfmdr1)
qqnorm(residuals(Mod2_pfmdr1))

#Interpretation:
#The effect of yearly lag in CQ use on the prevalence of pfmdr1-N86Y mutation over 3 years, suggests that the prevalence of pfmdr1-N86Y increased with an increase in CQ use, but this was not significant. Similar to pfcrt- K76T, Parasite rate also significantly (p<0.001) increased as pfmdr1-N86Y mutation increased. However, the interaction between Parasite rate and CQ use was negatively associated with pfmdr1-N86Y prevalence and this was not significant.

##Plots
mdr_na_rem<- mdr_na_rem%>% mutate(PR = PR * 100, Prev_DR = Prev_DR * 100)
#### 1. Plot of Pfmdr1_86Y Prevalence vs Mean_cq_All regions
ggplot(mdr_na_rem, aes(cut(mean_cq,10),Prev_DR,col=Region))+geom_boxplot()
ggsave("pfmdr1_86Y_p1.png") 

### 2. Plot of Pfmdr1_86Y Prevalence vs PR_All regions
ggplot(mdr_na_rem, aes(cut(PR,10),Prev_DR,col=Region))+geom_boxplot()
ggsave("pfmdr1_86Y_p2.png")

### 3. Plot of the relationship between Region and Year varied between pfmdr1 86Y Prevalence.
ggplot(data = mdr_na_rem %>% group_by(Region,year) %>% mutate(mean=mean(Prev_DR)), aes (y = Region, x = year,colour = mean)) + geom_point(aes(size=log(tested))) + scale_color_viridis_c(name="Pfmdr1 86Y Prevalence") + theme_classic() + xlab ("Year") + ylab ("Region")
ggsave("pfmdr1 86Y_Timep1.png")

### 4. Plot of the relationship between Region and Year varied between Parasite rate.
ggplot(data = mdr_na_rem %>% group_by(Region,year) %>% mutate(mean=mean(PR)), aes (y = Region, x = year,colour = mean)) + geom_point(aes(size=log(tested))) + scale_color_viridis_c(name="PR") + theme_classic() + xlab ("Year") + ylab ("Region")
ggsave("pfmdr1 86Y_Timep2.png")

### 5. Plot of the relationship between Region and Year varied between Mean_cq.
ggplot(data = mdr_na_rem %>% group_by(Region,year) %>% mutate(mean=mean(mean_cq)), aes (y = Region, x = year,colour = mean)) + geom_point(aes(size=log(tested))) + scale_color_viridis_c(name="Mean_cq") + theme_classic() + xlab ("Year") + ylab ("Region")
ggsave("pfmdr1 86Y_Timep3.png")


###PFDHPS A437G
#Prepare Database

sp<- read_csv("SP_imputed.csv")
#Set up lag years
sp<-sp%>%group_by(country)%>%mutate(drug_1yr_ago = lag(mean_sp,n=1),drug_2yr_ago = lag(mean_sp,n=2), drug_3yr_ago = lag(mean_sp,n=3),drug_4yr_ago = lag(mean_sp,n=4),drug_5yr_ago = lag(mean_sp,n=5))

# join original observation
dhps<-read_csv ("pfdhps_SNP.csv")%>%left_join(sp)
dhps_na_rem <- dhps %>% filter(!is.na(mean_sp))
dhps_na_rem<- dhps_na_rem%>% mutate(PR = PR / 100, Prev_DR = Prev_DR / 100)
write_csv(dhps_na_rem, file = "pfdhps_database.csv")

###Model
Mod1_pfdhps <- lmer(Prev_DR ~ PR +mean_sp +PR *mean_sp+drug_1yr_ago+drug_2yr_ago+drug_3yr_ago+ (1|country), weights=log(tested),data=dhps_na_rem)
summary(Mod1_pfdhps)
plot(Mod1_pfdhps)
qqnorm(residuals(Mod1_pfdhps))

#Interpretation:
#SP use was found to negatively influence pfdhps-A437G prevalence at lead times of up to 3 years, howbeit not significantly. Similarly, Parasite rate negatively correlated with pfdhps-A437G and this was not significant. The interaction between Parasite rate and SP use did not significantly influence pfdhps-A437G prevalence either.  

##Plots
dhps_na_rem<- dhps_na_rem%>% mutate(PR = PR * 100, Prev_DR = Prev_DR * 100)
#### 1. Plot of Pfdhps_437G Prevalence vs Mean_sp_All regions
ggplot(dhps_na_rem, aes(cut(mean_sp,10),Prev_DR,col=Region))+geom_boxplot()
ggsave("pfdhps_437G _p1.png")

### 2. Plot of Pfdhps_437G  Prevalence vs PR_All regions
ggplot(dhps_na_rem, aes(cut(PR,10),Prev_DR,col=Region))+geom_boxplot()
ggsave("pfdhps_437G _p2.png")

### 3. Plot of the relationship between Region and Year varied between pfdhps 437G Prevalence.
ggplot(data = dhps_na_rem %>% group_by(Region,year) %>% mutate(mean=mean(Prev_DR)), aes (y = Region, x = year,colour = mean)) + geom_point(aes(size=log(tested))) + scale_color_viridis_c(name="Pfdhps 437G Prevalence") + theme_classic() + xlab ("Year") + ylab ("Region")
ggsave("pfdhps 437G_Timep1.png")

### 4. Plot of the relationship between Region and Year varied between Parasite rate.
ggplot(data = dhps_na_rem %>% group_by(Region,year) %>% mutate(mean=mean(PR)), aes (y = Region, x = year,colour = mean)) + geom_point(aes(size=log(tested))) + scale_color_viridis_c(name="PR") + theme_classic() + xlab ("Year") + ylab ("Region")
ggsave("pfdhps 437G_Timep2.png")

### 5. Plot of the relationship between Region and Year varied between Mean_sp.
ggplot(data = dhps_na_rem %>% group_by(Region,year) %>% mutate(mean=mean(mean_sp)), aes (y = Region, x = year,colour = mean)) + geom_point(aes(size=log(tested))) + scale_color_viridis_c(name="Mean_sp") + theme_classic() + xlab ("Year") + ylab ("Region")
ggsave("pfdhps 437G_Timep3.png")


###PFDHFR S108N
#Prepare Database

# join original observation
dhfr<-read_csv ("pfdhfr_SNP.csv")%>%left_join(sp)
dhfr_na_rem <- dhfr %>% filter(!is.na(mean_sp))
dhfr_na_rem <- dhfr_na_rem%>% mutate(PR = PR / 100, Prev_DR = Prev_DR / 100)
write_csv(dhfr_na_rem, file = "pfdhfr_database.csv")

###Model
Mod1_pfdhfr <- lmer(Prev_DR ~ PR +mean_sp+PR *mean_sp+drug_1yr_ago+drug_2yr_ago+drug_3yr_ago+ (1|country), weights=log(tested),data=dhfr_na_rem)
summary(Mod1_pfdhfr)
plot(Mod1_pfdhfr)
qqnorm(residuals(Mod1_pfdhfr))

#Interpretation:
#Like in the case of pfdhps-A437G prevalence, the effect of SP use on the prevalence of pfdhfr-S108N mutation over 3 years lag period suggested that SP use negatively influenced pfdhfr-S108N mutation in 3years, but this was not found to be significant. On the other hand, Parasite rate, and the interaction between Parasite rate and SP use were significantly (p<0.001 and p<0.05) associated with pfdhfr-S108N prevalence, but with negative and positive estimates, respectively. 


###Plots
dhfr_na_rem <- dhfr_na_rem%>% mutate(PR = PR * 100, Prev_DR = Prev_DR * 100)
#### 1. Plots of Pfdhfr_108N Prevalence vs Mean_sp_All regions
ggplot(dhfr_na_rem, aes(cut(mean_sp,10),Prev_DR,col=Region))+geom_boxplot()
ggsave("pfdhfr_108N_p1.png") 

### 2. Plots of Pfdhfr_108N Prevalence vs PR_All regions
ggplot(dhfr_na_rem, aes(cut(PR,10),Prev_DR,col=Region))+geom_boxplot()
ggsave("pfdhfr_108N_p2.png")

### 3. Plot of the relationship between Region and Year varied between pfdhps 108N Prevalence.
ggplot(data = dhfr_na_rem %>% group_by(Region,year) %>% mutate(mean=mean(Prev_DR)), aes (y = Region, x = year,colour = mean)) + geom_point(aes(size=log(tested))) + scale_color_viridis_c(name="Pfdhfr 108N Prevalence") + theme_classic() + xlab ("Year") + ylab ("Region")
ggsave("pfdhfr 108N_Timep1.png")

### 4. Plot of the relationship between Region and Year varied between Parasite rate.
ggplot(data = dhfr_na_rem %>% group_by(Region,year) %>% mutate(mean=mean(PR)), aes (y = Region, x = year,colour = mean)) + geom_point(aes(size=log(tested))) + scale_color_viridis_c(name="PR") + theme_classic() + xlab ("Year") + ylab ("Region")
ggsave("pfdhfr 108N_Timep2.png")

### 5. Plot of the relationship between Region and Year varied between Mean_sp.
ggplot(data = dhfr_na_rem %>% group_by(Region,year) %>% mutate(mean=mean(mean_sp)), aes (y = Region, x = year,colour = mean)) + geom_point(aes(size=log(tested))) + scale_color_viridis_c(name="Mean_sp") + theme_classic() + xlab ("Year") + ylab ("Region")
ggsave("pfdhfr 108N_Timep3.png")


#PFK13 C580Y
act<- read_csv("ACT_imputed.csv")

#Set up lag years
act<-act%>%group_by(country)%>%mutate(drug_1yr_ago = lag(mean_act,n=1),drug_2yr_ago = lag(mean_act,n=2), drug_3yr_ago = lag(mean_act,n=3),drug_4yr_ago = lag(mean_act,n=4),drug_5yr_ago = lag(mean_act,n=5))

# join original observation
k13<-read_csv ("pfk13_SNP.csv")%>%left_join(act)
k13_na_rem <- k13 %>% filter(!is.na(mean_act))
k13_na_rem<- k13_na_rem %>% mutate(PR = PR / 100, Prev_DR = Prev_DR / 100)
write_csv(k13_na_rem, file = "pfk13_database.csv")

###Model
k13_na_rem <-k13_na_rem %>% filter(Prev_DR<1,Prev_DR>0)
Mod1_pfk13 <- lmer(Prev_DR ~ PR +mean_act+PR *mean_act +drug_1yr_ago+drug_2yr_ago+drug_3yr_ago+ (1|country), weights=log(tested),data=k13_na_rem)
summary(Mod1_pfk13)
plot(Mod1_pfk13)
qqnorm(residuals(Mod1_pfk13))

#Interpretation:
#ACT use was found to positively influence pfk13-C580Y prevalence at lead times of up to 3 years, although not significantly. However, Parasite rate negatively correlated with pfk13-C580Y mutation and this was also not significant. The interaction between Parasite rate and ACT use did not significantly influence pfk13-C580Y prevalence either. 

###Plots
k13_na_rem<- k13_na_rem %>% mutate(PR = PR * 100, Prev_DR = Prev_DR * 100)
#### 1. Plots of Pfk13_580Y Prevalence vs Mean_act_All regions
ggplot(k13_na_rem, aes(cut(mean_act,10),Prev_DR,col=country))+geom_boxplot()
ggsave("pfK13_580Y_p1.png") 

### 2. Plots of Pfdhfr_108N Prevalence vs PR_All regions
ggplot(k13_na_rem, aes(cut(PR,10),Prev_DR,col=country))+geom_boxplot()
ggsave("pfk13_580Y_p2.png")

### 3. Plot of the relationship between Region and Year varied between pfk13 580Y Prevalence.
ggplot(data = k13_na_rem %>% group_by(country,year) %>% mutate(mean=mean(Prev_DR)), aes (y = country, x = year,colour = mean)) + geom_point(aes(size=log(tested))) + scale_color_viridis_c(name="Pfk13 580Y Prevalence") + theme_classic() + xlab ("Year") + ylab ("Country")
ggsave("pfk13 580Y_Timep1.png")

### 4. Plot of the relationship between Region and Year varied between Parasite rate.
ggplot(data = k13_na_rem %>% group_by(country,year) %>% mutate(mean=mean(PR)), aes (y = country, x = year,colour = mean)) + geom_point(aes(size=log(tested))) + scale_color_viridis_c(name="PR") + theme_classic() + xlab ("Year") + ylab ("Country")
ggsave("pfk13 580Y_Timep2.png")

### 5. Plot of the relationship between Region and Year varied between Mean_act.
ggplot(data = k13_na_rem %>% group_by(country,year) %>% mutate(mean=mean(mean_act)), aes (y = country, x = year,colour = mean)) + geom_point(aes(size=log(tested))) + scale_color_viridis_c(name="Mean_act") + theme_classic() + xlab ("Year") + ylab ("Country")
ggsave("pfk13 580Y_Timep3.png")
