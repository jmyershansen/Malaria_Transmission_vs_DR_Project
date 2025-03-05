
library(tidyverse)
library(dplyr)
library(lme4)
library(lmerTest)
library(ggplot2)
install.packages ("forcats")
library(forcats)
#Summary:
#Multilevel mixed regression models coupled with time lag framework was used to access the effect of key malaria indicators: Parasite rate, PR (transmission intensity) and Drug use (Proportion of fevered children treated with CQ, SP and ACTs) on the reported prevalence of key drug (CQ, SP and ACTs) resistant markers: pfcrt-K76T, pfmdr1-N86Y, pfdhps-A437G, pfdhfr-S108N and pfk13-C580Y. independently.


###PFCRT K76T
#Prepare Database
cq<- read_csv("CQ_imputed.csv")
#Set up lag years
cq<-cq%>%group_by(country)%>%mutate(drug_1yr_ago = lag(mean_cq,n=1),drug_2yr_ago = lag(mean_cq,n=2), drug_3yr_ago = lag(mean_cq,n=3),drug_4yr_ago = lag(mean_cq,n=4),drug_5yr_ago = lag(mean_cq,n=5))

# join original observation
CRT<-read_csv ("pfcrt_SNP.csv")
CRT <- CRT[!is.na(CRT$PR),]
#write_csv(CRT, file = "inspect.csv")
#write_csv(CRT, file = "inspect.1.csv")
crt <- left_join(CRT,cq, by = c("country","year"))
#crt<-read_csv ("pfcrt_SNP.csv")%>%left_join(cq)
crt_na_rem <- crt %>% filter(!is.na(mean_cq))
#crt_na_rem <- crt_na_rem[!is.na(crt_na_rem$PR),]
#write_csv(crt_na_rem, file = "inspect.csv")
CRT$Region[is.na(CRT$Region)] <- "West Africa"
crt_na_rem$Region[is.na(crt_na_rem$Region)] <- "West Africa"
crt_na_rem$Region[crt_na_rem$country == "Madagascar"] <-"MAD"
crt_na_rem<- crt_na_rem%>% mutate(PR = PR / 100, Prev_DR = Prev_DR / 100)

crt_na_rem <- crt_na_rem$country %>% unique()
crt_na_rem <- crt_na_rem[!duplicated(crt_na_rem),]
Nstudies0 <- length(unique(crt_na_rem$`study Id`))
Nstudies1 <- length(unique(crt_na_rem$`PR`))
Nstudies2 <- length(unique(crt_na_rem$`Prev_DR`))
Nstudies3 <- length(unique(crt_na_rem$`mean_cq`))

write_csv(crt_na_rem, file = "inspect.csv")
write_csv(crt_na_rem, file = "pfcrt_database.csv")

###Model
Mod2_pfcrt <- lmer(Prev_DR ~ PR +mean_cq+PR *mean_cq +drug_1yr_ago+drug_2yr_ago+drug_3yr_ago+(1|Region), 
                   weights=log(tested),data=crt_na_rem)
summary(Mod2_pfcrt)
plot(Mod2_pfcrt)
qqnorm(residuals(Mod2_pfcrt))



### Recode regions

crt_na_rem<-crt_na_rem%>%mutate(TransZone = case_when(
  Region %in% c("West Africa","Central Africa") ~"High",
  Region %in% c("West Africa","Central Africa") ~"High",
  Region %in% c("East Africa","Southern Africa") ~"Med",
  Region %in% c("East Africa","Southern Africa") ~"Med",
  Region == "MAD" ~"MAD",
  Region %in% c("Asia", "SEA","South America")~"Low"))


###Model
Mod2_pfcrt <- lmer(Prev_DR ~ PR +mean_cq+PR *mean_cq +drug_1yr_ago+drug_2yr_ago+drug_3yr_ago+(1|TransZone), 
                   weights=log(tested),data=crt_na_rem)
summary(Mod2_pfcrt)
plot(Mod2_pfcrt)
qqnorm(residuals(Mod2_pfcrt))
###
###

Mod2_pfcrt <- crt_na_rem%>%filter(!is.na(drug_3yr_ago))
Mod2_pfcrt <-Mod2_pfcrt%>%mutate(cumDrug = ifelse(is.na(drug_3yr_ago),drug_2yr_ago,drug_3yr_ago)+
                                 drug_2yr_ago+drug_1yr_ago)
Mod2_pfcrt <- lmer(Prev_DR ~ PR +cumDrug+PR*cumDrug+(1|Region), 
                   weights=log(tested),data=Mod2_pfcrt)
summary(Mod2_pfcrt)
#plot(Mod2_pfcrt)
#qqnorm(residuals(Mod2_pfcrt))
###
###

#Modeling at High Transzone
highTrans <- crt_na_rem%>%filter(TransZone=="High")%>%filter(!is.na(drug_3yr_ago))

highTrans <-highTrans%>%mutate(cumDrug = ifelse(is.na(drug_3yr_ago),drug_2yr_ago,drug_3yr_ago)+
                                 drug_2yr_ago+drug_1yr_ago)

high_pfcrt <- lmer(Prev_DR ~ PR +cumDrug+PR*cumDrug+(1|Region), 
                   weights=log(tested),data=highTrans)
summary(high_pfcrt)
plot(high_pfcrt)
qqnorm(residuals(high_pfcrt))

high_pfcrt <- lm(Prev_DR ~ PR +cumDrug+PR*cumDrug, 
                   weights=log(tested),data=highTrans)
summary(high_pfcrt)

grid <- highTrans %>%  
  add_predictions(high_pfcrt)
ggplot(grid, aes(cumDrug,Prev_DR,col=Region))+
  geom_point(aes(color=country))+
  geom_line(aes(y=pred))+
  facet_grid(~Region)

ggplot(grid, aes(cumDrug,Prev_DR,col=Region))+
  geom_point(aes(color=country))+
  geom_line(aes(y=pred))+
  facet_grid(cut(PR,4)~Region)
ggplot(grid, aes(cumDrug,Prev_DR,col=Region))+
  geom_point(aes(color=country,size=log(tested)))+
  geom_line(aes(y=pred))

crt_na_rem <- crt_na_rem%>%add_predictions(Mod2_pfcrt)
ggplot(crt_na_rem, aes(drug_1yr_ago,Prev_DR,col=Region))+
  geom_point(aes(color=country))+
  geom_line(aes(y=pred))+
  facet_grid(cut(PR,4)~TransZone,scales="free_x")

ggplot(crt_na_rem, aes(drug_1yr_ago,PR,col=Region))+
  geom_point(aes(color=country))+
  geom_line(aes(y=pred))+
  facet_grid(cut(Prev_DR,4)~TransZone,scales="free_x")

library(modelr)

plot(high_pfcrt)

summary(high_pfcrt)
ggsave("pfcrt_76T_p1.png") 
###
###

#Modeling at Moderate Transzone

medTrans <- crt_na_rem%>%filter(TransZone=="Med")%>%filter(!is.na(drug_3yr_ago))

medTrans <-medTrans%>%mutate(cumDrug = ifelse(is.na(drug_3yr_ago),drug_2yr_ago,drug_3yr_ago)+
                                 drug_2yr_ago+drug_1yr_ago)

med_pfcrt <- lmer(Prev_DR ~ PR +cumDrug+PR*cumDrug+(1|Region), 
                   weights=log(tested),data=medTrans)
summary(med_pfcrt)
plot(med_pfcrt)
qqnorm(residuals(med_pfcrt))

med_pfcrt <- lm(Prev_DR ~ PR +cumDrug+PR*cumDrug, 
                 weights=log(tested),data=medTrans)
summary(med_pfcrt)


#Modeling at Low Transzone

lowTrans <- crt_na_rem%>%filter(TransZone=="Low")%>%filter(!is.na(drug_3yr_ago))

lowTrans <-lowTrans%>%mutate(cumDrug = ifelse(is.na(drug_3yr_ago),drug_2yr_ago,drug_3yr_ago)+
                                 drug_2yr_ago+drug_1yr_ago)

low_pfcrt <- lmer(Prev_DR ~ PR +cumDrug+PR*cumDrug+(1|Region), 
                   weights=log(tested),data=lowTrans)
summary(low_pfcrt)
plot(low_pfcrt)
qqnorm(residuals(low_pfcrt))

low_pfcrt <- lm(Prev_DR ~ PR +cumDrug+PR*cumDrug, 
                 weights=log(tested),data=lowTrans)
summary(low_pfcrt)

#Interpretation:
#The effect of CQ use on the prevalence of pfcrt-K76T mutation over 3 years lag period, suggested that CQ use negatively influenced pfcrt-K76T mutation in 3years, but this was not found to be significant. On the other hand, Parasite rate, and interaction between Parasite rate and CQ use were significantly (p<0.001 and p<0.05) associated with pfcrt-K76T mutatation, but with positive and negative estimates, respectively. 

#crt_na_rem<-read_csv ("pfcrt_database.csv")

##Plots
crt_na_rem<- crt_na_rem%>% mutate(PR = PR * 100, Prev_DR = Prev_DR * 100)
#### 1. Plot of Pfcrt_76T Prevalence vs Mean_cq_All regions
ggplot(crt_na_rem, aes(cut(mean_cq,10),Prev_DR,col=Region))+ geom_boxplot()+ xlab("CQ Usage") + ylab("Pfcrt 76T Prevalence")+ facet_grid(~TransZone)
#ggplot(crt_na_rem, aes(cut(drug_1yr_ago+drug_2yr_ago+drug_3yr_ago,10),Prev_DR,col=Region))+ geom_boxplot()+ xlab("Drug Usage") + ylab("Pfcrt 76T Prevalence")
#ggplot(crt_na_rem, aes(cut(drug_1yr_ago+drug_2yr_ago+drug_3yr_ago,10),Prev_DR,col=Region))+geom_boxplot()+ xlab("Drug Usage") + ylab("Pfcrt 76T Prevalence")+ 
#facet_grid(~TransZone)

crt_na_rem <- crt_na_rem%>%add_predictions(high_pfcrt)
ggplot(crt_na_rem, aes(cut(drug_1yr_ago+drug_2yr_ago+drug_3yr_ago,10),Prev_DR,col=Region))+
  geom_point(aes(color=country))+ geom_smooth(method = "lm", aes(group = country), se = F) +
  xlab("Drug Usage") + ylab("Pfcrt 76T Prevalence")+ 
  facet_grid(~TransZone)

ggplot(crt_na_rem, aes(cut(mean_cq,10),Prev_DR,col=Region))+
  geom_point(aes(color=country))+ geom_smooth(method = "lm", aes(group = country), se = F) + 
  xlab("CQ Usage") + ylab("Pfcrt 76T Prevalence")+ 
  facet_grid(~TransZone)


ggplot(crt_na_rem%>%filter(TransZone=="High"), aes(cut(drug_3yr_ago,10),Prev_DR,col=Region))+
  geom_boxplot(aes(color=country))+ xlab("Drug Usage") + ylab("Pfcrt 76T Prevalence")

ggplot(crt_na_rem%>%filter(TransZone=="High"), aes(cut(mean_cq,10),Prev_DR,col=Region))+
  geom_boxplot(aes(color=country))+ xlab("Drug Usage") + ylab("Pfcrt 76T Prevalence")

ggplot(crt_na_rem%>%filter(TransZone=="High"), aes(cut(mean_cq,10),Prev_DR,col=country))+
  geom_point(aes(color=country))+ geom_smooth(method = "lm", aes(group = country), se = F) + 
  xlab("Drug Usage") + ylab("Pfcrt 76T Prevalence")

###
###
#ggplot(crt_na_rem, aes(x = mean_cq, y = Prev_DR)) + 
          geom_point(aes(color=as.factor(Region), size=tested))+
        geom_smooth(method = "lm") + 
          xlab("CQ use") + ylab("Pfcrt 76T Prevalence")

ggplot(crt_na_rem%>%filter(TransZone=="High"), aes(x = mean_cq, y = Prev_DR)) + 
  geom_point(aes(color=as.factor(Region), size=tested))+
  geom_smooth(method = "lm") + 
  xlab("CQ use") + ylab("Pfcrt 76T Prevalence")
ggsave("pfcrt_76T_p1A.png")

ggplot(crt_na_rem%>%filter(TransZone=="Low"), aes(x = mean_cq, y = Prev_DR)) + 
  geom_point(aes(color=as.factor(Region), size=tested))+
  geom_smooth(method = "lm") + 
  xlab("CQ use") + ylab("Pfcrt 76T Prevalence")
ggsave("pfcrt_76T_p1B.png")


#ggplot(crt_na_rem%>%filter(country=="Burkina Faso"), aes(drug_3yr_ago,Prev_DR,col=Region))+
#geom_point(aes(color=PR))
#ggsave("pfcrt_76T_p1.png")

### 2. Plot of Pfcrt_76T Prevalence vs PR_All regions
ggplot(crt_na_rem, aes(cut(PR,10),Prev_DR,col=Region))+geom_boxplot()+ xlab("Parasite Rate") + ylab("Pfcrt 76T Prevalence")
ggplot(crt_na_rem, aes(cut(PR,10),Prev_DR,col=Region))+geom_boxplot()+ xlab("Parasite Rate") + ylab("Pfcrt 76T Prevalence")+ facet_grid(~TransZone)
ggsave("pfcrt_76T_p2.png")

###
###
#ggplot(crt_na_rem, aes(x = PR, y = Prev_DR)) + 
geom_point(aes(color=as.factor(Region), size=tested))+
  geom_smooth(method = "lm") + 
  xlab("Parasite rate") + ylab("Pfcrt 76T Prevalence")

ggplot(crt_na_rem%>%filter(TransZone=="High"), aes(x = PR, y = Prev_DR)) + 
  geom_point(aes(color=as.factor(Region), size=tested))+
  geom_smooth(method = "lm") + 
  xlab("Parasite rate") + ylab("Pfcrt 76T Prevalence")
ggsave("pfcrt_76T_p2A.png")

ggplot(crt_na_rem%>%filter(TransZone=="Low"), aes(x = PR, y = Prev_DR)) + 
  geom_point(aes(color=as.factor(Region), size=tested))+
  geom_smooth(method = "lm") + 
  xlab("Parasite rate") + ylab("Pfcrt 76T Prevalence")
ggsave("pfcrt_76T_p2B.png")


### 3. Plot of the relationship between Region and Year varied between pfcrt 76T Prevalence.
ggplot(data = crt_na_rem %>% group_by(Region,year) %>% mutate(mean=mean(Prev_DR)), aes (y = Region, x = year,colour = mean)) + geom_point(aes(size=log(tested))) + scale_color_viridis_c(name="Pfcrt 76T Prevalence") + theme_classic() + xlab ("Year") + ylab ("Region")
ggplot(data = CRT %>% group_by(Region,year) %>% mutate(mean=mean(Prev_DR)), aes (y = Region, x = year,colour = mean)) + geom_point(aes(size=log(tested))) + scale_color_viridis_c(name="Pfcrt 76T Prevalence") + theme_classic() + xlab ("Year") + ylab ("Region")+ ggtitle ("Time plot of Pfcrt_76T over time")

ggplot(data = CRT %>% group_by(Region,year) %>% mutate(Region = factor(Region, levels = c("West Africa", "Central Africa" , "East Africa", "Southern Africa", "North Africa", "Caribbean", "Central America", "South America", "SEA", "Asia", "Oceania")), mean=mean(Prev_DR))) + aes (y = Region, x = year,colour = mean) + geom_point(aes(size=log(tested))) + scale_color_viridis_c(name="Pfcrt 76T Prevalence") + theme_classic() + xlab ("Year") + ylab ("Region")+ ggtitle ("Time plot of Pfcrt_76T over time")
ggsave("pfcrt 76T_Timep1.png")
ggsave("pfcrt 76T_Timep1A.png")

### 4. Plot of the relationship between Region and Year varied between Parasite rate.
ggplot(data = crt_na_rem %>% group_by(Region,year) %>% mutate(mean=mean(PR)), aes (y = Region, x = year,colour = mean)) + geom_point(aes(size=log(tested))) + scale_color_viridis_c(name="PR") + theme_classic() + xlab ("Year") + ylab ("Region")
ggplot(data = CRT %>% group_by(Region,year) %>% mutate(mean=mean(PR)), aes (y = Region, x = year,colour = mean)) + geom_point(aes(size=log(tested))) + scale_color_viridis_c(name="Parasite Rate") + theme_classic() + xlab ("Year") + ylab ("Region") + ggtitle ("Time plot of PR over time")

ggplot(data = CRT %>% group_by(Region,year) %>% mutate(Region = factor(Region, levels = c("West Africa", "Central Africa" , "East Africa", "Southern Africa", "North Africa", "Caribbean", "Central America", "South America", "SEA", "Asia", "Oceania")), mean=mean(PR))) + aes (y = Region, x = year,colour = mean) + geom_point(aes(size=log(tested))) + scale_color_viridis_c(name="Parasite Rate") + theme_classic() + xlab ("Year") + ylab ("Region")+ ggtitle ("Time plot of PR over time")
ggsave("pfcrt 76T_Timep2.png")
ggsave("pfcrt 76T_Timep2A.png")


### 5. Plot of the relationship between Region and Year varied between Mean_cq.
ggplot(data = crt_na_rem %>% group_by(Region,year) %>% mutate(mean=mean(mean_cq)), aes (y = Region, x = year,colour = mean)) + geom_point(aes(size=log(tested))) + scale_color_viridis_c(name="CQ usage") + theme_classic() + xlab ("Year") + ylab ("Region")+ ggtitle ("Time plot of CQ usage over time")

ggplot(data = crt_na_rem %>% group_by(Region,year) %>% mutate(Region = factor(Region, levels = c("West Africa", "Central Africa" , "East Africa", "Southern Africa", "North Africa", "Caribbean", "Central America", "South America", "SEA", "Asia", "Oceania")), mean=mean(mean_cq))) + aes (y = Region, x = year,colour = mean) + geom_point(aes(size=log(tested))) + scale_color_viridis_c(name="CQ usage") + theme_classic() + xlab ("Year") + ylab ("Region")+ ggtitle ("Time plot of CQ usage over time")
ggsave("pfcrt 76T_Timep3.png")
ggsave("pfcrt 76T_Timep3A.png")


###PFMDR1 N86Y
#Prepare Database 
# join original observation
MDR<-read_csv ("pfmdr1_SNP.csv")
MDR <- MDR[!is.na(MDR$PR),]
mdr <- left_join(MDR,cq, by = c("country","year"))
#mdr<-read_csv ("pfmdr1_SNP.csv")%>%left_join(cq)
mdr_na_rem <- mdr %>% filter(!is.na(mean_cq))
#mdr_na_rem <- mdr_na_rem[!is.na(mdr_na_rem$PR),]
#write_csv(mdr_na_rem, file = "inspect.1.csv")
MDR$Region[is.na(MDR$Region)] <- "West Africa"
mdr_na_rem$Region[is.na(mdr_na_rem$Region)] <- "West Africa"
mdr_na_rem$Region[mdr_na_rem$country == "Madagascar"] <-"MAD"
mdr_na_rem<- mdr_na_rem%>% mutate(PR = PR / 100, Prev_DR = Prev_DR / 100)
write_csv(mdr_na_rem, file = "pfmdr1_database.csv")

###Model
Mod2_pfmdr1 <- lmer(Prev_DR ~ PR +mean_cq+PR *mean_cq+drug_1yr_ago+drug_2yr_ago+drug_3yr_ago+
                      (1+drug_1yr_ago|Region), weights=log(tested),data=mdr_na_rem)
Mod2_pfmdr1 <- lmer(Prev_DR ~ PR +drug_1yr_ago+
                      (1+drug_1yr_ago|Region), weights=log(tested),data=mdr_na_rem)
summary(Mod2_pfmdr1)
plot(Mod2_pfmdr1)
qqnorm(residuals(Mod2_pfmdr1))

mdr_na_rem<- mdr_na_rem%>%add_predictions(Mod2_pfmdr1)
ggplot(mdr_na_rem, aes(drug_2yr_ago+drug_1yr_ago,Prev_DR,col=Region))+
  geom_point(aes(color=country))+
  geom_line(aes(y=pred))+
  facet_grid(~Region)


### Recode regions

mdr_na_rem<-mdr_na_rem%>%mutate(TransZone = case_when(
  Region %in% c("West Africa","Central Africa") ~"High",
  Region %in% c("West Africa","Central Africa") ~"High",
  Region %in% c("East Africa","Southern Africa") ~"Med",
  Region %in% c("East Africa","Southern Africa") ~"Med",
  Region == "MAD" ~"MAD",
  Region %in% c("Asia", "SEA","South America")~"Low"))


###Model
Mod2_pfmdr1 <- lmer(Prev_DR ~ PR +mean_cq+PR *mean_cq +drug_1yr_ago+drug_2yr_ago+drug_3yr_ago+(1|TransZone), 
                   weights=log(tested),data=mdr_na_rem)
summary(Mod2_pfmdr1)
plot(Mod2_pfmdr1)
qqnorm(residuals(Mod2_pfmdr1))
###
###

Mod2_pfmdr1 <- mdr_na_rem%>%filter(!is.na(drug_3yr_ago))
Mod2_pfmdr1 <-Mod2_pfmdr1%>%mutate(cumDrug = ifelse(is.na(drug_3yr_ago),drug_2yr_ago,drug_3yr_ago)+
                                 drug_2yr_ago+drug_1yr_ago)
Mod2_pfmdr1 <- lmer(Prev_DR ~ PR +cumDrug+PR*cumDrug+(1|Region), 
                    weights=log(tested),data=Mod2_pfmdr1)
summary(Mod2_pfmdr1)
plot(Mod2_pfmdr1)
qqnorm(residuals(Mod2_pfmdr1))

###
###

#Modeling at High Transzone
highTrans <- mdr_na_rem%>%filter(TransZone=="High")%>%filter(!is.na(drug_3yr_ago))

highTrans <-highTrans%>%mutate(cumDrug = ifelse(is.na(drug_3yr_ago),drug_2yr_ago,drug_3yr_ago)+
                                 drug_2yr_ago+drug_1yr_ago)

high_pfmdr1 <- lmer(Prev_DR ~ PR +cumDrug+PR*cumDrug+(1|Region), 
                   weights=log(tested),data=highTrans)
summary(high_pfmdr1)
plot(high_pfmdr1)
qqnorm(residuals(high_pfmdr1))

high_pfmdr1 <- lm(Prev_DR ~ PR +cumDrug+PR*cumDrug, 
                 weights=log(tested),data=highTrans)
summary(high_pfmdr1)

#Modeling at Moderate Transzone

medTrans <- mdr_na_rem%>%filter(TransZone=="Med")%>%filter(!is.na(drug_3yr_ago))

medTrans <-medTrans%>%mutate(cumDrug = ifelse(is.na(drug_3yr_ago),drug_2yr_ago,drug_3yr_ago)+
                               drug_2yr_ago+drug_1yr_ago)

med_pfmdr1 <- lmer(Prev_DR ~ PR +cumDrug+PR*cumDrug+(1|Region), 
                  weights=log(tested),data=medTrans)
summary(med_pfmdr1)
plot(med_pfmdr1)
qqnorm(residuals(med_pfmdr1))

med_pfmdr1 <- lm(Prev_DR ~ PR +cumDrug+PR*cumDrug, 
                weights=log(tested),data=medTrans)
summary(med_pfmdr1)


#Modeling at Low Transzone

lowTrans <- mdr_na_rem%>%filter(TransZone=="Low")%>%filter(!is.na(drug_3yr_ago))

lowTrans <-lowTrans%>%mutate(cumDrug = ifelse(is.na(drug_3yr_ago),drug_2yr_ago,drug_3yr_ago)+
                               drug_2yr_ago+drug_1yr_ago)

low_pfmdr1 <- lmer(Prev_DR ~ PR +cumDrug+PR*cumDrug+(1|Region), 
                  weights=log(tested),data=lowTrans)
summary(low_pfmdr1)
plot(low_pfmdr1)
qqnorm(residuals(low_pfmdr1))

low_pfmdr1 <- lm(Prev_DR ~ PR +cumDrug+PR*cumDrug, 
                weights=log(tested),data=lowTrans)
summary(low_pfmdr1)

#Interpretation:
#The effect of yearly lag in CQ use on the prevalence of pfmdr1-N86Y mutation over 3 years, suggests that the prevalence of pfmdr1-N86Y increased with an increase in CQ use, but this was not significant. Similar to pfcrt- K76T, Parasite rate also significantly (p<0.001) increased as pfmdr1-N86Y mutation increased. However, the interaction between Parasite rate and CQ use was negatively associated with pfmdr1-N86Y prevalence and this was not significant.

##Plots
mdr_na_rem<- mdr_na_rem%>% mutate(PR = PR * 100, Prev_DR = Prev_DR * 100)
#### 1. Plot of Pfmdr1_86Y Prevalence vs Mean_cq_All regions
ggplot(mdr_na_rem, aes(cut(mean_cq,10),Prev_DR,col=Region))+geom_boxplot()+ xlab("Drug Usage") + ylab("Pfmdr1 86Y Prevalence")
ggplot(mdr_na_rem, aes(cut(mean_cq,10),Prev_DR,col=Region))+geom_boxplot()+ xlab("Drug Usage") + ylab("Pfmdr1 86Y Prevalence")+ facet_grid(~TransZone)
ggsave("pfmdr1_86Y_p1.png") 

### 2. Plot of Pfmdr1_86Y Prevalence vs PR_All regions
ggplot(mdr_na_rem, aes(cut(PR,10),Prev_DR,col=Region))+geom_boxplot()+ xlab("Parasite Rate") + ylab("Pfmdr1 86Y Prevalence")
ggplot(mdr_na_rem, aes(cut(PR,10),Prev_DR,col=Region))+geom_boxplot()+ xlab("Parasite Rate") + ylab("Pfmdr1 86Y Prevalence")+ facet_grid(~TransZone)
ggsave("pfmdr1_86Y_p2.png")

### 3. Plot of the relationship between Region and Year varied between pfmdr1 86Y Prevalence.
ggplot(data = mdr_na_rem %>% group_by(Region,year) %>% mutate(mean=mean(Prev_DR)), aes (y = Region, x = year,colour = mean)) + geom_point(aes(size=log(tested))) + scale_color_viridis_c(name="Pfmdr1 86Y Prevalence") + theme_classic() + xlab ("Year") + ylab ("Region")+ ggtitle ("Time plot of Pfmdr1 86Y over time")
ggplot(data = MDR %>% group_by(Region,year) %>% mutate(Region = factor(Region, levels = c("West Africa", "Central Africa" , "East Africa", "Southern Africa", "North Africa", "Caribbean", "Central America", "South America", "SEA", "Asia", "Oceania")), mean=mean(Prev_DR))) + aes (y = Region, x = year,colour = mean) + geom_point(aes(size=log(tested))) + scale_color_viridis_c(name="Pfmdr1 86Y Prevalence") + theme_classic() + xlab ("Year") + ylab ("Region")+ ggtitle ("Time plot of Pfmdr1 86Y over time")
ggsave("pfmdr1 86Y_Timep1.png")
ggsave("pfmdr1 86Y_Timep1A.png")

### 4. Plot of the relationship between Region and Year varied between Parasite rate.
ggplot(data = mdr_na_rem %>% group_by(Region,year) %>% mutate(mean=mean(PR)), aes (y = Region, x = year,colour = mean)) + geom_point(aes(size=log(tested))) + scale_color_viridis_c(name="PR") + theme_classic() + xlab ("Year") + ylab ("Region")
ggplot(data = MDR %>% group_by(Region,year) %>% mutate(Region = factor(Region, levels = c("West Africa", "Central Africa" , "East Africa", "Southern Africa", "North Africa", "Caribbean", "Central America", "South America", "SEA", "Asia", "Oceania")), mean=mean(PR))) + aes (y = Region, x = year,colour = mean) + geom_point(aes(size=log(tested))) + scale_color_viridis_c(name="Parasite Rate") + theme_classic() + xlab ("Year") + ylab ("Region")+ ggtitle ("Time plot of PR over time")
ggsave("pfmdr1 86Y_Timep2.png")
ggsave("pfmdr1 86Y_Timep2A.png")

### 5. Plot of the relationship between Region and Year varied between Mean_cq.
ggplot(data = mdr_na_rem %>% group_by(Region,year) %>% mutate(mean=mean(mean_cq)), aes (y = Region, x = year,colour = mean)) + geom_point(aes(size=log(tested))) + scale_color_viridis_c(name="Mean_cq") + theme_classic() + xlab ("Year") + ylab ("Region")
ggplot(data = mdr_na_rem %>% group_by(Region,year) %>% mutate(Region = factor(Region, levels = c("West Africa", "Central Africa" , "East Africa", "Southern Africa", "North Africa", "Caribbean", "Central America", "South America", "SEA", "Asia", "Oceania")), mean=mean(mean_cq))) + aes (y = Region, x = year,colour = mean) + geom_point(aes(size=log(tested))) + scale_color_viridis_c(name="CQ usage") + theme_classic() + xlab ("Year") + ylab ("Region")+ ggtitle ("Time plot of CQ usage over time")
ggsave("pfmdr1 86Y_Timep3.png")
ggsave("pfmdr1 86Y_Timep3A.png")


cq<- read_csv("CQ_imputed.csv")
#Set up lag years
cq<-cq%>%group_by(country)%>%mutate(drug_1yr_ago = lag(mean_cq,n=1),drug_2yr_ago = lag(mean_cq,n=2), drug_3yr_ago = lag(mean_cq,n=3),drug_4yr_ago = lag(mean_cq,n=4),drug_5yr_ago = lag(mean_cq,n=5))
# join original observation
CRT<-read_csv ("pfcrt_SNP.csv")
CRT <- CRT[!is.na(CRT$PR),]
#write_csv(CRT, file = "inspect.csv")
#write_csv(CRT, file = "inspect.1.csv")
crt <- left_join(CRT,cq, by = c("country","year"))
#crt<-read_csv ("pfcrt_SNP.csv")%>%left_join(cq)
crt_na_rem <- crt %>% filter(!is.na(mean_cq))
#crt_na_rem <- crt_na_rem[!is.na(crt_na_rem$PR),]
#write_csv(crt_na_rem, file = "inspect.csv")
CRT$Region[is.na(CRT$Region)] <- "West Africa"
crt_na_rem$Region[is.na(crt_na_rem$Region)] <- "West Africa"
crt_na_rem$Region[crt_na_rem$country == "Madagascar"] <-"MAD"


###PFDHPS A437G
#Prepare Database
sp<- read_csv("SP_imputed.csv")
#Set up lag years
sp<-sp%>%group_by(country)%>%mutate(drug_1yr_ago = lag(mean_sp,n=1),drug_2yr_ago = lag(mean_sp,n=2), drug_3yr_ago = lag(mean_sp,n=3),drug_4yr_ago = lag(mean_sp,n=4),drug_5yr_ago = lag(mean_sp,n=5))

# join original observation
DHPS<-read_csv ("pfdhps_SNP.csv")
DHPS <- DHPS[!is.na(DHPS$PR),]
dhps <- left_join(DHPS,sp, by = c("country","year"))
#dhps<-read_csv ("pfdhps_SNP.csv")%>%left_join(sp)
dhps_na_rem <- dhps %>% filter(!is.na(mean_sp))
#write_csv(dhps_na_rem, file = "inspect.1.csv")
DHPS$Region[is.na(DHPS$Region)] <- "West Africa"
dhps_na_rem$Region[is.na(dhps_na_rem$Region)] <- "West Africa"
dhps_na_rem$Region[dhps_na_rem$country == "Madagascar"] <-"MAD"
dhps_na_rem<- dhps_na_rem%>% mutate(PR = PR / 100, Prev_DR = Prev_DR / 100)
write_csv(dhps_na_rem, file = "pfdhps_database.csv")

###Model
Mod1_pfdhps <- lmer(Prev_DR ~ PR +mean_sp +PR *mean_sp+drug_1yr_ago+drug_2yr_ago+
                      drug_3yr_ago+ (1|Region), weights=log(tested),data=dhps_na_rem)
summary(Mod1_pfdhps)
plot(Mod1_pfdhps)
qqnorm(residuals(Mod1_pfdhps))
####
####

### Recode regions

dhps_na_rem<-dhps_na_rem%>%mutate(TransZone = case_when(
  Region %in% c("West Africa","Central Africa") ~"High",
  Region %in% c("West Africa","Central Africa") ~"High",
  Region %in% c("East Africa","Southern Africa") ~"Med",
  Region %in% c("East Africa","Southern Africa") ~"Med",
  Region == "MAD" ~"MAD",
  Region %in% c("Asia", "SEA","South America")~"Low"))

###Model
Mod2_pfdhps <- lmer(Prev_DR ~ PR +mean_sp+PR *mean_sp +drug_1yr_ago+drug_2yr_ago+drug_3yr_ago+(1|TransZone), 
                   weights=log(tested),data=dhps_na_rem)
summary(Mod2_pfdhps)
plot(Mod2_pfdhps)
qqnorm(residuals(Mod2_pfdhps))
###
###
Mod2_pfdhps <- dhps_na_rem%>%filter(!is.na(drug_3yr_ago))

Mod2_pfdhps <-Mod2_pfdhps%>%mutate(cumDrug = ifelse(is.na(drug_3yr_ago),drug_2yr_ago,drug_3yr_ago)+
                                 drug_2yr_ago+drug_1yr_ago)

Mod2_pfdhps <- lmer(Prev_DR ~ PR +cumDrug+PR*cumDrug+(1|'study id)', 
                    weights=log(tested),data=Mod2_pfdhps)
summary(Mod2_pfdhps)
plot(Mod2_pfdhps)
qqnorm(residuals(Mod2_pfdhps))
###

#Modeling at High Transzone
highTrans <- dhps_na_rem%>%filter(TransZone=="High")%>%filter(!is.na(drug_3yr_ago))

highTrans <-highTrans%>%mutate(cumDrug = ifelse(is.na(drug_3yr_ago),drug_2yr_ago,drug_3yr_ago)+
                                 drug_2yr_ago+drug_1yr_ago)

high_pfdhps <- lmer(Prev_DR ~ PR +cumDrug+PR*cumDrug+(1|Region), 
                   weights=log(tested),data=highTrans)
summary(high_pfdhps)
plot(high_pfdhps)
qqnorm(residuals(high_pfdhps))

high_pfdhps <- lm(Prev_DR ~ PR +cumDrug+PR*cumDrug, 
                 weights=log(tested),data=highTrans)
summary(high_pfdhps)
plot(high_pfdhps)
qqnorm(residuals(high_pfdhps))

#Modeling at Moderate Transzone
medTrans <- dhps_na_rem%>%filter(TransZone=="Med")%>%filter(!is.na(drug_3yr_ago))

medTrans <-medTrans%>%mutate(cumDrug = ifelse(is.na(drug_3yr_ago),drug_2yr_ago,drug_3yr_ago)+
                                 drug_2yr_ago+drug_1yr_ago)

med_pfdhps <- lmer(Prev_DR ~ PR +cumDrug+PR*cumDrug+(1|Region), 
                    weights=log(tested),data=medTrans)
summary(med_pfdhps)

med_pfdhps <- lm(Prev_DR ~ PR +cumDrug+PR*cumDrug, 
                  weights=log(tested),data=medTrans)
summary(med_pfdhps)


#Modeling at Low Transzone
lowTrans <- dhps_na_rem%>%filter(TransZone=="Low")%>%filter(!is.na(drug_3yr_ago))

lowTrans <-lowTrans%>%mutate(cumDrug = ifelse(is.na(drug_3yr_ago),drug_2yr_ago,drug_3yr_ago)+
                                 drug_2yr_ago+drug_1yr_ago)

low_pfdhps <- lmer(Prev_DR ~ PR +cumDrug+PR*cumDrug+(1|Region), 
                    weights=log(tested),data=lowTrans)
summary(low_pfdhps)

low_pfdhps <- lm(Prev_DR ~ PR +cumDrug+PR*cumDrug, 
                  weights=log(tested),data=lowTrans)
summary(low_pfdhps)


#Interpretation:
#SP use was found to negatively influence pfdhps-A437G prevalence at lead times of up to 3 years, howbeit not significantly. Similarly, Parasite rate negatively correlated with pfdhps-A437G and this was not significant. The interaction between Parasite rate and SP use did not significantly influence pfdhps-A437G prevalence either.  

##Plots
dhps_na_rem<- dhps_na_rem%>% mutate(PR = PR * 100, Prev_DR = Prev_DR * 100)
#### 1. Plot of Pfdhps_437G Prevalence vs Mean_sp_All regions
ggplot(dhps_na_rem, aes(cut(mean_sp,10),Prev_DR,col=Region))+geom_boxplot()+ xlab("SP Usage") + ylab("Pfdhps 437G Prevalence")
ggplot(dhps_na_rem, aes(cut(mean_sp,10),Prev_DR,col=Region))+geom_boxplot()+ xlab("SP Usage") + ylab("Pfdhps 437G Prevalence")+ facet_grid(~TransZone)
ggsave("pfdhps_437G _p1.png")

#ggplot(Mod2_pfdhps, aes(cut(cumDrug,10),Prev_DR,col=Region))+
geom_point(aes(color=country))+   geom_smooth(method = "lm", aes(group = TransZone), se = F) +
  xlab("Drug Usage") + ylab("Pfdhps 437G Prevalence")+ 
  facet_grid(~TransZone)

ggplot(Mod2_pfdhps, aes(cut(cumDrug,10),Prev_DR,col=Region))+geom_boxplot()+ xlab("Drug Usage") + 
  ylab("Pfdhps 437G Prevalence")+ facet_grid(~TransZone)

###
###
ggplot(dhps_na_rem, aes(x = mean_sp, y = Prev_DR)) + 
geom_point(aes(color=as.factor(Region), size=tested))+
  geom_smooth(method = "lm") + 
  xlab("SP usage") + ylab("Pfdhps 437G Prevalence")

ggplot(dhps_na_rem%>%filter(TransZone=="High"), aes(x = mean_sp, y = Prev_DR)) + 
  geom_point(aes(color=as.factor(Region), size=tested))+
  geom_smooth(method = "lm") + 
  xlab("SP usage") + ylab("Pfdhps 437G Prevalence")
ggsave("pfdhps_437G _p1A.png")

ggplot(dhps_na_rem%>%filter(TransZone=="Low"), aes(x = mean_sp, y = Prev_DR)) + 
  geom_point(aes(color=as.factor(Region), size=tested))+
  geom_smooth(method = "lm") + 
  xlab("SP usage") + ylab("Pfdhps 437G Prevalence")
ggsave("pfdhps_437G _p1B.png")


### 2. Plot of Pfdhps_437G  Prevalence vs PR_All regions
ggplot(dhps_na_rem, aes(cut(PR,10),Prev_DR,col=Region))+geom_boxplot()+ xlab("Parasite Rate") + ylab("Pfdhps 437G Prevalence")
ggplot(dhps_na_rem, aes(cut(PR,10),Prev_DR,col=Region))+geom_boxplot()+ xlab("Parasite Rate") + ylab("Pfdhps 437G Prevalence")+ facet_grid(~TransZone)
ggsave("pfdhps_437G _p2.png")

### 3. Plot of the relationship between Region and Year varied between pfdhps 437G Prevalence.
ggplot(data = dhps_na_rem %>% group_by(Region,year) %>% mutate(mean=mean(Prev_DR)), aes (y = Region, x = year,colour = mean)) + geom_point(aes(size=log(tested))) + scale_color_viridis_c(name="Pfdhps 437G Prevalence") + theme_classic() + xlab ("Year") + ylab ("Region")
ggplot(data = DHPS %>% group_by(Region,year) %>% mutate(mean=mean(Prev_DR)), aes (y = Region, x = year,colour = mean)) + geom_point(aes(size=log(tested))) + scale_color_viridis_c(name="Pfdhps 437G Prevalence") + theme_classic() + xlab ("Year") + ylab ("Region")+ ggtitle ("Time plot of Pfdhps_437G over time")

ggplot(data = DHPS %>% group_by(Region,year) %>% mutate(Region = factor(Region, levels = c("West Africa", "Central Africa" , "East Africa", "Southern Africa", "North Africa", "Caribbean", "Central America", "South America", "SEA", "Asia", "Oceania")), mean=mean(Prev_DR))) + aes (y = Region, x = year,colour = mean) + geom_point(aes(size=log(tested))) + scale_color_viridis_c(name="Pfdhps 437G Prevalence") + theme_classic() + xlab ("Year") + ylab ("Region")+ ggtitle ("Time plot of Pfdhps_437G over time")
ggsave("pfdhps 437G_Timep1.png")
ggsave("pfdhps 437G_Timep1A.png")

### 4. Plot of the relationship between Region and Year varied between Parasite rate.
ggplot(data = dhps_na_rem %>% group_by(Region,year) %>% mutate(mean=mean(PR)), aes (y = Region, x = year,colour = mean)) + geom_point(aes(size=log(tested))) + scale_color_viridis_c(name="PR") + theme_classic() + xlab ("Year") + ylab ("Region")
ggplot(data = DHPS %>% group_by(Region,year) %>% mutate(mean=mean(PR)), aes (y = Region, x = year,colour = mean)) + geom_point(aes(size=log(tested))) + scale_color_viridis_c(name="Parasite Rate") + theme_classic() + xlab ("Year") + ylab ("Region")+ ggtitle ("Time plot of PR over time")

ggplot(data = DHPS %>% group_by(Region,year) %>% mutate(Region = factor(Region, levels = c("West Africa", "Central Africa" , "East Africa", "Southern Africa", "North Africa", "Caribbean", "Central America", "South America", "SEA", "Asia", "Oceania")), mean=mean(PR))) + aes (y = Region, x = year,colour = mean) + geom_point(aes(size=log(tested))) + scale_color_viridis_c(name="Parasite Rate") + theme_classic() + xlab ("Year") + ylab ("Region")+ ggtitle ("Time plot of PR over time")
ggsave("pfdhps 437G_Timep2.png")
ggsave("pfdhps 437G_Timep2A.png")


### 5. Plot of the relationship between Region and Year varied between Mean_sp.
ggplot(data = dhps_na_rem %>% group_by(Region,year) %>% mutate(mean=mean(mean_sp)), aes (y = Region, x = year,colour = mean)) + geom_point(aes(size=log(tested))) + scale_color_viridis_c(name="SP usage") + theme_classic() + xlab ("Year") + ylab ("Region")+ ggtitle ("Time plot of SP usage over time")
#ggplot(data = dhps_na_rem %>% group_by(Region,year) %>% mutate(mean=mean(mean_sp)), aes (y = Region, x = year,colour = mean)) + geom_point(aes(size=log(tested))) + scale_color_gradient(low = "yellow", high = "red", name="Mean_sp") + theme_classic() + xlab ("Year") + ylab ("Region")

ggplot(data = dhps_na_rem %>% group_by(Region,year) %>% mutate(Region = factor(Region, levels = c("West Africa", "Central Africa" , "East Africa", "Southern Africa", "North Africa", "Caribbean", "Central America", "South America", "SEA", "Asia", "Oceania")), mean=mean(mean_sp))) + aes (y = Region, x = year,colour = mean) + geom_point(aes(size=log(tested))) + scale_color_viridis_c(name="SP usage") + theme_classic() + xlab ("Year") + ylab ("Region")+ ggtitle ("Time plot of SP usage over time")
ggsave("pfdhps 437G_Timep3.png")
ggsave("pfdhps 437G_Timep3A.png")


###PFDHFR S108N
#Prepare Database
# join original observation
DHFR<-read_csv ("pfdhfr_SNP.csv")
DHFR <- DHFR[!is.na(DHFR$PR),]
dhfr <- left_join(DHFR,sp, by = c("country","year"))
#dhfr<-read_csv ("pfdhfr_SNP.csv")%>%left_join(sp)
dhfr_na_rem <- dhfr %>% filter(!is.na(mean_sp))
#write_csv(dhfr_na_rem, file = "inspect.1.csv")
DHFR$Region[DHFR$country == "Senegal"] <-"West Africa"
DHFR$Region[DHFR$country == "Somalia"] <-"East Africa"
dhfr_na_rem$Region[dhfr_na_rem$country == "Senegal"] <-"West Africa"
dhfr_na_rem$Region[dhfr_na_rem$country == "Somalia"] <-"East Africa"
dhfr_na_rem$Region[dhfr_na_rem$country == "Madagascar"] <-"MAD"
dhfr_na_rem <- dhfr_na_rem%>% mutate(PR = PR / 100, Prev_DR = Prev_DR / 100)
write_csv(dhfr_na_rem, file = "pfdhfr_database.csv")

###Model
Mod1_pfdhfr <- lmer(Prev_DR ~ PR +mean_sp+PR *mean_sp+drug_1yr_ago+drug_2yr_ago+drug_3yr_ago+ (1|Region), weights=log(tested),data=dhfr_na_rem)

summary(Mod1_pfdhfr)
plot(Mod1_pfdhfr)
qqnorm(residuals(Mod1_pfdhfr))
###
###


### Recode regions

dhfr_na_rem<-dhfr_na_rem%>%mutate(TransZone = case_when(
  Region %in% c("West Africa","Central Africa") ~"High",
  Region %in% c("West Africa","Central Africa") ~"High",
  Region %in% c("East Africa","Southern Africa") ~"Med",
  Region %in% c("East Africa","Southern Africa") ~"Med",
  Region == "MAD" ~"MAD",
  Region %in% c("Asia", "SEA","South America")~"Low"))

###Model
Mod2_pfdhfr <- lmer(Prev_DR ~ PR +mean_sp+PR *mean_sp +drug_1yr_ago+drug_2yr_ago+drug_3yr_ago+(1|TransZone), 
                    weights=log(tested),data=dhfr_na_rem)
summary(Mod2_pfdhfr)
plot(Mod2_pfdhfr)
qqnorm(residuals(Mod2_pfdhfr))
###

Mod2_pfdhfr <- dhfr_na_rem%>%filter(!is.na(drug_3yr_ago))

Mod2_pfdhfr <-Mod2_pfdhfr%>%mutate(cumDrug = ifelse(is.na(drug_3yr_ago),drug_2yr_ago,drug_3yr_ago)+
                                     drug_2yr_ago+drug_1yr_ago)

Mod2_pfdhfr <- lmer(Prev_DR ~ PR +cumDrug+PR*cumDrug+(1|Region), 
                    weights=log(tested),data=Mod2_pfdhfr)
summary(Mod2_pfdhfr)
plot(Mod2_pfdhfr)
qqnorm(residuals(Mod2_pfdhfr))


#Modeling at High Transzone
highTrans <- dhfr_na_rem%>%filter(TransZone=="High")%>%filter(!is.na(drug_3yr_ago))

highTrans <-highTrans%>%mutate(cumDrug = ifelse(is.na(drug_3yr_ago),drug_2yr_ago,drug_3yr_ago)+
                                 drug_2yr_ago+drug_1yr_ago)

high_pfdhfr <- lmer(Prev_DR ~ PR +cumDrug+PR*cumDrug+ (1|Region), 
                   weights=log(tested),data=highTrans)
summary(high_pfdhfr)

high_pfdhfr  <- lm(Prev_DR ~ PR +cumDrug+PR*cumDrug, 
                 weights=log(tested),data=highTrans)
summary(high_pfdhfr)

#Modeling at Moderate Transzone
medTrans <- dhfr_na_rem%>%filter(TransZone=="Med")%>%filter(!is.na(drug_3yr_ago))

medTrans <-medTrans%>%mutate(cumDrug = ifelse(is.na(drug_3yr_ago),drug_2yr_ago,drug_3yr_ago)+
                                 drug_2yr_ago+drug_1yr_ago)

med_pfdhfr <- lmer(Prev_DR ~ PR +cumDrug+PR*cumDrug+ (1|Region), 
                    weights=log(tested),data=medTrans)
summary(med_pfdhfr)

med_pfdhfr  <- lm(Prev_DR ~ PR +cumDrug+PR*cumDrug, 
                   weights=log(tested),data=medTrans)
summary(med_pfdhfr)

#Modeling at Low Transzone
lowTrans <- dhfr_na_rem%>%filter(TransZone=="Low")%>%filter(!is.na(drug_3yr_ago))

lowTrans <-lowTrans%>%mutate(cumDrug = ifelse(is.na(drug_3yr_ago),drug_2yr_ago,drug_3yr_ago)+
                                 drug_2yr_ago+drug_1yr_ago)

low_pfdhfr <- lmer(Prev_DR ~ PR +cumDrug+PR*cumDrug+ (1|Region), 
                    weights=log(tested),data=lowTrans)
summary(low_pfdhfr)

low_pfdhfr  <- lm(Prev_DR ~ PR +cumDrug+PR*cumDrug, 
                   weights=log(tested),data=lowTrans)
summary(low_pfdhfr)



#Interpretation:
#Like in the case of pfdhps-A437G prevalence, the effect of SP use on the prevalence of pfdhfr-S108N mutation over 3 years lag period suggested that SP use negatively influenced pfdhfr-S108N mutation in 3years, but this was not found to be significant. On the other hand, Parasite rate, and the interaction between Parasite rate and SP use were significantly (p<0.001 and p<0.05) associated with pfdhfr-S108N prevalence, but with negative and positive estimates, respectively. 


###Plots
dhfr_na_rem <- dhfr_na_rem%>% mutate(PR = PR * 100, Prev_DR = Prev_DR * 100)
#### 1. Plots of Pfdhfr_108N Prevalence vs Mean_sp_All regions
ggplot(dhfr_na_rem, aes(cut(mean_sp,10),Prev_DR,col=Region))+geom_boxplot()+ xlab("Drug Usage") + ylab("Pfdhfr 108N Prevalence")
ggplot(dhfr_na_rem, aes(cut(mean_sp,10),Prev_DR,col=Region))+geom_boxplot()+ xlab("Drug Usage") + ylab("Pfdhfr 108N Prevalence")+ facet_grid(~TransZone)
ggsave("pfdhfr_108N_p1.png") 

### 2. Plots of Pfdhfr_108N Prevalence vs PR_All regions
ggplot(dhfr_na_rem, aes(cut(PR,10),Prev_DR,col=Region))+geom_boxplot()+ xlab("Parasite Rate") + ylab("Pfdhfr 108N Prevalence")
ggplot(dhfr_na_rem, aes(cut(PR,10),Prev_DR,col=Region))+geom_boxplot()+ xlab("Parasite Rate") + ylab("Pfdhfr 108N Prevalence")+ facet_grid(~TransZone)
ggsave("pfdhfr_108N_p2.png")

### 3. Plot of the relationship between Region and Year varied between pfdhps 108N Prevalence.
ggplot(data = dhfr_na_rem %>% group_by(Region,year) %>% mutate(mean=mean(Prev_DR)), aes (y = Region, x = year,colour = mean)) + geom_point(aes(size=log(tested))) + scale_color_viridis_c(name="Pfdhfr 108N Prevalence") + theme_classic() + xlab ("Year") + ylab ("Region")
ggplot(data = DHFR %>% group_by(Region,year) %>% mutate(mean=mean(Prev_DR)), aes (y = Region, x = year,colour = mean)) + geom_point(aes(size=log(tested))) + scale_color_viridis_c(name="Pfdhfr 108N Prevalence") + theme_classic() + xlab ("Year") + ylab ("Region")+ ggtitle ("Time plot of Pfdhfr_108N over time")

ggplot(data = DHFR %>% group_by(Region,year) %>% mutate(Region = factor(Region, levels = c("West Africa", "Central Africa" , "East Africa", "Southern Africa", "North Africa", "Caribbean", "Central America", "South America", "SEA", "Asia", "Oceania")), mean=mean(Prev_DR))) + aes (y = Region, x = year,colour = mean) + geom_point(aes(size=log(tested))) + scale_color_viridis_c(name="Pfdhfr 108N Prevalence") + theme_classic() + xlab ("Year") + ylab ("Region")+ ggtitle ("Time plot of Pfdhfr_108N over time")
ggsave("pfdhfr 108N_Timep1.png")
ggsave("pfdhfr 108N_Timep1A.png")


### 4. Plot of the relationship between Region and Year varied between Parasite rate.
ggplot(data = dhfr_na_rem %>% group_by(Region,year) %>% mutate(mean=mean(PR)), aes (y = Region, x = year,colour = mean)) + geom_point(aes(size=log(tested))) + scale_color_viridis_c(name="PR") + theme_classic() + xlab ("Year") + ylab ("Region")
ggplot(data = DHFR %>% group_by(Region,year) %>% mutate(mean=mean(PR)), aes (y = Region, x = year,colour = mean)) + geom_point(aes(size=log(tested))) + scale_color_viridis_c(name="Parasite Rate") + theme_classic() + xlab ("Year") + ylab ("Region")+ ggtitle ("Time plot of PR over time")

ggplot(data = DHFR %>% group_by(Region,year) %>% mutate(Region = factor(Region, levels = c("West Africa", "Central Africa" , "East Africa", "Southern Africa", "North Africa", "Caribbean", "Central America", "South America", "SEA", "Asia", "Oceania")), mean=mean(PR))) + aes (y = Region, x = year,colour = mean) + geom_point(aes(size=log(tested))) + scale_color_viridis_c(name="Parasite Rate") + theme_classic() + xlab ("Year") + ylab ("Region")+ ggtitle ("Time plot of PR over time")
ggsave("pfdhfr 108N_Timep2.png")
ggsave("pfdhfr 108N_Timep2A.png")

### 5. Plot of the relationship between Region and Year varied between Mean_sp.
ggplot(data = dhfr_na_rem %>% group_by(Region,year) %>% mutate(mean=mean(mean_sp)), aes (y = Region, x = year,colour = mean)) + geom_point(aes(size=log(tested))) + scale_color_viridis_c(name="SP usage") + theme_classic() + xlab ("Year") + ylab ("Region")+ ggtitle ("Time plot of SP usage over time")

ggplot(data = dhfr_na_rem %>% group_by(Region,year) %>% mutate(Region = factor(Region, levels = c("West Africa", "Central Africa" , "East Africa", "Southern Africa", "North Africa", "Caribbean", "Central America", "South America", "SEA", "Asia", "Oceania")), mean=mean(mean_sp))) + aes (y = Region, x = year,colour = mean) + geom_point(aes(size=log(tested))) + scale_color_viridis_c(name="SP usage") + theme_classic() + xlab ("Year") + ylab ("Region")+ ggtitle ("Time plot of SP usage over time")
ggsave("pfdhfr 108N_Timep3.png")
ggsave("pfdhfr 108N_Timep3A.png")


###PFK13 C580Y
act<- read_csv("ACT_imputed.csv")

#Set up lag years
act<-act%>%group_by(country)%>%mutate(drug_1yr_ago = lag(mean_act,n=1),drug_2yr_ago = lag(mean_act,n=2), drug_3yr_ago = lag(mean_act,n=3),drug_4yr_ago = lag(mean_act,n=4),drug_5yr_ago = lag(mean_act,n=5))

# join original observation
k13<-read_csv ("pfK13_SNP.csv")%>%left_join(act)
k13_na_rem <- k13 %>% filter(!is.na(mean_act))
#write_csv(k13_na_rem, file = "inspect.3.csv")
k13_na_rem<- k13_na_rem %>% mutate(PR = PR / 100, Prev_DR = Prev_DR / 100)
write_csv(k13_na_rem, file = "pfk13_database.csv")

###Model
k13_na_rem <-k13_na_rem %>% filter(Prev_DR<1,Prev_DR>0)
Mod1_pfk13 <- lmer(Prev_DR ~ PR +mean_act+PR *mean_act +drug_1yr_ago+drug_2yr_ago+drug_3yr_ago+ (1|country), weights=log(tested),data=k13_na_rem)
summary(Mod1_pfk13)
plot(Mod1_pfk13)
qqnorm(residuals(Mod1_pfk13))
###
###

### Recode regions

k13_na_rem <-k13_na_rem %>%mutate(TransZone = case_when(
  Region %in% c("West Africa","Central Africa") ~"High",
  Region %in% c("West Africa","Central Africa") ~"High",
  Region %in% c("East Africa","Southern Africa") ~"Med",
  Region %in% c("East Africa","Southern Africa") ~"Med",
  Region == "MAD" ~"MAD",
  Region %in% c("Asia", "SEA","South America")~"Low"))

###Model
Mod2_pfk13 <- lmer(Prev_DR ~ PR +mean_act+PR *mean_act +drug_1yr_ago+drug_2yr_ago+drug_3yr_ago+(1|country), 
                    weights=log(tested),data=k13_na_rem)
summary(Mod2_pfk13)
plot(Mod2_pfk13)
qqnorm(residuals(Mod2_pfk13))


#Modeling at Low Transzone
lowTrans <- k13_na_rem%>%filter(TransZone=="Low")%>%filter(!is.na(drug_3yr_ago))

lowTrans <-lowTrans%>%mutate(cumDrug = ifelse(is.na(drug_3yr_ago),drug_2yr_ago,drug_3yr_ago)+
                               drug_2yr_ago+drug_1yr_ago)

low_pfk13 <- lmer(Prev_DR ~ PR +cumDrug+PR*cumDrug+ (1|country), 
                   weights=log(tested),data=lowTrans)
summary(low_pfk13)

low_pfk13  <- lm(Prev_DR ~ PR +cumDrug+PR*cumDrug, 
                  weights=log(tested),data=lowTrans)
summary(low_pfk13)


#Interpretation:
#ACT use was found to positively influence pfk13-C580Y prevalence at lead times of up to 3 years, although not significantly. However, Parasite rate negatively correlated with pfk13-C580Y mutation and this was also not significant. The interaction between Parasite rate and ACT use did not significantly influence pfk13-C580Y prevalence either. 
ggplot(k13_na_rem%>%filter(tested>20), aes(drug_3yr_ago+drug_2yr_ago,Prev_DR,col=Region))+
  geom_point(aes(color=site,size=tested))+
  #geom_line(aes(y=pred))+
  facet_grid(~country,scales="free_x")

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



#####
####Maps
# work with maps package
install.packages("maps")
library(maps)
library(tidyverse)

###PFCRT K76T
#dataset = read_csv("pfcrt_SNP.csv")
#dataset <-crt_na_rem
dataset <-CRT
## Map of Pfcrt_76T Prevalence
ggplot() +
  # Draw the world map
  borders("world", colour = "gray85", fill = "gray80") +
  # Plot points using lat and long
  geom_point(data = CRT, aes(x = long, y = lat, color = Prev_DR), size = 2, alpha = 0.7) +
  # Customize the appearance
  theme_minimal() +
  labs(title = "Prevalence of Pfcrt_76T", x = "Longitude", y = "Latitude", color = "Pfcrt_76T") +
  theme(legend.position = "right")+
  scale_color_viridis_c()

## Map of PR
ggplot() +
  # Draw the world map
  borders("world", colour = "gray85", fill = "gray80") +
  # Plot points using lat and long
  geom_point(data = CRT, aes(x = long, y = lat, color = PR), size = 2, alpha = 0.7) +
  # Customize the appearance
  theme_minimal() +
  labs(title = "Prevalence of Pfcrt_76T", x = "Longitude", y = "Latitude", color = "Parasite Rate") +
  theme(legend.position = "right")+
  scale_color_viridis_c()

## Map of CQ_use
ggplot() +
  # Draw the world map
  borders("world", colour = "gray85", fill = "gray80") +
  # Plot points using lat and long
  geom_point(data = crt_na_rem, aes(x = long, y = lat, color = mean_cq), size = 2, alpha = 0.7) +
  # Customize the appearance
  theme_minimal() +
  labs(title = "Prevalence of Pfcrt_76T", x = "Longitude", y = "Latitude", color = "CQ usage") +
  theme(legend.position = "right")+
  scale_color_viridis_c()
####
####

## Map of Pfcrt_76T Prevalence
world_map <- map_data("world")
# Create a ggplot2 map
ggplot() +
  geom_polygon(data = world_map, aes(x = long, y = lat, group = group), 
               fill = "lightgray", color = "white") +
  geom_point(data = CRT, aes(x = long, y = lat, color = Prev_DR), 
             size = 3, alpha = 0.7) +
  scale_color_gradient(low = "yellow", high = "red") +
  theme_minimal() +
  labs(title = "Prevalence of Pfcrt_76T ", x = "Longitude", y = "Latitude", 
       color = "Pfcrt_76T")

## Map of PR
world_map <- map_data("world")
# Create a ggplot2 map
ggplot() +
  geom_polygon(data = world_map, aes(x = long, y = lat, group = group), 
               fill = "lightgray", color = "white") +
  geom_point(data = CRT, aes(x = long, y = lat, color = PR), 
             size = 3, alpha = 0.7) +
  scale_color_gradient(low = "yellow", high = "red") +
  theme_minimal() +
  labs(title = "Prevalence of Pfcrt_76T ", x = "Longitude", y = "Latitude", 
       color = "Parasite Rate")

## Map of CQ_use
world_map <- map_data("world")
# Create a ggplot2 map
ggplot() +
  geom_polygon(data = world_map, aes(x = long, y = lat, group = group), 
               fill = "lightgray", color = "white") +
  geom_point(data = crt_na_rem, aes(x = long, y = lat, color = mean_cq), 
             size = 3, alpha = 0.7) +
  scale_color_gradient(low = "yellow", high = "red") +
  theme_minimal() +
  labs(title = "Prevalence of Pfcrt_76T ", x = "Longitude", y = "Latitude", 
       color = "CQ usage")

##########
##########
## Map of Pfcrt_76T Prevalence
world_map <- map_data("world")
# Create a ggplot2 map
ggplot() +
  geom_polygon(data = world_map, aes(x = long, y = lat, group = group), 
               fill = "lightgray", color = "white") +
  geom_point(data = crt_na_rem, aes(x = long, y = lat, color = Prev_DR), 
             size = 3, alpha = 0.7) +
  scale_color_viridis_c()+
  theme_minimal() +
  labs(title = "Global distribution of Pfcrt_76T ", x = "Longitude", y = "Latitude", 
       color = "Pfcrt_76T")
ggsave("pfcrt 76T_Map1.png")

## Map of PR
world_map <- map_data("world")
# Create a ggplot2 map
ggplot() +
  geom_polygon(data = world_map, aes(x = long, y = lat, group = group), 
               fill = "lightgray", color = "white") +
  geom_point(data = crt_na_rem, aes(x = long, y = lat, color = PR), 
             size = 3, alpha = 0.7) +
  scale_color_viridis_c() +
  theme_minimal() +
  labs(title = "Global distribution of PR", x = "Longitude", y = "Latitude", 
       color = "Parasite Rate")
ggsave("pfcrt 76T_Map2.png")

## Map of CQ_use
world_map <- map_data("world")
# Create a ggplot2 map
ggplot() +
  geom_polygon(data = world_map, aes(x = long, y = lat, group = group), 
               fill = "lightgray", color = "white") +
  geom_point(data = crt_na_rem, aes(x = long, y = lat, color = mean_cq), 
             size = 3, alpha = 0.7) +
  scale_color_viridis_c() +
  theme_minimal() +
  labs(title = "Global distribution of CQ usage ", x = "Longitude", y = "Latitude", 
       color = "CQ usage")
ggsave("pfcrt 76T_Map3.png")


###PFMDR1 N86Y
#dataset = read_csv("pfmdr1_SNP.csv")
dataset <-mdr_na_rem
## Map of Pfmdr1_86Y Prevalence
ggplot() +
  # Draw the world map
  borders("world", colour = "gray85", fill = "gray80") +
  # Plot points using lat and long
  geom_point(data = dataset, aes(x = long, y = lat, color = Prev_DR), size = 2, alpha = 0.7) +
  # Customize the appearance
  theme_minimal() +
  labs(title = "Data Points on World Map", x = "Longitude", y = "Latitude", color = "Pfmdr1_86Y") +
  theme(legend.position = "right")+
  scale_color_viridis_c()

## Map of PR
ggplot() +
  # Draw the world map
  borders("world", colour = "gray85", fill = "gray80") +
  # Plot points using lat and long
  geom_point(data = dataset, aes(x = long, y = lat, color = PR), size = 2, alpha = 0.7) +
  # Customize the appearance
  theme_minimal() +
  labs(title = "Data Points on World Map", x = "Longitude", y = "Latitude", color = "Parasite Rate") +
  theme(legend.position = "right")+
  scale_color_viridis_c()

## Map of CQ_use
ggplot() +
  # Draw the world map
  borders("world", colour = "gray85", fill = "gray80") +
  # Plot points using lat and long
  geom_point(data = dataset, aes(x = long, y = lat, color = mean_cq), size = 2, alpha = 0.7) +
  # Customize the appearance
  theme_minimal() +
  labs(title = "Data Points on World Map", x = "Longitude", y = "Latitude", color = "CQ usage") +
  theme(legend.position = "right")+
  scale_color_viridis_c()
####
####

###PFDHPS A437G
#dataset = read_csv("pfdhps_SNP.csv")
dataset <-dhps_na_rem
## Map of Pfdhps_437G Prevalence
ggplot() +
  # Draw the world map
  borders("world", colour = "gray85", fill = "gray80") +
  # Plot points using lat and long
  geom_point(data = dataset, aes(x = long, y = lat, color = Prev_DR), size = 2, alpha = 0.7) +
  # Customize the appearance
  theme_minimal() +
  labs(title = "Data Points on World Map", x = "Longitude", y = "Latitude", color = "Pfdhps_437G") +
  theme(legend.position = "right")+
  scale_color_viridis_c()

## Map of PR
ggplot() +
  # Draw the world map
  borders("world", colour = "gray85", fill = "gray80") +
  # Plot points using lat and long
  geom_point(data = dataset, aes(x = long, y = lat, color = PR), size = 2, alpha = 0.7) +
  # Customize the appearance
  theme_minimal() +
  labs(title = "Data Points on World Map", x = "Longitude", y = "Latitude", color = "Parasite Rate") +
  theme(legend.position = "right")+
  scale_color_viridis_c()

## Map of SP_use
ggplot() +
  # Draw the world map
  borders("world", colour = "gray85", fill = "gray80") +
  # Plot points using lat and long
  geom_point(data = dataset, aes(x = long, y = lat, color = mean_sp), size = 2, alpha = 0.7) +
  # Customize the appearance
  theme_minimal() +
  labs(title = "Data Points on World Map", x = "Longitude", y = "Latitude", color = "SP usage") +
  theme(legend.position = "right")+
  scale_color_viridis_c()

##########
##########
## Map of Pfdhps_437G Prevalence
world_map <- map_data("world")
# Create a ggplot2 map
ggplot() +
  geom_polygon(data = world_map, aes(x = long, y = lat, group = group), 
               fill = "lightgray", color = "white") +
  geom_point(data = dhps_na_rem, aes(x = long, y = lat, color = Prev_DR), 
             size = 3, alpha = 0.7) +
  scale_color_viridis_c()+
  theme_minimal() +
  labs(title = "Global distribution of Pfdhps_437G ", x = "Longitude", y = "Latitude", 
       color = "Pfdhps_437G")
ggsave("pfdhps 437G_Map1.png")

## Map of PR
world_map <- map_data("world")
# Create a ggplot2 map
ggplot() +
  geom_polygon(data = world_map, aes(x = long, y = lat, group = group), 
               fill = "lightgray", color = "white") +
  geom_point(data = dhps_na_rem, aes(x = long, y = lat, color = PR), 
             size = 3, alpha = 0.7) +
  scale_color_viridis_c() +
  theme_minimal() +
  labs(title = "Global distribution of PR", x = "Longitude", y = "Latitude", 
       color = "Parasite Rate")
ggsave("pfdhps 437G_Map2.png")

## Map of SP_use
world_map <- map_data("world")
# Create a ggplot2 map
ggplot() +
  geom_polygon(data = world_map, aes(x = long, y = lat, group = group), 
               fill = "lightgray", color = "white") +
  geom_point(data = dhps_na_rem, aes(x = long, y = lat, color = mean_sp), 
             size = 3, alpha = 0.7) +
  scale_color_viridis_c() +
  theme_minimal() +
  labs(title = "Global distribution of SP usage ", x = "Longitude", y = "Latitude", 
       color = "SP usage")
ggsave("pfdhps 437G_Map3.png")


####
####
###PFDHFR S108N
#dataset = read_csv("pfdhfr_SNP.csv")
dataset <-dhfr_na_rem
## Map of Pfdhfr_108N Prevalence
ggplot() +
  # Draw the world map
  borders("world", colour = "gray85", fill = "gray80") +
  # Plot points using lat and long
  geom_point(data = dataset, aes(x = long, y = lat, color = Prev_DR), size = 2, alpha = 0.7) +
  # Customize the appearance
  theme_minimal() +
  labs(title = "Data Points on World Map", x = "Longitude", y = "Latitude", color = "Pfdhfr_108N") +
  theme(legend.position = "right")+
  scale_color_viridis_c()

## Map of PR
ggplot() +
  # Draw the world map
  borders("world", colour = "gray85", fill = "gray80") +
  # Plot points using lat and long
  geom_point(data = dataset, aes(x = long, y = lat, color = PR), size = 2, alpha = 0.7) +
  # Customize the appearance
  theme_minimal() +
  labs(title = "Data Points on World Map", x = "Longitude", y = "Latitude", color = "Parasite Rate") +
  theme(legend.position = "right")+
  scale_color_viridis_c()

## Map of SP_use
ggplot() +
  # Draw the world map
  borders("world", colour = "gray85", fill = "gray80") +
  # Plot points using lat and long
  geom_point(data = dataset, aes(x = long, y = lat, color = mean_sp), size = 2, alpha = 0.7) +
  # Customize the appearance
  theme_minimal() +
  labs(title = "Data Points on World Map", x = "Longitude", y = "Latitude", color = "SP usage") +
  theme(legend.position = "right")+
  scale_color_viridis_c()

##########
##########
## Map of Pfdhfr_108N  Prevalence
world_map <- map_data("world")
# Create a ggplot2 map
ggplot() +
  geom_polygon(data = world_map, aes(x = long, y = lat, group = group), 
               fill = "lightgray", color = "white") +
  geom_point(data = dhfr_na_rem, aes(x = long, y = lat, color = Prev_DR), 
             size = 3, alpha = 0.7) +
  scale_color_viridis_c()+
  theme_minimal() +
  labs(title = "Global distribution of Pfdhfr_108N ", x = "Longitude", y = "Latitude", 
       color = "Pfdhfr_108N")
ggsave("pfdhfr 108N_Map1.png")

## Map of PR
world_map <- map_data("world")
# Create a ggplot2 map
ggplot() +
  geom_polygon(data = world_map, aes(x = long, y = lat, group = group), 
               fill = "lightgray", color = "white") +
  geom_point(data = dhfr_na_rem, aes(x = long, y = lat, color = PR), 
             size = 3, alpha = 0.7) +
  scale_color_viridis_c() +
  theme_minimal() +
  labs(title = "Global distribution of PR", x = "Longitude", y = "Latitude", 
       color = "Parasite Rate")
ggsave("pfdhfr 108N_Map2.png")

## Map of CQ_use
world_map <- map_data("world")
# Create a ggplot2 map
ggplot() +
  geom_polygon(data = world_map, aes(x = long, y = lat, group = group), 
               fill = "lightgray", color = "white") +
  geom_point(data = dhfr_na_rem, aes(x = long, y = lat, color = mean_sp), 
             size = 3, alpha = 0.7) +
  scale_color_viridis_c() +
  theme_minimal() +
  labs(title = "Global distribution of SP usage ", x = "Longitude", y = "Latitude", 
       color = "SP usage")
ggsave("pfdhfr 108N_Map3.png")


####
####

###PFK13 C580Y
#dataset = read_csv("pfk13_SNP.csv")
dataset <-k13_na_rem
## Map of Pfk13_580Y Prevalence
ggplot() +
  # Draw the world map
  borders("world", colour = "gray85", fill = "gray80") +
  # Plot points using lat and long
  geom_point(data = dataset, aes(x = long, y = lat, color = Prev_DR), size = 2, alpha = 0.7) +
  # Customize the appearance
  theme_minimal() +
  labs(title = "Data Points on World Map", x = "Longitude", y = "Latitude", color = "Pfk13_580Y") +
  theme(legend.position = "right")+
  scale_color_viridis_c()

## Map of PR
ggplot() +
  # Draw the world map
  borders("world", colour = "gray85", fill = "gray80") +
  # Plot points using lat and long
  geom_point(data = dataset, aes(x = long, y = lat, color = PR), size = 2, alpha = 0.7) +
  # Customize the appearance
  theme_minimal() +
  labs(title = "Data Points on World Map", x = "Longitude", y = "Latitude", color = "Parasite Rate") +
  theme(legend.position = "right")+
  scale_color_viridis_c()

## Map of ACT_use
ggplot() +
  # Draw the world map
  borders("world", colour = "gray85", fill = "gray80") +
  # Plot points using lat and long
  geom_point(data = dataset, aes(x = long, y = lat, color = mean_act), size = 2, alpha = 0.7) +
  # Customize the appearance
  theme_minimal() +
  labs(title = "Data Points on World Map", x = "Longitude", y = "Latitude", color = "ACT usage") +
  theme(legend.position = "right")+
  scale_color_viridis_c()
####
####

