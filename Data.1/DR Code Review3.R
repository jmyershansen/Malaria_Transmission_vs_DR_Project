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
library(boot)

##CRT
cq<- read_csv("CQ_imputed.csv")
crt<-read_csv ("pfcrt_SNP.csv")

pfcrt_complete <- left_join(crt, cq, by = c("country","year"))
pfcrt_complete <- pfcrt_complete %>% arrange(country, year)
write_csv(pfcrt_complete, file = "pfcrt_prev&pr&cq.1.csv")

#Load data 
DU <- read_csv("pfcrt_prev&pr&cq.1.csv")
#To interpolate the missing drug use data, cut out "country, year and drug use" columns 
DU <- as.data.frame(cbind(DU$country,DU$year,DU$mean_cq))
colnames(DU)<-c("country","year","CQ")
DU <- DU %>% mutate(CQ = as.numeric(CQ))
DU <- DU %>% mutate(year = as.numeric(year))
#Summarize drug use data from the same 
DU <- DU %>% group_by(country, year) %>% summarise(mean_cq=mean(CQ, na.rm = T))
DU = DU %>% arrange(country, year)
#Expand drug use year from </= 2000 to 2022
c = DU$country %>% unique()
# loop through and find min and max year, if any years are missing for a country add that row with NAs
for(i in c){
  ro <- which(DU$country == i)
  mi <- min(DU$year[ro])
  mi<-ifelse(mi>2000,2000,mi)
  ma <- max(DU$year[ro])
  ma <- ifelse(ma<2022,2022,ma)
  my_range = mi:ma
  x <- which(!my_range %in% DU$year[ro])
  if(length(x)>0){
    new_ros = DU[1:length(x),]
    new_ros$country = i
    new_ros$year = my_range[x]
    new_ros[,3:ncol(new_ros)] = NA
    DU <- rbind(DU, new_ros)
  }
}
DU = DU %>% arrange(country, year)

# Say if data was observed or interpolated into column "du_obs" ie. NAs (data to be interpolated)=0, Observed data=1
DU$du_obs = ifelse(is.na(DU$mean_cq), 0, 1)
#Summarize observed data into column "num_obs"
DU<- DU %>% group_by(country) %>% mutate(num_obs=sum(du_obs))
#Filter out country with 3 or more data points for interpolation
c = DU %>% filter(num_obs>2) %>% 
  select(country) %>% unlist() %>% unique()
#Interpolate data
for(i in c){
  print(i)
  ro = which(DU$country == i)
  du = DU$mean_cq[ro]
  DU$mean_cq[ro] = na_kalman(du)
}
DU <- DU %>% filter(num_obs>2)
# correct for above 1 or below 0
DU$mean_cq = ifelse(DU$mean_cq > 1, 1, DU$mean_cq)
DU$mean_cq = ifelse(DU$mean_cq < 0, 0, DU$mean_cq)
#Filter out data from 2000 and above
DU <- DU %>% filter(year>1999)

cq1<- DU
crt<-read_csv ("pfcrt_SNP.csv")

#Set up lag years
cq1<-cq1%>%group_by(country)%>%mutate(drug_1yr_ago = lag(mean_cq,n=1),drug_2yr_ago = lag(mean_cq,n=2),
                                      drug_3yr_ago = lag(mean_cq,n=3),drug_4yr_ago = lag(mean_cq,n=4),drug_5yr_ago = lag(mean_cq,n=5))
write_csv(cq1, file = "CQ_lag.csv")
cq1 <- cq1 %>% mutate(year = as.numeric(year))
# Full join main data to drug use data
pfcrt_complete <- full_join(crt, cq1, by = c("country","year"))
pfcrt_complete <- pfcrt_complete %>% arrange(country, year)
write_csv(pfcrt_complete, file = "pfcrt_prev&pr&cq.2.csv")
#Load data
pfcrt_complete<-read_csv ("pfcrt_prev&pr&cq.2.csv")
pfcrt_complete2 <- pfcrt_complete$country %>% unique()
#Extract out NAs from drug use column and have a good lag year structure. NB.This is the same as right joining drug use data to main data. 
pfcrt_complete2 <- pfcrt_complete[!is.na(pfcrt_complete$mean_cq),]
#pfcrt_complete2 <- sum(duplicated(pfcrt_complete2))
pfcrt_complete2 <- pfcrt_complete2[!duplicated(pfcrt_complete2),]
write_csv(pfcrt_complete2, file = "pfcrt_prev&pr&cq.3.csv")


##MDR1
cq<- read_csv("CQ_imputed.csv")
mdr<-read_csv ("pfmdr1_SNP.csv")

pfmdr1_complete <- left_join(mdr, cq, by = c("country","year"))
pfmdr1_complete <- pfmdr1_complete %>% arrange(country, year)
write_csv(pfmdr1_complete, file = "pfmdr1_prev&pr&cq.1.csv")


#Load data 
DU <- read_csv("pfmdr1_prev&pr&cq.1.csv")
#To interpolate the missing drug use data, cut out "country, year and drug use" columns  
DU <- as.data.frame(cbind(DU$country,DU$year,DU$mean_cq))
colnames(DU)<-c("country","year","CQ")
DU <- DU %>% mutate(CQ = as.numeric(CQ))
DU <- DU %>% mutate(year = as.numeric(year))
#Summarize drug use data from the same 
DU <- DU %>% group_by(country, year) %>% summarise(mean_cq=mean(CQ, na.rm = T))
DU <- DU %>% arrange(country, year)
#Expand drug use year from </= 2000 to 2022
c <- DU$country %>% unique()
# loop through and find min and max year, if any years are missing for a country add that row with NAs
for(i in c){
  ro <- which(DU$country == i)
  mi <- min(DU$year[ro])
  mi<- ifelse(mi>2000,2000,mi)
  ma <- max(DU$year[ro])
  ma <- ifelse(ma<2022,2022,ma)
  my_range <- mi:ma
  x <- which(!my_range %in% DU$year[ro])
  if(length(x)>0){
    new_ros <- DU[1:length(x),]
    new_ros$country = i
    new_ros$year = my_range[x]
    new_ros[,3:ncol(new_ros)] = NA
    DU <- rbind(DU, new_ros)
  }
}
DU <- DU %>% arrange(country, year)

# Say if data was observed or interpolated into column "du_obs" ie. NAs (data to be interpolated)=0, Observed data=1
DU$du_obs <- ifelse(is.na(DU$mean_cq), 0, 1)
#Summarize observed data into column "num_obs"
DU<- DU %>% group_by(country) %>% mutate(num_obs=sum(du_obs))
#Filter out country with 3 or more data points for interpolation
c <- DU %>% filter(num_obs>2) %>% 
  select(country) %>% unlist() %>% unique()
#Interpolate data
for(i in c){
  print(i)
  ro <- which(DU$country == i)
  du <- DU$mean_cq[ro]
  DU$mean_cq[ro] <- na_kalman(du)
}
DU <- DU %>% filter(num_obs>2)
# correct for above 1 or below 0
DU$mean_cq <- ifelse(DU$mean_cq > 1, 1, DU$mean_cq)
DU$mean_cq <- ifelse(DU$mean_cq < 0, 0, DU$mean_cq)
#Filter out data from 2000 and above
DU <- DU %>% filter(year>1999)

cq1<- DU
mdr<-read_csv ("pfmdr1_SNP.csv")

#Set up lag years
cq1<-cq1%>%group_by(country)%>%mutate(drug_1yr_ago = lag(mean_cq,n=1),drug_2yr_ago = lag(mean_cq,n=2),
                                      drug_3yr_ago = lag(mean_cq,n=3),drug_4yr_ago = lag(mean_cq,n=4),drug_5yr_ago = lag(mean_cq,n=5))
write_csv(cq1, file = "CQ_lag.csv")
cq1 <- cq1 %>% mutate(year = as.numeric(year))
# Full join main data to drug use data
pfmdr1_complete <- full_join(mdr, cq1, by = c("country","year"))
pfmdr1_complete <- pfmdr1_complete %>% arrange(country, year)
write_csv(pfmdr1_complete, file = "pfmdr1_prev&pr&cq.2.csv")
#Load data
pfmdr1_complete<-read_csv ("pfmdr1_prev&pr&cq.2.csv")
pfmdr1_complete2 <- pfmdr1_complete$country %>% unique()
#Extract out NAs from drug use column and have a good lag year structure. NB.This is the same as right joining drug use data to main data (ie. mdr and cq1). 
pfmdr1_complete2 <- pfmdr1_complete[!is.na(pfmdr1_complete$mean_cq),]
#pfmdr1_complete2 <- sum(duplicated(pfmdr1_complete2))
pfmdr1_complete2 <- pfmdr1_complete2[!duplicated(pfmdr1_complete2),]
write_csv(pfmdr1_complete2, file = "pfmdr1_prev&pr&cq.3.csv")


##DHPS
sp<- read_csv("SP_imputed.csv")
dhps<-read_csv ("pfdhps_SNP.csv")

pfdhps_complete <- left_join(dhps, sp, by = c("country","year"))
pfdhps_complete <- pfdhps_complete %>% arrange(country, year)
write_csv(pfdhps_complete, file = "pfdhps_prev&pr&sp.1.csv")

#Load data 
DU <- read_csv("pfdhps_prev&pr&sp.1.csv")
#To interpolate the missing drug use data, cut out "country, year and drug use" columns  
DU <- as.data.frame(cbind(DU$country,DU$year,DU$mean_sp))
colnames(DU)<-c("country","year","SP")
DU <- DU %>% mutate(SP = as.numeric(SP))
DU <- DU %>% mutate(year = as.numeric(year))
#Summarize drug use data from the same 
DU <- DU %>% group_by(country, year) %>% summarise(mean_sp=mean(SP, na.rm = T))
DU <- DU %>% arrange(country, year)
#Expand drug use year from </= 2000 to 2022
c <- DU$country %>% unique()
# loop through and find min and max year, if any years are missing for a country add that row with NAs
for(i in c){
  ro <- which(DU$country == i)
  mi <- min(DU$year[ro])
  mi<- ifelse(mi>2000,2000,mi)
  ma <- max(DU$year[ro])
  ma <- ifelse(ma<2022,2022,ma)
  my_range <- mi:ma
  x <- which(!my_range %in% DU$year[ro])
  if(length(x)>0){
    new_ros <- DU[1:length(x),]
    new_ros$country = i
    new_ros$year = my_range[x]
    new_ros[,3:ncol(new_ros)] = NA
    DU <- rbind(DU, new_ros)
  }
}
DU <- DU %>% arrange(country, year)

# Say if data was observed or interpolated into column "du_obs" ie. NAs (data to be interpolated)=0, Observed data=1
DU$du_obs <- ifelse(is.na(DU$mean_sp), 0, 1)
#Summarize observed data into column "num_obs"
DU<- DU %>% group_by(country) %>% mutate(num_obs=sum(du_obs))
#Filter out country with 3 or more data points for interpolation
c <- DU %>% filter(num_obs>2) %>% 
  select(country) %>% unlist() %>% unique()
#Interpolate data
for(i in c){
  print(i)
  ro <- which(DU$country == i)
  du <- DU$mean_sp[ro]
  DU$mean_sp[ro] <- na_kalman(du)
}
DU <- DU %>% filter(num_obs>2)
# correct for above 1 or below 0
DU$mean_sp <- ifelse(DU$mean_sp > 1, 1, DU$mean_sp)
DU$mean_sp <- ifelse(DU$mean_sp < 0, 0, DU$mean_sp)
#Filter out data from 2000 and above
DU <- DU %>% filter(year>1999)

sp1<- DU
dhps<-read_csv ("pfdhps_SNP.csv")

#Set up lag years
sp1<-sp1%>%group_by(country)%>%mutate(drug_1yr_ago = lag(mean_sp,n=1),drug_2yr_ago = lag(mean_sp,n=2),
                                      drug_3yr_ago = lag(mean_sp,n=3),drug_4yr_ago = lag(mean_sp,n=4),drug_5yr_ago = lag(mean_sp,n=5))
write_csv(sp1, file = "SP_lag.csv")
sp1 <- sp1 %>% mutate(year = as.numeric(year))
# Full join main data to drug use data
pfdhps_complete <- full_join(dhps, sp1, by = c("country","year"))
pfdhps_complete <- pfdhps_complete %>% arrange(country, year)
write_csv(pfdhps_complete, file = "pfdhps_prev&pr&sp.2.csv")
#Load data
pfdhps_complete<-read_csv ("pfdhps_prev&pr&sp.2.csv")
pfdhps_complete2 <- pfdhps_complete$country %>% unique()
#Extract out NAs from drug use column and have a good lag year structure. NB.This is the same as right joining drug use data to main data (ie. mdr and cq1). 
pfdhps_complete2 <- pfdhps_complete[!is.na(pfdhps_complete$mean_sp),]
#pfdhps_complete2 <- sum(duplicated(pfdhps_complete2))
pfdhps_complete2 <- pfdhps_complete2[!duplicated(pfdhps_complete2),]
write_csv(pfdhps_complete2, file = "pfdhps_prev&pr&sp.3.csv")



##DHFR
sp<- read_csv("SP_imputed.csv")
dhfr<-read_csv ("pfdhfr_SNP.csv")

pfdhfr_complete <- left_join(dhfr, sp, by = c("country","year"))
pfdhfr_complete <- pfdhfr_complete %>% arrange(country, year)
write_csv(pfdhfr_complete, file = "pfdhfr_prev&pr&sp.1.csv")

#Load data 
DU <- read_csv("pfdhfr_prev&pr&sp.1.csv")
#To interpolate the missing drug use data, cut out "country, year and drug use" columns 
DU <- as.data.frame(cbind(DU$country,DU$year,DU$mean_sp))
colnames(DU)<-c("country","year","SP")
DU <- DU %>% mutate(SP = as.numeric(SP))
DU <- DU %>% mutate(year = as.numeric(year))
#Summarize drug use data from the same 
DU <- DU %>% group_by(country, year) %>% summarise(mean_sp=mean(SP, na.rm = T))
DU <- DU %>% arrange(country, year)
#Expand drug use year from </= 2000 to 2022
c <- DU$country %>% unique()
# loop through and find min and max year, if any years are missing for a country add that row with NAs
for(i in c){
  ro <- which(DU$country == i)
  mi <- min(DU$year[ro])
  mi<- ifelse(mi>2000,2000,mi)
  ma <- max(DU$year[ro])
  ma <- ifelse(ma<2022,2022,ma)
  my_range <- mi:ma
  x <- which(!my_range %in% DU$year[ro])
  if(length(x)>0){
    new_ros <- DU[1:length(x),]
    new_ros$country = i
    new_ros$year = my_range[x]
    new_ros[,3:ncol(new_ros)] = NA
    DU <- rbind(DU, new_ros)
  }
}
DU <- DU %>% arrange(country, year)

# Say if data was observed or interpolated into column "du_obs" ie. NAs (data to be interpolated)=0, Observed data=1
DU$du_obs <- ifelse(is.na(DU$mean_sp), 0, 1)
#Summarize observed data into column "num_obs"
DU<- DU %>% group_by(country) %>% mutate(num_obs=sum(du_obs))
#Filter out country with 3 or more data points for interpolation
c <- DU %>% filter(num_obs>2) %>% 
  select(country) %>% unlist() %>% unique()
#Interpolate data
for(i in c){
  print(i)
  ro <- which(DU$country == i)
  du <- DU$mean_sp[ro]
  DU$mean_sp[ro] <- na_kalman(du)
}
DU <- DU %>% filter(num_obs>2)
# correct for above 1 or below 0
DU$mean_sp <- ifelse(DU$mean_sp > 1, 1, DU$mean_sp)
DU$mean_sp <- ifelse(DU$mean_sp < 0, 0, DU$mean_sp)
#Filter out data from 2000 and above
DU <- DU %>% filter(year>1999)

sp1<- DU
dhfr<-read_csv ("pfdhfr_SNP.csv")

#Set up lag years
sp1<-sp1%>%group_by(country)%>%mutate(drug_1yr_ago = lag(mean_sp,n=1),drug_2yr_ago = lag(mean_sp,n=2),
                                      drug_3yr_ago = lag(mean_sp,n=3),drug_4yr_ago = lag(mean_sp,n=4),drug_5yr_ago = lag(mean_sp,n=5))
write_csv(sp1, file = "SP_lag.csv")
sp1 <- sp1 %>% mutate(year = as.numeric(year))
# Full join main data to drug use data
pfdhfr_complete <- full_join(dhfr, sp1, by = c("country","year"))
pfdhfr_complete <- pfdhfr_complete %>% arrange(country, year)
write_csv(pfdhfr_complete, file = "pfdhfr_prev&pr&sp.2.csv")
#Load data
pfdhfr_complete<-read_csv ("pfdhfr_prev&pr&sp.2.csv")
pfdhfr_complete2 <- pfdhfr_complete$country %>% unique()
#Extract out NAs from drug use column and have a good lag year structure. NB.This is the same as right joining drug use data to main data (ie. mdr and cq1). 
pfdhfr_complete2 <- pfdhfr_complete[!is.na(pfdhfr_complete$mean_sp),]
#pfdhfr_complete2 <- sum(duplicated(pfdhfr_complete2))
pfdhfr_complete2 <- pfdhfr_complete2[!duplicated(pfdhfr_complete2),]
write_csv(pfdhfr_complete2, file = "pfdhfr_prev&pr&sp.3.csv")


##K13
act<- read_csv("ACT_imputed.csv")
k13<-read_csv ("pfk13_SNP.csv")

pfk13_complete <- left_join(k13, act, by = c("country","year"))
pfk13_complete <- pfk13_complete %>% arrange(country, year)
write_csv(pfk13_complete, file = "pfk13_prev&pr&act.1.csv")

#Load data 
DU <- read_csv("pfk13_prev&pr&act.1.csv")
#To interpolate the missing drug use data, cut out "country, year and drug use" columns 
DU <- as.data.frame(cbind(DU$country,DU$year,DU$mean_act))
colnames(DU)<-c("country","year","ACT")
DU <- DU %>% mutate(ACT = as.numeric(ACT))
DU <- DU %>% mutate(year = as.numeric(year))
#Summarize drug use data from the same 
DU <- DU %>% group_by(country, year) %>% summarise(mean_act=mean(ACT, na.rm = T))
DU <- DU %>% arrange(country, year)
#Expand drug use year from </= 2000 to 2022
c <- DU$country %>% unique()
# loop through and find min and max year, if any years are missing for a country add that row with NAs
for(i in c){
  ro <- which(DU$country == i)
  mi <- min(DU$year[ro])
  mi<- ifelse(mi>2000,2000,mi)
  ma <- max(DU$year[ro])
  ma <- ifelse(ma<2022,2022,ma)
  my_range <- mi:ma
  x <- which(!my_range %in% DU$year[ro])
  if(length(x)>0){
    new_ros <- DU[1:length(x),]
    new_ros$country = i
    new_ros$year = my_range[x]
    new_ros[,3:ncol(new_ros)] = NA
    DU <- rbind(DU, new_ros)
  }
}
DU <- DU %>% arrange(country, year)

# Say if data was observed or interpolated into column "du_obs" ie. NAs (data to be interpolated)=0, Observed data=1
DU$du_obs <- ifelse(is.na(DU$mean_act), 0, 1)
#Summarize observed data into column "num_obs"
DU<- DU %>% group_by(country) %>% mutate(num_obs=sum(du_obs))
#Filter out country with 3 or more data points for interpolation
c <- DU %>% filter(num_obs>2) %>% 
  select(country) %>% unlist() %>% unique()
#Interpolate data
for(i in c){
  print(i)
  ro <- which(DU$country == i)
  du <- DU$mean_act[ro]
  DU$mean_act[ro] <- na_kalman(du)
}
DU <- DU %>% filter(num_obs>2)
# correct for above 1 or below 0
DU$mean_act <- ifelse(DU$mean_act > 1, 1, DU$mean_act)
DU$mean_act <- ifelse(DU$mean_act < 0, 0, DU$mean_act)
#Filter out data from 2000 and above
DU <- DU %>% filter(year>1999)

act1<- DU
k13<-read_csv ("pfk13_SNP.csv")

#Set up lag years
act1<-act1%>%group_by(country)%>%mutate(drug_1yr_ago = lag(mean_act,n=1),drug_2yr_ago = lag(mean_act,n=2),
                                      drug_3yr_ago = lag(mean_act,n=3),drug_4yr_ago = lag(mean_act,n=4),drug_5yr_ago = lag(mean_act,n=5))
write_csv(act1, file = "ACT_lag.csv")
act1 <- act1 %>% mutate(year = as.numeric(year))
# Full join main data to drug use data
pfk13_complete <- full_join(k13, act1, by = c("country","year"))
pfk13_complete <- pfk13_complete %>% arrange(country, year)
write_csv(pfk13_complete, file = "pfk13_prev&pr&act.2.csv")
#Load data
pfk13_complete<-read_csv ("pfk13_prev&pr&act.2.csv")
pfk13_complete2 <- pfk13_complete$country %>% unique()
#Extract out NAs from drug use column and have a good lag year structure. NB.This is the same as right joining drug use data to main data (ie. mdr and cq1). 
pfk13_complete2 <- pfk13_complete[!is.na(pfk13_complete$mean_act),]
#pfk13_complete2 <- sum(duplicated(pfk13_complete2))
pfk13_complete2 <- pfk13_complete2[!duplicated(pfk13_complete2),]
write_csv(pfk13_complete2, file = "pfk13_prev&pr&act.3.csv")
#####
#####
#####

##Models
###pfcrt-76T
pfcrt<-read_csv ("pfcrt_prev&pr&cq.3.csv")
pfcrt <- pfcrt%>% mutate(PR = PR / 100, Prev_DR = Prev_DR / 100)

#Mixed models
pfcrt <-pfcrt %>% filter(Prev_DR<1,Prev_DR>0, mean_cq<1,mean_cq>0, PR<1,PR>0,drug_1yr_ago<1,drug_1yr_ago>0,drug_2yr_ago<1,drug_2yr_ago>0,drug_3yr_ago<1,drug_3yr_ago>0,drug_4yr_ago<1,drug_4yr_ago>0,drug_5yr_ago<1,drug_5yr_ago>0)
Mod1_pfcrt <- lmer(logit(Prev_DR) ~  logit(PR) + logit(mean_cq) + logit(drug_1yr_ago) + logit(drug_2yr_ago) + logit(drug_3yr_ago) + logit(drug_4yr_ago) + logit(drug_5yr_ago) +(1|`study Id`), weights=tested, data=pfcrt)
summary(Mod1_pfcrt)
plot(Mod1_pfcrt)

Mod2_pfcrt <- lmer(Prev_DR ~ PR + mean_cq + drug_1yr_ago+drug_2yr_ago+drug_3yr_ago+drug_4yr_ago+drug_5yr_ago+ (1|`study Id`), data=pfcrt)
summary(Mod2_pfcrt)
plot(Mod2_pfcrt)

#Beta regtession
pfcrt <-pfcrt %>% filter(Prev_DR<1,Prev_DR>0)
df <- pfcrt
Mod3_pfcrt <- betareg(Prev_DR ~ PR + mean_cq+ drug_1yr_ago+drug_2yr_ago+drug_3yr_ago+drug_4yr_ago+drug_5yr_ago, data = df)
summary(Mod3_pfcrt)
plot(Mod3_pfcrt)

Mod4_pfcrt <- betareg(Prev_DR ~ PR + mean_cq +drug_1yr_ago+drug_2yr_ago+drug_3yr_ago+drug_4yr_ago+drug_5yr_ago, data = df, link = "loglog")
summary(Mod4_pfcrt)
plot(Mod4_pfcrt)


###pfmdr1-86Y
pfmdr1<-read_csv ("pfmdr1_prev&pr&cq.3.csv")
pfmdr1<- pfmdr1%>% mutate(PR = PR / 100, Prev_DR = Prev_DR / 100)

#Mixed models
pfmdr1 <-pfmdr1 %>% filter(Prev_DR<1,Prev_DR>0, mean_cq<1,mean_cq>0, PR<1,PR>0,drug_1yr_ago<1,drug_1yr_ago>0,drug_2yr_ago<1,drug_2yr_ago>0,drug_3yr_ago<1,drug_3yr_ago>0,drug_4yr_ago<1,drug_4yr_ago>0,drug_5yr_ago<1,drug_5yr_ago>0)
Mod1_pfmdr1 <- lmer(logit(Prev_DR) ~  logit(PR) + logit(mean_cq) + logit(drug_1yr_ago) + logit(drug_2yr_ago) + logit(drug_3yr_ago) + logit(drug_4yr_ago) + logit(drug_5yr_ago) +(1|`study Id`), weights=tested, data=pfmdr1)
summary(Mod1_pfmdr1)
plot(Mod1_pfmdr1)

Mod2_pfmdr1 <- lmer(Prev_DR ~ PR + mean_cq + drug_1yr_ago+drug_2yr_ago+drug_3yr_ago+drug_4yr_ago+drug_5yr_ago+ (1|`study Id`), data=pfmdr1)
summary(Mod2_pfmdr1)
plot(Mod2_pfmdr1)

#Beta regtession
pfmdr1 <-pfmdr1 %>% filter(Prev_DR<1,Prev_DR>0)
df <- pfmdr1
Mod3_pfmdr1 <- betareg(Prev_DR ~ PR + mean_cq+ drug_1yr_ago+drug_2yr_ago+drug_3yr_ago+drug_4yr_ago+drug_5yr_ago, data = df)
summary(Mod3_pfmdr1)
plot(Mod3_pfmdr1)

Mod4_pfmdr1 <- betareg(Prev_DR ~ PR + mean_cq +drug_1yr_ago+drug_2yr_ago+drug_3yr_ago+drug_4yr_ago+drug_5yr_ago, data = df, link = "loglog")
summary(Mod4_pfmdr1)
plot(Mod4_pfmdr1)

###pfdhps-108N
pfdhps<-read_csv ("pfdhps_prev&pr&sp.3.csv")
pfdhps<- pfdhps%>% mutate(PR = PR / 100, Prev_DR = Prev_DR / 100)

#Mixed models
pfdhps <-pfdhps %>% filter(Prev_DR<1,Prev_DR>0, mean_sp<1,mean_sp>0, PR<1,PR>0,drug_1yr_ago<1,drug_1yr_ago>0,drug_2yr_ago<1,drug_2yr_ago>0,drug_3yr_ago<1,drug_3yr_ago>0,drug_4yr_ago<1,drug_4yr_ago>0,drug_5yr_ago<1,drug_5yr_ago>0)
Mod1_pfdhps <- lmer(logit(Prev_DR) ~  logit(PR) + logit(mean_sp) + logit(drug_1yr_ago) + logit(drug_2yr_ago) + logit(drug_3yr_ago) + logit(drug_4yr_ago) + logit(drug_5yr_ago) +(1|`study id`), weights=tested, data=pfdhps)
summary(Mod1_pfdhps)
plot(Mod1_pfdhps)

Mod2_pfdhps <- lmer(Prev_DR ~ PR + mean_sp + drug_1yr_ago+drug_2yr_ago+drug_3yr_ago+drug_4yr_ago+drug_5yr_ago+ (1|`study id`), data=pfdhps)
summary(Mod2_pfdhps)
plot(Mod2_pfdhps)

#Beta regtession
pfdhps <-pfdhps %>% filter(Prev_DR<1,Prev_DR>0)
df <- pfdhps
Mod3_pfdhps <- betareg(Prev_DR ~ PR + mean_sp+ drug_1yr_ago+drug_2yr_ago+drug_3yr_ago+drug_4yr_ago+drug_5yr_ago, data = df)
summary(Mod3_pfdhps)
plot(Mod3_pfdhps)

Mod4_pfdhps <- betareg(Prev_DR ~ PR + mean_sp +drug_1yr_ago+drug_2yr_ago+drug_3yr_ago+drug_4yr_ago+drug_5yr_ago, data = df, link = "loglog")
summary(Mod4_pfdhps)
plot(Mod4_pfdhps)

###pfdhfr-437G
pfdhfr<-read_csv ("pfdhfr_prev&pr&sp.3.csv")
pfdhfr<- pfdhfr%>% mutate(PR = PR / 100, Prev_DR = Prev_DR / 100)

#Mixed models
pfdhfr <-pfdhfr %>% filter(Prev_DR<1,Prev_DR>0, mean_sp<1,mean_sp>0, PR<1,PR>0,drug_1yr_ago<1,drug_1yr_ago>0,drug_2yr_ago<1,drug_2yr_ago>0,drug_3yr_ago<1,drug_3yr_ago>0,drug_4yr_ago<1,drug_4yr_ago>0,drug_5yr_ago<1,drug_5yr_ago>0)
Mod1_pfdhfr <- lmer(logit(Prev_DR) ~  logit(PR) + logit(mean_sp) + logit(drug_1yr_ago) + logit(drug_2yr_ago) + logit(drug_3yr_ago) + logit(drug_4yr_ago) + logit(drug_5yr_ago) +(1|`study id`), weights=tested, data=pfdhfr)
summary(Mod1_pfdhfr)
plot(Mod1_pfdhfr)

Mod2_pfdhfr <- lmer(Prev_DR ~ PR + mean_sp + drug_1yr_ago+drug_2yr_ago+drug_3yr_ago+drug_4yr_ago+drug_5yr_ago+ (1|`study id`), data=pfdhfr)
summary(Mod2_pfdhfr)
plot(Mod2_pfdhfr)

#Beta regtession
pfdhfr <-pfdhfr %>% filter(Prev_DR<1,Prev_DR>0)
df <- pfdhfr
Mod3_pfdhfr <- betareg(Prev_DR ~ PR + mean_sp+ drug_1yr_ago+drug_2yr_ago+drug_3yr_ago+drug_4yr_ago+drug_5yr_ago, data = df)
summary(Mod3_pfdhfr)
plot(Mod3_pfdhfr)

Mod4_pfdhfr <- betareg(Prev_DR ~ PR + mean_sp +drug_1yr_ago+drug_2yr_ago+drug_3yr_ago+drug_4yr_ago+drug_5yr_ago, data = df, link = "loglog")
summary(Mod4_pfdhfr)
plot(Mod4_pfdhfr)


###pfk13-580Y
pfk13<-read_csv ("pfk13_prev&pr&act.3.csv")
pfk13<- pfk13%>% mutate(PR = PR / 100, Prev_DR = Prev_DR / 100)

#Mixed models
#pfk13 <-pfk13 %>% filter(Prev_DR<1,Prev_DR>0, mean_act<1,mean_act>0, PR<1,PR>0,drug_1yr_ago<1,drug_1yr_ago>0,drug_2yr_ago<1,drug_2yr_ago>0,drug_3yr_ago<1,drug_3yr_ago>0,drug_4yr_ago<1,drug_4yr_ago>0,drug_5yr_ago<1,drug_5yr_ago>0)
#Mod1_pfk13 <- lmer(logit(Prev_DR) ~  logit(PR) + logit(mean_act) + logit(drug_1yr_ago) + logit(drug_2yr_ago) + logit(drug_3yr_ago) + logit(drug_4yr_ago) + logit(drug_5yr_ago) +(1|sid), weights=tested, data=pfk13)
#summary(Mod1_pfk13)
#plot(Mod1_pfk13)

Mod2_pfk13 <- lmer(Prev_DR ~ PR + mean_act + drug_1yr_ago+drug_2yr_ago+drug_3yr_ago+drug_4yr_ago+drug_5yr_ago+ (1|sid), data=pfk13)
summary(Mod2_pfk13)
plot(Mod2_pfk13)

#Beta regtession
pfk13 <-pfk13 %>% filter(Prev_DR<1,Prev_DR>0)
df <- pfk13
Mod3_pfk13 <- betareg(Prev_DR ~ PR + mean_act + drug_1yr_ago+drug_2yr_ago+drug_3yr_ago+drug_4yr_ago+drug_5yr_ago, data = df)
summary(Mod3_pfk13)
plot(Mod3_pfk13)

Mod4_pfk13 <- betareg(Prev_DR ~ PR + mean_act +drug_1yr_ago+drug_2yr_ago+drug_3yr_ago+drug_4yr_ago+drug_5yr_ago, data = df, link = "loglog")
summary(Mod4_pfk13)
plot(Mod4_pfk13)


#####
####Maps
# work with maps package
#install.packages(maps)
library(maps)
library(tidyverse)

dataset = read_csv('~/../Downloads/pfmdr1_N86Y(PR.DR_2000_2021).csv')
dataset
ggplot() +
  # Draw the world map
  borders("world", colour = "gray85", fill = "gray80") +
  # Plot points using lat and long
  geom_point(data = dataset, aes(x = long, y = lat, color = CQ_use), size = 2, alpha = 0.7) +
  # Customize the appearance
  theme_minimal() +
  labs(title = "Data Points on World Map", x = "Longitude", y = "Latitude") +
  theme(legend.position = "right")+
  scale_color_viridis_c()











