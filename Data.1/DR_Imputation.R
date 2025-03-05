install.packages("imputeTS")
library(imputeTS)
library(tidyverse)
data(tsAirgap)
tsAirgap
imp <- na_kalman(tsAirgap)
ggplot_na_imputations(tsAirgap, imp, tsAirgapComplete)
ggplot_na_distribution(tsAirgap)
#View (tsAirgap)
#Structure(tsAirgap)          
data<-read_csv ("pfcrt_drug.years.csv")
imp1 <- na_kalman(tsAirgap)
imp1 <- na_kalman(DR)
#dataComplete<- imp1
#ggplot_na_imputations(data, imp, dataComplete)
ggplot_na_distribution(tsAirgap)
ggplot_na_distribution(DR)

DR <- as.data.frame(cbind(data$year,data$CQ_use, data$country))
colnames(DR)<-c("year","CQ_use","country")
####
####
data<-read_csv ("CQ_Kenya.csv")
data$imp = na_kalman(data$CQ_use)
#data$obs = !is.na(data$CQ_use)
#ggplot(data, aes(year, imp))+
  #geom_line()+
  #geom_point(aes(color = obs), size = 3)
ggplot_na_distribution(data$CQ_use)
write_csv(data,file="CQ_Kenya_imp.csv")
###
data<-read_csv ("CQ_Nigeria.csv")
data$imp<- na_kalman(data$CQ_use)
View(data)
ggplot_na_distribution(data$CQ_use)
write_csv(data,file="CQ_Nigeria_imp.csv")
###
data<-read_csv ("CQ_Tanzania.csv")
data$imp<- na_kalman(data$CQ_use)
View(data)
ggplot_na_distribution(data$CQ_use)
write_csv(data,file="CQ_Tanzania_imp.csv")
###
data<-read_csv ("CQ_BFaso.csv")
data$imp<- na_kalman(data$CQ_use)
View(data)
ggplot_na_distribution(data$CQ_use)
write_csv(data,file="CQ_BFaso_imp.csv")
###
data<-read_csv ("CQ_Ghana.csv")
data$imp<- na_kalman(data$CQ_use)
View(data)
ggplot_na_distribution(data$CQ_use)
write_csv(data,file="CQ_Ghana_imp.csv")
###
data<-read_csv ("CQ_Indonesia.csv")
data$imp<- na_kalman(data$CQ_use)
View(data)
ggplot_na_distribution(data$CQ_use)
write_csv(data,file="CQ_Indonesia_imp.csv")
###
data<-read_csv ("CQ_Malawi.csv")
data$imp<- na_kalman(data$CQ_use)
View(data)
ggplot_na_distribution(data$CQ_use)
write_csv(data,file="CQ_Malawi_imp.csv")
###
data<-read_csv ("CQ_Mali.csv")
data$imp<- na_kalman(data$CQ_use)
View(data)
ggplot_na_distribution(data$CQ_use)
write_csv(data,file="CQ_Mali_imp.csv")
###
data<-read_csv ("CQ_Mozambique.csv")
data$imp<- na_kalman(data$CQ_use)
View(data)
ggplot_na_distribution(data$CQ_use)
write_csv(data,file="CQ_Mozambique_imp.csv")
###
data<-read_csv ("CQ_Senegal.csv")
data$imp<- na_kalman(data$CQ_use)
View(data)
ggplot_na_distribution(data$CQ_use)
write_csv(data,file="CQ_Senegal_imp.csv")
###
data<-read_csv ("CQ_Uganda.csv")
data$imp<- na_kalman(data$CQ_use)
View(data)
ggplot_na_distribution(data$CQ_use)
write_csv(data,file="CQ_Uganda_imp.csv")
#####
#####
#####

# PFCRT Database prep
DU<-read_csv ("CQ_use.csv")
DU <- DU %>% group_by(country, year) %>% summarise(mean_cq=mean(CQ, na.rm = T))
#NCountry <- length(unique(DU$country))
c = DU$country %>% unique()
# loop through and find min and max year, if any years are missing for a country add that row with NAs
for(i in c){
  ro = which(DU$country == i)
  mi = min(DU$year[ro])
  mi= ifelse(mi>2000,2000,mi)
  ma = max(DU$year[ro])
  ma = ifelse(ma<2022,2022,ma)
  my_range = mi:ma
  x = which(!my_range %in% DU$year[ro])
  if(length(x)>0){
    new_ros = DU[1:length(x),]
    new_ros$country = i
    new_ros$year = my_range[x]
    new_ros[,3:ncol(new_ros)] = NA
    DU = rbind(DU, new_ros)
  }
}
DU = DU %>% arrange(country, year)

# say if data was observed or interpolated for each point ----
DU$du_obs = ifelse(is.na(DU$mean_cq), 0, 1)

DU<- DU %>% group_by(country) %>% mutate(num_obs=sum(du_obs))
c = DU %>% filter(num_obs>2) %>% 
  select(country) %>% unlist() %>% unique()

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
DU <- DU %>% filter(year>1999)
#DU <- DU %>% filter (year%in%(2000:2021))

write_csv(DU, file = "CQ_imputed.csv")
cq<- read_csv("CQ_imputed.csv")
cq2 <- cq$country %>% unique()
cq2 <- length(unique(cq$country))  

crt<-read_csv ("pfcrt_SNP.csv")
crt2 <- crt$country %>% unique()
crt2 <- length(unique(crt$country)) 
#crt <- crt %>% group_by(country, year) %>% mutate(mean_prev=mean(Prev_DR, na.rm = T))
#crt <- crt %>% mutate(PR = PR / 100, Prev_DR = Prev_DR / 100)

crt <- crt %>% group_by(country, year) %>% summarise(
  mean_Prev=mean(Prev_DR, na.rm = TRUE),
  mean_PR = mean(PR, na.rm = TRUE))

b = crt$country %>% unique()
# loop through and find min and max year, if any years are missing for a country add that row with NAs
for(i in b){
  ro = which(crt$country == i)
  mi = min(crt$year[ro])
  mi= ifelse(mi>2000,2000,mi)
  ma = max(crt$year[ro])
  ma = ifelse(ma<2022,2022,ma)
  my_range = mi:ma
  x = which(!my_range %in% crt$year[ro])
  if(length(x)>0){
    new_ros = crt[1:length(x),]
    new_ros$country = i
    new_ros$year = my_range[x]
    new_ros[,3:ncol(new_ros)] = NA
    crt = rbind(crt, new_ros)
  }
}
crt = crt %>% arrange(country, year)

  
pfcrt_complete <- right_join(crt, cq, by = c("country","year"))
write_csv(pfcrt_complete, file = "pfcrt_prev&pr&cq.csv")

pfcrt_complete<-read_csv ("pfcrt_prev&pr&cq.csv")
pfcrt_complete2 <- pfcrt_complete$country %>% unique()
pfcrt_complete2 <- sum(!is.na(pfcrt_complete$mean_Prev))
dt = pfcrt_complete[!is.na(pfcrt_complete$mean_Prev),]
dt = pfcrt_complete[!is.na(pfcrt_complete$mean_PR),]
dt = pfcrt_complete[!is.na(pfcrt_complete$mean_cq),]


#pfcrt_complete1 <- full_join(crt, cq, by = c("country","year"))
#write_csv(pfcrt_complete1, file = "pfcrt_prev&pr&cq1.csv")
#pfcrt_complete1 <- pfcrt_complete1$country %>% unique()
#####
#####
#####

# PFMDR1 Database prep
cq<- read_csv("CQ_imputed.csv")
cq2 <- cq$country %>% unique()
cq2 <- length(unique(cq$country))  

mdr<-read_csv ("pfmdr1_SNP.csv")
mdr2 <- mdr$country %>% unique()
mdr2 <- length(unique(mdr$country)) 

mdr <- mdr %>% group_by(country, year) %>% summarise(
  mean_Prev=mean(Prev_DR, na.rm = TRUE),
  mean_PR = mean(PR, na.rm = TRUE))

d = mdr$country %>% unique()
# loop through and find min and max year, if any years are missing for a country add that row with NAs
for(i in d){
  ro = which(mdr$country == i)
  mi = min(mdr$year[ro])
  mi= ifelse(mi>2000,2000,mi)
  ma = max(mdr$year[ro])
  ma = ifelse(ma<2022,2022,ma)
  my_range = mi:ma
  x = which(!my_range %in% mdr$year[ro])
  if(length(x)>0){
    new_ros = mdr[1:length(x),]
    new_ros$country = i
    new_ros$year = my_range[x]
    new_ros[,3:ncol(new_ros)] = NA
    mdr = rbind(mdr, new_ros)
  }
}
mdr = mdr %>% arrange(country, year)


pfmdr1_complete <- right_join(mdr, cq, by = c("country","year"))
write_csv(pfmdr1_complete, file = "pfmdr1_prev&pr&cq.csv")

pfmdr1_complete<-read_csv ("pfmdr1_prev&pr&cq.csv")
pfmdr1_complete2 <- pfmdr1_complete$country %>% unique()
pfmdr1_complete2 <- sum(!is.na(pfmdr1_complete$mean_Prev))
dt = pfmdr1_complete[!is.na(pfmdr1_complete$mean_Prev),]
dt = pfmdr1_complete[!is.na(pfmdr1_complete$mean_PR),]
dt = pfmdr1_complete[!is.na(pfmdr1_complete$mean_cq),]

#pfmdr1_complete1 <- full_join(mdr, cq, by = c("country","year"))
#write_csv(pfmdr1_complete1, file = "pfmdr1_prev&pr&cq.1.csv")
#pfmdr1_complete1 <- pfmdr1_complete1$country %>% unique()
#####
#####
#####

# PFDHPS Database prep
DU<-read_csv ("SP_use.csv")
DU <- DU %>% group_by(country, year) %>% summarise(mean_sp=mean(SP, na.rm = T))
#NCountry <- length(unique(DU$country))
c = DU$country %>% unique()
# loop through and find min and max year, if any years are missing for a country add that row with NAs
for(i in c){
  ro = which(DU$country == i)
  mi = min(DU$year[ro])
  mi= ifelse(mi>2000,2000,mi)
  ma = max(DU$year[ro])
  ma = ifelse(ma<2022,2022,ma)
  my_range = mi:ma
  x = which(!my_range %in% DU$year[ro])
  if(length(x)>0){
    new_ros = DU[1:length(x),]
    new_ros$country = i
    new_ros$year = my_range[x]
    new_ros[,3:ncol(new_ros)] = NA
    DU = rbind(DU, new_ros)
  }
}
DU = DU %>% arrange(country, year)

# say if data was observed or interpolated for each point ----
DU$du_obs = ifelse(is.na(DU$mean_sp), 0, 1)

DU<- DU %>% group_by(country) %>% mutate(num_obs=sum(du_obs))
c = DU %>% filter(num_obs>2) %>% 
  select(country) %>% unlist() %>% unique()

for(i in c){
  print(i)
  ro = which(DU$country == i)
  du = DU$mean_sp[ro]
  DU$mean_sp[ro] = na_kalman(du)
}

DU <- DU %>% filter(num_obs>2)
# correct for above 1 or below 0
DU$mean_sp = ifelse(DU$mean_sp > 1, 1, DU$mean_sp)
DU$mean_sp = ifelse(DU$mean_sp < 0, 0, DU$mean_sp)
DU <- DU %>% filter(year>1999)
#DU <- DU %>% filter (year%in%(2000:2021))

write_csv(DU, file = "SP_imputed.csv")
sp<- read_csv("SP_imputed.csv")
sp2 <- sp$country %>% unique()
sp2 <- length(unique(sp$country))  

dhps<-read_csv ("pfdhps_SNP.csv")
dhps2 <- dhps$country %>% unique()
dhps2 <- length(unique(dhps$country)) 

dhps <- dhps %>% group_by(country, year) %>% summarise(
  mean_Prev=mean(Prev_DR, na.rm = TRUE),
  mean_PR = mean(PR, na.rm = TRUE))

e = dhps$country %>% unique()
# loop through and find min and max year, if any years are missing for a country add that row with NAs
for(i in e){
  ro = which(dhps$country == i)
  mi = min(dhps$year[ro])
  mi= ifelse(mi>2000,2000,mi)
  ma = max(dhps$year[ro])
  ma = ifelse(ma<2022,2022,ma)
  my_range = mi:ma
  x = which(!my_range %in% dhps$year[ro])
  if(length(x)>0){
    new_ros = dhps[1:length(x),]
    new_ros$country = i
    new_ros$year = my_range[x]
    new_ros[,3:ncol(new_ros)] = NA
    dhps = rbind(dhps, new_ros)
  }
}
dhps = dhps %>% arrange(country, year)


pfdhps_complete <- right_join(dhps, sp, by = c("country","year"))
write_csv(pfdhps_complete, file = "pfdhps_prev&pr&sp.csv")

pfdhps_complete<-read_csv ("pfdhps_prev&pr&sp.csv")
pfdhps_complete2 <- pfdhps_complete$country %>% unique()
pfdhps_complete2 <- sum(!is.na(pfdhps_complete$mean_Prev))
dt = pfdhps_complete[!is.na(pfdhps_complete$mean_Prev),]
dt = pfdhps_complete[!is.na(pfdhps_complete$mean_PR),]
dt = pfdhps_complete[!is.na(pfdhps_complete$mean_sp),]

#pfdhps_complete1 <- full_join(dhps, sp, by = c("country","year"))
#write_csv(pfdhps_complete1, file = "pfdhps_prev&pr&sp.1.csv")
#pfdhps_complete1 <- pfdhps_complete1$country %>% unique()

#####
#####
#####

# PFDHFR Database prep
sp<- read_csv("SP_imputed.csv")
sp2 <- sp$country %>% unique()
sp2 <- length(unique(sp$country))  

dhfr<-read_csv ("pfdhfr_SNP.csv")
dhfr2 <- dhfr$country %>% unique()
dhfr2 <- length(unique(dhfr$country)) 

dhfr <- dhfr %>% group_by(country, year) %>% summarise(
  mean_Prev=mean(Prev_DR, na.rm = TRUE),
  mean_PR = mean(PR, na.rm = TRUE))

f = dhfr$country %>% unique()
# loop through and find min and max year, if any years are missing for a country add that row with NAs
for(i in f){
  ro = which(dhfr$country == i)
  mi = min(dhfr$year[ro])
  mi= ifelse(mi>2000,2000,mi)
  ma = max(dhfr$year[ro])
  ma = ifelse(ma<2022,2022,ma)
  my_range = mi:ma
  x = which(!my_range %in% dhfr$year[ro])
  if(length(x)>0){
    new_ros = dhfr[1:length(x),]
    new_ros$country = i
    new_ros$year = my_range[x]
    new_ros[,3:ncol(new_ros)] = NA
    dhfr = rbind(dhfr, new_ros)
  }
}
dhfr = dhfr %>% arrange(country, year)


pfdhfr_complete <- right_join(dhfr, sp, by = c("country","year"))
write_csv(pfdhfr_complete, file = "pfdhfr_prev&pr&sp.csv")

pfdhfr_complete<-read_csv ("pfdhfr_prev&pr&sp.csv")
pfdhfr_complete2 <- pfdhfr_complete$country %>% unique()
pfdhfr_complete2 <- sum(!is.na(pfdhfr_complete$mean_Prev))
dt = pfdhfr_complete[!is.na(pfdhfr_complete$mean_Prev),]
dt = pfdhfr_complete[!is.na(pfdhfr_complete$mean_PR),]
dt = pfdhfr_complete[!is.na(pfdhfr_complete$mean_sp),]

#pfdhfr_complete1 <- full_join(dhfr, sp, by = c("country","year"))
#write_csv(pfdhfr_complete1, file = "pfdhfr_prev&pr&sp.1.csv")
#pfdhfr_complete1 <- pfdhfr_complete1$country %>% unique()
#####
#####
#####

# PFK-13 Database prep
DU<-read_csv ("ACT_use.csv")
DU <- DU %>% group_by(country, year) %>% summarise(mean_act=mean(ACT, na.rm = T))
NCountry <- length(unique(DU$country))
c = DU$country %>% unique()
# loop through and find min and max year, if any years are missing for a country add that row with NAs
for(i in c){
  ro = which(DU$country == i)
  mi = min(DU$year[ro])
  mi= ifelse(mi>2000,2000,mi)
  ma = max(DU$year[ro])
  ma = ifelse(ma<2022,2022,ma)
  my_range = mi:ma
  x = which(!my_range %in% DU$year[ro])
  if(length(x)>0){
    new_ros = DU[1:length(x),]
    new_ros$country = i
    new_ros$year = my_range[x]
    new_ros[,3:ncol(new_ros)] = NA
    DU = rbind(DU, new_ros)
  }
}
DU = DU %>% arrange(country, year)

# say if data was observed or interpolated for each point ----
DU$du_obs = ifelse(is.na(DU$mean_act), 0, 1)

DU<- DU %>% group_by(country) %>% mutate(num_obs=sum(du_obs))
DU <- DU %>% group_by(country) %>% 
  mutate(sum_data = sum(mean_act, na.rm = T))

c = DU %>% filter(num_obs>2, sum_data > 0) %>% 
  select(country) %>% unlist() %>% unique()

for(i in c){
  print(i)
  ro = which(DU$country == i)
  du = DU$mean_act[ro]
  DU$mean_act[ro] = na_kalman(du)
}

DU <- DU %>% filter(country %in% c)

#DU <- DU %>% filter(num_obs>2)
# correct for above 1 or below 0
DU$mean_act = ifelse(DU$mean_act > 1, 1, DU$mean_act)
DU$mean_act = ifelse(DU$mean_act < 0, 0, DU$mean_act)
DU <- DU %>% filter(year>1999)
#DU <- DU %>% filter (year%in%(2000:2021))

write_csv(DU, file = "ACT_imputed.csv")
act<- read_csv("ACT_imputed.csv")
act2 <- act$country %>% unique()
act2 <- length(unique(act$country))  

k13<-read_csv ("pfK13_SNP.csv")
k13_2 <- k13$country %>% unique()
k13_2 <- length(unique(k13$country)) 

k13<- k13 %>% group_by(country, year) %>% summarise(
  mean_Prev=mean(Prev_DR, na.rm = TRUE),
  mean_PR = mean(PR, na.rm = TRUE))

e = k13$country %>% unique()
# loop through and find min and max year, if any years are missing for a country add that row with NAs
for(i in e){
  ro = which(k13$country == i)
  mi = min(k13$year[ro])
  mi= ifelse(mi>2000,2000,mi)
  ma = max(k13$year[ro])
  ma = ifelse(ma<2022,2022,ma)
  my_range = mi:ma
  x = which(!my_range %in% k13$year[ro])
  if(length(x)>0){
    new_ros = k13[1:length(x),]
    new_ros$country = i
    new_ros$year = my_range[x]
    new_ros[,3:ncol(new_ros)] = NA
    k13 = rbind(k13, new_ros)
  }
}
k13 = k13 %>% arrange(country, year)


pfk13_complete <- right_join(k13, act, by = c("country","year"))
pfk13_complete <- pfk13_complete %>% select(-sum_data)
write_csv(pfk13_complete, file = "pfk13_prev&pr&act.csv")

pfk13_complete<-read_csv ("pfk13_prev&pr&act.csv")
pfk13_complete2 <- pfk13_complete$country %>% unique()
pfk13_complete2 <- sum(!is.na(pfk13_complete$mean_Prev))
dt = pfk13_complete[!is.na(pfk13_complete$mean_Prev),]
dt = pfk13_complete[!is.na(pfk13_complete$mean_PR),]
dt = pfk13_complete[!is.na(pfk13_complete$mean_act),]

#pfk13_complete1 <- full_join(k13, act, by = c("country","year"))
#write_csv(pfk13_complete1, file = "pfk13_prev&pr&act.1.csv")
#pfk13_complete1 <- pfk13_complete1$country %>% unique()





####
####
####
WWARN <- read_csv("WWARN_ACT_Partner_Drug_Mol_Surveyor_Data(6).csv")
WWARN <- WWARN$country %>% unique()
WWARN1 <- read_csv("WWARN_pfdhfr_pfdhps_data.csv")
WWARN1 <- WWARN1$country %>% unique()
WWARN1 <- (unique(WWARN1$country))

crt<-read_csv ("pfcrt_SNP_complete.csv")
crt <- crt %>% group_by(country, year) %>% mutate(mean_prev=mean(Prev_DR, na.rm = T))
crt = crt$country %>% unique()
NCountry <- (unique(crt$country))
#NCountry <- length(unique(crt$country))
