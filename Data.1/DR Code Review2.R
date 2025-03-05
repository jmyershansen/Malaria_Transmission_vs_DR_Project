#WWARN_full <- read_csv("WWARN_ACT_Partner_Drug_Mol_Surveyor_Data(6).csv")
#pfcrt_SNP <- pfcrt_SNP %>% mutate(Prev_DR=(present/tested)*100)
#write_csv(pfcrt_SNP, file = "pfcrt_SNP.csv")

EXTR <- read.csv("Extracted-points-data.csv")
WWARN_full <- read_csv("WWARN_ACT_Partner_Drug_Mol_Surveyor_Data(6).csv")
WWARN_full <- WWARN_full %>% mutate(year=(`study start`+`study end`)/2)
WWARN_full$year <- ceiling (WWARN_full$year)

colnames(WWARN_full)[5] <- "long"

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
####

cq<- read_csv("CQ_imputed.csv")
#cq2 <- cq$country %>% unique()
#cq2 <- length(unique(cq$country)) 

crt<-read_csv ("pfcrt_SNP.csv")
#crt2 <- crt$country %>% unique()
#crt2 <- length(unique(crt$country)) 


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
crt <- crt %>% arrange(country, year)
#write_csv(crt, file = "pfcrt_SNP.csv")

cq<- read_csv("CQ_imputed.csv")
head(cq)
cq<-cq%>%group_by(country)%>%mutate(drug_1yr_ago = lag(mean_cq,n=1),drug_2yr_ago = lag(mean_cq,n=2),
                                    drug_3yr_ago = lag(mean_cq,n=3),drug_4yr_ago = lag(mean_cq,n=4),drug_5yr_ago = lag(mean_cq,n=5))
#write_csv(cq, file = "CQ_lag.csv")


pfcrt_complete <- right_join(crt, cq, by = c("country","year"))
pfcrt_complete <- pfcrt_complete %>% arrange(country, year)
#pfcrt_complete <- pfcrt_complete %>% arrange(year) %>% group_by(country, year)
write_csv(pfcrt_complete, file = "pfcrt_prev&pr&cq.csv")

pfcrt_complete<-read_csv ("pfcrt_prev&pr&cq.csv")
pfcrt_complete2 <- pfcrt_complete$country %>% unique()
pfcrt_complete2 <- sum(!is.na(pfcrt_complete$Prev_DR))
pfcrt_complete2 <- sum(!is.na(pfcrt_complete$PR))
pfcrt_complete2 <- sum(!is.na(pfcrt_complete$mean_cq))
pfcrt_complete2 <- sum(!duplicated(pfcrt_complete))

#dt = pfcrt_complete[!is.na(pfcrt_complete$Prev_DR),]
#dt = pfcrt_complete[!is.na(pfcrt_complete$PR),]
#dt = pfcrt_complete[!is.na(pfcrt_complete$mean_cq),]
dt = pfcrt_complete[!duplicated(pfcrt_complete),]
write_csv(dt, file = "pfcrt_prev&pr&cq.csv")
#####
#####

cq<- read_csv("CQ_imputed.csv")

pfcrt_complete <- left_join(crt, cq, by = c("country","year"))
pfcrt_complete <- pfcrt_complete %>% arrange(country, year)
#pfcrt_complete <- pfcrt_complete %>% arrange(year) %>% group_by(country, year)
write_csv(pfcrt_complete, file = "pfcrt_prev&pr&cq.1.csv")
#####
#####
#####




##MDR1
pfmdr1_SNP <- read_csv ("pfmdr1_SNP.csv")
pfmdr1_SNP <- pfmdr1_SNP %>% mutate(Prev_DR=(present/tested)*100)
write_csv(pfmdr1_SNP, file = "pfmdr1_SNP.csv")

cq<- read_csv("CQ_imputed.csv")
#cq2 <- cq$country %>% unique()
#cq2 <- length(unique(cq$country)) 

mdr<-read_csv ("pfmdr1_SNP.csv")
#mdr2 <- mdr$country %>% unique()
#mdr2 <- length(unique(mdr$country)) 


b =mdr$country %>% unique()
# loop through and find min and max year, if any years are missing for a country add that row with NAs
for(i in b){
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
#write_csv(mdr, file = "pfmdr1_SNP.csv")

cq<- read_csv("CQ_imputed.csv")
head(cq)
cq<-cq%>%group_by(country)%>%mutate(drug_1yr_ago = lag(mean_cq,n=1),drug_2yr_ago = lag(mean_cq,n=2),
                                    drug_3yr_ago = lag(mean_cq,n=3),drug_4yr_ago = lag(mean_cq,n=4),drug_5yr_ago = lag(mean_cq,n=5))

pfmdr1_complete <- right_join(mdr, cq, by = c("country","year"))
pfmdr1_complete <- pfcrt_complete %>% arrange(country, year)
write_csv(pfmdr1_complete, file = "pfmdr1_prev&pr&cq.csv")



pfmdr1_complete<-read_csv ("pfmdr1_prev&pr&cq.csv")
pfmdr1_complete2 <- pfmdr1_complete$country %>% unique()
pfmdr1_complete2 <- sum(!is.na(pfmdr1_complete$Prev_DR))
pfmdr1_complete2 <- sum(!is.na(pfmdr1_complete$PR))
pfmdr1_complete2 <- sum(!is.na(pfmdr1_complete$mean_cq))
pfmdr1_complete2 <- sum(!duplicated(pfmdr1_complete))
#dt = pfmdr1_complete[!is.na(pfmdr1_complete$Prev_DR),]
#dt = pfmdr1_complete[!is.na(pfmdr1_complete$PR),]
#dt = pfmdr1_complete[!is.na(pfmdr1_complete$mean_cq),]
dt = pfmdr1_complete[!duplicated(pfmdr1_complete),]
write_csv(dt, file = "pfmdr1_prev&pr&cq.csv")
#####
#####

cq<- read_csv("CQ_imputed.csv")

pfmdr1_complete <- left_join(mdr, cq, by = c("country","year"))
pfmdr1_complete <- pfmdr1_complete %>% arrange(country, year)
#pfmdr1_complete <- pfmdr1_complete %>% arrange(year) %>% group_by(country, year)
write_csv(pfmdr1_complete, file = "pfmdr1_prev&pr&cq.1.csv")
######
#####
#####

#DHPS
EXTR1 <- read.csv("Extracted-points-data.1.csv")
DHPS.DHFR <- read_csv("WWARN_pfdhfr_pfdhps_data.csv")
DHPS.DHFR <- DHPS.DHFR %>% mutate(year=(`study start year`+`study end year`)/2)
DHPS.DHFR$year <- ceiling (DHPS.DHFR$year)

colnames(DHPS.DHFR)[6] <- "long"

DHPS.DHFR <- DHPS.DHFR %>% arrange(year) %>% group_by(lat,long,year)

DHPS.DHFR$long <- round(DHPS.DHFR$long, digits = 2)
EXTR1$long <- round (EXTR1$long, digits =2)

DHPS.DHFR$lat <- round(DHPS.DHFR$lat, digits = 2)
EXTR1$lat <- round(EXTR1$lat, digits = 2)

Complete2 <- right_join(EXTR1,DHPS.DHFR,by = c("lat","long","year"))

Complete2 <- Complete2 %>% mutate(PR = as.numeric(value)*100)
View(Complete2)
#Numstudies7 <- length(unique(Complete3$`study id`))

Complete_pfdhps <- Complete2 %>% filter(`mutation`=="dhps 437G")
Complete_pfdhfr <- Complete2 %>% filter(`mutation`=="dhfr 108N")
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


Complete_pfdhfr_region <- right_join(Country_region, Complete_pfdhfr, by = "country")
#write_csv(Complete_pfdhfr_region, file = "pfdhfr_pm_region.csv")
write_csv(Complete_pfdhfr_region, file = "pfdhfr_SNP.csv")

pfdhfr_SNP <- read_csv ("pfdhfr_SNP.csv")
pfdhfr_SNP <- pfdhfr_SNP %>% mutate(Prev_DR=(present/tested)*100)
write_csv(pfdhfr_SNP, file = "pfdhfr_SNP.csv")
Numstudies0 <- length(unique(pfdhfr_SNP$`study id`))
Numstudies0 <- length(unique(pfdhfr_SNP$`PR`))
Numstudies0 <- length(unique(pfdhfr_SNP$`Prev_DR`))
######
#####
#####

sp<- read_csv("SP_imputed.csv")
#sq2 <- sp$country %>% unique()
#sq2 <- length(unique(sp$country)) 

dhps<-read_csv ("pfdhps_SNP.csv")
#dhps2 <- dhps$country %>% unique()
#dhps2 <- length(unique(dhps$country)) 


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
#write_csv(dhps, file = "pfdhps_SNP.csv")

sp<- read_csv("SP_imputed.csv")
head(sp)
sp<-sp%>%group_by(country)%>%mutate(drug_1yr_ago = lag(mean_sp,n=1),drug_2yr_ago = lag(mean_sp,n=2),
                                    drug_3yr_ago = lag(mean_sp,n=3),drug_4yr_ago = lag(mean_sp,n=4),drug_5yr_ago = lag(mean_sp,n=5))

pfdhps_complete <- right_join(dhps, sp, by = c("country","year"))
pfdhps_complete <- pfdhps_complete %>% arrange(country, year)
write_csv(pfdhps_complete, file = "pfdhps_prev&pr&sp.csv")



pfdhps_complete<-read_csv ("pfdhps_prev&pr&sp.csv")
pfdhps_complete2 <- pfdhps_complete$country %>% unique()
pfdhps_complete2 <- sum(!is.na(pfdhps_complete$Prev_DR))
pfdhps_complete2 <- sum(!is.na(pfdhps_complete$PR))
pfdhps_complete2 <- sum(!is.na(pfdhps_complete$mean_sp))
pfdhps_complete2 <- sum(!duplicated(pfdhps_complete))
#dt = pfdhps_complete[!is.na(pfdhps_complete$Prev_DR),]
#dt = pfdhps_complete[!is.na(pfdhps_complete$PR),]
#dt = pfdhps_complete[!is.na(pfdhps_complete$mean_sp),]
dt = pfdhps_complete[!duplicated(pfdhps_complete),]
dt1 = pfdhps_complete[duplicated(pfdhps_complete),]
write_csv(dt, file = "pfdhps_prev&pr&sp.csv")
#####
#####

sp<- read_csv("SP_imputed.csv")

pfdhps_complete <- left_join(dhps, sp, by = c("country","year"))
pfdhps_complete <- pfdhps_complete %>% arrange(country, year)
#pfdhps_complete <- pfdhps_complete %>% arrange(year) %>% group_by(country, year)
write_csv(pfdhps_complete, file = "pfdhps_prev&pr&sp.1.csv")
#####
#####
#####


##DHFR
sp<- read_csv("SP_imputed.csv")
#sp2 <- sp$country %>% unique()
#sp2 <- length(unique(sp$country)) 

dhfr<-read_csv ("pfdhfr_SNP.csv")
#dhfr2 <- dhfr$country %>% unique()
#dhfr2 <- length(unique(dhfr$country)) 


e = dhfr$country %>% unique()
# loop through and find min and max year, if any years are missing for a country add that row with NAs
for(i in e){
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
#write_csv(dhfr, file = "pfdhfr_SNP.csv")

sp<- read_csv("SP_imputed.csv")
head(sp)
sp<-sp%>%group_by(country)%>%mutate(drug_1yr_ago = lag(mean_sp,n=1),drug_2yr_ago = lag(mean_sp,n=2),
                                    drug_3yr_ago = lag(mean_sp,n=3),drug_4yr_ago = lag(mean_sp,n=4),drug_5yr_ago = lag(mean_sp,n=5))

pfdhfr_complete <- right_join(dhfr, sp, by = c("country","year"))
pfdhfr_complete <- pfdhfr_complete %>% arrange(country, year)
write_csv(pfdhfr_complete, file = "pfdhfr_prev&pr&sp.csv")



pfdhfr_complete<-read_csv ("pfdhfr_prev&pr&sp.csv")
pfdhfr_complete2 <- pfdhfr_complete$country %>% unique()
pfdhfr_complete2 <- sum(!is.na(pfdhfr_complete$Prev_DR))
pfdhfr_complete2 <- sum(!is.na(pfdhfr_complete$PR))
pfdhfr_complete2 <- sum(!is.na(pfdhfr_complete$mean_sp))
pfdhfr_complete2 <- sum(!duplicated(pfdhfr_complete))
#dt = pfdhfr_complete[!is.na(pfdhfr_complete$Prev_DR),]
#dt = pfdhfr_complete[!is.na(pfdhfr_complete$PR),]
#dt = pfdhfr_complete[!is.na(pfdhfr_complete$mean_sp),]
dt = pfdhfr_complete[!duplicated(pfdhfr_complete),]
write_csv(dt, file = "pfdhfr_prev&pr&sp.csv")
#####
#####

sp<- read_csv("SP_imputed.csv")

pfdhfr_complete <- left_join(dhfr, sp, by = c("country","year"))
pfdhfr_complete <- pfdhfr_complete %>% arrange(country, year)
#pfdhfr_complete <- pfdhfr_complete %>% arrange(year) %>% group_by(country, year)
write_csv(pfdhfr_complete, file = "pfdhfr_prev&pr&sp.1.csv")
#####
#####
#####

#K13
EXTR2 <- read.csv("Extracted-points-data.2.csv")
K13 <- read_csv("K13_surveyor_data.csv")

colnames(K13)[7] <- "long"

K13 <- K13 %>% arrange(year) %>% group_by(lat,long,year)

K13$long <- round(K13$long, digits = 2)
EXTR2$long <- round (EXTR2$long, digits =2)

K13$lat <- round(K13$lat, digits = 2)
EXTR2$lat <- round(EXTR2$lat, digits = 2)

Complete3 <- right_join(EXTR2,K13,by = c("lat","long","year"))

Complete3 <- Complete3 %>% mutate(PR = as.numeric(value)*100)
View(Complete3)
#Numstudies7 <- length(unique(Complete3$`sid`))

Complete_pfK13 <- Complete3 %>% filter(`mutation`=="C580Y")
#Numstudies8 <- length(unique(Complete_pfK13$`sid`))



Complete_pfK13_region <- right_join(Country_region, Complete_pfK13, by = "country")
write_csv(Complete_pfK13_region, file = "pfk13_SNP.csv")


pfK13_SNP <- read_csv ("pfk13_SNP.csv")
pfK13_SNP <- pfK13_SNP %>% mutate(Prev_DR=(present/tested)*100)
write_csv(pfK13_SNP, file = "pfk13_SNP.csv")
#Numstudies0 <- length(unique(pfK13_SNP$`sid`))
#Numstudies1 <- length(unique(pfK13_SNP$`PR`))
#Numstudies2 <- length(unique(pfK13_SNP$`Prev_DR`))
######
#####
#####

act<- read_csv("ACT_imputed.csv")
#act2 <- act$country %>% unique()
#act2 <- length(unique(act$country))

k13<-read_csv ("pfk13_SNP.csv")
#k13_2 <- k13$country %>% unique()
#k13_2  <- length(unique(k13$country)) 


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
#write_csv(k13, file = "pfk13_SNP.csv")

act<- read_csv("ACT_imputed.csv")
head(act)
act<-act%>%group_by(country)%>%mutate(drug_1yr_ago = lag(mean_act,n=1),drug_2yr_ago = lag(mean_act,n=2),
                                    drug_3yr_ago = lag(mean_act,n=3),drug_4yr_ago = lag(mean_act,n=4),drug_5yr_ago = lag(mean_act,n=5))

pfk13_complete <- right_join(k13, act, by = c("country","year"))
pfk13_complete <- pfk13_complete %>% arrange(country, year)
write_csv(pfk13_complete, file = "pfk13_prev&pr&act.csv")



pfk13_complete<-read_csv ("pfk13_prev&pr&act.csv")
pfk13_complete2 <- pfk13_complete$country %>% unique()
pfk13_complete2 <- sum(!is.na(pfk13_complete$Prev_DR))
pfk13_complete2 <- sum(!is.na(pfk13_complete$PR))
pfk13_complete2 <- sum(!is.na(pfk13_complete$mean_act))
pfk13_complete2 <- sum(!duplicated(pfk13_complete))
#dt = pfk13_complete[!is.na(pfk13_complete$Prev_DR),]
#dt = pfk13_complete[!is.na(pfk13_complete$PR),]
#dt = pfk13_complete[!is.na(pfk13_complete$mean_act),]
dt = pfk13_complete[!duplicated(pfk13_complete),]
write_csv(dt, file = "pfk13_prev&pr&act.csv")
#####
#####

act<- read_csv("ACT_imputed.csv")

pfk13_complete <- left_join(k13, act, by = c("country","year"))
pfk13_complete <- pfk13_complete %>% arrange(country, year)
#pfk13_complete <- pfk13_complete %>% arrange(year) %>% group_by(country, year)
write_csv(pfk13_complete, file = "pfk13_prev&pr&act.1.csv")
######
######
######
######

###CRT
DU <- read_csv("pfcrt_prev&pr&cq.1.csv")
#pfcrt_complete <- pfcrt_complete %>% distinct(country,year,mean_cq)
DU <- as.data.frame(cbind(DU$country,DU$year,DU$mean_cq))
colnames(DU)<-c("country","year","CQ")
                  
# PFCRT Database prep
#DU<-read_csv ("CQ_use.csv")
#DU <- DU %>% group_by(country, year) %>% summarise(mean_cq=mean(CQ, na.rm = T))
#NCountry <- length(unique(DU$country))
#DU <- DU %>% (DU$mean_cq, na.rm = T)
DU <- DU %>% mutate(CQ = as.numeric(CQ))
DU <- DU %>% mutate(year = as.numeric(year))
DU <- DU %>% group_by(country, year) %>% summarise(mean_cq=mean(CQ, na.rm = T))
#DU <- DU %>% group_by(country, year) %>% mutate(mean_cq=mean(mean_cq, na.rm = T))
#DU <- DU %>% group_by(country, year) %>% unique(country, year)
DU = DU %>% arrange(country, year)
#DU = DU[!duplicated(DU),]

# 
#DU %>% group_by(country) %>% 
  #summarise(count = c(), cq_data = sum(!is.na(mean_cq))) %>% view()

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
cq1<- DU
crt<-read_csv ("pfcrt_SNP.csv")

cq1<-cq1%>%group_by(country)%>%mutate(drug_1yr_ago = lag(mean_cq,n=1),drug_2yr_ago = lag(mean_cq,n=2),
                                    drug_3yr_ago = lag(mean_cq,n=3),drug_4yr_ago = lag(mean_cq,n=4),drug_5yr_ago = lag(mean_cq,n=5))

write_csv(cq1, file = "CQ_lag.csv")

#cq1$year = as.numeric(cq1$year)
cq1 <- cq1 %>% mutate(year = as.numeric(year))

pfcrt_complete <- full_join(crt, cq1, by = c("country","year"))
pfcrt_complete <- pfcrt_complete %>% arrange(country, year)
#pfcrt_complete <- pfcrt_complete %>% arrange(year) %>% group_by(country, year)
write_csv(pfcrt_complete, file = "pfcrt_prev&pr&cq.2.csv")

pfcrt_complete<-read_csv ("pfcrt_prev&pr&cq.2.csv")
pfcrt_complete2 <- pfcrt_complete$country %>% unique()
#pfcrt_complete2 <- sum(!is.na(pfcrt_complete$Prev_DR))
#pfcrt_complete2 <- sum(!is.na(pfcrt_complete$PR))
#pfcrt_complete2 <- sum(!is.na(pfcrt_complete$mean_cq))
#pfcrt_complete2 <- sum(!duplicated(pfcrt_complete))
#pfcrt_complete2 <- pfcrt_complete[!is.na(pfcrt_complete$mean_cq),]
#pfcrt_complete2 = pfcrt_complete[!duplicated(pfcrt_complete),]

pfcrt_complete2 <- pfcrt_complete[!is.na(pfcrt_complete$mean_cq),]
#pfmdr1_complete3 <- sum(!is.na(pfmdr1_complete2$Prev_DR))
#pfmdr1_complete2 <- pfmdr1_complete[!duplicated(pfmdr1_complete),]
write_csv(pfcrt_complete2, file = "pfcrt_prev&pr&cq.3.csv")



##MDR1
DU <- read_csv("pfmdr1_prev&pr&cq.1.csv")
DU <- as.data.frame(cbind(DU$country,DU$year,DU$mean_cq))
colnames(DU)<-c("country","year","CQ")

DU <- DU %>% mutate(CQ = as.numeric(CQ))
DU <- DU %>% mutate(year = as.numeric(year))
DU <- DU %>% group_by(country, year) %>% summarise(mean_cq=mean(CQ, na.rm = T))
DU = DU %>% arrange(country, year)

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
cq1<- DU

mdr<-read_csv ("pfmdr1_SNP.csv")

cq1<-cq1%>%group_by(country)%>%mutate(drug_1yr_ago = lag(mean_cq,n=1),drug_2yr_ago = lag(mean_cq,n=2),
                                      drug_3yr_ago = lag(mean_cq,n=3),drug_4yr_ago = lag(mean_cq,n=4),drug_5yr_ago = lag(mean_cq,n=5))
write_csv(cq1, file = "CQ_lag.csv")
cq1 <- cq1 %>% mutate(year = as.numeric(year))

pfmdr1_complete <- full_join(mdr, cq1, by = c("country","year"))
pfmdr1_complete <- pfmdr1_complete %>% arrange(country, year)
write_csv(pfmdr1_complete, file = "pfmdr1_prev&pr&cq.2.csv")

pfmdr1_complete<-read_csv ("pfmdr1_prev&pr&cq.2.csv")
pfmdr1_complete2 <- pfmdr1_complete$country %>% unique()
pfmdr1_complete2 <- sum(!is.na(pfmdr1_complete$Prev_DR))
pfmdr1_complete2 <- sum(!is.na(pfmdr1_complete$PR))
pfmdr1_complete2 <- sum(!is.na(pfmdr1_complete$mean_cq))
pfmdr1_complete2 <- sum(!duplicated(pfmdr1_complete))

pfmdr1_complete2 <- pfmdr1_complete[!is.na(pfmdr1_complete$mean_cq),]
pfmdr1_complete3 <- sum(!is.na(pfmdr1_complete2$Prev_DR))
#pfmdr1_complete2 = pfmdr1_complete[!duplicated(pfmdr1_complete),]
write_csv(pfmdr1_complete2, file = "pfmdr1_prev&pr&cq.3.csv")

##MDR1
DU <- read_csv("pfmdr1_prev&pr&cq.1.csv")
DU <- as.data.frame(cbind(DU$country,DU$year,DU$mean_cq))
colnames(DU)<-c("country","year","CQ")

DU <- DU %>% mutate(CQ = as.numeric(CQ))
DU <- DU %>% mutate(year = as.numeric(year))
DU <- DU %>% group_by(country, year) %>% summarise(mean_cq=mean(CQ, na.rm = T))
DU = DU %>% arrange(country, year)

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
cq1<- DU

mdr<-read_csv ("pfmdr1_SNP.csv")

cq1<-cq1%>%group_by(country)%>%mutate(drug_1yr_ago = lag(mean_cq,n=1),drug_2yr_ago = lag(mean_cq,n=2),
                                      drug_3yr_ago = lag(mean_cq,n=3),drug_4yr_ago = lag(mean_cq,n=4),drug_5yr_ago = lag(mean_cq,n=5))
write_csv(cq1, file = "CQ_lag.csv")
cq1 <- cq1 %>% mutate(year = as.numeric(year))

pfmdr1_complete <- full_join(mdr, cq1, by = c("country","year"))
pfmdr1_complete <- pfmdr1_complete %>% arrange(country, year)
write_csv(pfmdr1_complete, file = "pfmdr1_prev&pr&cq.2.csv")

pfmdr1_complete<-read_csv ("pfmdr1_prev&pr&cq.2.csv")
pfmdr1_complete2 <- pfmdr1_complete$country %>% unique()
pfmdr1_complete2 <- sum(!is.na(pfmdr1_complete$Prev_DR))
pfmdr1_complete2 <- sum(!is.na(pfmdr1_complete$PR))
pfmdr1_complete2 <- sum(!is.na(pfmdr1_complete$mean_cq))
pfmdr1_complete2 <- sum(!duplicated(pfmdr1_complete))

pfmdr1_complete2 <- pfmdr1_complete[!is.na(pfmdr1_complete$mean_cq),]
pfmdr1_complete3 <- sum(!is.na(pfmdr1_complete2$Prev_DR))
#pfmdr1_complete2 = pfmdr1_complete[!duplicated(pfmdr1_complete),]
write_csv(pfmdr1_complete2, file = "pfmdr1_prev&pr&cq.3.csv")
