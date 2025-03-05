rm(list = ls ())
#PFCRT
cq<- read_csv("CQ_imputed.csv")
#Set up lag years
cq<-cq%>%group_by(country)%>%mutate(drug_1yr_ago = lag(mean_cq,n=1),drug_2yr_ago = lag(mean_cq,n=2),
                                      drug_3yr_ago = lag(mean_cq,n=3),drug_4yr_ago = lag(mean_cq,n=4),drug_5yr_ago = lag(mean_cq,n=5))

# join original observation

crt<-read_csv ("pfcrt_SNP.csv")%>%left_join(cq)
crt_na_rem <- crt %>% filter(!is.na(mean_cq))
crt_na_rem$Region[is.na(crt_na_rem$Region)] <- "West Africa"
crt_na_rem<- crt_na_rem%>% mutate(PR = PR / 100, Prev_DR = Prev_DR / 100)
write_csv(crt_na_rem, file = "pfcrt_database.csv")


#crt_na_rem1 <- crt_na_rem[!is.na(crt_na_rem$mean_cq),]
#crt_na_rem2 <- crt_na_rem[!duplicated(crt_na_rem),]

#crt1<-read_csv ("pfcrt_SNP.csv")
#pfcrt_complete <- left_join(crt1, cq, by = c("country","year"))
#pfcrt_complete <- pfcrt_complete %>% arrange(country, year)
#pfcrt_complete <- pfcrt_complete[!is.na(pfcrt_complete$mean_cq),]


###Models
#Mod1_pfcrt <- lmer(Prev_DR ~ PR +mean_cq+PR *mean_cq +drug_1yr_ago+drug_2yr_ago+drug_3yr_ago+drug_4yr_ago+drug_5yr_ago+(1|Region), weights=log(tested),data=crt_na_rem)
#summary(Mod1_pfcrt)


Meta1_pfcrt <- glmer(Prev_DR ~ PR +mean_cq +PR *mean_cq+drug_1yr_ago+drug_2yr_ago+drug_3yr_ago+drug_4yr_ago+drug_5yr_ago+ (1|Region), weights=log(tested), data=crt_na_rem, family = binomial(link = "logit"))
summary(Meta1_pfcrt)
plot(Meta1_pfcrt)
qqnorm(residuals(Meta1_pfcrt))

crt_na_rem <-crt_na_rem %>% filter(Prev_DR<1,Prev_DR>0)
Mod2_pfcrt <- lmer(logit(Prev_DR) ~ PR +mean_cq+PR *mean_cq +drug_1yr_ago+drug_2yr_ago+drug_3yr_ago+drug_4yr_ago+drug_5yr_ago+(1|Region), weights=log(tested),data=crt_na_rem)
summary(Mod2_pfcrt)
# homogeneity of variance
plot(Mod2_pfcrt)
# normality of the residuals
qqnorm(residuals(Mod2_pfcrt))


Mod2_pfcrt <- lmer(Prev_DR ~ PR +mean_cq+PR *mean_cq +drug_1yr_ago+drug_2yr_ago+drug_3yr_ago+(1|country), weights=log(tested),data=crt_na_rem)
summary(Mod2_pfcrt)
# homogeneity of variance
plot(Mod2_pfcrt)
# normality of the residuals
qqnorm(residuals(Mod2_pfcrt))

# Check for additivity assumption
#ggplot(data = crt_na_rem, aes(x = PR, y = Prev_DR)) + geom_point() + geom_smooth(method = "lm") + scale_color_brewer(palette = "Dark2") + theme_classic()
#ggplot(data = crt_na_rem, aes(x = mean_cq, y = Prev_DR)) + geom_point() + geom_smooth(method = "lm") + scale_color_brewer(palette = "Dark2") + theme_classic()

# Check for collinearity assumption
#ind.cont <- crt_na_rem[c("PR", "mean_cq")]
#corr.test(ind.cont, use = "pairwise", method = "pearson", adjust = "none")

# Transform Prev_DR
crt_na_rem <- mutate(crt_na_rem, log_Prev_DR = log(Prev_DR + 1))
crt_na_rem <- mutate(crt_na_rem, sqrt_Prev_DR = sqrt(Prev_DR))
# Transform PR
crt_na_rem <- mutate(crt_na_rem, log_PR = log(PR+1))
crt_na_rem <- mutate(crt_na_rem, sqrt_Prev_DR = sqrt(PR))
# Transform cq
crt_na_rem <- mutate(crt_na_rem, log_mean_cq = log(mean_cq+1))
crt_na_rem <- mutate(crt_na_rem, sqrt_mean_cq = sqrt(mean_cq))

Mod3_pfcrt <- lmer(log_Prev_DR ~ log_PR+log_mean_cq+log_PR*log_mean_cq +drug_1yr_ago+drug_2yr_ago+drug_3yr_ago+drug_4yr_ago+drug_5yr_ago+(1|Region), weights=log(tested),data=crt_na_rem)
summary(Mod3_pfcrt)
# homogeneity of variance
plot(Mod3_pfcrt)
# normality of the residuals
qqnorm(residuals(Mod3_pfcrt))


ggplot(data = crt_na_rem, aes(x = PR, y = Prev_DR,col=Region)) + geom_point(aes(size=log(tested)) + geom_smooth(method = "lm") + scale_color_brewer(palette = "Dark2") + theme_classic())
ggplot(data = crt_na_rem, aes(x = drug_4yr_ago, y = Prev_DR,col=Region)) + geom_point(aes(size=log(tested)) + geom_smooth(method = "lm") + scale_color_brewer(palette = "Dark2") + theme_classic())
ggplot(crt_na_rem, aes(drug_4yr_ago,Prev_DR,col=Region))+geom_point(aes(size=log(tested)))
ggplot(crt_na_rem, aes(cut(drug_4yr_ago,10),Prev_DR,col=Region))+geom_boxplot()


ggplot(crt_na_rem%>%filter(Region=="East Africa"), 
       aes(drug_4yr_ago,Prev_DR,col=country))+geom_point(aes(size=log(tested)))
ggplot(crt_na_rem%>%filter(Region=="Central Africa"), 
       aes(drug_4yr_ago,Prev_DR,col=PR))+geom_point(aes(size=log(tested)))+
  scale_color_viridis_c()
ggplot(crt_na_rem%>%filter(Region=="West Africa"), 
       aes(drug_4yr_ago,Prev_DR,col=PR))+geom_point(aes(size=log(tested)))+
  scale_color_viridis_c()
ggplot(crt_na_rem%>%filter(Region %in% c("SEA","Asia","South America")), 
       aes(drug_4yr_ago,Prev_DR,col=country))+geom_point(aes(size=log(tested)))+
  scale_color_viridis_d()
####
####
####

crt_na_rem<- crt_na_rem%>% mutate(PR = PR * 100, Prev_DR = Prev_DR * 100)
#### Pfcrt_76T Prevalence vs Mean_cq_All regions
ggplot(data=crt_na_rem %>% filter(Region %in% c("West Africa","Central Africa", "East Africa", "Southern Africa", "SEA", "Asia", "South America")))+ aes(x = mean_cq, y = Prev_DR) + geom_point(aes(color=as.factor(Region), size=log(tested)))+ geom_smooth(method = "lm") + xlab("mean_cq") + ylab("Pfcrt 76T Prevalence") 
ggplot(crt_na_rem, aes(cut(mean_cq,10),Prev_DR,col=Region))+geom_boxplot()
ggsave("pfcrt_76T_p1.png") 

###Pfcrt_76T Prevalence vs PR_All regions
ggplot(data=crt_na_rem %>% filter(Region %in% c("West Africa","Central Africa", "East Africa", "Southern Africa", "SEA", "Asia", "South America")))+ aes(x = PR, y = Prev_DR) + geom_point(aes(color=as.factor(Region), size=tested))+ geom_smooth(method = "lm") + xlab("Parasite Rate") + ylab("Pfcrt 76T Prevalence")
ggplot(crt_na_rem, aes(cut(PR,10),Prev_DR,col=Region))+geom_boxplot()
ggsave("pfcrt_76T_p2.png")

### Plot of the relationship between Region and Year varied between pfcrt 76T Prevalence.
pfcrt_76T<-crt_na_rem%>%mutate(Prev_DRCat = cut(Prev_DR,breaks=10*(0:10),include.lowest=T))
ggplot(data = pfcrt_76T, aes (y = Region, x = year,colour = Prev_DRCat)) + geom_point(aes(size = 3)) + scale_color_viridis_d(name="Pfcrt 76T Prevalence") + theme_classic() + xlab ("Year") + ylab ("Region")

ggplot(data = crt_na_rem %>% group_by(Region,year) %>% mutate(mean=mean(Prev_DR)), aes (y = Region, x = year,colour = mean)) + geom_point(aes(size=log(tested))) + scale_color_viridis_c(name="Pfcrt 76T Prevalence") + theme_classic() + xlab ("Year") + ylab ("Region")
ggsave("pfcrt 76T_Timep1.png")

### Plot of the relationship between Region and Year varied between Parasite rate.
pfcrt_76T_I<-crt_na_rem%>%mutate(PRCat = cut(PR,breaks=10*(0:10),include.lowest=T))
ggplot(data = pfcrt_76T_I, aes (y = Region, x = year,colour = PRCat)) + geom_point(size = 3) + scale_color_viridis_d(name="PR") + theme_classic() + xlab ("Year") + ylab ("Region")

ggplot(data = crt_na_rem %>% group_by(Region,year) %>% mutate(mean=mean(PR)), aes (y = Region, x = year,colour = mean)) + geom_point(aes(size=log(tested))) + scale_color_viridis_c(name="PR") + theme_classic() + xlab ("Year") + ylab ("Region")
ggsave("pfcrt 76T_Timep2.png")


### Plot of the relationship between Region and Year varied between Mean_cq.
pfcrt_76T_II<-crt_na_rem%>%mutate(mean_cqCat = cut(mean_cq,breaks=10*(0:10),include.lowest=T))
ggplot(data = pfcrt_76T_II, aes (y = Region, x = year,colour = mean_cqCat)) + geom_point(size = 3) + scale_color_viridis_d(name="Mean_cq") + theme_classic() + xlab ("Year") + ylab ("Region")

ggplot(data = crt_na_rem %>% group_by(Region,year) %>% mutate(mean=mean(mean_cq)), aes (y = Region, x = year,colour = mean)) + geom_point(aes(size=log(tested))) + scale_color_viridis_c(name="Mean_cq") + theme_classic() + xlab ("Year") + ylab ("Region")
ggsave("pfcrt 76T_Timep3.png")


#PFMDR1
# join original observation

mdr<-read_csv ("pfmdr1_SNP.csv")%>%left_join(cq)
mdr_na_rem <- mdr %>% filter(!is.na(mean_cq))
mdr_na_rem$Region[is.na(mdr_na_rem$Region)] <- "West Africa"
mdr_na_rem<- mdr_na_rem%>% mutate(PR = PR / 100, Prev_DR = Prev_DR / 100)
write_csv(mdr_na_rem, file = "pfmdr1_database.csv")


#mdr_na_rem1 <- mdr_na_rem[!is.na(mdr_na_rem$mean_cq),]
#mdr_na_rem2 <- mdr_na_rem[!duplicated(mdr_na_rem),]

###Models
#Models with drug_3yr_ago seems to be better
#Mod1_pfmdr1 <- lmer(Prev_DR ~ PR +mean_cq+PR *mean_cq+drug_1yr_ago+drug_2yr_ago+drug_3yr_ago+drug_4yr_ago+drug_5yr_ago+ (1|Region), weights=log(tested),data=mdr_na_rem)
#summary(Mod1_pfmdr1)
#plot(Mod1_pfmdr1)
#qqnorm(residuals(Mod1_pfmdr1))


Meta1_pfmdr1 <- glmer(Prev_DR ~ PR +mean_cq +PR *mean_cq+drug_1yr_ago+drug_2yr_ago+drug_3yr_ago+drug_4yr_ago+drug_5yr_ago+ (1|Region), weights=log(tested), data=mdr_na_rem, family = binomial(link = "logit"))
summary(Meta1_pfmdr1)
plot(Meta1_pfmdr1)
qqnorm(residuals(Meta1_pfmdr1))

mdr_na_rem <-mdr_na_rem %>% filter(Prev_DR<1,Prev_DR>0)
Mod2_pfmdr1 <- lmer(logit(Prev_DR) ~ PR +mean_cq+PR *mean_cq+drug_1yr_ago+drug_2yr_ago+drug_3yr_ago+drug_4yr_ago+drug_5yr_ago+(1|Region), weights=log(tested),data=mdr_na_rem)
summary(Mod2_pfmdr1)
plot(Mod2_pfmdr1)
qqnorm(residuals(Mod2_pfmdr1))

Mod2_pfmdr1 <- lmer(Prev_DR ~ PR +mean_cq+PR *mean_cq+drug_1yr_ago+drug_2yr_ago+drug_3yr_ago+(1|country), weights=log(tested),data=mdr_na_rem)
summary(Mod2_pfmdr1)
plot(Mod2_pfmdr1)
qqnorm(residuals(Mod2_pfmdr1))

# Transform Prev_DR
mdr_na_rem <- mutate(mdr_na_rem, log_Prev_DR = log(Prev_DR + 1))
mdr_na_rem <- mutate(mdr_na_rem, sqrt_Prev_DR  = sqrt(Prev_DR ))
# Transform PR
mdr_na_rem <- mutate(mdr_na_rem, log_PR = log(PR+1))
mdr_na_rem <- mutate(mdr_na_rem, sqrt_PR = sqrt(PR))
# Transform cq
mdr_na_rem <- mutate(mdr_na_rem, log_mean_cq = log(mean_cq+1))
mdr_na_rem <- mutate(mdr_na_rem, sqrt_mean_cq = sqrt(mean_cq))

Mod3_pfmdr1 <- lmer(log_Prev_DR ~ log_PR+log_mean_cq+log_PR*log_mean_cq +drug_1yr_ago+drug_2yr_ago+drug_3yr_ago+drug_4yr_ago+drug_5yr_ago+(1|Region), weights=log(tested),data=mdr_na_rem)
summary(Mod3_pfmdr1)
# homogeneity of variance
plot(Mod3_pfmdr1)
# normality of the residuals
qqnorm(residuals(Mod3_pfmdr1))


ggplot(mdr_na_rem, aes(drug_3yr_ago,Prev_DR,col=Region))+geom_point(aes(size=log(tested)))
ggplot(mdr_na_rem, aes(drug_5yr_ago,Prev_DR,col=Region))+geom_point(aes(size=log(tested)))
ggplot(mdr_na_rem, aes(cut(drug_3yr_ago,10),Prev_DR,col=Region))+geom_boxplot()
ggplot(mdr_na_rem, aes(cut(drug_5yr_ago,10),Prev_DR,col=Region))+geom_boxplot()

ggplot(mdr_na_rem%>%filter(Region=="East Africa"), 
       aes(drug_3yr_ago,Prev_DR,col=PR))+geom_point(aes(size=log(tested)))+
  scale_color_viridis_c()
ggplot(mdr_na_rem%>%filter(Region=="Central Africa"), 
       aes(drug_3yr_ago,Prev_DR,col=PR))+geom_point(aes(size=log(tested)))+
  scale_color_viridis_c()
ggplot(mdr_na_rem%>%filter(Region=="West Africa"), 
       aes(drug_3yr_ago,Prev_DR,col=PR))+geom_point(aes(size=log(tested)))+
  scale_color_viridis_c()
ggplot(mdr_na_rem%>%filter(Region %in% c("SEA","Asia","South America")), 
       aes(drug_3yr_ago,Prev_DR,col=country))+geom_point(aes(size=log(tested)))+
  scale_color_viridis_d()
ggplot(mdr_na_rem%>%filter(Region=="South America"), 
       aes(drug_3yr_ago,Prev_DR,col=country))+geom_point(aes(size=log(tested)))+
  scale_color_viridis_d()
####
####
####

mdr_na_rem<- mdr_na_rem%>% mutate(PR = PR * 100, Prev_DR = Prev_DR * 100)
####Pfmdr1_86Y Prevalence vs Mean_cq_All regions
ggplot(data=mdr_na_rem %>% filter(Region %in% c("West Africa","Central Africa", "East Africa", "Southern Africa", "SEA", "Asia", "South America")))+ aes(x = mean_cq, y = Prev_DR) + geom_point(aes(color=as.factor(Region), size=log(tested)))+ geom_smooth(method = "lm") + xlab("mean_cq") + ylab("Pfmdr1 86Y Prevalence") 
ggplot(mdr_na_rem, aes(cut(mean_cq,10),Prev_DR,col=Region))+geom_boxplot()
ggsave("pfmdr1_86Y_p1.png") 

###Pfmdr1_86Y Prevalence vs PR_All regions
ggplot(data=mdr_na_rem %>% filter(Region %in% c("West Africa","Central Africa", "East Africa", "Southern Africa", "SEA", "Asia", "South America")))+ aes(x = PR, y = Prev_DR) + geom_point(aes(color=as.factor(Region), size=tested))+ geom_smooth(method = "lm") + xlab("Parasite Rate") + ylab("Pfmdr1 86Y Prevalence")
ggplot(mdr_na_rem, aes(cut(PR,10),Prev_DR,col=Region))+geom_boxplot()
ggsave("pfmdr1_86Y_p2.png")

### Plot of the relationship between Region and Year varied between pfmdr1 86Y Prevalence.
pfmdr1_86Y<-mdr_na_rem%>%mutate(Prev_DRCat = cut(Prev_DR,breaks=10*(0:10),include.lowest=T))
ggplot(data = pfmdr1_86Y, aes (y = Region, x = year,colour = Prev_DRCat)) + geom_point(aes(size = 3)) + scale_color_viridis_d(name="Pfmdr1 86Y Prevalence") + theme_classic() + xlab ("Year") + ylab ("Region")

ggplot(data = mdr_na_rem %>% group_by(Region,year) %>% mutate(mean=mean(Prev_DR)), aes (y = Region, x = year,colour = mean)) + geom_point(aes(size=log(tested))) + scale_color_viridis_c(name="Pfmdr1 86Y Prevalence") + theme_classic() + xlab ("Year") + ylab ("Region")
ggsave("pfmdr1 86Y_Timep1.png")

### Plot of the relationship between Region and Year varied between Parasite rate.
pfmdr1_86Y_I<-mdr_na_rem %>%mutate(PRCat = cut(PR,breaks=10*(0:10),include.lowest=T))
ggplot(data = pfmdr1_86Y_I, aes (y = Region, x = year,colour = PRCat)) + geom_point(size = 3) + scale_color_viridis_d(name="PR") + theme_classic() + xlab ("Year") + ylab ("Region")

ggplot(data = mdr_na_rem %>% group_by(Region,year) %>% mutate(mean=mean(PR)), aes (y = Region, x = year,colour = mean)) + geom_point(aes(size=log(tested))) + scale_color_viridis_c(name="PR") + theme_classic() + xlab ("Year") + ylab ("Region")
ggsave("pfmdr1 86Y_Timep2.png")


### Plot of the relationship between Region and Year varied between Mean_cq.
pfmdr1_86Y_II<-mdr_na_rem%>%mutate(mean_cqCat = cut(mean_cq,breaks=10*(0:10),include.lowest=T))
ggplot(data = pfmdr1_86Y_II, aes (y = Region, x = year,colour = mean_cqCat)) + geom_point(size = 3) + scale_color_viridis_d(name="Mean_cq") + theme_classic() + xlab ("Year") + ylab ("Region")

ggplot(data = mdr_na_rem %>% group_by(Region,year) %>% mutate(mean=mean(mean_cq)), aes (y = Region, x = year,colour = mean)) + geom_point(aes(size=log(tested))) + scale_color_viridis_c(name="Mean_cq") + theme_classic() + xlab ("Year") + ylab ("Region")
ggsave("pfmdr1 86Y_Timep3.png")


#PFDHPS
sp<- read_csv("SP_imputed.csv")
#Set up lag years
sp<-sp%>%group_by(country)%>%mutate(drug_1yr_ago = lag(mean_sp,n=1),drug_2yr_ago = lag(mean_sp,n=2),
                                    drug_3yr_ago = lag(mean_sp,n=3),drug_4yr_ago = lag(mean_sp,n=4),drug_5yr_ago = lag(mean_sp,n=5))

# join original observation

dhps<-read_csv ("pfdhps_SNP.csv")%>%left_join(sp)
dhps_na_rem <- dhps %>% filter(!is.na(mean_sp))
dhps_na_rem<- dhps_na_rem%>% mutate(PR = PR / 100, Prev_DR = Prev_DR / 100)
write_csv(dhps_na_rem, file = "pfdhps_database.csv")


#dhps_na_rem1 <- dhps_na_rem[!is.na(dhps_na_rem$mean_sp),]
#dhps_na_rem2 <- dhps_na_rem[!duplicated(dhps_na_rem),]


Mod1_pfdhps <- lmer(Prev_DR ~ PR +mean_sp +PR *mean_sp+drug_1yr_ago+drug_2yr_ago+drug_3yr_ago+ (1|country), weights=log(tested),data=dhps_na_rem)
summary(Mod1_pfdhps)
plot(Mod1_pfdhps)
qqnorm(residuals(Mod1_pfdhps))

Meta1_pfdhps <- glmer(Prev_DR ~ PR +mean_sp +PR *mean_sp+drug_1yr_ago+drug_2yr_ago+drug_3yr_ago+drug_4yr_ago+drug_5yr_ago+ (1|`study id`), weights=log(tested), data=dhps_na_rem, family = binomial(link = "logit"))
summary(Meta1_pfdhps)
plot(Meta1_pfdhps)
qqnorm(residuals(Meta1_pfdhps))

dhps_na_rem <-dhps_na_rem %>% filter(Prev_DR<1,Prev_DR>0)
Mod2_pfdhps <- lmer(logit(Prev_DR) ~ PR +mean_sp +PR *mean_sp+drug_1yr_ago+drug_2yr_ago+drug_3yr_ago+drug_4yr_ago+drug_5yr_ago+ (1|Region), weights=log(tested),data=dhps_na_rem)
summary(Mod2_pfdhps)
plot(Mod2_pfdhps)
qqnorm(residuals(Mod2_pfdhps))


# Transform Prev_DR
dhps_na_rem <- mutate(dhps_na_rem, log_Prev_DR = log(Prev_DR + 1))
dhps_na_rem <- mutate(dhps_na_rem, sqrt_Prev_DR  = sqrt(Prev_DR ))
# Transform PR
dhps_na_rem <- mutate(dhps_na_rem, log_PR = log(PR+1))
dhps_na_rem <- mutate(dhps_na_rem, sqrt_PR = sqrt(PR))
# Transform sp
dhps_na_rem <- mutate(dhps_na_rem, log_mean_sp = log(mean_sp+1))
dhps_na_rem <- mutate(dhps_na_rem, sqrt_mean_sp = sqrt(mean_sp))

Mod3_pfdhps <- lmer(log_Prev_DR ~ log_PR+log_mean_sp+log_PR*log_mean_sp +drug_1yr_ago+drug_2yr_ago+drug_3yr_ago+drug_4yr_ago+drug_5yr_ago+(1|Region), weights=log(tested),data=dhps_na_rem)
summary(Mod3_pfdhps)
# homogeneity of variance
plot(Mod3_pfdhps)
# normality of the residuals
qqnorm(residuals(Mod3_pfdhps))


ggplot(dhps_na_rem, aes(drug_1yr_ago,Prev_DR,col=Region))+geom_point(aes(size=log(tested)))
ggplot(dhps_na_rem, aes(drug_2yr_ago,Prev_DR,col=Region))+geom_point(aes(size=log(tested)))
ggplot(dhps_na_rem, aes(cut(drug_1yr_ago,10),Prev_DR,col=Region))+geom_boxplot()
ggplot(dhps_na_rem, aes(cut(drug_2yr_ago,10),Prev_DR,col=Region))+geom_boxplot()


ggplot(dhps_na_rem%>%filter(Region=="East Africa"), 
       aes(drug_1yr_ago,Prev_DR,col=country))+geom_point(aes(size=log(tested)))
ggplot(dhps_na_rem%>%filter(Region=="East Africa"), 
       aes(drug_1yr_ago,Prev_DR,col=PR))+geom_point(aes(size=log(tested)))+
  scale_color_viridis_c()
ggplot(dhps_na_rem%>%filter(Region=="Central Africa"), 
       aes(drug_1yr_ago,Prev_DR,col=PR))+geom_point(aes(size=log(tested)))+
  scale_color_viridis_c()
ggplot(dhps_na_rem%>%filter(Region=="West Africa"), 
       aes(drug_1yr_ago,Prev_DR,col=PR))+geom_point(aes(size=log(tested)))+
  scale_color_viridis_c()
ggplot(dhps_na_rem%>%filter(Region %in% c("SEA","Asia","South America")), 
       aes(drug_1yr_ago,Prev_DR,col=country))+geom_point(aes(size=log(tested)))+
  scale_color_viridis_d()
####
####
####

dhps_na_rem<- dhps_na_rem%>% mutate(PR = PR * 100, Prev_DR = Prev_DR * 100)

####Pfdhps_437G Prevalence vs Mean_sp_All regions
ggplot(data=dhps_na_rem %>% filter(Region %in% c("West Africa","Central Africa", "East Africa", "Southern Africa", "SEA", "Asia","South America")))+ aes(x = mean_sp, y = Prev_DR) + geom_point(aes(color=as.factor(Region), size=log(tested)))+ geom_smooth(method = "lm") + xlab("mean_sp") + ylab("Pfdhps 437G Prevalence") 
ggplot(dhps_na_rem, aes(cut(mean_sp,10),Prev_DR,col=Region))+geom_boxplot()
ggsave("pfdhps_437G _p1.png") 

###Pfdhps_437G  Prevalence vs PR_All regions
ggplot(data=dhps_na_rem %>% filter(Region %in% c("West Africa","Central Africa", "East Africa", "Southern Africa", "SEA", "Asia","South America")))+ aes(x = PR, y = Prev_DR) + geom_point(aes(color=as.factor(Region), size=log(tested)))+ geom_smooth(method = "lm") + xlab("Parasite Rate") + ylab("Pfdhps 437G Prevalence")
ggplot(dhps_na_rem, aes(cut(PR,10),Prev_DR,col=Region))+geom_boxplot()
ggsave("pfdhps_437G _p2.png")

### Plot of the relationship between Region and Year varied between pfdhps 437G Prevalence.
pfdhps_437G<-dhps_na_rem %>%mutate(Prev_DRCat = cut(Prev_DR,breaks=10*(0:10),include.lowest=T))
ggplot(data = pfdhps_437G, aes (y = Region, x = year,colour = Prev_DRCat)) + geom_point(aes(size = 3)) + scale_color_viridis_d(name="Pfdhps 437G Prevalence") + theme_classic() + xlab ("Year") + ylab ("Region")

ggplot(data = dhps_na_rem %>% group_by(Region,year) %>% mutate(mean=mean(Prev_DR)), aes (y = Region, x = year,colour = mean)) + geom_point(aes(size=log(tested))) + scale_color_viridis_c(name="Pfdhps 437G Prevalence") + theme_classic() + xlab ("Year") + ylab ("Region")
ggsave("pfdhps 437G_Timep1.png")

### Plot of the relationship between Region and Year varied between Parasite rate.
pfdhps_437G_I<-dhps_na_rem %>%mutate(PRCat = cut(PR,breaks=10*(0:10),include.lowest=T))
ggplot(data = pfdhps_437G_I, aes (y = Region, x = year,colour = PRCat)) + geom_point(size = 3) + scale_color_viridis_d(name="PR") + theme_classic() + xlab ("Year") + ylab ("Region")

ggplot(data = dhps_na_rem %>% group_by(Region,year) %>% mutate(mean=mean(PR)), aes (y = Region, x = year,colour = mean)) + geom_point(aes(size=log(tested))) + scale_color_viridis_c(name="PR") + theme_classic() + xlab ("Year") + ylab ("Region")
ggsave("pfdhps 437G_Timep2.png")


### Plot of the relationship between Region and Year varied between Mean_sp.
pfdhps_437G_II<-dhps_na_rem %>%mutate(mean_spCat = cut(mean_sp,breaks=10*(0:10),include.lowest=T))
ggplot(data = pfdhps_437G_II, aes (y = Region, x = year,colour = mean_spCat)) + geom_point(size = 3) + scale_color_viridis_d(name="Mean_sp") + theme_classic() + xlab ("Year") + ylab ("Region")

ggplot(data = dhps_na_rem %>% group_by(Region,year) %>% mutate(mean=mean(mean_sp)), aes (y = Region, x = year,colour = mean)) + geom_point(aes(size=log(tested))) + scale_color_viridis_c(name="Mean_sp") + theme_classic() + xlab ("Year") + ylab ("Region")
ggsave("pfdhps 437G_Timep3.png")


#PFDHFR
# join original observation

dhfr<-read_csv ("pfdhfr_SNP.csv")%>%left_join(sp)
dhfr_na_rem <- dhfr %>% filter(!is.na(mean_sp))
dhfr_na_rem <- dhfr_na_rem%>% mutate(PR = PR / 100, Prev_DR = Prev_DR / 100)
write_csv(dhfr_na_rem, file = "pfdhfr_database.csv")


#dhfr_na_rem1 <- dhfr_na_rem[!is.na(dhfr_na_rem$mean_sp),]
#dhfr_na_rem2 <- dhfr_na_rem[!duplicated(dhfr_na_rem),]


Mod1_pfdhfr <- lmer(Prev_DR ~ PR +mean_sp+PR *mean_sp+drug_1yr_ago+drug_2yr_ago+drug_3yr_ago+ (1|country), weights=log(tested),data=dhfr_na_rem)
summary(Mod1_pfdhfr)
plot(Mod1_pfdhfr)
qqnorm(residuals(Mod1_pfdhfr))

Meta1_pfdhfr <- glmer(Prev_DR ~ PR +mean_sp +PR *mean_sp+drug_1yr_ago+drug_2yr_ago+drug_3yr_ago+drug_4yr_ago+drug_5yr_ago+ (1|`study id`), weights=log(tested), data=dhfr_na_rem, family = binomial(link = "logit"))
summary(Meta1_pfdhfr)
plot(Meta1_pfdhfr)
qqnorm(residuals(Meta1_pfdhfr))


dhfr_na_rem <-dhfr_na_rem %>% filter(Prev_DR<1,Prev_DR>0)
Mod2_pfdhfr <- lmer(logit(Prev_DR) ~ PR +mean_sp +PR *mean_sp+drug_1yr_ago+drug_2yr_ago+drug_3yr_ago+drug_4yr_ago+drug_5yr_ago+ (1|Region), weights=log(tested),data=dhfr_na_rem)
summary(Mod2_pfdhfr)
plot(Mod2_pfdhfr)
qqnorm(residuals(Mod2_pfdhfr))


# Transform Prev_DR
dhfr_na_rem <- mutate(dhfr_na_rem, log_Prev_DR = log(Prev_DR + 1))
dhfr_na_rem <- mutate(dhfr_na_rem, sqrt_Prev_DR  = sqrt(Prev_DR ))
# Transform PR
dhfr_na_rem <- mutate(dhfr_na_rem, log_PR = log(PR+1))
dhfr_na_rem <- mutate(dhfr_na_rem, sqrt_PR = sqrt(PR))
# Transform sp
dhfr_na_rem <- mutate(dhfr_na_rem, log_mean_sp = log(mean_sp+1))
dhfr_na_rem <- mutate(dhfr_na_rem, sqrt_mean_sp = sqrt(mean_sp))

Mod3_pfdhfr <- lmer(log_Prev_DR ~ log_PR+log_mean_sp+log_PR*log_mean_sp +drug_1yr_ago+drug_2yr_ago+drug_3yr_ago+drug_4yr_ago+drug_5yr_ago+(1|Region), weights=log(tested),data=dhfr_na_rem)
summary(Mod3_pfdhfr)
# homogeneity of variance
plot(Mod3_pfdhfr)
# normality of the residuals
qqnorm(residuals(Mod3_pfdhfr))


ggplot(dhfr_na_rem, aes(drug_4yr_ago,Prev_DR,col=Region))+geom_point(aes(size=log(tested)))
ggplot(dhfr_na_rem, aes(drug_5yr_ago,Prev_DR,col=Region))+geom_point(aes(size=log(tested)))
ggplot(dhfr_na_rem, aes(cut(drug_4yr_ago,10),Prev_DR,col=Region))+geom_boxplot()
ggplot(dhfr_na_rem, aes(cut(drug_5yr_ago,10),Prev_DR,col=Region))+geom_boxplot()


ggplot(dhfr_na_rem%>%filter(Region=="East Africa"), 
       aes(drug_4yr_ago,Prev_DR,col=country))+geom_point(aes(size=log(tested)))
ggplot(dhfr_na_rem%>%filter(Region=="East Africa"), 
       aes(drug_4yr_ago,Prev_DR,col=PR))+geom_point(aes(size=log(tested)))+
  scale_color_viridis_c()
ggplot(dhfr_na_rem%>%filter(Region=="Central Africa"), 
       aes(drug_1yr_ago,Prev_DR,col=PR))+geom_point(aes(size=log(tested)))+
  scale_color_viridis_c()
ggplot(dhfr_na_rem%>%filter(Region=="West Africa"), 
       aes(drug_1yr_ago,Prev_DR,col=PR))+geom_point(aes(size=log(tested)))+
  scale_color_viridis_c()
ggplot(dhfr_na_rem%>%filter(Region %in% c("SEA","Asia","South America")), 
       aes(drug_1yr_ago,Prev_DR,col=PR))+geom_point(aes(size=log(tested)))+
  scale_color_viridis_c()
ggplot(dhfr_na_rem%>%filter(Region %in% c("SEA","Asia","South America")), 
       aes(drug_1yr_ago,Prev_DR,col=country))+geom_point(aes(size=log(tested)))+
  scale_color_viridis_d()
####
####
####

dhfr_na_rem <- dhfr_na_rem%>% mutate(PR = PR * 100, Prev_DR = Prev_DR * 100)

####Pfdhfr_108N Prevalence vs Mean_sp_All regions
ggplot(data=dhfr_na_rem %>% filter(Region %in% c("West Africa","Central Africa", "East Africa", "Southern Africa", "SEA", "Asia","South America")))+ aes(x = mean_sp, y = Prev_DR) + geom_point(aes(color=as.factor(Region), size=log(tested)))+ geom_smooth(method = "lm") + xlab("mean_sp") + ylab("Pfdhfr 108N Prevalence") 
ggplot(dhfr_na_rem, aes(cut(mean_sp,10),Prev_DR,col=Region))+geom_boxplot()
ggsave("pfdhfr_108N_p1.png") 

###Pfdhfr_108N Prevalence vs PR_All regions
ggplot(data=dhfr_na_rem %>% filter(Region %in% c("West Africa","Central Africa", "East Africa", "Southern Africa", "SEA", "Asia", "South America")))+ aes(x = PR, y = Prev_DR) + geom_point(aes(color=as.factor(Region), size=log(tested)))+ geom_smooth(method = "lm") + xlab("Parasite Rate") + ylab("Pfdhfr 108N Prevalence")
ggplot(dhfr_na_rem, aes(cut(PR,10),Prev_DR,col=Region))+geom_boxplot()
ggsave("pfdhfr_108N_p2.png")

### Plot of the relationship between Region and Year varied between pfdhps 108N Prevalence.
pfdhfr_108N<-dhfr_na_rem %>%mutate(Prev_DRCat = cut(Prev_DR,breaks=10*(0:10),include.lowest=T))
ggplot(data = pfdhfr_108N, aes (y = Region, x = year,colour = Prev_DRCat)) + geom_point(aes(size = 3)) + scale_color_viridis_d(name="Pfdhfr 108N Prevalence") + theme_classic() + xlab ("Year") + ylab ("Region")

ggplot(data = dhfr_na_rem %>% group_by(Region,year) %>% mutate(mean=mean(Prev_DR)), aes (y = Region, x = year,colour = mean)) + geom_point(aes(size=log(tested))) + scale_color_viridis_c(name="Pfdhfr 108N Prevalence") + theme_classic() + xlab ("Year") + ylab ("Region")
ggsave("pfdhfr 108N_Timep1.png")

### Plot of the relationship between Region and Year varied between Parasite rate.
pfdhfr_108N_I<-dhfr_na_rem %>%mutate(PRCat = cut(PR,breaks=10*(0:10),include.lowest=T))
ggplot(data = pfdhfr_108N_I, aes (y = Region, x = year,colour = PRCat)) + geom_point(size = 3) + scale_color_viridis_d(name="PR") + theme_classic() + xlab ("Year") + ylab ("Region")

ggplot(data = dhfr_na_rem %>% group_by(Region,year) %>% mutate(mean=mean(PR)), aes (y = Region, x = year,colour = mean)) + geom_point(aes(size=log(tested))) + scale_color_viridis_c(name="PR") + theme_classic() + xlab ("Year") + ylab ("Region")
ggsave("pfdhfr 108N_Timep2.png")


### Plot of the relationship between Region and Year varied between Mean_sp.
pfdhfr_108N_II<-dhfr_na_rem %>%mutate(mean_spCat = cut(mean_sp,breaks=10*(0:10),include.lowest=T))
ggplot(data = pfdhfr_108N_II, aes (y = Region, x = year,colour = mean_spCat)) + geom_point(size = 3) + scale_color_viridis_d(name="Mean_sp") + theme_classic() + xlab ("Year") + ylab ("Region")

ggplot(data = dhfr_na_rem %>% group_by(Region,year) %>% mutate(mean=mean(mean_sp)), aes (y = Region, x = year,colour = mean)) + geom_point(aes(size=log(tested))) + scale_color_viridis_c(name="Mean_sp") + theme_classic() + xlab ("Year") + ylab ("Region")
ggsave("pfdhfr 108N_Timep3.png")


#PFK13
act<- read_csv("ACT_imputed.csv")

#Set up lag years
act<-act%>%group_by(country)%>%mutate(drug_1yr_ago = lag(mean_act,n=1),drug_2yr_ago = lag(mean_act,n=2),
                                    drug_3yr_ago = lag(mean_act,n=3),drug_4yr_ago = lag(mean_act,n=4),drug_5yr_ago = lag(mean_act,n=5))

# join original observation

k13<-read_csv ("pfk13_SNP.csv")%>%left_join(act)
k13_na_rem <- k13 %>% filter(!is.na(mean_act))
k13_na_rem<- k13_na_rem %>% mutate(PR = PR / 100, Prev_DR = Prev_DR / 100)
write_csv(k13_na_rem, file = "pfk13_database.csv")


#k13_na_rem1 <- k13_na_rem[!is.na(k13_na_rem$mean_act),]
#k13_na_rem2 <- k13_na_rem[!duplicated(k13_na_rem),]


#This Model did not converge
Mod1_pfk13 <- lmer(Prev_DR ~ PR +mean_act+PR *mean_act +drug_1yr_ago+drug_2yr_ago+drug_3yr_ago+ (1|country), weights=log(tested),data=k13_na_rem)
summary(Mod1_pfk13)
plot(Mod1_pfk13)
qqnorm(residuals(Mod1_pfk13))

k13_na_rem <-k13_na_rem %>% filter(Prev_DR<1,Prev_DR>0)
Mod2_pfk13 <- lmer(logit(Prev_DR) ~ PR +mean_act+PR *mean_act +drug_1yr_ago+drug_2yr_ago+drug_3yr_ago+drug_4yr_ago+drug_5yr_ago+ (1|country), weights=log(tested),data=k13_na_rem)
summary(Mod2_pfk13)
plot(Mod2_pfk13)
qqnorm(residuals(Mod2_pfk13))

# Transform Prev_DR
k13_na_rem <- mutate(k13_na_rem, log_Prev_DR = log(Prev_DR + 1))
k13_na_rem <- mutate(k13_na_rem, sqrt_Prev_DR  = sqrt(Prev_DR ))
# Transform PR
k13_na_rem <- mutate(k13_na_rem, log_PR = log(PR+1))
k13_na_rem <- mutate(k13_na_rem, sqrt_PR = sqrt(PR))
# Transform act
k13_na_rem <- mutate(k13_na_rem, log_mean_act = log(mean_act+1))
k13_na_rem <- mutate(k13_na_rem, sqrt_mean_act = sqrt(mean_act))

Mod3_pfk13 <- lmer(log_Prev_DR ~ log_PR+log_mean_act+log_PR*log_mean_act +drug_1yr_ago+drug_2yr_ago+drug_3yr_ago+drug_4yr_ago+drug_5yr_ago+(1|country), weights=log(tested),data=k13_na_rem)
summary(Mod3_pfk13)
# homogeneity of variance
plot(Mod3_pfk13)
# normality of the residuals
qqnorm(residuals(Mod3_pfk13))


ggplot(k13_na_rem, aes(drug_4yr_ago,Prev_DR,col=country))+geom_point(aes(size=log(tested)))
ggplot(k13_na_rem, aes(drug_5yr_ago,Prev_DR,col=country))+geom_point(aes(size=log(tested)))
ggplot(k13_na_rem, aes(cut(drug_4yr_ago,10),Prev_DR,col=country))+geom_boxplot()
ggplot(k13_na_rem, aes(cut(drug_5yr_ago,10),Prev_DR,col=country))+geom_boxplot()

####
####
####

k13_na_rem<- k13_na_rem %>% mutate(PR = PR * 100, Prev_DR = Prev_DR * 100)

####Pfk13_580Y Prevalence vs Mean_act_All regions
ggplot(data=k13_na_rem, aes(x = mean_act, y = Prev_DR)) + geom_point(aes(color=as.factor(country), size=log(tested)))+ geom_smooth(method = "lm") + xlab("mean_act") + ylab("Pfk13 580Y Prevalence") 
ggplot(k13_na_rem, aes(cut(mean_act,10),Prev_DR,col=country))+geom_boxplot()
ggsave("pfK13_580Y_p1.png") 

###Pfdhfr_108N Prevalence vs PR_All regions
ggplot(data=k13_na_rem, aes(x = PR, y = Prev_DR)) + geom_point(aes(color=as.factor(country), size=log(tested)))+ geom_smooth(method = "lm") + xlab("Parasite Rate") + ylab("Pfk13 580Y Prevalence")
ggplot(k13_na_rem, aes(cut(PR,10),Prev_DR,col=country))+geom_boxplot()
ggsave("pfk13_580Y_p2.png")

### Plot of the relationship between Region and Year varied between pfk13 580Y Prevalence.
pfk13_580Y<-k13_na_rem %>%mutate(Prev_DRCat = cut(Prev_DR,breaks=10*(0:10),include.lowest=T))
ggplot(data = pfk13_580Y, aes (y = country, x = year,colour = Prev_DRCat)) + geom_point(aes(size = 3)) + scale_color_viridis_d(name="Pfk13 580Y Prevalence") + theme_classic() + xlab ("Year") + ylab ("Country")

ggplot(data = k13_na_rem %>% group_by(country,year) %>% mutate(mean=mean(Prev_DR)), aes (y = country, x = year,colour = mean)) + geom_point(aes(size=log(tested))) + scale_color_viridis_c(name="Pfk13 580Y Prevalence") + theme_classic() + xlab ("Year") + ylab ("Country")
ggsave("pfk13 580Y_Timep1.png")

### Plot of the relationship between Region and Year varied between Parasite rate.
pfk13_580Y_I<-k13_na_rem %>%mutate(PRCat = cut(PR,breaks=10*(0:10),include.lowest=T))
ggplot(data = pfk13_580Y_I, aes (y = country, x = year,colour = PRCat)) + geom_point(size = 3) + scale_color_viridis_d(name="PR") + theme_classic() + xlab ("Year") + ylab ("Country")

ggplot(data = k13_na_rem %>% group_by(country,year) %>% mutate(mean=mean(PR)), aes (y = country, x = year,colour = mean)) + geom_point(aes(size=log(tested))) + scale_color_viridis_c(name="PR") + theme_classic() + xlab ("Year") + ylab ("Country")
ggsave("pfk13 580Y_Timep2.png")

### Plot of the relationship between Region and Year varied between Mean_act.
pfk13_580Y_II<-k13_na_rem %>%mutate(mean_actCat = cut(mean_act,breaks=10*(0:10),include.lowest=T))
ggplot(data = pfk13_580Y_II, aes (y = country, x = year,colour = mean_actCat)) + geom_point(size = 3) + scale_color_viridis_d(name="Mean_act") + theme_classic() + xlab ("Year") + ylab ("Country")

ggplot(data = k13_na_rem %>% group_by(country,year) %>% mutate(mean=mean(mean_act)), aes (y = country, x = year,colour = mean)) + geom_point(aes(size=log(tested))) + scale_color_viridis_c(name="Mean_act") + theme_classic() + xlab ("Year") + ylab ("Country")
ggsave("pfk13 580Y_Timep3.png")
