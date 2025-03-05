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

#crt1<-read_csv ("pfcrt_SNP.csv")
#pfcrt_complete <- left_join(crt1, cq, by = c("country","year"))
#pfcrt_complete <- pfcrt_complete %>% arrange(country, year)
#pfcrt_complete <- pfcrt_complete[!is.na(pfcrt_complete$mean_cq),]



#Models with mean_cq, drug_1yr_ago and drug_4yr_ago gave stronger significance

#crt_na_rem <-crt_na_rem %>% filter(Prev_DR<1,Prev_DR>0)

Mod2_pfcrt <- lmer(Prev_DR ~ PR +mean_cq+PR *mean_cq +drug_1yr_ago+drug_2yr_ago+drug_3yr_ago+drug_4yr_ago+drug_5yr_ago+(1|Region), weights=log(tested),data=crt_na_rem)
summary(Mod2_pfcrt)

Mod2_pfcrt <- lmer(logit(Prev_DR) ~ PR +mean_cq+PR *mean_cq +drug_1yr_ago+drug_2yr_ago+drug_3yr_ago+drug_4yr_ago+drug_5yr_ago+(1|Region), weights=log(tested),data=crt_na_rem)
summary(Mod2_pfcrt)

#Anova(Mod2_pfcrt, type = 3)
#Anova(Mod2_pfcrt)
#anova(Mod2_pfcrt)
#anova_stats(Mod2_pfcrt)
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

Mod3_pfcrt <- lmer(log_Prev_DR ~ log_PR+mean_cq+log_PR*mean_cq +drug_1yr_ago+drug_2yr_ago+drug_3yr_ago+drug_4yr_ago+drug_5yr_ago+(1|Region), weights=log(tested),data=crt_na_rem)
summary(Mod3_pfcrt)
#Anova(Mod3_pfcrt, type = 3)
Anova(Mod3_pfcrt)
anova(Mod3_pfcrt)
# homogeneity of variance
plot(Mod3_pfcrt)
# normality of the residuals
qqnorm(residuals(Mod3_pfcrt))



ggplot(crt_na_rem, aes(drug_1yr_ago,Prev_DR,col=Region))+geom_point(aes(size=log(tested)))

ggplot(crt_na_rem, aes(drug_4yr_ago,Prev_DR,col=Region))+geom_point(aes(size=log(tested)))

ggplot(crt_na_rem, aes(cut(drug_1yr_ago,10),Prev_DR,col=Region))+geom_boxplot()

ggplot(crt_na_rem, aes(cut(drug_4yr_ago,10),Prev_DR,col=Region))+geom_boxplot()


ggplot(crt_na_rem%>%filter(Region=="East Africa"), 
       aes(drug_1yr_ago,Prev_DR,col=country))+geom_point(aes(size=log(tested)))

ggplot(crt_na_rem%>%filter(Region=="East Africa"), 
       aes(drug_4yr_ago,Prev_DR,col=country))+geom_point(aes(size=log(tested)))


ggplot(crt_na_rem%>%filter(Region=="Central Africa"), 
       aes(drug_2yr_ago,Prev_DR,col=PR))+geom_point(aes(size=log(tested)))+
  scale_color_viridis_c()

ggplot(crt_na_rem%>%filter(Region=="West Africa"), 
       aes(drug_4yr_ago,Prev_DR,col=PR))+geom_point(aes(size=log(tested)))+
  scale_color_viridis_c()

ggplot(crt_na_rem%>%filter(Region %in% c("SEA","Asia","South America")), 
       aes(drug_2yr_ago,Prev_DR,col=country))+geom_point(aes(size=log(tested)))+
  scale_color_viridis_d()

#PFMDR1
# join original observation

mdr<-read_csv ("pfmdr1_SNP.csv")%>%left_join(cq)
mdr_na_rem <- mdr %>% filter(!is.na(mean_cq))
mdr_na_rem$Region[is.na(mdr_na_rem$Region)] <- "West Africa"

#Models with drug_3yr_ago seems to be better
Mod2_pfmdr1 <- lmer(Prev_DR ~ PR +mean_cq+ (1|Region), weights=log(tested),data=mdr_na_rem)
summary(Mod2_pfmdr1)
Anova(Mod2_pfmdr1)
anova(Mod2_pfmdr1)
plot(Mod2_pfmdr1)
qqnorm(residuals(Mod2_pfmdr1))

# Transform Prev_DR
mdr_na_rem <- mutate(mdr_na_rem, log_Prev_DR = log(Prev_DR + 1))
# Transform PR
mdr_na_rem <- mutate(mdr_na_rem, log_PR = log(PR+1))

Mod3_pfmdr1 <- lmer(log_Prev_DR ~ log_PR+mean_cq+log_PR*mean_cq +(1|Region), weights=log(tested),data=mdr_na_rem)
summary(Mod3_pfmdr1)
#Anova(Mod3_pfmdr1, type = 3)
Anova(Mod3_pfmdr1)
anova(Mod3_pfmdr1)
# homogeneity of variance
plot(Mod3_pfmdr1)
# normality of the residuals
qqnorm(residuals(Mod3_pfmdr1))






Mod2_pfmdr1 <- lmer(Prev_DR ~ PR +drug_3yr_ago+ (1|Region), weights=log(tested),data=mdr_na_rem)
summary(Mod2_pfmdr1)
plot(Mod2_pfmdr1)
qqnorm(residuals(Mod2_pfmdr1))

Mod2_pfmdr1 <- lmer(Prev_DR ~ PR +drug_5yr_ago+ (1|Region), weights=log(tested),data=mdr_na_rem)
summary(Mod2_pfmdr1)
plot(Mod2_pfmdr1)
qqnorm(residuals(Mod2_pfmdr1))

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

ggplot(crt_na_rem%>%filter(Region %in% c("SEA","Asia","South America")), 
       aes(drug_3yr_ago,Prev_DR,col=country))+geom_point(aes(size=log(tested)))+
  scale_color_viridis_d()

ggplot(crt_na_rem%>%filter(Region=="South America"), 
       aes(drug_3yr_ago,Prev_DR,col=country))+geom_point(aes(size=log(tested)))+
  scale_color_viridis_d()


#PFDHPS
sp<- read_csv("SP_imputed.csv")
#Set up lag years
sp<-sp%>%group_by(country)%>%mutate(drug_1yr_ago = lag(mean_sp,n=1),drug_2yr_ago = lag(mean_sp,n=2),
                                    drug_3yr_ago = lag(mean_sp,n=3),drug_4yr_ago = lag(mean_sp,n=4),drug_5yr_ago = lag(mean_sp,n=5))

# join original observation

dhps<-read_csv ("pfdhps_SNP.csv")%>%left_join(sp)
dhps_na_rem <- dhps %>% filter(!is.na(mean_sp))

Mod2_pfdhps <- lmer(Prev_DR ~ PR +mean_sp+ (1|Region), weights=log(tested),data=dhps_na_rem)
summary(Mod2_pfdhps)
plot(Mod2_pfdhps)
qqnorm(residuals(Mod2_pfdhps))

Mod2_pfdhps <- lmer(Prev_DR ~ PR +drug_1yr_ago+ (1|Region), weights=log(tested),data=dhps_na_rem)
summary(Mod2_pfdhps)
plot(Mod2_pfdhps)
qqnorm(residuals(Mod2_pfdhps))

Mod2_pfdhps <- lmer(Prev_DR ~ PR +drug_2yr_ago+ (1|Region), weights=log(tested),data=dhps_na_rem)
summary(Mod2_pfdhps)
plot(Mod2_pfdhps)
qqnorm(residuals(Mod2_pfdhps))

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


#PFDHFR
# join original observation

dhfr<-read_csv ("pfdhfr_SNP.csv")%>%left_join(sp)
dhfr_na_rem <- dhfr %>% filter(!is.na(mean_sp))

Mod2_pfdhfr <- lmer(Prev_DR ~ PR +mean_sp+ (1|Region), weights=log(tested),data=dhfr_na_rem)
summary(Mod2_pfdhfr)
plot(Mod2_pfdhfr)
qqnorm(residuals(Mod2_pfdhfr))

Mod2_pfdhfr <- lmer(Prev_DR ~ PR +drug_4yr_ago+ (1|Region), weights=log(tested),data=dhfr_na_rem)
summary(Mod2_pfdhfr)
plot(Mod2_pfdhfr)
qqnorm(residuals(Mod2_pfdhfr))

Mod2_pfdhfr <- lmer(Prev_DR ~ PR +drug_5yr_ago+ (1|Region), weights=log(tested),data=dhfr_na_rem)
summary(Mod2_pfdhfr)
plot(Mod2_pfdhfr)
qqnorm(residuals(Mod2_pfdhfr))

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

#PFK13
act<- read_csv("ACT_imputed.csv")
#Set up lag years
act<-act%>%group_by(country)%>%mutate(drug_1yr_ago = lag(mean_act,n=1),drug_2yr_ago = lag(mean_act,n=2),
                                    drug_3yr_ago = lag(mean_act,n=3),drug_4yr_ago = lag(mean_act,n=4),drug_5yr_ago = lag(mean_act,n=5))

# join original observation

k13<-read_csv ("pfk13_SNP.csv")%>%left_join(act)
k13_na_rem <- k13 %>% filter(!is.na(mean_act))

#This Model did not converge
Mod2_pfk13 <- lmer(Prev_DR ~ PR +mean_act+drug_1yr_ago+ (1|country), weights=log(tested),data=k13_na_rem)
summary(Mod2_pfk13)
plot(Mod2_pfk13)
qqnorm(residuals(Mod2_pfk13))

Mod2_pfk13 <- lmer(Prev_DR ~ PR +drug_4yr_ago+ (1|country), weights=log(tested),data=k13_na_rem)
summary(Mod2_pfk13)
plot(Mod2_pfk13)
qqnorm(residuals(Mod2_pfk13))

Mod2_pfk13 <- lmer(Prev_DR ~ PR +drug_5yr_ago+ (1|country), weights=log(tested),data=k13_na_rem)
summary(Mod2_pfk13)
plot(Mod2_pfk13)
qqnorm(residuals(Mod2_pfk13))

ggplot(k13_na_rem, aes(drug_4yr_ago,Prev_DR,col=country))+geom_point(aes(size=log(tested)))

ggplot(k13_na_rem, aes(drug_5yr_ago,Prev_DR,col=country))+geom_point(aes(size=log(tested)))

ggplot(k13_na_rem, aes(cut(drug_4yr_ago,10),Prev_DR,col=country))+geom_boxplot()

ggplot(k13_na_rem, aes(cut(drug_5yr_ago,10),Prev_DR,col=country))+geom_boxplot()



