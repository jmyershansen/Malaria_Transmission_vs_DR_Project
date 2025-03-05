###PFCRT
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

###Model
Mod2_pfcrt <- lmer(Prev_DR ~ PR +mean_cq+PR *mean_cq +drug_1yr_ago+drug_2yr_ago+drug_3yr_ago+(1|country), weights=log(tested),data=crt_na_rem)
summary(Mod2_pfcrt)
# homogeneity of variance
plot(Mod2_pfcrt)
# normality of the residuals
qqnorm(residuals(Mod2_pfcrt))

###Interpretation


##Plots
crt_na_rem<- crt_na_rem%>% mutate(PR = PR * 100, Prev_DR = Prev_DR * 100)
#### 1. Pfcrt_76T Prevalence vs Mean_cq_All regions
ggplot(crt_na_rem, aes(cut(mean_cq,10),Prev_DR,col=Region))+geom_boxplot()
ggsave("pfcrt_76T_p1.png") 

### 2. Pfcrt_76T Prevalence vs PR_All regions
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


###PFMDR1
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

###Interpretation


##Plots
mdr_na_rem<- mdr_na_rem%>% mutate(PR = PR * 100, Prev_DR = Prev_DR * 100)
#### 1. Pfmdr1_86Y Prevalence vs Mean_cq_All regions
ggplot(mdr_na_rem, aes(cut(mean_cq,10),Prev_DR,col=Region))+geom_boxplot()
ggsave("pfmdr1_86Y_p1.png") 

### 2. Pfmdr1_86Y Prevalence vs PR_All regions
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


###PFDHPS
sp<- read_csv("SP_imputed.csv")
#Set up lag years
sp<-sp%>%group_by(country)%>%mutate(drug_1yr_ago = lag(mean_sp,n=1),drug_2yr_ago = lag(mean_sp,n=2),
                                    drug_3yr_ago = lag(mean_sp,n=3),drug_4yr_ago = lag(mean_sp,n=4),drug_5yr_ago = lag(mean_sp,n=5))

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

###Interpretation


##Plots
dhps_na_rem<- dhps_na_rem%>% mutate(PR = PR * 100, Prev_DR = Prev_DR * 100)
#### 1. Pfdhps_437G Prevalence vs Mean_sp_All regions
ggplot(dhps_na_rem, aes(cut(mean_sp,10),Prev_DR,col=Region))+geom_boxplot()
ggsave("pfdhps_437G _p1.png")

### 2. Pfdhps_437G  Prevalence vs PR_All regions
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


###PFDHFR
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

###Interpretation


###Plots
dhfr_na_rem <- dhfr_na_rem%>% mutate(PR = PR * 100, Prev_DR = Prev_DR * 100)
#### 1. Pfdhfr_108N Prevalence vs Mean_sp_All regions
ggplot(dhfr_na_rem, aes(cut(mean_sp,10),Prev_DR,col=Region))+geom_boxplot()
ggsave("pfdhfr_108N_p1.png") 

### 2. Pfdhfr_108N Prevalence vs PR_All regions
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

###Model
k13_na_rem <-k13_na_rem %>% filter(Prev_DR<1,Prev_DR>0)
Mod1_pfk13 <- lmer(Prev_DR ~ PR +mean_act+PR *mean_act +drug_1yr_ago+drug_2yr_ago+drug_3yr_ago+ (1|country), weights=log(tested),data=k13_na_rem)
summary(Mod1_pfk13)
plot(Mod1_pfk13)
qqnorm(residuals(Mod1_pfk13))

###Interpretation



###Plots
k13_na_rem<- k13_na_rem %>% mutate(PR = PR * 100, Prev_DR = Prev_DR * 100)
#### 1. Pfk13_580Y Prevalence vs Mean_act_All regions
ggplot(k13_na_rem, aes(cut(mean_act,10),Prev_DR,col=country))+geom_boxplot()
ggsave("pfK13_580Y_p1.png") 

### 2. Pfdhfr_108N Prevalence vs PR_All regions
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
