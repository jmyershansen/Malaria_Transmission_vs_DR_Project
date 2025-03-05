install.packages("dlnm")
install.packages("tidyverse")
install.packages("splines")
install.packages("lme4")
library(dlnm)
library(tidyverse)
library(splines)
library(lme4)


install.packages ("betareg")
library(betareg)
install.packages ("zoo")
library(zoo)
install.packages ("lmtest")
library(lmtest)

remove.packages("Matrix")
install.packages("Matrix")
remove.packages("Matrix")
remove.packages("lme4")
install.packages("lme4", type = "source")
library(lme4)
###1
df <- read_csv("pfcrt_lagdata.csv")

cq.cb <- crossbasis(df$CQ_use, lag=12, argvar=list(fun="lin"),
                    arglag=list(fun="poly",degree=2))

pr.cb <- crossbasis(df$PR, lag=12, argvar=list(fun="lin"),
                    arglag=list(fun="poly",degree=2))
#summary(cq.cb)
#df$`study Id` <- as.factor(df$`study Id`)

model1 <- glm(Prev_DR ~ cq.cb + pr.cb + ns(year, 1) + year,
              family=quasipoisson(), data=df)
summary(model1)
plot(model1)

model2 <- lmer(Prev_DR ~ cq.cb + pr.cb + (1 | `study Id`), data=df)
summary(model2)
plot(model2)

pred1.cq <- crosspred(cq.cb, model1, at=0:5, bylag=1, cumul=TRUE)
#pred1.pr <- crosspred(pr.cb, model1, at=0:5, bylag=1, cumul=TRUE)

plot(pred1.cq, "slices", var=2, col=3, ylab="Prev_DR", ci.arg=list(density=15,lwd=10),
     main="CQ use on Pfcrt-K76T")

plot(pred1.cq, "slices", var=2, col=2, cumul=TRUE, ylab="Cumulative Prev_DR",
     main="CQ use on Pfcrt-K76T@ 2-unit increase in CQ")

pred1.cq$allRRfit["2"]
cbind(pred1.cq$allRRlow, pred1.cq$allRRhigh)["2",]

#plot(pred1.pr, "slices", var=2, col=3, ylab="Prev_DR", ci.arg=list(density=15,lwd=2),
main="PR on Prev DR")

###2

df <- read_csv("pfmdr1_lagdata.csv")

cq.cb <- crossbasis(df$CQ_use, lag=10, argvar=list(fun="lin"),
                    arglag=list(fun="poly",degree=2))

pr.cb <- crossbasis(df$PR, lag=10, argvar=list(fun="lin"),
                    arglag=list(fun="poly",degree=2))

#summary(cq.cb)
#pfcrt_SNP1$study.Id <- as.factor(pfcrt_SNP1$study.Id)

model1 <- glm(Prev_DR ~ cq.cb + pr.cb + ns(year, 1) + year,
              family=quasipoisson(), data=df)
summary(model1)

#model2 <- lmer(Prev_DR ~ cq.cb + pr.cb + (1|study Id), data=dothisone)
#summary(model2)

pred1.cq <- crosspred(cq.cb, model1, at=0:5, bylag=1, cumul=TRUE)
#pred1.pr <- crosspred(pr.cb, model1, at=0:5, bylag=1, cumul=TRUE)

plot(pred1.cq, "slices", var=2, col=3, ylab="Prev_DR", ci.arg=list(density=15,lwd=10),
     main="CQ use on Pfmdr1-N86Y")

plot(pred1.cq, "slices", var=2, col=2, cumul=TRUE, ylab="Cumulative Prev_DR",
     main="CQ use on Pfmdr1-N86Y@ 2-unit increase in CQ")

pred1.cq$allRRfit["2"]
cbind(pred1.cq$allRRlow, pred1.cq$allRRhigh)["2",]

###3

df <- read_csv("pfdhps_lagdata.csv")

sp.cb <- crossbasis(df$SP_use, lag=8, argvar=list(fun="lin"),
                    arglag=list(fun="poly",degree=2))

pr.cb <- crossbasis(df$PR, lag=8, argvar=list(fun="lin"),
                    arglag=list(fun="poly",degree=2))

#summary(cq.cb)
#pfcrt_SNP1$study.Id <- as.factor(pfcrt_SNP1$study.Id)

model1 <- glm(Prev_DR ~ sp.cb + pr.cb + ns(year, 1) + year,
              family=quasipoisson(), data=df)
summary(model1)

#model2 <- lmer(Prev_DR ~  sp.cb + pr.cb + (1|study Id), data=dothisone)
#summary(model2)

pred1.sp <- crosspred(sp.cb, model1, at=0:5, bylag=1, cumul=TRUE)
#pred1.pr <- crosspred(pr.cb, model1, at=0:5, bylag=1, cumul=TRUE)

plot(pred1.sp, "slices", var=2, col=3, ylab="Prev_DR", ci.arg=list(density=15,lwd=10),
     main="SP use on Pfdhps-A437G")

plot(pred1.sp, "slices", var=2, col=2, cumul=TRUE, ylab="Cumulative Prev_DR",
     main="SP use on Pfdhps-A437G@ 2-unit increase in CQ")

pred1.cq$allRRfit["2"]
cbind(pred1.cq$allRRlow, pred1.cq$allRRhigh)["2",]

###4

df <- read_csv("pfdhfr_lagdata.csv")

sp.cb <- crossbasis(df$SP_use, lag=8, argvar=list(fun="lin"),
                    arglag=list(fun="poly",degree=2))

pr.cb <- crossbasis(df$PR, lag=8, argvar=list(fun="lin"),
                    arglag=list(fun="poly",degree=2))

#summary(cq.cb)
#pfcrt_SNP1$study.Id <- as.factor(pfcrt_SNP1$study.Id)

model1 <- glm(Prev_DR ~ sp.cb + pr.cb + ns(year, 1) + year,
              family=quasipoisson(), data=df)
summary(model1) nbinom()

#Model <- glm(Prev_DR ~ year + PR + CQ_use, data = pfcrt_SNP1, family = gaussian)
#summary(Model)

#model2 <- lmer(Prev_DR ~ sp.cb + pr.cb + (1|study Id), data=dothisone)
#summary(model2)

pred1.sp <- crosspred(sp.cb, model1, at=0:5, bylag=1, cumul=TRUE)
pred1.pr <- crosspred(pr.cb, model1, at=0:5, bylag=1, cumul=TRUE)

plot(pred1.sp, "slices", var=2, col=3, ylab="Prev_DR", ci.arg=list(density=15,lwd=10),
     main="SP use on Pfdhfr-S108N")

plot(pred1.sp, "slices", var=2, col=2, cumul=TRUE, ylab="Cumulative Prev_DR",
     main="SP use on Pfdhfr-S108N@ 2-unit increase in CQ")

pred1.cq$allRRfit["2"]
cbind(pred1.cq$allRRlow, pred1.cq$allRRhigh)["2",]
#####
#####
#####

df <- read_csv("pfcrt_prev&pr&cq.csv")
df <- df %>% mutate(mean_PR = mean_PR / 100, mean_Prev = mean_Prev / 100)
df$mean_PR <- round (df$mean_PR, digits =2)
df$mean_Prev <- round (df$mean_Prev, digits =2)
df$mean_cq <- round (df$mean_cq, digits =2)

df <- df[!is.na(df$mean_Prev),]
df <- df[!is.na(df$mean_PR),]


cq.cb <- crossbasis(df$mean_cq, lag=12, argvar=list(fun="lin"),
                    arglag=list(fun="poly",degree=2))

pr.cb <- crossbasis(df$mean_PR, lag=12, argvar=list(fun="lin"),
                    arglag=list(fun="poly",degree=2))

summary(cq.cb)
#df$`country` <- as.factor(df$`country`)

model1 <- glm(mean_Prev ~ cq.cb + pr.cb + ns(year, 1) + year,
              family=quasipoisson(), data=df)
summary(model1)
plot(model1)

model2 <- lmer(mean_Prev ~ cq.cb + pr.cb + (1 | `country`), data=df)
summary(model2)
#plot(model2)

install.packages ("betareg")
library(betareg)
install.packages ("zoo")
library(zoo)
install.packages ("lmtest")
library(lmtest)
data (GasolineYield)
head(GasolineYield,3)
View(GasolineYield)
##
data (FoodExpenditure)
head(FoodExpenditure,3)
View(FoodExpenditure)
#data("GasolineYield", package = "betareg")

###
df <- read_csv("pfcrt_prev&pr&cq.2.csv")
df <- df %>% mutate(PR = PR / 100, Prev_DR = Prev_DR / 100)
# changing 0 and 1
df2 = df
df2$Prev_DR = ifelse(df2$Prev_DR==1, 0.999, df2$Prev_DR)
df2$Prev_DR = ifelse(df2$Prev_DR==0, 0.0001, df2$Prev_DR)

gy_logit <- betareg(Prev_DR ~ mean_cq + PR, data = df2)
summary(gy_logit)
plot(gy_loglog)

gy_logit2 <- betareg(Prev_DR ~ mean_cq + PR | PR, data = df2)
summary(gy_logit2)

gy_loglog <- betareg(Prev_DR ~ mean_cq + PR, data = df2, link = "loglog")
summary(gy_loglog)
plot(gy_loglog)

###Pfcrt_76T Prev vs PR
#ggplot(data=df2, aes(x = mean_PR, y = mean_Prev) + geom_point(aes(color=as.factor(Country), size=tested))+ geom_smooth(method = "lm") + xlab("Parasite Rate") + ylab("Pfcrt 76T Prevalence")
#ggplot(data=df2 %>% filter(Region %in% c("West Africa","Central Africa", "East Africa", "Southern Africa", "SEA", "South America")))+ aes(x = PR, y = Prev_DR) + geom_point(aes(color=as.factor(Region), size=tested))+ geom_smooth(method = "lm") + xlab("Parasite Rate") + ylab("Pfcrt 76T Prevalence")
ggplot(data = df2, aes (y = Prev_DR, x = PR)) + geom_point(aes(size = 3))  + theme_classic() + geom_smooth(method = "lm") + xlab ("PR") + ylab ("Prevalence")


###
model2 <- lmer(Prev_DR ~ mean_cq + PR + (1 | `country`), data=df)
summary(model2)
plot(model2)

model2 <- glm(Prev_DR ~ mean_cq + PR, data=df)
summary(model2)
plot(model2)

pred1.cq <- crosspred(cq.cb, model1, at=0:5, bylag=1, cumul=TRUE)
#pred1.pr <- crosspred(pr.cb, model1, at=0:5, bylag=1, cumul=TRUE)

plot(pred1.cq, "slices", var=2, col=3, ylab="Prev_DR", ci.arg=list(density=15,lwd=10),
     main="CQ use on Pfcrt-K76T")

plot(pred1.cq, "slices", var=2, col=2, cumul=TRUE, ylab="Cumulative Prev_DR",
     main="CQ use on Pfcrt-K76T@ 2-unit increase in CQ")

pred1.cq$allRRfit["2"]
cbind(pred1.cq$allRRlow, pred1.cq$allRRhigh)["2",]
#####
#####
#####
df <- read_csv("pfmdr1_prev&pr&cq.csv")
df <- df %>% mutate(PR = PR / 100, Prev_DR = Prev_DR / 100)
# changing 0 and 1
df2 = df
df2$Prev_DR = ifelse(df2$Prev_DR==1, 0.999, df2$Prev_DR)
df2$Prev_DR = ifelse(df2$Prev_DR==0, 0.0001, df2$Prev_DR)

gy_logit <- betareg(Prev_DR ~ mean_cq + PR, data = df2)
summary(gy_logit)
plot(gy_logit)

#gy_loglog <- betareg(Prev_DR ~ mean_cq + PR, data = df2, link = "loglog")
#summary(gy_loglog)
#plot(gy_loglog)

model2 <- glm(Prev_DR ~ mean_cq + PR, data=df)
summary(model2)
plot(model2)

#####
#####
#####
df <- read_csv("pfdhps_prev&pr&sp.csv")
df <- df %>% mutate(PR = PR / 100, Prev_DR = Prev_DR / 100)
# changing 0 and 1
df2 = df
df2$Prev_DR = ifelse(df2$Prev_DR==1, 0.999, df2$Prev_DR)
df2$Prev_DR = ifelse(df2$Prev_DR==0, 0.0001, df2$Prev_DR)

gy_logit <- betareg(Prev_DR ~ mean_sp + PR, data = df2)
summary(gy_logit)
plot(gy_logit)

gy_loglog <- betareg(Prev_DR ~ mean_sp + PR, data = df2, link = "loglog")
summary(gy_loglog)
plot(gy_loglog)

#####
#####
#####
df <- read_csv("pfdhfr_prev&pr&sp.csv")
df <- df %>% mutate(PR = PR / 100, Prev_DR = Prev_DR / 100)
# changing 0 and 1
df2 = df
df2$Prev_DR = ifelse(df2$Prev_DR==1, 0.999, df2$Prev_DR)
df2$Prev_DR = ifelse(df2$Prev_DR==0, 0.0001, df2$Prev_DR)

gy_logit <- betareg(Prev_DR ~ mean_sp + PR, data = df2)
summary(gy_logit)
plot(gy_logit)

gy_loglog <- betareg(Prev_DR ~ mean_sp + PR, data = df2, link = "loglog")
summary(gy_loglog)
plot(gy_loglog)

#####
#####
#####
df <- read_csv("pfk13_prev&pr&act.csv")
df <- df %>% mutate(PR = PR / 100, Prev_DR = Prev_DR / 100)
# changing 0 and 1
df2 = df
df2$Prev_DR = ifelse(df2$Prev_DR==1, 0.999, df2$Prev_DR)
df2$Prev_DR = ifelse(df2$Prev_DR==0, 0.0001, df2$Prev_DR)

gy_logit <- betareg(Prev_DR ~ mean_act + PR, data = df2)
summary(gy_logit)
plot(gy_logit)

gy_loglog <- betareg(mean_Prev ~ mean_act + mean_PR, data = df2, link = "loglog")
summary(gy_loglog)
plot(gy_loglog)

model3 <- lmer(mean_Prev ~ mean_cq + mean_PR + (1 | `country`), data=df)
summary(model3)
plot(model3)


####
#####
#####
###Qixin's code

df <- read_csv("pfcrt_prev&pr&cq.csv")
df <- df %>% mutate(mean_PR = mean_PR / 100, mean_Prev = mean_Prev / 100)
#df$mean_PR <- round (df$mean_PR, digits =2)
#df$mean_Prev <- round (df$mean_Prev, digits =2)
#df$mean_cq <- round (df$mean_cq, digits =2)
#df <- df[!is.na(df$mean_Prev),]
#df <- df[!is.na(df$mean_PR),]

?lag

head(df)

df<-df%>%group_by(country)%>%mutate(drug_1yr_ago = lag(mean_cq,n=1))
df<-df%>%group_by(country)%>%mutate(drug_2yr_ago = lag(mean_cq,n=2),
                                    drug_3yr_ago = lag(mean_cq,n=3))
df<-df%>%group_by(country)%>%mutate(drug_4yr_ago = lag(mean_cq,n=4),
                                    drug_5yr_ago = lag(mean_cq,n=5))

df<- df %>% select(-mean_Prev, -mean_PR)
write_csv(df, file ="CQ_lag.csv")
  

# linear

model1<-lm(mean_Prev ~ mean_PR+mean_cq+drug_1yr_ago+drug_2yr_ago+drug_3yr_ago,data=df)
summary(model1)

model2<-lmer(mean_Prev ~ mean_PR+(mean_PR|country)+
             (drug_1yr_ago|country)+drug_2yr_ago+drug_3yr_ago,data=df)
summary(model2)

model3<-lmer(mean_Prev ~ drug_3yr_ago+(drug_2yr_ago|country)+
               (drug_3yr_ago|country),data=df)
summary(model3)

model4<-lmer(mean_Prev ~ drug_1yr_ago+drug_3yr_ago+drug_5yr_ago+
               (drug_3yr_ago|country)+(drug_5yr_ago|country),data=df)
summary(model4)
plot(model4)

ggplot(df, aes(drug_5yr_ago,mean_Prev))+geom_point()



#####
#####
#####
cq.cb <- crossbasis(df$mean_cq, lag=12, argvar=list(fun="lin"),
                    arglag=list(fun="poly",degree=2))

pr.cb <- crossbasis(df$mean_PR, lag=12, argvar=list(fun="lin"),
                    arglag=list(fun="poly",degree=2))

summary(cq.cb)
#df$`country` <- as.factor(df$`country`)

model1 <- glm(mean_Prev ~ cq.cb + pr.cb + ns(year, 1) + year,
              family=quasipoisson(), data=df)
summary(model1)
plot(model1)

model2 <- lmer(mean_Prev ~ cq.cb + pr.cb + (1 | `country`), data=df)
summary(model2)
#plot(model2)

####
#####
#####

df <- read_csv("pfmdr1_prev&pr&cq.csv")
df <- df %>% mutate(mean_PR = mean_PR / 100, mean_Prev = mean_Prev / 100)
#df$mean_PR <- round (df$mean_PR, digits =2)
#df$mean_Prev <- round (df$mean_Prev, digits =2)
#df$mean_cq <- round (df$mean_cq, digits =2)
df <- df[!is.na(df$mean_Prev),]
df <- df[!is.na(df$mean_PR),]


cq.cb <- crossbasis(df$mean_cq, lag=12, argvar=list(fun="lin"),
                    arglag=list(fun="poly",degree=2))

pr.cb <- crossbasis(df$mean_PR, lag=12, argvar=list(fun="lin"),
                    arglag=list(fun="poly",degree=2))

summary(cq.cb)
#df$`country` <- as.factor(df$`country`)

model1 <- glm(mean_Prev ~ cq.cb + pr.cb + ns(year, 1) + year,
              family=quasipoisson(), data=df)
summary(model1)
plot(model1)

model2 <- lmer(mean_Prev ~ cq.cb + pr.cb + (1 | `country`), data=df)
summary(model2)
#plot(model2)

df <- read_csv("pfdhps_prev&pr&sp.csv")
df <- df %>% mutate(mean_PR = mean_PR / 100, mean_Prev = mean_Prev / 100)
#df$mean_PR <- round (df$mean_PR, digits =2)
#df$mean_Prev <- round (df$mean_Prev, digits =2)
#df$mean_cq <- round (df$mean_cq, digits =2)
df <- df[!is.na(df$mean_Prev),]
df <- df[!is.na(df$mean_PR),]


sp.cb <- crossbasis(df$mean_sp, lag=12, argvar=list(fun="lin"),
                    arglag=list(fun="poly",degree=2))

pr.cb <- crossbasis(df$mean_PR, lag=12, argvar=list(fun="lin"),
                    arglag=list(fun="poly",degree=2))

summary(sp.cb)
#df$`country` <- as.factor(df$`country`)

model1 <- glm(mean_Prev ~ sp.cb + pr.cb + ns(year, 1) + year,
              family=quasipoisson(), data=df)
summary(model1)
plot(model1)

model2 <- lmer(mean_Prev ~ sp.cb + pr.cb + (1 | `country`), data=df)
summary(model2)
#plot(model2)

df <- read_csv("pfdhfr_prev&pr&sp.csv")
df <- df %>% mutate(mean_PR = mean_PR / 100, mean_Prev = mean_Prev / 100)
#df$mean_PR <- round (df$mean_PR, digits =2)
#df$mean_Prev <- round (df$mean_Prev, digits =2)
#df$mean_cq <- round (df$mean_cq, digits =2)
df <- df[!is.na(df$mean_Prev),]
df <- df[!is.na(df$mean_PR),]


sp.cb <- crossbasis(df$mean_sp, lag=12, argvar=list(fun="lin"),
                    arglag=list(fun="poly",degree=2))

pr.cb <- crossbasis(df$mean_PR, lag=12, argvar=list(fun="lin"),
                    arglag=list(fun="poly",degree=2))

summary(sp.cb)
#df$`country` <- as.factor(df$`country`)

model1 <- glm(mean_Prev ~ sp.cb + pr.cb + ns(year, 1) + year,
              family=quasipoisson(), data=df)
summary(model1)
plot(model1)

model2 <- lmer(mean_Prev ~ sp.cb + pr.cb + (1 | `country`), data=df)
summary(model2)
#plot(model2)