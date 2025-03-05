## looking at data distribution
library(tidyverse)

# load prfct_snp
#df = read_csv('../Data.1/pfcrt_SNP_complete.csv')
df <- read_csv("pfcrt_SNP_complete.csv")

# must have PR, CQ_use, must have Prev_DR
df = df[!is.na(df$PR),] # over 100 lack PR
df = df[!is.na(df$CQ_use),] # over 900 lack CQ_use
df = df[!is.na(df$Prev_DR),]

df = df[!duplicated(df),]

df <- df %>% distinct(`country`,`CQ_use`,`year`)
df  <- as.data.frame(cbind(df$country,df$CQ_use,df$year))
colnames(df)<-c("country","CQ_use","year")
df<- df %>% arrange(country,year)
NStudies <- length(unique(df$`country`))
write_csv(df,file="pfcrt_drug.years.csv")

#length(unique(df$`study Id`))
#df <- df %>% mutate(update_year=(`study start` + `study end`)/2)
#df$update_year <- round(df$update_year, digits = 0)
#df$update_year = (df$`study start` + df$`study end`) / 2
df <- df %>% select(-year)
df <- df %>% mutate(year=(update_year))

#DUdata <- read_csv ("FeveredChildrenDrugTreatProp_CQ.0.csv")
#df1 <- inner_join(DUdata,df, by = c("country","year"))

#df2 = df %>% group_by(country) %>% 
  #summarise(min_year = min(year, na.rm = T),
            #max_year = max(year, na.rm = T),
            #unique_year = length(unique(year)),
            #total_rows = n(),
            #middle = mean(c(min_year, max_year))) %>%
  #arrange(min_year)

# middle year
#df2$middle_year = 0

# this loop finds closest middle year
#for(i in 1:nrow(df2)){
  #ro = which(df$country == df2$country[i])
  #diffs = df$year[ro] - df2$middle[i]
  #diffs = sqrt(diffs ** 2)
  #x = which(diffs == min(diffs))
  #df2$middle_year[i] = df$year[ro[x]][1] %>% as.numeric()
}

#df2 = df2 %>% select(country, min_year, middle_year, max_year, everything())
#df2 = df2 %>% pivot_longer(cols = 2:4, names_to = "Year_Type", values_to = "Year")


### adding 11/1/23
# to make it I need at least two years and least two studies
madeit = df2 %>% filter(unique_year > 1, total_rows > 1)

by_country = list()
for(i in 1:nrow(madeit)){
  by_country[[i]] = df[df$country == madeit$country[i],]
}

## try this lag-model
library(dlnm)
library(splines)

dothisone = by_country %>% bind_rows() %>% arrange(country,year)%>% group_by(country,year) # gathers all the data 6:48
dothisone = df
#dothisone = by_country[[3]] ## choose a country
#length(unique(dothisone$`study Id`))
#length(unique(dothisone$`country`))

duplicated(df)
which(duplicated(df))
which(duplicated(df))->hi
df[hi,]
df = df[!duplicated(df),]
#pfcrt_clean <-df

duplicated(dothisone)
which(duplicated(dothisone))
which(duplicated(dothisone))->hi
dothisone[hi,]
dothisone = dothisone[!duplicated(dothisone),]

write_csv(dothisone, file = "pfcrt_lagdata.csv")
write_csv(dothisone, file = "pfcrt_lagdata.0.csv")

#joined <- right_join(dothisone,df,by = c("PR","Prev_DR","CQ_use"))
#joined <- right_join(dothisone,df,by = c("country","year"))

df <- read_csv("pfcrt_lagdata.csv")
pr.cb <- crossbasis(dothisone$PR, lag=3, argvar=list(fun="lin"),
                     arglag=list(fun="poly",degree=3))

cq.cb <- crossbasis(dothisone$CQ_use, lag=3, argvar=list(fun="lin"),
                    arglag=list(fun="poly",degree=3))
#summary(cq.cb)

pfcrt_SNP1$study.Id <- as.factor(pfcrt_SNP1$study.Id)

model1 <- glm(Prev_DR ~ pr.cb + cq.cb + ns(year, 1) + 0,
              family=quasipoisson(), dothisone)

model2 <- lmer(Prev_DR ~  pr.cb + cq.cb+ (1|study Id), data=dothisone)
summary(model1)

pred1.cq <- crosspred(cq.cb, model1, at=0:5, bylag=1, cumul=TRUE) # can't
pred1.pr <- crosspred(pr.cb, model1, at=0:5, bylag=1, cumul=TRUE)

plot(pred1.cq, "slices", var=5, col=3, ylab="Prev_DR", ci.arg=list(density=22,lwd=10),
   main="CQ use on Prev DR")

plot(pred1.pr, "slices", var=1, col=3, ylab="Prev_DR", ci.arg=list(density=15,lwd=2),
    main="PR on Prev DR")


# list of parameters we care about:
# plot(what is var =)
# crosspred(at =, bylag =, cumul =,)
# glm(ns, family)
# crossbasis(lag, fun='lin', fun = 'poly', degree)

# just for fun
ggplot(dothisone, aes(year))+
  geom_point(aes(y = CQ_use), color = 'green')+
  geom_smooth(aes(y = CQ_use), color = 'green', method = 'lm', se = F)+
  geom_point(aes(y = PR), color = 'red')+
  geom_smooth(aes(y = PR), color = 'red', method = 'lm', se = F)+
  geom_point(aes(y=Prev_DR), color = 'blue')+
  geom_smooth(aes(y=Prev_DR), color = 'blue', method = 'lm', se = F)+
  facet_wrap(~Region)

###

#write_csv(df2, file = "pfcrt_SNP_df.csv", row.names=F)
#unique(df$country) %>% sort()
####
####
#Check for studies

df <- read_csv("pfcrt_SNP_complete.csv")
df = df[!is.na(df$PR),] 
df = df[!is.na(df$CQ_use),] 
df = df[!is.na(df$Prev_DR),]
length(unique(df$`study Id`))
length(unique(df$`country`))
write_csv(df, file = "pfcrt_SNP_complete.csv")

df <- read_csv("pfmdr1_SNP_complete.csv")
df = df[!is.na(df$PR),] 
df = df[!is.na(df$CQ_use),] 
df = df[!is.na(df$Prev_DR),]
length(unique(df$`study Id`))
length(unique(df$`country`))
write_csv(df, file = "pfmdr1_SNP_complete.csv")

df <- read_csv("pfdhfr_SNP_complete.csv")
df = df[!is.na(df$PR),]
df = df[!is.na(df$SP_use),] 
df = df[!is.na(df$Prev_DR),]
length(unique(df$`study id`))
length(unique(df$`country`))
write_csv(df, file = "pfdhfr_SNP_complete.csv")

df <- read_csv("pfdhps_SNP_complete.csv")
df = df[!is.na(df$PR),]
df = df[!is.na(df$SP_use),] 
df = df[!is.na(df$Prev_DR),]
length(unique(df$`study id`))
length(unique(df$`country`))
write_csv(df, file = "pfdhps_SNP_complete.csv")

df <- read_csv("pfK13_SNP_complete.csv")
df = df[!is.na(df$PR),]
df = df[!is.na(df$ACT_use),] 
df = df[!is.na(df$Prev_DR),]
length(unique(df$`sid`))
length(unique(df$`country`))
write_csv(df, file = "pfK13_SNP_complete.csv")
#####
#####
df2$country = factor(df2$country, levels = unique(df2$country))

df2 %>% filter(unique_year > 1) %>% 
ggplot(., aes(x = min_year, xend = max_year, 
                y = country, yend = country, color = unique_year)) +
  geom_segment(size = 2) +
  geom_point()+
  labs(x = "Year Range", y = "Country")+
  scale_color_viridis_b()+
  theme_dark()

summary(df2$unique_year)
quantile(df2$unique_year, probs = seq(0,1,0.1))
table(df2$unique_year)

table(is.na(df$year), df$country)

table(df2$total_rows)
quantile(df2$total_rows, probs = seq(0,1,0.1))


## see pr vs dr_prev
ggplot(df, aes(PR, Prev_DR))+
  geom_point(show.legend = F)+facet_wrap(~year)+
  geom_smooth(se = F)

df %>% filter(country == "Mali") %>% 
ggplot(., aes(year, Prev_DR))+
  geom_smooth(color = 'blue', se = F)+
  geom_point(fill = 'blue', shape = 21)+
  facet_wrap(~country)+
  geom_smooth(aes(y = PR), color = 'red', se = F)+
  geom_point(aes(y = PR), fill = 'red', shape = 21)+
  geom_smooth(aes(y = CQ_use), color = 'green', se = F)+
  geom_point(aes(y = CQ_use), fill = 'green', shape = 21)

glm(Prev_DR ~ PR + CQ_use , data = df)

colnames(df)
####
####


## MDR1
df <- read_csv("pfmdr1_SNP_complete.csv")

# must have PR, CQ_use, must have Prev_DR
df = df[!is.na(df$PR),] # over 100 lack PR
df = df[!is.na(df$CQ_use),] # over 900 lack CQ_use
df = df[!is.na(df$Prev_DR),]

# filter out duplicates
df = df[!duplicated(df),]

duplicated(df)
which(duplicated(df))
which(duplicated(df))->hi
df[hi,]
df = df[!duplicated(df),]
pfmdr1_clean <-df

df <- df %>% distinct(`country`,`CQ_use`,`year`)
df  <- as.data.frame(cbind(df$country,df$CQ_use,df$year))
colnames(df)<-c("country","CQ_use","year")
df<- df %>% arrange(country,year)
NStudies <- length(unique(df$`country`))
write_csv(df,file="pfmdr1_drug.years.csv")


# count the number of years and studies per country
df2 = df %>% group_by(country) %>% 
  summarise(unique_year = length(unique(year)),
            total_rows = n())

### adding 11/1/23
# to make it I need at least two years and least two studies
madeit = df2 %>% filter(unique_year > 1, total_rows > 1)

# by country (list of dataframes)
by_country = list()
for(i in 1:nrow(madeit)){
  by_country[[i]] = df[df$country == madeit$country[i],]
}

# by regions (just know that this includes countries with only 1 study)
regions = df %>% select(Region) %>% unique() %>% unlist()
by_region = list()
for(i in regions){
  by_region[[i]] = df[df$Region == i,]
}

## try this lag-model
library(dlnm)
library(splines)

dothisone = by_country %>% bind_rows() %>% arrange(country,year)%>% group_by(country,year) # gathers all the data 6:48
dothisone = by_region %>% bind_rows() %>% arrange(country,year)%>% group_by(country,year) # gathers all the data 6:48

#unique(dothisone$Region)
#unique(df$Region)

write_csv(dothisone, file = "pfmdr1_lagdata.csv")
write_csv(dothisone, file = "pfmdr1_lagdata.0.csv")

####
####

## DHPS

df <- read_csv("pfdhps_SNP_complete.csv")

# must have PR, CQ_use, must have Prev_DR
df = df[!is.na(df$PR),] # over 100 lack PR
df = df[!is.na(df$SP_use),] # over 900 lack CQ_use
df = df[!is.na(df$Prev_DR),]

# filter out duplicates
df = df[!duplicated(df),]

duplicated(df)
which(duplicated(df))
which(duplicated(df))->hi
df[hi,]
df = df[!duplicated(df),]
pfdhps_clean <-df

df <- df %>% distinct(`country`,`SP_use`,`year`)
df  <- as.data.frame(cbind(df$country,df$SP_use,df$year))
colnames(df)<-c("country","SP_use","year")
df<- df %>% arrange(country,year)
NStudies <- length(unique(df$`country`))
write_csv(df,file="pfdhps_drug.years.csv")

# count the number of years and studies per country
df2 = df %>% group_by(country) %>% 
  summarise(unique_year = length(unique(year)),
            total_rows = n())

### adding 11/1/23
# to make it I need at least two years and least two studies
madeit = df2 %>% filter(unique_year > 1, total_rows > 1)


# by country (list of dataframes)
by_country = list()
for(i in 1:nrow(madeit)){
  by_country[[i]] = df[df$country == madeit$country[i],]
}

# by regions (just know that this includes countries with only 1 study)
regions = df %>% select(Region) %>% unique() %>% unlist()
by_region = list()
for(i in regions){
  by_region[[i]] = df[df$Region == i,]
}

## try this lag-model
library(dlnm)
library(splines)

dothisone = by_country %>% bind_rows() %>% arrange(country,year)%>% group_by(country,year) # gathers all the data 6:48
dothisone = by_region %>% bind_rows() %>% arrange(country,year)%>% group_by(country,year) # gathers all the data 6:48
dothisone = df
#unique(dothisone$Region)
#unique(df$Region)

write_csv(dothisone, file = "pfdhps_lagdata.csv")
write_csv(dothisone, file = "pfdhps_lagdata.0.csv")

###
###

## DHFR

df <- read_csv("pfdhfr_SNP_complete.csv")

# must have PR, CQ_use, must have Prev_DR
df = df[!is.na(df$PR),] # over 100 lack PR
df = df[!is.na(df$SP_use),] # over 900 lack CQ_use
df = df[!is.na(df$Prev_DR),]

# filter out duplicates
df = df[!duplicated(df),]

duplicated(df)
which(duplicated(df))
which(duplicated(df))->hi
df[hi,]
df = df[!duplicated(df),]
pfdhfr_clean <-df

df <- df %>% distinct(`country`,`SP_use`,`year`)
df  <- as.data.frame(cbind(df$country,df$SP_use,df$year))
colnames(df)<-c("country","SP_use","year")
df<- df %>% arrange(country,year)
NStudies <- length(unique(df$`country`))
write_csv(df,file="pfdhfr_drug.years.csv")


# count the number of years and studies per country
df2 = df %>% group_by(country) %>% 
  summarise(unique_year = length(unique(year)),
            total_rows = n())

### adding 11/1/23
# to make it I need at least two years and least two studies
madeit = df2 %>% filter(unique_year > 1, total_rows > 1)

# by country (list of dataframes)
by_country = list()
for(i in 1:nrow(madeit)){
  by_country[[i]] = df[df$country == madeit$country[i],]
}

# by regions (just know that this includes countries with only 1 study)
regions = df %>% select(Region) %>% unique() %>% unlist()
by_region = list()
for(i in regions){
  by_region[[i]] = df[df$Region == i,]
}

## try this lag-model
library(dlnm)
library(splines)

dothisone = by_country %>% bind_rows() %>% arrange(country,year)%>% group_by(country,year) # gathers all the data 6:48
dothisone = by_region %>% bind_rows() %>% arrange(country,year)%>% group_by(country,year) # gathers all the data 6:48

#unique(dothisone$Region)
#unique(df$Region)

write_csv(dothisone, file = "pfdhfr_lagdata.csv")
write_csv(dothisone, file = "pfdhfr_lagdata.0.csv")

###
###

## K13

df <- read_csv("pfK13_SNP_complete.csv")

# must have PR, CQ_use, must have Prev_DR
df = df[!is.na(df$PR),] # over 100 lack PR
df = df[!is.na(df$ACT_use),] # over 900 lack CQ_use
df = df[!is.na(df$Prev_DR),]

# filter out duplicates
df = df[!duplicated(df),]

duplicated(df)
which(duplicated(df))
which(duplicated(df))->hi
df[hi,]
df = df[!duplicated(df),]

# count the number of years and studies per country
df2 = df %>% group_by(country) %>% 
  summarise(unique_year = length(unique(year)),
            total_rows = n())

### adding 11/1/23
# to make it I need at least two years and least two studies
madeit = df2 %>% filter(unique_year > 1, total_rows > 1)

# by country (list of dataframes)
by_country = list()
for(i in 1:nrow(madeit)){
  by_country[[i]] = df[df$country == madeit$country[i],]
}

# by regions (just know that this includes countries with only 1 study)
regions = df %>% select(Region) %>% unique() %>% unlist()
by_region = list()
for(i in regions){
  by_region[[i]] = df[df$Region == i,]
}

## try this lag-model
library(dlnm)
library(splines)

dothisone = by_country %>% bind_rows() %>% arrange(country,year)%>% group_by(country,year) # gathers all the data 6:48
dothisone = by_region %>% bind_rows() %>% arrange(country,year)%>% group_by(country,year) # gathers all the data 6:48

#unique(dothisone$Region)
#unique(df$Region)

write_csv(dothisone, file = "pfK13_lagdata.csv")
write_csv(dothisone, file = "pfK13_lagdata.0.csv")

###
###
?logit
library(boot)
df2 <- df %>% filter(PR<100,PR>0,Prev_DR<100,Prev_DR>0, CQ_use<100,CQ_use>0)
model3 <- lmer(logit(Prev_DR/100) ~ logit(PR/100) + logit(CQ_use/100) + (1|`study Id`), weights=tested, data=df2)
summary(model3)
plot(model3)

####
####
> View(crt_mdr)
> table(pfcrt_clean$`study Id`)->x
> table(pfcrt_clean$`study Id`) %>% as_tibble()->x

> table(pfcrt_clean$`study Id`) %>% as.data.frame()->x
> View(x)
> table(pfmdr1_clean$`study Id`) %>% as.data.frame()->y
> View(y)

crt_mdr <- inner_join(pfcrt_clean, pfmdr1_clean, by = "study Id")
crt_mdr <- right_join(pfcrt_clean, pfmdr1_clean, by = c("study Id","year","country"))

# grabbing unique ids across markers
ids = c(pfmdr1_clean$`study Id`,
        pfcrt_clean$`study Id`, 
        pfdhfr_clean$`study id`,
        pfdhps_clean$`study id`)
length(ids)
ids %>% unique() %>% length()

crt_id <- pfcrt_clean$`study Id`
length(crt_id )
crt_id %>% unique() %>% length()

mdr1_id <- pfmdr1_clean$`study Id`
length(mdr1_id )
mdr1_id %>% unique() %>% length()
