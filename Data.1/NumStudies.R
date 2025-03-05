#########
WWARN_full <- read_csv("WWARN_ACT_Partner_Drug_Mol_Surveyor_Data(6).csv")
EXTR <- read.csv("Extracted-points-data.csv")

WWARN_full <- cbind(WWARN_full$`study end`,WWARN_full)

colnames(WWARN_full)[1] <- "year"

colnames(WWARN_full)[6] <- "long"

WWARN_full <- WWARN_full %>% arrange(year) %>% group_by(lat,long,year)

WWARN_full$long <- round(WWARN_full$long, digits = 2)
EXTR$long <- round (EXTR$long, digits =2)

WWARN_full$lat <- round(WWARN_full$lat, digits = 2)
EXTR$lat <- round(EXTR$lat, digits = 2)

Complete <- right_join(EXTR,WWARN_full,by = c("lat","long","year"))
Complete <- left_join(EXTR,WWARN_full,by = c("lat","long","year"))


Complete <- Complete %>% mutate(PR = as.numeric(value)*100)
View(Complete)
Complete_pfcrt <- Complete %>% filter(`marker group`=="pfcrt 76T")
Complete_pfmdr1 <- Complete %>% filter(`marker group`=="pfmdr1 86Y")

Complete_pfcrt <- Complete_pfcrt %>% mutate(Prev_DR=(present/tested)*100)
Numstudies1 <- length(unique(Complete_pfcrt$`study Id`))
Numstudies2 <- length((Complete_pfcrt$`Prev_DR`))

Complete_pfmdr1 <- Complete_pfmdr1 %>% mutate(Prev_DR=(present/tested)*100)
Numstudies1 <- length(unique(Complete_pfmdr1$`study Id`))
Numstudies2 <- length((Complete_pfmdr1$`Prev_DR`))
 
######
Old_data_pfcrt <- read_csv("WWARN_pfcrt.csv")
Country_region <- Old_data_pfcrt[c(3,5)] %>% distinct(country, Region) %>% arrange(country)
Complete_pfcrt_region <- right_join(Country_region, Complete_pfcrt, by = "country")



DHPS.DHFR <- read_csv("WWARN_pfdhfr_pfdhps_data.csv")
EXTR1 <- read.csv("Extracted-points-data.1.csv")

DHPS.DHFR <- cbind(DHPS.DHFR$`study end year`,DHPS.DHFR)

colnames(DHPS.DHFR)[1] <- "year"

colnames(DHPS.DHFR)[7] <- "long"

DHPS.DHFR <- DHPS.DHFR %>% arrange(year) %>% group_by(lat,long,year)

DHPS.DHFR$long <- round(DHPS.DHFR$long, digits = 2)
EXTR1$long <- round (EXTR1$long, digits =2)

DHPS.DHFR$lat <- round(DHPS.DHFR$lat, digits = 2)
EXTR1$lat <- round(EXTR1$lat, digits = 2)

Complete3 <- right_join(EXTR1,DHPS.DHFR,by = c("lat","long","year"))
Complete3 <- left_join(EXTR1,DHPS.DHFR,by = c("lat","long","year"))


Complete3 <- Complete3 %>% mutate(PR = as.numeric(value)*100)
#Numstudies7 <- length(unique(Complete3$`study id`))
Complete_pfdhps <- Complete3 %>% filter(`mutation`=="dhps 437G")
Complete_pfdhfr <- Complete3 %>% filter(`mutation`=="dhfr 108N")

Complete_pfdhps <- Complete_pfdhps %>% mutate(Prev_DR=(present/tested)*100)
Numstudies1 <- length(unique(Complete_pfdhps$`study id`))
Numstudies2 <- length((Complete_pfdhps$`Prev_DR`))

Complete_pfdhfr <- Complete_pfdhfr %>% mutate(Prev_DR=(present/tested)*100)
Numstudies1 <- length(unique(Complete_pfdhfr$`study id`))
Numstudies2 <- length((Complete_pfdhfr$`Prev_DR`))


# PFK13 POINT MUTATIONS
K13 <- read_csv("K13_surveyor_data.csv")
EXTR2 <- read.csv("Extracted-points-data.2.csv")

colnames(K13)[7] <- "long"

K13 <- K13 %>% arrange(year) %>% group_by(lat,long,year)

K13$long <- round(K13$long, digits = 2)
EXTR2$long <- round (EXTR2$long, digits =2)

K13$lat <- round(K13$lat, digits = 2)
EXTR2$lat <- round(EXTR2$lat, digits = 2)

Complete4 <- right_join(EXTR2,K13,by = c("lat","long","year"))
Complete4 <- left_join(EXTR2,K13,by = c("lat","long","year"))


Complete4 <- Complete4 %>% mutate(PR = as.numeric(value)*100)

Complete_pfK13 <- Complete4 %>% filter(`mutation`=="C580Y")
#Numstudies8 <- length(unique(Complete_pfK13$`sid`))

Complete_pfK13 <- Complete_pfK13 %>% mutate(Prev_DR=(present/tested)*100)
Numstudies1 <- length(unique(Complete_pfK13$`sid`))
Numstudies2 <- length((Complete_pfK13$`Prev_DR`))

