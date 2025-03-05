df <- read_csv("WWARN_ACT_Partner_Drug_Mol_Surveyor_Data(6).csv")
df1 <- df %>% filter(`marker group`=="pfcrt 76T")
df2 <- df %>% filter(`marker group`=="pfmdr1 86Y")
Numstudies <- length(unique(df$`study Id`))
Numstudies1 <- length(unique(df1$`study Id`))
Numstudies2 <- length(unique(df2$`study Id`))

df1 <- df1 %>% mutate(year=`study end`)
df2 <- df2 %>% mutate(year=`study end`)

duplicated(df1)
which(duplicated(df1))
which(duplicated(df1))->hi
df1[hi,]

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
df3 <- inner_join(EXTR,WWARN_full,by = c("lat","long","year"))
#Complete <- Complete %>% mutate(PR = as.numeric(value)*100)

Complete_pfcrt <- df3 %>% filter(`marker group`=="pfcrt 76T")
Complete_pfmdr1 <- df3 %>% filter(`marker group`=="pfmdr1 86Y")
Numstudies3 <- length(unique(Complete_pfcrt$`study Id`))
Numstudies4 <- length(unique(Complete_pfmdr1$`study Id`))

pfcrt_SNP <- read_csv ("pfcrt_SNP.csv")
pfmdr1_SNP <- read_csv ("pfmdr1_SNP.csv")
DUdata <- read_csv ("FeveredChildrenDrugTreatProp_CQ.0.csv")
#df3 <- length(unique(pfcrt_SNP$`study Id`))

pfcrt_complete <- inner_join(DUdata,pfcrt_SNP, by = c("country","year"))
pfcrt_complete0 <- inner_join(DUdata,df1, by = c("country","year"))
df4 <- length(unique(pfcrt_complete$`study Id`))
df4 <- length(unique(pfcrt_complete0$`study Id`)) 

pfcrt_complete1 <- inner_join(DUdata,Complete_pfcrt, by = c("country","year"))
df5 <- length(unique(pfcrt_complete1$`study Id`))

pfmdr_complete <- inner_join(DUdata,pfmdr1_SNP, by = c("country","year"))
pfmdr_complete0 <- inner_join(DUdata,df2, by = c("country","year"))
df6 <- length(unique(pfmdr_complete$`study Id`))
df6 <- length(unique(pfmdr_complete0$`study Id`))

pfmdr_complete1 <- inner_join(DUdata,Complete_pfmdr1, by = c("country","year"))
df7 <- length(unique(pfmdr_complete1$`study Id`))

###
###
df <- read_csv("WWARN_pfdhfr_pfdhps_data.csv")
df1 <- df %>% filter(`mutation`=="dhps 437G")
df2 <- df %>% filter(`mutation`=="dhfr 108N")
Numstudies <- length(unique(df$`study id`))
Numstudies1 <- length(unique(df1$`study id`))
Numstudies2 <- length(unique(df2$`study id`))

df1 <- df1 %>% mutate(year=`study end year`)
df2 <- df2 %>% mutate(year=`study end year`)

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
Complete3 <- inner_join(EXTR1,DHPS.DHFR,by = c("lat","long","year"))

Complete_pfdhps <- Complete3 %>% filter(`mutation`=="dhps 437G")
Complete_pfdhfr <- Complete3 %>% filter(`mutation`=="dhfr 108N")
Numstudies3 <- length(unique(Complete_pfdhps$`study id`))
Numstudies4 <- length(unique(Complete_pfdhfr$`study id`))

pfdhps_SNP <- read_csv ("pfdhps_SNP.csv")
pfdhfr_SNP <- read_csv ("pfdhfr_SNP.csv")
DUdata <- read_csv ("FeveredChildrenDrugTreatProp_SP.0.csv")
#df3 <- length(unique(pfdhps_SNP$`study id`))
#df3 <- length(unique(pfdhfr_SNP$`study id`))

pfdhps_complete <- inner_join(DUdata,pfdhps_SNP, by = c("country","year"))
pfdhps_complete0 <- inner_join(DUdata,df1, by = c("country","year"))
df4 <- length(unique(pfdhps_complete$`study id`))
df4 <- length(unique(pfdhps_complete0$`study id`))

pfdhps_complete1 <- inner_join(DUdata,Complete_pfdhps, by = c("country","year"))
df5 <- length(unique(pfdhps_complete1$`study id`))

pfdhfr_complete <- inner_join(DUdata,pfdhfr_SNP, by = c("country","year"))
pfdhfr_complete0 <- inner_join(DUdata,df2, by = c("country","year"))
df6 <- length(unique(pfdhfr_complete$`study id`))
df6 <- length(unique(pfdhfr_complete0$`study id`))

pfdhfr_complete1 <- inner_join(DUdata,Complete_pfdhfr, by = c("country","year"))
df7 <- length(unique(pfdhfr_complete1$`study id`))

###
###
df <- read_csv("K13_surveyor_data.csv")
df1 <- df %>% filter(`mutation`=="C580Y")
Numstudies <- length(unique(df$`sid`))
Numstudies1 <- length(unique(df1$`sid`))

K13 <- read_csv("K13_surveyor_data.csv")
EXTR2 <- read.csv("Extracted-points-data.2.csv")
colnames(K13)[7] <- "long"
K13 <- K13 %>% arrange(year) %>% group_by(lat,long,year)

K13$long <- round(K13$long, digits = 2)
EXTR2$long <- round (EXTR2$long, digits =2)
K13$lat <- round(K13$lat, digits = 2)
EXTR2$lat <- round(EXTR2$lat, digits = 2)

Complete4 <- inner_join(EXTR2,K13,by = c("lat","long","year"))
Complete_pfK13 <- Complete4 %>% filter(`mutation`=="C580Y")
Numstudies3 <- length(unique(Complete_pfK13$`sid`))

pfK13_SNP <- read_csv ("pfK13_SNP.csv")
DUdata <- read_csv ("ACT_use.0.csv")

pfK13_complete <- inner_join(DUdata,pfK13_SNP, by = c("country","year"))
pfK13_complete0 <- inner_join(DUdata,df1, by = c("country","year"))
df4 <- length(unique(pfK13_complete$`sid`))
df4 <- length(unique(pfK13_complete0$`sid`))

pfK13_complete1 <- inner_join(DUdata,Complete_pfK13, by = c("country","year"))
df5 <- length(unique(pfK13_complete1$`sid`))
