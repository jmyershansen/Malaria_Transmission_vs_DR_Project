library(tidyverse)

# load data
WWARN <- read_csv("WWARN_ACT_Partner_Drug_Mol_Surveyor_Data(6).csv")
WWARN1 <- read_csv("WWARN_pfdhfr_pfdhps_data.csv")
WWARN2 <- read_csv("K13_surveyor_data.csv")

df <- rbind(WWARN %>% select(country, lat, lon), WWARN1 %>% select(country, lat, lon), 
           WWARN2 %>% select(country, lat, lon))
df = df %>% distinct()
LATLONG <- cbind( df %>% select(lat, lon))
#LATLONG <- as.data.frame(cbind(df$lat,df$long))
#colnames(LATLONG)<-c("latitude","longitude")
write_csv(LATLONG,file="LatitudeLongitude.csv")

EXTR <- read.csv("Extracted-points-data.csv")
####
####
####


WWARN_full <- read_csv("WWARN_ACT_Partner_Drug_Mol_Surveyor_Data(6).csv")

WWARN_full <- cbind(WWARN_full$`study end`,WWARN_full)

colnames(WWARN_full)[1] <- "year"

colnames(WWARN_full)[6] <- "long"

WWARN_full <- WWARN_full %>% arrange(year) %>% group_by(lat,long,year)

WWARN_full$long <- round(WWARN_full$long, digits = 2)
EXTR$long <- round (EXTR$long, digits =2)

WWARN_full$lat <- round(WWARN_full$lat, digits = 2)
EXTR$lat <- round(EXTR$lat, digits = 2)

Complete <- full_join(EXTR,WWARN_full,by = c("lat","long","year"))

Complete <- Complete %>% mutate(PR = as.numeric(value)*100)
View(Complete)

Complete_pfcrt <- Complete %>% filter(`marker group`=="pfcrt 76T")

Complete_pfmdr1 <- Complete %>% filter(`marker group`=="pfmdr1 86Y")

#write_csv(Complete, file = "pfcrt_pfmdr1_pm.csv")
#write_csv(Complete_pfcrt, file = "pfcrt_pm.csv")
#write_csv(Complete_pfmdr1, file = "pfmdr1_pm.csv")

######
Old_data_pfcrt <- read_csv("WWARN_pfcrt.csv")

Country_region <- Old_data_pfcrt[c(3,5)] %>% distinct(country, Region) %>% arrange(country)

Complete_pfcrt_region <- right_join(Country_region, Complete_pfcrt, by = "country")
pfcrt_SNP<-Complete_pfcrt_region

write_csv(Complete_pfcrt_region, file = "pfcrt_SNP.csv")
######

Complete_pfmdr1_region <- right_join(Country_region, Complete_pfmdr1, by = "country")
pfmdr1_SNP<-Complete_pfmdr1_region 

write_csv(Complete_pfmdr1_region, file = "pfmdr1_SNP.csv")
######
#####
#####
#pfcrt_SNP <- read_csv ("pfcrt_SNP.csv")
pfcrt_SNP <- pfcrt_SNP %>% mutate(Prev_DR=(present/tested)*100)
pfmdr1_SNP <- pfmdr1_SNP %>% mutate(Prev_DR=(present/tested)*100)



write_csv(DU, file = "CQ_imputed.csv")
cq<- read_csv("CQ_imputed.csv")
cq2 <- cq$country %>% unique()
cq2 <- length(unique(cq$country))  

crt<-read_csv ("pfcrt_SNP.csv")
crt<-pfcrt_SNP
crt2 <- crt$country %>% unique()
crt2 <- length(unique(crt$country))
crt2 <- length(unique(crt$`study Id`))
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

mdr1<-read_csv ("pfmdr1_SNP.csv")
crt<-pfcrt_SNP
crt2 <- crt$country %>% unique()
crt2 <- length(unique(crt$country))
mdr12 <- length(unique(mdr1$`study Id`))

dhps<-read_csv ("pfdhps_SNP.csv")
crt<-pfcrt_SNP
crt2 <- crt$country %>% unique()
crt2 <- length(unique(crt$country))
dhps2 <- length(unique(dhps$`study id`))


dhfr<-read_csv ("pfdhfr_SNP.csv")
crt<-pfcrt_SNP
crt2 <- crt$country %>% unique()
crt2 <- length(unique(crt$country))
dhfr2 <- length(unique(dhfr$`study id`))

k13<-read_csv ("pfK13_SNP.csv")
crt<-pfcrt_SNP
crt2 <- crt$country %>% unique()
crt2 <- length(unique(crt$country))
k132 <- length(unique(k13$sid))


###Maps
###
###
install.packages("maps")
install.packages("rnaturalearth")
install.packages("rnaturalearthdata")
install.packages("sf")
library(maps)
library(tidyverse)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggplot2)
library(sf)

###Using map & ggplot2
dataset = read_csv("pfcrt_SNP.csv")
dataset <-crt_na_rem
# Convert map data to a data frame for ggplot2
world_map <- map_data("world")
# Create a ggplot2 map
ggplot() +
  geom_polygon(data = world_map, aes(x = long, y = lat, group = group), 
               fill = "lightgray", color = "white") +
  geom_point(data = crt_na_rem, aes(x = long, y = lat, color = Prev_DR), 
             size = 3, alpha = 0.7) +
  scale_color_gradient(low = "yellow", high = "red") +
  theme_minimal() +
  labs(title = "Prevalence of Pfcrt_76T ", x = "Longitude", y = "Latitude", 
       color = "Pfcrt_76T")

###Using sf & ggplot2
#data <- read.csv("pfcrt_SNP.csv")
data <-crt_na_rem
# Convert the data into a simple features (sf) object
data_sf <- st_as_sf(data, coords = c("long", "lat"), crs = 4326)
# Get the world map data
# Load world map as an sf object
world <- ne_countries(scale = "medium", returnclass = "sf") 
# Plot with ggplot2
ggplot(data = world) +
  geom_sf(fill = "gray90", color = "gray80") +               
  geom_sf(data = data_sf, aes(color = Prev_DR),           
          size = 2, alpha = 0.7) +                           
  scale_color_gradient(low = "yellow", high = "red") +       
  theme_minimal() +
  labs(title = "Malaria Prevalence Study Sites",
       color = "Prevalence") +
  theme(plot.title = element_text(hjust = 0.5)) 
###
###
