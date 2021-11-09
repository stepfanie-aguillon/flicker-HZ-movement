# This script contains the R code to map sampling localities.
#
# Authors: Aguillon SM, VG Rohwer
# Year: 2021
# Title: Revisiting a classic hybrid zone: rapid movement of the
#        northern flicker hybrid zone in contemporary times
# Journal Info: TBD
# bioRxiv DOI: 10.1101/2021.08.16.456504
#
# Edited date: June 2021
#
# Please cite the paper if you use these scripts
#


# load packages
library(tidyverse)
library(rgdal)
library(proj4)
library(maptools)
library(ggrepel)

# load df with lat/long points
localities <- read_tsv("./scoring_locality_distances.txt")
# which localities are sampled in both transects
localities <- localities %>%
  mutate(both_transects = ifelse(historic==contemporary,"both","one"))


# get state boundaries
CO_NE <- map_data("state") %>%
  filter(region=="colorado" | region=="nebraska")


# load river data
url.river_data <- url("http://sharpsightlabs.com/wp-content/datasets/usa_rivers.RData")
load(url.river_data)

# remove misc features
lines.rivers <- subset(lines.rivers, !(FEATURE %in% c("Shoreline"  ,"Shoreline Intermittent","Null","Closure Line","Apparent Limit")))

# subset to only Nebraska and Colorado
lines.rivers <- subset(lines.rivers, (STATE %in% c("NE","CO")))


# get Platte River
platteData = lines.rivers[grep("Platte",lines.rivers$NAME),]
df.platte_rivers <- fortify(platteData)
# inside the desired area
df.platte_rivers <- fortify(platteData)
platte_clean <- df.platte_rivers %>%
  filter(long>-104.9999) %>%
  filter(lat>40.1)

# get eastern rivers near Omaha, NE (also Platte but not picked up)
rivers_east <- fortify(lines.rivers)
rivers_east <- rivers_east %>%
  filter(id %in% c(48302, 48469, 47172, 47244, 47196))



# plot figure
ggplot() +
  geom_polygon(data=CO_NE,aes(x=long,y=lat,group=group),fill="gray",color="dark gray",alpha=0.3) +
  geom_path(data=platte_clean,aes(x=long,y=lat,group=group),color="deepskyblue2",size=0.5) +
  geom_path(data=rivers_east,aes(x=long,y=lat,group=group),color="deepskyblue2",size=0.5) +
  geom_point(data=localities,aes(x=longitude,y=latitude,size=both_transects)) +
  geom_text_repel(data=localities,aes(x=longitude,y=latitude,label=site_ID)) +
  scale_size_manual(values=c(2,1)) +
  theme_bw() +
  xlab("Longitude") +
  ylab("Latitude") +
  coord_map() +
  theme(legend.position="none",axis.title=element_text(face="bold",size=12), axis.text=element_text(size=10,color="black"))

# save plot as PDF with dimensions 6" x 4"
