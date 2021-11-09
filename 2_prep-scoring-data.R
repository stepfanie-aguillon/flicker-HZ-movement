# This script contains the R code to prepare data for nls clines.
#
# Authors: Aguillon SM, VG Rohwer
# Year: 2021
# Title: Revisiting a classic hybrid zone: rapid movement of the
#        northern flicker hybrid zone in contemporary times
# Journal Info: TBD
# bioRxiv DOI: 10.1101/2021.08.16.456504
#
# Edited date: Oct 2021
#
# Please cite the paper if you use these scripts
#


# load packages
library(tidyverse)
library(geosphere)

# read datafile
df <- read_tsv("./raw-data/flicker-scoring-HZonly.txt",col_names=TRUE)

# group dataset by site locality and historic/contemporary, summarize various aspects of each locale
sub_summary <- df %>%
  group_by(HZtransect,site_ID) %>%
  summarize(samples=n(),
            shaft_mean = mean(na.omit(shaft)),
            shaft_se = sd(na.omit(shaft))/samples,
            nuchal_mean = mean(na.omit(nuchal)),
            nuchal_se = sd(na.omit(nuchal))/samples,
            crown_mean = mean(na.omit(crown)),
            crown_se = sd(na.omit(crown))/samples,
            ears_mean = mean(na.omit(ear_coverts)),
            ears_se = sd(na.omit(ear_coverts))/samples,
            throat_mean = mean(na.omit(throat)),
            throat_se = sd(na.omit(throat))/samples,
            malar_mean = mean(na.omit(as.integer(malar))),
            malar_se = sd(na.omit(as.integer(malar)))/samples,
            standHI_mean = mean(na.omit(standHI)),
            standHI_se = sd(na.omit(standHI))/samples) %>%
  arrange(site_ID)
#View(sub_summary)


# getting info on localities for cline modelling

# group dataset by locality only (combined historic+contemporary), to get average lat/long for each locale
geog_summary <- df %>%
  group_by(site_ID) %>%
  summarize(samples=n(),latitude=mean(lat),longitude=mean(long)) %>%
  arrange(site_ID)

head(geog_summary)
#site_ID samples latitude longitude
#<dbl>   <int>    <dbl>     <dbl>
#  1       1       5     40.9     -106.
#2       2      15     40.3     -105.
#3       3      13     40.4     -104.
#4       4      10     40.4     -104.
#5       5      11     40.3     -104.
#6       6       4     42.0     -104.

mean(geog_summary$latitude) #41.03611
min(geog_summary$latitude) #40.26798
max(geog_summary$latitude) #41.95949

# distances between localities using the mean latitude to get straightline distance between different longitudes

# mean longitude values
mlat = mean(geog_summary$latitude)

# first locality
# set up in longitude, latitude
p1=c(geog_summary$longitude[1],mlat)

# for loop to work across all localities and calculate the distance
# the loop calculates the distance between longitude values from each site using a mean latitude value
# the loop then assigns this distance to the variable "distance"
distance <- NA
for(x in 1:length(geog_summary$site_ID)){
  distance[x] <- distm(p1,c(geog_summary$longitude[x],mlat))/1000
}

# takes the summary table geog_summary and adds the distance column
dist_table <- geog_summary %>%
  add_column(dist=distance)

# reads in a df with site_ID, site_name, and Y/N for presence of historic and contemporary samples at the site
locale_info <- read_tsv("./raw-data/locality_info.txt",col_names=TRUE)
final_table <- inner_join(dist_table,locale_info,by="site_ID") %>%
  dplyr::select(site_ID,site_name,samples:dist,historic,contemporary)

# writes the final dataset
write.table(final_table,"scoring_locality_distances.txt",quote=FALSE,sep="\t",row.names=FALSE)




# PREPARE DATASETS TO CREATE CLINES

# subset final_table dataset to only include locality ID, name, lat/long, and distance value
locale_df <- final_table %>%
  dplyr::select(site_ID,site_name,latitude,longitude,dist)


## CONTEMPORARY SAMPLES
# filter summary df from above to include only contemporary sites
contemporary_summary <- sub_summary %>%
  filter(HZtransect=="contemporary") %>%
  ungroup()
# merge the two datasets and reorder columns
contemporary_merge <- inner_join(contemporary_summary,locale_df,by="site_ID") %>%
  dplyr::select(site_ID,site_name,samples,latitude:dist,shaft_mean:standHI_se)
write.table(contemporary_merge,"scoring_contemporary_input.txt",quote=FALSE,sep="\t",row.names=FALSE)


## HISTORIC SAMPLES
# filter summary df from above to include only contemporary sites
historic_summary <- sub_summary %>%
  filter(HZtransect=="historic") %>%
  ungroup()
# merge the two datasets and reorder columns
historic_merge <- inner_join(historic_summary,locale_df,by="site_ID") %>%
  dplyr::select(site_ID,site_name,samples,latitude:dist,shaft_mean:standHI_se)
write.table(historic_merge,"scoring_historic_input.txt",quote=FALSE,sep="\t",row.names=FALSE)





### DATASETS FOR SEPARATE SEXES 

# MALES ONLY
sub_summary_males <- sub_df %>%
  filter(sex=="Male") %>%
  group_by(HZtransect,site_ID) %>%
  summarize(samples=n(),
            shaft_mean = mean(na.omit(shaft)),
            shaft_se = sd(na.omit(shaft))/samples,
            nuchal_mean = mean(na.omit(nuchal)),
            nuchal_se = sd(na.omit(nuchal))/samples, 
            crown_mean = mean(na.omit(crown)),
            crown_se = sd(na.omit(crown))/samples, 
            ears_mean = mean(na.omit(ear_coverts)),
            ears_se = sd(na.omit(ear_coverts))/samples,
            throat_mean = mean(na.omit(throat)),
            throat_se = sd(na.omit(throat))/samples,
            malar_mean = mean(na.omit(as.integer(malar))),
            malar_se = sd(na.omit(as.integer(malar)))/samples,
            standHI_mean = mean(na.omit(standHI)),
            standHI_se = sd(na.omit(standHI))/samples) %>%
  arrange(site_ID)

# contemporary samples
contemporary_males_summary <- sub_summary_males %>%
  filter(HZtransect=="contemporary") %>%
  ungroup()
# merge the two datasets and reorder columns
contemporary_males_merge <- inner_join(contemporary_males_summary,locale_df,by="site_ID") %>%
  dplyr::select(site_ID,site_name,samples,latitude:dist,shaft_mean:standHI_se)
write.table(contemporary_males_merge,"scoring_contemporary_input_males.txt",quote=FALSE,sep="\t",row.names=FALSE)

# historic samples
historic_males_summary <- sub_summary_males %>%
  filter(HZtransect=="historic") %>%
  ungroup()
# merge the two datasets and reorder columns
historic_males_merge <- inner_join(historic_males_summary,locale_df,by="site_ID") %>%
  dplyr::select(site_ID,site_name,samples,latitude:dist,shaft_mean:standHI_se)
write.table(historic_males_merge,"scoring_historic_input_males.txt",quote=FALSE,sep="\t",row.names=FALSE)




# FEMALES ONLY
sub_summary_females <- sub_df %>%
  filter(sex=="Female") %>%
  group_by(HZtransect,site_ID) %>%
  summarize(samples=n(),
            shaft_mean = mean(na.omit(shaft)),
            shaft_se = sd(na.omit(shaft))/samples,
            nuchal_mean = mean(na.omit(nuchal)),
            nuchal_se = sd(na.omit(nuchal))/samples, 
            crown_mean = mean(na.omit(crown)),
            crown_se = sd(na.omit(crown))/samples, 
            ears_mean = mean(na.omit(ear_coverts)),
            ears_se = sd(na.omit(ear_coverts))/samples,
            throat_mean = mean(na.omit(throat)),
            throat_se = sd(na.omit(throat))/samples,
            #malar_mean = mean(na.omit(as.integer(malar))),
            #malar_se = sd(na.omit(as.integer(malar)))/samples,
            standHI_mean = mean(na.omit(standHI)),
            standHI_se = sd(na.omit(standHI))/samples) %>%
  arrange(site_ID)

# contemporary samples
contemporary_females_summary <- sub_summary_females %>%
  filter(HZtransect=="contemporary") %>%
  ungroup()
# merge the two datasets and reorder columns
contemporary_females_merge <- inner_join(contemporary_females_summary,locale_df,by="site_ID") %>%
  dplyr::select(site_ID,site_name,samples,latitude:dist,shaft_mean:standHI_se)
write.table(contemporary_females_merge,"scoring_contemporary_input_females.txt",quote=FALSE,sep="\t",row.names=FALSE)

# historic samples
historic_females_summary <- sub_summary_females %>%
  filter(HZtransect=="historic") %>%
  ungroup()
# merge the two datasets and reorder columns
historic_females_merge <- inner_join(historic_females_summary,locale_df,by="site_ID") %>%
  dplyr::select(site_ID,site_name,samples,latitude:dist,shaft_mean:standHI_se)
write.table(historic_females_merge,"scoring_historic_input_females.txt",quote=FALSE,sep="\t",row.names=FALSE)
