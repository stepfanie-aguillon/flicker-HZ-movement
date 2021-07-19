# This script contains the R code to analyze photo data of contemporary samples.
#
# Authors: Aguillon SM, VG Rohwer
# Year: 2021
# Title: Revisiting a classic hybrid zone: rapid movement of the
#        northern flicker hybrid zone in contemporary times
# Journal Info: TBD
# bioRxiv DOI: TBD
#
# Edited date: June 2021
#
# Please cite the paper if you use these scripts
#


# load packages
library(tidyverse)
library(nlstools)
library(gridExtra)

# load dataset
photo_df <- read_tsv("./raw-data/photo_scores.txt",col_names=TRUE)


# STEPS FOR EACH PLUMAGE TRAIT
# 1. lm comparing distance and each photo channel
# 2. for significant channels, run clines using nls
# 3. run multivariate lms to compare categorical score with photo channels




########### STEP 1 ###########

# colorblind friendly plot colors
cbf_red <- "#DC3220"
cbf_green <- "#6EE26D"
cbf_blue <- "#0A58A6"
cbf_purple <- "#B43BEE"


### NUCHAL PATCH
nuchal_df <- photo_df %>%
  filter(trait=="nuchal") %>%
  distinct()

summary(lm(nuchal_df$area ~ nuchal_df$dist)) #significant

# for nuchal, area is significant

nuchal_plot <- ggplot(data=nuchal_df,aes(x=dist)) +
  geom_point(aes(y=area),color="black") +
  geom_smooth(aes(y=area),color="black",method=lm,se=FALSE) +
  ylab("Photo value") +
  xlab("Distance (km)") +
  ggtitle("Nuchal patch") +
  theme_bw() +
  theme(axis.title=element_text(face="bold",size=12), axis.text=element_text(size=10,color="black"))



### SHAFT
shaft_df <- photo_df %>%
  filter(trait=="vane") %>%
  distinct()

summary(lm(shaft_df$lwMean ~ shaft_df$dist)) #significant
summary(lm(shaft_df$mwMean ~ shaft_df$dist)) #significant
summary(lm(shaft_df$swMean ~ shaft_df$dist)) #NS
summary(lm(shaft_df$uvMean ~ shaft_df$dist)) #NS
summary(lm(shaft_df$lumMean ~ shaft_df$dist)) #significant

# for shaft, lwMean, mwMean, and lumMean are significant

shaft_plot <- ggplot(data=shaft_df,aes(x=dist)) +
  geom_point(aes(y=lwMean),color=cbf_red) +
  geom_smooth(aes(y=lwMean),color=cbf_red,method=lm,se=FALSE) +
  geom_point(aes(y=mwMean),color=cbf_green) +
  geom_smooth(aes(y=mwMean),color=cbf_green,method=lm,se=FALSE) +
  geom_point(aes(y=swMean),color=cbf_blue) +
  geom_smooth(aes(y=swMean),color=cbf_blue,method=lm,se=FALSE) +
  geom_point(aes(y=uvMean),color=cbf_purple) +
  geom_smooth(aes(y=uvMean),color=cbf_purple,method=lm,se=FALSE) +
  geom_point(aes(y=lumMean),color="gray55") +
  geom_smooth(aes(y=lumMean),color="gray55",method=lm,se=FALSE) +
  ylab("Photo value") +
  xlab("Distance (km)") +
  ggtitle("Shaft") +
  theme_bw() +
  theme(axis.title=element_text(face="bold",size=12), axis.text=element_text(size=10,color="black"))



### THROAT
throat_df <- photo_df %>%
  filter(trait=="throat") %>%
  distinct()

summary(lm(throat_df$lwMean ~ throat_df$dist)) #significant
summary(lm(throat_df$mwMean ~ throat_df$dist)) #marginally significant
summary(lm(throat_df$swMean ~ throat_df$dist)) #significant
summary(lm(throat_df$uvMean ~ throat_df$dist)) #marginally significant
summary(lm(throat_df$lumMean ~ throat_df$dist)) #significant

# for throat, lwMean, swMean, and lumMean are significant

throat_plot <- ggplot(data=throat_df,aes(x=dist)) +
  geom_point(aes(y=lwMean),color=cbf_red) +
  geom_smooth(aes(y=lwMean),color=cbf_red,method=lm,se=FALSE) +
  geom_point(aes(y=mwMean),color=cbf_green) +
  geom_smooth(aes(y=mwMean),color=cbf_green,method=lm,se=FALSE) +
  geom_point(aes(y=swMean),color=cbf_blue) +
  geom_smooth(aes(y=swMean),color=cbf_blue,method=lm,se=FALSE) +
  geom_point(aes(y=uvMean),color=cbf_purple) +
  geom_smooth(aes(y=uvMean),color=cbf_purple,method=lm,se=FALSE) +
  geom_point(aes(y=lumMean),color="gray55") +
  geom_smooth(aes(y=lumMean),color="gray55",method=lm,se=FALSE) +
  ylab("Photo value") +
  xlab("Distance (km)") +
  ggtitle("Throat") +
  theme_bw() +
  theme(axis.title=element_text(face="bold",size=12), axis.text=element_text(size=10,color="black"))





### CROWN
crown_df <- photo_df %>%
  filter(trait=="crown") %>%
  distinct()

summary(lm(crown_df$lwMean ~ crown_df$dist)) #NS
summary(lm(crown_df$mwMean ~ crown_df$dist)) #significant
summary(lm(crown_df$swMean ~ crown_df$dist)) #significant
summary(lm(crown_df$uvMean ~ crown_df$dist)) #significant
summary(lm(crown_df$lumMean ~ crown_df$dist)) #significant

# for crown, mwMean, swMean, uvMean, and lumMean are significant

crown_plot <- ggplot(data=crown_df,aes(x=dist)) +
  geom_point(aes(y=lwMean),color=cbf_red) +
  geom_smooth(aes(y=lwMean),color=cbf_red,method=lm,se=FALSE) +
  geom_point(aes(y=mwMean),color=cbf_green) +
  geom_smooth(aes(y=mwMean),color=cbf_green,method=lm,se=FALSE) +
  geom_point(aes(y=swMean),color=cbf_blue) +
  geom_smooth(aes(y=swMean),color=cbf_blue,method=lm,se=FALSE) +
  geom_point(aes(y=uvMean),color=cbf_purple) +
  geom_smooth(aes(y=uvMean),color=cbf_purple,method=lm,se=FALSE) +
  geom_point(aes(y=lumMean),color="gray55") +
  geom_smooth(aes(y=lumMean),color="gray55",method=lm,se=FALSE) +
  ylab("Photo value") +
  xlab("Distance (km)") +
  ggtitle("Crown") +
  theme_bw() +
  theme(axis.title=element_text(face="bold",size=12), axis.text=element_text(size=10,color="black"))





### EAR COVERTS
ear_df <- photo_df %>%
  filter(trait=="ear") %>%
  distinct()

summary(lm(ear_df$lwMean ~ ear_df$dist)) #significant
summary(lm(ear_df$mwMean ~ ear_df$dist)) #NS
summary(lm(ear_df$swMean ~ ear_df$dist)) #significant
summary(lm(ear_df$uvMean ~ ear_df$dist)) #NS
summary(lm(ear_df$lumMean ~ ear_df$dist)) #NS

# for ear coverts, lwMean and swMean are significant

ear_plot <- ggplot(data=ear_df,aes(x=dist)) +
  geom_point(aes(y=lwMean),color=cbf_red) +
  geom_smooth(aes(y=lwMean),color=cbf_red,method=lm,se=FALSE) +
  geom_point(aes(y=mwMean),color=cbf_green) +
  geom_smooth(aes(y=mwMean),color=cbf_green,method=lm,se=FALSE) +
  geom_point(aes(y=swMean),color=cbf_blue) +
  geom_smooth(aes(y=swMean),color=cbf_blue,method=lm,se=FALSE) +
  geom_point(aes(y=uvMean),color=cbf_purple) +
  geom_smooth(aes(y=uvMean),color=cbf_purple,method=lm,se=FALSE) +
  geom_point(aes(y=lumMean),color="gray55") +
  geom_smooth(aes(y=lumMean),color="gray55",method=lm,se=FALSE) +
  ylab("Photo value") +
  xlab("Distance (km)") +
  ggtitle("Ear coverts") +
  theme_bw() +
  theme(axis.title=element_text(face="bold",size=12), axis.text=element_text(size=10,color="black"))





### MALE MALAR
malar_df <- photo_df %>%
  filter(trait=="malar") %>%
  filter(sex=="Male") %>%
  distinct()

summary(lm(malar_df$lwMean ~ malar_df$dist)) #significant
summary(lm(malar_df$mwMean ~ malar_df$dist)) #significant
summary(lm(malar_df$swMean ~ malar_df$dist)) #significant
summary(lm(malar_df$uvMean ~ malar_df$dist)) #marginally significant
summary(lm(malar_df$lumMean ~ malar_df$dist)) #significant

# for malar, lwMean, mwMean, swMean, and lumMean are significant

malar_plot <- ggplot(data=malar_df,aes(x=dist)) +
  geom_point(aes(y=lwMean),color=cbf_red) +
  geom_smooth(aes(y=lwMean),color=cbf_red,method=lm,se=FALSE) +
  geom_point(aes(y=mwMean),color=cbf_green) +
  geom_smooth(aes(y=mwMean),color=cbf_green,method=lm,se=FALSE) +
  geom_point(aes(y=swMean),color=cbf_blue) +
  geom_smooth(aes(y=swMean),color=cbf_blue,method=lm,se=FALSE) +
  geom_point(aes(y=uvMean),color=cbf_purple) +
  geom_smooth(aes(y=uvMean),color=cbf_purple,method=lm,se=FALSE) +
  geom_point(aes(y=lumMean),color="gray55") +
  geom_smooth(aes(y=lumMean),color="gray55",method=lm,se=FALSE) +
  ylab("Photo value") +
  xlab("Distance (km)") +
  ggtitle("Male malar stripe") +
  theme_bw() +
  theme(axis.title=element_text(face="bold",size=12), axis.text=element_text(size=10,color="black"))


# plot lms for all six traits together
grid.arrange(crown_plot, ear_plot, malar_plot, nuchal_plot, shaft_plot, throat_plot, nrow=3)
# save plot as PDF with dimensions 6" x 7.75"







######### STEP 2 ###########

# obtain average and se for each site
# both for photo values and for categorical scoring
cline_df <- photo_df %>%
  group_by(site_ID,site_name,dist,trait) %>%
  summarize(samples=n(),
            # photo traits
            lw_mean = mean(na.omit(lwMean)),
            lw_se = sd(na.omit(lwMean))/samples,
            mw_mean = mean(na.omit(mwMean)),
            mw_se = sd(na.omit(mwMean))/samples,
            sw_mean = mean(na.omit(swMean)),
            sw_se = sd(na.omit(swMean))/samples,
            uv_mean = mean(na.omit(uvMean)),
            uv_se = sd(na.omit(uvMean))/samples,
            lum_mean = mean(na.omit(lumMean)),
            lum_se = sd(na.omit(lumMean))/samples,
            area_mean = mean(na.omit(area)),
            area_se = mean(na.omit(area))/samples) %>%
  arrange(site_ID)



### NUCHAL PATCH
# photo characters: area
nuchal_cline <- cline_df %>%
  filter(trait=="nuchal") %>%
  select(site_ID:samples,starts_with("area"))

nuchal_area_sc <- mean(filter(nuchal_cline,site_ID>25)$area_mean)

nuchal_cline_plot <- ggplot() +
  geom_point(data=nuchal_cline,aes(x=dist,y=area_mean),color="black") +
  geom_smooth(data=nuchal_cline,aes(x=dist,y=area_mean),method="nls",formula=y~nuchal_area_sc/(1+exp(-(4*(x-c))/w)),se=FALSE,method.args=list(start=list(c=490,w=490)),color="black") +
  xlab("Distance (km)") +
  ylab("Photo value") +
  ggtitle("Nuchal patch") +
  theme_bw() +
  theme(axis.title=element_text(face="bold",size=12), axis.text=element_text(size=10,color="black"))



### SHAFT
# photo characters: lwMean, mwMean, lumMean
shaft_cline <- cline_df %>%
  filter(trait=="vane")

shaft_lw_sc <- mean(filter(shaft_cline,site_ID>25)$lw_mean)
shaft_mw_sc <- mean(filter(shaft_cline,site_ID>25)$mw_mean)
#shaft_sw_sc <- mean(filter(shaft_cline,site_ID>25)$sw_mean)
#shaft_uv_sc <- mean(filter(shaft_cline,site_ID>25)$uv_mean)
shaft_lum_sc <- mean(filter(shaft_cline,site_ID>25)$lum_mean)

shaft_cline_plot <- ggplot() +
  geom_point(data=shaft_cline,aes(x=dist,y=lw_mean),color=cbf_red) +
  geom_smooth(data=shaft_cline,aes(x=dist,y=lw_mean),method="nls",formula=y~shaft_lw_sc/(1+exp(-(4*(x-c))/w)),se=FALSE,method.args=list(start=list(c=490,w=490)),color=cbf_red) +
  geom_point(data=shaft_cline,aes(x=dist,y=mw_mean),color=cbf_green) +
  geom_smooth(data=shaft_cline,aes(x=dist,y=mw_mean),method="nls",formula=y~shaft_mw_sc/(1+exp(-(4*(x-c))/w)),se=FALSE,method.args=list(start=list(c=490,w=490)),color=cbf_green) +
  #geom_point(data=shaft_cline,aes(x=dist,y=sw_mean),color=cbf_blue) +
  #geom_smooth(data=shaft_cline,aes(x=dist,y=sw_mean),method="nls",formula=y~shaft_sw_sc/(1+exp(-(4*(x-c))/w)),se=FALSE,method.args=list(start=list(c=490,w=490)),color=cbf_blue) +
  #geom_point(data=shaft_cline,aes(x=dist,y=uv_mean),color=cbf_purple) +
  #geom_smooth(data=shaft_cline,aes(x=dist,y=uv_mean),method="nls",formula=y~shaft_uv_sc/(1+exp(-(4*(x-c))/w)),se=FALSE,method.args=list(start=list(c=490,w=490)),color=cbf_purple) +
  geom_point(data=shaft_cline,aes(x=dist,y=lum_mean),color="gray44") +
  geom_smooth(data=shaft_cline,aes(x=dist,y=lum_mean),method="nls",formula=y~shaft_lum_sc/(1+exp(-(4*(x-c))/w)),se=FALSE,method.args=list(start=list(c=490,w=490)),color="gray44") +
  xlab("Distance (km)") +
  ylab("Photo value") +
  ggtitle("Shaft") +
  theme_bw() +
  theme(axis.title=element_text(face="bold",size=12), axis.text=element_text(size=10,color="black"))




### THROAT
# photo characters: lwMean, swMean, and lumMean
throat_cline <- cline_df %>%
  filter(trait=="throat")

throat_lw_sc <- mean(filter(throat_cline,site_ID>25)$lw_mean)
#throat_mw_sc <- mean(filter(throat_cline,site_ID>25)$mw_mean)
throat_sw_sc <- mean(filter(throat_cline,site_ID>25)$sw_mean)
#throat_uv_sc <- mean(filter(throat_cline,site_ID>25)$uv_mean)
throat_lum_sc <- mean(filter(throat_cline,site_ID>25)$lum_mean)

throat_cline_plot <- ggplot() +
  geom_point(data=throat_cline,aes(x=dist,y=lw_mean),color=cbf_red) +
  geom_smooth(data=throat_cline,aes(x=dist,y=lw_mean),method="nls",formula=y~throat_lw_sc/(1+exp(-(4*(x-c))/w)),se=FALSE,method.args=list(start=list(c=490,w=490)),color=cbf_red) +
  #geom_point(data=throat_cline,aes(x=dist,y=mw_mean),color=cbf_green) +
  #geom_smooth(data=throat_cline,aes(x=dist,y=mw_mean),method="nls",formula=y~throat_mw_sc/(1+exp(-(4*(x-c))/w)),se=FALSE,method.args=list(start=list(c=490,w=490)),color=cbf_green) +
  geom_point(data=throat_cline,aes(x=dist,y=sw_mean),color=cbf_blue) +
  geom_smooth(data=throat_cline,aes(x=dist,y=sw_mean),method="nls",formula=y~throat_sw_sc/(1+exp(-(4*(x-c))/w)),se=FALSE,method.args=list(start=list(c=490,w=490)),color=cbf_blue) +
  #geom_point(data=throat_cline,aes(x=dist,y=uv_mean),color=cbf_purple) +
  #geom_smooth(data=throat_cline,aes(x=dist,y=uv_mean),method="nls",formula=y~throat_uv_sc/(1+exp(-(4*(x-c))/w)),se=FALSE,method.args=list(start=list(c=490,w=490)),color=cbf_purple) +
  geom_point(data=throat_cline,aes(x=dist,y=lum_mean),color="gray44") +
  geom_smooth(data=throat_cline,aes(x=dist,y=lum_mean),method="nls",formula=y~throat_lum_sc/(1+exp(-(4*(x-c))/w)),se=FALSE,method.args=list(start=list(c=490,w=490)),color="gray44") +
  xlab("Distance (km)") +
  ylab("Photo value") +
  ggtitle("Throat") +
  theme_bw() +
  theme(axis.title=element_text(face="bold",size=12), axis.text=element_text(size=10,color="black"))




### CROWN
# photo characters: mwMean, swMean, uvMean, and lumMean are significant
crown_cline <- cline_df %>%
  filter(trait=="crown")

#crown_lw_sc <- mean(filter(crown_cline,site_ID>25)$lw_mean)
crown_mw_sc <- mean(filter(crown_cline,site_ID>25)$mw_mean)
crown_sw_sc <- mean(filter(crown_cline,site_ID>25)$sw_mean)
crown_uv_sc <- mean(filter(crown_cline,site_ID>25)$uv_mean)
crown_lum_sc <- mean(filter(crown_cline,site_ID>25)$lum_mean)

crown_cline_plot <- ggplot() +
  #geom_point(data=crown_cline,aes(x=dist,y=lw_mean),color=cbf_red) +
  #geom_smooth(data=crown_cline,aes(x=dist,y=lw_mean),method="nls",formula=y~crown_lw_sc/(1+exp(-(4*(x-c))/w)),se=FALSE,method.args=list(start=list(c=490,w=490)),color=cbf_red) +
  geom_point(data=crown_cline,aes(x=dist,y=mw_mean),color=cbf_green) +
  geom_smooth(data=crown_cline,aes(x=dist,y=mw_mean),method="nls",formula=y~crown_mw_sc/(1+exp(-(4*(x-c))/w)),se=FALSE,method.args=list(start=list(c=490,w=490)),color=cbf_green) +
  geom_point(data=crown_cline,aes(x=dist,y=sw_mean),color=cbf_blue) +
  geom_smooth(data=crown_cline,aes(x=dist,y=sw_mean),method="nls",formula=y~crown_sw_sc/(1+exp(-(4*(x-c))/w)),se=FALSE,method.args=list(start=list(c=490,w=490)),color=cbf_blue) +
  geom_point(data=crown_cline,aes(x=dist,y=uv_mean),color=cbf_purple) +
  geom_smooth(data=crown_cline,aes(x=dist,y=uv_mean),method="nls",formula=y~crown_uv_sc/(1+exp(-(4*(x-c))/w)),se=FALSE,method.args=list(start=list(c=490,w=490)),color=cbf_purple) +
  geom_point(data=crown_cline,aes(x=dist,y=lum_mean),color="gray44") +
  geom_smooth(data=crown_cline,aes(x=dist,y=lum_mean),method="nls",formula=y~crown_lum_sc/(1+exp(-(4*(x-c))/w)),se=FALSE,method.args=list(start=list(c=490,w=490)),color="gray44") +
  xlab("Distance (km)") +
  ylab("Photo value") +
  ggtitle("Crown") +
  theme_bw() +
  theme(axis.title=element_text(face="bold",size=12), axis.text=element_text(size=10,color="black"))






### EAR COVERTS
# photo traits: lwMean and swMean are significant
ear_cline <- cline_df %>%
  filter(trait=="ear")

ear_lw_sc <- mean(filter(ear_cline,site_ID>25)$lw_mean)
#ear_mw_sc <- mean(filter(ear_cline,site_ID>25)$mw_mean)
ear_sw_sc <- mean(filter(ear_cline,site_ID>25)$sw_mean)
#ear_uv_sc <- mean(filter(ear_cline,site_ID>25)$uv_mean)
#ear_lum_sc <- mean(filter(ear_cline,site_ID>25)$lum_mean)

ear_cline_plot <- ggplot() +
  geom_point(data=ear_cline,aes(x=dist,y=lw_mean),color=cbf_red) +
  geom_smooth(data=ear_cline,aes(x=dist,y=lw_mean),method="nls",formula=y~ear_lw_sc/(1+exp(-(4*(x-c))/w)),se=FALSE,method.args=list(start=list(c=490,w=490)),color=cbf_red) +
  #geom_point(data=ear_cline,aes(x=dist,y=mw_mean),color=cbf_green) +
  #geom_smooth(data=ear_cline,aes(x=dist,y=mw_mean),method="nls",formula=y~ear_mw_sc/(1+exp(-(4*(x-c))/w)),se=FALSE,method.args=list(start=list(c=490,w=490)),color=cbf_green) +
  geom_point(data=ear_cline,aes(x=dist,y=sw_mean),color=cbf_blue) +
  geom_smooth(data=ear_cline,aes(x=dist,y=sw_mean),method="nls",formula=y~ear_sw_sc/(1+exp(-(4*(x-c))/w)),se=FALSE,method.args=list(start=list(c=490,w=490)),color=cbf_blue) +
  #geom_point(data=ear_cline,aes(x=dist,y=uv_mean),color=cbf_purple) +
  #geom_smooth(data=ear_cline,aes(x=dist,y=uv_mean),method="nls",formula=y~ear_uv_sc/(1+exp(-(4*(x-c))/w)),se=FALSE,method.args=list(start=list(c=490,w=490)),color=cbf_purple) +
  #geom_point(data=ear_cline,aes(x=dist,y=lum_mean),color="gray44") +
  #geom_smooth(data=ear_cline,aes(x=dist,y=lum_mean),method="nls",formula=y~ear_lum_sc/(1+exp(-(4*(x-c))/w)),se=FALSE,method.args=list(start=list(c=490,w=490)),color="gray44") +
  xlab("Distance (km)") +
  ylab("Photo value") +
  ggtitle("Ear coverts") +
  theme_bw() +
  theme(axis.title=element_text(face="bold",size=12), axis.text=element_text(size=10,color="black"))





### MALE MALAR STRIPE

malar_cline <- photo_df %>%
  filter(sex=="Male") %>%
  filter(trait=="malar") %>%
  group_by(site_ID,site_name,dist,trait) %>%
  summarize(samples=n(),
            # photo traits
            lw_mean = mean(na.omit(lwMean)),
            lw_se = sd(na.omit(lwMean))/samples,
            mw_mean = mean(na.omit(mwMean)),
            mw_se = sd(na.omit(mwMean))/samples,
            sw_mean = mean(na.omit(swMean)),
            sw_se = sd(na.omit(swMean))/samples,
            uv_mean = mean(na.omit(uvMean)),
            uv_se = sd(na.omit(uvMean))/samples,
            lum_mean = mean(na.omit(lumMean)),
            lum_se = sd(na.omit(lumMean))/samples) %>%
  arrange(site_ID)

#photo characters: lwMean, mwMean, swMean, and lumMean are significant

malar_lw_sc <- mean(filter(malar_cline,site_ID<3)$lw_mean)
malar_mw_sc <- mean(filter(malar_cline,site_ID<3)$mw_mean)
malar_sw_sc <- mean(filter(malar_cline,site_ID<3)$sw_mean)
#malar_uv_sc <- mean(filter(malar_cline,site_ID<3)$uv_mean)
malar_lum_sc <- mean(filter(malar_cline,site_ID<3)$lum_mean)

malar_cline_plot <- ggplot() +
  geom_point(data=malar_cline,aes(x=dist,y=lw_mean),color=cbf_red) +
  geom_smooth(data=malar_cline,aes(x=dist,y=lw_mean),method="nls",formula=y~malar_lw_sc/(1+exp((4*(x-c))/w)),se=FALSE,method.args=list(start=list(c=490,w=490)),color=cbf_red) +
  geom_point(data=malar_cline,aes(x=dist,y=mw_mean),color=cbf_green) +
  geom_smooth(data=malar_cline,aes(x=dist,y=mw_mean),method="nls",formula=y~malar_mw_sc/(1+exp((4*(x-c))/w)),se=FALSE,method.args=list(start=list(c=490,w=490)),color=cbf_green) +
  geom_point(data=malar_cline,aes(x=dist,y=sw_mean),color=cbf_blue) +
  geom_smooth(data=malar_cline,aes(x=dist,y=sw_mean),method="nls",formula=y~malar_sw_sc/(1+exp((4*(x-c))/w)),se=FALSE,method.args=list(start=list(c=490,w=490)),color=cbf_blue) +
  #geom_point(data=malar_cline,aes(x=dist,y=uv_mean),color=cbf_purple) +
  #geom_smooth(data=malar_cline,aes(x=dist,y=uv_mean),method="nls",formula=y~malar_uv_sc/(1+exp((4*(x-c))/w)),se=FALSE,method.args=list(start=list(c=490,w=490)),color=cbf_purple) +
  geom_point(data=malar_cline,aes(x=dist,y=lum_mean),color="gray44") +
  geom_smooth(data=malar_cline,aes(x=dist,y=lum_mean),method="nls",formula=y~malar_lum_sc/(1+exp((4*(x-c))/w)),se=FALSE,method.args=list(start=list(c=490,w=490)),color="gray44") +
  xlab("Distance (km)") +
  ylab("Photo value") +
  ggtitle("Male malar stripe") +
  theme_bw() +
  theme(axis.title=element_text(face="bold",size=12), axis.text=element_text(size=10,color="black"))


# plot clines for all six traits
grid.arrange(crown_cline_plot,ear_cline_plot,malar_cline_plot,nuchal_cline_plot,shaft_cline_plot,throat_cline_plot,nrow=3)
# save plot as PDF with dimensions 6" x 7.75"




### run full nls modelling for two nice clines

# set up some parameters
# cline equation
rhs <- function(sc, x, c, w) {
  sc/(1+exp(-(4*(x-c))/w))
}
# cline equation 2
rhs2 <- function(sc, x, c, w) {
  sc/(1+exp((4*(x-c))/w))
}
# for plotting bootstrap
s <- seq(-100,880,length=100)



### NUCHAL PATCH
# photo characters: area
nuchal_cline <- cline_df %>%
  filter(trait=="nuchal") %>%
  select(site_ID:samples,starts_with("area"))

# scaling factor for end of the cline
# averages across the last three sites of YSFLs
nuchal_area_sc <- mean(filter(nuchal_cline,site_ID>25)$area_mean)


# modelling the cline
m_nuchal_area <- nls(area_mean ~ rhs(nuchal_area_sc,dist,center,width),
                     data=nuchal_cline,
                     start=list(center=490,width=490),
                     trace=T)

# summarizing the output
summary(m_nuchal_area)
coef(m_nuchal_area)

# calculating a confidence interval
CI_nuchal_area <- confint(m_nuchal_area,parm=c("center","width"))
CI_nuchal_area

# bootstrapping the output
bootstrap_nuchal_area <- nlsBoot(m_nuchal_area,niter=999)
summary(bootstrap_nuchal_area)
boot_coeff_nuchal_area <- bootstrap_nuchal_area$coefboot

# bootstrapping for plotting
y_nuchal_area <- matrix(NA,nrow=length(bootstrap_nuchal_area),ncol=length(s))
for(i in 1:length(bootstrap_nuchal_area)){
  y_nuchal_area[i,]=rhs(nuchal_area_sc,s,boot_coeff_nuchal_area[i,1],boot_coeff_nuchal_area[i,2])
}
boot_ci_nuchal_area <- apply(y_nuchal_area,2,quantile,c(0.025,0.975))
boot_df_nuchal_area <- data.frame(s,t(boot_ci_nuchal_area))
names(boot_df_nuchal_area) = c("dist","lower","upper")

sink("./cline-output-photos/nuchal_area_output.txt")
print('##### MODEL SUMMARY #####')
summary(m_nuchal_area)
print('##### COEFF FROM MODEL #####')
coef(m_nuchal_area)
print('##### CI FROM MODEL #####')
CI_nuchal_area
print('##### BOOTSTRAP SUMMARY #####')
summary(bootstrap_nuchal_area)
sink()

write.table(boot_df_nuchal_area,"./cline-output-photos/bootstrap-CI_nuchal_area.txt",row.names=FALSE,sep="\t",quote=FALSE)

# cline plot
nuchal_cline <- ggplot() +
  geom_ribbon(data=boot_df_nuchal_area,aes(x=dist,ymin=lower,ymax=upper),alpha=0.5,fill="gray") +
  geom_smooth(data=nuchal_cline,aes(x=dist,y=area_mean),method="nls",formula=y~nuchal_area_sc/(1+exp(-(4*(x-c))/w)),se=FALSE,method.args=list(start=list(c=490,w=490)), color="black") +
  geom_point(data=nuchal_cline,aes(x=dist,y=area_mean)) +
  xlab("Distance (km)") +
  ylab("Area") +
  ggtitle("Nuchal patch") +
  theme_bw() +
  theme(axis.title=element_text(face="bold",size=12), axis.text=element_text(size=10,color="black"))






### MALE MALAR STRIPE
# photo characters: lwmean
malar_cline <- photo_df %>%
  filter(sex=="Male") %>%
  filter(trait=="malar") %>%
  group_by(site_ID,site_name,dist,trait) %>%
  summarize(samples=n(),
            lw_mean = mean(na.omit(lwMean)),
            lw_se = sd(na.omit(lwMean))/samples) %>%
  arrange(site_ID)

# scaling factor for end of the cline
# averages across the last three sites of YSFLs
malar_lw_sc <- mean(filter(malar_cline,site_ID<3)$lw_mean)


# modelling the cline
m_malar_lw <- nls(lw_mean ~ rhs2(malar_lw_sc,dist,center,width),
                     data=malar_cline,
                     start=list(center=490,width=490),
                     trace=T)

# summarizing the output
summary(m_malar_lw)
coef(m_malar_lw)

# calculating a confidence interval
CI_malar_lw <- confint(m_malar_lw,parm=c("center","width"))
CI_malar_lw

# bootstrapping the output
bootstrap_malar_lw <- nlsBoot(m_malar_lw,niter=999)
summary(bootstrap_malar_lw)
boot_coeff_malar_lw <- bootstrap_malar_lw$coefboot

# bootstrapping for plotting
y_malar_lw <- matrix(NA,nrow=length(bootstrap_malar_lw),ncol=length(s))
for(i in 1:length(bootstrap_malar_lw)){
  y_malar_lw[i,]=rhs2(malar_lw_sc,s,boot_coeff_malar_lw[i,1],boot_coeff_malar_lw[i,2])
}
boot_ci_malar_lw <- apply(y_malar_lw,2,quantile,c(0.025,0.975))
boot_df_malar_lw <- data.frame(s,t(boot_ci_malar_lw))
names(boot_df_malar_lw) = c("dist","lower","upper")

sink("./cline-output-photos/malar_lw_output.txt")
print('##### MODEL SUMMARY #####')
summary(m_malar_lw)
print('##### COEFF FROM MODEL #####')
coef(m_malar_lw)
print('##### CI FROM MODEL #####')
CI_malar_lw
print('##### BOOTSTRAP SUMMARY #####')
summary(bootstrap_malar_lw)
sink()

write.table(boot_df_malar_lw,"./cline-output-photos/bootstrap-CI_malar_lw.txt",row.names=FALSE,sep="\t",quote=FALSE)

# cline plot
malar_cline <- ggplot() +
  geom_ribbon(data=boot_df_malar_lw,aes(x=dist,ymin=lower,ymax=upper),alpha=0.5,fill="gray") +
  geom_smooth(data=malar_cline,aes(x=dist,y=lw_mean),method="nls",formula=y~malar_lw_sc/(1+exp((4*(x-c))/w)),se=FALSE,method.args=list(start=list(c=490,w=490)), color=cbf_red) +
  geom_point(data=malar_cline,aes(x=dist,y=lw_mean),color=cbf_red) +
  xlab("Distance (km)") +
  ylab("lwMean") +
  ggtitle("Male malar stripe") +
  theme_bw() +
  theme(axis.title=element_text(face="bold",size=12), axis.text=element_text(size=10,color="black"))


# plot two clines together
grid.arrange(nuchal_cline, malar_cline, nrow=2)
# save plot as PDF with dimensions 4" x 8"



#### COMPARE PHOTO CLINES WITH SCORING CLINES

comparison_df <- read_tsv("./cline-output-photos/model-output-comparisons.txt", col_names=TRUE)

ggplot(data=comparison_df, aes(x=trait,y=center,color=data_source)) +
  geom_pointrange(aes(ymin=center_low,ymax=center_high),position=position_dodge(width=0.5)) +
  scale_color_manual(values=c(cbf_red,"black")) +
  coord_flip() +
  ylim(c(0,800)) +
  xlab("Trait") +
  ylab("Distance (km)") +
  labs(color="Data Source") +
  theme_bw() +
  theme(axis.title=element_text(face="bold",size=12), axis.text=element_text(size=10,color="black"),
        legend.position=c(0.8, 0.8), legend.background=element_rect(linetype="solid",color="black"))
# save plot as PDF with dimensions 3.5" x 4"








########### STEP 3 ###########


# NUCHAL PATCH
nuchal_df <- photo_df %>%
  filter(trait=="nuchal")

# regression
summary(lm(nuchal ~ area,data=nuchal_df))
#Call:
#  lm(formula = nuchal ~ area, data = nuchal_df)
#
#Residuals:
#  Min       1Q   Median       3Q      Max
#-2.86025 -0.65480 -0.03915  0.84335  2.39893
#
#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)
#(Intercept)  3.157e+00  1.563e-01    20.2   <2e-16 ***
#  area        -2.562e-04  1.612e-05   -15.9   <2e-16 ***
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#Residual standard error: 0.883 on 100 degrees of freedom
#Multiple R-squared:  0.7165,	Adjusted R-squared:  0.7136
#F-statistic: 252.7 on 1 and 100 DF,  p-value: < 2.2e-16




# SHAFT
shaft_df <- photo_df %>%
  filter(trait=="vane")

# regression
summary(lm(shaft ~ lwMean + mwMean + swMean + uvMean + lumMean, data=shaft_df))
#Call:
#  lm(formula = shaft ~ lwMean + mwMean + swMean + uvMean + lumMean,
#     data = shaft_df)
#
#Residuals:
#  Min       1Q   Median       3Q      Max
#-0.99768 -0.24735  0.03285  0.25263  1.86768
#
#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)
#(Intercept) -5.524e-02  1.556e+00  -0.036 0.971749
#lwMean       4.434e+02  2.274e+02   1.950 0.054149 .
#mwMean       7.549e+02  4.608e+02   1.638 0.104644
#swMean       3.153e+02  1.442e+02   2.186 0.031237 *
#  uvMean      -2.978e+01  7.549e+00  -3.945 0.000152 ***
#  lumMean     -1.472e+03  8.218e+02  -1.792 0.076337 .
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#Residual standard error: 0.452 on 96 degrees of freedom
#Multiple R-squared:  0.9242,	Adjusted R-squared:  0.9203
#F-statistic: 234.2 on 5 and 96 DF,  p-value: < 2.2e-16




# THROAT
throat_df <- photo_df %>%
  filter(trait=="throat")

# regression
summary(lm(throat ~ lwMean + mwMean + swMean + uvMean + lumMean, data=throat_df))
#Call:
#  lm(formula = throat ~ lwMean + mwMean + swMean + uvMean + lumMean,
#     data = throat_df)
#
#Residuals:
#  Min       1Q   Median       3Q      Max
#-2.18396 -0.37330  0.02787  0.35565  2.17901
#
#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)
#(Intercept)    4.658      3.018   1.544    0.126
#lwMean      -474.791    449.144  -1.057    0.293
#mwMean      -626.573    925.556  -0.677    0.500
#swMean      -282.521    294.289  -0.960    0.339
#uvMean        42.794     29.539   1.449    0.151
#lumMean     1344.392   1632.641   0.823    0.412
#
#Residual standard error: 0.6098 on 96 degrees of freedom
#Multiple R-squared:  0.8784,	Adjusted R-squared:  0.872
#F-statistic: 138.7 on 5 and 96 DF,  p-value: < 2.2e-16




# CROWN
crown_df <- photo_df %>%
  filter(trait=="crown")

# regression
summary(lm(crown ~ lwMean + mwMean + swMean + uvMean + lumMean, data=crown_df))
#Call:
#  lm(formula = crown ~ lwMean + mwMean + swMean + uvMean + lumMean,
#     data = crown_df)
#
#Residuals:
#  Min      1Q  Median      3Q     Max
#-2.8623 -0.8427 -0.2223  0.7704  3.2035
#
#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)
#(Intercept)    21.878      9.707   2.254   0.0265 *
#  lwMean      -1683.623   1395.911  -1.206   0.2307
#mwMean      -4335.483   2954.701  -1.467   0.1456
#swMean      -1144.436    901.443  -1.270   0.2073
#uvMean         52.657     84.197   0.625   0.5332
#lumMean      7014.261   5141.696   1.364   0.1757
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#Residual standard error: 1.226 on 96 degrees of freedom
#Multiple R-squared:  0.4773,	Adjusted R-squared:  0.4501
#F-statistic: 17.53 on 5 and 96 DF,  p-value: 2.653e-12




# EAR COVERTS
ear_df <- photo_df %>%
  filter(trait=="ear") %>%
  distinct()

# regression
summary(lm(ear_coverts ~ lwMean + mwMean + swMean + uvMean + lumMean, data=ear_df))
#Call:
#  lm(formula = ear_coverts ~ lwMean + mwMean + swMean + uvMean +
#       lumMean, data = ear_df)
#
#Residuals:
#  Min       1Q   Median       3Q      Max
#-2.21057 -0.39170  0.06766  0.46375  2.02944
#
#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)
#(Intercept)     5.493      5.585   0.984    0.328
#lwMean       -575.062    867.247  -0.663    0.509
#mwMean      -1055.446   1776.178  -0.594    0.554
#swMean       -228.953    562.000  -0.407    0.685
#uvMean         -3.106     49.746  -0.062    0.950
#lumMean      1868.908   3144.036   0.594    0.554
#
#Residual standard error: 0.8745 on 95 degrees of freedom
#Multiple R-squared:  0.7454,	Adjusted R-squared:  0.732
#F-statistic: 55.62 on 5 and 95 DF,  p-value: < 2.2e-16




# MALE MALAR
malar_df <- photo_df %>%
  filter(trait=="malar") %>%
  filter(sex=="Male") %>%
  distinct()

# regression
summary(lm(malar ~ lwMean + mwMean + swMean + uvMean + lumMean, data=malar_df))
#Call:
#  lm(formula = malar ~ lwMean + mwMean + swMean + uvMean + lumMean,
#     data = malar_df)
#
#Residuals:
#  Min      1Q  Median      3Q     Max
#-1.1716 -0.4353 -0.1010  0.4375  1.4345
#
#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)
#(Intercept)     9.571      6.369   1.503    0.138
#lwMean      -1300.276   1015.513  -1.280    0.205
#mwMean      -2728.244   2059.983  -1.324    0.190
#swMean       -837.926    646.235  -1.297    0.200
#uvMean         42.027     52.438   0.801    0.426
#lumMean      4801.124   3657.776   1.313    0.194
#
#Residual standard error: 0.5855 on 60 degrees of freedom
#Multiple R-squared:  0.8907,	Adjusted R-squared:  0.8815
#F-statistic: 97.75 on 5 and 60 DF,  p-value: < 2.2e-16
