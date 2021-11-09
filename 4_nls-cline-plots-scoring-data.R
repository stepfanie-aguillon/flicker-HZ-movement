# This script contains the R code to plot nls cline output.
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
library(ghibli)
library(gridExtra)

# load datasets
contemporary <- read_tsv("scoring_contemporary_input.txt",col_names=TRUE)
historic <- read_tsv("scoring_historic_input.txt",col_names=TRUE)

# load boostrapping results
bootstrap_c <- read_tsv("./cline-output-scores/bootstrap-CI_standHI_c.txt",col_names=TRUE)
bootstrap_h <- read_tsv("./cline-output-scores/bootstrap-CI_standHI_h.txt",col_names=TRUE)

# limit bootstrap for historic (to match observations)
bootstrap_h <- bootstrap_h %>%
  filter(dist>56)

# load individual trait datasets
models <- read_tsv("./cline-output-scores/model-output-summary.txt",col_names=TRUE)
models <- models %>%
  group_by(transect)


# historic vs. contemporary cline comparison (in text)
# standardized hybrid index

fig2a <- ggplot() +
  # historic, standard HI, GRAY
  geom_ribbon(data=bootstrap_h,aes(x=dist,ymin=lower,ymax=upper),alpha=0.25,fill="gray55") +
  geom_smooth(data=historic,aes(x=dist,y=standHI_mean),method="nls",formula=y~1/(1+exp((4*(x-c))/w)),se=FALSE,method.args=list(start=list(c=360,w=360)),color="gray55",size=1) +
  geom_pointrange(data=historic,aes(x=dist,y=standHI_mean,ymin=standHI_mean-standHI_se,ymax=standHI_mean+standHI_se),color="gray55",shape=17) +

  # contemporary, standard HI, BLACK
  geom_ribbon(data=bootstrap_c,aes(x=dist,ymin=lower,ymax=upper),alpha=0.25,fill="black") +
  geom_smooth(data=contemporary,aes(x=dist,y=standHI_mean),method="nls",formula=y~1/(1+exp((4*(x-c))/w)),se=FALSE,method.args=list(start=list(c=360,w=360)),color="black",size=1) +
  geom_pointrange(data=contemporary,aes(x=dist,y=standHI_mean,ymin=standHI_mean-standHI_se,ymax=standHI_mean+standHI_se),color="black",shape=16) +

  ylim(c(-0.01,1.01)) +
  xlim(c(0,800)) +
  xlab("Distance (km)") +
  ylab("Hybrid index") +
  theme_bw() +
  theme(axis.title=element_text(face="bold",size=12), axis.text=element_text(size=10,color="black"))




# historic vs. contemporary cline centers (in text)
# standardized hybrid index AND individual traits

# order traits for plotting
models$trait <- factor(models$trait, levels=c("throat","shaft","nuchal","malar","ear coverts","crown","hybrid index"))

fig2b <- ggplot(data=models, aes(x=trait,y=center, color=transect,shape=transect)) +
  geom_pointrange(aes(ymin=center_low,ymax=center_high),position=position_dodge(width=0.5)) +
  scale_color_manual(values=c("black","gray55")) +
  scale_shape_manual(values=c(16,17)) +
  coord_flip() +
  ylim(c(0,400)) +
  xlab("Trait") +
  ylab("Distance (km)") +
  theme_bw() +
  theme(legend.position="none", axis.title=element_text(face="bold",size=12), axis.text=element_text(size=10,color="black"))


grid.arrange(fig2a,fig2b,nrow=1)
# save plot as PDF with dimensions 8" x 4"






# individual trait clines (supplement)
# historic transect
historic_clines <- ggplot() +
  # shaft, GREEN
  geom_smooth(data=historic,aes(x=dist,y=shaft_mean),method="nls",formula=y~4/(1+exp((4*(x-c))/w)),se=FALSE,method.args=list(start=list(c=360,w=360)),color=ghibli_palettes$MarnieMedium2[6],size=1) +
  geom_jitter(data=historic,aes(x=dist,y=shaft_mean),color=ghibli_palettes$MarnieMedium2[6],shape=17,size=2) +

  # throat, BLUE
  geom_smooth(data=historic,aes(x=dist,y=throat_mean),method="nls",formula=y~4/(1+exp((4*(x-c))/w)),se=FALSE,method.args=list(start=list(c=360,w=360)),color=ghibli_palettes$PonyoMedium[3],size=1) +
  geom_jitter(data=historic,aes(x=dist,y=throat_mean),color=ghibli_palettes$PonyoMedium[3],shape=17,size=2) +

  # nuchal patch, RED
  geom_smooth(data=historic,aes(x=dist,y=nuchal_mean),method="nls",formula=y~4/(1+exp((4*(x-c))/w)),se=FALSE,method.args=list(start=list(c=360,w=360)),color=ghibli_palettes$KikiMedium[3],size=1) +
  geom_jitter(data=historic,aes(x=dist,y=nuchal_mean),color=ghibli_palettes$KikiMedium[3],shape=17,size=2) +

  # crown, NAVY
  geom_smooth(data=historic,aes(x=dist,y=crown_mean),method="nls",formula=y~4/(1+exp((4*(x-c))/w)),se=FALSE,method.args=list(start=list(c=360,w=360)),color=ghibli_palettes$PonyoMedium[2],size=1) +
  geom_jitter(data=historic,aes(x=dist,y=crown_mean),color=ghibli_palettes$PonyoMedium[2],shape=17,size=2) +

  # ear covert, PEACH
  geom_smooth(data=historic,aes(x=dist,y=ears_mean),method="nls",formula=y~4/(1+exp((4*(x-c))/w)),se=FALSE,method.args=list(start=list(c=360,w=360)),color=ghibli_palettes$PonyoMedium[5],size=1) +
  geom_jitter(data=historic,aes(x=dist,y=ears_mean),color=ghibli_palettes$PonyoMedium[5],shape=17,size=2) +

  # male malar, YELLOW
  geom_smooth(data=historic,aes(x=dist,y=malar_mean),method="nls",formula=y~4/(1+exp((4*(x-c))/w)),se=FALSE,method.args=list(start=list(c=360,w=360)),color=ghibli_palettes$PonyoMedium[6],size=1) +
  geom_jitter(data=historic,aes(x=dist,y=malar_mean),color=ghibli_palettes$PonyoMedium[6],shape=17,size=2) +

  ylim(c(-0.1,4.1)) +
  xlim(c(0,800)) +
  xlab("Distance (km)") +
  ylab("Trait score") +
  theme_bw() +
  theme(axis.title=element_text(face="bold",size=12), axis.text=element_text(size=10,color="black"))


# individual trait clines (supplement)
# contemporary transect
contemporary_clines <- ggplot() +
  # shaft, GREEN
  geom_smooth(data=contemporary,aes(x=dist,y=shaft_mean),method="nls",formula=y~4/(1+exp((4*(x-c))/w)),se=FALSE,method.args=list(start=list(c=360,w=360)),color=ghibli_palettes$MarnieMedium2[6],size=1) +
  geom_jitter(data=contemporary,aes(x=dist,y=shaft_mean),color=ghibli_palettes$MarnieMedium2[6],shape=16,size=2) +

  # throat, BLUE
  geom_smooth(data=contemporary,aes(x=dist,y=throat_mean),method="nls",formula=y~4/(1+exp((4*(x-c))/w)),se=FALSE,method.args=list(start=list(c=360,w=360)),color=ghibli_palettes$PonyoMedium[3],size=1) +
  geom_jitter(data=contemporary,aes(x=dist,y=throat_mean),color=ghibli_palettes$PonyoMedium[3],shape=16,size=2) +

  # nuchal patch, RED
  geom_smooth(data=contemporary,aes(x=dist,y=nuchal_mean),method="nls",formula=y~4/(1+exp((4*(x-c))/w)),se=FALSE,method.args=list(start=list(c=360,w=360)),color=ghibli_palettes$KikiMedium[3],size=1) +
  geom_jitter(data=contemporary,aes(x=dist,y=nuchal_mean),color=ghibli_palettes$KikiMedium[3],shape=16,size=2) +

  # crown, NAVY
  geom_smooth(data=contemporary,aes(x=dist,y=crown_mean),method="nls",formula=y~4/(1+exp((4*(x-c))/w)),se=FALSE,method.args=list(start=list(c=360,w=360)),color=ghibli_palettes$PonyoMedium[2],size=1) +
  geom_jitter(data=contemporary,aes(x=dist,y=crown_mean),color=ghibli_palettes$PonyoMedium[2],shape=16,size=2) +

  # ear covert, PEACH
  geom_smooth(data=contemporary,aes(x=dist,y=ears_mean),method="nls",formula=y~4/(1+exp((4*(x-c))/w)),se=FALSE,method.args=list(start=list(c=360,w=360)),color=ghibli_palettes$PonyoMedium[5],size=1) +
  geom_jitter(data=contemporary,aes(x=dist,y=ears_mean),color=ghibli_palettes$PonyoMedium[5],shape=16,size=2) +

  # male malar, YELLOW
  geom_smooth(data=contemporary,aes(x=dist,y=malar_mean),method="nls",formula=y~4/(1+exp((4*(x-c))/w)),se=FALSE,method.args=list(start=list(c=360,w=360)),color=ghibli_palettes$PonyoMedium[6],size=1) +
  geom_jitter(data=contemporary,aes(x=dist,y=malar_mean),color=ghibli_palettes$PonyoMedium[6],shape=16,size=2) +

  ylim(c(-0.1,4.1)) +
  xlim(c(0,800)) +
  xlab("Distance (km)") +
  ylab("Trait score") +
  theme_bw() +
  theme(axis.title=element_text(face="bold",size=12), axis.text=element_text(size=10,color="black"))

grid.arrange(historic_clines,contemporary_clines,nrow=2)
# save plot as PDF with dimensions 5" x 8"






# influence of including samples on the North Platte (supplement)

# load dataset, edit to add in info on North v South Platte
contemporary <- read_tsv("scoring_contemporary_input.txt",col_names=TRUE)
contemporary <- contemporary %>%
  mutate(Platte = ifelse(site_ID %in% c(6,9,11,13,16,17),"North Platte","South Platte"))
# filter dataset to remove samples on the North Platte
S_contemporary <- contemporary %>%
  filter(!site_ID %in% c(6,9,11,13,16,17))

# load boostrapping results
bootstrap_c <- read_tsv("./cline-output-scores/bootstrap-CI_standHI_c.txt",col_names=TRUE)
bootstrap_cS <- read_tsv("./cline-output-scores/bootstrap-CI_standHI_SPlatte_c.txt",col_names=TRUE)

# standardized hybrid index
ggplot() +
  # contemporary, standard HI, BLACK
  geom_ribbon(data=bootstrap_c,aes(x=dist,ymin=lower,ymax=upper),alpha=0.25,fill="gray44") +
  geom_smooth(data=contemporary,aes(x=dist,y=standHI_mean),method="nls",formula=y~1/(1+exp((4*(x-c))/w)),se=FALSE,method.args=list(start=list(c=360,w=360)),color="black",size=2) +
  geom_pointrange(data=contemporary,aes(x=dist,y=standHI_mean,ymin=standHI_mean-standHI_se,ymax=standHI_mean+standHI_se,shape=Platte),color="black") +
  
  # make the North Platte localities open circles
  scale_shape_manual(values=c(21,16)) +
  
  # contemporary without North Platte, standard HI, RED
  geom_ribbon(data=bootstrap_cS,aes(x=dist,ymin=lower,ymax=upper),alpha=0.25,fill="gray44") +
  geom_smooth(data=S_contemporary,aes(x=dist,y=standHI_mean),method="nls",formula=y~1/(1+exp((4*(x-c))/w)),se=FALSE,method.args=list(start=list(c=360,w=360)),color="#DC3220",size=1,linetype="longdash") +

  ylim(c(-0.01,1.01)) +
  xlim(c(-10,810)) +
  xlab("Distance (km)") +
  ylab("Hybrid index") +
  theme_bw() +
  theme(axis.title=element_text(face="bold",size=12), axis.text=element_text(size=10,color="black"),legend.position="none")
# save plot as PDF with dimensions 4" x 5"
