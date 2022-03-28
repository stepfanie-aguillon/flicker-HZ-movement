# This script contains the R code to compare only sites with repeat sampling efforts
#
# Authors: Aguillon SM, VG Rohwer
# Year: 2022
# Title: Revisiting a classic hybrid zone: movement of the
#        northern flicker hybrid zone in contemporary times
# Journal Info: Evolution
# DOI: 10.1111/evo.14474
#
# Edited date: Oct 2021
#
# Please cite the paper if you use these scripts
#


# load packages 
library(tidyverse)
library(ghibli)
library(nlstools)
library(gridExtra)



##### CLINES WITH ONLY THE SITES WITH REPEAT SAMPLING


# load datasets
contemporary <- read_tsv("scoring_contemporary_input.txt",col_names=TRUE)
historic <- read_tsv("scoring_historic_input.txt",col_names=TRUE)

# filter dataset to keep only sampling sites with repeat sampling
localities <- read_tsv("./scoring_locality_distances.txt") %>%
  mutate(both_transects = ifelse(historic==contemporary,"both","one")) %>%
  filter(both_transects == "both")

repeat_sites <- localities$site_ID

contemporary_repeat <- contemporary %>%
  filter(site_ID %in% repeat_sites)
historic_repeat <- historic %>%
  filter(site_ID %in% repeat_sites)



# set up some parameters
# cline for HI
rhs <- function(x, c, w) {
  1/(1+exp((4*(x-c))/w))
}
# for plotting bootstrap
s <- seq(0,800,length=100)



### CONTEMPORARY HI
# modelling the cline
standHI_c <- nls(standHI_mean ~ rhs(dist, center, width), 
                       data=contemporary_repeat,
                       start=list(center=390,width=390),
                       trace=T)

# summarizing the output
summary(standHI_c)
coef(standHI_c)

# calculating a confidence interval
CI_standHI_c <- confint(standHI_c,parm=c("center","width"))
CI_standHI_c

# bootstrapping the output
bootstrap_standHI_c <- nlsBoot(standHI_c,niter=999)
summary(bootstrap_standHI_c)
boot_coeff_standHI_c <- bootstrap_standHI_c$coefboot

# bootstrapping for plotting
y_standHI_c <- matrix(NA,nrow=length(bootstrap_standHI_c),ncol=length(s))
for(i in 1:length(bootstrap_standHI_c)){
  y_standHI_c[i,]=rhs(s,boot_coeff_standHI_c[i,1],boot_coeff_standHI_c[i,2])
}
boot_ci_standHI_c <- apply(y_standHI_c,2,quantile,c(0.025,0.975))
boot_df_standHI_c <- data.frame(s,t(boot_ci_standHI_c))
names(boot_df_standHI_c) = c("dist","lower","upper")

sink("./cline-output-scores/standHI_c_repeatsites_output.txt")
print('##### MODEL SUMMARY #####')
summary(standHI_c)
print('##### COEFF FROM MODEL #####')
coef(standHI_c)
print('##### CI FROM MODEL #####')
CI_standHI_c
print('##### BOOTSTRAP SUMMARY #####')
summary(bootstrap_standHI_c)
sink()

write.table(boot_df_standHI_c,"./cline-output-scores/bootstrap-CI_standHI_repeatsites_c.txt",row.names=FALSE,sep="\t",quote=FALSE)



### HISTORIC HI
# modelling the cline
standHI_h <- nls(standHI_mean ~ rhs(dist, center, width), 
                       data=historic_repeat,
                       start=list(center=390,width=390),
                       trace=T)

# summarizing the output
summary(standHI_h)
coef(standHI_h)

# calculating a confidence interval
CI_standHI_h <- confint(standHI_h,parm=c("center","width"))
CI_standHI_h

# bootstrapping the output
bootstrap_standHI_h <- nlsBoot(standHI_h,niter=999)
summary(bootstrap_standHI_h)
boot_hoeff_standHI_h <- bootstrap_standHI_h$coefboot

# bootstrapping for plotting
y_standHI_h <- matrix(NA,nrow=length(bootstrap_standHI_h),ncol=length(s))
for(i in 1:length(bootstrap_standHI_h)){
  y_standHI_h[i,]=rhs(s,boot_hoeff_standHI_h[i,1],boot_hoeff_standHI_h[i,2])
}
boot_hi_standHI_h <- apply(y_standHI_h,2,quantile,c(0.025,0.975))
boot_df_standHI_h <- data.frame(s,t(boot_hi_standHI_h))
names(boot_df_standHI_h) = c("dist","lower","upper")

sink("./cline-output-scores/standHI_h_repeatsites_output.txt")
print('##### MODEL SUMMARY #####')
summary(standHI_h)
print('##### COEFF FROM MODEL #####')
coef(standHI_h)
print('##### CI FROM MODEL #####')
CI_standHI_h
print('##### BOOTSTRAP SUMMARY #####')
summary(bootstrap_standHI_h)
sink()

write.table(boot_df_standHI_h,"./cline-output-scores/bootstrap-CI_standHI_repeatsites_h.txt",row.names=FALSE,sep="\t",quote=FALSE)



##### PLOTTING THE CLINES

# load bootstrap output from original modelling
bootstrap_c <- read_tsv("./cline-output-scores/bootstrap-CI_standHI_c.txt",col_names=TRUE)
bootstrap_h <- read_tsv("./cline-output-scores/bootstrap-CI_standHI_h.txt",col_names=TRUE)

# load bootstrap output from repeat sampling sites modelling (from above)
bootstrap_c_repeat <- read_tsv("./cline-output-scores/bootstrap-CI_standHI_repeatsites_c.txt",col_names=TRUE)
bootstrap_h_repeat <- read_tsv("./cline-output-scores/bootstrap-CI_standHI_repeatsites_h.txt",col_names=TRUE)


ggplot() +
  ### HISTORIC, grayed colors
  # full dataset, gray
  geom_ribbon(data=bootstrap_h,aes(x=dist,ymin=lower,ymax=upper),alpha=0.25,fill="gray44") +
  geom_smooth(data=historic,aes(x=dist,y=standHI_mean),method="nls",formula=y~1/(1+exp((4*(x-c))/w)),se=FALSE,method.args=list(start=list(c=360,w=360)),color="gray55",size=2) +

  # only repeat sampling sites, pink
  geom_ribbon(data=bootstrap_h_repeat,aes(x=dist,ymin=lower,ymax=upper),alpha=0.25,fill="gray44") +
  geom_smooth(data=historic_repeat,aes(x=dist,y=standHI_mean),method="nls",formula=y~1/(1+exp((4*(x-c))/w)),se=FALSE,method.args=list(start=list(c=360,w=360)),color="#e86e61",size=1,linetype="longdash") +
  
  ### CONTEMPORARY, dark colors
  # full dataset, black
  geom_ribbon(data=bootstrap_c,aes(x=dist,ymin=lower,ymax=upper),alpha=0.25,fill="gray44") +
  geom_smooth(data=contemporary,aes(x=dist,y=standHI_mean),method="nls",formula=y~1/(1+exp((4*(x-c))/w)),se=FALSE,method.args=list(start=list(c=360,w=360)),color="black",size=2) +
  
  # only repeat sampling sites, red
  geom_ribbon(data=bootstrap_c_repeat,aes(x=dist,ymin=lower,ymax=upper),alpha=0.25,fill="gray44") +
  geom_smooth(data=contemporary_repeat,aes(x=dist,y=standHI_mean),method="nls",formula=y~1/(1+exp((4*(x-c))/w)),se=FALSE,method.args=list(start=list(c=360,w=360)),color="#DC3220",size=1,linetype="longdash") +
  
  ylim(c(-0.01,1.01)) +
  xlim(c(-10,810)) +
  xlab("Distance (km)") +
  ylab("Hybrid index") +
  theme_bw() +
  theme(axis.title=element_text(face="bold",size=12), axis.text=element_text(size=10,color="black"))
# save plot as PDF with dimensions 4" x 5"






#### WILCOXON TESTS TO DIRECTLY COMPARE SITES WITH REPEAT SAMPLING

# load in dataset
scoring_df <- read_tsv("./raw-data/flicker-scoring-HZonly.txt",col_names=TRUE)

# load dataset with info on samples present in each transect
localities <- read_tsv("./scoring_locality_distances.txt") %>%
  mutate(both_transects = ifelse(historic==contemporary,"both","one")) %>%
  filter(both_transects == "both")

repeat_sites <- localities$site_ID


# plotting comparisons

# subset dataset to include only desired site localities
plotting_df <- scoring_df %>%
  filter(site_ID %in% repeat_sites)

# re-order factors for correct plotting
plotting_df$HZtransect <- factor(plotting_df$HZtransect, levels=c("historic","contemporary"))

# make plot
ggplot() + 
  geom_boxplot(data=plotting_df,aes(x=HZtransect,y=standHI,color=HZtransect),outlier.shape=NA) +
  geom_jitter(data=plotting_df,aes(x=HZtransect,y=standHI,color=HZtransect),width=0.25) +
  facet_wrap(~site_ID, nrow=2) +
  scale_color_manual(values=c("gray55","black")) +
  ylim(c(-0.01,1.01)) +  
  #xlab("Transect") +
  ylab("Hybrid index") +
  theme_bw() +
  theme(legend.position="none", axis.title=element_text(face="bold",size=11), axis.title.x=element_blank(), 
        axis.text=element_text(size=9,color="black"), strip.background=element_rect(fill="white"), 
        strip.text.x=element_text(size=9,color="black"))
# save plot as PDF with dimensions 5" x 8"




### WILCOXON TESTS

# site 4
site4 <- filter(scoring_df,site_ID==4) # n = 9
wilcox.test(filter(site4,HZtransect=="contemporary")$standHI,filter(site4,HZtransect=="historic")$standHI)
#data:  filter(site4, HZtransect == "contemporary")$standHI and filter(site4, HZtransect == "historic")$standHI
#W = 2, p-value = 0.184

# site 12
site12 <- filter(scoring_df,site_ID==12) # n = 24
wilcox.test(filter(site12,HZtransect=="contemporary")$standHI,filter(site12,HZtransect=="historic")$standHI)
#data:  filter(site12, HZtransect == "contemporary")$standHI and filter(site12, HZtransect == "historic")$standHI
#W = 27, p-value = 0.0164

# site 15
site15 <- filter(scoring_df,site_ID==15) # n = 9
wilcox.test(filter(site15,HZtransect=="contemporary")$standHI,filter(site15,HZtransect=="historic")$standHI)
#data:  filter(site15, HZtransect == "contemporary")$standHI and filter(site15, HZtransect == "historic")$standHI
#W = 2.5, p-value = 0.2356

# site 18
site18 <- filter(scoring_df,site_ID==18) # n = 32
wilcox.test(filter(site18,HZtransect=="contemporary")$standHI,filter(site18,HZtransect=="historic")$standHI)
#data:  filter(site18, HZtransect == "contemporary")$standHI and filter(site18, HZtransect == "historic")$standHI
#W = 53.5, p-value = 0.4295

# site 21
site21 <- filter(scoring_df,site_ID==21) # n = 38
wilcox.test(filter(site21,HZtransect=="contemporary")$standHI,filter(site21,HZtransect=="historic")$standHI)
#data:  filter(site21, HZtransect == "contemporary")$standHI and filter(site21, HZtransect == "historic")$standHI
#W = 32.5, p-value = 0.02334

# site 22
site22 <- filter(scoring_df,site_ID==22) # n = 28
wilcox.test(filter(site22,HZtransect=="contemporary")$standHI,filter(site22,HZtransect=="historic")$standHI)
#data:  filter(site22, HZtransect == "contemporary")$standHI and filter(site22, HZtransect == "historic")$standHI
#W = 51.5, p-value = 0.7242

# site 24
site24 <- filter(scoring_df,site_ID==24) # n = 43
wilcox.test(filter(site24,HZtransect=="contemporary")$standHI,filter(site24,HZtransect=="historic")$standHI)
#data:  filter(site24, HZtransect == "contemporary")$standHI and filter(site24, HZtransect == "historic")$standHI
#W = 8, p-value = 0.2952

# site 25
site25 <- filter(scoring_df,site_ID==25) # n = 31
wilcox.test(filter(site25,HZtransect=="contemporary")$standHI,filter(site25,HZtransect=="historic")$standHI)
#data:  filter(site25, HZtransect == "contemporary")$standHI and filter(site25, HZtransect == "historic")$standHI
#W = 28, p-value = 0.9648

# site 26
site26 <- filter(scoring_df,site_ID==26) # n = 40
wilcox.test(filter(site26,HZtransect=="contemporary")$standHI,filter(site26,HZtransect=="historic")$standHI)
#data:  filter(site26, HZtransect == "contemporary")$standHI and filter(site26, HZtransect == "historic")$standHI
#W = 46, p-value = 0.1656

# site 27
site27 <- filter(scoring_df,site_ID==27) # n = 5
wilcox.test(filter(site27,HZtransect=="contemporary")$standHI,filter(site27,HZtransect=="historic")$standHI)
#data:  filter(site27, HZtransect == "contemporary")$standHI and filter(site27, HZtransect == "historic")$standHI
#W = 2.5, p-value = 1
