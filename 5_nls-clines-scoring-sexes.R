# This script contains the R code to compare male and female clines
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
library(nlstools)
library(gridExtra)

# load full datasets (sexes combined)
contemporary <- read_tsv("scoring_contemporary_input.txt",col_names=TRUE)
historic <- read_tsv("scoring_historic_input.txt",col_names=TRUE)

# load female only datasets
contemporary_female <- read_tsv("scoring_contemporary_input_females.txt",col_names=TRUE)
historic_female <- read_tsv("scoring_historic_input_females.txt",col_names=TRUE)

# load male only datasets
contemporary_male <- read_tsv("scoring_contemporary_input_males.txt",col_names=TRUE)
historic_male <- read_tsv("scoring_historic_input_males.txt",col_names=TRUE)




# set up some parameters
# cline for HI
rhs <- function(x, c, w) {
  1/(1+exp((4*(x-c))/w))
}
# for plotting bootstrap
s <- seq(0,800,length=100)



### MALES, CONTEMPORARY, HI
# modelling the cline
males_standHI_c <- nls(standHI_mean ~ rhs(dist, center, width), 
                   data=contemporary_male,
                   start=list(center=390,width=390),
                   trace=T)

# summarizing the output
summary(males_standHI_c)
coef(males_standHI_c)

# calculating a confidence interval
CI_standHI_c <- confint(males_standHI_c,parm=c("center","width"))
CI_standHI_c

# bootstrapping the output
bootstrap_standHI_c <- nlsBoot(males_standHI_c,niter=999)
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

sink("./cline-output-scores/standHI_c_malesonly_output.txt")
print('##### MODEL SUMMARY #####')
summary(males_standHI_c)
print('##### COEFF FROM MODEL #####')
coef(males_standHI_c)
print('##### CI FROM MODEL #####')
CI_standHI_c
print('##### BOOTSTRAP SUMMARY #####')
summary(bootstrap_standHI_c)
sink()

write.table(boot_df_standHI_c,"./cline-output-scores/bootstrap-CI_standHI_malesonly_c.txt",row.names=FALSE,sep="\t",quote=FALSE)



### FEMALES, CONTEMPORARY, HI
# modelling the cline
females_standHI_c <- nls(standHI_mean ~ rhs(dist, center, width), 
                       data=contemporary_female,
                       start=list(center=390,width=390),
                       trace=T)

# summarizing the output
summary(females_standHI_c)
coef(females_standHI_c)

# calculating a confidence interval
CI_standHI_c <- confint(females_standHI_c,parm=c("center","width"))
CI_standHI_c

# bootstrapping the output
bootstrap_standHI_c <- nlsBoot(females_standHI_c,niter=999)
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

sink("./cline-output-scores/standHI_c_femalesonly_output.txt")
print('##### MODEL SUMMARY #####')
summary(females_standHI_c)
print('##### COEFF FROM MODEL #####')
coef(females_standHI_c)
print('##### CI FROM MODEL #####')
CI_standHI_c
print('##### BOOTSTRAP SUMMARY #####')
summary(bootstrap_standHI_c)
sink()

write.table(boot_df_standHI_c,"./cline-output-scores/bootstrap-CI_standHI_femalesonly_c.txt",row.names=FALSE,sep="\t",quote=FALSE)




### MALES, HISTORIC, HI
# modelling the cline
males_standHI_h <- nls(standHI_mean ~ rhs(dist, center, width), 
                       data=historic_male,
                       start=list(center=390,width=390),
                       trace=T)

# summarizing the output
summary(males_standHI_h)
coef(males_standHI_h)

# calculating a confidence interval
CI_standHI_h <- confint(males_standHI_h,parm=c("center","width"))
CI_standHI_h

# bootstrapping the output
bootstrap_standHI_h <- nlsBoot(males_standHI_h,niter=999)
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

sink("./cline-output-scores/standHI_h_malesonly_output.txt")
print('##### MODEL SUMMARY #####')
summary(males_standHI_h)
print('##### COEFF FROM MODEL #####')
coef(males_standHI_h)
print('##### CI FROM MODEL #####')
CI_standHI_h
print('##### BOOTSTRAP SUMMARY #####')
summary(bootstrap_standHI_h)
sink()

write.table(boot_df_standHI_h,"./cline-output-scores/bootstrap-CI_standHI_malesonly_h.txt",row.names=FALSE,sep="\t",quote=FALSE)


### FEMALES, HISTORIC, HI
# modelling the cline
females_standHI_h <- nls(standHI_mean ~ rhs(dist, center, width), 
                         data=historic_female,
                         start=list(center=390,width=390),
                         trace=T)

# summarizing the output
summary(females_standHI_h)
coef(females_standHI_h)

# calculating a confidence interval
CI_standHI_h <- confint(females_standHI_h,parm=c("center","width"))
CI_standHI_h

# bootstrapping the output
bootstrap_standHI_h <- nlsBoot(females_standHI_h,niter=999)
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

sink("./cline-output-scores/standHI_h_femalesonly_output.txt")
print('##### MODEL SUMMARY #####')
summary(females_standHI_h)
print('##### COEFF FROM MODEL #####')
coef(females_standHI_h)
print('##### CI FROM MODEL #####')
CI_standHI_h
print('##### BOOTSTRAP SUMMARY #####')
summary(bootstrap_standHI_h)
sink()

write.table(boot_df_standHI_h,"./cline-output-scores/bootstrap-CI_standHI_femalesonly_h.txt",row.names=FALSE,sep="\t",quote=FALSE)






##### PLOTTING THE CLINES

# load bootstrap output original from modelling
bootstrap_c <- read_tsv("./cline-output-scores/bootstrap-CI_standHI_c.txt",col_names=TRUE)
bootstrap_h <- read_tsv("./cline-output-scores/bootstrap-CI_standHI_h.txt",col_names=TRUE)

# load bootstrap output from single sex modelling
boot_male_c <- read_tsv("./cline-output-scores/bootstrap-CI_standHI_malesonly_c.txt",col_names=TRUE)
boot_male_h <- read_tsv("./cline-output-scores/bootstrap-CI_standHI_malesonly_h.txt",col_names=TRUE)
boot_female_c <- read_tsv("./cline-output-scores/bootstrap-CI_standHI_femalesonly_c.txt",col_names=TRUE)
boot_female_h <- read_tsv("./cline-output-scores/bootstrap-CI_standHI_femalesonly_h.txt",col_names=TRUE)


#### HYBRID INDEX, HISTORIC

historic_sex <- ggplot() +
  # both sexes
  geom_ribbon(data=bootstrap_h,aes(x=dist,ymin=lower,ymax=upper),alpha=0.15,fill="gray44") +
  geom_smooth(data=historic,aes(x=dist,y=standHI_mean),method="nls",formula=y~1/(1+exp((4*(x-c))/w)),se=FALSE,method.args=list(start=list(c=360,w=360)),color="gray55",size=2) +

  # females, red
  geom_ribbon(data=boot_female_h,aes(x=dist,ymin=lower,ymax=upper),alpha=0.15,fill="gray44") +
  geom_smooth(data=historic_female,aes(x=dist,y=standHI_mean),method="nls",formula=y~1/(1+exp((4*(x-c))/w)),se=FALSE,method.args=list(start=list(c=360,w=360)),color="#e86e61",size=1,linetype="longdash") +

  # males, blue
  geom_ribbon(data=boot_male_h,aes(x=dist,ymin=lower,ymax=upper),alpha=0.15,fill="gray44") +
  geom_smooth(data=historic_male,aes(x=dist,y=standHI_mean),method="nls",formula=y~1/(1+exp((4*(x-c))/w)),se=FALSE,method.args=list(start=list(c=360,w=360)),color="#0d72d6",size=1,linetype="longdash") +

  ylim(c(-0.01,1.01)) +
  xlim(c(-10,810)) +
  xlab("Distance (km)") +
  ylab("Hybrid index") +
  theme_bw() +
  theme(axis.title=element_text(face="bold",size=12), axis.text=element_text(size=10,color="black"))


#### HYBRID INDEX, CONTEMPORARY

contemporary_sex <- ggplot() +
  # both sexes
  geom_ribbon(data=bootstrap_c,aes(x=dist,ymin=lower,ymax=upper),alpha=0.15,fill="gray44") +
  geom_smooth(data=contemporary,aes(x=dist,y=standHI_mean),method="nls",formula=y~1/(1+exp((4*(x-c))/w)),se=FALSE,method.args=list(start=list(c=360,w=360)),color="black",size=2) +

  # females, red
  geom_ribbon(data=boot_female_c,aes(x=dist,ymin=lower,ymax=upper),alpha=0.15,fill="gray44") +
  geom_smooth(data=contemporary_female,aes(x=dist,y=standHI_mean),method="nls",formula=y~1/(1+exp((4*(x-c))/w)),se=FALSE,method.args=list(start=list(c=360,w=360)),color="#DC3220",size=1,linetype="longdash") +

  # males, blue
  geom_ribbon(data=boot_male_c,aes(x=dist,ymin=lower,ymax=upper),alpha=0.15,fill="gray44") +
  geom_smooth(data=contemporary_male,aes(x=dist,y=standHI_mean),method="nls",formula=y~1/(1+exp((4*(x-c))/w)),se=FALSE,method.args=list(start=list(c=360,w=360)),color="#0A58A6",size=1,linetype="longdash") +

  ylim(c(-0.01,1.01)) +
  xlim(c(-10,810)) +
  xlab("Distance (km)") +
  ylab("Hybrid index") +
  theme_bw() +
  theme(axis.title=element_text(face="bold",size=12), axis.text=element_text(size=10,color="black"))


grid.arrange(historic_sex,contemporary_sex,nrow=2)
# save plot as PDF with dimensions 5" x 8"
