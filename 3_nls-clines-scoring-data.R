# This script contains the R code for nls cline fitting.
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
library(nlstools)


# load datasets
contemporary <- read_tsv("scoring_contemporary_input.txt",col_names=TRUE)
historic <- read_tsv("scoring_historic_input.txt",col_names=TRUE)


### CONTEMPORARY ANALYSES

# set up parameters
# cline for HI
rhs <- function(x, c, w) {
  1/(1+exp((4*(x-c))/w))
}
# cline for individual traits
rhs2 <- function(x, c, w) {
  4/(1+exp((4*(x-c))/w))
}
# for plotting bootstrap
s <- seq(0,800,length=100)



### HYBRID INDEX
# modelling the cline
m_standHI_c <- nls(standHI_mean ~ rhs(dist, center, width),
                   data=contemporary,
                   start=list(center=390,width=390),
                   trace=T)

# summarizing the output
summary(m_standHI_c)
coef(m_standHI_c)

# calculating a confidence interval
CI_standHI_c <- confint(m_standHI_c,parm=c("center","width"))
CI_standHI_c

# bootstrapping the output
bootstrap_standHI_c <- nlsBoot(m_standHI_c,niter=999)
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

sink("./cline-output-scores/standHI_c_output.txt")
print('##### MODEL SUMMARY #####')
summary(m_standHI_c)
print('##### COEFF FROM MODEL #####')
coef(m_standHI_c)
print('##### CI FROM MODEL #####')
CI_standHI_c
print('##### BOOTSTRAP SUMMARY #####')
summary(bootstrap_standHI_c)
sink()

write.table(boot_df_standHI_c,"./cline-output-scores/bootstrap-CI_standHI_c.txt",row.names=FALSE,sep="\t",quote=FALSE)


### SHAFT
# modelling the cline
m_shaft_c <- nls(shaft_mean ~ rhs2(dist, center, width),
                 data=contemporary,
                 start=list(center=390,width=390),
                 trace=T)

# summarizing the output
summary(m_shaft_c)
coef(m_shaft_c)

# calculating a confidence interval
CI_shaft_c <- confint(m_shaft_c,parm=c("center","width"))
CI_shaft_c

# bootstrapping the output
bootstrap_shaft_c <- nlsBoot(m_shaft_c,niter=999)
summary(bootstrap_shaft_c)
boot_coeff_shaft_c <- bootstrap_shaft_c$coefboot

# bootstrapping for plotting
y_shaft_c <- matrix(NA,nrow=length(bootstrap_shaft_c),ncol=length(s))
for(i in 1:length(bootstrap_shaft_c)){
  y_shaft_c[i,]=rhs2(s,boot_coeff_shaft_c[i,1],boot_coeff_shaft_c[i,2])
}
boot_ci_shaft_c <- apply(y_shaft_c,2,quantile,c(0.025,0.975))
boot_df_shaft_c <- data.frame(s,t(boot_ci_shaft_c))
names(boot_df_shaft_c) = c("dist","lower","upper")

sink("./cline-output-scores/shaft_c_output.txt")
print('##### MODEL SUMMARY #####')
summary(m_shaft_c)
print('##### COEFF FROM MODEL #####')
coef(m_shaft_c)
print('##### CI FROM MODEL #####')
CI_shaft_c
print('##### BOOTSTRAP SUMMARY #####')
summary(bootstrap_shaft_c)
sink()

write.table(boot_df_shaft_c,"./cline-output-scores/bootstrap-CI_shaft_c.txt",row.names=FALSE,sep="\t",quote=FALSE)


### CROWN
# modelling the cline
m_crown_c <- nls(crown_mean ~ rhs2(dist, center, width),
                 data=contemporary,
                 start=list(center=390,width=390),
                 trace=T)

# summarizing the output
summary(m_crown_c)
coef(m_crown_c)

# calculating a confidence interval
CI_crown_c <- confint(m_crown_c,parm=c("center","width"))
CI_crown_c

# bootstrapping the output
bootstrap_crown_c <- nlsBoot(m_crown_c,niter=999)
summary(bootstrap_crown_c)
boot_coeff_crown_c <- bootstrap_crown_c$coefboot

# bootstrapping for plotting
y_crown_c <- matrix(NA,nrow=length(bootstrap_crown_c),ncol=length(s))
for(i in 1:length(bootstrap_crown_c)){
  y_crown_c[i,]=rhs2(s,boot_coeff_crown_c[i,1],boot_coeff_crown_c[i,2])
}
boot_ci_crown_c <- apply(y_crown_c,2,quantile,c(0.025,0.975))
boot_df_crown_c <- data.frame(s,t(boot_ci_crown_c))
names(boot_df_crown_c) = c("dist","lower","upper")

sink("./cline-output-scores/crown_c_output.txt")
print('##### MODEL SUMMARY #####')
summary(m_crown_c)
print('##### COEFF FROM MODEL #####')
coef(m_crown_c)
print('##### CI FROM MODEL #####')
CI_crown_c
print('##### BOOTSTRAP SUMMARY #####')
summary(bootstrap_crown_c)
sink()

write.table(boot_df_crown_c,"./cline-output-scores/bootstrap-CI_crown_c.txt",row.names=FALSE,sep="\t",quote=FALSE)


### THROAT
# modelling the cline
m_throat_c <- nls(throat_mean ~ rhs2(dist, center, width),
                  data=contemporary,
                  start=list(center=390,width=390),
                  trace=T)

# summarizing the output
summary(m_throat_c)
coef(m_throat_c)

# calculating a confidence interval
CI_throat_c <- confint(m_throat_c,parm=c("center","width"))
CI_throat_c

# bootstrapping the output
bootstrap_throat_c <- nlsBoot(m_throat_c,niter=999)
summary(bootstrap_throat_c)
boot_coeff_throat_c <- bootstrap_throat_c$coefboot

# bootstrapping for plotting
y_throat_c <- matrix(NA,nrow=length(bootstrap_throat_c),ncol=length(s))
for(i in 1:length(bootstrap_throat_c)){
  y_throat_c[i,]=rhs2(s,boot_coeff_throat_c[i,1],boot_coeff_throat_c[i,2])
}
boot_ci_throat_c <- apply(y_throat_c,2,quantile,c(0.025,0.975))
boot_df_throat_c <- data.frame(s,t(boot_ci_throat_c))
names(boot_df_throat_c) = c("dist","lower","upper")

sink("./cline-output-scores/throat_c_output.txt")
print('##### MODEL SUMMARY #####')
summary(m_throat_c)
print('##### COEFF FROM MODEL #####')
coef(m_throat_c)
print('##### CI FROM MODEL #####')
CI_throat_c
print('##### BOOTSTRAP SUMMARY #####')
summary(bootstrap_throat_c)
sink()

write.table(boot_df_throat_c,"./cline-output-scores/bootstrap-CI_throat_c.txt",row.names=FALSE,sep="\t",quote=FALSE)


### NUCHAL
# modelling the cline
m_nuchal_c <- nls(nuchal_mean ~ rhs2(dist, center, width),
                  data=contemporary,
                  start=list(center=390,width=390),
                  trace=T)

# summarizing the output
summary(m_nuchal_c)
coef(m_nuchal_c)

# calculating a confidence interval
CI_nuchal_c <- confint(m_nuchal_c,parm=c("center","width"))
CI_nuchal_c

# bootstrapping the output
bootstrap_nuchal_c <- nlsBoot(m_nuchal_c,niter=999)
summary(bootstrap_nuchal_c)
boot_coeff_nuchal_c <- bootstrap_nuchal_c$coefboot

# bootstrapping for plotting
y_nuchal_c <- matrix(NA,nrow=length(bootstrap_nuchal_c),ncol=length(s))
for(i in 1:length(bootstrap_nuchal_c)){
  y_nuchal_c[i,]=rhs2(s,boot_coeff_nuchal_c[i,1],boot_coeff_nuchal_c[i,2])
}
boot_ci_nuchal_c <- apply(y_nuchal_c,2,quantile,c(0.025,0.975))
boot_df_nuchal_c <- data.frame(s,t(boot_ci_nuchal_c))
names(boot_df_nuchal_c) = c("dist","lower","upper")

sink("./cline-output-scores/nuchal_c_output.txt")
print('##### MODEL SUMMARY #####')
summary(m_nuchal_c)
print('##### COEFF FROM MODEL #####')
coef(m_nuchal_c)
print('##### CI FROM MODEL #####')
CI_nuchal_c
print('##### BOOTSTRAP SUMMARY #####')
summary(bootstrap_nuchal_c)
sink()

write.table(boot_df_nuchal_c,"./cline-output-scores/bootstrap-CI_nuchal_c.txt",row.names=FALSE,sep="\t",quote=FALSE)


### EAR COVERTS
# modelling the cline
m_ears_c <- nls(ears_mean ~ rhs2(dist, center, width),
                data=contemporary,
                start=list(center=390,width=390),
                trace=T)

# summarizing the output
summary(m_ears_c)
coef(m_ears_c)

# calculating a confidence interval
CI_ears_c <- confint(m_ears_c,parm=c("center","width"))
CI_ears_c

# bootstrapping the output
bootstrap_ears_c <- nlsBoot(m_ears_c,niter=999)
summary(bootstrap_ears_c)
boot_coeff_ears_c <- bootstrap_ears_c$coefboot

# bootstrapping for plotting
y_ears_c <- matrix(NA,nrow=length(bootstrap_ears_c),ncol=length(s))
for(i in 1:length(bootstrap_ears_c)){
  y_ears_c[i,]=rhs2(s,boot_coeff_ears_c[i,1],boot_coeff_ears_c[i,2])
}
boot_ci_ears_c <- apply(y_ears_c,2,quantile,c(0.025,0.975))
boot_df_ears_c <- data.frame(s,t(boot_ci_ears_c))
names(boot_df_ears_c) = c("dist","lower","upper")

sink("./cline-output-scores/ears_c_output.txt")
print('##### MODEL SUMMARY #####')
summary(m_ears_c)
print('##### COEFF FROM MODEL #####')
coef(m_ears_c)
print('##### CI FROM MODEL #####')
CI_ears_c
print('##### BOOTSTRAP SUMMARY #####')
summary(bootstrap_ears_c)
sink()

write.table(boot_df_ears_c,"./cline-output-scores/bootstrap-CI_ears_c.txt",row.names=FALSE,sep="\t",quote=FALSE)


### MALAR
# modelling the cline
m_malar_c <- nls(malar_mean ~ rhs2(dist, center, width),
                 data=contemporary,
                 start=list(center=390,width=390),
                 trace=T)

# summarizing the output
summary(m_malar_c)
coef(m_malar_c)

# calculating a confidence interval
CI_malar_c <- confint(m_malar_c,parm=c("center","width"))
CI_malar_c

# bootstrapping the output
bootstrap_malar_c <- nlsBoot(m_malar_c,niter=999)
summary(bootstrap_malar_c)
boot_coeff_malar_c <- bootstrap_malar_c$coefboot

# bootstrapping for plotting
y_malar_c <- matrix(NA,nrow=length(bootstrap_malar_c),ncol=length(s))
for(i in 1:length(bootstrap_malar_c)){
  y_malar_c[i,]=rhs2(s,boot_coeff_malar_c[i,1],boot_coeff_malar_c[i,2])
}
boot_ci_malar_c <- apply(y_malar_c,2,quantile,c(0.025,0.975))
boot_df_malar_c <- data.frame(s,t(boot_ci_malar_c))
names(boot_df_malar_c) = c("dist","lower","upper")

sink("./cline-output-scores/malar_c_output.txt")
print('##### MODEL SUMMARY #####')
summary(m_malar_c)
print('##### COEFF FROM MODEL #####')
coef(m_malar_c)
print('##### CI FROM MODEL #####')
CI_malar_c
print('##### BOOTSTRAP SUMMARY #####')
summary(bootstrap_malar_c)
sink()

write.table(boot_df_malar_c,"./cline-output-scores/bootstrap-CI_malar_c.txt",row.names=FALSE,sep="\t",quote=FALSE)







### HISTORIC ANALYSES

# set up some parameters
# note: these are the same parameters created above
# cline for HI
rhs <- function(x, c, w) {
  1/(1+exp((4*(x-c))/w))
}
# cline for individual traits
rhs2 <- function(x, c, w) {
  4/(1+exp((4*(x-c))/w))
}
# for plotting bootstrap
s <- seq(0,800,length=100)



### HYBRID INDEX
# modelling the cline
m_standHI_h <- nls(standHI_mean ~ rhs(dist, center, width),
                   data=historic,
                   start=list(center=390,width=390),
                   trace=T)

# summarizing the output
summary(m_standHI_h)
coef(m_standHI_h)

# calculating a confidence interval
CI_standHI_h <- confint(m_standHI_h,parm=c("center","width"))
CI_standHI_h

# bootstrapping the output
bootstrap_standHI_h <- nlsBoot(m_standHI_h,niter=999)
summary(bootstrap_standHI_h)
boot_coeff_standHI_h <- bootstrap_standHI_h$coefboot

# bootstrapping for plotting
y_standHI_h <- matrix(NA,nrow=length(bootstrap_standHI_h),ncol=length(s))
for(i in 1:length(bootstrap_standHI_h)){
  y_standHI_h[i,]=rhs(s,boot_coeff_standHI_h[i,1],boot_coeff_standHI_h[i,2])
}
boot_hi_standHI_h <- apply(y_standHI_h,2,quantile,c(0.025,0.975))
boot_df_standHI_h <- data.frame(s,t(boot_hi_standHI_h))
names(boot_df_standHI_h) = c("dist","lower","upper")

sink("./cline-output-scores/standHI_h_output.txt")
print('##### MODEL SUMMARY #####')
summary(m_standHI_h)
print('##### COEFF FROM MODEL #####')
coef(m_standHI_h)
print('##### CI FROM MODEL #####')
CI_standHI_h
print('##### BOOTSTRAP SUMMARY #####')
summary(bootstrap_standHI_h)
sink()

write.table(boot_df_standHI_h,"./cline-output-scores/bootstrap-CI_standHI_h.txt",row.names=FALSE,sep="\t",quote=FALSE)


### SHAFT
# modelling the cline
m_shaft_h <- nls(shaft_mean ~ rhs2(dist, center, width),
                 data=historic,
                 start=list(center=390,width=390),
                 trace=T)

# summarizing the output
summary(m_shaft_h)
coef(m_shaft_h)

# calculating a confidence interval
CI_shaft_h <- confint(m_shaft_h,parm=c("center","width"))
CI_shaft_h

# bootstrapping the output
bootstrap_shaft_h <- nlsBoot(m_shaft_h,niter=999)
summary(bootstrap_shaft_h)
boot_coeff_shaft_h <- bootstrap_shaft_h$coefboot

# bootstrapping for plotting
y_shaft_h <- matrix(NA,nrow=length(bootstrap_shaft_h),ncol=length(s))
for(i in 1:length(bootstrap_shaft_h)){
  y_shaft_h[i,]=rhs2(s,boot_coeff_shaft_h[i,1],boot_coeff_shaft_h[i,2])
}
boot_hi_shaft_h <- apply(y_shaft_h,2,quantile,c(0.025,0.975))
boot_df_shaft_h <- data.frame(s,t(boot_hi_shaft_h))
names(boot_df_shaft_h) = c("dist","lower","upper")

sink("./cline-output-scores/shaft_h_output.txt")
print('##### MODEL SUMMARY #####')
summary(m_shaft_h)
print('##### COEFF FROM MODEL #####')
coef(m_shaft_h)
print('##### CI FROM MODEL #####')
CI_shaft_h
print('##### BOOTSTRAP SUMMARY #####')
summary(bootstrap_shaft_h)
sink()

write.table(boot_df_shaft_h,"./cline-output-scores/bootstrap-CI_shaft_h.txt",row.names=FALSE,sep="\t",quote=FALSE)


### CROWN
# modelling the cline
m_crown_h <- nls(crown_mean ~ rhs2(dist, center, width),
                 data=historic,
                 start=list(center=390,width=390),
                 trace=T)

# summarizing the output
summary(m_crown_h)
coef(m_crown_h)

# calculating a confidence interval
CI_crown_h <- confint(m_crown_h,parm=c("center","width"))
CI_crown_h

# bootstrapping the output
bootstrap_crown_h <- nlsBoot(m_crown_h,niter=999)
summary(bootstrap_crown_h)
boot_coeff_crown_h <- bootstrap_crown_h$coefboot

# bootstrapping for plotting
y_crown_h <- matrix(NA,nrow=length(bootstrap_crown_h),ncol=length(s))
for(i in 1:length(bootstrap_crown_h)){
  y_crown_h[i,]=rhs2(s,boot_coeff_crown_h[i,1],boot_coeff_crown_h[i,2])
}
boot_hi_crown_h <- apply(y_crown_h,2,quantile,c(0.025,0.975))
boot_df_crown_h <- data.frame(s,t(boot_hi_crown_h))
names(boot_df_crown_h) = c("dist","lower","upper")

sink("./cline-output-scores/crown_h_output.txt")
print('##### MODEL SUMMARY #####')
summary(m_crown_h)
print('##### COEFF FROM MODEL #####')
coef(m_crown_h)
print('##### CI FROM MODEL #####')
CI_crown_h
print('##### BOOTSTRAP SUMMARY #####')
summary(bootstrap_crown_h)
sink()

write.table(boot_df_crown_h,"./cline-output-scores/bootstrap-CI_crown_h.txt",row.names=FALSE,sep="\t",quote=FALSE)


### THROAT
# modelling the cline
m_throat_h <- nls(throat_mean ~ rhs2(dist, center, width),
                  data=historic,
                  start=list(center=390,width=390),
                  trace=T)

# summarizing the output
summary(m_throat_h)
coef(m_throat_h)

# calculating a confidence interval
CI_throat_h <- confint(m_throat_h,parm=c("center","width"))
CI_throat_h

# bootstrapping the output
bootstrap_throat_h <- nlsBoot(m_throat_h,niter=999)
summary(bootstrap_throat_h)
boot_coeff_throat_h <- bootstrap_throat_h$coefboot

# bootstrapping for plotting
y_throat_h <- matrix(NA,nrow=length(bootstrap_throat_h),ncol=length(s))
for(i in 1:length(bootstrap_throat_h)){
  y_throat_h[i,]=rhs2(s,boot_coeff_throat_h[i,1],boot_coeff_throat_h[i,2])
}
boot_hi_throat_h <- apply(y_throat_h,2,quantile,c(0.025,0.975))
boot_df_throat_h <- data.frame(s,t(boot_hi_throat_h))
names(boot_df_throat_h) = c("dist","lower","upper")

sink("./cline-output-scores/throat_h_output.txt")
print('##### MODEL SUMMARY #####')
summary(m_throat_h)
print('##### COEFF FROM MODEL #####')
coef(m_throat_h)
print('##### CI FROM MODEL #####')
CI_throat_h
print('##### BOOTSTRAP SUMMARY #####')
summary(bootstrap_throat_h)
sink()

write.table(boot_df_throat_h,"./cline-output-scores/bootstrap-CI_throat_h.txt",row.names=FALSE,sep="\t",quote=FALSE)


### NUCHAL
# modelling the cline
m_nuchal_h <- nls(nuchal_mean ~ rhs2(dist, center, width),
                  data=historic,
                  start=list(center=390,width=390),
                  trace=T)

# summarizing the output
summary(m_nuchal_h)
coef(m_nuchal_h)

# calculating a confidence interval
CI_nuchal_h <- confint(m_nuchal_h,parm=c("center","width"))
CI_nuchal_h

# bootstrapping the output
bootstrap_nuchal_h <- nlsBoot(m_nuchal_h,niter=999)
summary(bootstrap_nuchal_h)
boot_coeff_nuchal_h <- bootstrap_nuchal_h$coefboot

# bootstrapping for plotting
y_nuchal_h <- matrix(NA,nrow=length(bootstrap_nuchal_h),ncol=length(s))
for(i in 1:length(bootstrap_nuchal_h)){
  y_nuchal_h[i,]=rhs2(s,boot_coeff_nuchal_h[i,1],boot_coeff_nuchal_h[i,2])
}
boot_hi_nuchal_h <- apply(y_nuchal_h,2,quantile,c(0.025,0.975))
boot_df_nuchal_h <- data.frame(s,t(boot_hi_nuchal_h))
names(boot_df_nuchal_h) = c("dist","lower","upper")

sink("./cline-output-scores/nuchal_h_output.txt")
print('##### MODEL SUMMARY #####')
summary(m_nuchal_h)
print('##### COEFF FROM MODEL #####')
coef(m_nuchal_h)
print('##### CI FROM MODEL #####')
CI_nuchal_h
print('##### BOOTSTRAP SUMMARY #####')
summary(bootstrap_nuchal_h)
sink()

write.table(boot_df_nuchal_h,"./cline-output-scores/bootstrap-CI_nuchal_h.txt",row.names=FALSE,sep="\t",quote=FALSE)


### EAR COVERTS
# modelling the cline
m_ears_h <- nls(ears_mean ~ rhs2(dist, center, width),
                data=historic,
                start=list(center=390,width=390),
                trace=T)

# summarizing the output
summary(m_ears_h)
coef(m_ears_h)

# calculating a confidence interval
CI_ears_h <- confint(m_ears_h,parm=c("center","width"))
CI_ears_h

# bootstrapping the output
bootstrap_ears_h <- nlsBoot(m_ears_h,niter=999)
summary(bootstrap_ears_h)
boot_coeff_ears_h <- bootstrap_ears_h$coefboot

# bootstrapping for plotting
y_ears_h <- matrix(NA,nrow=length(bootstrap_ears_h),ncol=length(s))
for(i in 1:length(bootstrap_ears_h)){
  y_ears_h[i,]=rhs2(s,boot_coeff_ears_h[i,1],boot_coeff_ears_h[i,2])
}
boot_hi_ears_h <- apply(y_ears_h,2,quantile,c(0.025,0.975))
boot_df_ears_h <- data.frame(s,t(boot_hi_ears_h))
names(boot_df_ears_h) = c("dist","lower","upper")

sink("./cline-output-scores/ears_h_output.txt")
print('##### MODEL SUMMARY #####')
summary(m_ears_h)
print('##### COEFF FROM MODEL #####')
coef(m_ears_h)
print('##### CI FROM MODEL #####')
CI_ears_h
print('##### BOOTSTRAP SUMMARY #####')
summary(bootstrap_ears_h)
sink()

write.table(boot_df_ears_h,"./cline-output-scores/bootstrap-CI_ears_h.txt",row.names=FALSE,sep="\t",quote=FALSE)


### MALAR
# modelling the cline
m_malar_h <- nls(malar_mean ~ rhs2(dist, center, width),
                 data=historic,
                 start=list(center=390,width=390),
                 trace=T)

# summarizing the output
summary(m_malar_h)
coef(m_malar_h)

# calculating a confidence interval
CI_malar_h <- confint(m_malar_h,parm=c("center","width"))
CI_malar_h

# bootstrapping the output
bootstrap_malar_h <- nlsBoot(m_malar_h,niter=999)
summary(bootstrap_malar_h)
boot_coeff_malar_h <- bootstrap_malar_h$coefboot

# bootstrapping for plotting
y_malar_h <- matrix(NA,nrow=length(bootstrap_malar_h),ncol=length(s))
for(i in 1:length(bootstrap_malar_h)){
  y_malar_h[i,]=rhs2(s,boot_coeff_malar_h[i,1],boot_coeff_malar_h[i,2])
}
boot_hi_malar_h <- apply(y_malar_h,2,quantile,c(0.025,0.975))
boot_df_malar_h <- data.frame(s,t(boot_hi_malar_h))
names(boot_df_malar_h) = c("dist","lower","upper")

sink("./cline-output-scores/malar_h_output.txt")
print('##### MODEL SUMMARY #####')
summary(m_malar_h)
print('##### COEFF FROM MODEL #####')
coef(m_malar_h)
print('##### CI FROM MODEL #####')
CI_malar_h
print('##### BOOTSTRAP SUMMARY #####')
summary(bootstrap_malar_h)
sink()

write.table(boot_df_malar_h,"./cline-output-scores/bootstrap-CI_malar_h.txt",row.names=FALSE,sep="\t",quote=FALSE)









### TESTING INFLUENCE OF INCLUDING SAMPLES ON THE NORTH PLATTE

# load datasets
contemporary <- read_tsv("scoring_contemporary_input.txt",col_names=TRUE)

# filter dataset to remove samples on the North Platte
S_contemporary <- contemporary %>%
  filter(!site_ID %in% c(6,9,11,13,16,17))


# set up parameters
# cline for HI
rhs <- function(x, c, w) {
  1/(1+exp((4*(x-c))/w))
}

# for plotting bootstrap
s <- seq(0,800,length=100)



### HYBRID INDEX
# modelling the cline
m_standHI_SPlatte_c <- nls(standHI_mean ~ rhs(dist, center, width),
                   data=S_contemporary,
                   start=list(center=390,width=390),
                   trace=T)

# summarizing the output
summary(m_standHI_SPlatte_c)
coef(m_standHI_SPlatte_c)

# calculating a confidence interval
CI_standHI_SPlatte_c <- confint(m_standHI_SPlatte_c,parm=c("center","width"))
CI_standHI_SPlatte_c

# bootstrapping the output
bootstrap_standHI_SPlatte_c <- nlsBoot(m_standHI_SPlatte_c,niter=999)
summary(bootstrap_standHI_SPlatte_c)
boot_coeff_standHI_SPlatte_c <- bootstrap_standHI_SPlatte_c$coefboot

# bootstrapping for plotting
y_standHI_SPlatte_c <- matrix(NA,nrow=length(bootstrap_standHI_SPlatte_c),ncol=length(s))
for(i in 1:length(bootstrap_standHI_SPlatte_c)){
  y_standHI_SPlatte_c[i,]=rhs(s,boot_coeff_standHI_SPlatte_c[i,1],boot_coeff_standHI_SPlatte_c[i,2])
}
boot_ci_standHI_SPlatte_c <- apply(y_standHI_SPlatte_c,2,quantile,c(0.025,0.975))
boot_df_standHI_SPlatte_c <- data.frame(s,t(boot_ci_standHI_SPlatte_c))
names(boot_df_standHI_SPlatte_c) = c("dist","lower","upper")

sink("./cline-output-scores/standHI_SPlatte_c_output.txt")
print('##### MODEL SUMMARY #####')
summary(m_standHI_SPlatte_c)
print('##### COEFF FROM MODEL #####')
coef(m_standHI_SPlatte_c)
print('##### CI FROM MODEL #####')
CI_standHI_SPlatte_c
print('##### BOOTSTRAP SUMMARY #####')
summary(bootstrap_standHI_SPlatte_c)
sink()

write.table(boot_df_standHI_SPlatte_c,"./cline-output-scores/bootstrap-CI_standHI_SPlatte_c.txt",row.names=FALSE,sep="\t",quote=FALSE)
