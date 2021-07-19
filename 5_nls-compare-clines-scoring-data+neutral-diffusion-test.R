# This script contains the R code to compare historic and contemporary
# HI clines with AIC and test for neutral diffusion.
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


# load datasets
contemporary <- read_tsv("scoring_contemporary_input.txt",col_names=TRUE)
historic <- read_tsv("scoring_historic_input.txt",col_names=TRUE)

# set up some parameters
# cline for HI
rhs <- function(x, c, w) {
  1/(1+exp((4*(x-c))/w))
}



# REGULAR MODEL FOR THE CONTEMPORARY HI
m_standHI_c <- nls(standHI_mean ~ rhs(dist, center, width),
                   data=contemporary,
                   start=list(center=390,width=390),
                   trace=T)

# summarizing the output
coef(m_standHI_c)
#center    width
#135.6074 274.5806


# REGULAR MODEL FOR THE HISTORIC HI (TO GET PARAMETERS)
m_standHI_h <- nls(standHI_mean ~ rhs(dist, center, width),
                   data=historic,
                   start=list(center=390,width=390),
                   trace=T)

# summarizing the output
coef(m_standHI_h)
#center    width
#208.3763 251.1938


# CONTEMPORARY CLINE WITH THE CENTER FIXED TO THE HISTORIC LOCATION
m_standHI_fixed_c <- nls(standHI_mean ~ rhs(dist, 208.3763, width),
                         data=contemporary,
                         start=list(width=390),
                         trace=T)

# summarizing the output
coef(m_standHI_fixed_c)
#width
#366.8351



### testing AIC between the models

AIC(m_standHI_c,m_standHI_fixed_c)
#df        AIC
#m_standHI_c        3 -32.650649
#m_standHI_fixed_c  2  -6.597632

summary(AIC(m_standHI_c,m_standHI_fixed_c))
#df            AIC
#Min.   :2.00   Min.   :-32.651
#1st Qu.:2.25   1st Qu.:-26.137
#Median :2.50   Median :-19.624
#Mean   :2.50   Mean   :-19.624
#3rd Qu.:2.75   3rd Qu.:-13.111
#Max.   :3.00   Max.   : -6.598



# plot output to compare the regular cline (black) with the center fixed to the historic location cline (gray, dashed)
fixed_cline <- ggplot() +
  geom_pointrange(data=contemporary,aes(x=dist,y=standHI_mean,ymin=standHI_mean-standHI_se,ymax=standHI_mean+standHI_se),color="black",shape=16) +
  geom_smooth(data=contemporary,aes(x=dist,y=standHI_mean),method="nls",formula=y~1/(1+exp((4*(x-208.3763))/w)),se=FALSE,method.args=list(start=list(w=360)),color="gray55",linetype="dashed") +
  geom_smooth(data=contemporary,aes(x=dist,y=standHI_mean),method="nls",formula=y~1/(1+exp((4*(x-c))/w)),se=FALSE,method.args=list(start=list(c=360,w=360)),color="black") +
  ylim(c(-0.01,1.01)) +
  xlim(c(0,800)) +
  xlab("Distance (km)") +
  ylab("Hybrid index") +
  theme_bw() +
  theme(axis.title=element_text(face="bold",size=12), axis.text=element_text(size=10,color="black"))






#### TEST NEUTRAL DIFFUSION ####

# for plotting bootstrap
s <- seq(0,800,length=100)

# bootstrapping the output from the contemporary nls cline
bootstrap_standHI_c <- nlsBoot(m_standHI_c,niter=100000)
boot_coeff_standHI_c <- as.data.frame(bootstrap_standHI_c$coefboot)
bootstrap_contemporary <- boot_coeff_standHI_c %>%
  rename(width_c = width)

# bootstrapping the output from the historic nls cline
bootstrap_standHI_h <- nlsBoot(m_standHI_h,niter=100000)
boot_coeff_standHI_h <- as.data.frame(bootstrap_standHI_h$coefboot)
bootstrap_historic <- boot_coeff_standHI_h %>%
  rename(width_h = width)


### WANG ET. AL 2019 ANALYSIS (EQUATION 3)

# differences in width of observed data (lhs of equation 3)
width_df <- as.data.frame(cbind(bootstrap_contemporary$width_c,bootstrap_historic$width_h))
width_df <- width_df %>%
  rename(width_c = V1, width_h = V2)

width_df <- width_df %>%
  mutate(sq_width_c = width_c^2) %>%
  mutate(sq_width_h = width_h^2) %>%
  mutate(diff = sq_width_c - sq_width_h)

# 95% confidence interval of the output
quantile(width_df$diff,c(0.025,0.975))
#2.5%     97.5%
#  -37666.95  90968.12

#write.table(width_df,"./neutral-diffusion-test_bootstrap-hist_063021.txt",row.names=FALSE,sep="\t",quote=FALSE)


# estimation of neutral expectation (rhs of equation3)
# 2*pi*10^2*12

# possible sigma values
sigma <- 100.7   #RMS calculated in Moore and Buchanan 1985 from banding data
sigma2 <- 15*2   #KL Wiebe estimates natal dispersal to be "typically greater than 15 km, and probably much greater" (BotW)

# possible delta t values
#1955-1957 = historic
#2016-2018 = contemporary
delta_t <- 60/1.0   #62 years between time points divided by 1.0 years/generation
delta_t2 <- 60/1.8   #using same estimate as Wang et al. (2019); originally from Mila, Smith, & Wayne (2007)

# calculating the equation
# realistic
neutral <- 2*pi*(sigma)^2*(delta_t)   #3822875
# very conservative
neutral2 <- 2*pi*(sigma2)^2*(delta_t2)   #188495.6




# plot output
plot1 <- ggplot(data=width_df) +
  geom_histogram(aes(x=diff),fill="black",bins=500) +
  geom_vline(xintercept=neutral,color="red") +
  xlab(expression(w[c]^2-w[h]^2)) +
  ylab("Count") +
  theme_bw() +
  theme(axis.title=element_text(face="bold",size=12), axis.text=element_text(size=10,color="black"))

plot2 <- ggplot(data=width_df) +
  geom_histogram(aes(x=diff),fill="black",bins=500) +
  geom_vline(xintercept=neutral2,color="red") +
  #xlab(expression(w[c]^2-w[h]^2)) +
  #ylab("Count") +
  scale_x_continuous(breaks=c(0e+00, 2e+05, 4e+05)) +
  theme_bw() +
  theme(axis.title=element_text(face="bold",size=10), axis.text=element_text(size=8,color="black"), axis.title.x=element_blank(), axis.title.y=element_blank())

inset_plot <- plot1 + annotation_custom(ggplotGrob(plot2), xmin=1e+06, xmax=3.5e+06, ymin=4000, ymax=10500)







#### arrange the plots together
grid.arrange(fixed_cline, inset_plot, nrow=2)
# save plot as PDF with dimensions 8" x 4"
