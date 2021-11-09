# This script contains the R code to analyze photo data of contemporary samples.
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
library(gridExtra)
library(ghibli)

# load dataset
photo_df <- read_tsv("./raw-data/photo_scores_HZonly.txt",col_names=TRUE)


# NUCHAL PATCH
nuchal_df <- photo_df %>%
  filter(trait=="nuchal") %>%
  distinct()

summary(lm(nuchal ~ area,data=nuchal_df))
#Call:
#  lm(formula = nuchal ~ area, data = nuchal_df)
#
#Residuals:
#  Min       1Q   Median       3Q      Max 
#-2.63734 -0.60240 -0.07847  0.77006  2.45033 
#
#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)    
#(Intercept)  2.9112264  0.1810593   16.08   <2e-16 ***
#  area        -0.0002367  0.0000182  -13.00   <2e-16 ***
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#Residual standard error: 0.8795 on 88 degrees of freedom
#Multiple R-squared:  0.6577,	Adjusted R-squared:  0.6538 
#F-statistic: 169.1 on 1 and 88 DF,  p-value: < 2.2e-16

nuchal_plot <- ggplot(data=nuchal_df,aes(x=factor(nuchal),y=area)) +
  geom_boxplot(outlier.shape=NA,color="gray44") +
  coord_flip() +
  geom_point(color=ghibli_palettes$KikiMedium[3]) +
  geom_smooth(aes(group=1),color=ghibli_palettes$KikiMedium[3],method=lm) +
  xlab("Nuchal score") +
  theme_bw() +
  theme(axis.title=element_text(face="bold",size=10), axis.text=element_text(size=9,color="black"), 
        axis.title.x=element_blank())
nuchal_plot



# SHAFT
shaft_df <- photo_df %>%
  filter(trait=="vane")

summary(lm(shaft ~ lwMean + mwMean + swMean + uvMean + lumMean, data=shaft_df))
#Call:
#  lm(formula = shaft ~ lwMean + mwMean + swMean + uvMean + lumMean, 
#     data = shaft_df)
#
#Residuals:
#  Min       1Q   Median       3Q      Max 
#-0.81793 -0.21217 -0.03536  0.20575  1.92759 
#
#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)    
#(Intercept)    -3.201      1.759  -1.820  0.07232 .  
#lwMean        798.968    251.952   3.171  0.00212 ** 
#mwMean       1450.125    508.487   2.852  0.00547 ** 
#swMean        535.961    159.384   3.363  0.00116 ** 
#uvMean        -34.756      7.799  -4.457 2.55e-05 ***
#lumMean     -2725.386    908.219  -3.001  0.00354 ** 
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#Residual standard error: 0.4325 on 84 degrees of freedom
#Multiple R-squared:  0.915,	Adjusted R-squared:  0.9099 
#F-statistic: 180.7 on 5 and 84 DF,  p-value: < 2.2e-16

shaft_plot <- ggplot(data=shaft_df,aes(x=factor(shaft),y=lwMean+mwMean+swMean+uvMean+lumMean)) +
  geom_boxplot(outlier.shape=NA,color="gray44") +
  coord_flip() +
  geom_point(color=ghibli_palettes$MarnieMedium2[6]) +
  geom_smooth(aes(group=1),color=ghibli_palettes$MarnieMedium2[6],method=lm) +
  xlab("Shaft score") +
  ylab("Sum of photo values") +
  theme_bw() +
  theme(axis.title=element_text(face="bold",size=10), axis.text=element_text(size=9,color="black"))
shaft_plot



# THROAT
throat_df <- photo_df %>%
  filter(trait=="throat") %>%
  distinct()

summary(lm(throat ~ lwMean + mwMean + swMean + uvMean + lumMean, data=throat_df))
#Call:
#  lm(formula = throat ~ lwMean + mwMean + swMean + uvMean + lumMean, 
#     data = throat_df)
#
#Residuals:
#  Min      1Q  Median      3Q     Max 
#-2.0786 -0.3771 -0.0006  0.3505  2.2378 
#
#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)  
#(Intercept)     7.413      3.535   2.097   0.0390 *
#lwMean       -930.082    535.093  -1.738   0.0858 .
#mwMean      -1550.254   1096.899  -1.413   0.1613  
#swMean       -588.426    353.356  -1.665   0.0996 .
#uvMean         71.327     35.489   2.010   0.0477 *
#lumMean      2994.036   1941.640   1.542   0.1268  
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#Residual standard error: 0.6178 on 84 degrees of freedom
#Multiple R-squared:  0.8623,	Adjusted R-squared:  0.8541 
#F-statistic: 105.2 on 5 and 84 DF,  p-value: < 2.2e-16

throat_plot <- ggplot(data=throat_df,aes(x=factor(throat),y=lwMean+mwMean+swMean+uvMean+lumMean)) +
  geom_boxplot(outlier.shape=NA,color="gray44") +
  coord_flip() +
  geom_point(color=ghibli_palettes$PonyoMedium[3]) +
  geom_smooth(aes(group=1),color=ghibli_palettes$PonyoMedium[3],method=lm) +
  xlab("Throat score") +
  ylab("Sum of photo values") +
  theme_bw() +
  theme(axis.title=element_text(face="bold",size=10), axis.text=element_text(size=9,color="black"))
throat_plot



# CROWN
crown_df <- photo_df %>%
  filter(trait=="crown") %>%
  distinct()

summary(lm(crown ~ lwMean + mwMean + swMean + uvMean + lumMean, data=crown_df))
#Call:
#  lm(formula = crown ~ lwMean + mwMean + swMean + uvMean + lumMean, 
#     data = crown_df)
#
#Residuals:
#  Min      1Q  Median      3Q     Max 
#-3.1689 -0.7811 -0.2990  0.6205  3.2926 
#
#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)  
#(Intercept)    20.88      10.74   1.944   0.0553 .
#lwMean      -1448.00    1515.57  -0.955   0.3421  
#mwMean      -3935.41    3219.34  -1.222   0.2250  
#swMean      -1021.14     978.95  -1.043   0.2999  
#uvMean         34.08      91.56   0.372   0.7107  
#lumMean      6273.81    5591.50   1.122   0.2650  
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#Residual standard error: 1.175 on 84 degrees of freedom
#Multiple R-squared:  0.4318,	Adjusted R-squared:  0.398 
#F-statistic: 12.77 on 5 and 84 DF,  p-value: 3.109e-09

crown_plot <- ggplot(data=crown_df,aes(x=factor(crown),y=lwMean+mwMean+swMean+uvMean+lumMean)) +
  geom_boxplot(outlier.shape=NA,color="gray44") +
  coord_flip() +
  geom_point(color=ghibli_palettes$PonyoMedium[2]) +
  geom_smooth(aes(group=1),color=ghibli_palettes$PonyoMedium[2],method=lm) +
  xlab("Crown score") +
  theme_bw() +
  theme(axis.title=element_text(face="bold",size=10), axis.text=element_text(size=9,color="black"),
        axis.title.x=element_blank())
crown_plot



# EAR COVERTS
ear_df <- photo_df %>%
  filter(trait=="ear") %>%
  distinct()

summary(lm(ear_coverts ~ lwMean + mwMean + swMean + uvMean + lumMean, data=ear_df))
#Call:
#  lm(formula = ear_coverts ~ lwMean + mwMean + swMean + uvMean + 
#       lumMean, data = ear_df)
#
#Residuals:
#  Min       1Q   Median       3Q      Max 
#-2.17868 -0.47363  0.04957  0.39108  2.07186 
#
#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)
#(Intercept)    10.656      7.518   1.417    0.160
#lwMean      -1405.151   1162.717  -1.209    0.230
#mwMean      -2749.376   2376.036  -1.157    0.251
#swMean       -776.488    757.357  -1.025    0.308
#uvMean         40.741     65.097   0.626    0.533
#lumMean      4882.607   4213.413   1.159    0.250
#
#Residual standard error: 0.8873 on 84 degrees of freedom
#Multiple R-squared:  0.6991,	Adjusted R-squared:  0.6812 
#F-statistic: 39.04 on 5 and 84 DF,  p-value: < 2.2e-16

ear_plot <- ggplot(data=ear_df,aes(x=factor(ear_coverts),y=lwMean+mwMean+swMean+uvMean+lumMean)) +
  geom_boxplot(outlier.shape=NA,color="gray44") +
  coord_flip() +
  geom_point(color=ghibli_palettes$PonyoMedium[5]) +
  geom_smooth(aes(group=1),color=ghibli_palettes$PonyoMedium[5],method=lm) +
  xlab("Ear covert score") +
  theme_bw() +
  theme(axis.title=element_text(face="bold",size=10), axis.text=element_text(size=9,color="black"),
        axis.title.x=element_blank())
ear_plot



# MALE MALAR
malar_df <- photo_df %>%
  filter(trait=="malar") %>%
  filter(sex=="Male") %>%
  distinct()

summary(lm(malar ~ lwMean + mwMean + swMean + uvMean + lumMean, data=malar_df))
#Call:
#  lm(formula = malar ~ lwMean + mwMean + swMean + uvMean + lumMean, 
#     data = malar_df)
#
#Residuals:
#  Min      1Q  Median      3Q     Max 
#-1.0090 -0.4473 -0.1487  0.4507  1.4299 
#
#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)
#(Intercept)    10.625      9.194   1.156    0.253
#lwMean      -1418.970   1434.256  -0.989    0.327
#mwMean      -2958.755   2890.221  -1.024    0.311
#swMean       -929.032    933.326  -0.995    0.324
#uvMean         44.261     75.155   0.589    0.558
#lumMean      5228.492   5160.379   1.013    0.316
#
#Residual standard error: 0.5868 on 53 degrees of freedom
#Multiple R-squared:  0.8828,	Adjusted R-squared:  0.8718 
#F-statistic: 79.88 on 5 and 53 DF,  p-value: < 2.2e-16

malar_plot <- ggplot(data=malar_df,aes(x=factor(malar),y=lwMean+mwMean+swMean+uvMean+lumMean)) +
  geom_boxplot(outlier.shape=NA,color="gray44") +
  coord_flip() +
  geom_point(color=ghibli_palettes$PonyoMedium[6]) +
  geom_smooth(aes(group=1),color=ghibli_palettes$PonyoMedium[6],method=lm) +
  xlab("Malar stripe score") +
  theme_bw() +
  theme(axis.title=element_text(face="bold",size=10), axis.text=element_text(size=9,color="black"),
        axis.title.x=element_blank())
malar_plot



# arrange the plots together
grid.arrange(crown_plot, ear_plot, malar_plot, nuchal_plot, shaft_plot, throat_plot, nrow=3,heights=c(0.9,0.9,1))
# save plot as PDF with dimensions 6" x 7.75"
