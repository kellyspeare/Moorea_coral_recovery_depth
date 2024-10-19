library(glmmTMB)
library(DHARMa)
library(car)
library(emmeans)
library(ggplot2)
library(tidyverse)

# 4_juvenile_coral_mortality
# script analyzes and plots data on the mortality of juvenile corals

#------------------------------------------------------------------------------------------------------------#
#--------------------------------------------- Data ---------------------------------------------------------
#------------------------------------------------------------------------------------------------------------#

# data collected from time-series photoquadrats
# data from MCR LTER: Coral Reef: Early life stage bottleneck determines rates of coral recovery following severe disturbance; Data for Speare et al., 2024, Ecology
# https://doi.org/10.6073/pasta/f006b56623d3dd61d689237b048b53d8

juvs_mort_first_year<-read.csv("data_clean/juvenile_coral_mortality.csv", header=TRUE, stringsAsFactors = TRUE)

#------------------------------------------------------------------------------------------------------------#
#--------------- Does the mortality of juvenile corals differ across depths?----------------------------------
#------------------------------------------------------------------------------------------------------------#

juv_mort_glmm <- glmmTMB(Dead_year_2 ~ Depth + timeframe + Depth:timeframe + (1|Site/Site_Depth_Pole_Quad), 
                         data=juvs_mort_first_year, family=binomial(link = "logit"))

Anova(juv_mort_glmm)
car::Anova(juv_mort_glmm, type=2) ## model results with Wald chi squared tests
summary(juv_mort_glmm)
plot(allEffects(juv_mort_glmm))
check_overdispersion(juv_mort_glmm)
plot(simulateResiduals(juv_mort_glmm, quantileFunction = qnorm))

juv_mort_glmm_predicted<-data.frame(emmeans(juv_mort_glmm, type="response", ~Depth*timeframe))

juv_mortality_plot<-ggplot()+
  geom_jitter(data=juvs_mort_first_year, aes(x=timeframe, y=Dead_year_2, group=Depth, color=Depth), width=0.3, height=0.05, alpha=(0.5), shape=1, size=2)+
  geom_errorbar(data=juv_mort_glmm_predicted, aes(x=timeframe, ymin=prob-SE, ymax=prob+SE, width=0), color="black")+
  geom_line(data=juv_mort_glmm_predicted, aes(x=timeframe, y=prob, group=Depth, color=Depth))+
  geom_point(data=juv_mort_glmm_predicted, aes(x=timeframe, y=prob, group=Depth, fill=Depth), pch=21, size=4)+
  scale_color_manual(values=c("#6DDED2", "#4A918A"), na.translate = F)+
  scale_fill_manual(values=c("#6DDED2", "#4A918A"), na.translate = F)+
  theme_classic()+
  ylab(expression(paste("Prob. of mortaltity for juvenile ", italic("Pocillopora"))))+
  xlab("Year")+
  #scale_y_continuous(breaks=c(0,2,4,6,8))+
  theme(axis.text.x = element_text(colour="black", angle=-45, size=14, hjust=0), 
        axis.text.y = element_text(colour="black", size=16))+
  theme(axis.title.x=element_text(size=16), axis.title.y=element_text(size=16))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.ticks = element_line(color="black"))+
  theme(legend.title = element_blank(), legend.text=element_text(size=16))+
  theme(legend.position="none")+
  theme(aspect.ratio = 4/4)
                                       