library(glmmTMB)
library(DHARMa)
library(car)
library(emmeans)
library(ggplot2)
library(tidyverse)

# 3_juvenile_coral_density
# script analyzes and plots data on the density of juvenile corals

#------------------------------------------------------------------------------------------------------------#
#--------------------------------------------- Data ---------------------------------------------------------
#------------------------------------------------------------------------------------------------------------#

# data collected from time-series photoquadrats
# data from MCR LTER: Coral Reef: Early life stage bottleneck determines rates of coral recovery following severe disturbance; Data for Speare et al., 2024, Ecology
# https://doi.org/10.6073/pasta/f006b56623d3dd61d689237b048b53d8

juv_density<-read.csv("data/juvenile_coral_density.csv", header=TRUE, stringsAsFactors = TRUE)

#------------------------------------------------------------------------------------------------------------#
#---------------- Does the density of juvenile corals differ across depths?----------------------------------
#------------------------------------------------------------------------------------------------------------#

#glmm with negative binomial distribution
# including quadrat nested within site
juv_dens_glmm2 <- glmmTMB(recruits_m2 ~ Depth + YearFirstAppeared + Depth:YearFirstAppeared + (1|Site/Site_Depth_Pole_Quad), data=recruits_per_quad_year, family=nbinom2(link = "log"))

car::Anova(juv_dens_glmm2, type=2) ## model results with Wald chi squared tests
check_overdispersion(juv_dens_glmm2)
plot(simulateResiduals(fittedModel=juv_dens_glmm2, quantileFunction = qnorm))
residuals(simulateResiduals(fittedModel=juv_dens_glmm2, quantileFunction = qnorm))
plotResiduals(simulateResiduals(fittedModel=juv_dens_glmm2, quantileFunction = qnorm), form=recruits_per_quad_year$YearFirstAppeared)
plotResiduals(simulateResiduals(fittedModel=juv_dens_glmm2, quantileFunction = qnorm), form=recruits_per_quad_year$Depth)

# model predicted means
juv_dens_glmm_predicted<-data.frame(emmeans(juv_dens_glmm2, type="response", ~Depth*YearFirstAppeared))

# plotting it
juv_density_plot<-ggplot()+
  geom_jitter(data=subset(recruits_per_quad_year, YearFirstAppeared!="2010"), aes(x=YearFirstAppeared, y=recruits_m2, group=Depth, color=Depth), width=0.3, alpha=(0.5), shape=1, size=2)+
  geom_errorbar(data=subset(juv_dens_glmm_predicted, YearFirstAppeared!="2010"), aes(x=YearFirstAppeared, ymin=response-SE, ymax=response+SE, width=0), color="black")+
  geom_line(data=subset(juv_dens_glmm_predicted, YearFirstAppeared!="2010"), aes(x=YearFirstAppeared, y=response, color=Depth, group=Depth))+
  geom_point(data=subset(juv_dens_glmm_predicted, YearFirstAppeared!="2010"), aes(x=YearFirstAppeared, y=response, fill=Depth, group=Depth), pch=21, size=4)+
  scale_color_manual(values=c("#6DDED2", "#4A918A"), na.translate = F)+
  scale_fill_manual(values=c("#6DDED2", "#4A918A"), na.translate = F)+
  theme_classic()+
  ylab(expression(paste("New ", italic("Pocillopora"), " juveniles m"^-2)))+
  xlab("Year")+
  #scale_y_continuous(breaks=c(0,2,4,6,8))+
  theme(axis.text.x = element_text(colour="black", angle=-45, size=14, hjust=0), 
        axis.text.y = element_text(colour="black", size=16))+
  theme(axis.title.x=element_blank(), axis.title.y=element_text(size=16))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.ticks = element_line(color="black"))+
  theme(legend.title = element_blank(), legend.text=element_text(size=16))+
  theme(legend.position="none")+
  theme(aspect.ratio = 4/4)

