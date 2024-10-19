
library(glmmTMB)
library(Matrix)
library(DHARMa)
library(emmeans)
library(ggplot2)
library(car)
library(effects)

# 5_juvenile_coral_growth
# script analyzes and plots data on the growth of juvenile corals
# combines all plots for figure 2

#------------------------------------------------------------------------------------------------------------#
#--------------------------------------------- Data ---------------------------------------------------------
#------------------------------------------------------------------------------------------------------------#

# data collected from time-series photoquadrats
# data from MCR LTER: Coral Reef: Early life stage bottleneck determines rates of coral recovery following severe disturbance; Data for Speare et al., 2024, Ecology
# https://doi.org/10.6073/pasta/f006b56623d3dd61d689237b048b53d8

juv_growth<-read.csv("data_clean/juvenile_coral_growth.csv", header=TRUE, stringsAsFactors = TRUE)

#------------------------------------------------------------------------------------------------------------#
#--------------------- Does the growth of juvenile corals differ across depths? ------------------------------
#------------------------------------------------------------------------------------------------------------#

juv_growth$Site_Depth_Pole_Quad<-as.factor(paste(juv_growth$Site, juv_growth$Depth, juv_growth$Pole, juv_growth$Quad, sep="_"))
juv_growth$year1<-as.factor(juv_growth$year1)

juv_growth$log_ch_area<-log(juv_growth$ch_area) # logging for future use. add a number because there are no zeros
juv_growth$log1_ch_area<-log(1+(juv_growth$ch_area))
# you want to use the log+1 data because otherwise you end up with negative values

#---------------------------------- gaussian. log transformed (+1) 

juv_growth_glmm1 <- glmmTMB(log1_ch_area ~ Depth*year1 + (1|Site/Site_Depth_Pole_Quad), 
                            data=juv_growth, family=gaussian(link = "identity"))

car::Anova(juv_growth_glmm1, type=2) ## model results with Wald chi squared tests
summary(juv_growth_glmm1)
plot(allEffects(juv_growth_glmm1))
plot(simulateResiduals(juv_growth_glmm1, quantileFunction = qnorm))

juv_growth_glmm1_predicted<-data.frame(emmeans(juv_growth_glmm1, type="response", ~Depth*year1))

# Plot it 

juv_growth_plot<-ggplot()+
  geom_jitter(data=juv_growth, aes(x=year1, y=log1_ch_area, group=Depth, color=Depth), width=0.3, alpha=(0.4), shape=1, size=2)+
  geom_errorbar(data=juv_growth_glmm1_predicted, aes(x=year1, ymin=emmean-SE, ymax=emmean+SE, width=0) , color="black")+
  geom_line(data=juv_growth_glmm1_predicted, aes(x=year1, y=emmean, group=Depth, color=Depth))+
  geom_point(data=juv_growth_glmm1_predicted, aes(x=year1, y=emmean, group=Depth, fill=Depth), colour="black",pch=21, size=4)+
  scale_color_manual(values=c("#6DDED2", "#4A918A"), na.translate = F)+
  scale_fill_manual(values=c("#6DDED2", "#4A918A"), na.translate = F)+
  theme_classic()+
  ylab(expression ("log (change in area + 1) ("~cm^2 ~yr^-1~")"))+
  xlab("Year")+
  scale_x_discrete(labels=c('2011-2012', '2012-2013', '2013-2014', '2014-2015', '2015-2016'))+
  theme(axis.text.x = element_text(colour="black", angle=-45, size=14, hjust=0), 
        axis.text.y = element_text(colour="black", size=16))+
  theme(axis.title.x=element_text(size=16), axis.title.y=element_text(size=16))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.ticks = element_line(color="black"))+
  theme(legend.title = element_blank(), legend.text=element_text(size=16))+
  theme(legend.position="none")+
  theme(aspect.ratio = 4/4)

#getting legend

figlegend<-ggplot()+
  geom_line(data=juv_growth_glmm1_predicted, aes(x=year1, y=emmean, group=Depth, color=Depth))+
  geom_errorbar(data=juv_growth_glmm1_predicted, aes(x=year1, ymin=emmean-SE, ymax=emmean+SE, width=0), color="black")+
  geom_point(data=juv_growth_glmm1_predicted, aes(x=year1, y=emmean, group=Depth, fill=Depth), colour="black", pch=21, size=4)+
  scale_color_manual(values=c("#6DDED2", "#4A918A"), na.translate = F)+
  scale_fill_manual(values=c("#6DDED2", "#4A918A"), na.translate = F)+
  theme_minimal()+theme(legend.title = element_blank(), legend.text =element_text(size=12) )
legend_fig2 <- get_legend(figlegend + theme(legend.box.margin = margin(0, 0, 0, 12)))

#---------------------------------- gaussian. log transformed (+1) response- including size in model

#including size in year1 in the model
juv_growth_glmm2 <- glmmTMB(log1_ch_area ~ Depth*year1 + diameter_year1 + (1|Site/Site_Depth_Pole_Quad), 
                            data=juv_growth, family=gaussian(link = "identity"))
car::Anova(juv_growth_glmm2, type=2) ## model results with Wald chi squared tests
summary(juv_growth_glmm2)
plot(simulateResiduals(juv_growth_glmm2, quantileFunction = qnorm))

# sometimes getting a warning message that appears to be a bug with matrix package 
#https://stackoverflow.com/questions/77466641/error-in-fittmbtmbstruc-negative-log-likelihood-is-nan-at-starting-parameter
#https://stackoverflow.com/questions/77533611/error-in-glmmtmb-after-matrix-update-negative-log-likelihood-is-nan-at-starting

juv_growth_glmm2_predicted<-data.frame(emmeans(juv_growth_glmm2, type="response", ~Depth*year1))

juv_growth_size_plot<-ggplot()+
  geom_jitter(data=juv_growth, aes(x=year1, y=log1_ch_area, group=Depth, color=Depth), width=0.3, alpha=(0.4), shape=1, size=2)+
  geom_errorbar(data=juv_growth_glmm2_predicted, aes(x=year1, ymin=emmean-SE, ymax=emmean+SE, width=0) , color="black")+
  geom_line(data=juv_growth_glmm2_predicted, aes(x=year1, y=emmean, group=Depth, color=Depth))+
  geom_point(data=juv_growth_glmm2_predicted, aes(x=year1, y=emmean, group=Depth, fill=Depth), colour="black",pch=21, size=4)+
  scale_color_manual(values=c("#6DDED2", "#4A918A"), na.translate = F)+
  scale_fill_manual(values=c("#6DDED2", "#4A918A"), na.translate = F)+
  theme_classic()+
  ylab(expression ("log (change in area + 1) ("~cm^2 ~yr^-1~")"))+
  xlab("Year")+
  scale_x_discrete(labels=c('2011-2012', '2012-2013', '2013-2014', '2014-2015', '2015-2016'))+
  theme(axis.text.x = element_text(colour="black", angle=-45, size=14, hjust=0), 
        axis.text.y = element_text(colour="black", size=16))+
  theme(axis.title.x=element_text(size=16), axis.title.y=element_text(size=16))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.ticks = element_line(color="black"))+
  theme(legend.title = element_blank(), legend.text=element_text(size=16))+
  #theme(legend.position="none")+
  theme(aspect.ratio = 4/4)
ggsave("figures/growth_size_supp.pdf", width=5, height=5, units="in")

#------------------------------------------------------------------------------------------------------------#
#------------------------------------- Combine plots for Figure 2  -------------------------------------------
#------------------------------------------------------------------------------------------------------------#

### -------------- combining all plots for figure 2

recruits_juvs_mortality_growth<-cowplot::plot_grid(recruit_tile_plot, juv_density_plot, juv_mortality_plot, juv_growth_plot,
                                                   nrow = 2, align = "vh", labels = c("(a)", "(b)", "(c)", "(d)") )
recruits_juvs_mortality_growth_leg<-plot_grid(recruits_juvs_mortality_growth, legend_fig2, nrow = 1, align = "vh", rel_widths = c(0.93,0.07))

ggsave("figures/recruits_juvs_mortality_growth_leg.pdf", width=11.5, height=10, units="in")


   