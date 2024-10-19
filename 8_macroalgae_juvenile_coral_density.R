library(glmmTMB)
library(DHARMa)
library(car)
library(emmeans)
library(ggplot2)
library(tidyverse)

# 8_macroalgae_juvenile_coral_density
# script analyzes data the abundance of macroalgae and relationship to juvenile coral density

#------------------------------------------------------------------------------------------------------------#
#--------------------------------------------- Data ---------------------------------------------------------
#------------------------------------------------------------------------------------------------------------#

# requires benthic_cover dataframe summarized in script "1_time_series_coral_cover.R"

## benthic cover data from photoquad timeseries in script 1_time_series_coral_cover.R
benthic_cover<-read.csv("data/benthic_cover.csv", header=TRUE, stringsAsFactors = TRUE)

#------------------------------------------------------------------------------------------------------------#
#--------------------------------------- Percent cover of Macroalgae  ---------------------------------------#
#------------------------------------------------------------------------------------------------------------#
## percent cover of macroalgae
macro<-benthic_cover[c(1,2,37:40,16)]

macro$Site<-revalue(macro$Site, c("LTER 1"="LTER1", "LTER 2"="LTER2"))

macro$Site_Depth_Pole_Quad<-as.factor(paste(macro$Site, macro$Depth, macro$Transect, macro$Quadrat, sep="_"))

macro<-subset(macro, Site=="LTER1" | Site=="LTER2")

macro<-subset(macro, Year!="2005" & Year!="2006" & Year!="2007" & Year!="2008" & Year!="2009")


macro_means<-ddply(macro, .(Site, Depth, Year), summarize,
                   mean=mean(Macroalgae), 
                   SE=sd(Macroalgae)/sqrt(length(Macroalgae)))

macro.cover.mod<-lmer(Macroalgae ~ Depth + Year + Depth:Year+ (1|Site), data=macro)
plot(macro.cover.mod)

macro.cover.mod.1<-lmer(Macroalgae ~ Depth + Year + (1|Site), data=macro)
anova(macro.cover.mod, macro.cover.mod.1) #significant interaction of depth and year

macro_glmm<-glmmTMB(Macroalgae ~ Depth + Year + (1|Site), data=macro, family=gaussian(link = "identity"))
car::Anova(macro_glmm, type=2) ## model results with Wald chi squared tests
summary(macro_glmm)
plot(simulateResiduals(macro_glmm, quantileFunction = qnorm))
plot(allEffects(macro_glmm))

predicted_macro_cover<-data.frame(emmeans(macro.cover.mod, ~ Depth*Year))

macroalgae_plot<-ggplot()+
  geom_jitter(data=macro, aes(x=Year, y=Macroalgae, group=Depth, color=Depth), width=0.3, alpha=(0.3), shape=1, size=2)+ 
  geom_point(data=predicted_macro_cover, aes(x=Year, y=emmean, group=Depth, color=Depth),shape=16, size=3)+
  geom_errorbar(data=predicted_macro_cover, aes(x=Year,ymin=emmean-SE, ymax=emmean+SE, width=0, color=Depth))+
  geom_line(data=predicted_macro_cover, aes(x=Year, y=emmean, group=Depth, color=Depth))+
  scale_color_manual(values=c("#6DDED2", "#4A918A"), na.translate = F)+
  theme_classic()+
  xlab("Year")+
  ylab("% Cover of fleshy macroalgae")+
  #ylim(c(0,100))+
  labs(color="")+
  theme(axis.text.x = element_text(colour="black", size=12, hjust=-0.2, angle=-45), 
        axis.text.y = element_text(colour="black", size=12))+
  theme(axis.title.x=element_text(size=16), axis.title.y=element_text(size=16))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.ticks = element_line(color="black"))+
  theme(legend.title=element_blank(), legend.text = element_text(size=16))+
  theme(strip.text.x = element_text(size = 20),strip.background = element_blank())+
  theme(aspect.ratio = 4/4)
ggsave("figures/macroalgae_plot.pdf", width=5,height=5)



#------------------------------------------------------------------------------------------------------------#
#------------- Relationship between macroalgae cover and the density of juvenile pocillopora ----------------
#------------------------------------------------------------------------------------------------------------#

# summarizing abundance of macroalgae in each quadrat
macro_cover_quad_means<-ddply(macro, .(Site, Depth, Year, Site_Depth_Pole_Quad), summarize,
                   mean_macroalgae=mean(Macroalgae))

#subsetting years of interest 
macro_sub<-subset(macro, Year=="2011" |Year=="2012"|Year=="2013")
# calculating mean macroalgae cover across those 3 years
macro_cover_quad_means<-ddply(macro_sub, .(Site, Depth,Site_Depth_Pole_Quad), summarize, mean_macro=mean(Macroalgae))

macro_juvs<-left_join(recruits_per_quad_2011_2013, macro_cover_quad_means, by=c("Site_Depth_Pole_Quad"))
macro_juvs$log_mean_macro<-log(macro_juvs$mean_macro + 1)

#adding factors for site and depth 
macro_juvs$Site_Depth_Pole_Quad2<-macro_juvs$Site_Depth_Pole_Quad
macro_juvs<-macro_juvs%>% separate(Site_Depth_Pole_Quad2, c("Site", "Depth", "Pole", "Quad"), convert=TRUE)
macro_juvs$Site<-as.factor(macro_juvs$Site)
macro_juvs$Depth<-as.factor(macro_juvs$Depth)


## glmm 
macro_juvs_glmm<-glmmTMB(recruits_m2_2011_2013 ~ log_mean_macro*Depth + (1|Site), 
                         data=macro_juvs, family=gaussian(link = "identity"))
car::Anova(macro_juvs_glmm, type=2) ## model results with Wald chi squared tests
summary(macro_juvs_glmm)
plot(simulateResiduals(macro_juvs_glmm, quantileFunction = qnorm))
plot(allEffects(recov_glmm2))


