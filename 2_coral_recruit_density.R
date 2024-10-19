library(ggplot2)
library(tidyverse)
library(cowplot)
library(ggpubr)
library(plyr)
library(dplyr)
library(splitstackshape)
library(tidyr)
library(stats)
library(lme4)
library(data.table)
library(effects)
library(lattice)
library(pbkrtest)
library(emmeans)
library(car)
library(lmerTest)
library(glmmTMB)
library(DHARMa)
library(performance)

# 2_coral_recruit_density
# script analyzes and plots data on the density of coral recruits on tiles

#-------------------------------------------------------------------------------------------------------#
### ---------------------------------------- Data -------------------------------------------------------
#-------------------------------------------------------------------------------------------------------#

# data on coral recruit density on coral settlement tiles
# data from MCR LTER: Coral Reef: Coral Community Dynamics: Coral Recruitment
# https://doi.org/10.6073/pasta/57b3d0d926dd643b471d261a5984f078

recs<-read.csv("data/coral_recruit_tile_spat_counts_2006-2016_20180626.csv", header=TRUE, stringsAsFactors = TRUE)

###----------------------------------- data wrangling  --------------------------------------------------

recs$nominal_year<-as.factor(recs$nominal_year)
recs$tile_id<-as.factor(recs$tile_id)
recs$season<-as.factor(recs$season)

names(recs)[names(recs) == 'nominal_year'] <- "year"

recs<-subset(recs, year!="2006" & year!="2007" & year!="2008" & year!="2009" & year!="2016") #subsetting out years before 2010 and after 2015
recs<-subset(recs, habitat=="Forereef") #only interested in the forereef habitat
recs<-subset(recs, family=="Pocilloporidae") #only considering Pocilloporidae coral recruits

recs$location<-factor(recs$location)
#renaming the locations
recs$location<-revalue(recs$location, c("LTER 1 Outer 10m coral Recruit Tiles"="LTER 1 10 m", "LTER 2 Outer 10m coral Recruit Tiles"="LTER 2 10 m",
                                        "LTER 1 Outer 17m coral Recruit Tiles"="LTER 1 17 m", "LTER 2 Outer 17m coral Recruit Tiles"="LTER 2 17 m"))
#creating a separate factor for depth
recs$depth<-ifelse(recs$location=="LTER 1 10 m" | recs$location=="LTER 2 10 m", "10 m", "17 m")
recs$depth<-as.factor(recs$depth)

#creating a separate factor for site
recs$site<-ifelse(recs$location=="LTER 1 10 m" | recs$location=="LTER 1 17 m", "LTER1", "LTER2")
recs$depth<-as.factor(recs$depth)

#------------------------------------ summarizing data ----------------------------------------------

#tiles are scored and replaced twice annually. Summing the two numbers for each year for total annual recruitment
recs_tile_season<-ddply(recs, .(year, tile_id, depth, site, season), summarize, poc_sum=sum(count))
recs_tile_ann<-ddply(recs_tile_season, .(year, tile_id, depth, site), summarize, poc_ann_sum=sum(poc_sum))

recs_site_ann_means<-ddply(recs_tile_ann, .(year, depth, site), summarize, poc_means=mean(poc_ann_sum))


recs_depth_ann_means<-ddply(recs_site_ann_means, .(year, depth), summarize, 
                            mean_poc=mean(poc_means),
                            se_poc=sd(poc_means)/sqrt(length(poc_means)))


#---------------------------------- GLMM analysis ---------------------------------------------

glmm1 <- glmer.nb(poc_ann_sum ~ depth*year + (1|site), data=recs_tile_ann, family="nbinom2")
Anova(glmm1)
car::Anova(glmm1, type=2) ## model results with Wald chi squared tests
summary(glmm1)
plot(allEffects(glmm1))
plot(predictorEffect("depth", mod=glmm1))
plot(simulateResiduals(glmm1, quantileFunction = qnorm))

glmm1_tukey<-emmeans(glmm1, type="response", ~depth*year)
glmm1_tukey
cld(glmm1_tukey, Letters=letters)

glmm1_tukey2<-emmeans(glmm1, type="response", ~depth)
cld(glmm1_tukey2, Letters=letters)

# pulling out model fitted values for each level depth, year
# making a unique dataframe to predict over
recruit_glmm_predicted<-data.frame(emmeans(glmm1, type="response", ~depth*year))

recruit_tile_plot<-ggplot()+
  geom_jitter(data=recs_tile_ann, aes(x=year, y=poc_ann_sum, group=depth, color=depth), width=0.3, alpha=(0.5), shape=1, size=2)+
  geom_errorbar(data=recruit_glmm_predicted, aes(x=year, ymin=response-SE, ymax=response+SE, width=0), color="black")+
  geom_line(data=recruit_glmm_predicted, aes(x=year, y=response, group=depth, color=depth))+
  geom_point(data=recruit_glmm_predicted, aes(x=year, y=response, group=depth, fill=depth), pch=21, size=4)+
  scale_color_manual(values=c("#6DDED2", "#4A918A"), na.translate = F)+
  scale_fill_manual(values=c("#6DDED2", "#4A918A"), na.translate = F)+
  theme_classic()+
  ylab(expression(paste("Pocilloporidae recruits ", tile^-1)))+
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


