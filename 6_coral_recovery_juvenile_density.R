library(glmmTMB)
library(Matrix)
library(DHARMa)
library(emmeans)
library(ggplot2)
library(car)
library(effects)

# 6_coral_recovery_juvenile_density
# script analyzes relationship between density of juvenile corals and recovery rate of corals
# requires dataframes generated in script 1_time_series_coral_cover.R

#------------------------------------------------------------------------------------------------------------#
#--------------------------------------------- Data ---------------------------------------------------------
#------------------------------------------------------------------------------------------------------------#

# data collected from time-series photoquadrats
# data from MCR LTER: Coral Reef: Early life stage bottleneck determines rates of coral recovery following severe disturbance; Data for Speare et al., 2024, Ecology
# https://doi.org/10.6073/pasta/f006b56623d3dd61d689237b048b53d8

# coral cover for each quadrat
#reading coral cover data generated from script 1
coral_cover_quad_means<-read.csv("data/coral_cover_quad_means.csv", header = TRUE)

# density of juvenile corals
juvs_per_quad_year<-read.csv("data_clean/juvenile_coral_density.csv", header=TRUE, stringsAsFactors = TRUE)

#------------------------------------------------------------------------------------------------------------#
#--------------Relationship between the density of juvenile pocs and the rate of recovery---------------------
#------------------------------------------------------------------------------------------------------------#


# the peak years of new juvenile corals were from 2011-2013. summing the number of unique juvenile corals per quadrat for those years

juvs_per_quad_2011_2013<-ddply(subset(juvs_per_quad_year, YearFirstAppeared=="2011" | YearFirstAppeared=="2012" | YearFirstAppeared=="2013"),
                                   .(Site_Depth_Pole_Quad), summarize,
                                   count_juvs_2011_2013=sum(count_juvs),
                                   juvs_m2_2011_2013=sum(juvs_m2), 
                                   n_years_data=length(count_juvs))

#####---------------------------- Recovery of Pocillopora ----------------------------------------------

# nadir = minimum coral cover
nadir_pocillopora<-coral_cover_quad_means
nadir_pocillopora$Mean_Stony_Coral<-NULL
nadir_pocillopora<-subset(nadir_pocillopora, Year!="2005" & Year!="2006" & Year!="2007" & Year!="2008" & Year!="2009")
nadir_pocillopora<-nadir_pocillopora[nadir_pocillopora$Mean_Pocillopora == ave(nadir_pocillopora$Mean_Pocillopora, nadir_pocillopora$Site_Depth_Pole_Quad, FUN=min), ]

# some quadrats had the same minimum coral cover for 2 years in a row
# used the first year that a quadrat reached the minimum as the nadir
# I filtered them by year and used the first year that a quadrat reached minimum coral cover 
nadir_pocillopora %>% arrange(Mean_Pocillopora) %>% distinct(Site_Depth_Pole_Quad, .keep_all = TRUE)

nadir_pocillopora<-nadir_pocillopora %>% group_by(Site_Depth_Pole_Quad) %>% arrange(Year)%>% filter(row_number()==1)
nadir_pocillopora$Site_Depth_Pole_Quad<-as.factor(nadir_pocillopora$Site_Depth_Pole_Quad)
names(nadir_pocillopora)[names(nadir_pocillopora) == 'Mean_Pocillopora'] <- "Nadir_Pocillopora"
nadir_pocillopora$Year<-as.numeric(as.character(nadir_pocillopora$Year))
names(nadir_pocillopora)[names(nadir_pocillopora)=="Year"]<-"Year_nadir"

#then coral cover in 2018 (max coral cover)
pocillopora_2018<-subset(coral_cover_quad_means, Year=="2018")
names(pocillopora_2018)[names(pocillopora_2018)=="Mean_Pocillopora"]<-"Pocillopora_2018"
pocillopora_2018$Year<-as.numeric(as.character(pocillopora_2018$Year))
names(pocillopora_2018)[names(pocillopora_2018)=="Year"]<-"Year_2018"
pocillopora_2018<-pocillopora_2018[c(1,4:6)]

# merging dataframes for nadier (minimum) pocillopora, and pocillopora cover in 2018
recovery_pocillopora<-merge(nadir_pocillopora, pocillopora_2018, by="Site_Depth_Pole_Quad")

# number of years from nadir to 2018
recovery_pocillopora$Recovery_Years<-recovery_pocillopora$Year_2018 - recovery_pocillopora$Year_nadir
# recovery rate is the change in pocillopora cover from nadier to 2018, divided by the number of years
recovery_pocillopora$Recovery_Rate_Pocillopora<-(recovery_pocillopora$Pocillopora_2018 - recovery_pocillopora$Nadir_Pocillopora)/recovery_pocillopora$Recovery_Years

#####---------------------------- Rate of recovery vs Density of Juvneile corals  ----------------------------------------

#combining data on recovery rate of pocillopora and density of juvs in quadrats
recovery_pocillopora<-merge(juvs_per_quad_2011_2013, recovery_pocillopora, by="Site_Depth_Pole_Quad")
recovery_pocillopora$Site_Depth<-as.factor(paste(recovery_pocillopora$Site, recovery_pocillopora$Depth, sep="_"))
recovery_pocillopora$Site_Depth<-factor(recovery_pocillopora$Site_Depth, levels=c("LTER1_10m", "LTER2_10m","LTER1_17m", "LTER2_17m"))

# quick plot: recruit density v rate of recovery. diff colors by depth
ggplot(recovery_pocillopora, aes(x=juvs_m2_2011_2013, y=Recovery_Rate_Pocillopora, color=Depth))+
  geom_point(aes(color=Depth), size=2.5)+
  geom_smooth(method=lm, se=TRUE)+
  theme_classic()+
  scale_color_manual(values=c("#6DDED2", "#4A918A"))+
  labs(x = expression ("juvs"~m^-2~"(2011-2013)"),
       y = expression ("Rate of recovery (%"~yr^-1~")"))+
  theme(axis.text.x = element_text(colour="black", size=12), 
        axis.text.y = element_text(colour="black", size=14))+
  theme(axis.title.x=element_text(size=14), axis.title.y=element_text(size=16))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.ticks = element_line(color="black"))+
  theme(legend.title = element_blank(), legend.text=element_text(size=12))+
  theme(aspect.ratio = 4/4)


#####---------------------------- GLMM analysis ----------------------------------------------

## density of juveniles vs recovery 

### model that includes depth as a fixed effect
recov_glmm1<-glmmTMB(Recovery_Rate_Pocillopora ~  juvs_m2_2011_2013*Depth +(1|Site), 
                    data=recovery_pocillopora, family=gaussian(link = "identity"))
car::Anova(recov_glmm1, type=2) ## model results with Wald chi squared tests
summary(recov_glmm1)
plot(simulateResiduals(recov_glmm1, quantileFunction = qnorm))
plot(allEffects(recov_glmm1))


emtrends(recov_glmm1, pairwise ~ Depth, var = "juvs_m2_2011_2013")
# slopes of lines are not statistically different

### predicted recovery
recov_ref<-ref_grid(recov_glmm2, at=list(juvs_m2_2011_2013=c(0,5,10,20,45,60,80,95,100,110,120,130,140)))
recov_pred<-data.frame(emmip(recov_ref, Depth ~ juvs_m2_2011_2013, cov.reduce=range, plotit=FALSE)) #look at them visually
#subsetting the predicted data to only include the range of observed data for each depth
recov_pred10<-subset(recov_pred, Depth=="10m" & juvs_m2_2011_2013>40) #10m: 45-136
recov_pred17<-subset(recov_pred, Depth=="17m" & juvs_m2_2011_2013<100) #17m: 4:92
recov_pred_sub<-rbind(recov_pred10, recov_pred17)


# plot figure

ggplot()+
  geom_point(data=recovery_pocillopora, aes(x=juvs_m2_2011_2013, y=Recovery_Rate_Pocillopora, color=Depth), alpha=0.8, shape=1, size=2, stroke=1)+ #colored points are actual data
  geom_ribbon(data=recov_pred_sub, aes(x=juvs_m2_2011_2013, ymin=yvar-SE, ymax=yvar+SE, fill=Depth),alpha= 0.5) +
  geom_line(data=recov_pred_sub, aes(x=juvs_m2_2011_2013, y = yvar, color=Depth), size = 0.9)+
  theme_classic()+
  scale_y_continuous(limits=c(0,11.5), breaks=c(0,2,4,6,8,10))+
  scale_color_manual(values=c("#6DDED2", "#4A918A"))+
  scale_fill_manual(values = c("#6DDED2", "#4A918A"))+
  labs(x = expression ("Juvenile corals"~m^-2~"(2011-2013)"),
       y = expression ("Rate of recovery (%"~yr^-1~")"))+
  theme(axis.text.x = element_text(colour="black", size=12), 
        axis.text.y = element_text(colour="black", size=14))+
  theme(axis.title.x=element_text(size=14), axis.title.y=element_text(size=16))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.ticks = element_line(color="black"))+
  #theme(legend.position = "none")+
  theme(aspect.ratio = 4/4)


