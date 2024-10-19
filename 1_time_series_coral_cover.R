library(ggplot2)
library(cowplot)
library(ggpubr)
library(tidyverse)
library(plyr)
library(dplyr)
library(splitstackshape)
library(tidyr)
library(stats)
library(lme4)
library(emmeans)

# 1_time_series_coral_cover
# plots patterns of coral decline and recovery at 10 and 17m at 6 MCR LTER sites

# data on coral population and community dynamics from time series photoquadrats 
# data from MCR LTER: Coral Reef: Long-term Population and Community Dynamics: Corals, ongoing since 2005
# https://doi.org/10.6073/pasta/15d5120fb4f7b79811b16287eae15a35

#------------------------------------------------------------------------------------------------------------#
#--------------------------------------------------- Data ----------------------------------------------------
#------------------------------------------------------------------------------------------------------------#

## photoquadrat metadata
photoquad_meta<-read.csv("data/knb-lter-mcr.4_1_20191119.csv", header=TRUE, stringsAsFactors = TRUE)
photoquad_meta<-photoquad_meta[c(2:6)]
photoquad_meta$Location<-as.factor(photoquad_meta$Location)
photoquad_meta<-distinct(photoquad_meta, Location, .keep_all=TRUE)


## benthic cover data from photoquad timeseries
benthic_cover<-read.csv("data/knb-lter-mcr.4_2_20190321.csv", header=TRUE, stringsAsFactors = TRUE)

benthic_cover<-left_join(benthic_cover, photoquad_meta, by="Location") #merging with metadata

benthic_cover$Habitat<-as.factor(benthic_cover$Habitat)                    
benthic_cover<-subset(benthic_cover, Habitat=="Outer 10" | Habitat=="Outer 17")
benthic_cover$Habitat<-revalue(benthic_cover$Habitat, c("Outer 10"="10m", "Outer 17"="17m"))
names(benthic_cover)[names(benthic_cover)=="Habitat"]<-"Depth"

names(benthic_cover)[names(benthic_cover)=="Date"]<-"Year"
benthic_cover$Year<-as.factor(benthic_cover$Year)
benthic_cover$Year<-revalue(benthic_cover$Year, c("2005-05"="2005","2006-04"="2006","2007-04"="2007",
                                              "2008-04"="2008", "2009-04"="2009", "2010-04"="2010", "2011-04"="2011",
                                              "2012-04"="2012", "2013-04"="2013", "2014-04"="2014", "2015-04"="2015",
                                              "2016-04"="2016", "2017-04"="2017", "2018-04"="2018"))
# save for use in script 8
# write.csv(benthic_cover, "data/benthic_cover.csv", row.names=FALSE)

# making a dataframe for coral cover. including columns for Acropora, Pocillopora, Porites, and total coral cover
coral_cover<-benthic_cover[c(1,2,37:40, 4,23,24:27,36)]
coral_cover$Location<-as.factor(coral_cover$Location)

#combining all the porites species
coral_cover$Porites_sum<-coral_cover$Porites + coral_cover$Porites_irregularis + coral_cover$Porites_rus + coral_cover$Porites_spp_Massive
coral_cover$Porites_irregularis<-NULL
coral_cover$Porites_rus<-NULL
coral_cover$Porites_spp_Massive<-NULL
coral_cover$Porites<-NULL

names(coral_cover)[names(coral_cover)=="Porites_sum"]<-"Porites"

#subsetting out the df with all sites for plotting later
coral_cover_all_sites<-coral_cover

coral_cover_all_sites_means<-ddply(coral_cover_all_sites, .(Site, Depth, Year), summarize,
                                 mean_stony_coral=mean(Stony_Coral), 
                                 se_stony_coral=sd(Stony_Coral)/sqrt(length(Stony_Coral)))

#------------------------------------------------------------------------------------------------------------#
#-------------------------- Patterns of coral cover at LTER 1 and LTER 2 -------------------------------------
#------------------------------------------------------------------------------------------------------------#

# subesetting to only look at LTER 1 and 2
coral_cover<-subset(coral_cover, Site=="LTER 1" | Site=="LTER 2")
coral_cover$Site<-revalue(coral_cover$Site, c("LTER 1"="LTER1", "LTER 2"="LTER2"))

coral_cover$Site_Depth_Pole_Quad<-as.factor(paste(coral_cover$Site, coral_cover$Depth, coral_cover$Transect, coral_cover$Quadrat, sep="_"))

#setting asside this coral cover dataframe for plotting
coral_cover_for_plotting<-coral_cover

coral_cover_quad_means<-ddply(coral_cover, .(Site_Depth_Pole_Quad, Depth, Site, Year), summarize,
                         Mean_Stony_Coral=mean(Stony_Coral),
                         Mean_Pocillopora=mean(Pocillopora))

# need this dataframe for script "6_coral_recovery_juv_density.R" 
#saving this dataframe as a csv file
#write.csv(coral_cover_quad_means, "coral_cover_quad_means.csv", row.names = FALSE)

                        
## ------------ mixed effects model 

coral.cover.mod<-lmer(Mean_Stony_Coral ~ Depth + Year + Depth:Year+ (1|Site), data=coral_cover_quad_means)

predicted_coral_cover<-data.frame(emmeans(coral.cover.mod, ~ Depth*Year))

stony_coral_fig<-ggplot(predicted_coral_cover, aes(x=Year, y=emmean, color=Depth, group=Depth))+
  geom_point(shape=16, size=4)+
  geom_line()+
  geom_errorbar(aes(ymin=emmean-SE, ymax=emmean+SE, width=0))+
  scale_color_manual(values=c("#6DDED2", "#4A918A"), na.translate = F)+
  theme_classic()+
  ylim(c(-3,70))+
  ylab("% Cover of stony corals")+
  xlab("Year")+
  theme(axis.text.x = element_text(colour="black", angle=-45, size=14, hjust=0), 
        axis.text.y = element_text(colour="black", size=16))+
  theme(axis.title.x=element_text(size=16), axis.title.y=element_text(size=16))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.ticks = element_line(color="black"))+
  theme(legend.title = element_blank(), legend.text=element_text(size=16))+
  theme(aspect.ratio = 4/4)


#------------------------------------------------------------------------------------------------------------#
#-------------------------------- Absolute abundance of coral genera -----------------------------------------
#------------------------------------------------------------------------------------------------------------#

coral_cover_l<-coral_cover_for_plotting %>% pivot_longer(cols=c(7:10),
                                                         names_to = "Genus",
                                                         values_to = "Percent")
coral_cover_l$Genus<-as.factor(coral_cover_l$Genus)
coral_cover_l$Site<-factor(coral_cover_l$Site)

coral_cover_l_means<- ddply(coral_cover_l, .(Year, Depth, Genus), summarize,
                            Mean_Percent=mean(Percent))

ggplot(coral_cover_l_means, aes(x=Year, y=Mean_Percent, color=Genus, group=Genus))+
  geom_point(shape=16, size=4)+
  geom_line()+
 #geom_errorbar(aes(ymin=emmean-SE, ymax=emmean+SE, width=0))+
  #scale_color_manual(values=c("#6DDED2", "#4A918A"), na.translate = F)+
  facet_wrap(~Depth)+
  theme_classic()+
  #ylim(c(-3,70))+
  ylab("% Cover of stony corals")+
  xlab("Year")+
  theme(axis.text.x = element_text(colour="black", angle=-45, size=14, hjust=0), 
        axis.text.y = element_text(colour="black", size=16))+
  theme(axis.title.x=element_text(size=16), axis.title.y=element_text(size=16))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.ticks = element_line(color="black"))+
  theme(legend.title = element_blank(), legend.text=element_text(size=16))+
  theme(aspect.ratio = 4/4)  

#------------------------------------------------------------------------------------------------------------#
#-------------------------------- Relative abundance of coral genera -----------------------------------------
#------------------------------------------------------------------------------------------------------------#

coral_rel_abundance<-coral_cover_for_plotting
#calculating the relative abundance of Poc, Acr, and Porites
# doing the division with an ifelse() function so that it returns 0 instead of NaN when percent cover of corals is 0
coral_rel_abundance$Acropora_rel<-ifelse(coral_rel_abundance$Acropora==0, 0, coral_rel_abundance$Acropora / coral_rel_abundance$Stony_Coral)*100
coral_rel_abundance$Pocillopora_rel<-ifelse(coral_rel_abundance$Pocillopora==0, 0, coral_rel_abundance$Pocillopora / coral_rel_abundance$Stony_Coral)*100
coral_rel_abundance$Porites_rel<-ifelse(coral_rel_abundance$Porites==0, 0, coral_rel_abundance$Porites / coral_rel_abundance$Stony_Coral)*100

names(coral_rel_abundance)[names(coral_rel_abundance)=="Habitat"]<-"Depth"

coral_rel_abundance_means<-ddply(coral_rel_abundance, .(Depth, Site, Year), summarize,
             Mean_Pocillopora_rel=mean(Pocillopora_rel),
             Mean_Acropora_rel=mean(Acropora_rel),
             Mean_Porites_rel=mean(Porites_rel))

coral_rel_abundance_depth_means<-ddply(coral_rel_abundance_means, .(Depth, Year), summarize,
                                       Pocillopora=mean(Mean_Pocillopora_rel),
                                       SE_Poc=sd(Mean_Pocillopora_rel)/sqrt(length(Mean_Pocillopora_rel)),
                                       Acropora=mean(Mean_Acropora_rel),
                                       SE_Acr=sd(Mean_Acropora_rel)/sqrt(length(Mean_Acropora_rel)),
                                       Porites=mean(Mean_Porites_rel),
                                       SE_Por=sd(Mean_Porites_rel)/sqrt(length(Mean_Porites_rel)),)

coral_rel_abundance_depth_means$Other<-100-(coral_rel_abundance_depth_means$Pocillopora + coral_rel_abundance_depth_means$Acropora + coral_rel_abundance_depth_means$Porites)

coral_rel_abundance_depth_means_l<-gather(coral_rel_abundance_depth_means, Genus, Relative_Abundance, 3:6)
coral_rel_abundance_depth_means_l$Depth<-factor(coral_rel_abundance_depth_means_l$Depth)
coral_rel_abundance_depth_means_l$Genus<-factor(coral_rel_abundance_depth_means_l$Genus, levels=c("Other", "Acropora", "Porites", "Pocillopora"))


Relative_abundance_fig<-ggplot(coral_rel_abundance_depth_means_l, aes(x=Year, y=Relative_Abundance, fill=Genus))+
  geom_bar(position="stack", stat="identity")+
  facet_wrap(~Depth, ncol=1)+
  scale_fill_manual(values=c("#EA7580","#F8CD9C", "#1BB6AF","#172869"))+
  theme_classic()+
  xlab("Year")+
  ylab("Relative abundance of coral genera")+
  labs(color="")+
  theme(axis.text.x = element_text(colour="black", size=12, hjust=-0.2, angle=-45), 
        axis.text.y = element_text(colour="black", size=14))+
  theme(axis.title.x=element_text(size=14), axis.title.y=element_text(size=16))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.ticks = element_line(color="black"))+
  theme(legend.title=element_blank(), legend.text = element_text(size=16))+
  theme(strip.text.x = element_text(size = 20), strip.background = element_blank())+
  theme(aspect.ratio = 3/5)

#coral_covere_and_rel_abundance_fig<-cowplot::plot_grid(stony_coral_fig, Relative_abundance_fig, align = "vh", ncol = 2, axis = "bt", rel_widths = c(1.125, 1), labels = c("A","B"))
#ggsave("coral_covere_and_rel_abundance_fig.pdf", height=6, width=12, units="in")


#------------------------------------------------------------------------------------------------------------#
#-------------------------- Patterns of coral decline and recovery at 6 LTER sites  -------------------------#
#------------------------------------------------------------------------------------------------------------#

LTER12_color<-ggplot()+
  geom_line(data=subset(coral_cover_all_sites_means, Site=="LTER 1" | Site=="LTER 2"), aes(x=Year, y=mean_stony_coral, group=Depth, color=Depth))+
  geom_errorbar(data=subset(coral_cover_all_sites_means, Site=="LTER 1" | Site=="LTER 2"), aes(x=Year, ymin=mean_stony_coral-se_stony_coral, ymax=mean_stony_coral+se_stony_coral), width=0)+
  geom_point(data=subset(coral_cover_all_sites_means, Site=="LTER 1" | Site=="LTER 2"), aes(x=Year, y=mean_stony_coral, group=Depth, fill=Depth), pch=21, size=3)+
  facet_wrap(~Site, nrow=1)+
  scale_color_manual(values=c("#6DDED2", "#4A918A"), na.translate = F, guide="none")+
  scale_fill_manual(values=c("#6DDED2", "#4A918A"), na.translate = F, labels=c("10 m", "17 m"))+
  theme_classic()+
  #annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, size=1)+
  xlab("Year")+
  ylab("% Cover of stony corals")+
  ylim(c(0,85))+
  labs(color="")+
  theme(axis.text.x = element_text(colour="black", size=8, hjust=-0.2, angle=-45), 
        axis.text.y = element_text(colour="black", size=12))+
  theme(axis.title.x=element_text(size=12), axis.title.y=element_text(size=16))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.ticks = element_line(color="black"))+
  theme(legend.title=element_blank(), legend.text = element_text(size=16))+
  theme(strip.text.x = element_text(size = 20),strip.background = element_blank())+
  theme(aspect.ratio = 4/4)
#ggsave("figures/revision/LTER12_color.pdf", width=6,height=4)

LTER34_color<-ggplot()+
  geom_line(data=subset(coral_cover_all_sites_means, Site=="LTER 3" | Site=="LTER 4"), aes(x=Year, y=mean_stony_coral, group=Depth, color=Depth))+
  geom_errorbar(data=subset(coral_cover_all_sites_means, Site=="LTER 3" | Site=="LTER 4"), aes(x=Year, ymin=mean_stony_coral-se_stony_coral, ymax=mean_stony_coral+se_stony_coral), width=0)+
  geom_point(data=subset(coral_cover_all_sites_means, Site=="LTER 3" | Site=="LTER 4"), aes(x=Year, y=mean_stony_coral, group=Depth, fill=Depth), pch=21, size=3)+
  facet_wrap(~Site, nrow=2)+
  scale_color_manual(values=c("#6DDED2", "#4A918A"), na.translate = F, guide="none")+
  scale_fill_manual(values=c("#6DDED2", "#4A918A"), na.translate = F, labels=c("10 m", "17 m"))+
  theme_classic()+
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, size=1)+
  xlab("Year")+
  ylab("% Cover of stony corals")+
  ylim(c(0,85))+
  labs(color="")+
  theme(axis.text.x = element_text(colour="black", size=8, hjust=-0.2, angle=-45), 
        axis.text.y = element_text(colour="black", size=12))+
  theme(axis.title.x=element_text(size=12), axis.title.y=element_text(size=16))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.ticks = element_line(color="black"))+
  theme(legend.title=element_blank(), legend.text = element_text(size=16))+
  theme(strip.text.x = element_text(size = 20),strip.background = element_blank())+
  theme(aspect.ratio = 4/4)
#ggsave("figures/LTER34_color.pdf", width=4,height=6)

## LTER 3 LTER 4
LTER34_color<-geom_line(data=subset(coral_cover_all_sites_means, Site=="LTER 1" | Site=="LTER 2"), aes(x=Year, y=mean_stony_coral, group=Depth, color=Depth))+
#ggsave("figures/LTER34_color.pdf", width=4,height=6)

LTER34<-geom_line(data=subset(coral_cover_all_sites_means, Site=="LTER 1" | Site=="LTER 2"), aes(x=Year, y=mean_stony_coral, group=Depth, color=Depth))+
#ggsave("figures/LTER34.pdf", width=4,height=6)


coral_cover_all_sites_means_56<-subset(coral_cover_all_sites_means, Site=="LTER 5" | Site=="LTER 6")
coral_cover_all_sites_means_56$Site<-factor(coral_cover_all_sites_means_56$Site, levels=c("LTER 6", "LTER 5"))

## LTER 5 LTER 6
LTER56_color<-ggplot()+
  geom_line(data=subset(coral_cover_all_sites_means, Site=="LTER 5" | Site=="LTER 6"), aes(x=Year, y=mean_stony_coral, group=Depth, color=Depth))+
  geom_errorbar(data=subset(coral_cover_all_sites_means, Site=="LTER 5" | Site=="LTER 6"), aes(x=Year, ymin=mean_stony_coral-se_stony_coral, ymax=mean_stony_coral+se_stony_coral), width=0)+
  geom_point(data=subset(coral_cover_all_sites_means, Site=="LTER 5" | Site=="LTER 6"), aes(x=Year, y=mean_stony_coral, group=Depth, fill=Depth), pch=21, size=3)+
  facet_wrap(~factor(Site, levels=c("LTER 6", "LTER 5")), nrow=2)+
  scale_color_manual(values=c("#6DDED2", "#4A918A"), na.translate = F, guide="none")+
  scale_fill_manual(values=c("#6DDED2", "#4A918A"), na.translate = F, labels=c("10 m", "17 m"))+
  theme_classic()+
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, size=1)+
  xlab("Year")+
  ylab("% Cover of stony corals")+
  ylim(c(0,85))+
  labs(color="")+
  theme(axis.text.x = element_text(colour="black", size=8, hjust=-0.2, angle=-45), 
        axis.text.y = element_text(colour="black", size=12))+
  theme(axis.title.x=element_text(size=12), axis.title.y=element_text(size=16))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.ticks = element_line(color="black"))+
  theme(legend.title=element_blank(), legend.text = element_text(size=16))+
  theme(strip.text.x = element_text(size = 20),strip.background = element_blank())+
  theme(aspect.ratio = 4/4)
#ggsave("figures/LTER56_color.pdf", width=4,height=6)




