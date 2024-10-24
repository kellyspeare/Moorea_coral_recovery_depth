library(emmeans)
library(multcomp)
library(tidyverse)
library(car)
library(ggplot2)
library(DHARMa)

# 7_substratum_coral_mortality
# script analyzes data on substratum availability and juvenile coral mortality

#------------------------------------------------------------------------------------------------------------#
#--------------------------------------------- Data ---------------------------------------------------------
#------------------------------------------------------------------------------------------------------------#

# data from MCR LTER: Coral Reef: Early life stage bottleneck determines rates of coral recovery following severe disturbance; Data for Speare et al., 2024, Ecology
# https://doi.org/10.6073/pasta/f006b56623d3dd61d689237b048b53d8

subst_avail<-read.csv("data_clean/substratum_availability.csv", header = TRUE, stringsAsFactors = TRUE)

#means by depth
subst_depth_means<-ddply(subst_avail, .(Depth, Label), summarise, 
                         mean_abund=mean(percent),
                         se_abund=sd(percent)/sqrt(length(mean_abund)))

# calculating percent from proportion
subst_avail$percent<-subst_avail$Proportion*100

#------------------------------------------------------------------------------------------------------------#
#-------------------- does the abundance (percent cover of substratum type)  differ by depth?-----------------
#------------------------------------------------------------------------------------------------------------#
# glmm 

subst_glmm<-glmmTMB(percent ~ Depth*Label + (1|Site), data=subst_avail,
                     family=gaussian(link="identity"))

car::Anova(subst_glmm, type=2) ## model results with Wald chi squared tests
plot(allEffects(subst_glmm))
plot(simulateResiduals(subst_glmm, quantileFunction = qnorm))
testCategorical(subst_glmm, subst_avail$Label)
testCategorical(subst_glmm, subst_avail$Depth)
testResiduals(subst_glmm)

# Tukey Post hoc
subst_glmm_tukey<-emmeans(subst_glmm, type="response", ~Depth*Label)
cld(subst_glmm_tukey, Letters=letters)

subst_glmm_predicted<-data.frame(emmeans(subst_glmm, type="response", ~Depth*Label))

subst_names <- c(
  `skeleton` = "Dead skeletons",
  `hard` = "Hard substratum",
  `rubble` = "Rubble"
)

substrate_abund_plot<-ggplot()+
  geom_jitter(data=subst_avail, aes(x=Depth, y=percent, group=Depth, color=Depth), width=0.2, alpha=(0.5), shape=1, size=2)+
  geom_errorbar(data=subst_glmm_predicted, aes(x=Depth, ymin=emmean-SE, ymax=emmean+SE, width=0.3))+
  geom_point(data=subst_glmm_predicted, aes(x=Depth, y=emmean, group=Depth, fill=Depth), colour="black",pch=21, size=4)+
  scale_color_manual(values=c("#6DDED2", "#4A918A"), na.translate = F)+
  scale_fill_manual(values=c("#6DDED2", "#4A918A"), na.translate = F)+
  theme_classic()+
  ylab("Availability of substratum (%)")+
  facet_wrap(~Label, labeller = as_labeller(subst_names))+
  xlab("Depth")+
  theme_classic()+
  theme(strip.text = element_text(size = 14, color = "black", face = "bold"))+ #facet text
  theme(axis.text.x = element_text(colour="black", size=14), 
        axis.text.y = element_text(colour="black", size=14))+
  theme(axis.title.x=element_blank(), axis.title.y=element_text(size=16))+
  theme(strip.background = element_rect(colour="white", fill="white"))+
  theme(axis.ticks = element_line(color="black"))+  
  theme(legend.position="none")+
  theme(aspect.ratio = 4/4)
  
#------------------------------------------------------------------------------------------------------------#
#------------------- Proportion of juvenile corals on each substratum type -----------------------------------
#------------------------------------------------------------------------------------------------------------#

# reading in the juvenile coral mortality data 
juv_subst<-read.csv("data_clean/juvenile_coral_mortality.csv", header=TRUE, stringsAsFactors = TRUE)

#---------------- GLMMs 
# for each substratum category, askng whether the probability of being on that substratum type differs by depth

##--------------- Rubble 
juv_subst$BinaryRubble<-ifelse(juv_subst$Substrate=="rubble", 1, 0)

juv_rub_glmm <- glmmTMB(BinaryRubble ~ Depth + (1|Site/Site_Depth_Pole_Quad), 
                         data=juv_subst, family=binomial(link = "logit"))

car::Anova(juv_rub_glmm, type=2) ## model results with Wald chi squared tests
plot(allEffects(juv_rub_glmm))
plot(simulateResiduals(juv_rub_glmm, quantileFunction = qnorm))

juv_rubble_glmm_predicted<-data.frame(emmeans(juv_rub_glmm, type="response", ~Depth))
juv_rubble_glmm_predicted$Label<-as.factor("rubble")

##--------------- Hard substrate
juv_subst$BinaryHard<-ifelse(juv_subst$Substrate=="hard", 1, 0)

juv_hard_glmm <- glmmTMB(BinaryHard ~ Depth + (1|Site/Site_Depth_Pole_Quad), 
                         data=juv_subst, family=binomial(link = "logit"))

car::Anova(juv_hard_glmm, type=2) ## model results with Wald chi squared tests
plot(allEffects(juv_hard_glmm))
plot(simulateResiduals(juv_hard_glmm, quantileFunction = qnorm))

juv_hard_glmm_predicted<-data.frame(emmeans(juv_hard_glmm, type="response", ~Depth))
juv_hard_glmm_predicted$Label<-as.factor("hard")

##--------------- dead Skeleton substrate
juv_subst$BinarySkeleton<-ifelse(juv_subst$Substrate=="deadSkeleton", 1, 0)

juv_skeleton_glmm <- glmmTMB(BinarySkeleton ~ Depth + (1|Site/Site_Depth_Pole_Quad), 
                         data=juv_subst, family=binomial(link = "logit"))

car::Anova(juv_skeleton_glmm, type=2) ## model results with Wald chi squared tests
plot(allEffects(juv_skeleton_glmm))
plot(simulateResiduals(juv_skeleton_glmm, quantileFunction = qnorm))

juv_skeleton_glmm_predicted<-data.frame(emmeans(juv_skeleton_glmm, type="response", ~Depth))
juv_skeleton_glmm_predicted$Label<-as.factor("deadSkeleton")

prob_predicted<-rbind(juv_skeleton_glmm_predicted, juv_rubble_glmm_predicted, juv_hard_glmm_predicted)
prob_predicted$Label <- factor(prob_predicted$Label, levels = c("deadSkeleton", "hard", "rubble"))

## --------------- combining the data for plotting 
# will plot the probability of being on each substratum type
# extract those probabilities from each model

rubble_binary<-juv_subst %>% dplyr::select(Site, Depth, Site_Depth_Pole_Quad, Coral, timeframe, BinaryRubble)
rubble_binary$Label<-as.factor("rubble")
names(rubble_binary)[names(rubble_binary) == "BinaryRubble"] <- "Binary"

hard_binary<-juv_subst %>% dplyr::select(Site, Depth, Site_Depth_Pole_Quad, Coral, timeframe, BinaryHard)
hard_binary$Label<-as.factor("hard")
names(hard_binary)[names(hard_binary) == "BinaryHard"] <- "Binary"

skeleton_binary<-juv_subst %>% dplyr::select(Site, Depth, Site_Depth_Pole_Quad, Coral, timeframe, BinarySkeleton)
skeleton_binary$Label<-as.factor("deadSkeleton")
names(skeleton_binary)[names(skeleton_binary) == "BinarySkeleton"] <- "Binary"

binary_plotting<-rbind(hard_binary, skeleton_binary, rubble_binary)
binary_plotting$Label <- factor(binary_plotting$Label, levels = c("deadSkeleton", "hard", "rubble"))

label_names <- c(
  `deadSkeleton` = "Dead skeletons",
  `hard` = "Hard substrate",
  `rubble` = "Rubble"
)

prob_substrate_plot<-ggplot()+
  geom_jitter(data=binary_plotting, aes(x=Depth, y=Binary, group=Depth, color=Depth), width=0.2,height=0.1, alpha=(0.5), shape=1, size=2)+
  geom_errorbar(data=prob_predicted, aes(x=Depth, ymin=prob-SE, ymax=prob+SE, fill=Depth), width=0.2, position=position_dodge(0.9))+
  geom_point(data=prob_predicted, aes(x=Depth, y=prob, fill=Depth), pch=21,size=4, position = position_dodge(width = 0.9))+
  scale_fill_discrete(guide = "none")+
  scale_color_manual(values=c("#6DDED2", "#4A918A"), labels = c("10 m", "17 m"))+
  scale_fill_manual(values=c("#6DDED2", "#4A918A"), labels = c("10 m", "17 m"))+
  facet_wrap(~Label, labeller=as_labeller(label_names))+
  theme_classic()+
  ylab("Probability of being \non substratum type")+
  xlab("Depth")+
  theme_classic()+
  scale_y_continuous(breaks=seq(0.0, 1.0, 0.2))+
  #ylim(c(-0.3,1.3), breaks=c(0.0, 0.2, 0.4, 0.6, 0.8, 1))+
  theme(strip.text = element_text(size = 14, color = "white", face = "bold"))+ #facet text
  theme(axis.text.x = element_text(colour="black", size=14), 
        axis.text.y = element_text(colour="black", size=14))+
  theme(axis.title.x=element_blank(), axis.title.y=element_text(size=16))+
  theme(strip.background = element_rect(colour="white", fill="white"))+
  theme(axis.ticks = element_line(color="black"))+  
  theme(legend.position="none")+
  theme(aspect.ratio = 4/4)

#------------------------------------------------------------------------------------------------------------#
#------------------- Mortality of juvenile corals on each substratum type -----------------------------------
#------------------------------------------------------------------------------------------------------------#

### probability of mortality on each subsratum type 
### GLMM
#juv_subst$Dead_year_2<-ifelse(juv_subst$Alive_year_2==1, 0, 1)
juv_subst$Substrate<-factor(juv_subst$Substrate)
juv_subst_mort_glmm<-glmmTMB(Dead_year_2 ~ Substrate*Depth + (1|Site/Site_Depth_Pole_Quad), 
                             data=juv_subst, family=binomial(link = "logit"))

Anova(juv_subst_mort_glmm)
car::Anova(juv_subst_mort_glmm, type=2) ## model results with Wald chi squared tests
plot(allEffects(juv_subst_mort_glmm))
check_overdispersion(juv_subst_mort_glmm)
plot(simulateResiduals(juv_subst_mort_glmm, quantileFunction = qnorm))

#tukey post hoc
juv_subst_mort_glmm_tukey<-emmeans(juv_subst_mort_glmm, type="response", ~Depth*Substrate)
#pairs(juv_subst_mort_glmm_tukey, adjust="tukey")
cld(juv_subst_mort_glmm_tukey, Letters=letters)

#predicted data frame
juv_subst_mort_glmm_predicted<-data.frame(emmeans(juv_subst_mort_glmm, type="response", ~Depth*Substrate))

# label names
subst_names2 <- c(
  `deadSkeleton` = "Dead skeletons",
  `hard` = "Hard substrate",
  `rubble` = "Rubble"
)

mort_on_substrate_plot<-ggplot()+
  geom_jitter(data=juv_subst, aes(x=Depth, y=Dead_year_2, group=Depth, color=Depth), width=0.2,height=0.1, alpha=(0.5), shape=1, size=2)+
  geom_errorbar(data=juv_subst_mort_glmm_predicted, aes(x=Depth, ymin=prob-SE, ymax=prob+SE, width=0.3))+
  geom_point(data=juv_subst_mort_glmm_predicted, aes(x=Depth, y=prob, group=Depth, fill=Depth), colour="black",pch=21, size=4)+
  scale_color_manual(values=c("#6DDED2", "#4A918A"), na.translate = F)+
  scale_fill_manual(values=c("#6DDED2", "#4A918A"), na.translate = F)+
  theme_classic()+
  ylab("Probability of mortality")+
  facet_wrap(~Substrate, labeller = as_labeller(subst_names2))+
  xlab("Depth")+
  theme_classic()+
  scale_y_continuous(breaks=seq(0.0, 1.0, 0.2))+
  theme(strip.text = element_text(size = 14, color = "white", face = "bold"))+ #facet text
  theme(axis.text.x = element_text(colour="black", size=14), 
        axis.text.y = element_text(colour="black", size=14))+
  theme(axis.title.x=element_blank(), axis.title.y=element_text(size=16))+
  theme(strip.background = element_rect(colour="white", fill="white"))+
  theme(axis.ticks = element_line(color="black"))+  
  theme(legend.position="none")+
  theme(aspect.ratio = 4/4)

mort_on_substrate_plot_legend<-ggplot()+
  geom_jitter(data=juv_subst, aes(x=Depth, y=Dead_year_2, group=Depth, color=Depth), width=0.3,height=0.1, alpha=(0.5), shape=1, size=2)+
  geom_errorbar(data=juv_subst_mort_glmm_predicted, aes(x=Depth, ymin=prob-SE, ymax=prob+SE, width=0.3))+
  geom_point(data=juv_subst_mort_glmm_predicted, aes(x=Depth, y=prob, group=Depth, fill=Depth), colour="black",pch=21, size=4)+
  scale_color_manual(values=c("#6DDED2", "#4A918A"), na.translate = F)+
  scale_fill_manual(values=c("#6DDED2", "#4A918A"), na.translate = F)+
  theme_classic()+
  ylab("Probability of mortality")+
  facet_wrap(~Substrate, labeller = as_labeller(subst_names2))+
  xlab("Depth")+
  theme_classic()+
  scale_y_continuous(breaks=seq(0.0, 1.0, 0.2))+
  theme(strip.text = element_text(size = 14, color = "white", face = "bold"))+ #facet text
  theme(axis.text.x = element_text(colour="black", size=14, hjust=0), 
        axis.text.y = element_text(colour="black", size=14))+
  theme(axis.title.x=element_blank(), axis.title.y=element_text(size=16))+
  theme(strip.background = element_rect(colour="white", fill="white"))+
  theme(axis.ticks = element_line(color="black"))+  
  #theme(legend.position="none")+
  theme(aspect.ratio = 4/4)

substrate_legend <- cowplot::get_legend(mort_on_substrate_plot_legend + theme(legend.box.margin = margin(0, 0, 0, 12)))


#### ------------------------ Plot Figure ------------------------

substrate_fig<-cowplot::plot_grid(substrate_abund_plot, prob_substrate_plot, 
                   mort_on_substrate_plot, nrow = 3, align="vh",labels=c("(a)","(b)","(c)"))

Figure4<-cowplot::plot_grid(substrate_fig, substrate_legend, rel_widths = c(9,1))

ggsave("figures/revision/Figure4.pdf", width=9, height=8.5, units="in")

#ggsave("figures/legend_fig4.pdf", width=4, height=4, units="in")
#Figure4<-cowplot::plot_grid(juvenile_mortality_fig_v2, legend_fig4, rel_widths = c(3.5, 0.5))
#ggsave("figures/Figure4_v2.pdf", width=15.5, height=4.2, units="in")


#------------------------------------------------------------------------------------------------------------#
# ------------------------- Relative contribution to mortality calculations ---------------------------------
#------------------------------------------------------------------------------------------------------------#

# determine the relative contributions of differences in the availability of substratum types 
# and differences in mortality rates on dead coral skeletons to overall mortality patterns,

juv_subst %>% group_by(Depth) %>% summarize(total=sum(BinaryRubble)) 

#----- 10m -----
# calculate number of individuals on each substrate type:
# dead skeletons
juv_subst %>% filter(Depth=="10m") %>% summarize(total=sum(BinarySkeleton)) #305
# hard
juv_subst %>% filter(Depth=="10m") %>% summarize(total=sum(BinaryHard)) #389
# rubble
juv_subst %>% filter(Depth=="10m") %>% summarize(total=sum(BinaryRubble)) #16
# total 
juv_subst %>% filter(Depth=="10m") %>% summarize(total=length(BinaryRubble)) #710

mort_10<-filter(juv_subst_mort_glmm_predicted[,c(1:3)], Depth=="10m") #mortality rate from data
names(mort_10)[names(mort_10)=="prob"]<-"mort_rate"
mort_10$no_juvs<-c(305, 389, 16) #number of juveniles
mort_10$prop_juvs<-mort_10$no_juvs/710 #proportion of total juveniles 
mort_10$no_survived<-(mort_10$no_juvs)*(1-mort_10$mort_rate)# num of indiv surviving = (1-mort_rate) * number of individuals

#calc total prob of mortality at depth
mort_at_10<-(sum(mort_10$no_juvs) - sum(mort_10$no_survived))/sum(mort_10$no_juvs) #0.2044324 is overall mortality rate

#----- 17m -----
# calculate number of individuals on each substrate type:
# dead skeletons
juv_subst %>% filter(Depth=="17m") %>% summarize(total=sum(BinarySkeleton)) #169
# hard
juv_subst %>% filter(Depth=="17m") %>% summarize(total=sum(BinaryHard)) #458
# rubble
juv_subst %>% filter(Depth=="17m") %>% summarize(total=sum(BinaryRubble)) #116
# total 
juv_subst %>% filter(Depth=="17m") %>% summarize(total=length(BinaryRubble)) #743

mort_17<-filter(juv_subst_mort_glmm_predicted[,c(1:3)], Depth=="17m") #mortality rate from data
names(mort_17)[names(mort_17)=="prob"]<-"mort_rate"
mort_17$no_juvs<-c(169, 458, 116) #number of juveniles
mort_17$prop_juvs<-mort_17$no_juvs/743 #proportion of total juveniles 
mort_17$no_survived<-(mort_17$no_juvs)*(1-mort_17$mort_rate)# num of indiv surviving = (1-mort_rate) * number of individuals

#calc total prob of mortality at depth
mort_at_17<-(sum(mort_17$no_juvs) - sum(mort_17$no_survived))/sum(mort_17$no_juvs) #0.362642 is overall mortality rate

#----- projected at 17m ----- 
# if the proportion of juveniles on each substrate type at 17m was the same as observed at 10m
mort_17_proj<-mort_10[,c(2,5)] #proportion of juveniles across substrate types at 10m
mort_17_proj$no_juvs<-mort_17_proj$prop_juvs*743 #times total number of juveniles at 17m
mort_17_proj<-left_join(mort_17_proj, mort_17[,c(2,3)], by="Substrate")# probability of mortality on each substrate type at 17m
mort_17_proj$no_surviving<-(mort_17_proj$no_juvs)*(1-mort_17_proj$prob)
mort_at_17_projected <- (743 - sum(mort_17_proj$no_surviving))/743 #total number of juvs at 17m -  number that survived divided by total number
#projected mortality rate of 0.332629

# how much higher is mortality at 17m than 10m
depth_diff_mort <- mort_at_17 - mort_at_10 #0.1582095

# how much of this difference in mortality can be attributed to higher mortality rate
mort_at_17_projected - mort_at_17
