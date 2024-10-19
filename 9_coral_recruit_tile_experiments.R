
library(ggplot2)
library(ggthemes)
library(plyr)
library(dplyr)
library(tidyr)
library(lme4)
library(glmmTMB)

# 9_coral_recruit_tile_experiments
# script analyzes data from experiments with coral settlement tiles

# data from MCR LTER: Coral Reef: Early life stage bottleneck determines rates of coral recovery following severe disturbance; Data for Speare et al., 2024, Ecology
# https://doi.org/10.6073/pasta/f006b56623d3dd61d689237b048b53d8

#------------------------------------------------------------------------------------------------------------#
#-------------------------- Data - recruitment tile caging experiment-----------------------------------------
#------------------------------------------------------------------------------------------------------------#


caging_expt<-read.csv("data_clean/coral_recruit_caging_experiment.csv", header=TRUE, stringsAsFactors = TRUE)

caging_expt$Tile<-as.factor(caging_expt$Tile)

# summarizing the number of recruits per tile
recs_per_tile<-ddply(caging_expt, .(Tile, Depth, Caging), summarize,
                     no_recruits=length(Tile))

#mean number of settlers on tiles that were caged vs uncaged
aggregate(no_recruits ~ Caging, FUN=mean, data=recs_per_tile) #Caged 5.9, Uncaged 2.4

#mean number of settlers on tiles that were at 10m vs 17m
aggregate(no_recruits ~ Depth, FUN=mean, data=recs_per_tile) #10m 4.1, 17m 4.6

# do the number of recruits per tile differ between depths for caged/uncaged groups
t.test(no_recruits ~ Depth, data=subset(recs_per_tile, Caging=="Uncaged")) # no difference p-value = 0.4904
t.test(no_recruits ~ Depth, data=subset(recs_per_tile, Caging=="Caged")) #no difference p-value = 0.6972

#### ------------------------- Recruit tile caging experiment - GLMM ------------------------------

#is there a relationship between depth and caging and # of recruits on a tile?
rec_glm<-glmmTMB(no_recruits ~ Depth + Caging + Depth*Caging, data=recs_per_tile, family = nbinom2(link="log"))
car::Anova(rec_glm, type=2) ## significant effect of caging
plot(simulateResiduals(fittedModel=rec_glm, quantileFunction = qnorm)) # all looks good

recruit_dens_glm_predicted<-data.frame(emmeans(rec_glm, type="response", ~Depth*Caging))

#Figure 5A
recruits_per_tile_plot<-ggplot()+
  geom_point(data=recs_per_tile, aes(x=Depth, y=no_recruits, shape=Caging, color=Depth), size=3, alpha=0.5, position=position_jitterdodge(jitter.width = 0.35, dodge.width=0.9))+
  # scale_shape_manual(values=c(2,1))+
  scale_color_manual(values=c("#6DDED2", "#4A918A"))+
  geom_errorbar(data=recruit_dens_glm_predicted, aes(x=Depth, ymin=response-SE, ymax=response+SE, group=Caging), width=0.35, position=position_dodge(0.9))+
  geom_point(data=recruit_dens_glm_predicted, aes(x=Depth, y=response, fill=Depth, shape=Caging),size=3, position = position_dodge(width = 0.9), inherit.aes=FALSE)+
  scale_shape_manual(values=c(24, 21))+
  scale_fill_manual(values=c("#6DDED2", "#4A918A"))+
  xlab("Depth")+
  ylab("Recruits per tile")+
  ggtitle("")+
  theme_classic()+
  theme(legend.text=element_text(size=16), legend.title=element_blank())+
  theme(axis.title=element_text(size=16))+
  theme(axis.text.x=element_text(size=16, color="black"), axis.text.y=element_text(size=16, color="black"))+
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),axis.title.y=element_text(vjust=1.2), axis.title.x=element_text(vjust=-0.4))+
  theme(panel.border=element_blank(), axis.line=element_line()) + 
  theme(panel.grid.minor=element_blank())+
  #theme(legend.position = "none")+
  theme(aspect.ratio = 4/4)

ggplot()+
  geom_point(data=recruit_dens_glm_predicted, aes(x=Depth, y=response, fill=Depth, shape=Caging),size=3, position = position_dodge(width = 0.9))+
  scale_shape_manual(values=c(24, 21))+
  scale_fill_manual(values=c("#6DDED2", "#4A918A"))+
  xlab("Depth")+
  ylab("Recruits per tile")+
  ggtitle("")+
  theme_classic()+
  theme(legend.text=element_text(size=16), legend.title=element_blank())+
  theme(axis.title=element_text(size=16))+
  theme(axis.text.x=element_text(size=16, color="black"), axis.text.y=element_text(size=16, color="black"))+
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),axis.title.y=element_text(vjust=1.2), axis.title.x=element_text(vjust=-0.4))+
  theme(panel.border=element_blank(), axis.line=element_line()) + 
  theme(panel.grid.minor=element_blank())+
  #theme(legend.position = "none")+
  theme(aspect.ratio = 4/4)


shape_legend<-ggplot()+
  geom_point(data=recruit_dens_glm_predicted, aes(x=Depth, y=response, color=Depth, shape=Caging),size=3, position = position_dodge(width = 0.9))+
  scale_shape_manual(values=c(24, 21))+
  scale_color_manual(values=c("#6DDED2", "#4A918A"))+
  xlab("Depth")+
  ylab("Recruits per tile")+
  ggtitle("")+
  theme_classic()+
  theme(legend.text=element_text(size=16), legend.title=element_blank())+
  theme(axis.title=element_text(size=16))+
  theme(axis.text.x=element_text(size=16, color="black"), axis.text.y=element_text(size=16, color="black"))+
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),axis.title.y=element_text(vjust=1.2), axis.title.x=element_text(vjust=-0.4))+
  theme(panel.border=element_blank(), axis.line=element_line()) + 
  theme(panel.grid.minor=element_blank())+
  #theme(legend.position = "none")+
  theme(aspect.ratio = 4/4)


#------------------------------------------------------------------------------------------------------------#
#-------------------------- Data - recruitment tile transplant experiment-----------------------------------------
#------------------------------------------------------------------------------------------------------------#

transplant_expt<-read.csv("data/recruit_transplant_experiment.csv", header=TRUE, stringsAsFactors = TRUE)

# calculate percent mortality by tile 
#first need a count of the number of recruits on each tile that were then transplanted
transplant_expt$Tile<-as.factor(transplant_expt$Tile)

transp_recs<-ddply(transplant_expt, .(Tile, OriginalDepth, TransplantedTo), summarize, initial_recs=length(CoralNumber))

surviving_recs<-ddply(subset(transplant_expt, BinarySurvival=="Alive"), .(Tile), summarize, surviving_recs=length(CoralNumber))

dead_recs<-ddply(subset(transplant_expt, BinarySurvival!="Alive"), .(Tile), summarize, dead_recs=length(CoralNumber))

transp_recs<-left_join(transp_recs, dead_recs, by="Tile")
transp_recs<-left_join(transp_recs, surviving_recs, by="Tile")

transp_recs[is.na(transp_recs)]<-0 

#simple stats from start of experiment
mean(transp_recs$initial_recs) #mean recruits per tile: 6.2
median(transp_recs$initial_recs) #median recruits per tile: 5.5

#calculating percent mortality
transp_recs$mortality<-(transp_recs$dead_recs / transp_recs$initial_recs)*100

#creating new factor for transplant treatment
transp_recs$transp_trt<-as.factor(ifelse(transp_recs$OriginalDepth=="17m" & transp_recs$TransplantedTo=="17m", "17_to_17",
                                         ifelse(transp_recs$OriginalDepth=="17m" & transp_recs$TransplantedTo=="10m", "17_to_10",
                                                ifelse(transp_recs$OriginalDepth=="10m" & transp_recs$TransplantedTo=="10m", "10_to_10", "10_to_17"))))

#renaming original depth and transplant depth factors
transp_recs$OriginalDepth<-revalue(transp_recs$OriginalDepth, c("10m"="from10","17m"="from17"))
transp_recs$TransplantedTo<-revalue(transp_recs$TransplantedTo, c("10m"="to10","17m"="to17"))

#creating new factor for plotting (same as "TransplantedTo" just renamed)
transp_recs$Depth<-as.factor(ifelse(transp_recs$TransplantedTo=="to10", "10m", "17m"))

# t-tests comparing the initial densities of recruits on tiles
t.test(initial_recs ~ TransplantedTo, data=transp_recs) #p-value = 0.4856
t.test(initial_recs ~ OriginalDepth, data=transp_recs) #p-value = 0.5918

#### ------------------------- Recruit tile transplant experiment - GLMM ------------------------------

## does mortality differ by depth? 
transplant_expt$BinaryMortality<-ifelse(transplant_expt$BinarySurvival=="Alive", 0, 1)

rec_mort_glmm <- glmmTMB(BinaryMortality ~ TransplantedTo + (1|Tile) + (1|OriginalDepth), 
                         data=transplant_expt, family=binomial(link = "logit"))

car::Anova(rec_mort_glmm, type=2) ## model results with Wald chi squared tests
plot(allEffects(rec_mort_glmm))
check_overdispersion(rec_mort_glmm)
plot(simulateResiduals(rec_mort_glmm, quantileFunction = qnorm))
# Depth P=0.002332 **, model checks look good

recruit_mort_glmm_predicted<-data.frame(emmeans(rec_mort_glmm, type="response", ~TransplantedTo))

recruit_mortality_plot<-ggplot()+
  geom_point(data=transplant_expt, aes(x=TransplantedTo, y=BinaryMortality, color=TransplantedTo), 
             size=3, alpha=0.4, shape=1, position=position_jitterdodge(jitter.width = 0.7,jitter.height = 0.05, dodge.width=0.9))+
  geom_errorbar(data=recruit_mort_glmm_predicted, aes(x=TransplantedTo, ymin=prob-SE, ymax=prob+SE), position=position_dodge(0.64),color="black", width=0.35)+
  geom_point(data=recruit_mort_glmm_predicted, aes(x=TransplantedTo, y=prob, fill=TransplantedTo), shape=21, size=3, position = position_dodge(width = 0.9))+
  scale_color_manual(values=c("#6DDED2", "#4A918A"))+
  scale_fill_manual(values=c("#6DDED2", "#4A918A"))+
  xlab("Depth")+
  ylab("Probability of mortality")+
  ggtitle("")+
  theme_classic()+
  theme(legend.text=element_text(size=16), legend.title=element_blank())+
  theme(axis.title=element_text(size=16))+
  theme(axis.text.x=element_text(size=16, color="black"), axis.text.y=element_text(size=16, color="black"))+
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),axis.title.y=element_text(vjust=1.2), axis.title.x=element_text(vjust=-0.4))+
  theme(panel.border=element_blank(), axis.line=element_line()) + 
  theme(panel.grid.minor=element_blank())+
  #theme(legend.position = "none")+
  theme(aspect.ratio = 4/4)


# ------------------------- recruitment tile transplant expt - dead recruit skeletal condition  -------------------------

transplant_expt_dead<-subset(transplant_expt, BinarySurvival=="Dead")

## GLMM
transplant_expt_dead$BinaryDamagedRemoved<-ifelse(transplant_expt_dead$SkeletonStructure=="Intact", 0, 1)
# does the probability of being damaged or removed differ by depth -- not how we asked this question originally
rec_skel_glmm <- glmmTMB(BinaryDamagedRemoved ~ TransplantedTo + (1|Tile) + (1|OriginalDepth), 
                         data=transplant_expt_dead, family=binomial(link = "logit"))


car::Anova(rec_skel_glmm, type=2) ## model results with Wald chi squared tests
plot(allEffects(rec_skel_glmm))
check_overdispersion(rec_skel_glmm)
plot(simulateResiduals(rec_skel_glmm, quantileFunction = qnorm))

recruit_skel_glmm_predicted<-data.frame(emmeans(rec_skel_glmm, type="response", ~TransplantedTo))

recruit_skeleton_plot<-ggplot()+
  geom_point(data=transplant_expt_dead, aes(x=TransplantedTo, y=BinaryDamagedRemoved, color=TransplantedTo), shape=1,
             size=3, alpha=0.4, position=position_jitterdodge(jitter.width = 0.7,jitter.height = 0.05, dodge.width=0.9))+
  geom_errorbar(data=recruit_skel_glmm_predicted, aes(x=TransplantedTo, ymin=prob-SE, ymax=prob+SE, fill=TransplantedTo), width=0.35, position=position_dodge(0.9))+
  geom_point(data=recruit_skel_glmm_predicted, aes(x=TransplantedTo, y=prob, fill=TransplantedTo), shape=21,size=3, position = position_dodge(width = 0.9))+
  scale_fill_manual(values=c("#6DDED2", "#4A918A"), labels = c("10 m", "17 m"))+
  scale_color_manual(values=c("#6DDED2", "#4A918A"), labels = c("10 m", "17 m"))+
  xlab("Depth")+
  ylab("Probability of being damaged or removed")+
  ggtitle("")+
  theme_classic()+
  theme(legend.text=element_text(size=16), legend.title=element_blank())+
  theme(axis.title=element_text(size=16))+
  theme(axis.text.x=element_text(size=16, color="black"), axis.text.y=element_text(size=16, color="black"))+
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),axis.title.y=element_text(vjust=1.2), axis.title.x=element_text(vjust=-0.4))+
  theme(panel.border=element_blank(), axis.line=element_line()) + 
  theme(panel.grid.minor=element_blank())+
  theme(legend.position = "none")+
  theme(aspect.ratio = 4/4)

color_legend<-ggplot()+
  geom_point(data=recruit_skel_glmm_predicted, aes(x=TransplantedTo, y=prob, fill=TransplantedTo), shape=21,size=3, position = position_dodge(width = 0.9))+
  scale_fill_manual(values=c("#6DDED2", "#4A918A"), labels = c("10 m", "17 m"))+
  xlab("Depth")+
  ylab("Probability of being damaged or removed")+
  ggtitle("")+
  theme_classic()+
  theme(legend.text=element_text(size=16), legend.title=element_blank())+
  theme(axis.title=element_text(size=16))+
  theme(axis.text.x=element_text(size=16, color="black"), axis.text.y=element_text(size=16, color="black"))+
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),axis.title.y=element_text(vjust=1.2), axis.title.x=element_text(vjust=-0.4))+
  theme(panel.border=element_blank(), axis.line=element_line()) + 
  theme(panel.grid.minor=element_blank())+
  theme(aspect.ratio = 4/4)
color_legend2 <- get_legend(color_legend + theme(legend.box.margin = margin(0, 0, 0, 0)))

## t test: are more dead recruits damaged or removed than intact?
dead_rec_tile_sums2<-dead_rec_tile_sums

dead_rec_tile_sums2$SkeletonStructure2<-ifelse(dead_rec_tile_sums2$SkeletonStructure=="Intact", "Intact", "Damaged or removed")

dead_rec_tile_sums2<-ddply(dead_rec_tile_sums2, .(Tile, TransplantedTo, SkeletonStructure2), summarize,
                           sum_percent=sum(percent))

t.test(sum_percent ~ SkeletonStructure2, data=dead_rec_tile_sums2) #P<0.0001

t.test(sum_percent ~ TransplantedTo, data=subset(dead_rec_tile_sums2, SkeletonStructure2=="Damaged or removed"))



# ------------------------------------------ plot figure 5  ---------------------------------------#

Figure5<-cowplot::plot_grid(recruits_per_tile_plot, recruit_mortality_plot, skeleton_plot, nrow = 1, align="vh", labels=c("A","B","C"))
ggsave("figures/Figure5.pdf", width=14.5, height=4, units="in")

Figure5_v2<-cowplot::plot_grid(recruits_per_tile_plot_v2, recruit_mortality_plot_v2, skeleton_plot_v2, nrow = 1, align="vh", labels=c("(a)","(b)","(c)"))
ggsave("figures/Figure5_v2.pdf", width=16, height=4, units="in")


















