### Loading Packages
library(openxlsx)
library(phytools)
library(ape)
library(ggplot2)
library(lsr)
library(adephylo)
library(missForest)
library(treeplyr)
library(Rmisc)
library(lme4qtl)
library(car)
library(ggpubr)
library(patchwork)
library(tidyr)
library(BestFitM)
library(ggtrendline)
library(pscl)
library(smatr)
library(deming)
library(AICcmodavg)
library(ggsignif)
library(lme4)
library(lmerTest)
setwd(getwd())

### Loading additional informations
### Common species
Common_sp_list = read.xlsx("Data/Common_species_list.xlsx", sheet = "Common_sp_list", colNames = TRUE, rowNames = FALSE)
Common_sp_list_SLA = c(na.omit(Common_sp_list$Species_SLA))
Common_sp_list_AGB = c(na.omit(Common_sp_list$Species_AGB))
Common_sp_list$Origin = factor(Common_sp_list$Origin, levels = c("Native","Exotic"))
#knitr::kable(head(Common_sp_list))
### pot experiment taits data 
pot_trait = read.xlsx("Data/Pot_traits_mean.xlsx", sheet = "Pot_means", colNames = TRUE, rowNames = FALSE)
colnames(pot_trait) <- paste0(colnames(pot_trait), "_green")
colnames(pot_trait)[1] <- "Species"
head(pot_trait)

################################################################################
############################ Fig 3 & Fig S4 ####################################
################################################################################
mytheme = theme(panel.background = element_rect(fill='white', colour='black'),
                panel.grid=element_blank(), 
                legend.position = "none",
                axis.title = element_text(color='black',size=24),
                axis.ticks.length = unit(0.4,"lines"), axis.ticks = element_line(color='black'),
                axis.line = element_line(colour = "black"), 
                axis.title.x=element_text(colour='black', size=13,vjust = 1),
                axis.title.y=element_text(colour='black', size=13,vjust = 1),
                axis.text=element_text(colour='black',size=11))+
  theme(plot.title = element_text(color = "black", size = 13, hjust = 0.5))

### The relationship between characters and relative biomass in the second year
composition_change <- read.xlsx("Data/Field_composition_database.xlsx",sheet = "Field_composition", rowNames = FALSE, colNames = TRUE)
unique(composition_change$Species)
#head(composition_change)
composition_change$Origin = factor(composition_change$Origin, levels = c("Native","Exotic"))

#########################################################################################
### Select common species
composition_change_common <- composition_change[(composition_change$Species %in% Common_sp_list_AGB), ]
length(unique(composition_change_common$Species))
##
composition_change_common$rebio2020_100 = sqrt(composition_change_common$rebio2020*100)
composition_change_common$Block = as.factor(composition_change_common$Block)
composition_change_common$Plot_num = as.factor(composition_change_common$Plot_num)
###
composition_change_common = composition_change_common %>% drop_na(Field_SLA_imp)
length(unique(composition_change_common$Species))
### Remove the missing data of relative biomass in the second year
Com_relative_bio = composition_change_common %>% drop_na(rebio2020)
length(unique(Com_relative_bio$Species))

### Comparison of relative Biomass between exotic and Native plants 
### in Field experiment in the second year
Com_relative_bio2 = Com_relative_bio
mod12 <- lmer(rebio2020_100 ~ Origin +  (1|Block/Plot_num), data=Com_relative_bio2)
summary(mod12)
anova(mod12)
Anova(mod12, type="II", test.statistic ="Chisq")

### Fig S4b
P1 = ggplot(Com_relative_bio2,aes(x=Origin,y=rebio2020_100))+
  geom_violin(data=Com_relative_bio2, aes(y=rebio2020_100,x=Origin,fill=Origin,color=Origin),trim=T,
              scale = "width",position = position_dodge(0.5),width=0.6,alpha = 0.8, size = 1.5)+
  geom_boxplot(data=Com_relative_bio2, aes(y=rebio2020_100,x=Origin),fill="white",color="black",
               position = position_dodge(0.5),width=0.15,outlier.size=0.8,outlier.shape = 1)+
  stat_summary(data=Com_relative_bio2, aes(y=rebio2020_100,x=Origin),fill="white",color="black",
               fun="mean",position = position_dodge(0.5),geom="point",shape=21, size=1.5)+
  scale_color_manual(values = c('#60A7A6','#FEA6A6'))+
  scale_fill_manual(values = c('#60A7A6','#FEA6A6'))+
  mytheme+
  labs(x = NULL, y = 'Relative abundance within the \n plot in second year (%, sqrt) ' ,title = NULL) + 
  scale_y_continuous(position = "left",labels = scales::label_comma(accuracy =0.01))+
  geom_signif(comparisons = list(c(1,2)),test="t.test", annotations='p = 0.018',tip_length = 0.02,size = 0.5,
              textsize = 4,y_position = 8.8); P1


### Relationship between specific leaf area and relative biomass of plants 
### in the second year in field experiments
### (Inferred trait data)
Com_relative_bio2 = Com_relative_bio %>% drop_na(Field_SLA_imp)
length(unique(Com_relative_bio2$Species))
Com_relative_bio2$Field_SLA_imp_log = log10(Com_relative_bio2$Field_SLA_imp)

mod12 <- lmer(rebio2020_100 ~ (Field_SLA_imp_log) +  (1|Block/Plot_num), data=Com_relative_bio2)
summary(mod12)
anova(mod12)
MuMIn::r.squaredGLMM(mod12)

Com_relative_bio2$F0 = predictSE(mod12, Com_relative_bio2, level = 0)$fit
Com_relative_bio2$SE <- predictSE(mod12, Com_relative_bio2, level = 0)$se.fit

SLA_lmer1 = ggplot(Com_relative_bio2, aes(x=(Field_SLA_imp_log), y=rebio2020_100)) +
  geom_point(aes(fill = Origin),size = 2.2, pch = 21) + # fill = "#00000022", color = "#0A0A0A"
  geom_line(aes(y=F0), size=1, linetype = 2) + #color = "#8B0000"
  labs(x = ' \n ',
       y = 'Relative abundance within the \n plot in second year (%, sqrt)',
       title = "Specific leaf area") +
  mytheme + 
  scale_fill_manual(values = c('#60A7A6','#FEA6A6'))+
  scale_x_continuous(labels = scales::label_comma(accuracy =0.01)) +
  scale_y_continuous(labels = scales::label_comma(accuracy =0.01)) + 
  geom_line(aes(y =  F0 - 1.96 * SE), color = "grey40", linetype = "dashed", size=0.6) +
  geom_line(aes(y =  F0 + 1.96 * SE), color = "grey40", linetype = "dashed", size=0.6) +
  annotate('text',x = 2.2,y=7.5,label=('R2 = 0.004 \n p = 0.560'),size=4,color='black'); SLA_lmer1

### Remove perennials
Com_relative_bio3 = subset(Com_relative_bio2, Lifeform != "Perennial")
mod12 <- lmer(rebio2020_100 ~ Field_SLA_imp_log +  (1|Block/Plot_num), data=Com_relative_bio3)
summary(mod12)
anova(mod12)

MuMIn::r.squaredGLMM(mod12)
Com_relative_bio3$F0 = predictSE(mod12, Com_relative_bio3, level = 0)$fit
Com_relative_bio3$SE <- predictSE(mod12, Com_relative_bio3, level = 0)$se.fit

SLA_lmer2 = ggplot(Com_relative_bio3, aes(x=(Field_SLA_imp_log), y=rebio2020_100)) +
  geom_point(aes(fill = Origin),size = 2.2, pch = 21) + # fill = "#00000022", color = "#0A0A0A"
  geom_line(aes(y=F0), size=1, linetype = 2) + #color = "#8B0000"
  labs(x = ' \n ',
       y = 'Relative abundance within the \n plot in second year (%, sqrt)',
       title = "Specific leaf area") +
  mytheme + 
  scale_fill_manual(values = c('#60A7A6','#FEA6A6'))+
  scale_x_continuous(labels = scales::label_comma(accuracy =0.01)) +
  scale_y_continuous(labels = scales::label_comma(accuracy =0.01)) + 
  geom_line(aes(y =  F0 - 1.96 * SE), color = "grey40", linetype = "dashed", size=0.6) +
  geom_line(aes(y =  F0 + 1.96 * SE), color = "grey40", linetype = "dashed", size=0.6) +
  annotate('text',x = 2.2,y=7.5,label=('R2 = 0.003 \n p = 0.634'),size=4,color='black'); SLA_lmer2


################################################################################
### Relationship between maximum Plant height and relative Biomass 
### in the second year in Field experiment
Com_relative_bio2 = composition_change_common %>% drop_na(rebio2020)
length(unique(Com_relative_bio2$Species))

Com_relative_bio2$Field_Hmax_log = log10(Com_relative_bio2$Field_Hmax)
###
Com_relative_bio2$Block = as.factor(Com_relative_bio2$Block)
Com_relative_bio2$Plot_num = as.factor(Com_relative_bio2$Plot_num)

mod12 <- lmer(rebio2020_100 ~ Field_Hmax_log +  (1|Block/Plot_num), data=Com_relative_bio2)
anova(mod12)

MuMIn::r.squaredGLMM(mod12)
Com_relative_bio2$F0 = predictSE(mod12, Com_relative_bio2, level = 0)$fit
Com_relative_bio2$SE <- predictSE(mod12, Com_relative_bio2, level = 0)$se.fit

Hmax_lmer1 = ggplot(Com_relative_bio2, aes(x=(Field_Hmax_log), y=rebio2020_100)) +
  geom_point(aes(fill = Origin),size = 2.2, pch = 21) + # fill = "#00000022", color = "#0A0A0A"
  geom_line(aes(y=F0), size=1, linetype = 1) + #color = "#8B0000"
  labs(x = 'Traits value (log10) estimated in field experiment',
       y = ' \n ',
       title = "Maximum height") +
  mytheme + 
  scale_fill_manual(values = c('#60A7A6','#FEA6A6'))+
  scale_x_continuous(labels = scales::label_comma(accuracy =0.01)) +
  scale_y_continuous(labels = scales::label_comma(accuracy =0.01)) + 
  geom_line(aes(y =  F0 - 1.96 * SE), color = "grey40", linetype = "dashed", size=0.6) +
  geom_line(aes(y =  F0 + 1.96 * SE), color = "grey40", linetype = "dashed", size=0.6) +
  annotate('text',x = 1.7,y=7.5,label=('R2 = 0.215 \n p < 0.001'),size=4,color='black'); Hmax_lmer1

### Remove perennials
Com_relative_bio3 = subset(Com_relative_bio2, Lifeform != "Perennial")
mod12 <- lmer(rebio2020_100 ~ Field_Hmax_log +  (1|Block/Plot_num), data=Com_relative_bio3)
anova(mod12)

MuMIn::r.squaredGLMM(mod12)
Com_relative_bio3$F0 = predictSE(mod12, Com_relative_bio3, level = 0)$fit
Com_relative_bio3$SE <- predictSE(mod12, Com_relative_bio3, level = 0)$se.fit

Hmax_lmer2 = ggplot(Com_relative_bio3, aes(x=(Field_Hmax_log), y=rebio2020_100)) +
  geom_point(aes(fill = Origin),size = 2.2, pch = 21) + # fill = "#00000022", color = "#0A0A0A"
  geom_line(aes(y=F0), size=1, linetype = 1) + #color = "#8B0000"
  labs(x = 'Traits value (log10) estimated in field experiment',
       y = ' \n ',
       title = "Maximum height") +
  mytheme + 
  scale_fill_manual(values = c('#60A7A6','#FEA6A6'))+
  scale_x_continuous(labels = scales::label_comma(accuracy =0.01)) +
  scale_y_continuous(labels = scales::label_comma(accuracy =0.01)) + 
  geom_line(aes(y =  F0 - 1.96 * SE), color = "grey40", linetype = "dashed", size=0.6) +
  geom_line(aes(y =  F0 + 1.96 * SE), color = "grey40", linetype = "dashed", size=0.6) +
  annotate('text',x = 1.7,y=7.5,label=('R2 = 0.149 \n p < 0.001'),size=4,color='black')

Hmax_lmer2

### Relationship between Aboveground biomass and relative Biomass 
### in the second year in Field experiment

Com_relative_bio2 = composition_change_common %>% drop_na(rebio2020)
length(unique(Com_relative_bio2$Species))

Com_relative_bio2$Field_AGB_log = log10(Com_relative_bio2$Field_AGB)
###
Com_relative_bio2$Block = as.factor(Com_relative_bio2$Block)
Com_relative_bio2$Plot_num = as.factor(Com_relative_bio2$Plot_num)

mod12 <- lmer((rebio2020_100) ~ Field_AGB_log +  (1|Block/Plot_num), data=Com_relative_bio2)
anova(mod12)

MuMIn::r.squaredGLMM(mod12)
Com_relative_bio2$F0 = predictSE(mod12, Com_relative_bio2, level = 0)$fit
Com_relative_bio2$SE <- predictSE(mod12, Com_relative_bio2, level = 0)$se.fit

AGB_lmer1 = ggplot(Com_relative_bio2, aes(x=(Field_AGB_log), y=(rebio2020_100))) +
  geom_point(aes(fill = Origin),size = 2.2, pch = 21) + # fill = "#00000022", color = "#0A0A0A"
  geom_line(aes(y=F0), size=1, linetype = 1) + #color = "#8B0000"
  scale_fill_manual(values = c('#60A7A6','#FEA6A6'))+
  scale_x_continuous(labels = scales::label_comma(accuracy =0.01)) +
  scale_y_continuous(labels = scales::label_comma(accuracy =0.01)) + 
  labs(x = ' \n ',
       y = ' \n ',
       title = "Aboveground biomass") +
  mytheme + 
  geom_line(aes(y =  F0 - 1.96 * SE), color = "grey40", linetype = "dashed", size=0.6) +
  geom_line(aes(y =  F0 + 1.96 * SE), color = "grey40", linetype = "dashed", size=0.6) +
  annotate('text',x = 1.3,y=7.5,label=('R2 = 0.320 \n p < 0.001'),size=4,color='black');AGB_lmer1

### Remove perennials
Com_relative_bio3 = subset(Com_relative_bio2, Lifeform != "Perennial")
mod12 <- lmer((rebio2020_100) ~ Field_AGB_log +  (1|Block/Plot_num), data=Com_relative_bio3)
anova(mod12)

MuMIn::r.squaredGLMM(mod12)
Com_relative_bio3$F0 = predictSE(mod12, Com_relative_bio3, level = 0)$fit
Com_relative_bio3$SE <- predictSE(mod12, Com_relative_bio3, level = 0)$se.fit

AGB_lmer2 = ggplot(Com_relative_bio3, aes(x=(Field_AGB_log), y=(rebio2020_100))) +
  geom_point(aes(fill = Origin),size = 2.2, pch = 21) + # fill = "#00000022", color = "#0A0A0A"
  geom_line(aes(y=F0), size=1, linetype = 1) + #color = "#8B0000"
  scale_fill_manual(values = c('#60A7A6','#FEA6A6'))+
  scale_x_continuous(labels = scales::label_comma(accuracy =0.01)) +
  scale_y_continuous(labels = scales::label_comma(accuracy =0.01)) + 
  labs(x = ' \n ',
       y = ' \n ',
       title = "Aboveground biomass") +
  mytheme + 
  geom_line(aes(y =  F0 - 1.96 * SE), color = "grey40", linetype = "dashed", size=0.6) +
  geom_line(aes(y =  F0 + 1.96 * SE), color = "grey40", linetype = "dashed", size=0.6) +
  annotate('text',x = 1.3,y=7.5,label=('R2 = 0.270 \n p < 0.001'),size=4,color='black');AGB_lmer2


### Evaluating the relationship between the persistence of species and the 
### functional traits of plants in the field
composition_change$exist_prob <- ifelse(!is.na(composition_change$rebio2020) & composition_change$rebio2020 > 0, 1, 0)

### Select common species
composition_change_common <- composition_change[(composition_change$Species %in% Common_sp_list_AGB), ]
length(unique(composition_change_common$Species))

###
data_a = nrow(subset(composition_change_common, Origin == "Native" & exist_prob == "1"))
data_b = nrow(subset(composition_change_common, Origin == "Native" & exist_prob == "0"))

data_c = nrow(subset(composition_change_common, Origin == "Exotic" & exist_prob == "1"))
data_d = nrow(subset(composition_change_common, Origin == "Exotic" & exist_prob == "0"))

Native = as.data.frame(c(data_a, data_b))
Exotic = as.data.frame(c(data_c, data_d))

data_per = data.frame(Native = Native, Exotic = Exotic)
data_per$Type = c("Persistence","Disappear") 
colnames(data_per)[c(1,2)] = c("Native","Exotic")

dat <- reshape2::melt(data_per, id = 'Type')

### Fig S4a
P2 <- ggplot(dat, aes(variable, value, fill = Type)) +
  geom_col(width = 0.7, position = 'fill', color = "black") +
  scale_fill_manual(values = c('#EAA37C', '#86C0DD'),
                    limits = c('Disappear', 'Persistence')) +
  labs(x = NULL, y = 'Proportion', fill = 'Status') +
  mytheme +
  theme(legend.position = "right"); P2

prop.test(c(40, 54), c(181, 223))  #Persistence
prop.test(c(141, 169), c(181, 223))  #Disappear

###
composition_change_common_exist = composition_change_common
composition_change_common_exist$Block = as.factor(composition_change_common_exist$Block)
composition_change_common_exist$Plot_num = as.factor(composition_change_common_exist$Plot_num)
composition_change_common_exist2 = subset(composition_change_common_exist, Lifeform != "Perennial")
unique(composition_change_common_exist$Species)
### 
composition_change_common_exist$Field_SLA_imp_log = log10(composition_change_common_exist$Field_SLA_imp)
composition_change_common_exist$Field_Hmax_log = log10(composition_change_common_exist$Field_Hmax)
composition_change_common_exist$Field_AGB_log = log10(composition_change_common_exist$Field_AGB)

### SLA
mod12 <- glmer(exist_prob ~ (Field_SLA_imp_log) + (1|Block/Plot_num), na.action=na.omit, family=binomial, data=composition_change_common_exist)
table(composition_change_common_exist$Origin)
Anova(mod12)
summary(mod12)
MuMIn::r.squaredGLMM(mod12)
composition_change_common_exist$F0 = predictSE(mod12, composition_change_common_exist, level = 0)$fit
composition_change_common_exist$SE <- predictSE(mod12, composition_change_common_exist, level = 0)$se.fit

SLA_glmer1 = ggplot(composition_change_common_exist, aes(x=(Field_SLA_imp_log), y=exist_prob)) +
  geom_point(aes(fill = Origin),size = 2.2, pch = 21) + # fill = "#00000022", color = "#0A0A0A"
  geom_line(aes(y=F0), size=1, linetype = 2) + #color = "#8B0000"
  labs(x = ' \n ',
       y = 'Probability of persistence',
       title = "Specific leaf area") +
  mytheme + 
  scale_fill_manual(values = c('#60A7A6','#FEA6A6'))+
  scale_x_continuous(labels = scales::label_comma(accuracy =0.01)) +
  scale_y_continuous(labels = scales::label_comma(accuracy =0.01)) + 
  geom_line(aes(y =  F0 - 1.96 * SE), color = "grey40", linetype = "dashed", size=0.6) +
  geom_line(aes(y =  F0 + 1.96 * SE), color = "grey40", linetype = "dashed", size=0.6) +
  annotate('text',x = 2.2,y=0.75,label=('Z = -1.00 \n R2 = 0.004 \n p = 0.318'),size=4,color='black');SLA_glmer1


### Remove perennials
composition_change_common_exist2 = subset(composition_change_common_exist, Lifeform != "Perennial")
unique(composition_change_common_exist2$Species)
mod12 <- glmer(exist_prob ~ (Field_SLA_imp_log) + (1|Block/Plot_num), na.action=na.omit, family=binomial, data=composition_change_common_exist2)
Anova(mod12)
summary(mod12)
MuMIn::r.squaredGLMM(mod12)
composition_change_common_exist2$F0 = predictSE(mod12, composition_change_common_exist2, level = 0)$fit
composition_change_common_exist2$SE <- predictSE(mod12, composition_change_common_exist2, level = 0)$se.fit

SLA_glmer2 = ggplot(composition_change_common_exist2, aes(x=(Field_SLA_imp_log), y=exist_prob)) +
  geom_point(aes(fill = Origin),size = 2.2, pch = 21) + # fill = "#00000022", color = "#0A0A0A"
  geom_line(aes(y=F0), size=1, linetype = 2) + #color = "#8B0000"
  labs(x = ' \n ',
       y = 'Probability of persistence',
       title = "Specific leaf area") +
  mytheme + 
  scale_fill_manual(values = c('#60A7A6','#FEA6A6'))+
  scale_x_continuous(labels = scales::label_comma(accuracy =0.01)) +
  scale_y_continuous(labels = scales::label_comma(accuracy =0.01)) + 
  geom_line(aes(y =  F0 - 1.96 * SE), color = "grey40", linetype = "dashed", size=0.6) +
  geom_line(aes(y =  F0 + 1.96 * SE), color = "grey40", linetype = "dashed", size=0.6) +
  annotate('text',x = 2.2,y=0.75,label=('Z = -0.76 \n R2 = 0.003 \n p = 0.450'),size=4,color='black'); SLA_glmer2


################################################################################
mod12 <- glmer(exist_prob ~ (Field_Hmax_log) + (1|Block/Plot_num), na.action=na.omit, family=binomial, data=composition_change_common_exist)
Anova(mod12)
summary(mod12)
MuMIn::r.squaredGLMM(mod12)
composition_change_common_exist$F0 = predictSE(mod12, composition_change_common_exist, level = 0)$fit
composition_change_common_exist$SE <- predictSE(mod12, composition_change_common_exist, level = 0)$se.fit

Hmax_glmer1 = ggplot(composition_change_common_exist, aes(x=(Field_Hmax_log), y=exist_prob)) +
  geom_point(aes(fill = Origin),size = 2.2, pch = 21) + # fill = "#00000022", color = "#0A0A0A"
  geom_line(aes(y=F0), size=1, linetype = 1) + #color = "#8B0000"
  labs(x = ' \n ',
       y = ' \n ',
       title = "Maximum height") +
  mytheme +   
  scale_fill_manual(values = c('#60A7A6','#FEA6A6'))+
  scale_x_continuous(labels = scales::label_comma(accuracy =0.01)) +
  scale_y_continuous(labels = scales::label_comma(accuracy =0.01)) + 
  geom_line(aes(y =  F0 - 1.96 * SE), color = "grey40", linetype = "dashed", size=0.6) +
  geom_line(aes(y =  F0 + 1.96 * SE), color = "grey40", linetype = "dashed", size=0.6) +
  annotate('text',x = 2.2,y=0.75,label=('Z = 3.86 \n R2 = 0.076 \n p < 0.001'),size=4,color='black')

Hmax_glmer1

### Remove perennials
mod12 <- glmer(exist_prob ~ (Field_Hmax_log) + (1|Block/Plot_num), na.action=na.omit, family=binomial, data=composition_change_common_exist2)
Anova(mod12)
summary(mod12)
MuMIn::r.squaredGLMM(mod12)
composition_change_common_exist2$F0 = predictSE(mod12, composition_change_common_exist2, level = 0)$fit
composition_change_common_exist2$SE <- predictSE(mod12, composition_change_common_exist2, level = 0)$se.fit

Hmax_glmer2 = ggplot(composition_change_common_exist2, aes(x=(Field_Hmax_log), y=exist_prob)) +
  geom_point(aes(fill = Origin),size = 2.2, pch = 21) + # fill = "#00000022", color = "#0A0A0A"
  geom_line(aes(y=F0), size=1, linetype = 1) + #color = "#8B0000"
  labs(x = ' \n ',
       y = ' \n ',
       title = "Maximum height") +
  mytheme +   
  scale_fill_manual(values = c('#60A7A6','#FEA6A6'))+
  scale_x_continuous(labels = scales::label_comma(accuracy =0.01)) +
  scale_y_continuous(labels = scales::label_comma(accuracy =0.01)) + 
  geom_line(aes(y =  F0 - 1.96 * SE), color = "grey40", linetype = "dashed", size=0.6) +
  geom_line(aes(y =  F0 + 1.96 * SE), color = "grey40", linetype = "dashed", size=0.6) +
  annotate('text',x = 2.2,y=0.75,label=('Z = 5.63 \n R2 = 0.237 \n p < 0.001'),size=4,color='black'); Hmax_glmer2

### AGB
mod12 <- glmer(exist_prob ~ (Field_AGB_log) + (1|Block/Plot_num), na.action=na.omit, family=binomial, data=composition_change_common_exist)
Anova(mod12)
summary(mod12)
MuMIn::r.squaredGLMM(mod12)
composition_change_common_exist$F0 = predictSE(mod12, composition_change_common_exist, level = 0)$fit
composition_change_common_exist$SE <- predictSE(mod12, composition_change_common_exist, level = 0)$se.fit

AGB_glmer1 = ggplot(composition_change_common_exist, aes(x=(Field_AGB_log), y=exist_prob)) +
  geom_point(aes(fill = Origin),size = 2.2, pch = 21) + # fill = "#00000022", color = "#0A0A0A"
  geom_line(aes(y=F0), size=1, linetype = 1) + #color = "#8B0000"
  labs(x = ' \n ',
       y = ' \n ',
       title = "Aboveground biomass") +
  mytheme +
  scale_fill_manual(values = c('#60A7A6','#FEA6A6'))+
  scale_x_continuous(labels = scales::label_comma(accuracy =0.01)) +
  scale_y_continuous(labels = scales::label_comma(accuracy =0.01)) + 
  geom_line(aes(y =  F0 - 1.96 * SE), color = "grey40", linetype = "dashed", size=0.6) +
  geom_line(aes(y =  F0 + 1.96 * SE), color = "grey40", linetype = "dashed", size=0.6) +
  annotate('text',x = 2.2,y=0.75,label=('Z = 4.09 \n R2 = 0.074 \n p < 0.001'),size=4,color='black');AGB_glmer1


### Remove perennials
mod12 <- glmer(exist_prob ~ (Field_AGB_log) + (1|Block/Plot_num), na.action=na.omit, family=binomial, data=composition_change_common_exist2)
Anova(mod12)
summary(mod12)
MuMIn::r.squaredGLMM(mod12)
composition_change_common_exist2$F0 = predictSE(mod12, composition_change_common_exist2, level = 0)$fit
composition_change_common_exist2$SE <- predictSE(mod12, composition_change_common_exist2, level = 0)$se.fit

AGB_glmer2 = ggplot(composition_change_common_exist2, aes(x=(Field_AGB_log), y=exist_prob)) +
  geom_point(aes(fill = Origin),size = 2.2, pch = 21) + # fill = "#00000022", color = "#0A0A0A"
  geom_line(aes(y=F0), size=1, linetype = 1) + #color = "#8B0000"
  labs(x = ' \n ',
       y = ' \n ',
       title = "Aboveground biomass") +
  mytheme +
  scale_fill_manual(values = c('#60A7A6','#FEA6A6'))+
  scale_x_continuous(labels = scales::label_comma(accuracy =0.01)) +
  scale_y_continuous(labels = scales::label_comma(accuracy =0.01)) + 
  geom_line(aes(y =  F0 - 1.96 * SE), color = "grey40", linetype = "dashed", size=0.6) +
  geom_line(aes(y =  F0 + 1.96 * SE), color = "grey40", linetype = "dashed", size=0.6) +
  annotate('text',x = 2.2,y=0.75,label=('Z = 6.34 \n R2 = 0.247 \n p < 0.001'),size=4,color='black')

AGB_glmer2

### Fig S4
P2 + P1 + plot_annotation(tag_levels = "a")

### Fig 3
SLA_glmer1+Hmax_glmer1+AGB_glmer1+SLA_lmer1+Hmax_lmer1+AGB_lmer1+
  plot_layout(ncol = 3, nrow = 2) +
  plot_annotation(tag_levels = 'a')

### Fig 3 (Remove perennials)
SLA_glmer2+Hmax_glmer2+AGB_glmer2+SLA_lmer2+Hmax_lmer2+AGB_lmer2+
  plot_layout(ncol = 3, nrow = 2) +
  plot_annotation(tag_levels = 'a')

################################################################################
##################################   Fig S3   ##################################
################################################################################
### Traits database based on pot experiment
### Select common species
composition_change_common <- composition_change[(composition_change$Species %in% Common_sp_list_AGB), ]
length(unique(composition_change_common$Species))
composition_change_common$rebio2020_100 = sqrt(composition_change_common$rebio2020*100)
composition_change_common$Block = as.factor(composition_change_common$Block)
composition_change_common$Plot_num = as.factor(composition_change_common$Plot_num)
### Remove the missing data of relative biomass in the second year
Com_relative_bio = composition_change_common %>% drop_na(rebio2020) %>% 
  left_join(pot_trait, by = "Species")

### Relationship between Specific leaf area and relative biomass 
### in the second year in pot experiment
### Inferred data
Com_relative_bio2 = Com_relative_bio %>% drop_na(SLA_green)
length(unique(Com_relative_bio2$Species))
Com_relative_bio2$SLA_green_log = log10(Com_relative_bio2$SLA_green)
###
mod12 <- lmer(rebio2020_100 ~ SLA_green_log +  (1|Block/Plot_num), data=Com_relative_bio2)
summary(mod12)
anova(mod12)
MuMIn::r.squaredGLMM(mod12)

Com_relative_bio2$F0 = predictSE(mod12, Com_relative_bio2, level = 0)$fit
Com_relative_bio2$SE <- predictSE(mod12, Com_relative_bio2, level = 0)$se.fit

SLA_lmer1 = ggplot(Com_relative_bio2, aes(x=(SLA_green_log), y=rebio2020_100)) +
  geom_point(aes(fill = Origin),size = 2.2, pch = 21) + # fill = "#00000022", color = "#0A0A0A"
  geom_line(aes(y=F0), size=1, linetype = 2) + #color = "#8B0000"
  labs(x = ' \n ',
       y = 'Relative abundance within the \n plot in second year (%, sqrt)',
       title = "Specific leaf area") +
  mytheme + 
  scale_fill_manual(values = c('#60A7A6','#FEA6A6'))+
  scale_x_continuous(labels = scales::label_comma(accuracy =0.01)) +
  scale_y_continuous(labels = scales::label_comma(accuracy =0.01)) + 
  geom_line(aes(y =  F0 - 1.96 * SE), color = "grey40", linetype = "dashed", size=0.6) +
  geom_line(aes(y =  F0 + 1.96 * SE), color = "grey40", linetype = "dashed", size=0.6) +
  annotate('text',x = 2.2,y=7.5,label=('R2 < 0.001 \n p = 0.784'),size=4,color='black'); SLA_lmer1

### Remove perennials
Com_relative_bio3 = subset(Com_relative_bio2, Lifeform != "Perennial")
mod12 <- lmer(rebio2020_100 ~ SLA_green_log +  (1|Block/Plot_num), data=Com_relative_bio3)
summary(mod12)
anova(mod12)

MuMIn::r.squaredGLMM(mod12)
Com_relative_bio3$F0 = predictSE(mod12, Com_relative_bio3, level = 0)$fit
Com_relative_bio3$SE <- predictSE(mod12, Com_relative_bio3, level = 0)$se.fit

SLA_lmer2 = ggplot(Com_relative_bio3, aes(x=(SLA_green_log), y=rebio2020_100)) +
  geom_point(aes(fill = Origin),size = 2.2, pch = 21) + # fill = "#00000022", color = "#0A0A0A"
  geom_line(aes(y=F0), size=1, linetype = 2) + #color = "#8B0000"
  labs(x = ' \n ',
       y = 'Relative abundance within the \n plot in second year (%, sqrt)',
       title = "Specific leaf area") +
  mytheme + 
  scale_fill_manual(values = c('#60A7A6','#FEA6A6'))+
  scale_x_continuous(labels = scales::label_comma(accuracy =0.01)) +
  scale_y_continuous(labels = scales::label_comma(accuracy =0.01)) + 
  geom_line(aes(y =  F0 - 1.96 * SE), color = "grey40", linetype = "dashed", size=0.6) +
  geom_line(aes(y =  F0 + 1.96 * SE), color = "grey40", linetype = "dashed", size=0.6) +
  annotate('text',x = 2.2,y=7.5,label=('R2 = 0.021 \n p = 0.215'),size=4,color='black'); SLA_lmer2

SLA_lmer2

################################################################################
### Relationship between maximum plant height and relative biomass 
### in the second year in pot experiment
Com_relative_bio2 = Com_relative_bio %>% drop_na(rebio2020)
length(unique(Com_relative_bio2$Species))
Com_relative_bio2$Hmax_green_log = log10(Com_relative_bio2$Hmax_green)

###
Com_relative_bio2$Block = as.factor(Com_relative_bio2$Block)
Com_relative_bio2$Plot_num = as.factor(Com_relative_bio2$Plot_num)

mod12 <- lmer(rebio2020_100 ~ Hmax_green_log +  (1|Block/Plot_num), data=Com_relative_bio2)
anova(mod12)

MuMIn::r.squaredGLMM(mod12)
Com_relative_bio2$F0 = predictSE(mod12, Com_relative_bio2, level = 0)$fit
Com_relative_bio2$SE <- predictSE(mod12, Com_relative_bio2, level = 0)$se.fit

Hmax_lmer1 = ggplot(Com_relative_bio2, aes(x=(Hmax_green_log), y=rebio2020_100)) +
  geom_point(aes(fill = Origin),size = 2.2, pch = 21) + # fill = "#00000022", color = "#0A0A0A"
  geom_line(aes(y=F0), size=1, linetype = 1) + #color = "#8B0000"
  labs(x = 'Traits value (log10) estimated in pot experiment',
       y = ' \n ',
       title = "Maximum height") +
  mytheme + 
  scale_fill_manual(values = c('#60A7A6','#FEA6A6'))+
  scale_x_continuous(labels = scales::label_comma(accuracy =0.01)) +
  scale_y_continuous(labels = scales::label_comma(accuracy =0.01)) + 
  geom_line(aes(y =  F0 - 1.96 * SE), color = "grey40", linetype = "dashed", size=0.6) +
  geom_line(aes(y =  F0 + 1.96 * SE), color = "grey40", linetype = "dashed", size=0.6) +
  annotate('text',x = 1.5,y=7.5,label=('R2 = 0.051 \n p = 0.024'),size=4,color='black');Hmax_lmer1

### Remove perennials
Com_relative_bio3 = subset(Com_relative_bio2, Lifeform != "Perennial")
mod12 <- lmer(rebio2020_100 ~ Hmax_green_log +  (1|Block/Plot_num), data=Com_relative_bio3)
anova(mod12)

MuMIn::r.squaredGLMM(mod12)
Com_relative_bio3$F0 = predictSE(mod12, Com_relative_bio3, level = 0)$fit
Com_relative_bio3$SE <- predictSE(mod12, Com_relative_bio3, level = 0)$se.fit

Hmax_lmer2 = ggplot(Com_relative_bio3, aes(x=(Hmax_green_log), y=rebio2020_100)) +
  geom_point(aes(fill = Origin),size = 2.2, pch = 21) + # fill = "#00000022", color = "#0A0A0A"
  geom_line(aes(y=F0), size=1, linetype = 2) + #color = "#8B0000"
  labs(x = 'Traits value (log10) estimated in pot experiment',
       y = ' \n ',
       title = "Maximum height") +
  mytheme + 
  scale_fill_manual(values = c('#60A7A6','#FEA6A6'))+
  scale_x_continuous(labels = scales::label_comma(accuracy =0.01)) +
  scale_y_continuous(labels = scales::label_comma(accuracy =0.01)) + 
  geom_line(aes(y =  F0 - 1.96 * SE), color = "grey40", linetype = "dashed", size=0.6) +
  geom_line(aes(y =  F0 + 1.96 * SE), color = "grey40", linetype = "dashed", size=0.6) +
  annotate('text',x = 1.5,y=7.5,label=('R2 = 0.019 \n p = 0.226'),size=4,color='black'); Hmax_lmer2

### Relationship between Aboveground biomass and relative biomass 
### in the second year in pot experiment
Com_relative_bio2$AGB_green_log = log10(Com_relative_bio2$AGB_green)
unique(Com_relative_bio2$Species)
###
Com_relative_bio2$Block = as.factor(Com_relative_bio2$Block)
Com_relative_bio2$Plot_num = as.factor(Com_relative_bio2$Plot_num)

mod12 <- lmer((rebio2020_100) ~ AGB_green_log +  (1|Block/Plot_num), data=Com_relative_bio2)
anova(mod12)
MuMIn::r.squaredGLMM(mod12)

Com_relative_bio2$F0 = predictSE(mod12, Com_relative_bio2, level = 0)$fit
Com_relative_bio2$SE <- predictSE(mod12, Com_relative_bio2, level = 0)$se.fit

AGB_lmer1 = ggplot(Com_relative_bio2, aes(x=(AGB_green_log), y=(rebio2020_100))) +
  geom_point(aes(fill = Origin),size = 2.2, pch = 21) + # fill = "#00000022", color = "#0A0A0A"
  geom_line(aes(y=F0), size=1, linetype = 2) + #color = "#8B0000"
  scale_fill_manual(values = c('#60A7A6','#FEA6A6'))+
  scale_x_continuous(labels = scales::label_comma(accuracy =0.01)) +
  scale_y_continuous(labels = scales::label_comma(accuracy =0.01)) + 
  labs(x = ' \n ',
       y = ' \n ',
       title = "Aboveground biomass") +
  mytheme + 
  geom_line(aes(y =  F0 - 1.96 * SE), color = "grey40", linetype = "dashed", size=0.6) +
  geom_line(aes(y =  F0 + 1.96 * SE), color = "grey40", linetype = "dashed", size=0.6) +
  annotate('text',x = 0,y=7.5,label=('R2 = 0.010 \n p = 0.325'),size=4,color='black'); AGB_lmer1

### Remove perennials
Com_relative_bio3 = subset(Com_relative_bio2, Lifeform != "Perennial")
mod12 <- lmer((rebio2020_100) ~ AGB_green_log +  (1|Block/Plot_num), data=Com_relative_bio3)
anova(mod12)
MuMIn::r.squaredGLMM(mod12)

Com_relative_bio3$F0 = predictSE(mod12, Com_relative_bio3, level = 0)$fit
Com_relative_bio3$SE <- predictSE(mod12, Com_relative_bio3, level = 0)$se.fit

AGB_lmer2 = ggplot(Com_relative_bio3, aes(x=(AGB_green_log), y=(rebio2020_100))) +
  geom_point(aes(fill = Origin),size = 2.2, pch = 21) + # fill = "#00000022", color = "#0A0A0A"
  geom_line(aes(y=F0), size=1, linetype = 2) + #color = "#8B0000"
  scale_fill_manual(values = c('#60A7A6','#FEA6A6'))+
  scale_x_continuous(labels = scales::label_comma(accuracy =0.01)) +
  scale_y_continuous(labels = scales::label_comma(accuracy =0.01)) + 
  labs(x = ' \n ',
       y = ' \n ',
       title = "Aboveground biomass") +
  mytheme + 
  geom_line(aes(y =  F0 - 1.96 * SE), color = "grey40", linetype = "dashed", size=0.6) +
  geom_line(aes(y =  F0 + 1.96 * SE), color = "grey40", linetype = "dashed", size=0.6) +
  annotate('text',x = 0,y=7.5,label=('R2 = 0.003 \n p = 0.647'),size=4,color='black'); AGB_lmer2

################################################################################
# Evaluate the relationship between the persistence of species and plant functional traits
composition_change_common$exist_prob <- ifelse(!is.na(composition_change_common$rebio2020) & composition_change_common$rebio2020 > 0, 1, 0)
composition_change_common_exist = composition_change_common %>% left_join(pot_trait, by = "Species") %>% drop_na(SLA_green)
length(unique(composition_change_common_exist$Species))

composition_change_common_exist$SLA_green_log = log10(composition_change_common_exist$SLA_green)
composition_change_common_exist$Hmax_green_log = log10(composition_change_common_exist$Hmax_green)
composition_change_common_exist$AGB_green_log = log10(composition_change_common_exist$AGB_green)

### SLA
library(lme4)
library(AICcmodavg)
mod12 <- glmer(exist_prob ~ (SLA_green_log) + (1|Block/Plot_num), na.action=na.omit, family=binomial, data=composition_change_common_exist)
table(composition_change_common_exist$Origin)
Anova(mod12)
summary(mod12)
MuMIn::r.squaredGLMM(mod12)
composition_change_common_exist$F0 = predictSE(mod12, composition_change_common_exist, level = 0)$fit
composition_change_common_exist$SE <- predictSE(mod12, composition_change_common_exist, level = 0)$se.fit

SLA_glmer1 = ggplot(composition_change_common_exist, aes(x=(SLA_green_log), y=exist_prob)) +
  geom_point(aes(fill = Origin),size = 2.2, pch = 21) + # fill = "#00000022", color = "#0A0A0A"
  geom_line(aes(y=F0), size=1, linetype = 1) + #color = "#8B0000"
  labs(x = ' \n ',
       y = 'Probability of persistence',
       title = "Specific leaf area") +
  mytheme + 
  scale_fill_manual(values = c('#60A7A6','#FEA6A6'))+
  scale_x_continuous(labels = scales::label_comma(accuracy =0.01)) +
  scale_y_continuous(labels = scales::label_comma(accuracy =0.01)) + 
  geom_line(aes(y =  F0 - 1.96 * SE), color = "grey40", linetype = "dashed", size=0.6) +
  geom_line(aes(y =  F0 + 1.96 * SE), color = "grey40", linetype = "dashed", size=0.6) +
  annotate('text',x = 2.2,y=0.75,label=('Z = -2.68 \n R2 = 0.027 \n p = 0.007'),size=4,color='black')

SLA_glmer1

### Remove perennials
composition_change_common_exist2 = subset(composition_change_common_exist, Lifeform != "Perennial")
mod12 <- glmer(exist_prob ~ (SLA_green_log) + (1|Block/Plot_num), na.action=na.omit, family=binomial, data=composition_change_common_exist2)
Anova(mod12)
summary(mod12)
MuMIn::r.squaredGLMM(mod12)
composition_change_common_exist2$F0 = predictSE(mod12, composition_change_common_exist2, level = 0)$fit
composition_change_common_exist2$SE <- predictSE(mod12, composition_change_common_exist2, level = 0)$se.fit

SLA_glmer2 = ggplot(composition_change_common_exist2, aes(x=(SLA_green_log), y=exist_prob)) +
  geom_point(aes(fill = Origin),size = 2.2, pch = 21) + # fill = "#00000022", color = "#0A0A0A"
  geom_line(aes(y=F0), size=1, linetype = 1) + #color = "#8B0000"
  labs(x = ' \n ',
       y = 'Probability of persistence',
       title = "Specific leaf area") +
  mytheme + 
  scale_fill_manual(values = c('#60A7A6','#FEA6A6'))+
  scale_x_continuous(labels = scales::label_comma(accuracy =0.01)) +
  scale_y_continuous(labels = scales::label_comma(accuracy =0.01)) + 
  geom_line(aes(y =  F0 - 1.96 * SE), color = "grey40", linetype = "dashed", size=0.6) +
  geom_line(aes(y =  F0 + 1.96 * SE), color = "grey40", linetype = "dashed", size=0.6) +
  annotate('text',x = 2.2,y=0.75,label=('Z = -3.50 \n R2 = 0.052 \n p < 0.001'),size=4,color='black')

SLA_glmer2

### Hmax
mod12 <- glmer(exist_prob ~ (Hmax_green_log) + (1|Block/Plot_num), na.action=na.omit, family=binomial, data=composition_change_common_exist)
Anova(mod12)
summary(mod12)
MuMIn::r.squaredGLMM(mod12)
composition_change_common_exist$F0 = predictSE(mod12, composition_change_common_exist, level = 0)$fit
composition_change_common_exist$SE <- predictSE(mod12, composition_change_common_exist, level = 0)$se.fit

Hmax_glmer1 = ggplot(composition_change_common_exist, aes(x=(Hmax_green_log), y=exist_prob)) +
  geom_point(aes(fill = Origin),size = 2.2, pch = 21) + # fill = "#00000022", color = "#0A0A0A"
  geom_line(aes(y=F0), size=1, linetype = 2) + #color = "#8B0000"
  labs(x = ' \n ',
       y = ' \n ',
       title = "Maximum height") +
  mytheme +   
  scale_fill_manual(values = c('#60A7A6','#FEA6A6'))+
  scale_x_continuous(labels = scales::label_comma(accuracy =0.01)) +
  scale_y_continuous(labels = scales::label_comma(accuracy =0.01)) + 
  geom_line(aes(y =  F0 - 1.96 * SE), color = "grey40", linetype = "dashed", size=0.6) +
  geom_line(aes(y =  F0 + 1.96 * SE), color = "grey40", linetype = "dashed", size=0.6) +
  annotate('text',x = 1.5,y=0.75,label=('Z = 1.06 \n R2 = 0.005 \n p = 0.290'),size=4,color='black')

Hmax_glmer1

### Remove perennials
mod12 <- glmer(exist_prob ~ (Hmax_green_log) + (1|Block/Plot_num), na.action=na.omit, family=binomial, data=composition_change_common_exist2)
Anova(mod12)
summary(mod12)
MuMIn::r.squaredGLMM(mod12)
composition_change_common_exist2$F0 = predictSE(mod12, composition_change_common_exist2, level = 0)$fit
composition_change_common_exist2$SE <- predictSE(mod12, composition_change_common_exist2, level = 0)$se.fit

Hmax_glmer2 = ggplot(composition_change_common_exist2, aes(x=(Hmax_green_log), y=exist_prob)) +
  geom_point(aes(fill = Origin),size = 2.2, pch = 21) + # fill = "#00000022", color = "#0A0A0A"
  geom_line(aes(y=F0), size=1, linetype = 2) + #color = "#8B0000"
  labs(x = ' \n ',
       y = ' \n ',
       title = "Maximum height") +
  mytheme +   
  scale_fill_manual(values = c('#60A7A6','#FEA6A6'))+
  scale_x_continuous(labels = scales::label_comma(accuracy =0.01)) +
  scale_y_continuous(labels = scales::label_comma(accuracy =0.01)) + 
  geom_line(aes(y =  F0 - 1.96 * SE), color = "grey40", linetype = "dashed", size=0.6) +
  geom_line(aes(y =  F0 + 1.96 * SE), color = "grey40", linetype = "dashed", size=0.6) +
  annotate('text',x = 1.2,y=0.75,label=('Z = 1.93 \n R2 = 0.022 \n p = 0.053'),size=4,color='black')

Hmax_glmer2

### AGB
unique(composition_change_common_exist$Species)
mod12 <- glmer(exist_prob ~ (AGB_green_log) + (1|Block/Plot_num), na.action=na.omit, family=binomial, data=composition_change_common_exist)
Anova(mod12)
summary(mod12)
MuMIn::r.squaredGLMM(mod12)
composition_change_common_exist$F0 = predictSE(mod12, composition_change_common_exist, level = 0)$fit
composition_change_common_exist$SE <- predictSE(mod12, composition_change_common_exist, level = 0)$se.fit

AGB_glmer1 = ggplot(composition_change_common_exist, aes(x=(AGB_green_log), y=exist_prob)) +
  geom_point(aes(fill = Origin),size = 2.2, pch = 21) + # fill = "#00000022", color = "#0A0A0A"
  geom_line(aes(y=F0), size=1, linetype = 1) + #color = "#8B0000"
  labs(x = ' \n ',
       y = ' \n ',
       title = "Aboveground biomass") +
  mytheme +
  scale_fill_manual(values = c('#60A7A6','#FEA6A6'))+
  scale_x_continuous(labels = scales::label_comma(accuracy =0.01)) +
  scale_y_continuous(labels = scales::label_comma(accuracy =0.01)) + 
  geom_line(aes(y =  F0 - 1.96 * SE), color = "grey40", linetype = "dashed", size=0.6) +
  geom_line(aes(y =  F0 + 1.96 * SE), color = "grey40", linetype = "dashed", size=0.6) +
  annotate('text',x = 0,y=0.75,label=('Z = 2.12 \n R2 = 0.021 \n p = 0.034'),size=4,color='black')

AGB_glmer1

### Remove perennials
mod12 <- glmer(exist_prob ~ (AGB_green_log) + (1|Block/Plot_num), na.action=na.omit, family=binomial, data=composition_change_common_exist2)
Anova(mod12)
summary(mod12)
MuMIn::r.squaredGLMM(mod12)
composition_change_common_exist2$F0 = predictSE(mod12, composition_change_common_exist2, level = 0)$fit
composition_change_common_exist2$SE <- predictSE(mod12, composition_change_common_exist2, level = 0)$se.fit

AGB_glmer2 = ggplot(composition_change_common_exist2, aes(x=(AGB_green_log), y=exist_prob)) +
  geom_point(aes(fill = Origin),size = 2.2, pch = 21) + # fill = "#00000022", color = "#0A0A0A"
  geom_line(aes(y=F0), size=1, linetype = 1) + #color = "#8B0000"
  labs(x = ' \n ',
       y = ' \n ',
       title = "Aboveground biomass") +
  mytheme +
  scale_fill_manual(values = c('#60A7A6','#FEA6A6'))+
  scale_x_continuous(labels = scales::label_comma(accuracy =0.01)) +
  scale_y_continuous(labels = scales::label_comma(accuracy =0.01)) + 
  geom_line(aes(y =  F0 - 1.96 * SE), color = "grey40", linetype = "dashed", size=0.6) +
  geom_line(aes(y =  F0 + 1.96 * SE), color = "grey40", linetype = "dashed", size=0.6) +
  annotate('text',x = 0,y=0.75,label=('Z = 3.15 \n R2 = 0.067 \n p = 0.002'),size=4,color='black')

AGB_glmer2

### Fig S3
SLA_glmer1+Hmax_glmer1+AGB_glmer1+SLA_lmer1+Hmax_lmer1+AGB_lmer1+
  plot_layout(ncol = 3, nrow = 2) +
  plot_annotation(tag_levels = 'a')

### Fig S3 (Remove perennials)
SLA_glmer2+Hmax_glmer2+AGB_glmer2+SLA_lmer2+Hmax_lmer2+AGB_lmer2+
  plot_layout(ncol = 3, nrow = 2) +
  plot_annotation(tag_levels = 'a')

################################################################################
#################################   Fig S1   ###################################
################################################################################
### Sensitivity analysis
### Traits deviation estimation
### Loading pot experiment database
pot_trait = read.xlsx("Data/Pot_traits_mean.xlsx", sheet = "Pot_means", colNames = TRUE, rowNames = FALSE)
colnames(pot_trait) <- paste0(colnames(pot_trait), "_green")
colnames(pot_trait)[1] <- "Species"
pot_trait = pot_trait[pot_trait$Species %in% unique(Common_sp_list_SLA), ]
unique(pot_trait$Species)

### Loading field experiment database
trait_data = read.xlsx("Data/Field_traits_mean.xlsx", sheet = "Field_means", colNames = TRUE, rowNames = FALSE)
colnames(trait_data) <- paste0(colnames(trait_data), "_field")
colnames(trait_data)[1] <- "Species"
trait_data = trait_data[trait_data$Species %in% unique(Common_sp_list_SLA), ]
unique(trait_data$Species)

trait_data = trait_data %>% left_join(pot_trait, by = "Species") %>%
  left_join(Common_sp_list[,c(1,2,5)], by= "Species")
length(unique(trait_data$Species))

trait_data$Origin = factor(trait_data$Origin, levels = c("Native", "Exotic"))
summary(lm(log10(SLA_field)~ log10(SLA_green), data = trait_data))
ma.test <- sma(SLA_field ~ SLA_green, log='xy', slope.test=1, data = trait_data, type = "shift", na.action = na.omit)
ma.test$n
summary(ma.test)

df2 = (trait_data %>% arrange(desc(SLA_field)))[c(1:6),]
df3 = (trait_data %>% arrange((SLA_field)))[c(1:6),]

p1 = ggplot(trait_data, mapping = aes(x = log10(SLA_green) , y = log10(SLA_field))) + 
  geom_point(trait_data, mapping = aes(x = log10(SLA_green) , y = log10(SLA_field), shape = Origin, fill = Origin),size=2.2, color = "black")+
  scale_shape_manual(values = c(21,21))+
  scale_color_manual(values=c('#60A7A6','#FEA6A6'))+
  scale_fill_manual(values=c('#60A7A6','#FEA6A6'))+
  scale_x_continuous(labels = scales::label_comma(accuracy =0.1)) +
  scale_y_continuous(labels = scales::label_comma(accuracy =0.1)) + 
  labs(x = 'Specific leaf area (cm2/g, log10) \nestimated in pot experiment',
       y = 'Specific leaf area (cm2/g, log10) \nestimated in field experiment') + 
  geom_abline(intercept=1.353949 ,slope=0.4013532 ,size=.8, linetype = 1)+
  mytheme + 
  annotate('text',x = 1.8,y=2.2,label=('R2 = 0.083 \n p = 0.048'),size=4,color='black') +
  ggrepel::geom_text_repel(aes(label=Species), df2,size = 2.5,segment.color = "black", color = "black",direction = "both",box.padding = 0.7,
                           max.overlaps = getOption("ggrepel.max.overlaps", default = 100)) +
  ggrepel::geom_text_repel(aes(label=Species), df3,size = 2.5,segment.color = "black", color = "black",direction = "both",box.padding = 0.7,
                           max.overlaps = getOption("ggrepel.max.overlaps", default = 100))


p1

### DI
library(funrar)
standr = function(x){(x-min(x))/(max(x)-min(x))} 
r2 <- function(x) {x$`Sum Sq`[1]/ (x$`Sum Sq`[1] + x$`Sum Sq`[2])}
trait_data = read.xlsx("Data/Field_traits_mean.xlsx", sheet = "Field_means", colNames = TRUE, rowNames = FALSE)
colnames(trait_data) <- paste0(colnames(trait_data), "_field")
colnames(trait_data)[1] <- "Species"

trait_data = trait_data %>% left_join(pot_trait, by = "Species") %>% drop_na(SLA_green) %>% drop_na(SLA_field)
rownames(trait_data) = trait_data$Species
length(unique(trait_data$Species))

trait_data$SLA_field_log = log10(trait_data$SLA_field)
trait_data$SLA_green_log = log10(trait_data$SLA_green)

### Pot experiment SLA differences
dist_mat <- compute_dist_matrix(trait_data['SLA_green_log'], metric = 'euclidean',center = TRUE, scale = TRUE)
dist_mat <- standr(dist_mat)
diag(dist_mat) <- NA
### Calculate the mean of each row (excluding diagonals)
mean_vals <- apply(dist_mat, 1, mean, na.rm = TRUE)
### Calculate the standard deviation for each row (excluding diagonals)
sd_vals <- apply(dist_mat, 1, sd, na.rm = TRUE)
### Calculate the standard error for each line (removing diagonal lines)
se_vals <- sd_vals / sqrt(ncol(dist_mat) - 1)
SLA_green_dis <- data.frame(SLA_green_mean = mean_vals, SLA_green_sd = sd_vals, SLA_green_se = se_vals)

### Field SLA differences
dist_mat <- compute_dist_matrix(trait_data['SLA_field_log'], metric = 'euclidean',center = TRUE, scale = TRUE)#gower
dist_mat <- standr(dist_mat)
diag(dist_mat) <- NA
### Calculate the mean of each row (excluding diagonals)
mean_vals <- apply(dist_mat, 1, mean, na.rm = TRUE)
### Calculate the standard deviation for each row (excluding diagonals) 
sd_vals <- apply(dist_mat, 1, sd, na.rm = TRUE)
### Calculate the standard error for each line (removing diagonal lines) 
se_vals <- sd_vals / sqrt(ncol(dist_mat) - 1)
SLA_field_dis <- data.frame(SLA_field_mean = mean_vals, SLA_field_sd = sd_vals, SLA_field_se = se_vals)
SLA_total_dis = cbind(SLA_green_dis, SLA_field_dis)
SLA_total_dis$Species = rownames(SLA_total_dis)
SLA_total_dis = SLA_total_dis %>% left_join(Common_sp_list[,c(1,2,5)], by = "Species")
#knitr::kable(head(SLA_total_dis))

#### Deming regression
fit <- deming(SLA_field_mean ~ SLA_green_mean, ystd=SLA_field_sd, xstd=SLA_green_sd, data=SLA_total_dis)
print(fit)
SLA_total_dis$Origin = factor(SLA_total_dis$Origin, levels = c("Native","Exotic"))

### 
phylogenyAux = read.tree("Data/All_species.newick")
to_drop = phylogenyAux$tip.label[!phylogenyAux$tip.label %in% trait_data$Species]
tree <- drop.tip(as.phylo(phylogenyAux), to_drop) 

Test_total_dis = SLA_total_dis
rownames(Test_total_dis) = Test_total_dis$Species
###
SLA_pot_di<-Test_total_dis$SLA_green_mean
names(SLA_pot_di)<- rownames(Test_total_dis)
phytools::phylosig(tree, SLA_pot_di, method = "K", test = TRUE, nsim =  1000)

SLA_field_di<-Test_total_dis$SLA_field_mean
names(SLA_field_di)<- rownames(Test_total_dis)
phytools::phylosig(tree, SLA_field_di, method = "K", test = TRUE, nsim =  1000)

###
Origin<-Test_total_dis$Origin
names(Origin)<- rownames(Test_total_dis)

### phylANOVA
aov_di <- phylANOVA(tree, Origin, SLA_pot_di, nsim = 1000, posthoc=TRUE, p.adj = 'bonferroni')
aov_di
r2(aov_di)
aov_di <- phylANOVA(tree, Origin, SLA_field_di, nsim = 1000, posthoc=TRUE, p.adj = 'bonferroni')
aov_di
r2(aov_di)

### 
data1 = as.data.frame(SLA_pot_di)
data1$Type = "pot"
data1$Species = rownames(data1)
colnames(data1)[1] = "Di"
rownames(data1) = NULL

data2 = as.data.frame(SLA_field_di)
data2$Type = "field"
data2$Species = rownames(data2)
colnames(data2)[1] = "Di"
rownames(data2) = NULL
##
data_all_di = rbind(data1,data2) 
data_all_di = data_all_di %>% left_join(Common_sp_list[,c(1,2,5)], by = "Species")
data_all_di$Origin = factor(data_all_di$Origin, levels = c("Native", "Exotic"))
data_all_di$Type = factor(data_all_di$Type, levels = c("pot", "field"))
data_all_di$Species = as.factor(data_all_di$Species)
###
phyloMat = vcv.phylo(tree)
phyloMat = phyloMat / max(phyloMat)
dim(phyloMat)
###
mod = relmatLmer(Di ~ Origin*Type + (1|Species), data = data_all_di, relmat = list(Species=phyloMat))
Anova(mod, type="II", test.statistic ="F")

###
colnames(SLA_total_dis)
dat1 = Rmisc::summarySE(SLA_total_dis, groupvars = c("Origin"), measurevar = c("SLA_green_mean"))[,c(1,3,5)]
colnames(dat1)[3] = "pot_se"
dat2 = Rmisc::summarySE(SLA_total_dis, groupvars = c("Origin"), measurevar = c("SLA_field_mean"))[,c(3,5)]
colnames(dat2)[2] = "field_se"
dattt = cbind(dat1, dat2)

df2 = (SLA_total_dis %>% arrange(desc(SLA_field_mean)))[c(1:12),]

p2 = ggplot(SLA_total_dis,aes(x=SLA_green_mean,y=SLA_field_mean,fill = Origin, color = Origin, shape = Origin))+
  labs(x = "Specific leaf area(cm2/g) \ninterspecific differences in monoculture",
       y="Specific leaf area(cm2/g) \nInterspecific differences in field")+
  geom_point(size=2.2, pch = 21, color = "black")+
  geom_errorbar(data = dattt,mapping = aes(ymax = SLA_field_mean+field_se, ymin=SLA_field_mean-field_se),width=0.01,size=0.5,alpha = 1, color = "black")+#
  geom_errorbarh(data = dattt,mapping = aes(xmax=SLA_green_mean+pot_se,xmin=SLA_green_mean-pot_se),height=0.01,size=0.5,alpha = 1, color = "black")+#
  geom_point(data = dattt,mapping = aes(x = SLA_green_mean, y = SLA_field_mean),size=3.8, pch = 21, color = "black")+
  scale_shape_manual(values = c(21,21))+
  scale_color_manual(values=c('#60A7A6','#FEA6A6'))+
  scale_fill_manual(values=c('#60A7A6','#FEA6A6'))+
  scale_x_continuous(labels = scales::label_comma(accuracy =0.1)) +
  scale_y_continuous(labels = scales::label_comma(accuracy =0.1)) + 
  guides(col = guide_legend(ncol = 1))+
  #geom_abline(intercept=0.97990531, slope=0.06054834 ,size=1, linetype = 2)+
  geom_abline(intercept=-2.72938, slope=15.70254 ,size=1, linetype = 2)+
  geom_abline(intercept=0,slope=1 ,size=1, linetype = 1, color = "#95373B")+
  mytheme + 
  #annotate('text',x = 3,y=1,label=('Slope = 0.06054834 , \n95% CI [-55.23847, 55.35957]'),size=4,color='black')+
  ggrepel::geom_text_repel(aes(label=Species), df2,size = 2.5,segment.color = "black", color = "black",direction = "both",box.padding = 0.6,
                           max.overlaps = getOption("ggrepel.max.overlaps", default = 100))

p2

### Relative biomass in the second year
composition_change_common <- composition_change[(composition_change$Species %in% Common_sp_list_AGB), ]
length(unique(composition_change_common$Species))
composition_change_common$rebio2020_100 = sqrt(composition_change_common$rebio2020*100)
### Remove missing field specific leaf area data
composition_change_common = composition_change_common %>% drop_na(Field_SLA)
length(unique(composition_change_common$Species))
### Remove the missing data of relative biomass in the second year
Com_relative_bio = composition_change_common %>% drop_na(rebio2020)
length(unique(Com_relative_bio$Species))
### Relationship between specific leaf area and relative biomass of plants 
### in the second year in field experiments
### row_data
Com_relative_bio2 = Com_relative_bio %>% drop_na(Field_SLA)
length(unique(Com_relative_bio2$Species))
Com_relative_bio2$Field_SLA = log10(Com_relative_bio2$Field_SLA)

###
mod12 <- lmer((rebio2020_100) ~ Field_SLA +  (1|Block/Plot_num), data=Com_relative_bio2)
anova(mod12)
MuMIn::r.squaredGLMM(mod12)

Com_relative_bio2$F0 = predictSE(mod12, Com_relative_bio2, level = 0)$fit
Com_relative_bio2$SE <- predictSE(mod12, Com_relative_bio2, level = 0)$se.fit

p33 = ggplot(Com_relative_bio2, aes(x=(Field_SLA), y=(rebio2020_100))) +
  geom_point(aes(fill = Origin),size = 2.2, pch = 21) + # fill = "#00000022", color = "#0A0A0A"
  geom_line(aes(y=F0), size=1, linetype = 2) + #color = "#8B0000"
  mytheme + 
  scale_fill_manual(values=c('#60A7A6','#FEA6A6'))+
  scale_x_continuous(labels = scales::label_comma(accuracy =0.1)) +
  scale_y_continuous(labels = scales::label_comma(accuracy =0.1)) + 
  ylab('Relative biomass within the \n plot in second year (%, sqrt)') + 
  xlab('Specific leaf area (cm2/g, log10) \n measured in field experiment') +
  geom_line(aes(y =  F0 - 1.96 * SE), color = "grey40", linetype = "dashed", size=0.6) +
  geom_line(aes(y =  F0 + 1.96 * SE), color = "grey40", linetype = "dashed", size=0.6) +
  annotate('text',x = 2.2,y=7.5,label=('R2 = 0.002 \n p = 0.665'),size=4,color='black')

p33

### Evaluating the relationship between the persistence of species and 
### the functional traits of plants in the field
composition_change$exist_prob <- ifelse(!is.na(composition_change$rebio2020) & composition_change$rebio2020 > 0, 1, 0)

### Select common species
composition_change_common_exist <- composition_change[(composition_change$Species %in% Common_sp_list_AGB), ]
composition_change_common_exist <- composition_change_common_exist %>% drop_na(Field_SLA)
length(unique(composition_change_common_exist$Species))

### 
composition_change_common_exist$Field_SLA_log = log10(composition_change_common_exist$Field_SLA)
### SLA
composition_change_common_exist$Block = as.factor(composition_change_common_exist$Block)
composition_change_common_exist$Plot_num = as.factor(composition_change_common_exist$Plot_num)

mod12 <- glmer(exist_prob ~ (Field_SLA_log) + (1|Block/Plot_num), na.action=na.omit, family=binomial, data=composition_change_common_exist)
table(composition_change_common_exist$Origin)
Anova(mod12)
summary(mod12)
MuMIn::r.squaredGLMM(mod12)
composition_change_common_exist$F0 = predictSE(mod12, composition_change_common_exist, level = 0)$fit
composition_change_common_exist$SE <- predictSE(mod12, composition_change_common_exist, level = 0)$se.fit

p44 = ggplot(composition_change_common_exist, aes(x=(Field_SLA_log), y=exist_prob)) +
  geom_point(aes(fill = Origin),size = 2.2, pch = 21) + # fill = "#00000022", color = "#0A0A0A"
  geom_line(aes(y=F0), size=1, linetype = 2) + #color = "#8B0000"
  ylab('Probability of exist \n in second year') + 
  xlab('Specific leaf area (cm2/g, log10) \n in field experiment') +
  mytheme +   
  scale_fill_manual(values=c('#60A7A6','#FEA6A6'))+
  scale_x_continuous(labels = scales::label_comma(accuracy =0.01)) + 
  geom_line(aes(y =  F0 - 1.96 * SE), color = "grey40", linetype = "dashed", size=0.6) +
  geom_line(aes(y =  F0 + 1.96 * SE), color = "grey40", linetype = "dashed", size=0.6) +
  annotate('text',x = 2.2,y=0.75,label=('Z = -1.25 \n R2 = 0.008 \n p = 0.212'),size=4,color='black')

p44

### Fig S1
p1+p2+p44+p33+
  plot_layout(ncol = 2, nrow = 2) +
  plot_annotation(tag_levels = 'a')

### Notice that,
### For more picture details, we have further adjusted it in Adobe illustrator.
