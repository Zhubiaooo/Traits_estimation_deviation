###Loading Packages
library(openxlsx)
library(phytools)
library(ape)
library(ggplot2)
library(lsr)
library(adephylo)
library(multcomp)
library(treeplyr)
library(Rmisc)
library(lme4qtl)
library(car)
library(ggpubr)
library(patchwork)
library(tidyr)
library(lmerTest)
library(Rmisc)

### Analyzing species lists (SLA = 48, Hmax&AGB = 64)
Common_sp_list = read.xlsx("Data/Common_species_list.xlsx", sheet = "Common_sp_list", colNames = TRUE, rowNames = FALSE)
Common_sp_list_SLA = c(na.omit(Common_sp_list$Species_SLA))
Common_sp_list_AGB = c(na.omit(Common_sp_list$Species_AGB))

### Reading into the plant phylogenetic tree
phylogenyAux = read.tree("Data/iq_tree.treefile")

to_drop = phylogenyAux$tip.label[!phylogenyAux$tip.label %in% Common_sp_list$Species]
phylogenyAux <- drop.tip(as.phylo(phylogenyAux), to_drop) 
### Species phylogenetic correlation matrix
phyloMat = vcv.phylo(phylogenyAux)
phyloMat = phyloMat / max(phyloMat)
dim(phyloMat)

### Set the picture theme
mytheme = theme_bw()+
  theme( panel.background = element_rect(fill='white', colour='black'),
         panel.grid=element_blank(), 
         axis.ticks.length = unit(0.4,"lines"), 
         axis.ticks = element_line(color='black'),
         axis.line = element_line(colour = "black"), 
         axis.title.x=element_text(colour='black', size=13,vjust = 1),
         axis.title.y=element_text(colour='black', size=13,vjust = 1),
         axis.text=element_text(colour='black',size=11))+
  theme(legend.position = "none")+
  theme(plot.title = element_text(color = "black", size = 13, hjust = 0.5))

### Comparison of the first and second measurements of plant height in pot experiment
### first measurements
trait_data <- read.xlsx("Data/Pot_traits_database0831.xlsx", sheet = "Pot_data", colNames = TRUE, rowNames = FALSE)
trait_data <- trait_data[trait_data$Species %in% unique(Common_sp_list_AGB), ]
trait_data <- trait_data %>% drop_na(Hmax_time1)
length(unique(trait_data$Species))
colnames(trait_data)
#shapiro.test(trait_data$Hmax)
shapiro.test(sqrt(trait_data$Hmax_time1))
shapiro.test(sqrt(trait_data$Hmax_time2))
trait_data$Origin = factor(trait_data$Origin, levels = c("Native", "Exotic"))
colnames(trait_data)
#data11 = Rmisc::summarySE(trait_data, groupvars = c("Species"), measurevar = "Hmax_time1")
#data22 = Rmisc::summarySE(trait_data, groupvars = c("Species"), measurevar = "Hmax_time2")
#data_all = data11[c(1,3)] %>% left_join(data22[,c(1,3)])
#write.csv(data_all, "data_all.csv")
mod <- relmatLmer(sqrt(Hmax_time1) ~ Origin + (1|Species), data = trait_data, relmat = list(Species=phyloMat))
Anova(mod, type="II", test.statistic ="Chisq")
shapiro.test(residuals(mod))
#ggpubr::ggqqplot(residuals(mod))
hist(residuals(mod))

mod <- relmatLmer(sqrt(Hmax_time2) ~ Origin + (1|Species), data = trait_data, relmat = list(Species=phyloMat))
Anova(mod, type="II", test.statistic ="Chisq")
shapiro.test(residuals(mod))
#ggpubr::ggqqplot(residuals(mod))
hist(residuals(mod))

tuk_quasipoisson <- glht(mod, alternative = 'two.sided', linfct = mcp(Origin = 'Tukey'))
summary(tuk_quasipoisson) 
tuk_quasipoisson.cld <- cld(tuk_quasipoisson, level = 0.05, decreasing = TRUE)
plot(tuk_quasipoisson.cld, col = c('#60A7A6','#FEA6A6'))

###Percentage increase
summarySE(trait_data, measurevar = c("Hmax_time1"), groupvars = c("Origin"))
summarySE(trait_data, measurevar = c("Hmax_time2"), groupvars = c("Origin"))
summarySE(trait_data, measurevar = c("Hmax"), groupvars = c("Origin"))

###Phylogenetic signal test of plant traits
summary_data <- summarySE(trait_data, measurevar = "Hmax", groupvars=c("Species"))
test <- make.treedata(tree = phylogenyAux, data = summary_data[,c(1,3)], name_column = "Species")
myTree.withBrLe <- compute.brlen(test$phy)
myProx <- vcv.phylo(myTree.withBrLe)
set.seed(4321)
Physignal <- abouheif.moran(test$dat, W = myProx, method = "oriAbouheif", nrepet = 999)
Physignal
plot(Physignal)

#### plot
trait_data$Hmax = sqrt(trait_data$Hmax)
trait_data$Hmax_time1 = sqrt(trait_data$Hmax_time1)
trait_data$Hmax_time2 = sqrt(trait_data$Hmax_time2)

#### Plant height of first measurements (cm, sqrt)
ggplot(trait_data,aes(x=Origin,y=Hmax_time1))+
  geom_violin(data=trait_data, aes(y=Hmax_time1,x=Origin,fill=Origin,color=Origin),trim=T,
              scale = "width",position = position_dodge(0.5),width=0.6,alpha = 0.8, size = 1.5)+
  geom_boxplot(data=trait_data, aes(y=Hmax_time1,x=Origin),fill="white",color="black",
               position = position_dodge(0.5),width=0.15,outlier.size=0.8,outlier.shape = 1)+
  stat_summary(data=trait_data, aes(y=Hmax_time1,x=Origin),fill="white",color="black",
               fun="mean",position = position_dodge(0.5),geom="point",shape=21, size=1.5)+
  scale_color_manual(values = c('#60A7A6','#FEA6A6'))+
  scale_fill_manual(values = c('#60A7A6','#FEA6A6'))+
  mytheme+
  theme(axis.text.x = element_text(angle = 25,size=11,colour = "black", hjust = 1,vjust = 1)) + 
  labs(x = NULL, y = 'Plant height of\nfirst measurements (cm, sqrt)' ,title = NULL) + 
  scale_y_continuous(position = "left",labels = scales::label_comma(accuracy =0.1))+
  geom_signif(comparisons = list(c(1,2)),test="t.test", annotations='ns.',tip_length = 0.02,size = 0.5,
              textsize = 4,y_position = 7.5) -> p1a; p1a

####
#mod12 = relmatLmer(Hmax ~ Hmax_time1 + (1|Species), data = trait_data, relmat = list(Species=phyloMat))
mod12 <- lmer((Hmax) ~ Hmax_time1 + (1|Species), data=trait_data)
anova(mod12)
summary(mod12)
MuMIn::r.squaredGLMM(mod12)

summary(lm(Hmax ~ Hmax_time1, data = trait_data))
ggplot(trait_data, mapping = aes(x = Hmax_time1 , y = Hmax)) + 
  geom_smooth(method = lm, se = T, color = "black" , fill = "#EBEBEB", alpha = 0.6) + 
  geom_point(trait_data, mapping = aes(x = Hmax_time1 , y = Hmax, shape = Origin, fill = Origin),size=2.2, color = "black")+
  scale_shape_manual(values = c(21,21))+
  scale_color_manual(values=c('#60A7A6','#FEA6A6'))+
  scale_fill_manual(values=c('#60A7A6','#FEA6A6'))+
  labs(x = 'Plant height of first measurements (cm, sqrt)',
       y = 'Maximum height (cm, sqrt)') + 
  mytheme + 
  scale_x_continuous(labels = scales::label_comma(accuracy =0.1)) +
  scale_y_continuous(labels = scales::label_comma(accuracy =0.1)) + 
  annotate('text', label = expression(italic(R)^2~"="~0.427*","~italic(p)~"<"~0.001), x = 4.3, y = 2.5, size = 4, parse = TRUE) -> p1b; p1b

library(patchwork)
a = (p1a+p1b) + plot_layout(widths = c(0.3,0.7)) + plot_annotation(tag_levels = 'a')

#### Plant height of second measurements (cm, sqrt)
ggplot(trait_data,aes(x=Origin,y=Hmax_time2))+
  geom_violin(data=trait_data, aes(y=Hmax_time2,x=Origin,fill=Origin,color=Origin),trim=T,
              scale = "width",position = position_dodge(0.5),width=0.6,alpha = 0.8, size = 1.5)+
  geom_boxplot(data=trait_data, aes(y=Hmax_time2,x=Origin),fill="white",color="black",
               position = position_dodge(0.5),width=0.15,outlier.size=0.8,outlier.shape = 1)+
  stat_summary(data=trait_data, aes(y=Hmax_time2,x=Origin),fill="white",color="black",
               fun="mean",position = position_dodge(0.5),geom="point",shape=21, size=1.5)+
  scale_color_manual(values = c('#60A7A6','#FEA6A6'))+
  scale_fill_manual(values = c('#60A7A6','#FEA6A6'))+
  mytheme+
  theme(axis.text.x = element_text(angle = 25,size=11,colour = "black", hjust = 1,vjust = 1)) + 
  labs(x = NULL, y = 'Plant height of\nsecond measurements (cm, sqrt)' ,title = NULL) + 
  scale_y_continuous(position = "left",labels = scales::label_comma(accuracy =0.1))+
  geom_signif(comparisons = list(c(1,2)),test="t.test", annotations='ns.',tip_length = 0.02,size = 0.5,
              textsize = 4,y_position = 10.5) -> p2a; p2a


summary(lm(Hmax ~ Hmax_time2, data = trait_data))
ggplot(trait_data, mapping = aes(x = Hmax_time2 , y = Hmax)) + 
  geom_smooth(method = lm, se = T, color = "black" , fill = "#EBEBEB", alpha = 0.6) + 
  geom_point(trait_data, mapping = aes(x = Hmax_time2 , y = Hmax, shape = Origin, fill = Origin),size=2.2, color = "black")+
  scale_shape_manual(values = c(21,21))+
  scale_color_manual(values=c('#60A7A6','#FEA6A6'))+
  scale_fill_manual(values=c('#60A7A6','#FEA6A6'))+
  labs(x = 'Plant height of second measurements (cm, sqrt)',
       y = 'Maximum height (cm, sqrt)') + 
  mytheme + 
  scale_x_continuous(labels = scales::label_comma(accuracy =0.1)) +
  scale_y_continuous(labels = scales::label_comma(accuracy =0.1)) + 
  #geom_abline(intercept=0.3272 ,slope=1.0290 ,size=.8, linetype = 1)+
  annotate('text', label = expression(italic(R)^2~"="~0.865*","~italic(p)~"<"~0.001), 
           x = 6.3, y = 2.5, size = 4, parse = TRUE) -> p2b; p2b


library(patchwork)
b = (p2a+p2b) + plot_layout(widths = c(0.3,0.7)) + plot_annotation(tag_levels = 'a')
(a/b) + plot_annotation(tag_levels = "a")

