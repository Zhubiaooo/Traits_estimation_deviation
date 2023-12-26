################################################################################
##################################   Fig 1a   ##################################
################################################################################
### This part mainly visualizes the phylogenetic relationship of species used in pot and field experiments.

### Loading Packages
library(openxlsx)
library(V.PhyloMaker)
library(plantlist)
library(ape)
library(tidyverse)
library(treeio)
library(ggtree)
library(rphylopic)
library(patchwork)
setwd(getwd())

### The mega-tree “GBOTB.extended.tre” with the V. PhyloMaker R package 
tree_data <- read.xlsx("Data/Total_species_list.xlsx", sheet = "Species_tree", colNames = T, rowNames = F)
#knitr::kable(head(tree_data))
# Build a phylogeny using phylo.maker
species <- gsub("_", " ", as.character(tree_data$Species))
sp.list <- subset(TPL(gsub("_", " ", species)), select = c("YOUR_SEARCH", "POSSIBLE_GENUS", "FAMILY"))
colnames(sp.list) <- c("species", "genus", "family")
#knitr::kable(head(sp.list))sp.list

phyloPlants <- phylo.maker(sp.list = sp.list, tree = GBOTB.extended,
                           nodes = nodes.info.1,
                           output.sp.list = TRUE, output.tree = FALSE,
                           scenarios = "S1")
phylogenyAux <- phyloPlants$scenario.1
plot(phylogenyAux)

### Loading tree file
tree <- read.tree("Data/All_species.newick") |> multi2di()
plot(tree)
mf.clrs <- c('#60A7A6','#FEA6A6')
### Read additional data, including families, genera, species, etc.
dat <- read.xlsx("Data/Common_species_list.xlsx", sheet = "Common_sp_list", colNames = T, rowNames = F)
head(dat)

### Delete the dataset where there is no species (delete a branch length) to get the common species information
to_drop <- tree$tip.label[!tree$tip.label %in% dat$Species]
tree <- drop.tip(as.phylo(tree), to_drop) 

### Add groups to tree
# setNames(to.new.name, from.old.name)
dat$Species <- str_replace_all(dat$Species, "_", " ")
mf.species.replacement <- setNames(dat$Origin, nm = dat$Species)
Origin <- split(tree$tip.label, recode(str_replace_all(tree$tip.label, "_", " "), !!!mf.species.replacement))
tree <- groupOTU(tree, Origin, group_name = "Origin") 

tree_df_raw <- fortify(tree)
tree_df_raw$label <- str_replace_all(tree_df_raw$label, "_", " ")
tree_df <- full_join(tree_df_raw, dat[,c(-2)], by = c("label" = "Species"))
tree_df$Origin <- factor(tree_df$Origin, levels = c("Native","Exotic"))

###
ordater.species.replacement <- setNames(dat$Family, nm = dat$Species)
tree_df2 <- tree_df |>
  mutate(label = str_replace_all(label, "_", " "),order = case_when(
    isTip ~ recode(label, !!!ordater.species.replacement), TRUE ~ "")) # 只替换 tip 的 label

## Maximum time span
max.geo.time <- max(tree_df2$x, na.rm = TRUE)
max.geo.time

p.tree2 <- tree_df2 |> 
  ggtree::ggtree(size = 0.4, color = "black") +
  geom_tiplab(aes(label = sub("_", " ", label)),size = 2, offset=3, fontface = "italic") +
  geom_tippoint(aes(color = Origin),size = 1) +
  scale_color_manual(name = NULL, values = mf.clrs) +
  scale_x_continuous(breaks = seq(max.geo.time - 140, max.geo.time, 20),
                     labels = rev(seq(0, 140, 20))) + # Add differentiation Timeline
  coord_cartesian(xlim = c(max.geo.time - 140, max.geo.time),ylim = c(0, 65),expand = FALSE,clip = "off") +
  theme(plot.background = element_blank(),
        #panel.grid.major.x = element_line(color = "grey", linewidth = 0.2, linetype = 2),
        plot.margin = margin(0.5,4,0.5,4, unit = "cm"),
        strip.background = element_blank(),
        strip.text = element_text(size = 10),
        axis.text.x = element_text(size = 8),
        axis.line.x = element_line(),
        axis.ticks.x.bottom = element_line(),
        legend.position = c(0.2,0.1), 
        #legend.background = element_rect(fill=rgb(1,1,1,alpha=0.001), color = "grey"),
        legend.title = element_blank()) ;p.tree2
### rotate tree
#p.tree2+geom_nodelab(aes(label=node))
ggtree::rotate(p.tree2,node = 66) -> p.tree3 ;p.tree3

###
right.df <- tree_df2[c(1:65),] |> 
  summarise(y1 = min(y), y2 = max(y), y = mean(y),.by = order) |> 
  na.omit()
# right.df |> View()

### Set the color corresponding to each Family
order.clrs = c("Poaceae" = "#1F5346", "Malvaceae" = "#277567", "Fabaceae" = "#549D94", "Polygonaceae" = "#82C4B8", 
               "Nyctaginaceae" = "#B9E1D8","Amaranthaceae" = "#EAD09D", "Plantaginaceae" = "#E4CC90", 
               "Lamiaceae" = "#BF8D46", "Solanaceae" = "#905A1C", "Apiaceae" = "#55330E", "Asteraceae" = "#44310D")

p.right <- tree_df2 |> 
  filter(isTip) |> 
  ggplot(aes(x = -0.2, y, color = order, fill = order)) +
  ggrepel::geom_text_repel(data = right.df,aes(x = -0.2, label = order),size = 2.7,
                           max.overlaps = 100,direction = "y",xlim = c(1, 3),segment.size = 0.2,hjust = 0, color = "black") +
  geom_point(data = right.df,aes(x = -0.2),size = 1) +
  geom_segment(data = right.df,aes(x = -0.2, xend = -0.2, y = y1, yend = y2),lwd = 6) +
  scale_color_manual(name = NULL, values = order.clrs) +
  scale_fill_manual(name = NULL, values = order.clrs) +
  #facet_grid(. ~ "Order") +
  scale_y_reverse() +
  coord_cartesian(xlim = c(-0.5, 5),ylim = c(65.5, 1),expand = FALSE,clip = "off") +
  theme_tree() +
  theme(plot.background = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 10),
        legend.position = "none"); p.right

### manually set the x-axis position of PhyloPic
order = c("Asteraceae","Apiaceae","Solanaceae","Lamiaceae","Plantaginaceae","Amaranthaceae",
          "Nyctaginaceae","Polygonaceae","Fabaceae","Malvaceae","Poaceae")
right.df2 = data.frame(order = order,
                       y = c(seq(0, 50, by = 5)),
                       pic.x.pos = c(rep(c(1,1.8), 5), 1))

### Add plant silhouette
p.pic <- right.df2 |> 
  ggplot(aes(x = 0, y, color = order, fill = order)) +
  lapply( # add PhyloPic in batch
    1:nrow(right.df2),
    function(i) {
      order.id <- unique(right.df2$order)
      temp.df <- subset(as.data.frame(right.df2), order == order.id[i])
      rphylopic::add_phylopic(x = temp.df$pic.x.pos, y = temp.df$y, ysize = 5,
                              name = temp.df$order,color = order.clrs[temp.df$order])
    }) + # add unloaded pictures manually
  # ggrepel::geom_text_repel(
  #   data = right.df2,
  #   aes(x = pic.x.pos, label = order),
  #   size = 2.5,
  #   min.segment.length = unit(0, "pt")
  # ) +
  shadowtext::geom_shadowtext(
    data = right.df2,
    aes(x = pic.x.pos, label = order),size = 4,bg.color = "white") +
  scale_color_manual(name = NULL, values = order.clrs) +
  scale_fill_manual(name = NULL, values = order.clrs) +
  coord_cartesian(xlim = c(0.5, 4),ylim = c(-2, 52),expand = FALSE,clip = "off") +
  theme_tree() +
  theme(plot.background = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 10),
        legend.position = "none")
p.pic

### Stitching pictures (Fig 1a)
tree_all = p.tree3+p.right+p.pic+plot_layout(ncol = 3,widths = c(1,2,2)); tree_all

##################################################################################
#############################   Fig 1b,1c,1d   ###################################
##################################################################################
### This section of analysis mainly focuses on comparing and analyzing 
### the common species of traits in pot experiments and field experiments

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

### Reading into the plant phylogenetic tree
phylogenyAux <- read.tree("Data/All_species.newick")
### Species phylogenetic correlation matrix
phyloMat <- vcv.phylo( phylogenyAux )
phyloMat <- phyloMat / max( phyloMat )
dim(phyloMat)

### Analyzing species lists (SLA = 48, Hmax&AGB = 64)
trait_data <- read.xlsx("Data/Common_species_list.xlsx", sheet = "Common_sp_list", colNames = TRUE, rowNames = FALSE)
Common_sp_list_SLA <- c(na.omit(trait_data$Species_SLA))
Common_sp_list_AGB <- c(na.omit(trait_data$Species_AGB))
par(mar = c(5,9,5,9))

### Specific leaf area of common species in pot experiment
trait_data <- read.xlsx("Data/Pot_traits_database.xlsx", sheet = "Pot_data", colNames = TRUE, rowNames = FALSE)
trait_data <- trait_data[trait_data$Species %in% unique(Common_sp_list_SLA), ]
trait_data <- trait_data %>% drop_na(SLA)
length(unique(trait_data$Species))

#shapiro.test(trait_data$SLA)
#shapiro.test(sqrt(trait_data$SLA))
shapiro.test(log(trait_data$SLA, 10))
trait_data$Origin <- factor(trait_data$Origin, levels = c("Native", "Exotic"))
trait_data$Species <- as.factor(trait_data$Species)
mod <- relmatLmer(log(SLA, 10) ~ Origin + (1|Species), data = trait_data, relmat = list(Species=phyloMat))
summary(mod)
Anova(mod, type="II", test.statistic ="F")
shapiro.test(residuals(mod))
#ggpubr::ggqqplot(residuals(mod))
#hist(residuals(mod))
tuk_quasipoisson <- glht(mod, alternative = 'two.sided', linfct = mcp(Origin = 'Tukey'))
summary(tuk_quasipoisson) 
tuk_quasipoisson.cld <- cld(tuk_quasipoisson, level = 0.05, decreasing = TRUE)
plot(tuk_quasipoisson.cld, col = c('#60A7A6','#FEA6A6'))

###Percentage increase
summarySE(trait_data, measurevar = c("SLA"), groupvars = c("Origin"))
(184.8034-168.9243)/168.9243

###Phylogenetic signal test of plant traits
summary_data <- summarySE(trait_data, measurevar="SLA", groupvars=c("Species"))
test = make.treedata(tree = phylogenyAux,  data = summary_data[,c(1,3)], name_column = "Species")
myTree.withBrLe <- compute.brlen(test$phy)
myProx <- vcv.phylo(myTree.withBrLe)
set.seed(4321)
Physignal <- abouheif.moran(test$dat, W = myProx, method = "oriAbouheif", nrepet = 999)
Physignal
plot(Physignal)

###Maximum height of common species in pot experiment
trait_data <- read.xlsx("Data/Pot_traits_database.xlsx", sheet = "Pot_data", colNames = TRUE, rowNames = FALSE)
trait_data <- trait_data[trait_data$Species %in% unique(Common_sp_list_AGB), ]
trait_data <- trait_data %>% drop_na(Hmax)
length(unique(trait_data$Species))

#shapiro.test(trait_data$Hmax)
shapiro.test(sqrt(trait_data$Hmax))
#shapiro.test(log(trait_data$Hmax, 10))
trait_data$Origin = factor(trait_data$Origin, levels = c("Native", "Exotic"))

mod <- relmatLmer(sqrt(Hmax) ~ Origin + (1|Species), data = trait_data, relmat = list(Species=phyloMat))
summary(mod)
Anova(mod, type="II", test.statistic ="F")
shapiro.test(residuals(mod))
#ggpubr::ggqqplot(residuals(mod))
hist(residuals(mod))

tuk_quasipoisson <- glht(mod, alternative = 'two.sided', linfct = mcp(Origin = 'Tukey'))
summary(tuk_quasipoisson) 
tuk_quasipoisson.cld <- cld(tuk_quasipoisson, level = 0.05, decreasing = TRUE)
plot(tuk_quasipoisson.cld, col = c('#60A7A6','#FEA6A6'))

###Percentage increase
summarySE(trait_data, measurevar = c("Hmax"), groupvars = c("Origin"))
(26.19310-25.91875)/26.19310

###Phylogenetic signal test of plant traits
summary_data <- summarySE(trait_data, measurevar = "Hmax", groupvars=c("Species"))
test <- make.treedata(tree = phylogenyAux, data = summary_data[,c(1,3)], name_column = "Species")
myTree.withBrLe <- compute.brlen(test$phy)
myProx <- vcv.phylo(myTree.withBrLe)
set.seed(4321)
Physignal <- abouheif.moran(test$dat, W = myProx, method = "oriAbouheif", nrepet = 999)
Physignal
plot(Physignal)

####Aboveground biomass of common species in pot experiment
trait_data <- read.xlsx("Data/Pot_traits_database.xlsx", sheet = "Pot_data", colNames = TRUE, rowNames = FALSE)
trait_data <- trait_data[trait_data$Species %in% unique(Common_sp_list_AGB), ]
trait_data <- trait_data %>% drop_na(AGB)
length(unique(trait_data$Species))

#shapiro.test(trait_data$AGB)
#shapiro.test(sqrt(trait_data$AGB))
shapiro.test(log(trait_data$AGB, 10))
trait_data$Origin <- factor(trait_data$Origin, levels = c("Native", "Exotic"))

mod <- relmatLmer(log10(AGB) ~ Origin + (1|Species), data = trait_data, relmat = list(Species=phyloMat))
summary(mod)
Anova(mod, type="II", test.statistic ="F")
shapiro.test(residuals(mod))
#ggpubr::ggqqplot(residuals(mod))
hist(residuals(mod))

tuk_quasipoisson <- glht(mod, alternative = 'two.sided', linfct = mcp(Origin = 'Tukey'))
summary(tuk_quasipoisson) 
tuk_quasipoisson.cld <- cld(tuk_quasipoisson, level = 0.05, decreasing = TRUE)
plot(tuk_quasipoisson.cld, col = c('#60A7A6','#FEA6A6'))

###Percentage increase
summarySE(trait_data, measurevar = c("AGB"), groupvars = c("Origin"))
(1.644797-1.537843)/1.537843

###Phylogenetic signal test of plant traits
summary_data <- summarySE(trait_data, measurevar="AGB", groupvars=c("Species"))
test <- make.treedata(tree = phylogenyAux,  data = summary_data[,c(1,3)], name_column = "Species")
myTree.withBrLe <- compute.brlen(test$phy)
myProx <- vcv.phylo(myTree.withBrLe)
set.seed(4321)
Physignal <- abouheif.moran(test$dat, W = myProx, method = "oriAbouheif", nrepet = 999)
Physignal
plot(Physignal)

###Specific leaf area of common species in field experiment
trait_data <- read.xlsx("Data/Field_traits_database(SLA).xlsx", sheet = "Field_SLA", colNames = TRUE, rowNames = FALSE)
trait_data <- trait_data[trait_data$Species %in% unique(Common_sp_list_SLA), ]
trait_data <- trait_data %>% drop_na(SLA)
length(unique(trait_data$Species))

#shapiro.test(trait_data$SLA)
#shapiro.test(sqrt(trait_data$SLA))
shapiro.test(log(trait_data$SLA, 10))
trait_data$Origin <- factor(trait_data$Origin, levels = c("Native", "Exotic"))
mod <- relmatLmer(log(SLA, 10) ~ Origin +(1|Species), data = trait_data, relmat = list(Species=phyloMat))
summary(mod)
Anova(mod, type="II", test.statistic ="F")
shapiro.test(residuals(mod))
#ggpubr::ggqqplot(residuals(mod))
#hist(residuals(mod))

tuk_quasipoisson <- glht(mod, alternative = 'two.sided', linfct = mcp(Origin = 'Tukey'))
summary(tuk_quasipoisson) 
tuk_quasipoisson.cld = cld(tuk_quasipoisson, level = 0.05, decreasing = TRUE)
plot(tuk_quasipoisson.cld, col = c('#60A7A6','#FEA6A6'))

###Percentage increase
summarySE(trait_data, measurevar = c("SLA"), groupvars = c("Origin"))
(178.6640-171.0985)/171.0985	

###Phylogenetic signal test of plant traits
summary_data <- summarySE(trait_data, measurevar="SLA", groupvars=c("Species"))
test = make.treedata(tree = phylogenyAux,  data = summary_data[,c(1,3)], name_column = "Species")
myTree.withBrLe <- compute.brlen(test$phy)
myProx <- vcv.phylo(myTree.withBrLe)
set.seed(4321)
Physignal = abouheif.moran(test$dat, W = myProx, method = "oriAbouheif", nrepet = 999)
Physignal
plot(Physignal)

###Maximum height of common species in field garden experiment
trait_data <- read.xlsx("Data/Field_traits_database(Hmax&AGB).xlsx", sheet = "Field_Hmax&AGB", colNames = TRUE, rowNames = FALSE)
trait_data <- trait_data[trait_data$Species %in% unique(Common_sp_list_AGB), ]
trait_data <- trait_data %>% drop_na(Field_AGB)
length(unique(trait_data$Species))

#shapiro.test(trait_data$Field_Hmax)
shapiro.test(sqrt(trait_data$Field_Hmax))
#shapiro.test(log(trait_data$Field_Hmax, 10))
trait_data$Origin <- factor(trait_data$Origin, levels = c("Native", "Exotic"))
trait_data$Block <- as.factor(trait_data$Block)
trait_data$Plot_num <- as.factor(trait_data$Plot_num)

mod <- relmatLmer(sqrt(Field_Hmax) ~ Origin + (1|Block/Plot_num) + (1|Species) , data = trait_data, relmat = list(Species=phyloMat))
summary(mod)
Anova(mod, type="II", test.statistic ="F")
shapiro.test(residuals(mod))
#ggpubr::ggqqplot(residuals(mod))
#hist(residuals(mod))

tuk_quasipoisson <- glht(mod, alternative = 'two.sided', linfct = mcp(Origin = 'Tukey'))
summary(tuk_quasipoisson) 
tuk_quasipoisson.cld <- cld(tuk_quasipoisson, level = 0.05, decreasing = TRUE)
plot(tuk_quasipoisson.cld, col = c('#60A7A6','#FEA6A6'))

###Percentage increase
summarySE(trait_data, measurevar = c("Field_Hmax"), groupvars = c("Origin"))
(123.78475-97.18785)/97.18785	

###Phylogenetic signal test of plant traits
summary_data <- summarySE(trait_data, measurevar="Field_Hmax", groupvars=c("Species"))
test <- make.treedata(tree = phylogenyAux,  data = summary_data[,c(1,3)], name_column = "Species")
myTree.withBrLe <- compute.brlen(test$phy)
myProx <- vcv.phylo(myTree.withBrLe)
set.seed(4321)
Physignal <- abouheif.moran(test$dat, W = myProx, method = "oriAbouheif", nrepet = 999)
Physignal
plot(Physignal)

###Aboveground biomass of common species in field experiment
trait_data <- read.xlsx("Data/Field_traits_database(Hmax&AGB).xlsx", sheet = "Field_Hmax&AGB", colNames = TRUE, rowNames = FALSE)
trait_data <- trait_data[trait_data$Species %in% unique(Common_sp_list_AGB), ]
trait_data <- trait_data %>% drop_na(Field_AGB)
length(unique(trait_data$Species))

#shapiro.test(trait_data$Field_AGB)
#shapiro.test(sqrt(trait_data$Field_AGB))
shapiro.test(log(trait_data$Field_AGB, 10))
trait_data$Origin <- factor(trait_data$Origin, levels = c("Native", "Exotic"))
trait_data$Block <- as.factor(trait_data$Block)
trait_data$Plot_num <- as.factor(trait_data$Plot_num)

mod <- relmatLmer(log(Field_AGB, 10) ~ Origin + (1|Block/Plot_num) + (1|Species)  , data = trait_data, relmat = list(Species=phyloMat))
summary(mod)
Anova(mod, type="II", test.statistic ="F")
shapiro.test(residuals(mod))
#ggpubr::ggqqplot(residuals(mod))
#hist(residuals(mod))

tuk_quasipoisson <- glht(mod, alternative = 'two.sided', linfct = mcp(Origin = 'Tukey'))
summary(tuk_quasipoisson) 
tuk_quasipoisson.cld <- cld(tuk_quasipoisson, level = 0.05, decreasing = TRUE)
plot(tuk_quasipoisson.cld, col = c('#60A7A6','#FEA6A6'))

### Percentage increase
summarySE(trait_data, measurevar = c("Field_AGB"), groupvars = c("Origin"))
(151.77568-50.41337)/50.41337

### Phylogenetic signal test of plant traits
summary_data <- summarySE(trait_data, measurevar = "Field_AGB", groupvars=c("Species"))
test <- make.treedata(tree = phylogenyAux, data = summary_data[,c(1,3)], name_column = "Species")
myTree.withBrLe <- compute.brlen(test$phy)
myProx <- vcv.phylo(myTree.withBrLe)
set.seed(4321)
Physignal <- abouheif.moran(test$dat, W = myProx, method = "oriAbouheif", nrepet = 999)
Physignal
plot(Physignal)

### Comparison and visualization of common species traits in two experiments
### Set the picture theme
mytheme = theme_classic()+
  theme(panel.grid = element_blank(), panel.background = element_blank(), plot.title = element_text(size = 14), 
        axis.ticks.y = element_line(), axis.line.y = element_line(),axis.line.x = element_line(color = 'black'),
        axis.text.y = element_text(size = 11,angle = 0, colour = "black"),
        axis.text.x = element_text(angle = 25,size=11,colour = "black", hjust = 1,vjust = 1),
        axis.title.y=element_text(size = 13),
        axis.title.x=element_text(size = 13),
        plot.margin = unit(c(1, 3, 1, 3), "lines"),
        legend.position = "none")+
  theme(plot.title = element_text(color = "black", size = 13, hjust = 0.5))

### Pot experiment
trait_data = read.xlsx("Data/Pot_traits_database.xlsx", sheet = "Pot_data", colNames = TRUE, rowNames = FALSE)
trait_data = trait_data[trait_data$Species %in% unique(Common_sp_list_SLA), ]
trait_data = trait_data %>% drop_na(SLA)
trait_data$Origin=factor(trait_data$Origin, levels = c("Native","Exotic"))
length(unique(trait_data$Species))
trait_data$SLA = log10(trait_data$SLA)

p1 = ggplot(trait_data,aes(x=Origin,y=SLA))+
  geom_violin(data=trait_data, aes(y=SLA,x=Origin,fill=Origin,color=Origin),trim=T,scale = "width",
              position = position_dodge(0.5),width=0.6,alpha = 0.8,linewidth = 1.5)+
  geom_boxplot(data=trait_data, aes(y=SLA,x=Origin),fill="white",color="black",
               position = position_dodge(0.5),width=0.15,outlier.size=0.8,outlier.shape = 1)+
  stat_summary(data=trait_data, aes(y=SLA,x=Origin),fill="white",color="black",
               fun="mean",position = position_dodge(0.5),geom="point",shape=21, size=1.5)+
  scale_color_manual(values = c('#60A7A6','#FEA6A6'))+
  scale_fill_manual(values = c('#60A7A6','#FEA6A6'))+
  mytheme+
  labs(x = NULL, y = expression('Specific leaf area (cm'^ 2*'/g, log10)'),  title = NULL) + 
  scale_y_continuous(position = "left",labels = scales::label_comma(accuracy =0.1), limits = c(1,3.7))+
  geom_signif(comparisons = list(c(1,2)),test="t.test", annotations='ns.',tip_length = 0.02,size = 0.5,
              textsize = 4,y_position = 3.4); p1


###
trait_data = read.xlsx("Data/Pot_traits_database.xlsx", sheet = "Pot_data", colNames = TRUE, rowNames = FALSE)
trait_data = trait_data[trait_data$Species %in% unique(Common_sp_list_AGB), ]
trait_data = trait_data %>% drop_na(Hmax)
length(unique(trait_data$Species))
trait_data$Origin=factor(trait_data$Origin, levels = c("Native","Exotic"))
trait_data$Hmax = sqrt(trait_data$Hmax)

p2 = ggplot(trait_data,aes(x=Origin,y=Hmax))+
  geom_violin(data=trait_data, aes(y=Hmax,x=Origin,fill=Origin,color=Origin),trim=T,
              scale = "width",position = position_dodge(0.5),width=0.6,alpha = 0.8, size = 1.5)+
  geom_boxplot(data=trait_data, aes(y=Hmax,x=Origin),fill="white",color="black",
               position = position_dodge(0.5),width=0.15,outlier.size=0.8,outlier.shape = 1)+
  stat_summary(data=trait_data, aes(y=Hmax,x=Origin),fill="white",color="black",
               fun="mean",position = position_dodge(0.5),geom="point",shape=21, size=1.5)+
  scale_color_manual(values = c('#60A7A6','#FEA6A6'))+
  scale_fill_manual(values = c('#60A7A6','#FEA6A6'))+
  mytheme+
  labs(x = NULL, y = 'Maximum height (cm, sqrt)' ,title = NULL) + 
  scale_y_continuous(position = "left",labels = scales::label_comma(accuracy =0.1))+
  geom_signif(comparisons = list(c(1,2)),test="t.test", annotations='ns.',tip_length = 0.02,size = 0.5,
              textsize = 4,y_position = 10) ; p2


###
trait_data = read.xlsx("Data/Pot_traits_database.xlsx", sheet = "Pot_data", colNames = TRUE, rowNames = FALSE)
trait_data = trait_data[trait_data$Species %in% unique(Common_sp_list_AGB), ]
trait_data = trait_data %>% drop_na(AGB)
length(unique(trait_data$Species))
##AA %in% BB
trait_data$Origin=factor(trait_data$Origin, levels = c("Native","Exotic"))
trait_data$AGB = log10(trait_data$AGB)

p3 = ggplot(trait_data,aes(x=Origin,y=AGB))+
  geom_violin(data=trait_data, aes(y=AGB,x=Origin,fill=Origin,color=Origin),trim=T,scale = "width",
              position = position_dodge(0.5),width=0.6,alpha = 0.8,size = 1.5)+
  geom_boxplot(data=trait_data, aes(y=AGB,x=Origin),fill="white",color="black",
               position = position_dodge(0.5),width=0.15,outlier.size=0.8,outlier.shape = 1)+
  stat_summary(data=trait_data, aes(y=AGB,x=Origin),fill="white",color="black",
               fun="mean",position = position_dodge(0.5),geom="point",shape=21, size=1.5)+
  scale_color_manual(values = c('#60A7A6','#FEA6A6'))+
  scale_fill_manual(values = c('#60A7A6','#FEA6A6'))+
  mytheme+
  labs(x = NULL, y = 'Aboveround biomass (g, log10)' ,title = NULL) + 
  scale_y_continuous(position = "left",labels = scales::label_comma(accuracy =0.1), limits = c(-1.7,1.5))+
  geom_signif(comparisons = list(c(1,2)),test="t.test", annotations='ns.',tip_length = 0.02,size = 0.5,
              textsize = 4,y_position = 1.2); p3



### Traits measured in Field experiment
trait_data = read.xlsx("Data/Field_traits_database(SLA).xlsx", sheet = "Field_SLA", colNames = TRUE, rowNames = FALSE)
trait_data = trait_data[trait_data$Species %in% unique(Common_sp_list_SLA), ] %>% drop_na(SLA)
length(unique(trait_data$Species))
trait_data$SLA = log10(trait_data$SLA)
trait_data$Origin = factor(trait_data$Origin, levels = c("Native", "Exotic"))

p4 = ggplot(trait_data,aes(x=Origin,y=SLA))+
  geom_violin(data=trait_data, aes(y=SLA,x=Origin,fill=Origin,color=Origin),trim=T,scale = "width",
              position = position_dodge(0.5),width=0.6,alpha = 0.8,size = 1.5)+
  geom_boxplot(data=trait_data, aes(y=SLA,x=Origin),fill="white",color="black",
               position = position_dodge(0.5),width=0.15,outlier.size=0.8,outlier.shape = 1)+
  stat_summary(data=trait_data, aes(y=SLA,x=Origin),fill="white",color="black",
               fun="mean",position = position_dodge(0.5),geom="point",shape=21, size=1.5)+
  scale_color_manual(values = c('#60A7A6','#FEA6A6'))+
  scale_fill_manual(values = c('#60A7A6','#FEA6A6'))+
  labs(x = NULL, y = expression('Specific leaf area (cm'^ 2*'/g, log10)'),  title = NULL) + 
  mytheme+
  scale_y_continuous(position = "right",labels = scales::label_comma(accuracy =0.1))+
  geom_signif(comparisons = list(c(1,2)),test="t.test", annotations='ns.',tip_length = 0.02,size = 0.5,
              textsize = 4,y_position = 2.7); p4

###
trait_data = read.xlsx("Data/Field_traits_database(Hmax&AGB).xlsx", sheet = "Field_Hmax&AGB", colNames = TRUE, rowNames = FALSE)
trait_data = trait_data[trait_data$Species %in% unique(Common_sp_list_AGB), ] %>% drop_na(Field_Hmax)
length(unique(trait_data$Species))
trait_data$Origin=factor(trait_data$Origin, levels = c("Native","Exotic"))
trait_data$Field_Hmax = sqrt(trait_data$Field_Hmax)

p5 = ggplot(trait_data,aes(x=Origin,y=Field_Hmax))+
  geom_violin(data=trait_data, aes(y=Field_Hmax,x=Origin,fill=Origin,color=Origin),trim=T,scale = "width",
              position = position_dodge(0.5),width=0.6,alpha = 0.8,size = 1.5)+
  geom_boxplot(data=trait_data, aes(y=Field_Hmax,x=Origin),fill="white",color="black",
               position = position_dodge(0.5),width=0.15,outlier.size=0.8,outlier.shape = 1)+
  stat_summary(data=trait_data, aes(y=Field_Hmax,x=Origin),fill="white",color="black",
               fun="mean",position = position_dodge(0.5),geom="point",shape=21, size=1.5)+
  scale_color_manual(values = c('#60A7A6','#FEA6A6'))+
  scale_fill_manual(values = c('#60A7A6','#FEA6A6'))+
  mytheme+
  labs(x = NULL, y = 'Maximum height (cm, sqrt)' ,title = NULL) + 
  scale_y_continuous(position = "right",labels = scales::label_comma(accuracy =1))+
  geom_signif(comparisons = list(c(1,2)),test="t.test", annotations='ns.',tip_length = 0.02,size = 0.5,
              textsize = 4,y_position = 18); p5


###
trait_data = read.xlsx("Data/Field_traits_database(Hmax&AGB).xlsx", sheet = "Field_Hmax&AGB", colNames = TRUE, rowNames = FALSE)
trait_data = trait_data[trait_data$Species %in% unique(Common_sp_list_AGB), ] %>% drop_na(Field_AGB)
length(unique(trait_data$Species))
trait_data$Origin=factor(trait_data$Origin, levels = c("Native","Exotic"))
trait_data$Field_AGB = log10(trait_data$Field_AGB)

p6 = ggplot(trait_data,aes(x=Origin,y=Field_AGB))+
  geom_violin(data=trait_data, aes(y=Field_AGB,x=Origin,fill=Origin,color=Origin),trim=T,scale = "width",
              position = position_dodge(0.5),width=0.6,alpha = 0.8,size = 1.5)+
  geom_boxplot(data=trait_data, aes(y=Field_AGB,x=Origin),fill="white",color="black",size = 0.5,
               position = position_dodge(0.5),width=0.15,outlier.size=0.8,outlier.shape = 1)+
  stat_summary(data=trait_data, aes(y=Field_AGB,x=Origin),fill="white",color="black",
               fun="mean",position = position_dodge(0.5),geom="point",shape=21, size=1.5)+
  scale_color_manual(values = c('#60A7A6','#FEA6A6'))+
  scale_fill_manual(values = c('#60A7A6','#FEA6A6'))+
  mytheme+
  labs(x = NULL, y = 'Aboveround biomass (g, log10)' ,title = NULL) + 
  scale_y_continuous(position = "right",labels = scales::label_comma(accuracy =0.1))+
  geom_signif(comparisons = list(c(1,2)),test="t.test", annotations='*',tip_length = 0.02,size = 0.5,
              textsize = 4,y_position = 3.5); p6


(p1+p4+p2+p5+p3+p6) + plot_layout(ncol = 2, nrow = 3)

### Notice that,
### For more picture details, we have further adjusted it in Adobe illustrator.