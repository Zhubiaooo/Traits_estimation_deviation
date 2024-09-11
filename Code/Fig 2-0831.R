### Performance of the same set of species in the pot and field experiments, 
### Correlation between pot and field traits

################################################################################
##################################   Fig 2   ###################################
################################################################################

### Loading Packages
library(openxlsx)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(missForest)
library(adephylo)
library(ape)
library(smatr)
library(deming)
library(funrar)
library(BestFitM)
library(ggtrendline)
library(AICcmodavg)
library(pscl)
library(vegan)
library(ggfortify)
library(ggExtra)
library(phytools)
library(lme4qtl)
library(car)


### Analyzing species lists
Common_sp_list = read.xlsx("Data/Common_species_list.xlsx", sheet = "Common_sp_list", colNames = TRUE, rowNames = FALSE)
Common_sp_list_SLA = c(na.omit(Common_sp_list$Species_SLA)); length(Common_sp_list_SLA)
Common_sp_list_AGB = c(na.omit(Common_sp_list$Species_AGB)); length(Common_sp_list_AGB)

### Custom themes and functions
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

r2 <- function(x) {x$`Sum Sq`[1]/ (x$`Sum Sq`[1] + x$`Sum Sq`[2])}
standr = function(x){(x-min(x))/(max(x)-min(x))} 

### Loading tree
phylogenyAux = read.tree("Data/iq_tree.treefile")
to_drop = phylogenyAux$tip.label[!phylogenyAux$tip.label %in% Common_sp_list$Species]
tree <- drop.tip(as.phylo(phylogenyAux), to_drop) 
phyloMat = vcv.phylo(tree)
phyloMat = phyloMat / max(phyloMat)
dim(phyloMat)
plot(tree)

### Loading pot experiment database (mean value of traits)
pot_trait = read.xlsx("Data/Pot_traits_mean0831.xlsx", sheet = "Pot_means", colNames = TRUE, rowNames = FALSE)
colnames(pot_trait) <- paste0(colnames(pot_trait), "_pot")
colnames(pot_trait)[1] <- "Species"
nrow(pot_trait)

#################################### Fig 2A ####################################
### Loading field experiment database (mean value of traits)
### Deming regressive analysis
### Specific leaf area
trait_data = read.xlsx("Data/Field_traits_mean0831_2.xlsx", sheet = "Field_means", colNames = TRUE, rowNames = FALSE)
colnames(trait_data) <- paste0(colnames(trait_data), "_field")
colnames(trait_data)[1] <- "Species"
nrow(trait_data)
head(trait_data)
trait_data = trait_data %>% left_join(pot_trait, by = "Species") %>% 
  left_join(Common_sp_list[,c(1,2,5)], by = "Species") %>%
  drop_na(SLA_pot) %>% drop_na(SLA_imp_field)
rownames(trait_data) = trait_data$Species
length(unique(trait_data$Species))

trait_data$SLA_imp_field_sqrt = sqrt(trait_data$SLA_imp_field)
trait_data$SLA_pot_sqrt = sqrt(trait_data$SLA_pot)

### Pot experiment SLA differences
dist_mat <- compute_dist_matrix(trait_data['SLA_pot_sqrt'], metric = 'euclidean',center = TRUE, scale = TRUE)
dist_mat <- standr(dist_mat)
diag(dist_mat) <- NA

### Calculate the mean of each row (excluding diagonals)
mean_vals <- apply(dist_mat, 1, mean, na.rm = TRUE)
### Calculate the standard deviation for each row (excluding diagonals)
sd_vals <- apply(dist_mat, 1, sd, na.rm = TRUE)
### Calculate the standard error for each line (removing diagonal lines)
se_vals <- sd_vals / sqrt(ncol(dist_mat) - 1)
SLA_pot_dis <- data.frame(SLA_pot_mean = mean_vals, SLA_pot_sd = sd_vals, SLA_pot_se = se_vals)

### Field SLA differences
dist_mat <- compute_dist_matrix(trait_data['SLA_imp_field_sqrt'], metric = 'euclidean',center = TRUE, scale = TRUE)#gower
dist_mat <- standr(dist_mat)
diag(dist_mat) <- NA
### Calculate the mean of each row (excluding diagonals)
mean_vals <- apply(dist_mat, 1, mean, na.rm = TRUE)
### Calculate the standard deviation for each row (excluding diagonals) 
sd_vals <- apply(dist_mat, 1, sd, na.rm = TRUE)
### Calculate the standard error for each line (removing diagonal lines) 
se_vals <- sd_vals / sqrt(ncol(dist_mat) - 1)
SLA_field_dis <- data.frame(SLA_Field_means = mean_vals, SLA_field_sd = sd_vals, SLA_field_se = se_vals)
SLA_total_dis = cbind(SLA_pot_dis, SLA_field_dis)
SLA_total_dis$Species = rownames(SLA_total_dis)
SLA_total_dis = SLA_total_dis %>% left_join(Common_sp_list[,c(1,2,5)], by = "Species")
#knitr::kable(head(SLA_total_dis))

###
Test_total_dis = SLA_total_dis
rownames(Test_total_dis) = Test_total_dis$Species
###
SLA_pot_di<-Test_total_dis$SLA_pot_mean
names(SLA_pot_di)<- rownames(Test_total_dis)
###
SLA_field_di<-Test_total_dis$SLA_Field_means
names(SLA_field_di)<- rownames(Test_total_dis)

### The effect of origin and types of experiments on DI
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
###
data_all_di = rbind(data1,data2) 
data_all_di = data_all_di %>% left_join(Common_sp_list[,c(1,2,5)], by = "Species")
###
shapiro.test(data_all_di$Di)
data_all_di$Origin = factor(data_all_di$Origin, levels = c("Native", "Exotic"))
data_all_di$Type = factor(data_all_di$Type, levels = c("pot", "field"))
data_all_di$Species = as.factor(data_all_di$Species)

mod = relmatLmer(Di ~ Origin*Type + (1|Species), data = data_all_di, relmat = list(Species=phyloMat))
Anova(mod, type="II", test.statistic ="Chisq")
Rmisc::summarySE(data_all_di, measurevar = "Di", groupvars = c("Type"))
####
library(emmeans)
emmeans(mod, specs = pairwise ~ Type|Origin, type = 'response', adjust = 'none')
emmeans(mod, specs = pairwise ~ Type, type = 'response', adjust = 'none')

fit <- deming(SLA_Field_means ~ SLA_pot_mean, ystd=SLA_field_sd, xstd=SLA_pot_sd, data=SLA_total_dis)
print(fit)
SLA_total_dis$Origin = factor(SLA_total_dis$Origin, levels = c("Native","Exotic"))

first_char <- substr(SLA_total_dis$Species, 1, 1)
sub_str <- gsub(".*_", "", SLA_total_dis$Species)
Latin_name <- paste(first_char, sub_str, sep = ". ")
SLA_total_dis$Latin_name = Latin_name

df2 = (SLA_total_dis %>% arrange(desc(SLA_Field_means)))[c(1:8),]
df3 = (SLA_total_dis %>% arrange(desc(SLA_pot_mean)))[c(1:8),]

dat1 = Rmisc::summarySE(SLA_total_dis, groupvars = c("Origin"), measurevar = c("SLA_pot_mean"))[,c(1,3,5)]
colnames(dat1)[3] = "pot_se"
dat2 = Rmisc::summarySE(SLA_total_dis, groupvars = c("Origin"), measurevar = c("SLA_Field_means"))[,c(3,5)]
colnames(dat2)[2] = "field_se"
dattt = cbind(dat1, dat2)

p1 = ggplot()+
  labs(x = " \n ",
       y= "Functional distinctiveness \n estimated in field experiment",
       title = "Specific leaf area")+
  #geom_smooth(data = SLA_total_dis, mapping = aes(x=SLA_pot_mean,y=SLA_Field_means),method = "lm", se = T, color = "black", fill = "grey80", linetype = 2) + 
  geom_point(SLA_total_dis, mapping = aes(x=SLA_pot_mean,y=SLA_Field_means, fill = Origin, color = Origin), size=2.2, pch = 21, color = "black")+
  geom_errorbar(data = dattt,mapping = aes(x = SLA_pot_mean,ymax = SLA_Field_means+field_se, ymin=SLA_Field_means-field_se),width=0.01,size=0.5,alpha = 1, color = "black")+#
  geom_errorbarh(data = dattt,mapping = aes(y = SLA_Field_means,xmax=SLA_pot_mean+pot_se,xmin=SLA_pot_mean-pot_se),height=0.01,size=0.5,alpha = 1, color = "black")+#
  geom_point(data = dattt,mapping = aes(x = SLA_pot_mean, y = SLA_Field_means, fill = Origin, color = Origin),size=3.8, pch = 21, color = "black")+
  scale_color_manual(values=c('#60A7A6','#FEA6A6'))+
  scale_fill_manual(values=c('#60A7A6','#FEA6A6'))+
  guides(col = guide_legend(ncol = 1))+
  geom_abline(intercept=-0.1589007, slope=1.5577529, size=1, linetype = 2)+
  geom_abline(intercept=0,slope=1 ,size=1, linetype = 1, color = "#8B0000")+
  scale_x_continuous(labels = scales::label_comma(accuracy =0.01), limits = c(0.13,0.68)) +
  scale_y_continuous(labels = scales::label_comma(accuracy =0.01), limits = c(0.13,0.68)) + 
  mytheme 
  #ggrepel::geom_text_repel(mapping = aes(x=SLA_pot_mean,y=SLA_Field_means,label=Latin_name), data = df3,size = 2.8,segment.color = "black", 
  #                         color = "black",direction = "both",box.padding = 0.6,
  #                         max.overlaps = getOption("ggrepel.max.overlaps", default = 25), fontface = "italic") +
  #ggrepel::geom_text_repel(mapping = aes(x=SLA_pot_mean,y=SLA_Field_means,label=Latin_name), data = df2,size = 2.8,segment.color = "black", 
  #                         color = "black",direction = "both",box.padding = 0.6,
  #                         max.overlaps = getOption("ggrepel.max.overlaps", default = 25), fontface = "italic") ;p1

#################################### Fig 2B ####################################
### Hmax
trait_data = read.xlsx("Data/Field_traits_mean0831_2.xlsx", sheet = "Field_means", colNames = TRUE, rowNames = FALSE)
colnames(trait_data) <- paste0(colnames(trait_data), "_field")
colnames(trait_data)[1] <- "Species"

trait_data = trait_data %>% left_join(pot_trait, by = "Species") %>% 
  left_join(Common_sp_list[,c(1,2,5)], by = "Species") %>%
  drop_na(Hmax_field) %>% drop_na(Hmax_pot)
rownames(trait_data) = trait_data$Species
length(unique(trait_data$Species))

trait_data$Hmax_field_sqrt = sqrt(trait_data$Hmax_field)
trait_data$Hmax_pot_sqrt = sqrt(trait_data$Hmax_pot)

### Hmax difference in Pot experiment 
dist_mat <- compute_dist_matrix(trait_data['Hmax_pot_sqrt'], metric = 'euclidean',center = TRUE, scale = TRUE)#gower
dist_mat <- standr(dist_mat)
diag(dist_mat) <- NA
mean_vals <- apply(dist_mat, 1, mean, na.rm = TRUE)
sd_vals <- apply(dist_mat, 1, sd, na.rm = TRUE)
se_vals <- sd_vals / sqrt(ncol(dist_mat) - 1)
Hmax_pot_dis <- data.frame(Hmax_pot_mean = mean_vals, Hmax_pot_sd = sd_vals, Hmax_pot_se = se_vals)

### Hmax difference in Field
dist_mat <- compute_dist_matrix(trait_data['Hmax_field_sqrt'], metric = 'euclidean',center = TRUE, scale = TRUE)#gower
dist_mat <- standr(dist_mat)
diag(dist_mat) <- NA
mean_vals <- apply(dist_mat, 1, mean, na.rm = TRUE)
sd_vals <- apply(dist_mat, 1, sd, na.rm = TRUE)
se_vals <- sd_vals / sqrt(ncol(dist_mat) - 1)
Hmax_field_dis <- data.frame(Hmax_Field_means = mean_vals, Hmax_field_sd = sd_vals, Hmax_field_se = se_vals)

Hmax_total_dis = cbind(Hmax_pot_dis, Hmax_field_dis)
Hmax_total_dis$Species = rownames(Hmax_total_dis)
Hmax_total_dis = Hmax_total_dis %>% left_join(Common_sp_list[,c(1,2,5)], by = "Species")
#knitr::kable(head(Hmax_total_dis))

### The effect of origin and types of experiments on DI
Test_total_dis = Hmax_total_dis
rownames(Test_total_dis) = Test_total_dis$Species
###
Hmax_pot_di<-Test_total_dis$Hmax_pot_mean
names(Hmax_pot_di)<- rownames(Test_total_dis)
phytools::phylosig(tree, Hmax_pot_di, method = "K", test = TRUE, nsim =  1000)

Hmax_field_di<-Test_total_dis$Hmax_Field_means
names(Hmax_field_di)<- rownames(Test_total_dis)
phytools::phylosig(tree, Hmax_field_di, method = "K", test = TRUE, nsim =  1000)

###
Origin <-Test_total_dis$Origin
names(Origin)<- rownames(Test_total_dis)

### phylANOVA
aov_di <- phylANOVA(tree, Origin, Hmax_pot_di, nsim = 1000, posthoc=TRUE, p.adj = 'bonferroni')
aov_di
r2(aov_di)
aov_di <- phylANOVA(tree, Origin, Hmax_field_di, nsim = 1000, posthoc=TRUE, p.adj = 'bonferroni')
aov_di
r2(aov_di)

### 
data1 = as.data.frame(Hmax_pot_di)
data1$Type = "pot"
data1$Species = rownames(data1)
colnames(data1)[1] = "Di"
rownames(data1) = NULL

data2 = as.data.frame(Hmax_field_di)
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
mod = relmatLmer(Di ~ Origin*Type + (1|Species), data = data_all_di, relmat = list(Species=phyloMat))
Anova(mod, type="II", test.statistic ="Chisq")
###
fit <- deming(Hmax_Field_means ~ Hmax_pot_mean, ystd=Hmax_field_sd, xstd=Hmax_pot_sd, data=Hmax_total_dis)
print(fit)
Hmax_total_dis$Origin = factor(Hmax_total_dis$Origin, levels = c("Native","Exotic"))

first_char <- substr(Hmax_total_dis$Species, 1, 1)
sub_str <- gsub(".*_", "", Hmax_total_dis$Species)
Latin_name <- paste(first_char, sub_str, sep = ". ")
Hmax_total_dis$Latin_name = Latin_name

df2 = (Hmax_total_dis %>% arrange(desc(Hmax_Field_means)))[c(1:8),]
df3 = (Hmax_total_dis %>% arrange(desc(Hmax_pot_mean)))[c(1:8),]

dat1 = Rmisc::summarySE(Hmax_total_dis, groupvars = c("Origin"), measurevar = c("Hmax_pot_mean"))[,c(1,3,5)]
colnames(dat1)[3] = "pot_se"
dat2 = Rmisc::summarySE(Hmax_total_dis, groupvars = c("Origin"), measurevar = c("Hmax_Field_means"))[,c(3,5)]
colnames(dat2)[2] = "field_se"
dattt = cbind(dat1, dat2)

p2 = ggplot()+
  labs(x = " \n ",
       y= " \n ",
       title = "Maximum height")+
  #geom_smooth(data = Hmax_total_dis, mapping = aes(x=Hmax_pot_mean,y=Hmax_Field_means),method = "lm", se = T, color = "black", fill = "grey80", linetype = 1) + 
  geom_point(data = Hmax_total_dis, mapping = aes(x=Hmax_pot_mean,y=Hmax_Field_means,fill = Origin, color = Origin), size=2.2, pch = 21, color = "black")+
  geom_errorbar(data = dattt,mapping = aes(x = Hmax_pot_mean, ymax = Hmax_Field_means+field_se, ymin=Hmax_Field_means-field_se),width=0.01,size=0.5,alpha = 1, color = "black")+#
  geom_errorbarh(data = dattt,mapping = aes(y = Hmax_Field_means ,xmax=Hmax_pot_mean+pot_se,xmin=Hmax_pot_mean-pot_se),height=0.01,size=0.5,alpha = 1, color = "black")+#
  geom_point(data = dattt,mapping = aes(x = Hmax_pot_mean, y = Hmax_Field_means,fill = Origin, color = Origin),size=3.8, pch = 21, color = "black")+
  scale_color_manual(values=c('#60A7A6','#FEA6A6'))+
  scale_fill_manual(values=c('#60A7A6','#FEA6A6'))+
  guides(col = guide_legend(ncol = 1))+
  geom_abline(intercept=-0.002267806, slope=0.941337257, size=1, linetype = 2)+
  geom_abline(intercept=0,slope=1 ,size=1, linetype = 1, color = "#8B0000")+
  mytheme + 
  scale_x_continuous(labels = scales::label_comma(accuracy =0.01), limits = c(0.16,0.65))+
  scale_y_continuous(labels = scales::label_comma(accuracy =0.01), limits = c(0.16,0.65))
  #ggrepel::geom_text_repel(mapping = aes(x=Hmax_pot_mean,y=Hmax_Field_means,label=Latin_name), data = df2,size = 2.8,segment.color = "black", color = "black",direction = "both",box.padding = 0.6,
  #                         max.overlaps = getOption("ggrepel.max.overlaps", default = 25), fontface = "italic") + 
  #ggrepel::geom_text_repel(mapping = aes(x=Hmax_pot_mean,y=Hmax_Field_means,label=Latin_name), df3,size = 2.8,segment.color = "black", color = "black",direction = "both",box.padding = 0.6,
  #                         max.overlaps = getOption("ggrepel.max.overlaps", default = 25), fontface = "italic");p2


### AGB
trait_data = read.xlsx("Data/Field_traits_mean0831_2.xlsx", sheet = "Field_means", colNames = TRUE, rowNames = FALSE)
colnames(trait_data) <- paste0(colnames(trait_data), "_field")
colnames(trait_data)[1] <- "Species"

trait_data = trait_data %>% left_join(pot_trait, by = "Species") %>%
  left_join(Common_sp_list[,c(1,2,5)], by = "Species") %>%
  drop_na(AGB_field) %>% drop_na(AGB_pot)
rownames(trait_data) = trait_data$Species
length(unique(trait_data$Species))

trait_data$AGB_field_log = log10(trait_data$AGB_field)
trait_data$AGB_pot_log = log10(trait_data$AGB_pot)
### AGB differences in pot experiment
dist_mat <- compute_dist_matrix(trait_data['AGB_pot_log'], metric = 'euclidean',center = TRUE, scale = TRUE)#gower
dist_mat <- standr(dist_mat)
diag(dist_mat) <- NA
mean_vals <- apply(dist_mat, 1, mean, na.rm = TRUE)
sd_vals <- apply(dist_mat, 1, sd, na.rm = TRUE)
se_vals <- sd_vals / sqrt(ncol(dist_mat) - 1)
AGB_pot_dis <- data.frame(AGB_pot_mean = mean_vals, AGB_pot_sd = sd_vals, AGB_pot_se = se_vals)

### AGB differences in Field
dist_mat <- compute_dist_matrix(trait_data['AGB_field_log'], metric = 'euclidean',center = TRUE, scale = TRUE)#gower
dist_mat <- standr(dist_mat)
diag(dist_mat) <- NA
mean_vals <- apply(dist_mat, 1, mean, na.rm = TRUE)
sd_vals <- apply(dist_mat, 1, sd, na.rm = TRUE)
se_vals <- sd_vals / sqrt(ncol(dist_mat) - 1)
AGB_field_dis <- data.frame(AGB_Field_means = mean_vals, AGB_field_sd = sd_vals, AGB_field_se = se_vals)

AGB_total_dis = cbind(AGB_pot_dis, AGB_field_dis)
AGB_total_dis$Species = rownames(AGB_total_dis)
AGB_total_dis = AGB_total_dis %>% left_join(Common_sp_list[,c(1,2,5)], by = "Species")
#knitr::kable(head(AGB_total_dis))

### The effect of origin and types of experiments on DI
Test_total_dis = AGB_total_dis
rownames(Test_total_dis) = Test_total_dis$Species
###
AGB_pot_di<-Test_total_dis$AGB_pot_mean
names(AGB_pot_di)<- rownames(Test_total_dis)
phytools::phylosig(tree, AGB_pot_di, method = "K", test = TRUE, nsim =  1000)

AGB_field_di<-Test_total_dis$AGB_Field_means
names(AGB_field_di)<- rownames(Test_total_dis)
phytools::phylosig(tree, AGB_field_di, method = "K", test = TRUE, nsim =  1000)

###
Origin<-Test_total_dis$Origin
names(Origin)<- rownames(Test_total_dis)

### phylANOVA
aov_di <- phylANOVA(tree, Origin, AGB_pot_di, nsim = 1000, posthoc=TRUE, p.adj = 'bonferroni')
aov_di
r2(aov_di)
aov_di <- phylANOVA(tree, Origin, AGB_field_di, nsim = 1000, posthoc=TRUE, p.adj = 'bonferroni')
aov_di
r2(aov_di)

### 
data1 = as.data.frame(AGB_pot_di)
data1$Type = "pot"
data1$Species = rownames(data1)
colnames(data1)[1] = "Di"
rownames(data1) = NULL

data2 = as.data.frame(AGB_field_di)
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

mod = relmatLmer(Di ~ Origin*Type + (1|Species), data = data_all_di, relmat = list(Species=phyloMat))
Anova(mod, type="II", test.statistic ="Chisq")
###
fit <- deming(AGB_Field_means ~ AGB_pot_mean, ystd=AGB_field_sd, xstd=AGB_pot_sd, data=AGB_total_dis)
print(fit)
AGB_total_dis$Origin = factor(AGB_total_dis$Origin, levels = c("Native","Exotic"))

first_char <- substr(AGB_total_dis$Species, 1, 1)
sub_str <- gsub(".*_", "", AGB_total_dis$Species)
Latin_name <- paste(first_char, sub_str, sep = ". ")
AGB_total_dis$Latin_name = Latin_name

df2 = (AGB_total_dis %>% arrange(desc(AGB_Field_means)))[c(1:8),]
df3 = (AGB_total_dis %>% arrange(desc(AGB_pot_mean)))[c(1:8),]

dat1 = Rmisc::summarySE(AGB_total_dis, groupvars = c("Origin"), measurevar = c("AGB_pot_mean"))[,c(1,3,5)]
colnames(dat1)[3] = "pot_se"
dat2 = Rmisc::summarySE(AGB_total_dis, groupvars = c("Origin"), measurevar = c("AGB_Field_means"))[,c(3,5)]
colnames(dat2)[2] = "field_se"
dattt = cbind(dat1, dat2)

#AGB_total_dis = AGB_total_dis[!AGB_total_dis$Species %in% df2$Species,]
p3 = ggplot()+
  labs(x = "Functional distinctiveness \n estimated in pot experiment",
       y = "Functional distinctiveness \n estimated in field experiment",
       title = "Aboveground biomass")+
  #geom_smooth(data = AGB_total_dis, mapping = aes(x=AGB_pot_mean,y=AGB_Field_means),method = "lm", se = T, color = "black", fill = "grey80", linetype = 2) + 
  geom_errorbar(data = dattt,mapping = aes(x = AGB_pot_mean,ymax = AGB_Field_means+field_se, ymin=AGB_Field_means-field_se),width=0.01,size=0.5,alpha = 1, color = "black")+#
  geom_errorbarh(data = dattt,mapping = aes(y = AGB_Field_means,xmax=AGB_pot_mean+pot_se,xmin=AGB_pot_mean-pot_se),height=0.01,size=0.5,alpha = 1, color = "black")+#
  geom_point(data = AGB_total_dis, mapping = aes(x=AGB_pot_mean,y=AGB_Field_means,fill = Origin, color = Origin),size=2.2, pch = 21, color = "black")+
  geom_point(data = dattt,mapping = aes(x = AGB_pot_mean, y = AGB_Field_means,fill = Origin, color = Origin),size=3.8, pch = 21, color = "black")+
  scale_color_manual(values=c('#60A7A6','#FEA6A6'))+
  scale_fill_manual(values=c('#60A7A6','#FEA6A6'))+
  geom_abline(intercept=0,slope=1 ,size=1, linetype = 1, color = "#8B0000")+
  geom_abline(intercept=0.1471211,slope=0.4719499 ,size=1, linetype = 2)+
  mytheme + 
  scale_x_continuous(labels = scales::label_comma(accuracy =0.01), limits = c(0.16,0.67))+
  scale_y_continuous(labels = scales::label_comma(accuracy =0.01), limits = c(0.16,0.67))
  #ggrepel::geom_text_repel(aes(x=AGB_pot_mean,y=AGB_Field_means,label=Latin_name), data = df2,size = 2.8,segment.color = "black", color = "black",direction = "both",box.padding = 0.6,
  #                         max.overlaps = getOption("ggrepel.max.overlaps", default = 25), fontface = "italic") + 
  #ggrepel::geom_text_repel(aes(x=AGB_pot_mean,y=AGB_Field_means,label=Latin_name), data = df3,size = 2.8,segment.color = "black", color = "black",direction = "both",box.padding = 0.6,
  #                         max.overlaps = getOption("ggrepel.max.overlaps", default = 25), fontface = "italic")

p3


### all traits
trait_data = read.xlsx("Data/Field_traits_mean0831_2.xlsx", sheet = "Field_means", colNames = TRUE, rowNames = FALSE)
colnames(trait_data) <- paste0(colnames(trait_data), "_field")
colnames(trait_data)[1] <- "Species"

trait_data = trait_data %>% left_join(pot_trait, by = "Species") %>% 
  left_join(Common_sp_list[,c(1,2,5)], by = "Species") %>%
  drop_na(AGB_field) %>% drop_na(AGB_pot)
rownames(trait_data) = trait_data$Species
length(unique(trait_data$Species))

###
trait_data$SLA_pot = sqrt(trait_data$SLA_pot)
trait_data$Hmax_pot = sqrt(trait_data$Hmax_pot)
trait_data$AGB_pot = log10(trait_data$AGB_pot)
trait_data$SLA_imp_field = sqrt(trait_data$SLA_imp_field)
trait_data$Hmax_field = sqrt(trait_data$Hmax_field)
trait_data$AGB_field = log10(trait_data$AGB_field)
### 
dist_mat <- compute_dist_matrix(trait_data[ ,c('SLA_pot', 'Hmax_pot', 'AGB_pot')], metric = 'euclidean',center = TRUE, scale = TRUE)#gower
dist_mat <- standr(dist_mat)
diag(dist_mat) <- NA
mean_vals <- apply(dist_mat, 1, mean, na.rm = TRUE)
sd_vals <- apply(dist_mat, 1, sd, na.rm = TRUE)
se_vals <- sd_vals / sqrt(ncol(dist_mat) - 1)
All_pot_dis <- data.frame(All_pot_mean = mean_vals, All_pot_sd = sd_vals, All_pot_se = se_vals)

### 
dist_mat <- compute_dist_matrix(trait_data[ ,c('SLA_imp_field','Hmax_field','AGB_field')], metric = 'euclidean',center = TRUE, scale = TRUE)#gower
dist_mat <- standr(dist_mat)
diag(dist_mat) <- NA
mean_vals <- apply(dist_mat, 1, mean, na.rm = TRUE)
sd_vals <- apply(dist_mat, 1, sd, na.rm = TRUE)
se_vals <- sd_vals / sqrt(ncol(dist_mat) - 1)
All_field_dis <- data.frame(All_Field_means = mean_vals, All_field_sd = sd_vals, All_field_se = se_vals)

All_total_dis = cbind(All_pot_dis, All_field_dis)
All_total_dis$Species = rownames(All_total_dis)
All_total_dis = All_total_dis %>% left_join(Common_sp_list[,c(1,2,5)], by = "Species")
#knitr::kable(head(All_total_dis))

### The effect of origin and types of experiments on DI
Test_total_dis = All_total_dis
rownames(Test_total_dis) = Test_total_dis$Species
###
All_pot_di<-Test_total_dis$All_pot_mean
names(All_pot_di)<- rownames(Test_total_dis)
phytools::phylosig(tree, All_pot_di, method = "K", test = TRUE, nsim =  1000)

All_field_di<-Test_total_dis$All_Field_means
names(All_field_di)<- rownames(Test_total_dis)
phytools::phylosig(tree, All_field_di, method = "K", test = TRUE, nsim =  1000)

### 
data1 = as.data.frame(All_pot_di)
data1$Type = "pot"
data1$Species = rownames(data1)
colnames(data1)[1] = "Di"
rownames(data1) = NULL

data2 = as.data.frame(All_field_di)
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

mod = relmatLmer(Di ~ Origin*Type + (1|Species), data = data_all_di, relmat = list(Species=phyloMat))
Anova(mod, type="II", test.statistic ="Chisq")
###
fit <- deming(All_Field_means ~ All_pot_mean, ystd=All_field_sd, xstd=All_pot_sd, data=All_total_dis)
print(fit)
All_total_dis$Origin = factor(All_total_dis$Origin, levels = c("Native","Exotic"))

first_char <- substr(All_total_dis$Species, 1, 1)
sub_str <- gsub(".*_", "", All_total_dis$Species)
Latin_name <- paste(first_char, sub_str, sep = ". ")
All_total_dis$Latin_name = Latin_name

df2 = (All_total_dis %>% arrange(desc(All_Field_means)))[c(1:8),]
df3 = (All_total_dis %>% arrange(desc(All_pot_mean)))[c(1:8),]
#All_total_dis = All_total_dis[!All_total_dis$Species %in% df2$Species,]

dat1 = Rmisc::summarySE(All_total_dis, groupvars = c("Origin"), measurevar = c("All_pot_mean"))[,c(1,3,5)]
colnames(dat1)[3] = "pot_se"
dat2 = Rmisc::summarySE(All_total_dis, groupvars = c("Origin"), measurevar = c("All_Field_means"))[,c(3,5)]
colnames(dat2)[2] = "field_se"
dattt = cbind(dat1, dat2)

p_all = ggplot()+
  labs(y = " \n ",
       x = "Functional distinctiveness \n estimated in pot experiment",
       title = "The three traits")+
  #geom_smooth(data = All_total_dis, mapping = aes(x=All_pot_mean,y=All_Field_means),method = "lm", se = T, color = "black", fill = "grey80", linetype = 2) + 
  geom_errorbar(data = dattt,mapping = aes(x = All_pot_mean,ymax = All_Field_means+field_se, ymin=All_Field_means-field_se),width=0.01,size=0.5,alpha = 1, color = "black")+#
  geom_errorbarh(data = dattt,mapping = aes(y=All_Field_means,xmax=All_pot_mean+pot_se,xmin=All_pot_mean-pot_se),height=0.01,size=0.5,alpha = 1, color = "black")+#
  geom_point(data = All_total_dis, mapping = aes(x=All_pot_mean,y=All_Field_means,fill = Origin, color = Origin),size=2.2, pch = 21, color = "black")+
  geom_point(data = dattt,mapping = aes(x = All_pot_mean, y = All_Field_means,fill = Origin, color = Origin),size=3.8, pch = 21, color = "black")+
  scale_color_manual(values=c('#60A7A6','#FEA6A6'))+
  scale_fill_manual(values=c('#60A7A6','#FEA6A6'))+
  guides(col = guide_legend(ncol = 1))+
  geom_abline(intercept=0,slope=1 ,size=1, linetype = 1, color = "#8B0000")+
  geom_abline(intercept=0.0605427 ,slope=0.8563063,size=1, linetype = 1)+
  mytheme + 
  scale_y_continuous(labels = scales::label_comma(accuracy =0.01), limits = c(0.2, 0.65))+
  scale_x_continuous(labels = scales::label_comma(accuracy =0.01), limits = c(0.2, 0.65))
  #ggrepel::geom_text_repel(mapping = aes(x=All_pot_mean,y=All_Field_means,label=Latin_name), data = df2,size = 2.8,segment.color = "black", color = "black",direction = "both",box.padding = 0.6,
  #                         max.overlaps = getOption("ggrepel.max.overlaps", default = 25), fontface = "italic") +
  #ggrepel::geom_text_repel(mapping = aes(x=All_pot_mean,y=All_Field_means,label=Latin_name), data = df3,size = 2.8,segment.color = "black", color = "black",direction = "both",box.padding = 0.6,
  #                         max.overlaps = getOption("ggrepel.max.overlaps", default = 25), fontface = "italic")

p_all

### Fig 2
library(cowplot)
plot_grid(p1,p2,p3,p_all,
          labels = c('(a)', '(b)', '(c)', '(d)'))

### Notice that,
### For more picture details, we have further adjusted it in Adobe illustrator.
