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
setwd(getwd())

### Analyzing species lists
Common_sp_list = read.xlsx("Data/Common_species_list.xlsx", sheet = "Common_sp_list", colNames = TRUE, rowNames = FALSE)
Common_sp_list_SLA = c(na.omit(Common_sp_list$Species_SLA))
Common_sp_list_AGB = c(na.omit(Common_sp_list$Species_AGB))

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
phylogenyAux=read.tree("Data/All_species.newick")
to_drop = phylogenyAux$tip.label[!phylogenyAux$tip.label %in% Common_sp_list$Species]
tree <- drop.tip(as.phylo(phylogenyAux), to_drop) 
phyloMat = vcv.phylo(tree)
phyloMat = phyloMat / max(phyloMat)
dim(phyloMat)

### Loading pot experiment database (mean value of traits)
pot_trait = read.xlsx("Data/Pot_traits_mean.xlsx", sheet = "Pot_means", colNames = TRUE, rowNames = FALSE)
colnames(pot_trait) <- paste0(colnames(pot_trait), "_pot")
colnames(pot_trait)[1] <- "Species"
nrow(pot_trait)

### Loading field experiment database (mean value of traits)
### Deming regressive analysis
### Specific leaf area
trait_data = read.xlsx("Data/Field_traits_mean.xlsx", sheet = "Field_means", colNames = TRUE, rowNames = FALSE)
colnames(trait_data) <- paste0(colnames(trait_data), "_field")
colnames(trait_data)[1] <- "Species"
nrow(trait_data)
head(trait_data)
trait_data = trait_data %>% left_join(pot_trait, by = "Species") %>% 
  left_join(Common_sp_list[,c(1,2,5)], by = "Species") %>%
  drop_na(SLA_pot) %>% drop_na(SLA_imp_field)
rownames(trait_data) = trait_data$Species
length(unique(trait_data$Species))

trait_data$SLA_imp_field_log = log10(trait_data$SLA_imp_field)
trait_data$SLA_pot_log = log10(trait_data$SLA_pot)

### Pot experiment SLA differences
dist_mat <- compute_dist_matrix(trait_data['SLA_pot_log'], metric = 'euclidean',center = TRUE, scale = TRUE)
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
dist_mat <- compute_dist_matrix(trait_data['SLA_imp_field_log'], metric = 'euclidean',center = TRUE, scale = TRUE)#gower
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
phytools::phylosig(tree, SLA_pot_di, method = "K", test = TRUE, nsim =  1000)

SLA_field_di<-Test_total_dis$SLA_Field_means
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
summary(mod)
Anova(mod, type="II", test.statistic ="F")

### Deming regressive
fit <- deming(SLA_Field_means ~ SLA_pot_mean, ystd=SLA_field_sd, xstd=SLA_pot_sd, data=SLA_total_dis)
print(fit)
SLA_total_dis$Origin = factor(SLA_total_dis$Origin, levels = c("Native","Exotic"))

first_char <- substr(SLA_total_dis$Species, 1, 1)
sub_str <- gsub(".*_", "", SLA_total_dis$Species)
Latin_name <- paste(first_char, sub_str, sep = ". ")
SLA_total_dis$Latin_name = Latin_name

df2 = (SLA_total_dis %>% arrange(desc(SLA_Field_means)))[c(1:12),]

dat1 = Rmisc::summarySE(SLA_total_dis, groupvars = c("Origin"), measurevar = c("SLA_pot_mean"))[,c(1,3,5)]
colnames(dat1)[3] = "pot_se"
dat2 = Rmisc::summarySE(SLA_total_dis, groupvars = c("Origin"), measurevar = c("SLA_Field_means"))[,c(3,5)]
colnames(dat2)[2] = "field_se"
dattt = cbind(dat1, dat2)

p1 = ggplot(SLA_total_dis,aes(x=SLA_pot_mean,y=SLA_Field_means, fill = Origin, color = Origin))+
  labs(x = " \n ",
       y= "Functional distinctiveness \n estimated in field experiment",
       title = "Specific leaf area")+
  geom_point(size=2.2, pch = 21, color = "black")+
  geom_errorbar(data = dattt,mapping = aes(ymax = SLA_Field_means+field_se, ymin=SLA_Field_means-field_se),width=0.01,size=0.5,alpha = 1, color = "black")+#
  geom_errorbarh(data = dattt,mapping = aes(xmax=SLA_pot_mean+pot_se,xmin=SLA_pot_mean-pot_se),height=0.01,size=0.5,alpha = 1, color = "black")+#
  geom_point(data = dattt,mapping = aes(x = SLA_pot_mean, y = SLA_Field_means),size=3.8, pch = 21, color = "black")+
  scale_color_manual(values=c('#60A7A6','#FEA6A6'))+
  scale_fill_manual(values=c('#60A7A6','#FEA6A6'))+
  guides(col = guide_legend(ncol = 1))+
  geom_abline(intercept=0.1590778, slope=0.2311216, size=1, linetype = 2)+
  geom_abline(intercept=0,slope=1 ,size=1, linetype = 1, color = "#8B0000")+
  scale_x_continuous(labels = scales::label_comma(accuracy =0.01)) +
  scale_y_continuous(labels = scales::label_comma(accuracy =0.01), 
                     limits = c(ggplot_build(p1b)$layout$panel_scales_y[[1]]$range$range)) + 
  mytheme + 
  ggrepel::geom_text_repel(aes(label=Latin_name), df2,size = 2.8,segment.color = "black", 
                           color = "black",direction = "both",box.padding = 0.6,
                           max.overlaps = getOption("ggrepel.max.overlaps", default = 25), 
                           fontface = "italic") ;p1


### Hmax
trait_data = read.xlsx("Data/Field_traits_mean.xlsx", sheet = "Field_means", colNames = TRUE, rowNames = FALSE)
colnames(trait_data) <- paste0(colnames(trait_data), "_field")
colnames(trait_data)[1] <- "Species"

trait_data = trait_data %>% left_join(pot_trait, by = "Species") %>% 
  left_join(Common_sp_list[,c(1,2,5)], by = "Species") %>%
  drop_na(Hmax_field) %>% drop_na(Hmax_pot)
rownames(trait_data) = trait_data$Species
length(unique(trait_data$Species))

trait_data$Hmax_field_log = log10(trait_data$Hmax_field)
trait_data$Hmax_pot_log = log10(trait_data$Hmax_pot)

### Hmax difference in Pot experiment 
dist_mat <- compute_dist_matrix(trait_data['Hmax_pot_log'], metric = 'euclidean',center = TRUE, scale = TRUE)#gower
dist_mat <- standr(dist_mat)
diag(dist_mat) <- NA
mean_vals <- apply(dist_mat, 1, mean, na.rm = TRUE)
sd_vals <- apply(dist_mat, 1, sd, na.rm = TRUE)
se_vals <- sd_vals / sqrt(ncol(dist_mat) - 1)
Hmax_pot_dis <- data.frame(Hmax_pot_mean = mean_vals, Hmax_pot_sd = sd_vals, Hmax_pot_se = se_vals)

### Hmax difference in Field
dist_mat <- compute_dist_matrix(trait_data['Hmax_field_log'], metric = 'euclidean',center = TRUE, scale = TRUE)#gower
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
Origin<-Test_total_dis$Origin
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
mod = relmatLmer(Di ~ Type*Origin + (1|Species), data = data_all_di, relmat = list(Species=phyloMat))
Anova(mod, type="II", test.statistic ="F")

###
fit <- deming(Hmax_Field_means ~ Hmax_pot_mean, ystd=Hmax_field_sd, xstd=Hmax_pot_sd, data=Hmax_total_dis)
print(fit)
Hmax_total_dis$Origin = factor(Hmax_total_dis$Origin, levels = c("Native","Exotic"))

first_char <- substr(Hmax_total_dis$Species, 1, 1)
sub_str <- gsub(".*_", "", Hmax_total_dis$Species)
Latin_name <- paste(first_char, sub_str, sep = ". ")
Hmax_total_dis$Latin_name = Latin_name

dat1 = Rmisc::summarySE(Hmax_total_dis, groupvars = c("Origin"), measurevar = c("Hmax_pot_mean"))[,c(1,3,5)]
colnames(dat1)[3] = "pot_se"
dat2 = Rmisc::summarySE(Hmax_total_dis, groupvars = c("Origin"), measurevar = c("Hmax_Field_means"))[,c(3,5)]
colnames(dat2)[2] = "field_se"
dattt = cbind(dat1, dat2)

df2 = (Hmax_total_dis %>% arrange(desc(Hmax_Field_means)))[c(1:12),]

p2 = ggplot(Hmax_total_dis,aes(x=Hmax_pot_mean,y=Hmax_Field_means,fill = Origin, color = Origin))+
  labs(x = " \n ",
       y= " \n ",
       title = "Maximum height")+
  geom_point(size=2.2, pch = 21, color = "black")+
  geom_errorbar(data = dattt,mapping = aes(ymax = Hmax_Field_means+field_se, ymin=Hmax_Field_means-field_se),width=0.01,size=0.5,alpha = 1, color = "black")+#
  geom_errorbarh(data = dattt,mapping = aes(xmax=Hmax_pot_mean+pot_se,xmin=Hmax_pot_mean-pot_se),height=0.01,size=0.5,alpha = 1, color = "black")+#
  geom_point(data = dattt,mapping = aes(x = Hmax_pot_mean, y = Hmax_Field_means),size=3.8, pch = 21, color = "black")+
  scale_color_manual(values=c('#60A7A6','#FEA6A6'))+
  scale_fill_manual(values=c('#60A7A6','#FEA6A6'))+
  guides(col = guide_legend(ncol = 1))+
  geom_abline(intercept=0.02865052, slope=0.83819137, size=1, linetype = 2)+
  geom_abline(intercept=0,slope=1 ,size=1, linetype = 1, color = "#8B0000")+
  mytheme + 
  scale_x_continuous(labels = scales::label_comma(accuracy =0.01))+
  scale_y_continuous(labels = scales::label_comma(accuracy =0.01),
                     limits = c(ggplot_build(p2b)$layout$panel_scales_y[[1]]$range$range))+
  ggrepel::geom_text_repel(aes(label=Latin_name), df2,size = 2.8,segment.color = "black", color = "black",direction = "both",box.padding = 0.6,
                           max.overlaps = getOption("ggrepel.max.overlaps", default = 25), fontface = "italic") ;p2


### AGB
trait_data = read.xlsx("Data/Field_traits_mean.xlsx", sheet = "Field_means", colNames = TRUE, rowNames = FALSE)
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
Anova(mod, type="II", test.statistic ="F")

###
fit <- deming(AGB_Field_means ~ AGB_pot_mean, ystd=AGB_field_sd, xstd=AGB_pot_sd, data=AGB_total_dis)
print(fit)
AGB_total_dis$Origin = factor(AGB_total_dis$Origin, levels = c("Native","Exotic"))

first_char <- substr(AGB_total_dis$Species, 1, 1)
sub_str <- gsub(".*_", "", AGB_total_dis$Species)
Latin_name <- paste(first_char, sub_str, sep = ". ")
AGB_total_dis$Latin_name = Latin_name

df2 = (AGB_total_dis %>% arrange(desc(AGB_Field_means)))[c(1:12),]

dat1 = Rmisc::summarySE(AGB_total_dis, groupvars = c("Origin"), measurevar = c("AGB_pot_mean"))[,c(1,3,5)]
colnames(dat1)[3] = "pot_se"
dat2 = Rmisc::summarySE(AGB_total_dis, groupvars = c("Origin"), measurevar = c("AGB_Field_means"))[,c(3,5)]
colnames(dat2)[2] = "field_se"
dattt = cbind(dat1, dat2)

#AGB_total_dis = AGB_total_dis[!AGB_total_dis$Species %in% df2$Species,]
p3 = ggplot(AGB_total_dis,aes(x=AGB_pot_mean,y=AGB_Field_means,fill = Origin, color = Origin))+
  labs(x = "Functional distinctiveness \n estimated in pot experiment",
       y = "Functional distinctiveness \n estimated in field experiment",
       title = "Aboveground biomass")+
  geom_errorbar(data = dattt,mapping = aes(ymax = AGB_Field_means+field_se, ymin=AGB_Field_means-field_se),width=0.01,size=0.5,alpha = 1, color = "black")+#
  geom_errorbarh(data = dattt,mapping = aes(xmax=AGB_pot_mean+pot_se,xmin=AGB_pot_mean-pot_se),height=0.01,size=0.5,alpha = 1, color = "black")+#
  geom_point(size=2.2, pch = 21, color = "black")+
  geom_point(data = dattt,mapping = aes(x = AGB_pot_mean, y = AGB_Field_means),size=3.8, pch = 21, color = "black")+
  scale_color_manual(values=c('#60A7A6','#FEA6A6'))+
  scale_fill_manual(values=c('#60A7A6','#FEA6A6'))+
  geom_abline(intercept=0,slope=1 ,size=1, linetype = 1, color = "#8B0000")+
  geom_abline(intercept=0.2205223,slope=0.1847623 ,size=1, linetype = 2)+
  mytheme + 
  scale_x_continuous(labels = scales::label_comma(accuracy =0.01))+
  scale_y_continuous(labels = scales::label_comma(accuracy =0.01))+
  ggrepel::geom_text_repel(aes(label=Latin_name), df2,size = 2.8,segment.color = "black", color = "black",direction = "both",box.padding = 0.6,
                           max.overlaps = getOption("ggrepel.max.overlaps", default = 25), fontface = "italic")

p3


### all traits
trait_data = read.xlsx("Data/Field_traits_mean.xlsx", sheet = "Field_means", colNames = TRUE, rowNames = FALSE)
colnames(trait_data) <- paste0(colnames(trait_data), "_field")
colnames(trait_data)[1] <- "Species"

trait_data = trait_data %>% left_join(pot_trait, by = "Species") %>% 
  left_join(Common_sp_list[,c(1,2,5)], by = "Species") %>%
  drop_na(AGB_field) %>% drop_na(AGB_pot)
rownames(trait_data) = trait_data$Species
length(unique(trait_data$Species))

###
trait_data$SLA_pot = log10(trait_data$SLA_pot)
trait_data$Hmax_pot = log10(trait_data$Hmax_pot)
trait_data$AGB_pot = log10(trait_data$AGB_pot)
trait_data$SLA_imp_field = log10(trait_data$SLA_imp_field)
trait_data$Hmax_field = log10(trait_data$Hmax_field)
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
Anova(mod, type="II", test.statistic ="F")

###
fit <- deming(All_Field_means ~ All_pot_mean, ystd=All_field_sd, xstd=All_pot_sd, data=All_total_dis)
print(fit)
All_total_dis$Origin = factor(All_total_dis$Origin, levels = c("Native","Exotic"))

first_char <- substr(All_total_dis$Species, 1, 1)
sub_str <- gsub(".*_", "", All_total_dis$Species)
Latin_name <- paste(first_char, sub_str, sep = ". ")
All_total_dis$Latin_name = Latin_name

df2 = (All_total_dis %>% arrange(desc(All_Field_means)))[c(1:12),]
#All_total_dis = All_total_dis[!All_total_dis$Species %in% df2$Species,]

dat1 = Rmisc::summarySE(All_total_dis, groupvars = c("Origin"), measurevar = c("All_pot_mean"))[,c(1,3,5)]
colnames(dat1)[3] = "pot_se"
dat2 = Rmisc::summarySE(All_total_dis, groupvars = c("Origin"), measurevar = c("All_Field_means"))[,c(3,5)]
colnames(dat2)[2] = "field_se"
dattt = cbind(dat1, dat2)

p_all = ggplot(All_total_dis,aes(x=All_pot_mean,y=All_Field_means,fill = Origin, color = Origin))+
  labs(y = " \n ",
       x = "Functional distinctiveness \n estimated in pot experiment",
       title = "The three traits")+
  geom_errorbar(data = dattt,mapping = aes(ymax = All_Field_means+field_se, ymin=All_Field_means-field_se),width=0.01,size=0.5,alpha = 1, color = "black")+#
  geom_errorbarh(data = dattt,mapping = aes(xmax=All_pot_mean+pot_se,xmin=All_pot_mean-pot_se),height=0.01,size=0.5,alpha = 1, color = "black")+#
  geom_point(size=2.2, pch = 21, color = "black")+
  geom_point(data = dattt,mapping = aes(x = All_pot_mean, y = All_Field_means),size=3.8, pch = 21, color = "black")+
  scale_color_manual(values=c('#60A7A6','#FEA6A6'))+
  scale_fill_manual(values=c('#60A7A6','#FEA6A6'))+
  guides(col = guide_legend(ncol = 1))+
  geom_abline(intercept=0,slope=1 ,size=1, linetype = 1, color = "#8B0000")+
  geom_abline(intercept=0.1072445 ,slope=0.9052349 ,size=1, linetype = 2)+
  mytheme + 
  scale_y_continuous(labels = scales::label_comma(accuracy =0.01))+
  scale_x_continuous(labels = scales::label_comma(accuracy =0.01))+
  ggrepel::geom_text_repel(aes(label=Latin_name), df2,size = 2.8,segment.color = "black", color = "black",direction = "both",box.padding = 0.6,
                           max.overlaps = getOption("ggrepel.max.overlaps", default = 25), fontface = "italic")

p_all

### Fig 2
library(cowplot)
plot_grid(p1,p2,p3,p_all,
          labels = c('(a)', '(b)', '(c)', '(d)'))

##################################################################################
###################################   Fig 4   ####################################
##################################################################################
### Relationship between functional distinctiveness (muti-traits), persistence 
### and relative abundance within the plot in second year

composition_change <- read.xlsx("Data/Field_composition_database.xlsx",sheet = "Field_composition", rowNames = FALSE, colNames = TRUE)
composition_change$exist_prob <- ifelse(!is.na(composition_change$rebio2020) & composition_change$rebio2020 > 0, 1, 0)
#head(composition_change)
composition_change$Origin = factor(composition_change$Origin, levels = c("Native","Exotic"))
### 
composition_change_common <- composition_change[(composition_change$Species %in% Common_sp_list_AGB), ]
composition_change_common$rebio2020_100 = sqrt(composition_change_common$rebio2020*100)
composition_change_common$Block = as.factor(composition_change_common$Block)
composition_change_common$Plot_num = as.factor(composition_change_common$Plot_num)
length(unique(composition_change_common$Species))
All_total_dis2 = All_total_dis[ ,c("Species","All_pot_mean","All_Field_means")]
names(All_total_dis2)
composition_change_common = composition_change_common %>% left_join(All_total_dis2, by = "Species")

### Remove the missing data of relative biomass in the second year
Com_relative_bio = composition_change_common %>% drop_na(rebio2020)
length(unique(Com_relative_bio$Species))

### based on pot experiment
Com_relative_bio2 = Com_relative_bio
length(unique(Com_relative_bio2$Species))
Com_relative_bio2$Block = as.factor(Com_relative_bio2$Block)
Com_relative_bio2$Plot_num = as.factor(Com_relative_bio2$Plot_num)

library(lme4)
library(lmerTest)
mod12 <- lmer((rebio2020_100) ~ All_pot_mean +  (1|Block/Plot_num), data=Com_relative_bio2)
anova(mod12)
summary(mod12)
MuMIn::r.squaredGLMM(mod12)
Com_relative_bio2$F0 <- predictSE(mod12, Com_relative_bio2, level = 0)$fit
Com_relative_bio2$SE <- predictSE(mod12, Com_relative_bio2, level = 0)$se.fit

p11 = ggplot(Com_relative_bio2, aes(x=(All_pot_mean), y=(rebio2020_100))) +
  geom_point(aes(fill = Origin),size = 2.2, pch = 21) + # fill = "#00000022", color = "#0A0A0A"
  geom_line(aes(y=F0), size=1, linetype = 1) + #color = "#8B0000"
  mytheme + 
  scale_fill_manual(values=c('#60A7A6','#FEA6A6'))+
  scale_x_continuous(labels = scales::label_comma(accuracy =0.01))+
  scale_y_continuous(labels = scales::label_comma(accuracy =0.01))+
  ylab('Relative abundance within the \n plot in second year (%, sqrt)') + 
  xlab('Functional distinctiveness \n measured in pot experiment') +
  geom_line(aes(y =  F0 - 1.96 * SE), color = "grey40", linetype = "dashed", size=0.6) +
  geom_line(aes(y =  F0 + 1.96 * SE), color = "grey40", linetype = "dashed", size=0.6) +
  annotate('text',x = 0.5,y=7.5,label=('R2 = 0.078 \n p = 0.006'),size=4,color='black')

p11

###
library(lme4)
library(AICcmodavg)
mod12 <- glmer(exist_prob ~ (All_pot_mean) + (1|Block/Plot_num), na.action=na.omit, family=binomial, data=composition_change_common)
table(composition_change_common$Origin)
Anova(mod12)
summary(mod12)
MuMIn::r.squaredGLMM(mod12)
composition_change_common$F0 = predictSE(mod12, composition_change_common, level = 0)$fit
composition_change_common$SE <- predictSE(mod12, composition_change_common, level = 0)$se.fit

p22 = ggplot(composition_change_common, aes(x=(All_pot_mean), y=exist_prob)) +
  geom_point(aes(fill = Origin),size = 2.2, pch = 21) + # fill = "#00000022", color = "#0A0A0A"
  geom_line(aes(y=F0), size=1, linetype = 2) + #color = "#8B0000"
  ylab('Probability of persistence') + 
  xlab('Functional distinctiveness \n in pot experiment') +
  scale_x_continuous(labels = scales::label_comma(accuracy =0.01))+
  scale_y_continuous(labels = scales::label_comma(accuracy =0.01))+
  mytheme +   
  scale_fill_manual(values=c('#60A7A6','#FEA6A6'))+
  geom_line(aes(y =  F0 - 1.96 * SE), color = "grey40", linetype = "dashed", size=0.6) +
  geom_line(aes(y =  F0 + 1.96 * SE), color = "grey40", linetype = "dashed", size=0.6) +
  annotate('text',x = 0.5,y=0.75,label=('Z = 0.99 \n R2 = 0.004 \n p = 0.321'),size=4,color='black')

p22

### based on field experiment
Com_relative_bio = composition_change_common %>% drop_na(rebio2020)
length(unique(Com_relative_bio$Species))

### 田间实验
Com_relative_bio2 = Com_relative_bio
length(unique(Com_relative_bio2$Species))

mod12 <- lmer((rebio2020_100) ~ All_Field_means +  (1|Block/Plot_num), data=Com_relative_bio2)
anova(mod12)
MuMIn::r.squaredGLMM(mod12)

Com_relative_bio2$F0 = predictSE(mod12, Com_relative_bio2, level = 0)$fit
Com_relative_bio2$SE <- predictSE(mod12, Com_relative_bio2, level = 0)$se.fit

p33 = ggplot(Com_relative_bio2, aes(x=(All_Field_means), y=(rebio2020_100))) +
  geom_point(aes(fill = Origin),size = 2.2, pch = 21) + # fill = "#00000022", color = "#0A0A0A"
  geom_line(aes(y=F0), size=1, linetype = 1) + #color = "#8B0000"
  mytheme + 
  scale_fill_manual(values=c('#60A7A6','#FEA6A6'))+
  scale_x_continuous(labels = scales::label_comma(accuracy =0.01))+
  scale_y_continuous(labels = scales::label_comma(accuracy =0.01),
                     limits = c(ggplot_build(p11)$layout$panel_scales_y[[1]]$range$range))+
  ylab(' \n ') + 
  xlab('Functional distinctiveness \n measured in field experiment') +
  geom_line(aes(y =  F0 - 1.96 * SE), color = "grey40", linetype = "dashed", size=0.6) +
  geom_line(aes(y =  F0 + 1.96 * SE), color = "grey40", linetype = "dashed", size=0.6) +
  annotate('text',x = 0.35,y=7.5,label=('R2 = 0.089 \n p = 0.003'),size=4,color='black')

p33

###
mod12 <- glmer(exist_prob ~ (All_Field_means) + (1|Block/Plot_num), na.action=na.omit, family=binomial, data=composition_change_common,
               control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
table(composition_change_common$Origin)
Anova(mod12)
summary(mod12)
MuMIn::r.squaredGLMM(mod12)
composition_change_common$F0 = predictSE(mod12, composition_change_common, level = 0)$fit
composition_change_common$SE <- predictSE(mod12, composition_change_common, level = 0)$se.fit

p44 = ggplot(composition_change_common, aes(x=(All_Field_means), y=exist_prob)) +
  geom_point(aes(fill = Origin),size = 2.2, pch = 21) + # fill = "#00000022", color = "#0A0A0A"
  geom_line(aes(y=F0), size=1, linetype = 1) + #color = "#8B0000"
  ylab(' \n ') + 
  xlab('Functional distinctiveness \n in field experiment') +
  mytheme +   
  scale_x_continuous(labels = scales::label_comma(accuracy =0.01))+
  scale_y_continuous(labels = scales::label_comma(accuracy =0.01))+
  scale_fill_manual(values=c('#60A7A6','#FEA6A6'))+
  geom_line(aes(y =  F0 - 1.96 * SE), color = "grey40", linetype = "dashed", size=0.6) +
  geom_line(aes(y =  F0 + 1.96 * SE), color = "grey40", linetype = "dashed", size=0.6) +
  annotate('text',x = 0.4,y=0.75,label=('Z = 4.56 \n R2 = 0.097 \n p < 0.001'),size=4,color='black')

p44

### Fig 4
p22+p44+p11+p33 +
  plot_layout(ncol = 2, nrow = 2) +
  plot_annotation(tag_levels = 'a')


################################################################################
################################## Fig S2 ######################################
################################################################################
### Correlation of measurement values under two planting scenarios for the same trait
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

### Loading pot experiment database 
pot_trait = read.xlsx("Data/Pot_traits_mean.xlsx", sheet = "Pot_means", colNames = TRUE, rowNames = FALSE)
colnames(pot_trait) <- paste0(colnames(pot_trait), "_pot")
colnames(pot_trait)[1] <- "Species"

### Loading field experiment database 
### AGB
trait_data = read.xlsx("Data/Field_traits_mean.xlsx", sheet = "Field_means", colNames = TRUE, rowNames = FALSE)
colnames(trait_data) <- paste0(colnames(trait_data), "_field")
colnames(trait_data)[1] <- "Species"

trait_data = trait_data %>% left_join(pot_trait, by = "Species") %>% 
  left_join(Common_sp_list[,c(1,2,5)], by = "Species") %>%
  drop_na(AGB_field) %>% drop_na(AGB_pot)
length(unique(trait_data$Species))

trait_data$Origin = factor(trait_data$Origin, levels = c("Native", "Exotic"))
ma.test <- sma(AGB_field ~ AGB_pot, log='xy', slope.test=1, data = trait_data, type = "shift", na.action = na.omit)
ma.test$n
summary(ma.test)

first_char <- substr(trait_data$Species, 1, 1)
sub_str <- gsub(".*_", "", trait_data$Species)
Latin_name <- paste(first_char, sub_str, sep = ". ")
trait_data$Latin_name = Latin_name

df2 = (trait_data %>% arrange(desc(AGB_field)))[c(1:6),]
df3 = (trait_data %>% arrange((AGB_field)))[c(1:6),]

P_AGB = ggplot(trait_data, mapping = aes(x = log10(AGB_pot) , y = log10(AGB_field))) + 
  geom_point(trait_data, mapping = aes(x = log10(AGB_pot) , y = log10(AGB_field), shape = Origin, fill = Origin),size=2.2, color = "black")+
  scale_shape_manual(values = c(21,21))+
  scale_color_manual(values=c('#60A7A6','#FEA6A6'))+
  scale_fill_manual(values=c('#60A7A6','#FEA6A6'))+
  labs(x = ' \n ',
       y = ' \n ',
       title = "Aboveground biomass") + 
  mytheme + 
  scale_x_continuous(labels = scales::label_comma(accuracy =0.1)) +
  scale_y_continuous(labels = scales::label_comma(accuracy =0.1)) + 
  geom_abline(intercept=1.427132 ,slope=1.623992 ,size=.8, linetype = 1)+
  ggrepel::geom_text_repel(aes(label=Latin_name), df2,size = 2.5,segment.color = "black", 
                           color = "black",direction = "both",box.padding = 0.7, fontface = "italic",
                           max.overlaps = getOption("ggrepel.max.overlaps", default = 20))+
  ggrepel::geom_text_repel(aes(label=Latin_name), df3,size = 2.5,segment.color = "black", 
                           color = "black",direction = "both",box.padding = 0.7, fontface = "italic",
                           max.overlaps = getOption("ggrepel.max.overlaps", default = 20))+
  annotate('text', label = expression(italic(R)^2~"="~0.198*","~italic(p)~"<"~0.001), x = -0.6, y = 2.5, size = 3, parse = TRUE)

P_AGB

### Hmax
trait_data = read.xlsx("Data/Field_traits_mean.xlsx", sheet = "Field_means", colNames = TRUE, rowNames = FALSE)
colnames(trait_data) <- paste0(colnames(trait_data), "_field")
colnames(trait_data)[1] <- "Species"

trait_data = trait_data %>% left_join(pot_trait, by = "Species") %>% 
  left_join(Common_sp_list[,c(1,2,5)], by = "Species") %>%
  drop_na(Hmax_field) %>% drop_na(Hmax_pot)
length(unique(trait_data$Species))

trait_data$Origin = factor(trait_data$Origin, levels = c("Native", "Exotic"))
summary(lm(log10(Hmax_field)~ log10(Hmax_pot), data = trait_data))
ma.test <- sma(Hmax_field ~ Hmax_pot, log='xy', slope.test=1, data = trait_data, type = "shift", na.action = na.omit)
ma.test$n
summary(ma.test)
#plot(ma.test)
first_char <- substr(trait_data$Species, 1, 1)
sub_str <- gsub(".*_", "", trait_data$Species)
Latin_name <- paste(first_char, sub_str, sep = ". ")
trait_data$Latin_name = Latin_name

df2 = (trait_data %>% arrange(desc(Hmax_field)))[c(1:6),]
df3 = (trait_data %>% arrange((Hmax_field)))[c(1:6),]

P_Hmax = ggplot(trait_data, mapping = aes(x = log10(Hmax_pot) , y = log10(Hmax_field))) + 
  geom_point(aes(shape = Origin, fill = Origin),size=2.2, color = "black")+
  scale_shape_manual(values = c(21,21))+
  scale_color_manual(values=c('#60A7A6','#FEA6A6'))+
  scale_fill_manual(values=c('#60A7A6','#FEA6A6'))+
  labs(y = ' \n ',
       x = bquote(atop('Traits value'~(log[10]), 'estimated in pot experiment')),
       title = "Maximum height") + 
  geom_abline(intercept=1.1498689 ,slope=0.6356759 ,size=.8, linetype = 1)+
  mytheme + 
  ggrepel::geom_text_repel(aes(label=Latin_name), df2,size = 2.5,segment.color = "black",
                           color = "black",direction = "both",box.padding = 0.7, fontface = "italic",
                           max.overlaps = getOption("ggrepel.max.overlaps", default = 15))+
  ggrepel::geom_text_repel(aes(label=Latin_name), df3,size = 2.5,segment.color = "black",
                           color = "black",direction = "both",box.padding = 0.7, fontface = "italic",
                           max.overlaps = getOption("ggrepel.max.overlaps", default = 15))+
  annotate('text', label = expression(italic(R)^2~"="~0.338*","~italic(p)~"<"~0.001), x = 0.95, y = 2.2, size = 3, parse = TRUE)


P_Hmax

### SLA (row data)
trait_data = read.xlsx("Data/Field_traits_mean.xlsx", sheet = "Field_means", colNames = TRUE, rowNames = FALSE)
colnames(trait_data) <- paste0(colnames(trait_data), "_field")
colnames(trait_data)[1] <- "Species"

trait_data = trait_data %>% left_join(pot_trait, by = "Species") %>% 
  left_join(Common_sp_list[,c(1,2,5)], by = "Species") %>%
  drop_na(SLA_imp_field) %>% drop_na(SLA_pot) 
length(unique(trait_data$Species))

trait_data$Origin = factor(trait_data$Origin, levels = c("Native", "Exotic"))
summary(lm(log10(SLA_field)~ log10(SLA_pot), data = trait_data))
ma.test <- sma(SLA_field ~ SLA_pot, log='xy', slope.test=1, data = trait_data, type = "shift", na.action = na.omit)
ma.test$n
summary(ma.test)
#plot(ma.test)
first_char <- substr(trait_data$Species, 1, 1)
sub_str <- gsub(".*_", "", trait_data$Species)
Latin_name <- paste(first_char, sub_str, sep = ". ")
trait_data$Latin_name = Latin_name

df2 = (trait_data %>% arrange(desc(SLA_field)))[c(1:6),]
df3 = (trait_data %>% arrange((SLA_field)))[c(1:6),]

P_SLA = ggplot(trait_data, mapping = aes(x = log10(SLA_pot) , y = log10(SLA_field))) + 
  geom_point(trait_data, mapping = aes(x = log10(SLA_pot) , y = log10(SLA_field), shape = Origin, fill = Origin),size=2.2, color = "black")+
  scale_shape_manual(values = c(21,21))+
  scale_color_manual(values=c('#60A7A6','#FEA6A6'))+
  scale_fill_manual(values=c('#60A7A6','#FEA6A6'))+
  labs(x = ' \n ',
       y = bquote(atop('Traits value'~(log[10]), 'estimated in field experiment')),
       title = "Specific leaf area") + 
  geom_abline(intercept=1.353949 ,slope=0.4013532 ,size=.8, linetype = 1)+
  mytheme + 
  ggrepel::geom_text_repel(aes(label=Latin_name), df2,size = 2.5,segment.color = "black", 
                           color = "black",direction = "both",box.padding = 0.7, fontface = "italic",
                           max.overlaps = getOption("ggrepel.max.overlaps", default = 15)) +
  ggrepel::geom_text_repel(aes(label=Latin_name), df3,size = 2.5,segment.color = "black", 
                           color = "black",direction = "both",box.padding = 0.7,fontface = "italic",
                           max.overlaps = getOption("ggrepel.max.overlaps", default = 15)) +
  annotate('text', label = expression(italic(R)^2~"="~0.083*","~italic(p)~"="~0.048), x = 1.8, y = 2.2, size = 3, parse = TRUE)



P_SLA

### (Inferred data)
summary(lm(log10(SLA_imp_field)~ log10(SLA_pot), data = trait_data))
ma.test <- sma(SLA_imp_field ~ SLA_pot, log='xy', slope.test=1, data = trait_data, type = "shift", na.action = na.omit)
ma.test$n
summary(ma.test)
#plot(ma.test)

df2 = (trait_data %>% arrange(desc(SLA_imp_field)))[c(1:6),]
df3 = (trait_data %>% arrange((SLA_imp_field)))[c(1:6),]

P_SLA = ggplot(trait_data, mapping = aes(x = log10(SLA_pot) , y = log10(SLA_imp_field))) + 
  geom_point(trait_data, mapping = aes(x = log10(SLA_pot) , y = log10(SLA_imp_field), 
                                       shape = Origin, fill = Origin),size=2.2, color = "black")+
  scale_shape_manual(values = c(21,21))+
  scale_color_manual(values=c('#60A7A6','#FEA6A6'))+
  scale_fill_manual(values=c('#60A7A6','#FEA6A6'))+
  labs(x = ' \n ',
       y = bquote(atop('Traits value'~(log[10]), 'estimated in field experiment')),
       title = "Specific leaf area") + 
  geom_abline(intercept=1.542867 ,slope=0.3206351 ,size=.8, linetype = 1)+
  mytheme + 
  ggrepel::geom_text_repel(aes(label=Latin_name), df2,size = 2.5,segment.color = "black", 
                           color = "black",direction = "both",box.padding = 0.7, fontface = "italic",
                           max.overlaps = getOption("ggrepel.max.overlaps", default = 100)) +
  ggrepel::geom_text_repel(aes(label=Latin_name), df3,size = 2.5,segment.color = "black",
                           color = "black",direction = "both",box.padding = 0.7, fontface = "italic",
                           max.overlaps = getOption("ggrepel.max.overlaps", default = 100)) + 
  annotate('text', label = expression(italic(R)^2~"="~0.098*","~italic(p)~"="~0.012), x = 1.5, y = 2.3, size = 3, parse = TRUE)


P_SLA

### Fig S2
P_SLA+P_Hmax+P_AGB+plot_annotation(tag_levels = 'a')


### Other related analysis
################################################################################
### Comparison of survival rate between exotic and native plants 
### in the first year of field experiment
suv <- read.xlsx("Data/Survival_database_in_field.xlsx",sheet = "Survival", rowNames = FALSE, colNames = TRUE)
rownames(suv) = suv$Species
phylogenyAux <- read.tree("Data/All_species.newick")
to_drop = setdiff(as.vector(phylogenyAux$tip.label),as.vector((suv$Species)))
phy_tree <- drop.tip(as.phylo(phylogenyAux), to_drop) 

###
suv2<-suv$Survival
names(suv2)<- rownames(suv)
phytools::phylosig(phy_tree, suv2, method = "K", test = TRUE, nsim =  1000)
###
Origin<-suv$Origin
names(Origin)<- rownames(suv)

### phylANOVA
aov_di <- phylANOVA(phy_tree, Origin, suv2, nsim = 1000, posthoc=TRUE, p.adj = 'bonferroni')
aov_di
Rmisc::summarySE(suv, measurevar = "Survival", groupvars = "Origin")

### Notice that,
### For more picture details, we have further adjusted it in Adobe illustrator.