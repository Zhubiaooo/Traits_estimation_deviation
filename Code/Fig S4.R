################################################################################
#################################   Fig S4   ###################################
################################################################################
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
summary(lm(sqrt(SLA_field)~ sqrt(SLA_green), data = trait_data))

first_char <- substr(trait_data$Species, 1, 1)
sub_str <- gsub(".*_", "", trait_data$Species)
Latin_name <- paste(first_char, sub_str, sep = ". ")
trait_data$Latin_name = Latin_name

df2 = (trait_data %>% arrange(desc(SLA_field)))[c(1:6),]
df3 = (trait_data %>% arrange((SLA_field)))[c(1:6),]

p1 = ggplot(trait_data, mapping = aes(x = sqrt(SLA_green) , y = sqrt(SLA_field))) + 
  geom_smooth(data = trait_data, mapping = aes(x=sqrt(SLA_green),y=sqrt(SLA_field)),method = "lm", se = T, color = "black", fill = "grey80", linetype = 1) + 
  geom_point(trait_data, mapping = aes(x = sqrt(SLA_green) , y = sqrt(SLA_field), shape = Origin, fill = Origin),size=2.2, color = "black")+
  scale_shape_manual(values = c(21,21))+
  scale_color_manual(values=c('#005097','#FDB435'))+
  scale_fill_manual(values=c('#005097','#FDB435'))+
  scale_x_continuous(labels = scales::label_comma(accuracy =0.1)) +
  scale_y_continuous(labels = scales::label_comma(accuracy =0.1)) + 
  labs(x = 'Specific leaf area (cm2/g, sqrt) \nestimated in pot experiment',
       y = 'Specific leaf area (cm2/g, sqrt) \nestimated in field experiment') + 
  #geom_abline(intercept=11.63547 ,slope=0.12705 ,size=.8, linetype = 1)+
  mytheme + 
  annotate('text',x = 20,y=10.5,label=('R2 = 0.229 \n p < 0.001'),size=4,color='black') +
  ggrepel::geom_text_repel(aes(label=Latin_name), df2,size = 2.5,segment.color = "black", color = "black",direction = "both",box.padding = 0.7,
                           max.overlaps = getOption("ggrepel.max.overlaps", default = 100)) +
  ggrepel::geom_text_repel(aes(label=Latin_name), df3,size = 2.5,segment.color = "black", color = "black",direction = "both",box.padding = 0.7,
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

trait_data$SLA_field_log = sqrt(trait_data$SLA_field)
trait_data$SLA_green_log = sqrt(trait_data$SLA_green)

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
phylogenyAux = read.tree("Data/iq_tree.treefile")
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

mod = relmatLmer(Di ~ Origin*Type + (1|Species), data = data_all_di, relmat = list(Species=phyloMat))
Anova(mod, type="II", test.statistic ="Chisq")
Rmisc::summarySE(data_all_di, measurevar = "Di", groupvars = c("Type"))

###
colnames(SLA_total_dis)
dat1 = Rmisc::summarySE(SLA_total_dis, groupvars = c("Origin"), measurevar = c("SLA_green_mean"))[,c(1,3,5)]
colnames(dat1)[3] = "pot_se"
dat2 = Rmisc::summarySE(SLA_total_dis, groupvars = c("Origin"), measurevar = c("SLA_field_mean"))[,c(3,5)]
colnames(dat2)[2] = "field_se"
dattt = cbind(dat1, dat2)

first_char <- substr(SLA_total_dis$Species, 1, 1)
sub_str <- gsub(".*_", "", SLA_total_dis$Species)
Latin_name <- paste(first_char, sub_str, sep = ". ")
SLA_total_dis$Latin_name = Latin_name

df2 = (SLA_total_dis %>% arrange(desc(SLA_field_mean)))[c(1:10),]
df3 = (SLA_total_dis %>% arrange(desc(SLA_green_mean)))[c(1:10),]

p2 = ggplot()+
  #geom_smooth(data = SLA_total_dis, mapping = aes(x=SLA_green_mean,y=SLA_field_mean),method = "lm", se = T, color = "black", fill = "grey80", linetype = 2) + 
  labs(x = "Functional distinctiveness\nestimated in pot experiment",
       y="Functional distinctiveness\nestimated in field experiment")+
  geom_point(SLA_total_dis,mapping = aes(x=SLA_green_mean,y=SLA_field_mean,fill = Origin, color = Origin, shape = Origin),size=2.2, pch = 21, color = "black")+
  geom_errorbar(data = dattt,mapping = aes(x = SLA_green_mean, ymax = SLA_field_mean+field_se, ymin=SLA_field_mean-field_se),width=0.01,size=0.5,alpha = 1, color = "black")+#
  geom_errorbarh(data = dattt,mapping = aes(y = SLA_field_mean,xmax=SLA_green_mean+pot_se,xmin=SLA_green_mean-pot_se),height=0.01,size=0.5,alpha = 1, color = "black")+#
  geom_point(data = dattt,mapping = aes(x = SLA_green_mean, y = SLA_field_mean,fill = Origin, color = Origin),size=3.8, pch = 21, color = "black")+
  scale_shape_manual(values = c(21,21))+
  scale_color_manual(values=c('#005097','#FDB435'))+
  scale_fill_manual(values=c('#005097','#FDB435'))+
  scale_x_continuous(labels = scales::label_comma(accuracy =0.1), limits = c(0.15,0.69)) +
  scale_y_continuous(labels = scales::label_comma(accuracy =0.1), limits = c(0.15,0.69)) + 
  guides(col = guide_legend(ncol = 1))+
  geom_abline(intercept=0.04950864, slope=0.81094810 ,size=1, linetype = 2)+
  geom_abline(intercept=0,slope=1 ,size=1, linetype = 1, color = "#95373B")+
  mytheme + 
  #annotate('text',x = 3,y=1,label=('Slope = 0.06054834 , \n95% CI [-55.23847, 55.35957]'),size=4,color='black')+
  ggrepel::geom_text_repel(aes(x=SLA_green_mean,y=SLA_field_mean,label=Latin_name), df2,size = 2.5,segment.color = "black", color = "black",direction = "both",box.padding = 0.6,
                           max.overlaps = getOption("ggrepel.max.overlaps", default = 100)) +
  ggrepel::geom_text_repel(aes(x=SLA_green_mean,y=SLA_field_mean,label=Latin_name), df3,size = 2.5,segment.color = "black", color = "black",direction = "both",box.padding = 0.6,
                           max.overlaps = getOption("ggrepel.max.overlaps", default = 100))

p2

library(openxlsx)
library(lmerTest)
library(lme4)
library(MuMIn)
library(car)
library(AICcmodavg)
library(ggplot2)
field_com_data = read.xlsx("Data/all_row_data.xlsx", sheet = "field_data_mean", rowNames = F, colNames = T)
head(field_com_data)

### 
Field_trait = read.xlsx("Data/Field_traits_mean.xlsx", sheet = "Field_means", colNames = TRUE, rowNames = FALSE)
colnames(Field_trait) <- paste0( "Field_", colnames(Field_trait))
colnames(Field_trait)[1] <- "Species"

### 
Pot_trait = read.xlsx("Data/Pot_traits_mean.xlsx", sheet = "Pot_means", colNames = TRUE, rowNames = FALSE)
colnames(Pot_trait) <- paste0("Pot_", colnames(Pot_trait))
colnames(Pot_trait)[1] <- "Species"
nrow(Pot_trait)

### 
field_com_data = field_com_data %>% left_join(Field_trait[,c(1:5)], by = "Species") %>%
  left_join(Pot_trait, by = "Species")
colnames(field_com_data)
field_com_data$rebio2022_100 = log10(field_com_data$rebio2022*100)
field_com_data$Block = as.factor(field_com_data$Block)
length(unique(field_com_data$Species))
field_com_data$Seed_source <- factor(field_com_data$Seed_source, levels = rev(c("Guangdong","Guangxi","Hunan","Hubei","Henan","Shandong")))

#### 
field_com_data$Field_SLA_sqrt = sqrt(field_com_data$Field_SLA)
field_com_data$Field_SLA_imp_sqrt = sqrt(field_com_data$Field_SLA_imp)
field_com_data$Field_Hmax_sqrt = sqrt(field_com_data$Field_Hmax)
field_com_data$Field_AGB_lg = log10(field_com_data$Field_AGB)
##
field_com_data$Pot_SLA_sqrt = sqrt(field_com_data$Pot_SLA)
field_com_data$Pot_Hmax_sqrt = sqrt(field_com_data$Pot_Hmax)
field_com_data$Pot_AGB_lg = log10(field_com_data$Pot_AGB)
##
field_com_data$exist_prob <- ifelse(!is.na(field_com_data$rebio2022) & field_com_data$rebio2022 > 0, 1, 0)
field_com_data$Origin = factor(field_com_data$Origin, levels = c("Native","Exotic"))

#### 
field_com_data = field_com_data %>% tidyr::drop_na(Field_SLA)
length(unique(field_com_data$Species))

################################################################################
##### odds of persistence
Seed_source = unique(field_com_data$Seed_source)
all_summary_glmer = NULL
for (i in Seed_source) {
  select_data = subset(field_com_data, Seed_source == i)
  mod1 <- glmer(exist_prob ~ Field_SLA_sqrt + (1|Block) , na.action=na.omit, family=binomial, data=select_data)
  mod1_summary = data.frame(Pop = i, Traits = "SLA", N = nrow(select_data), sp_num = length(unique(select_data$Species)),
                            Chisq = as.data.frame(Anova(mod1))[,1], p_value = as.data.frame(Anova(mod1))[,3],
                            R2m = r.squaredGLMM(mod1)[1,1], R2c = r.squaredGLMM(mod1)[1,2])
  all_summary_glmer = rbind(all_summary_glmer, mod1_summary)
}

##### 
colnames(field_com_data)
mod1 = glmer(exist_prob ~ Field_SLA_sqrt + (1|Block/Plot_num), na.action=na.omit, family=binomial, data=field_com_data)
summary(mod1)
Anova(mod1); r.squaredGLMM(mod1)[1,]

nrow(field_com_data); length(unique(field_com_data$Species))

##### relative abundance of species in second year
Seed_source = unique(field_com_data$Seed_source)
all_summary_lmer = NULL
for (i in Seed_source) {
  select_data = subset(field_com_data, Seed_source == i) %>% tidyr::drop_na(rebio2022_100)
  mod1 <- lmer(rebio2022_100 ~ Field_SLA_sqrt + (1|Block) , na.action=na.omit, data=select_data)
  mod1_summary = data.frame(Pop = i, Traits = "SLA", N = nrow(select_data), sp_num = length(unique(select_data$Species)),
                            Chisq = as.data.frame(Anova(mod1))[,1], p_value = as.data.frame(Anova(mod1))[,3],
                            R2m = r.squaredGLMM(mod1)[1,1], R2c = r.squaredGLMM(mod1)[1,2])
  all_summary_lmer = rbind(all_summary_lmer, mod1_summary)
}
all_summary_lmer$p_value = round(all_summary_lmer$p_value, 3)

################################################################################
colnames(field_com_data)
mod1 = lmer(rebio2022_100 ~ Field_SLA_sqrt + (1|Block/Plot_num) , na.action=na.omit, data=field_com_data)
Anova(mod1); r.squaredGLMM(mod1)[1,]

select_data = field_com_data %>% tidyr::drop_na(rebio2022_100)
nrow(select_data); length(unique(select_data$Species))


######### Visualization
length(unique(field_com_data$Species))
unique(field_com_data$Seed_source)

#### SLA
mod12 <- glmer(exist_prob ~ Field_SLA_sqrt*Seed_source + (1|Block) , na.action=na.omit, family=binomial, data=field_com_data)
Anova(mod12)
summary(mod12)
MuMIn::r.squaredGLMM(mod12)

field_com_data$F0 = predictSE(mod12, field_com_data, level = 0)$fit
field_com_data$SE <- predictSE(mod12, field_com_data, level = 0)$se.fit

####
field_com_data_all = field_com_data
mod12 <- glmer(exist_prob ~ Field_SLA_sqrt + (1|Block/Plot_num) , na.action=na.omit, family=binomial, data=field_com_data_all)
Anova(mod12)
r.squaredGLMM(mod12)
field_com_data_all$F0 = predictSE(mod12, field_com_data_all, level = 0)$fit
field_com_data_all$SE <- predictSE(mod12, field_com_data_all, level = 0)$se.fit

ggplot(field_com_data, aes(x=Field_SLA_sqrt, y=exist_prob)) +
  #geom_ribbon(aes(x=Field_SLA_sqrt,ymin = F0 - 1.96 * SE,ymax = F0 + 1.96 * SE, fill = Seed_source), alpha = I(0.1)) +
  geom_line(data = subset(field_com_data, Seed_source != "Hubei"), mapping = aes(x=Field_SLA_sqrt, y=F0, color = Seed_source), size=1, linetype = 2) +
  geom_line(data = subset(field_com_data, Seed_source == "Hubei"), mapping = aes(x=Field_SLA_sqrt, y=F0, color = Seed_source), size=1, linetype = 1) +
  geom_line(data = field_com_data_all, mapping = aes(x=Field_SLA_sqrt, y=F0), size=1.5, linetype = 1, color = "black") +
  geom_ribbon(data = field_com_data_all, mapping = aes(x=Field_SLA_sqrt, ymin = F0 - 1.96 * SE,ymax = F0 + 1.96 * SE, fill = "#EBEBEB"), alpha = I(0.1)) +
  geom_point(aes(fill = Seed_source),size = 2.2, pch = 21) + 
  labs(x = "Specific leaf area (cm2/g, log10)\nestimated in field experiment", y = 'Odds of persistence', title = NULL) +
  mytheme + theme(legend.position = "none") + 
  scale_fill_manual(values=c("Shandong" = "#E69F00", "Henan" = "#57B4E9", "Hubei" = "#019E73",
                             "Hunan" = "#F0E442","Guangxi" = "#0072B2","Guangdong" = "#D55E00"))+
  scale_color_manual(values=c("Shandong" = "#E69F00", "Henan" = "#57B4E9", "Hubei" = "#019E73",
                              "Hunan" = "#F0E442","Guangxi" = "#0072B2","Guangdong" = "#D55E00"))+
  scale_x_continuous(labels = scales::label_comma(accuracy =0.1)) +
  scale_y_continuous(labels = scales::label_comma(accuracy =0.01))-> p3; p3

################################################################################
unique(field_com_data$Seed_source)
field_com_data$Seed_source = as.factor(field_com_data$Seed_source)
field_com_data_no = field_com_data %>% tidyr::drop_na(rebio2022_100)
unique(field_com_data_no$Species)

################################################################################
sp = unique(field_com_data_no$Species)
final_selected_row = NULL
for (i in sp) {
  data_sub = subset(field_com_data_no, Species == i)
  max_index <- which.max(data_sub$rebio2022)
  selected_row <- data_sub[max_index, ]
  final_selected_row = rbind(final_selected_row, selected_row)
}
head(final_selected_row)

first_char <- substr(final_selected_row$Species, 1, 1)
sub_str <- gsub(".*_", "", final_selected_row$Species)
Latin_name <- paste(first_char, sub_str, sep = ". ")
final_selected_row$Latin_name = Latin_name
df2 = (final_selected_row %>% arrange(desc(rebio2022_100)))[c(1:8),]


################################################################################
#### SLA
mod12 <- lmer(rebio2022_100 ~ Field_SLA_sqrt*Seed_source + (1|Block), data = field_com_data_no)
summary(mod12)
anova(mod12)
MuMIn::r.squaredGLMM(mod12)
field_com_data_no$F0 = predictSE(mod12, field_com_data_no, level = 0)$fit
field_com_data_no$SE <- predictSE(mod12, field_com_data_no, level = 0)$se.fit

####
field_com_data_all = field_com_data_no
mod12 <- lmer(rebio2022_100 ~ Field_SLA_sqrt + (1|Block/Plot_num), data = field_com_data_all)
Anova(mod12)
field_com_data_all$F0 = predictSE(mod12, field_com_data_all, level = 0)$fit
field_com_data_all$SE <- predictSE(mod12, field_com_data_all, level = 0)$se.fit
summary(field_com_data_all$rebio2022_100)
summary(field_com_data_all$Field_SLA_sqrt)

ggplot(field_com_data_no, aes(x=Field_SLA_sqrt, y=rebio2022_100)) +
  #geom_ribbon(aes(x=Field_SLA_sqrt,ymin = F0 - 1.96 * SE,ymax = F0 + 1.96 * SE, fill = Seed_source), alpha = I(0.1)) +
  geom_line(data = subset(field_com_data_no, Seed_source != "Hubei" & Seed_source != "Hunan"), mapping = aes(x=Field_SLA_sqrt, y=F0, color = Seed_source), size=1, linetype = 1) +
  geom_line(data = subset(field_com_data_no, Seed_source == "Hubei"), mapping = aes(x=Field_SLA_sqrt, y=F0, color = Seed_source), size=1, linetype = 2) +
  geom_line(data = subset(field_com_data_no, Seed_source == "Hunan"), mapping = aes(x=Field_SLA_sqrt, y=F0, color = Seed_source), size=1, linetype = 2) +
  geom_line(data = field_com_data_all, mapping = aes(x=Field_SLA_sqrt, y=F0), size=1.5, linetype = 2, color = "black") +
  geom_ribbon(data = field_com_data_all, mapping = aes(x=Field_SLA_sqrt, ymin = F0 - 1.96 * SE,ymax = F0 + 1.96 * SE, fill = "#EBEBEB"), alpha = I(0.1)) +
  geom_point(aes(fill = Seed_source),size = 2.2, color = "black", pch = 21) + # fill = "#00000022", color = "#0A0A0A"
  scale_fill_manual(values=c("Shandong" = "#E69F00", "Henan" = "#57B4E9", "Hubei" = "#019E73",
                             "Hunan" = "#F0E442","Guangxi" = "#0072B2","Guangdong" = "#D55E00"))+
  scale_color_manual(values=c("Shandong" = "#E69F00", "Henan" = "#57B4E9", "Hubei" = "#019E73",
                              "Hunan" = "#F0E442","Guangxi" = "#0072B2","Guangdong" = "#D55E00"))+
  ggrepel::geom_text_repel(mapping = aes(x=Field_SLA_sqrt,y=rebio2022_100,label=Latin_name), data = df2,size = 2.8,segment.color = "black", color = "black",direction = "both",box.padding = 0.6,
                           max.overlaps = getOption("ggrepel.max.overlaps", default = 25), fontface = "italic") +
  scale_x_continuous(labels = scales::label_comma(accuracy =0.1)) +
  scale_y_continuous(labels = scales::label_comma(accuracy =0.01)) +
  labs(x = "Specific leaf area (cm2/g, log10)\nestimated in field experiment",y = 'Relative abundance within the \n plot in second year (%, log10)',title = NULL) +
  mytheme -> p4; p4 

(p1+p2)/(p3+p4)
################################################################################
### Fig S4
p1+p2+p3+p4+
  plot_layout(ncol = 2, nrow = 2) +
  plot_annotation(tag_levels = 'a')

### Notice that,
### For more picture details, we have further adjusted it in Adobe illustrator.
