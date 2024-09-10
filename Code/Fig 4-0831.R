##################################################################################
###################################   Fig 4   ####################################
##################################################################################
### Relationship between functional distinctiveness (muti-traits), persistence 
### and relative abundance within the plot in second year
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
phylogenyAux = read.tree("Data/iq_tree.treefile")
plot(phylogenyAux)
to_drop = phylogenyAux$tip.label[!phylogenyAux$tip.label %in% Common_sp_list$Species]
tree <- drop.tip(as.phylo(phylogenyAux), to_drop) 
phyloMat = vcv.phylo(tree)
phyloMat = phyloMat / max(phyloMat)
dim(phyloMat)

### Loading pot experiment database (mean value of traits)
pot_trait = read.xlsx("Data/Pot_traits_mean0831.xlsx", sheet = "Pot_means", colNames = TRUE, rowNames = FALSE)
colnames(pot_trait) <- paste0(colnames(pot_trait), "_pot")
colnames(pot_trait)[1] <- "Species"
nrow(pot_trait)

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
All_total_dis2 = All_total_dis[ ,c("Species","All_pot_mean","All_Field_means")]
names(All_total_dis2)

library(openxlsx)
library(lmerTest)
library(lme4)
library(MuMIn)
library(car)
library(AICcmodavg)
library(ggplot2)
field_com_data = read.xlsx("Data/all_row_data0829.xlsx", sheet = "field_data_mean", rowNames = F, colNames = T)
head(field_com_data)

### Add uniqueness parameter
field_com_data = field_com_data %>% left_join(All_total_dis2, by = "Species")
colnames(field_com_data)
field_com_data$rebio2022_100 = log10(field_com_data$rebio2022*100)
field_com_data$Block = as.factor(field_com_data$Block)
length(unique(field_com_data$Species))
field_com_data$Seed_source <- factor(field_com_data$Seed_source, levels = rev(c("Guangdong","Guangxi","Hunan","Hubei","Henan","Shandong")))
field_com_data$exist_prob <- ifelse(!is.na(field_com_data$rebio2022) & field_com_data$rebio2022 > 0, 1, 0)

##### Based on field measurement trait data
##### odds of persistence
Seed_source = unique(field_com_data$Seed_source)
all_summary_glmer = NULL
for (i in Seed_source) {
  select_data = subset(field_com_data, Seed_source == i)
  mod1 <- glmer(exist_prob ~ All_pot_mean + (1|Block) , na.action=na.omit, family=binomial, data=select_data)
  mod1_summary = data.frame(Pop = i, Traits = "All_pot_mean", N = nrow(select_data), sp_num = length(unique(select_data$Species)),
                            Chisq = as.data.frame(Anova(mod1))[,1], p_value = as.data.frame(Anova(mod1))[,3],
                            R2m = r.squaredGLMM(mod1)[1,1], R2c = r.squaredGLMM(mod1)[1,2])
  mod2 <- glmer(exist_prob ~ All_Field_means + (1|Block) , na.action=na.omit, family=binomial, data=select_data)
  mod2_summary = data.frame(Pop = i, Traits = "All_Field_means", N = nrow(select_data), sp_num = length(unique(select_data$Species)),
                            Chisq = as.data.frame(Anova(mod2))[,1], p_value = as.data.frame(Anova(mod2))[,3],
                            R2m = r.squaredGLMM(mod2)[1,1], R2c = r.squaredGLMM(mod2)[1,2])
  Total_summary = rbind(mod1_summary,mod2_summary)
  all_summary_glmer = rbind(all_summary_glmer, Total_summary)
}
all_summary_glmer$R2m = round(all_summary_glmer$R2m, 3)

##### Overall field data analysis
colnames(field_com_data)
mod1 = glmer(exist_prob ~ All_pot_mean + (1|Block/Plot_num) , na.action=na.omit, family=binomial, data=field_com_data)
Anova(mod1); r.squaredGLMM(mod1)[1,]
summary(mod1)

mod2 = glmer(exist_prob ~ All_Field_means + (1|Block/Plot_num) , na.action=na.omit, family=binomial, data=field_com_data)
Anova(mod2); r.squaredGLMM(mod2)[1,]
summary(mod2)

nrow(field_com_data); length(unique(field_com_data$Species))

##### relative abundance of species in second year
Seed_source = unique(field_com_data$Seed_source)
all_summary_lmer = NULL
for (i in Seed_source) {
  select_data = subset(field_com_data, Seed_source == i) %>% tidyr::drop_na(rebio2022_100)
  mod1 <- lmer(rebio2022_100 ~ All_pot_mean + (1|Block) , na.action=na.omit, data=select_data)
  mod1_summary = data.frame(Pop = i, Traits = "All_pot_mean", N = nrow(select_data), sp_num = length(unique(select_data$Species)),
                            Chisq = as.data.frame(Anova(mod1))[,1], p_value = as.data.frame(Anova(mod1))[,3],
                            R2m = r.squaredGLMM(mod1)[1,1], R2c = r.squaredGLMM(mod1)[1,2])
  mod2 <- lmer(rebio2022_100 ~ All_Field_means + (1|Block) , na.action=na.omit, data=select_data)
  mod2_summary = data.frame(Pop = i, Traits = "All_Field_means", N = nrow(select_data), sp_num = length(unique(select_data$Species)),
                            Chisq = as.data.frame(Anova(mod2))[,1], p_value = as.data.frame(Anova(mod2))[,3],
                            R2m = r.squaredGLMM(mod2)[1,1], R2c = r.squaredGLMM(mod2)[1,2])
  Total_summary = rbind(mod1_summary,mod2_summary)
  all_summary_lmer = rbind(all_summary_lmer, Total_summary)
}
all_summary_lmer$p_value = round(all_summary_lmer$p_value, 3)

select_data = subset(field_com_data, Seed_source == "Guangdong") %>% tidyr::drop_na(rebio2022_100)
summary(select_data$All_Field_means)
summary(select_data$rebio2022_100)
mod1 <- lmer(rebio2022_100 ~ All_pot_mean + (1|Block) , na.action=na.omit, data=select_data)
Anova(mod1)

##### Overall data analysis
colnames(field_com_data)
select_data = field_com_data %>% tidyr::drop_na(rebio2022_100)
nrow(select_data); length(unique(select_data$Species))

mod1 = lmer(rebio2022_100 ~ All_pot_mean + (1|Block/Plot_num) , na.action=na.omit, data=field_com_data)
mod1 = lmer(rebio2022_100 ~ All_pot_mean + (1|Block/Plot_num) , data=select_data)
Anova(mod1); r.squaredGLMM(mod1)[1,]

mod2 = lmer(rebio2022_100 ~ All_Field_means + (1|Block/Plot_num) , na.action=na.omit, data=field_com_data)
Anova(mod2); r.squaredGLMM(mod2)[1,]


################################################################################
######### Visualization
length(unique(field_com_data$Species))
unique(field_com_data$Seed_source)

#### All_pot_mean
mod12 <- glmer(exist_prob ~ All_pot_mean*Seed_source + (1|Block) , na.action=na.omit, family=binomial, data=field_com_data)
Anova(mod12)
summary(mod12)
MuMIn::r.squaredGLMM(mod12)

field_com_data$F0 = predictSE(mod12, field_com_data, level = 0)$fit
field_com_data$SE <- predictSE(mod12, field_com_data, level = 0)$se.fit

####
field_com_data_all = field_com_data
mod12 <- glmer(exist_prob ~ All_pot_mean + (1|Block/Plot_num) , na.action=na.omit, family=binomial, data=field_com_data_all)
Anova(mod12)
field_com_data_all$F0 = predictSE(mod12, field_com_data_all, level = 0)$fit
field_com_data_all$SE <- predictSE(mod12, field_com_data_all, level = 0)$se.fit

ggplot(field_com_data, aes(x=All_pot_mean, y=exist_prob)) +
  #geom_ribbon(aes(x=All_pot_mean,ymin = F0 - 1.96 * SE,ymax = F0 + 1.96 * SE, fill = Seed_source), alpha = I(0.1)) +
  geom_line(data = subset(field_com_data, Seed_source != "Hubei" & Seed_source != "Shandong"), mapping = aes(x=All_pot_mean, y=F0, color = Seed_source), size=1, linetype = 2) +
  geom_line(data = subset(field_com_data, Seed_source == "Hubei"), mapping = aes(x=All_pot_mean, y=F0, color = Seed_source), size=1, linetype = 1) +
  geom_line(data = subset(field_com_data, Seed_source == "Shandong"), mapping = aes(x=All_pot_mean, y=F0, color = Seed_source), size=1, linetype = 1) +
  geom_line(data = field_com_data_all, mapping = aes(x=All_pot_mean, y=F0), size=1.5, linetype = 1, color = "black") +
  geom_ribbon(data = field_com_data_all, mapping = aes(x=All_pot_mean,ymin = F0 - 1.96 * SE,ymax = F0 + 1.96 * SE, fill = "#EBEBEB"), alpha = I(0.1)) +
  geom_point(aes(fill = Seed_source),size = 2.2, pch = 21) + 
  ylab('Odds of persistence') + 
  xlab('Functional distinctiveness\nestimated in pot experiment') +
  mytheme +
  scale_fill_manual(values=c("Shandong" = "#376694", "Henan" = "#73BCD5", "Hubei" = "#ABDCE0",
                             "Hunan" = "#EFBA55","Guangxi" = "#E68C51","Guangdong" = "#994240"))+
  scale_color_manual(values=c("Shandong" = "#376694", "Henan" = "#73BCD5", "Hubei" = "#ABDCE0",
                              "Hunan" = "#EFBA55","Guangxi" = "#E68C51","Guangdong" = "#994240"))+
  scale_x_continuous(labels = scales::label_comma(accuracy =0.1)) +
  scale_y_continuous(labels = scales::label_comma(accuracy =0.01))-> p1; p1

################################################################################
#### All_Field_means
mod12 <- glmer(exist_prob ~ All_Field_means*Seed_source + (1|Block) , na.action=na.omit, family=binomial, data=field_com_data)
Anova(mod12)
summary(mod12)
MuMIn::r.squaredGLMM(mod12)

field_com_data$F0 = predictSE(mod12, field_com_data, level = 0)$fit
field_com_data$SE <- predictSE(mod12, field_com_data, level = 0)$se.fit

####
field_com_data_all = field_com_data
mod12 <- glmer(exist_prob ~ All_Field_means + (1|Block/Plot_num) , na.action=na.omit, family=binomial, data=field_com_data_all)
Anova(mod12)
field_com_data_all$F0 = predictSE(mod12, field_com_data_all, level = 0)$fit
field_com_data_all$SE <- predictSE(mod12, field_com_data_all, level = 0)$se.fit

ggplot(field_com_data, aes(x=(All_Field_means), y=exist_prob)) +
  #geom_ribbon(aes(x=All_Field_means,ymin = F0 - 1.96 * SE,ymax = F0 + 1.96 * SE, fill = Seed_source), alpha = I(0.1)) +
  geom_line(data = subset(field_com_data, Seed_source != "Guangdong" & Seed_source != "Hubei" & Seed_source != "Hunan"), mapping = aes(x=All_Field_means, y=F0, color = Seed_source), size=1, linetype = 2) +
  geom_line(data = subset(field_com_data, Seed_source == "Guangdong"), mapping = aes(x=All_Field_means, y=F0, color = Seed_source), size=1, linetype = 1) +
  geom_line(data = subset(field_com_data, Seed_source == "Hubei"), mapping = aes(x=All_Field_means, y=F0, color = Seed_source), size=1, linetype = 1) +
  geom_line(data = subset(field_com_data, Seed_source == "Hunan"), mapping = aes(x=All_Field_means, y=F0, color = Seed_source), size=1, linetype = 1) +
  geom_line(data = field_com_data_all, mapping = aes(x=All_Field_means, y=F0), size=1.5, linetype = 1, color = "black") +
  geom_ribbon(data = field_com_data_all, mapping = aes(x=All_Field_means,ymin = F0 - 1.96 * SE,ymax = F0 + 1.96 * SE, fill = "#EBEBEB"), alpha = I(0.1)) +
  geom_point(aes(fill = Seed_source),size = 2.2, pch = 21) + # fill = "#00000022", color = "#0A0A0A"
  ylab(' \n ') + 
  xlab('Functional distinctiveness\nestimated in field experiment') +
  mytheme + theme(legend.position = "none") + 
  #scale_fill_manual(values = c('#60A7A6','#FEA6A6'))+
  scale_fill_manual(values=c("Shandong" = "#376694", "Henan" = "#73BCD5", "Hubei" = "#ABDCE0",
                             "Hunan" = "#EFBA55","Guangxi" = "#E68C51","Guangdong" = "#994240"))+
  scale_color_manual(values=c("Shandong" = "#376694", "Henan" = "#73BCD5", "Hubei" = "#ABDCE0",
                              "Hunan" = "#EFBA55","Guangxi" = "#E68C51","Guangdong" = "#994240"))+
  scale_x_continuous(labels = scales::label_comma(accuracy =0.1)) +
  scale_y_continuous(labels = scales::label_comma(accuracy =0.01)) ->p2; p2

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
#### All_pot_mean
mod12 <- lmer(rebio2022_100 ~ All_pot_mean*Seed_source + (1|Block), data = field_com_data_no)
summary(mod12)
anova(mod12)
MuMIn::r.squaredGLMM(mod12)
field_com_data_no$F0 = predictSE(mod12, field_com_data_no, level = 0)$fit
field_com_data_no$SE <- predictSE(mod12, field_com_data_no, level = 0)$se.fit

####
field_com_data_all = field_com_data_no
mod12 <- lmer(rebio2022_100 ~ All_pot_mean + (1|Block/Plot_num), data = field_com_data_all)
Anova(mod12)
field_com_data_all$F0 = predictSE(mod12, field_com_data_all, level = 0)$fit
field_com_data_all$SE <- predictSE(mod12, field_com_data_all, level = 0)$se.fit

ggplot(field_com_data_no, aes(x=All_pot_mean, y=rebio2022_100)) +
  #geom_ribbon(aes(x=All_pot_mean,ymin = F0 - 1.96 * SE,ymax = F0 + 1.96 * SE, fill = Seed_source), alpha = I(0.1)) +
  geom_line(data = subset(field_com_data_no, Seed_source != "Guangxi" & Seed_source != "Henan" & Seed_source != "Shandong"), mapping = aes(x=All_pot_mean, y=F0, color = Seed_source), size=1, linetype = 1) +
  geom_line(data = subset(field_com_data_no, Seed_source == "Guangxi"), mapping = aes(x=All_pot_mean, y=F0, color = Seed_source), size=1, linetype = 2) +
  geom_line(data = subset(field_com_data_no, Seed_source == "Henan"), mapping = aes(x=All_pot_mean, y=F0, color = Seed_source), size=1, linetype = 2) +
  geom_line(data = subset(field_com_data_no, Seed_source == "Shandong"), mapping = aes(x=All_pot_mean, y=F0, color = Seed_source), size=1, linetype = 2) +
  geom_line(data = field_com_data_all, mapping = aes(x=All_pot_mean, y=F0), size=1.5, linetype = 1, color = "black") +
  geom_ribbon(data = field_com_data_all, mapping = aes(x=All_pot_mean,ymin = F0 - 1.96 * SE,ymax = F0 + 1.96 * SE, fill = "#EBEBEB"), alpha = I(0.1)) +
  geom_point(aes(fill = Seed_source),size = 2.2, pch = 21) + # fill = "#00000022", color = "#0A0A0A"
  scale_fill_manual(values=c("Shandong" = "#376694", "Henan" = "#73BCD5", "Hubei" = "#ABDCE0",
                             "Hunan" = "#EFBA55","Guangxi" = "#E68C51","Guangdong" = "#994240"))+
  scale_color_manual(values=c("Shandong" = "#376694", "Henan" = "#73BCD5", "Hubei" = "#ABDCE0",
                              "Hunan" = "#EFBA55","Guangxi" = "#E68C51","Guangdong" = "#994240"))+
  #ggrepel::geom_text_repel(mapping = aes(x=All_pot_mean,y=rebio2022_100,label=Latin_name), data = df2,size = 2.8,segment.color = "black", color = "black",direction = "both",box.padding = 0.6,
  #                         max.overlaps = getOption("ggrepel.max.overlaps", default = 25), fontface = "italic") +
  scale_x_continuous(labels = scales::label_comma(accuracy =0.1)) +
  scale_y_continuous(labels = scales::label_comma(accuracy =0.01)) +
  ylab('Relative abundance within the\nplot in second year (%, log10)') + 
  xlab('Functional distinctiveness\nestimated in pot experiment') +
  mytheme -> p3; p3 


#### All_Field_means
mod12 <- lmer(rebio2022_100 ~ All_Field_means*Seed_source + (1|Block), data = field_com_data_no)
summary(mod12)
anova(mod12)
MuMIn::r.squaredGLMM(mod12)
field_com_data_no$F0 = predictSE(mod12, field_com_data_no, level = 0)$fit
field_com_data_no$SE <- predictSE(mod12, field_com_data_no, level = 0)$se.fit

####
field_com_data_all = field_com_data_no
mod12 <- lmer(rebio2022_100 ~ All_Field_means + (1|Block/Plot_num), data = field_com_data_all)
Anova(mod12)
field_com_data_all$F0 = predictSE(mod12, field_com_data_all, level = 0)$fit
field_com_data_all$SE <- predictSE(mod12, field_com_data_all, level = 0)$se.fit

ggplot(field_com_data_no, aes(x=All_Field_means, y=rebio2022_100)) +
  #geom_ribbon(aes(x=All_Field_means,ymin = F0 - 1.96 * SE,ymax = F0 + 1.96 * SE, fill = Seed_source), alpha = I(0.1)) +
  geom_line(data = subset(field_com_data_no, Seed_source != "Guangdong" & Seed_source != "Henan"), mapping = aes(x=All_Field_means, y=F0, color = Seed_source), size=1, linetype = 1) +
  geom_line(data = subset(field_com_data_no, Seed_source == "Guangdong"), mapping = aes(x=All_Field_means, y=F0, color = Seed_source), size=1, linetype = 2) +
  geom_line(data = subset(field_com_data_no, Seed_source == "Henan"), mapping = aes(x=All_Field_means, y=F0, color = Seed_source), size=1, linetype = 2) +
  geom_line(data = field_com_data_all, mapping = aes(x=All_Field_means, y=F0), size=1.5, linetype = 1, color = "black") +
  geom_ribbon(data = field_com_data_all, mapping = aes(x=All_Field_means,ymin = F0 - 1.96 * SE,ymax = F0 + 1.96 * SE, fill = "#EBEBEB"), alpha = I(0.1)) +
  geom_point(aes(fill = Seed_source),size = 2.2, pch = 21) + # fill = "#00000022", color = "#0A0A0A"
  scale_fill_manual(values=c("Shandong" = "#376694", "Henan" = "#73BCD5", "Hubei" = "#ABDCE0",
                             "Hunan" = "#EFBA55","Guangxi" = "#E68C51","Guangdong" = "#994240"))+
  scale_color_manual(values=c("Shandong" = "#376694", "Henan" = "#73BCD5", "Hubei" = "#ABDCE0",
                              "Hunan" = "#EFBA55","Guangxi" = "#E68C51","Guangdong" = "#994240"))+
  #ggrepel::geom_text_repel(mapping = aes(x=All_Field_means,y=rebio2022_100,label=Latin_name), data = df2,size = 2.8,segment.color = "black", color = "black",direction = "both",box.padding = 0.6,
  #                         max.overlaps = getOption("ggrepel.max.overlaps", default = 25), fontface = "italic") +
  scale_x_continuous(labels = scales::label_comma(accuracy =0.1)) +
  scale_y_continuous(labels = scales::label_comma(accuracy =0.01)) +
  ylab(' \n ') + 
  xlab('Functional distinctiveness\nestimated in field experiment') +
  mytheme -> p4; p4


#####
(p1|p2)/(p3|p4)

