############packages##########
library(openxlsx)
library(tidyverse)
library(multcomp)
library(emmeans)
library(tidyr)
library(pheatmap)
load("heatmap.Rdata")###load plot.Rdata

###
Common_sp_list = read.xlsx("Data/Common_species_list.xlsx", sheet = "Common_sp_list", colNames = TRUE, rowNames = FALSE)
Common_sp_list_SLA = c(na.omit(Common_sp_list$Species_SLA))
Common_sp_list_AGB = c(na.omit(Common_sp_list$Species_AGB))
###

trait_data_total <- read.xlsx("Data/all_row_data0829.xlsx", sheet = "Field_data_row_mean", rowNames = F, colNames = T)
trait_data_total <- trait_data_total[trait_data_total$Species %in% Common_sp_list_AGB, ]
unique(trait_data_total$Species)

####SLA
trait_data <- subset(trait_data_total, SLA != "NA")
unique(trait_data$Species)
####
summary_sp_data <- (trait_data[,c("Species","Seed_source")])
summary_sp_data$count = 1
summary_sp_data <- Rmisc::summarySE(summary_sp_data, measurevar = "count", groupvars = c("Species","Seed_source"))
summary_sp_data <- summary_sp_data[summary_sp_data$N >2, ]
unique(summary_sp_data$Species)

summary_sp_data <- summary_sp_data %>% group_by(Species) %>% summarise(sum_count = sum(count))
select_sp_SLA <- unique(summary_sp_data[summary_sp_data$sum_count > 1, ]$Species)

final_result_mod = NULL
for (i in select_sp_SLA) {
  select_data <- trait_data[trait_data$Species %in% i,]
  colnames(select_data)
  mod1 <- aov(sqrt(SLA) ~ Seed_source, data = select_data)
  summary_mod1 <- summary(mod1)
  result_mod1 <- data.frame(Species = i, Traits = "SLA", F_value = summary_mod1[[1]]$`F value`[1], p_value = summary_mod1[[1]]$"Pr(>F)"[1])
  final_result_mod = rbind(final_result_mod,result_mod1)
}

final_result_mod$p_value = round(final_result_mod$p_value, 3)
final_result_mod

######Species with remarkable differences
sig_SLA_spe <- subset(final_result_mod, p_value < 0.05)$Species
cld_result_mod = NULL
for (i in sig_SLA_spe) {
  select_data <- trait_data[trait_data$Species %in% i,]
  colnames(select_data)
  mod1 <- aov(sqrt(SLA) ~ Seed_source, data = select_data)
  model_means <- emmeans(mod1, specs = ~ Seed_source)
  model_means_cld <- cld(object = model_means,
                         adjust = "none",
                         Letters = letters,
                         alpha = 0.05)
  result_cld1 = data.frame(Species = i , Seed_source = model_means_cld$Seed_source, group = model_means_cld$.group)
  cld_result_mod = rbind(cld_result_mod,result_cld1)
}
cld_result_mod
unique(cld_result_mod$Species)

###Hmax
trait_data <- subset(trait_data_total, Hmax != "NA")
unique(trait_data$Species)
####
summary_sp_data <- (trait_data[,c("Species","Seed_source")])
summary_sp_data$count = 1
summary_sp_data <- Rmisc::summarySE(summary_sp_data, measurevar = "count", groupvars = c("Species","Seed_source"))
summary_sp_data <- summary_sp_data[summary_sp_data$N >2, ]
unique(summary_sp_data$Species)
unique(summary_sp_data$Species)
summary_sp_data <- summary_sp_data %>% group_by(Species) %>% summarise(sum_count = sum(count))
select_sp_Hmax <- unique(summary_sp_data[summary_sp_data$sum_count >1, ]$Species)

final_result_mod = NULL
for (i in select_sp_Hmax) {
  select_data <- trait_data[trait_data$Species %in% i,]
  colnames(select_data)
  mod2 <- aov(sqrt(Hmax) ~ Seed_source, data = select_data)
  summary_mod2 <- summary(mod2)
  result_mod2 <- data.frame(Species = i, Traits = "Hmax", F_value = summary_mod2[[1]]$`F value`[1], p_value = summary_mod2[[1]]$"Pr(>F)"[1])
  final_result_mod = rbind(final_result_mod,result_mod2)
}

final_result_mod$p_value = round(final_result_mod$p_value, 3)
final_result_mod

######Species with remarkable differences
sig_Hmax_spe <- subset(final_result_mod, p_value < 0.05)$Species
cld_result_mod = NULL
for (i in sig_Hmax_spe) {
  select_data <- trait_data[trait_data$Species %in% i,]
  colnames(select_data)
  mod1 <- aov(sqrt(Hmax) ~ Seed_source, data = select_data)
  model_means <- emmeans(mod1, specs = ~ Seed_source)
  model_means_cld <- cld(object = model_means,
                         adjust = "none",
                         Letters = letters,
                         alpha = 0.05)
  result_cld1 = data.frame(Species = i , Seed_source = model_means_cld$Seed_source, group = model_means_cld$.group)
  cld_result_mod = rbind(cld_result_mod,result_cld1)
}
cld_result_mod
unique(cld_result_mod$Species)

###AGB
trait_data <- subset(trait_data_total, AGB != "NA")
unique(trait_data$Species)
####
summary_sp_data <- (trait_data[,c("Species","Seed_source")])
summary_sp_data$count = 1
summary_sp_data <- Rmisc::summarySE(summary_sp_data, measurevar = "count", groupvars = c("Species","Seed_source"))
summary_sp_data <- summary_sp_data[summary_sp_data$N >2, ]
unique(summary_sp_data$Species)

summary_sp_data <- summary_sp_data %>% group_by(Species) %>% summarise(sum_count = sum(count))
select_sp_AGB <- unique(summary_sp_data[summary_sp_data$sum_count >1, ]$Species)

###
final_result_mod = NULL
for (i in select_sp_AGB) {
  select_data <- trait_data[trait_data$Species %in% i,]
  colnames(select_data)
  mod3 <- aov(log10(AGB) ~ Seed_source, data = select_data)
  summary_mod3 <- summary(mod3)
  result_mod3 <- data.frame(Species = i, Traits = "AGB",F_value = summary_mod3[[1]]$`F value`[1], p_value = summary_mod3[[1]]$"Pr(>F)"[1])
  result_mod <- rbind(result_mod2,result_mod3)
  ###
  final_result_mod = rbind(final_result_mod,result_mod3)
}

final_result_mod$p_value = round(final_result_mod$p_value, 3)
final_result_mod

######Species with remarkable differences
sig_AGB_spe <- subset(final_result_mod, p_value < 0.05)$Species
cld_result_mod = NULL
for (i in sig_AGB_spe) {
  select_data <- trait_data[trait_data$Species %in% i,]
  colnames(select_data)
  mod1 <- aov(log10(AGB) ~ Seed_source, data = select_data)
  model_means <- emmeans(mod1, specs = ~ Seed_source)
  model_means_cld <- cld(object = model_means,
                         adjust = "none",
                         Letters = letters,
                         alpha = 0.05)
  result_cld1 = data.frame(Species = i , Seed_source = model_means_cld$Seed_source, group = model_means_cld$.group)
  cld_result_mod = rbind(cld_result_mod,result_cld1)
}
cld_result_mod

###Integrated species data
sp_select_SLA_Hmax = union(select_sp_SLA,select_sp_Hmax)
sp_select = union(sp_select_SLA_Hmax,select_sp_AGB)

spe_select = (trait_data[trait_data$Species %in% sp_select,])
spe_select$SLA <- sqrt(spe_select$SLA)
spe_select$Hmax <- sqrt(spe_select$Hmax)
spe_select$AGB <- log10(spe_select$AGB)
spe_select_mean <- aggregate(Hmax ~ Seed_source + Species, data = spe_select, FUN = mean) 

spe_select_mean$Seed_source <- factor(spe_select_mean$Seed_source, levels = c("Shandong","Henan","Hubei","Hunan","Guangxi","Guangdong"))
spe_select_mean_trans <- spread(spe_select_mean, key = "Seed_source",
                                value = "Hmax")
spe_select_mean_trans[is.na(spe_select_mean_trans)] = NA
rownames(spe_select_mean_trans) <- spe_select_mean_trans[, 1]
spe_select_mean_trans <- spe_select_mean_trans[-1]

spe_select_mean_trans$latin = rownames(spe_select_mean_trans)
spe_select_mean_trans$latin <- gsub("_", " ", spe_select_mean_trans$latin)
rownames(spe_select_mean_trans) = spe_select_mean_trans$latin
spe_select_mean_trans$latin = NULL

pheatmap(spe_select_mean_trans,
         annotation_row = NULL,
         annotation_names_row = TRUE,
         labels_row = rownames(spe_select_mean_trans),
         legend = TRUE,
         show_rownames = TRUE, 
         angle_col = 90,
         fontsize_number = textSize,
         border_color = "#C1C1C1",
         na_col = "white",
         fontsize = legendTextSize,
         treeheight_row = 0, 
         treeheight_col = 0,
         legend_labels = "sig*log(p-value)",
         color = colorRampPalette(colors = c("#F0F6E4","#E29C35"))(40),###You can change the heatmap by changing colors 
         fontsize_col = colTextSize,
         fontsize_row = rowTextSize,
         cellwidth = cellWidth,
         cellheight = cellHeight,
         filename = NA,
         cluster_cols = FALSE,
         cluster_rows = FALSE)
##color we selected: SLA:#D9352A    Hmax:#E29C35   AGB:#414887

