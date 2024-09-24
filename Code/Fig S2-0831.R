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

trait_data_total <- read.xlsx("Data/all_row_data.xlsx", sheet = "field_data_mean2", rowNames = F, colNames = T)
trait_data_total <- trait_data_total[trait_data_total$Species %in% Common_sp_list_AGB, ]
unique(trait_data_total$Species)

###Hmax
trait_data <- subset(trait_data_total, Hmax_obs != "NA")
unique(trait_data$Species)
summary(trait_data$Hmax); dim(trait_data)

####
summary_sp_data <- (trait_data[,c("Species","Seed_source")])
summary_sp_data$count = 1
summary_sp_data <- Rmisc::summarySE(summary_sp_data, measurevar = "count", groupvars = c("Species","Seed_source"))
summary_sp_data <- summary_sp_data[summary_sp_data$N >2, ]
unique(summary_sp_data$Species)

trait_data = trait_data %>% left_join(summary_sp_data[,c("Species","Seed_source","N")], by = c("Species","Seed_source"))
trait_data = trait_data %>% tidyr::drop_na(N)

summary_sp_data <- summary_sp_data %>% group_by(Species) %>% summarise(sum_count = sum(count))
select_sp_Hmax <- unique(summary_sp_data[summary_sp_data$sum_count >1, ]$Species)

final_result_mod = NULL
for (i in select_sp_Hmax) {
  select_data <- trait_data[trait_data$Species %in% i,]
  colnames(select_data)
  mod2 <- aov(sqrt(Hmax_obs) ~ Seed_source, data = select_data)
  summary_mod2 <- summary(mod2)
  result_mod2 <- data.frame(Species = i, Traits = "Hmax", F_value = summary_mod2[[1]]$`F value`[1], p_value = summary_mod2[[1]]$"Pr(>F)"[1])
  final_result_mod = rbind(final_result_mod,result_mod2)
}

final_result_mod$p_value = round(final_result_mod$p_value, 3)
View(subset(final_result_mod, p_value < 0.05))
dim(subset(final_result_mod, p_value < 0.05))

######Species with remarkable differences
sig_Hmax_spe <- subset(final_result_mod, p_value < 0.05)$Species
cld_result_mod = NULL
for (i in sig_Hmax_spe) {
  select_data <- trait_data[trait_data$Species %in% i,]
  colnames(select_data)
  mod1 <- aov(sqrt(Hmax_obs) ~ Seed_source, data = select_data)
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


### plot
spe_select = (trait_data[trait_data$Species %in% select_sp_Hmax,])
spe_select$Hmax <- sqrt(spe_select$Hmax_obs)
spe_select_mean <- aggregate(Hmax ~ Seed_source + Species, data = spe_select, FUN = mean)

spe_select_mean$Seed_source <- factor(spe_select_mean$Seed_source, levels = c("Shandong","Henan","Hubei","Hunan","Guangxi","Guangdong"))
spe_select_mean_trans <- spread(spe_select_mean, key = "Seed_source",value = "Hmax")

spe_select_mean_trans$latin <- gsub("_", " ", spe_select_mean_trans$Species)
rownames(spe_select_mean_trans) = spe_select_mean_trans$latin
spe_select_mean_trans$latin = NULL; spe_select_mean_trans$Species = NULL
spe_select_mean_trans$p_value = NA
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


###AGB
trait_data <- subset(trait_data_total, AGB_obs != "NA")
unique(trait_data$Species); dim(trait_data)
####
summary_sp_data <- (trait_data[,c("Species","Seed_source")])
summary_sp_data$count = 1
summary_sp_data <- Rmisc::summarySE(summary_sp_data, measurevar = "count", groupvars = c("Species","Seed_source"))
summary_sp_data <- summary_sp_data[summary_sp_data$N >2, ]
unique(summary_sp_data$Species)

trait_data = trait_data %>% left_join(summary_sp_data[,c("Species","Seed_source","N")], by = c("Species","Seed_source"))
trait_data = trait_data %>% tidyr::drop_na(N)
unique(trait_data$Species)

summary_sp_data <- summary_sp_data %>% group_by(Species) %>% summarise(sum_count = sum(count))
select_sp_AGB <- unique(summary_sp_data[summary_sp_data$sum_count >1, ]$Species)

###
final_result_mod = NULL
for (i in select_sp_AGB) {
  select_data <- trait_data[trait_data$Species %in% i,]
  colnames(select_data)
  mod3 <- aov(log10(AGB_obs) ~ Seed_source, data = select_data)
  summary_mod3 <- summary(mod3)
  result_mod3 <- data.frame(Species = i, Traits = "AGB",F_value = summary_mod3[[1]]$`F value`[1], p_value = summary_mod3[[1]]$"Pr(>F)"[1])
  result_mod <- rbind(result_mod2,result_mod3)
  ###
  final_result_mod = rbind(final_result_mod,result_mod3)
}

final_result_mod$p_value = round(final_result_mod$p_value, 3)
View(subset(final_result_mod, p_value < 0.05))
dim(subset(final_result_mod, p_value < 0.05))


######Species with remarkable differences
sig_AGB_spe <- subset(final_result_mod, p_value < 0.05)$Species
cld_result_mod = NULL
for (i in sig_AGB_spe) {
  select_data <- trait_data[trait_data$Species %in% i,]
  colnames(select_data)
  mod1 <- aov(log10(AGB_obs) ~ Seed_source, data = select_data)
  model_means <- emmeans(mod1, specs = ~ Seed_source)
  model_means_cld <- cld(object = model_means,
                         adjust = "none",
                         Letters = letters,
                         alpha = 0.05)
  result_cld1 = data.frame(Species = i , Seed_source = model_means_cld$Seed_source, group = model_means_cld$.group)
  cld_result_mod = rbind(cld_result_mod,result_cld1)
}
cld_result_mod

### plot
spe_select = (trait_data[trait_data$Species %in% select_sp_AGB,])
spe_select$AGB <- log10(spe_select$AGB_obs)
spe_select_mean <- aggregate(AGB ~ Seed_source + Species, data = spe_select, FUN = mean)

spe_select_mean$Seed_source <- factor(spe_select_mean$Seed_source, levels = c("Shandong","Henan","Hubei","Hunan","Guangxi","Guangdong"))
spe_select_mean_AGB <- spread(spe_select_mean, key = "Seed_source",value = "AGB")

spe_select_mean_AGB$latin <- gsub("_", " ", spe_select_mean_AGB$Species)
rownames(spe_select_mean_AGB) = spe_select_mean_AGB$latin
spe_select_mean_AGB$latin = NULL; spe_select_mean_AGB$Species = NULL
spe_select_mean_AGB$p_value = NA
pheatmap(spe_select_mean_AGB,
         annotation_row = NULL,
         annotation_names_row = TRUE,
         labels_row = rownames(spe_select_mean_AGB),
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
         color = colorRampPalette(colors = c("#F0F6E4","#414887"))(40),###You can change the heatmap by changing colors 
         fontsize_col = colTextSize,
         fontsize_row = rowTextSize,
         cellwidth = cellWidth,
         cellheight = cellHeight,
         filename = NA,
         cluster_cols = FALSE,
         cluster_rows = FALSE)


################################################################################
####SLA
trait_data <- subset(trait_data_total, SLA_obs != "NA")
unique(trait_data$Species)

#view(trait_data[trait_data$Species %in% "Mosla_dianthera", ])
####
summary_sp_data <- (trait_data[,c("Species","Seed_source")])
summary_sp_data$count = 1
summary_sp_data <- Rmisc::summarySE(summary_sp_data, measurevar = "count", groupvars = c("Species","Seed_source"))
summary_sp_data <- summary_sp_data[summary_sp_data$N >2, ]
unique(summary_sp_data$Species)

trait_data = trait_data %>% left_join(summary_sp_data[,c("Species","Seed_source","N")], by = c("Species","Seed_source"))
trait_data = trait_data %>% tidyr::drop_na(N)
unique(trait_data$Species)

summary_sp_data <- summary_sp_data %>% group_by(Species) %>% summarise(sum_count = sum(count))
select_sp_SLA <- unique(summary_sp_data[summary_sp_data$sum_count > 1, ]$Species)

final_result_mod = NULL
for (i in select_sp_SLA) {
  select_data <- trait_data[trait_data$Species %in% i,]
  colnames(select_data)
  mod1 <- aov(sqrt(SLA_obs) ~ Seed_source, data = select_data)
  summary_mod1 <- summary(mod1)
  result_mod1 <- data.frame(Species = i, Traits = "SLA", F_value = summary_mod1[[1]]$`F value`[1], p_value = summary_mod1[[1]]$"Pr(>F)"[1])
  final_result_mod = rbind(final_result_mod,result_mod1)
}

final_result_mod$p_value = round(final_result_mod$p_value, 4)
View(subset(final_result_mod, p_value < 0.05))


######Species with remarkable differences
sig_SLA_spe <- subset(final_result_mod, p_value < 0.05)$Species
cld_result_mod = NULL
for (i in sig_SLA_spe) {
  select_data <- trait_data[trait_data$Species %in% i,]
  colnames(select_data)
  mod1 <- aov(sqrt(SLA_obs) ~ Seed_source, data = select_data)
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

### plot
spe_select = (trait_data[trait_data$Species %in% select_sp_SLA,])
spe_select$SLA <- sqrt(spe_select$SLA_obs)
spe_select_mean <- aggregate(SLA ~ Seed_source + Species, data = spe_select, FUN = mean)

spe_select_mean$Seed_source <- factor(spe_select_mean$Seed_source, levels = c("Shandong","Henan","Hubei","Hunan","Guangxi","Guangdong"))
spe_select_mean_trans <- spread(spe_select_mean, key = "Seed_source",value = "SLA")

### 
spe_select_mean_trans2 = data.frame(Species = gsub(" ", "_", rownames(spe_select_mean_AGB)))
spe_select_mean_trans2 = spe_select_mean_trans2 %>% left_join(spe_select_mean_trans, by = "Species")
#spe_select_mean_trans$Species %in% gsub(" ", "_", rownames(spe_select_mean_AGB))

spe_select_mean_trans2$latin <- gsub("_", " ", spe_select_mean_trans2$Species)
rownames(spe_select_mean_trans2) = spe_select_mean_trans2$latin
spe_select_mean_trans2$latin = NULL; spe_select_mean_trans2$Species = NULL
spe_select_mean_trans2$p_value = NA
pheatmap(spe_select_mean_trans2,
         annotation_row = NULL,
         annotation_names_row = TRUE,
         labels_row = rownames(spe_select_mean_trans2),
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
         color = colorRampPalette(colors = c("#F0F6E4","#D9352A"))(40),###You can change the heatmap by changing colors 
         fontsize_col = colTextSize,
         fontsize_row = rowTextSize,
         cellwidth = cellWidth,
         cellheight = cellHeight,
         filename = NA,
         cluster_cols = FALSE,
         cluster_rows = FALSE)


##color we selected: SLA:#D9352A    Hmax:#E29C35   AGB:#414887

### bar plot
bar_data = data.frame(var = c("SLA", "Hmax", "AGB", "Not_sig"),
                      count_sum = c(5,14,12,87))
library(ggplot2)
bar_data$var = factor(bar_data$var, levels = c("SLA","Hmax","AGB","Not_sig"))
ggplot(bar_data, aes( x = "", y = count_sum,fill = var))+
  geom_col(position = 'stack', width = 0.6, color = "black") + 
  scale_fill_manual(values = c("#DC5246", "#E2A03D", "#525990","#998675")) +
  theme_classic() + 
  theme(axis.title = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(), legend.position = "none")
