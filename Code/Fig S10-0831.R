library(openxlsx)
library(lmerTest)
library(lme4)
library(MuMIn)
library(car)
library(AICcmodavg)
library(ggplot2)
field_com_data = read.xlsx("Data/all_row_data0829.xlsx", sheet = "field_data_mean", rowNames = F, colNames = T)
head(field_com_data)

### 
Field_trait = read.xlsx("Data/Field_traits_mean0831_2.xlsx", sheet = "Field_means", colNames = TRUE, rowNames = FALSE)
colnames(Field_trait) <- paste0( "Field_", colnames(Field_trait))
colnames(Field_trait)[1] <- "Species"

### 
Pot_trait = read.xlsx("Data/Pot_traits_mean0831.xlsx", sheet = "Pot_means", colNames = TRUE, rowNames = FALSE)
colnames(Pot_trait) <- paste0("Pot_", colnames(Pot_trait))
colnames(Pot_trait)[1] <- "Species"
nrow(Pot_trait)
summary(Pot_trait)
### 
field_com_data = field_com_data %>% left_join(Field_trait[,c(1:5)], by = "Species") %>%
  left_join(Pot_trait, by = "Species")
colnames(field_com_data)
field_com_data$rebio2022_100 = log10(field_com_data$rebio2022*100)
field_com_data$Block = as.factor(field_com_data$Block)
length(unique(field_com_data$Species))
field_com_data$Seed_source <- factor(field_com_data$Seed_source, levels = rev(c("Guangdong","Guangxi","Hunan","Hubei","Henan","Shandong")))

#### 
field_com_data$Field_SLA_imp_sqrt = sqrt(field_com_data$Field_SLA_imp)
field_com_data$Field_Hmax_sqrt = sqrt(field_com_data$Field_Hmax)
field_com_data$Field_AGB_lg = log10(field_com_data$Field_AGB)
##
field_com_data$Pot_SLA_sqrt = sqrt(field_com_data$Pot_SLA)
field_com_data$Pot_Hmax_sqrt = sqrt(field_com_data$Pot_Hmax)
field_com_data$Pot_Hmax_time1_sqrt = log10(field_com_data$Pot_Hmax_time1)
field_com_data$Pot_Hmax_time2_sqrt = log10(field_com_data$Pot_Hmax_time2)
field_com_data$Pot_AGB_lg = log10(field_com_data$Pot_AGB)

##
field_com_data$exist_prob <- ifelse(!is.na(field_com_data$rebio2022) & field_com_data$rebio2022 > 0, 1, 0)

################################################################################
#### odds of persistence
Seed_source = unique(field_com_data$Seed_source)
all_summary_glmer = NULL
for (i in Seed_source) {
  select_data = subset(field_com_data, Seed_source == i)
  mod2 <- glmer(exist_prob ~ Pot_Hmax_time1_sqrt + (1|Block) , na.action=na.omit, family=binomial, data=select_data)
  mod2_summary = data.frame(Pop = i, Traits = "Hmax", N = nrow(select_data), sp_num = length(unique(select_data$Species)),
                            Chisq = as.data.frame(Anova(mod2))[,1], p_value = as.data.frame(Anova(mod2))[,3],
                            R2m = r.squaredGLMM(mod2)[1,1], R2c = r.squaredGLMM(mod2)[1,2])
  all_summary_glmer = rbind(all_summary_glmer, mod2_summary)
}
all_summary_glmer$R2m = round(all_summary_glmer$R2m, 3)

##### 
mod2 = glmer(exist_prob ~ Pot_Hmax_time1_sqrt + (1|Block/Plot_num) , na.action=na.omit, family=binomial, data=field_com_data)
Anova(mod2); r.squaredGLMM(mod2)[1,]

nrow(field_com_data); length(unique(field_com_data$Species))

##### relative abundance of species in second year
Seed_source = unique(field_com_data$Seed_source)
all_summary_lmer = NULL
for (i in Seed_source) {
  select_data = subset(field_com_data, Seed_source == i) %>% tidyr::drop_na(rebio2020_100)
  mod2 <- lmer(rebio2022_100 ~ Pot_Hmax_time1_sqrt + (1|Block) , na.action=na.omit, data=select_data)
  mod2_summary = data.frame(Pop = i, Traits = "Hmax", N = nrow(select_data), sp_num = length(unique(select_data$Species)),
                            Chisq = as.data.frame(Anova(mod2))[,1], p_value = as.data.frame(Anova(mod2))[,3],
                            R2m = r.squaredGLMM(mod2)[1,1], R2c = r.squaredGLMM(mod2)[1,2])
  all_summary_lmer = rbind(all_summary_lmer, mod2_summary)
}
all_summary_lmer$R2m = round(all_summary_lmer$R2m, 3)
##### 
str(select_data)
mod2 = lmer(rebio2022_100 ~ Pot_Hmax_time1_sqrt + (1|Block/Plot_num), data=subset(field_com_data, rebio2022_100 != "NA"))
Anova(mod2); r.squaredGLMM(mod2)[1,]

select_data = field_com_data %>% tidyr::drop_na(rebio2022_100)
nrow(select_data); length(unique(select_data$Species))

######### 
length(unique(field_com_data$Species))
unique(field_com_data$Seed_source)

#### Hmax
mod12 <- glmer(exist_prob ~ Pot_Hmax_time1_sqrt*Seed_source + (1|Block) , na.action=na.omit, family=binomial, data=field_com_data)
Anova(mod12)
summary(mod12)
MuMIn::r.squaredGLMM(mod12)

field_com_data$F0 = predictSE(mod12, field_com_data, level = 0)$fit
field_com_data$SE <- predictSE(mod12, field_com_data, level = 0)$se.fit

####
field_com_data_all = field_com_data
mod12 <- glmer(exist_prob ~ Pot_Hmax_time1_sqrt + (1|Block/Plot_num) , na.action=na.omit, family=binomial, data=field_com_data_all)
Anova(mod12)
field_com_data_all$F0 = predictSE(mod12, field_com_data_all, level = 0)$fit
field_com_data_all$SE <- predictSE(mod12, field_com_data_all, level = 0)$se.fit

ggplot(field_com_data, aes(x=(Pot_Hmax_time1_sqrt), y=exist_prob)) +
  #geom_ribbon(aes(x=Pot_Hmax_time1_sqrt,ymin = F0 - 1.96 * SE,ymax = F0 + 1.96 * SE, fill = Seed_source), alpha = I(0.1)) +
  geom_line(data = subset(field_com_data, Seed_source != "Hunan"), mapping = aes(x=Pot_Hmax_time1_sqrt,y=F0, color = Seed_source), size=1, linetype = 2) + #color = "#8B0000"
  geom_line(data = subset(field_com_data, Seed_source == "Hunan"), mapping = aes(x=Pot_Hmax_time1_sqrt,y=F0, color = Seed_source), size=1, linetype = 1) + #color = "#8B0000"
  geom_line(data = field_com_data_all, mapping = aes(x=Pot_Hmax_time1_sqrt, y=F0), size=1.5, linetype = 2, color = "black") +
  geom_ribbon(data = field_com_data_all, mapping = aes(x=Pot_Hmax_time1_sqrt, ymin = F0 - 1.96 * SE,ymax = F0 + 1.96 * SE, fill = "#EBEBEB"), alpha = I(0.1)) +
  geom_point(aes(fill = Seed_source),size = 2.2, pch = 21) + # fill = "#00000022", color = "#0A0A0A"
  labs(x = 'Maximum height estimated\nin first time of pot experiment',y = 'Odds of persistence',title = NULL) +  
  mytheme + theme(legend.position = "none") + 
  scale_fill_manual(values=c("Shandong" = "#376694", "Henan" = "#73BCD5", "Hubei" = "#ABDCE0",
                             "Hunan" = "#EFBA55","Guangxi" = "#E68C51","Guangdong" = "#994240"))+
  scale_color_manual(values=c("Shandong" = "#376694", "Henan" = "#73BCD5", "Hubei" = "#ABDCE0",
                              "Hunan" = "#EFBA55","Guangxi" = "#E68C51","Guangdong" = "#994240"))+
  scale_x_continuous(labels = scales::label_comma(accuracy =0.1)) +
  scale_y_continuous(labels = scales::label_comma(accuracy =0.01)) -> p2; p2

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

#### Hmax
mod12 <- lmer(rebio2022_100 ~ Pot_Hmax_time1_sqrt*Seed_source + (1|Block), data = field_com_data_no)
summary(mod12)
anova(mod12)
MuMIn::r.squaredGLMM(mod12)
field_com_data_no$F0 = predictSE(mod12, field_com_data_no, level = 0)$fit
field_com_data_no$SE <- predictSE(mod12, field_com_data_no, level = 0)$se.fit

####
field_com_data_all = field_com_data_no
mod12 <- lmer(rebio2022_100 ~ Pot_Hmax_time1_sqrt + (1|Block/Plot_num), data = field_com_data_all)
Anova(mod12)
field_com_data_all$F0 = predictSE(mod12, field_com_data_all, level = 0)$fit
field_com_data_all$SE <- predictSE(mod12, field_com_data_all, level = 0)$se.fit

ggplot(field_com_data_no, aes(x=Pot_Hmax_time1_sqrt, y=rebio2022_100)) +
  #geom_ribbon(aes(x=Pot_Hmax_time1_sqrt,ymin = F0 - 1.96 * SE,ymax = F0 + 1.96 * SE, fill = Seed_source), alpha = I(0.1)) +
  geom_line(data = subset(field_com_data_no, Seed_source != "Guangdong" & Seed_source != "Hunan"), mapping = aes(x=Pot_Hmax_time1_sqrt, y=F0, color = Seed_source), size=1, linetype = 2) +
  geom_line(data = subset(field_com_data_no, Seed_source == "Guangdong"), mapping = aes(x=Pot_Hmax_time1_sqrt, y=F0, color = Seed_source), size=1, linetype = 1) +
  geom_line(data = subset(field_com_data_no, Seed_source == "Hunan"), mapping = aes(x=Pot_Hmax_time1_sqrt, y=F0, color = Seed_source), size=1, linetype = 1) +
  geom_line(data = field_com_data_all, mapping = aes(x=Pot_Hmax_time1_sqrt, y=F0), size=1.5, linetype = 1, color = "black") +
  geom_ribbon(data = field_com_data_all, mapping = aes(x=Pot_Hmax_time1_sqrt, ymin = F0 - 1.96 * SE,ymax = F0 + 1.96 * SE, fill = "#EBEBEB"), alpha = I(0.1)) +
  geom_point(aes(fill = Seed_source),size = 2.2, color = "black", pch = 21) + # fill = "#00000022", color = "#0A0A0A"
  scale_fill_manual(values=c("Shandong" = "#376694", "Henan" = "#73BCD5", "Hubei" = "#ABDCE0",
                             "Hunan" = "#EFBA55","Guangxi" = "#E68C51","Guangdong" = "#994240"))+
  scale_color_manual(values=c("Shandong" = "#376694", "Henan" = "#73BCD5", "Hubei" = "#ABDCE0",
                              "Hunan" = "#EFBA55","Guangxi" = "#E68C51","Guangdong" = "#994240"))+
  ggrepel::geom_text_repel(mapping = aes(x=Pot_Hmax_time1_sqrt,y=rebio2022_100,label=Latin_name), data = df2,size = 2.8,segment.color = "black", color = "black",direction = "both",box.padding = 0.6,
                           max.overlaps = getOption("ggrepel.max.overlaps", default = 25), fontface = "italic") +
  scale_x_continuous(labels = scales::label_comma(accuracy =0.1)) +
  scale_y_continuous(labels = scales::label_comma(accuracy =0.01)) +
  labs(x = 'Maximum height estimated\nin first time of pot experiment',
       y = 'Relative abundance within the\nplot in second year (%, log10)',title = NULL) +  
  mytheme -> p5; p5

(p2|p5)/(p2|p5)
