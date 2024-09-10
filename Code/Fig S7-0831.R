################################################################################
################################## Fig S7 ######################################
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
pot_trait = read.xlsx("Data/Pot_traits_mean0831.xlsx", sheet = "Pot_means", colNames = TRUE, rowNames = FALSE)
colnames(pot_trait) <- paste0(colnames(pot_trait), "_pot")
colnames(pot_trait)[1] <- "Species"

### Loading field experiment database 
### AGB
trait_data = read.xlsx("Data/Field_traits_mean0831_2.xlsx", sheet = "Field_means", colNames = TRUE, rowNames = FALSE)
colnames(trait_data) <- paste0(colnames(trait_data), "_field")
colnames(trait_data)[1] <- "Species"

trait_data = trait_data %>% left_join(pot_trait, by = "Species") %>% 
  left_join(Common_sp_list[,c(1,2,5)], by = "Species") %>%
  drop_na(AGB_field) %>% drop_na(AGB_pot)
length(unique(trait_data$Species))
trait_data$Origin = factor(trait_data$Origin, levels = c("Native", "Exotic"))

first_char <- substr(trait_data$Species, 1, 1)
sub_str <- gsub(".*_", "", trait_data$Species)
Latin_name <- paste(first_char, sub_str, sep = ". ")
trait_data$Latin_name = Latin_name

summary(lm(log10(AGB_field)~ log10(AGB_pot), data = trait_data))

df2 = (trait_data %>% arrange(desc(AGB_field)))[c(1:6),]
df3 = (trait_data %>% arrange((AGB_field)))[c(1:6),]

P_AGB = ggplot(trait_data, mapping = aes(x = log10(AGB_pot) , y = log10(AGB_field))) + 
  geom_smooth(data = trait_data, mapping = aes(x=log10(AGB_pot),y=log10(AGB_field)),method = "lm", se = T, color = "black", fill = "grey80", linetype = 1) + 
  geom_point(trait_data, mapping = aes(x = log10(AGB_pot) , y = log10(AGB_field), shape = Origin, fill = Origin),size=2.2, color = "black")+
  scale_shape_manual(values = c(21,21))+
  scale_color_manual(values=c('#60A7A6','#FEA6A6'))+
  scale_fill_manual(values=c('#60A7A6','#FEA6A6'))+
  labs(x = NULL,
       y = NULL,
       title = "Aboveground biomass (g, log10)") + 
  mytheme + 
  scale_x_continuous(labels = scales::label_comma(accuracy =0.1)) +
  scale_y_continuous(labels = scales::label_comma(accuracy =0.1)) + 
  #geom_abline(intercept= 1.4322 ,slope=0.7259 ,size=.8, linetype = 1)+
  ggrepel::geom_text_repel(aes(label=Latin_name), df2,size = 2.5,segment.color = "black", 
                           color = "black",direction = "both",box.padding = 0.7, fontface = "italic",
                           max.overlaps = getOption("ggrepel.max.overlaps", default = 20))+
  ggrepel::geom_text_repel(aes(label=Latin_name), df3,size = 2.5,segment.color = "black", 
                           color = "black",direction = "both",box.padding = 0.7, fontface = "italic",
                           max.overlaps = getOption("ggrepel.max.overlaps", default = 20))
P_AGB

### Hmax
trait_data = read.xlsx("Data/Field_traits_mean0831_2.xlsx", sheet = "Field_means", colNames = TRUE, rowNames = FALSE)
colnames(trait_data) <- paste0(colnames(trait_data), "_field")
colnames(trait_data)[1] <- "Species"

trait_data = trait_data %>% left_join(pot_trait, by = "Species") %>% 
  left_join(Common_sp_list[,c(1,2,5)], by = "Species") %>%
  drop_na(Hmax_field) %>% drop_na(Hmax_pot)
length(unique(trait_data$Species))

colnames(trait_data)
trait_data$Origin = factor(trait_data$Origin, levels = c("Native", "Exotic"))
summary(lm(sqrt(Hmax_field)~ sqrt(Hmax_pot), data = trait_data))
summary(lm(sqrt(Hmax_field)~ sqrt(Hmax_time1_pot), data = trait_data))
#plot(ma.test)
first_char <- substr(trait_data$Species, 1, 1)
sub_str <- gsub(".*_", "", trait_data$Species)
Latin_name <- paste(first_char, sub_str, sep = ". ")
trait_data$Latin_name = Latin_name

df2 = (trait_data %>% arrange(desc(Hmax_field)))[c(1:6),]
df3 = (trait_data %>% arrange((Hmax_field)))[c(1:6),]

####
P_Hmax = ggplot(trait_data, mapping = aes(x = sqrt(Hmax_pot) , y = sqrt(Hmax_field))) + 
  geom_smooth(data = trait_data, mapping = aes(x=sqrt(Hmax_pot),y=sqrt(Hmax_field)),method = "lm", se = T, color = "black", fill = "grey80", linetype = 1) + 
  geom_point(aes(shape = Origin, fill = Origin),size=2.2, color = "black")+
  scale_shape_manual(values = c(21,21))+
  scale_color_manual(values=c('#60A7A6','#FEA6A6'))+
  scale_fill_manual(values=c('#60A7A6','#FEA6A6'))+
  labs(x = 'Traits value estimated in pot experiment', y = NULL, title = "Maximum height") + 
  mytheme + 
  scale_x_continuous(labels = scales::label_comma(accuracy =0.1)) +
  ggrepel::geom_text_repel(aes(label=Latin_name), df2,size = 2.5,segment.color = "black",
                           color = "black",direction = "both",box.padding = 0.7, fontface = "italic",
                           max.overlaps = getOption("ggrepel.max.overlaps", default = 15))+
  ggrepel::geom_text_repel(aes(label=Latin_name), df3,size = 2.5,segment.color = "black",
                           color = "black",direction = "both",box.padding = 0.7, fontface = "italic",
                           max.overlaps = getOption("ggrepel.max.overlaps", default = 15))
P_Hmax

###
summary(lm(sqrt(Hmax_field)~ sqrt(Hmax_time1_pot), data = trait_data))

P_Hmax1 = ggplot(trait_data, mapping = aes(x = sqrt(Hmax_time1_pot) , y = sqrt(Hmax_field))) + 
  geom_smooth(data = trait_data, mapping = aes(x=sqrt(Hmax_time1_pot),y=sqrt(Hmax_field)),method = "lm", se = T, color = "black", fill = "grey80", linetype = 1) + 
  geom_point(aes(shape = Origin, fill = Origin),size=2.2, color = "black")+
  scale_shape_manual(values = c(21,21))+
  scale_color_manual(values=c('#60A7A6','#FEA6A6'))+
  scale_fill_manual(values=c('#60A7A6','#FEA6A6'))+
  labs(y = 'Maximum height estimated in field experiment',
       x = 'Maximum height estimated in first time of pot experiment',
       title = "Maximum height") + 
  #geom_abline(intercept=6.1373 ,slope=0.7809 ,size=.8, linetype = 1)+
  mytheme + 
  scale_x_continuous(labels = scales::label_comma(accuracy =0.1)) +
  ggrepel::geom_text_repel(aes(label=Latin_name), df2,size = 2.5,segment.color = "black",
                           color = "black",direction = "both",box.padding = 0.7, fontface = "italic",
                           max.overlaps = getOption("ggrepel.max.overlaps", default = 15))+
  ggrepel::geom_text_repel(aes(label=Latin_name), df3,size = 2.5,segment.color = "black",
                           color = "black",direction = "both",box.padding = 0.7, fontface = "italic",
                           max.overlaps = getOption("ggrepel.max.overlaps", default = 15))+
  annotate('text', label = expression(italic(R)^2~"="~0.277*","~italic(p)~"<"~0.001), x = 4, y = 13, size = 3, parse = TRUE)


P_Hmax1

summary(lm(sqrt(Hmax_field)~ sqrt(Hmax_time2_pot), data = trait_data))

P_Hmax2 = ggplot(trait_data, mapping = aes(x = sqrt(Hmax_time2_pot) , y = sqrt(Hmax_field))) + 
  geom_smooth(data = trait_data, mapping = aes(x=sqrt(Hmax_time2_pot),y=sqrt(Hmax_field)),method = "lm", se = T, color = "black", fill = "grey80", linetype = 1) + 
  geom_point(aes(shape = Origin, fill = Origin),size=2.2, color = "black")+
  scale_shape_manual(values = c(21,21))+
  scale_color_manual(values=c('#60A7A6','#FEA6A6'))+
  scale_fill_manual(values=c('#60A7A6','#FEA6A6'))+
  labs(y = NULL,
       x = 'Traits value estimated in pot experiment (second)',
       title = "Maximum height") + 
  #geom_abline(intercept=1.1498689 ,slope=0.6356759 ,size=.8, linetype = 1)+
  mytheme + 
  scale_x_continuous(labels = scales::label_comma(accuracy =0.1)) +
  scale_y_continuous(position = "left",labels = scales::label_comma(accuracy =0.1), 
                     limits = c(ggplot_build(P_Hmax1)$layout$panel_scales_y[[1]]$range$range)) +
  ggrepel::geom_text_repel(aes(label=Latin_name), df2,size = 2.5,segment.color = "black",
                           color = "black",direction = "both",box.padding = 0.7, fontface = "italic",
                           max.overlaps = getOption("ggrepel.max.overlaps", default = 15))+
  ggrepel::geom_text_repel(aes(label=Latin_name), df3,size = 2.5,segment.color = "black",
                           color = "black",direction = "both",box.padding = 0.7, fontface = "italic",
                           max.overlaps = getOption("ggrepel.max.overlaps", default = 15))+
  annotate('text', label = expression(italic(R)^2~"="~0.203*","~italic(p)~"<"~0.001), x = 3, y = 13, size = 3, parse = TRUE)


P_Hmax2

P_Hmax1|P_Hmax2

### SLA (row data)
trait_data = read.xlsx("Data/Field_traits_mean0831_2.xlsx", sheet = "Field_means", colNames = TRUE, rowNames = FALSE)
colnames(trait_data) <- paste0(colnames(trait_data), "_field")
colnames(trait_data)[1] <- "Species"

trait_data = trait_data %>% left_join(pot_trait, by = "Species") %>% 
  left_join(Common_sp_list[,c(1,2,5)], by = "Species") %>%
  drop_na(SLA_imp_field) %>% drop_na(SLA_pot) 
length(unique(trait_data$Species))
trait_data$Origin = factor(trait_data$Origin, levels = c("Native", "Exotic"))

#plot(ma.test)
first_char <- substr(trait_data$Species, 1, 1)
sub_str <- gsub(".*_", "", trait_data$Species)
Latin_name <- paste(first_char, sub_str, sep = ". ")
trait_data$Latin_name = Latin_name

df2 = (trait_data %>% arrange(desc(SLA_field)))[c(1:6),]
df3 = (trait_data %>% arrange((SLA_field)))[c(1:6),]

summary(lm(sqrt(SLA_field)~ sqrt(SLA_pot), data = trait_data))

P_SLA = ggplot(trait_data, mapping = aes(x = sqrt(SLA_pot) , y = sqrt(SLA_field))) + 
  geom_smooth(data = trait_data, mapping = aes(x=sqrt(SLA_pot),y=sqrt(SLA_field)),method = "lm", se = T, color = "black", fill = "grey80", linetype = 1) + 
  geom_point(trait_data, mapping = aes(x = sqrt(SLA_pot) , y = sqrt(SLA_field), shape = Origin, fill = Origin),size=2.2, color = "black")+
  scale_shape_manual(values = c(21,21))+
  scale_color_manual(values=c('#60A7A6','#FEA6A6'))+
  scale_fill_manual(values=c('#60A7A6','#FEA6A6'))+
  labs(x = NULL,
       y = bquote(atop('Traits value'~(log[10]), 'estimated in field experiment')),
       title = "Specific leaf area") + 
  mytheme + 
  ggrepel::geom_text_repel(aes(label=Latin_name), df2,size = 2.5,segment.color = "black", 
                           color = "black",direction = "both",box.padding = 0.7, fontface = "italic",
                           max.overlaps = getOption("ggrepel.max.overlaps", default = 15)) +
  ggrepel::geom_text_repel(aes(label=Latin_name), df3,size = 2.5,segment.color = "black", 
                           color = "black",direction = "both",box.padding = 0.7,fontface = "italic",
                           max.overlaps = getOption("ggrepel.max.overlaps", default = 15))
P_SLA

### (Inferred data)
summary(lm(sqrt(SLA_imp_field)~ sqrt(SLA_pot), data = trait_data))

df2 = (trait_data %>% arrange(desc(SLA_imp_field)))[c(1:6),]
df3 = (trait_data %>% arrange((SLA_imp_field)))[c(1:6),]

P_SLA = ggplot(trait_data, mapping = aes(x = sqrt(SLA_pot) , y = sqrt(SLA_imp_field))) + 
  geom_smooth(data = trait_data, mapping = aes(x=sqrt(SLA_pot),y=sqrt(SLA_imp_field)),method = "lm", se = T, color = "black", fill = "grey80", linetype = 1) + 
  geom_point(trait_data, mapping = aes(x = sqrt(SLA_pot) , y = sqrt(SLA_imp_field), 
                                       shape = Origin, fill = Origin),size=2.2, color = "black")+
  scale_shape_manual(values = c(21,21))+
  scale_color_manual(values=c('#60A7A6','#FEA6A6'))+
  scale_fill_manual(values=c('#60A7A6','#FEA6A6'))+
  labs(x = NULL,
       y = "Traits value estimated in field experiment",
       title = "Specific leaf area (cm2/g, sqrt)") + 
  mytheme + 
  ggrepel::geom_text_repel(aes(label=Latin_name), df2,size = 2.5,segment.color = "black", 
                           color = "black",direction = "both",box.padding = 0.7, fontface = "italic",
                           max.overlaps = getOption("ggrepel.max.overlaps", default = 100)) +
  ggrepel::geom_text_repel(aes(label=Latin_name), df3,size = 2.5,segment.color = "black",
                           color = "black",direction = "both",box.padding = 0.7, fontface = "italic",
                           max.overlaps = getOption("ggrepel.max.overlaps", default = 100))

P_SLA

### Fig S2
P_SLA+P_Hmax+P_AGB+plot_annotation(tag_levels = 'a')

