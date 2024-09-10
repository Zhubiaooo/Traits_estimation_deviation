################################################################################
### Loading pot experiment database 
pot_trait = read.xlsx("Data/Pot_traits_mean0831.xlsx", sheet = "Pot_means", colNames = TRUE, rowNames = FALSE)
#pot_trait = read.xlsx("Data/Pot_traits_database.xlsx", sheet = "Pot_data", colNames = TRUE, rowNames = FALSE)[,c("Species","SLA")]
pot_trait = subset(pot_trait, SLA != "NA")
colnames(pot_trait) <- paste0(colnames(pot_trait), "_pot")
colnames(pot_trait)[1] <- "Species"
#rownames(pot_trait) <- pot_trait$Species

### Loading Field experiment database 
field_trait = read.xlsx("Data/Field_traits_mean0831_2.xlsx", sheet = "Field_means", colNames = TRUE, rowNames = FALSE)
#field_trait = read.xlsx("Data/Field_traits_database.xlsx", sheet = "Field_data", colNames = TRUE, rowNames = FALSE)[,c("Species","SLA")]
#field_trait = subset(field_trait, SLA != "NA")
colnames(field_trait) <- paste0(colnames(field_trait), "_field")
colnames(field_trait)[1] <- "Species"
#rownames(field_trait) <- field_trait$Species

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

### Loading tree
phylogenyAux = read.tree("Data/iq_tree.treefile")
plot(phylogenyAux)

to_drop = phylogenyAux$tip.label[!phylogenyAux$tip.label %in% Common_sp_list$Species]
tree <- drop.tip(as.phylo(phylogenyAux), to_drop) 
phyloMat = vcv.phylo(tree)
phyloMat = phyloMat / max(phyloMat)
dim(phyloMat)
plot(tree)
####
aa = pot_trait[,c("Species","SLA_pot","Hmax_pot","AGB_pot")] %>%
  left_join(Common_sp_list[c("Species","Origin")], by = "Species")
colnames(aa)[c(2,3,4)] = c("SLA","Hmax","AGB")
aa$type = "pot"

bb = field_trait[,c("Species","SLA_imp_field","Hmax_field","AGB_field")] %>%
  left_join(Common_sp_list[c("Species","Origin")], by = "Species")
colnames(bb)[c(2,3,4)] = c("SLA","Hmax","AGB")
bb$type = "field"

aaa = rbind(aa, bb)
###
aaa = aaa %>% left_join(Common_sp_list[,c("Species","Origin")])
aaa$Origin = factor(aaa$Origin, levels = c("Native", "Exotic"))
mod <- relmatLmer(sqrt(SLA) ~ Origin*type + (1|Species), data = aaa, relmat = list(Species=phyloMat))
Anova(mod, type="II", test.statistic ="Chisq")
shapiro.test(residuals(mod))
aaa$type = factor(aaa$type, levels = c("pot","field"))
aaa$SLA_sqrt = sqrt(aaa$SLA)
### plot
ggplot(aaa,aes(x=type,y=sqrt(SLA)))+
  geom_violin(data=aaa, aes(x=type,y=sqrt(SLA),fill=type,color=type),trim=T,scale = "width",
              position = position_dodge(0.7),width=0.6,alpha = 0.8,linewidth = 1.5)+
  scale_color_manual(values = c("#A0A1A1","#A0A1A1")) + 
  scale_fill_manual(values = c("#A0A1A1","#A0A1A1")) + 
  ggnewscale::new_scale_fill() +
  ggnewscale::new_scale_color() + 
  geom_boxplot(data=aaa, aes(x=type,y=sqrt(SLA),fill=type,color=type ),
               position = position_dodge(0.7),width=0.15,outlier.size=0.8,outlier.shape = 1)+
  stat_summary(data=aaa, aes(x=type,y=sqrt(SLA),fill=type,color=type),
               fun="mean",position = position_dodge(0.7),geom="point",shape=21, size=1.5)+
  scale_color_manual(values = c('black','black'))+
  scale_fill_manual(values = c('white','white'))+
  mytheme+
  scale_x_discrete(labels = c("pot" = "Pot", "field" = "Field")) + 
  geom_signif(comparisons = list(c(1,2)),test="t.test", annotations='p < 0.001',tip_length = 0.02,size = 0.5,
              textsize = 4,y_position = 22) + 
  labs(x = NULL, y = expression('Specific leaf area (cm'^ 2*'/g, sqrt)'),  title = NULL)

### Percentage increase
Rmisc::summarySE(aaa, measurevar = c("SLA"), groupvars = c("type"))
(242.0664 - 175.6049)/175.6049 
