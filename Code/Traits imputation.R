### Inferring and interpolating missing data on specific leaf area
#Imputation is based on correlations among traits and species phylogenetic relatedness
### Field trait mean
library(openxlsx)
library(ape)
library(adephylo)
library(missForest)
library(ggplot2)
library(dplyr)
library(tidyverse)

###
mytheme = theme_bw()+
  theme( panel.background = element_rect(fill='white', colour='black'),
         panel.grid=element_blank(), 
         legend.position = "none",
         axis.ticks.length = unit(0.4,"lines"), axis.ticks = element_line(color='black'),
         axis.line = element_line(colour = "black"), 
         axis.title.x=element_text(colour='black', size=13,vjust = 1),
         axis.title.y=element_text(colour='black', size=13,vjust = 1),
         axis.text=element_text(colour='black',size=11))

Field_traits <- read.xlsx("Data/Imputation_data0831_2.xlsx",sheet = "Imputation_data", rowNames = FALSE, colNames = TRUE)
rownames(Field_traits) = Field_traits$Species
### Analyzing species lists (SLA = 49, Hmax&AGB = 64)
Common_sp_list = read.xlsx("Data/Common_species_list.xlsx", sheet = "Common_sp_list", colNames = TRUE, rowNames = FALSE)
Common_sp_list_SLA = c(na.omit(Common_sp_list$Species_SLA))
Common_sp_list_AGB = c(na.omit(Common_sp_list$Species_AGB))
Field_traits <- Field_traits[Field_traits$Species %in% Common_sp_list_AGB, ]

head(Field_traits)
summary(Field_traits)
### N = 88 (Species)
phylogenyAux <- read.tree("Data/iq_tree.treefile")
to_drop = setdiff(as.vector(phylogenyAux$tip.label),as.vector((Field_traits$Species)))
phy_tree <- drop.tip(as.phylo(phylogenyAux), to_drop) 
plot(phy_tree)
#reorder species according to phylogenetic tree
Field_traits <- Field_traits[phy_tree$tip.label, ]
nrow(Field_traits)

#create phylogenetic proximity table
prox.Ab.all <- proxTips(phy_tree, method = "Abouheif", normalize="none")
prox.Ab.all[1:5,1:5]
dim(prox.Ab.all)

#Calculate Moran eigenvetors to filter out phylogenetic autocorrelation
prox <- prop.table(prox.Ab.all, 1) #standardize by row
prox <- 0.5 * (prox + t(prox)) #make matrix symetric
prox[1:5,1:5]

ME <- me.phylo(prox = prox)
#??adephylo::me.phylo
ME <- ME[rownames(Field_traits),]
dim(ME)
head(ME)

###
#pvr_tree <- PVRdecomp(phy_tree, scale = TRUE)
#eigenvec <- data.frame(Species=phy_tree$tip.label, pvr_tree@Eigen$vectors)
#rownames(eigenvec) = eigenvec$Species
#eigenvec$Species = NULL
#eigenvec <- eigenvec[rownames(Field_traits),]

##impute missing trait values using missForest function with default settings
set.seed(1234)
trait_miss <- cbind(log10(Field_traits[,2:4]), ME[,1:30])
#trait_miss <- cbind(log10(Field_traits[,2:4]), eigenvec[,1:30])
trait.imp <- missForest(trait_miss)
#trait.imp <- missForest(cbind(log10(Field_traits[,2:4]), ME[,2:31]))
(imp.error <- trait.imp$OOBerror)#check OOB rerror
trait.imp <- as.data.frame(10^trait.imp$ximp[,1:3])
trait.imp <- trait.imp[phy_tree$tip.label, ]
trait.imp$Species = rownames(trait.imp)
colnames(trait.imp)[3] = "Field_SLA_imp"
# write.csv(trait.imp,"trait.imp0901.csv")


################################################################################
Field_traits <- read.xlsx("Data/Imputation_data0831_2.xlsx",sheet = "Imputation_data", rowNames = FALSE, colNames = TRUE)
rownames(Field_traits) = Field_traits$Species
Field_traits <- Field_traits %>% drop_na(SLA); nrow(Field_traits)

phylogenyAux <- read.tree("Data/iq_tree.treefile")
to_drop = setdiff(as.vector(phylogenyAux$tip.label),as.vector((Field_traits$Species)))
phy_tree <- drop.tip(as.phylo(phylogenyAux), to_drop) 
plot(phy_tree)
#reorder species according to phylogenetic tree
Field_traits <- Field_traits[phy_tree$tip.label, ]
nrow(Field_traits)

#create phylogenetic proximity table
prox.Ab.all <- proxTips(phy_tree, method = "Abouheif", normalize="none")
prox.Ab.all[1:5,1:5]
dim(prox.Ab.all)

#Calculate Moran eigenvetors to filter out phylogenetic autocorrelation
prox <- prop.table(prox.Ab.all, 1) #standardize by row
prox <- 0.5 * (prox + t(prox)) #make matrix symetric
prox[1:5,1:5]

ME <- me.phylo(prox = prox)
#??adephylo::me.phylo
ME <- ME[rownames(Field_traits),]
dim(ME)
head(ME)

###
#pvr_tree <- PVRdecomp(phy_tree, scale = TRUE)
#eigenvec <- data.frame(Species=phy_tree$tip.label, pvr_tree@Eigen$vectors)
#rownames(eigenvec) = eigenvec$Species
#eigenvec$Species = NULL
#eigenvec <- eigenvec[rownames(Field_traits),]

##impute missing trait values using missForest function with default settings
trait_miss <- cbind(log10(Field_traits[,2:4]), ME[,1:30])

### use 64 species to compare the methods
t<-na.omit(trait_miss)
t$Species <- rownames(t)
dim(t)
#write.csv(t, "t.csv")

set.seed(1234)
res<-vector("list",100)

for (ii in 1:100) {
  set350<-t[sample(1:dim(t)[1],64),]
  set350<-set350[order(set350$Species),]
  tr<-set350
  na<-15 ### use average NA percent for each trait
  for(i in 1:3)
  {tr[sample(1:64,na),i]<-NA}
  # check  NA %|
  N1<-apply(tr,2,function(x) length(na.omit(x)))
  1-sum(N1[1:3])/64/3
  1-N1[1:3]/64
  
  # use misforest fill the NA
  rff<-tr
  N <- ncol(rff)
  rf<-missForest::missForest(rff[,-N]) ### 
  rf_350 <- as.data.frame(10^rf$ximp[,1:3])
  
  # use mean
  Mean_350<-tr ## for using mean values
  M<-apply(tr[,1:3],2,function(x) mean(na.omit(x)))
  
  ###
  for(i in 1:3){
    Mean_350[which(is.na(Mean_350[,i])),i]<-M[i]}
  
  ###
  yy<-NULL
  for(i in 1:3){
    x<-cor.test(set350[,i],rf_350[,i],method="spearman")
    z<-cor.test(set350[,i],Mean_350[,i],method="spearman")
    y<-c(x$estimate,z$estimate)#y$estimate,
    yy<-rbind(yy,y)
  }
  ###
  p <- as.data.frame(yy)
  colnames(p)<-c("RandomForest","Average")
  p$name<-colnames(tr)[1:3]
  res[[ii]]<- p
  print(paste0("Sloop", ii))
}

res

combined_df <- do.call(rbind, res)
rownames(combined_df) = NULL
long_df <- combined_df %>% gather(key = "Methods", value = "rho", -name)
long_df$name = str_replace_all(long_df$name, "_", " ")

long_df$name<-factor(long_df$name,levels = c("SLA","Hmax","AGB"))
long_df$Methods<-factor(long_df$Methods,levels=c("RandomForest","Average"))

ggplot(long_df,aes(x=name,y=rho))+
  geom_violin(data=long_df, aes(x=name,y=rho,fill=Methods,color=Methods),trim=T,scale = "width",
              position = position_dodge(0.7),width=0.6,alpha = 0.8,linewidth = 1.5)+
  scale_color_manual(values = c("#EFF0F0","#A0A1A1")) + 
  scale_fill_manual(values = c("#EFF0F0","#A0A1A1")) + 
  ggnewscale::new_scale_fill() +
  ggnewscale::new_scale_color() + 
  geom_boxplot(data=long_df, aes(x=name,y=rho,fill=Methods,color=Methods ),
               position = position_dodge(0.7),width=0.15,outlier.size=0.8,outlier.shape = 1)+
  stat_summary(data=long_df, aes(x=name,y=rho,fill=Methods,color=Methods),
               fun="mean",position = position_dodge(0.7),geom="point",shape=21, size=1.5)+
  scale_color_manual(values = c('black','black'))+
  scale_fill_manual(values = c('white','white'))+
  mytheme+
  xlab("")+ylab("Spearman rho")+
  scale_x_discrete(labels = c("SLA" = "Specific leaf area", "Hmax" = "Maximum height", "AGB" = "Aboveground biomass")) + 
  theme(legend.position = "right",
        axis.text.x = element_text(colour='black', size=11, angle = 35, hjust = 1, vjust = 1))

rf<-apply(sapply(res,function(x) x$RandomForest),1,mean)
m<-apply(sapply(res,function(x) x$Average),1,mean)
rf.sd<-apply(sapply(res,function(x) x$RandomForest),1,function(y) sd(y)/sqrt(length(y)))
m.sd<-apply(sapply(res,function(x) x$Average),1,function(y) sd(y)/sqrt(length(y)))

# data
dat<-data.frame(name=rep(res[[1]]$name,2),group=rep(c("rf","mean"),each=3),mean=c(rf,m),sd=c(rf.sd,m.sd))
dat$name<-factor(dat$name,levels = c("SLA","Hmax","AGB"))
dat$group<-factor(dat$group,levels=c("rf","mean"))

## plot
tapply(dat$mean,dat$group,range)

ff<-ggplot(data=dat,aes(x=name,y=mean,ymin = mean-sd, ymax = mean+sd,fill=group))+
  geom_bar(stat="identity",position="dodge")+
  geom_errorbar(position=position_dodge(width=0.9), width = 0.2)+
  coord_cartesian(ylim = c(0.7,0.95))+
  scale_y_continuous(expand = c(0,0))+
  xlab("")+ylab("Spearman rho")+
  mytheme + 
  scale_fill_discrete(name="Methods",labels=c('RandomForest','Average'))+
  scale_x_discrete(labels = c("SLA" = "Specific leaf area", "Hmax" = "Maximum height", "AGB" = "Aboveground biomass")) + 
  theme(legend.position = c(0.15,0.85),
        axis.text.x = element_text(colour='black', size=11, angle = 35, hjust = 1, vjust = 1)) ; ff






