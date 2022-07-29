library(readxl)
library(dplyr)

# Figure 1a ------------------------------------------
# Figure 1a. taxonomy #######
rm(list = ls())
dat.taxa <- read_excel("Fig 1 Source Data.xlsx", sheet = "Fig 1a taxonomy PCA")
dat.funct <- read_excel("Fig 1 Source Data.xlsx", sheet = "Fig 1a metagenome PCA")

pca <- dat.taxa
Fig1a.taxa.pca <- ggplot(pca,aes(PC1,PC2))+
  geom_point(size=2,aes(col=Disease,shape=Cohort))+ 
  scale_color_manual(values=c("#f0999f","#46bbc0"))+
  scale_shape_manual(values=c(16,15)) +
  theme_bw()+theme(axis.line = element_line(colour = "black"),
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   panel.background = element_blank())

pca$Group = paste(pca$Disease,"|", pca$Cohort,sep = "")
Fig1a.taxa.pca

Fig1a.taxa.pc1.density <-
  ggplot(pca) +
  geom_density(aes(x=PC1, group=Group, fill=Disease, linetype=Cohort),
               color="black", alpha=0.6, position = 'identity') +
  scale_fill_manual(values=c("#f0999f","#46bbc0")) +
  theme_bw() +
  scale_linetype_manual(values = c("solid","dashed"))+
  labs(fill="")
Fig1a.taxa.pc1.density

Fig1a.taxa.pc2.density <-
  ggplot(pca) +
  geom_density(aes(x=PC2, group=Group, fill=Disease, linetype=Cohort),
               color="black", alpha=0.6, position = 'identity') +
  scale_fill_manual(values=c("#f0999f","#46bbc0")) +
  theme_bw() +
  scale_linetype_manual(values = c("solid","dashed"))+
  labs(fill="") + 
  coord_flip()
Fig1a.taxa.pc2.density

# Figure 1a. function #######

pca <- dat.funct
Fig1a.function.pca <- ggplot(pca,aes(PC1,PC2))+
  geom_point(size=2, aes(col=Disease,shape=Cohort))+
  scale_color_manual(values=c("#f0999f","#46bbc0"))+
  scale_shape_manual(values=c(16,15)) +
  theme_bw()+theme(axis.line = element_line(colour = "black"),
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   panel.background = element_blank())

pca$Group = paste(pca$Disease,"|", pca$Cohort,sep = "")
Fig1a.function.pca

Fig1a.function.pc1.density <-
  ggplot(pca) +
  geom_density(aes(x=PC1, group=Group, fill=Disease, linetype=Cohort),
               color="black", alpha=0.6, position = 'identity') +
  scale_fill_manual(values=c("#f0999f","#46bbc0")) +
  theme_bw() +
  scale_linetype_manual(values = c("solid","dashed"))+
  labs(fill="")
Fig1a.function.pc1.density

Fig1a.function.pc2.density <-
  ggplot(pca) +
  geom_density(aes(x=PC2, group=Group, fill=Disease, linetype=Cohort),
               color="black", alpha=0.6, position = 'identity') +
  scale_fill_manual(values=c("#f0999f","#46bbc0")) +
  theme_bw() +
  scale_linetype_manual(values = c("solid","dashed"))+
  labs(fill="") + 
  coord_flip()
Fig1a.function.pc2.density



# Figure 1b. --------------------------------------------
rm(list = ls())
library(reshape2)
dat <- read_excel("Fig 1 Source Data.xlsx", sheet ="Fig 1b species abund")
data2<-melt(dat,id.vals=c("NAME","Disease")) 
species=factor(data2$variable,
               levels=c("Moraxella_catarrhalis","Pseudomonas_aeruginosa","Haemophilus_parahaemolyticus","Stenotrophomonas_maltophilia","Streptococcus_intermedius","Acinetobacter_johnsonii","Neisseria_subflava","Prevotella_melaninogenica","Prevotella_intermedia","Haemophilus_parainfluenzae","Mogibacterium_diversum","Fusobacterium_pseudoperio.","Fusobacterium_nucleatum","Neisseria_flavescens","Porphyromonas_gingivalis","Veillonella_parvula","Lactobacillus_oris","Prevotella_scopos","Neisseria_meningitidis","Campylobacter_concisus","Prevotella_sp_oral_taxon_299","Neisseria_elongata","Filifactor_alocis","Prevotella_denticola","Parvimonas_micra","Treponema_sp_OMZ_838","Tannerella_sp_oral_taxon_286","Prevotella_fusca","Treponema_denticola","Tannerella_forsythia","Capnocytophaga_gingivalis"))
Fig1b <- ggplot(data2,aes(x=species,y=data2$value),colour=factor(data2$Disease))+
  geom_boxplot(aes(fill=data2$Disease),outlier.colour=NULL,outlier.shape=21,outlier.size=2)+
  ylim(0,0.5)+
  scale_fill_manual(values=c("#f0999f","#46bbc0"))+
  xlab("") + ylab("Normalized relative abundance") +
  theme_bw()+theme(axis.line = element_line(colour = "black"),
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   #  panel.border = element_blank(),
                   panel.background = element_blank(),
                   axis.text.x=element_text(angle=90, hjust=1, vjust=0.3, color="black",size=11),
                   axis.text.y=element_text(color="black",size=11),
                   axis.title.x=element_text(),
                   axis.title.y=element_text())

Fig1b
