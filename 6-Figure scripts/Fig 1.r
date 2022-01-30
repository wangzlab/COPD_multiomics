list.files()

library(readxl)
library(dplyr)

excel_sheets("Fig 1 data.xlsx")


#Figure 1a --------------------------------------
# Figure 1a. genes ####
dat.genes <- read_excel("Fig 1 data.xlsx", sheet = "Fig 1a rarefaction genes")
dat.KO <- read_excel("Fig 1 data.xlsx", sheet = "Fig 1a rarefaction KOs")


data1<-dat.genes[(dat.genes$Group=="COPD"),]
data2<-dat.genes[(dat.genes$Group=="Health"),]
dat_combined <- merge(data1,data2, by = "Samples",all = T)


library(ggplot2)
p1<-ggplot(dat_combined,aes(x=factor(Samples),y=Genes.x/1000000))+
  geom_boxplot(colour="#f0999f",outlier.colour=NULL,outlier.shape=19,outlier.size=1)+
  #scale_y_continuous(limits = c(0,2))  +
  theme_bw() + theme(panel.grid.minor =element_blank(), panel.grid.major = element_blank()) 
p1

library(ggpubr)
library(cowplot)
p2<-ggplot(dat_combined,aes(x=factor(Samples),y=Genes.y/1000000)) + 
  geom_boxplot(colour="#46bbc0",outlier.colour=NULL,outlier.shape=19,outlier.size=1)+
  scale_y_continuous(limits = c(0,1.25), position = "right",
                     breaks = seq(0,1.25, by=0.25))  +
  theme_half_open(11, rel_small = 1) + 
  rremove("x.axis")+
  rremove("xlab") +
  rremove("x.text") +
  rremove("x.ticks") +
  rremove("legend") 

p2

aligned_plots <- align_plots(p1, p2, align="hv", axis="tblr")
Fig1a.genes <- ggdraw(aligned_plots[[1]]) + draw_plot(aligned_plots[[2]])
Fig1a.genes


# Figure 1a. KO ####

data1<-dat.KO[(dat.KO$Group=="COPD"),]
data2<-dat.KO[(dat.KO$Group=="Health"),]
dat_combined <- merge(data1,data2, by = "Samples",all = T)

p1<-ggplot(dat_combined,aes(x=factor(Samples),y=Genes.x/1000))+
  geom_boxplot(colour="#f0999f",outlier.colour=NULL,outlier.shape=19,outlier.size=1)+
  scale_y_continuous(limits = c(2.5,7))  +
  theme_bw() + theme(panel.grid.minor =element_blank(), panel.grid.major = element_blank()) 
p1


p2<-ggplot(dat_combined,aes(x=factor(Samples),y=Genes.y/1000)) + 
  geom_boxplot(colour="#46bbc0",outlier.colour=NULL,outlier.shape=19,outlier.size=1)+
  scale_y_continuous(limits =  c(2.5,7), position = "right",
                     breaks = seq(3,7, by=1))  +
  theme_half_open(11, rel_small = 1) + 
  rremove("x.axis")+
  rremove("xlab") +
  rremove("x.text") +
  rremove("x.ticks") +
  rremove("legend") 

p2

aligned_plots <- align_plots(p1, p2, align="hv", axis="tblr")
Fig1a.KO <- ggdraw(aligned_plots[[1]]) + draw_plot(aligned_plots[[2]])
Fig1a.KO

# Figure 1b ------------------------------------------
# Figure 1b. taxonomy #######
rm(list = ls())
dat.taxa <- read_excel("Fig 1 data.xlsx", sheet = "Fig 1b taxonomy PCA")
dat.funct <- read_excel("Fig 1 data.xlsx", sheet = "Fig 1b metagenome PCA")

pca <- dat.taxa
Fig1b.taxa.pca <- ggplot(pca,aes(PC1,PC2))+
  geom_point(size=2,aes(col=Disease,shape=Cohort))+ 
  scale_color_manual(values=c("#f0999f","#46bbc0"))+
  scale_shape_manual(values=c(16,15)) +
  theme_bw()+theme(axis.line = element_line(colour = "black"),
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   panel.background = element_blank())

pca$Group = paste(pca$Disease,"|", pca$Cohort,sep = "")
Fig1b.taxa.pca

Fig1b.taxa.pc1.density <-
  ggplot(pca) +
  geom_density(aes(x=PC1, group=Group, fill=Disease, linetype=Cohort),
                 color="black", alpha=0.6, position = 'identity') +
  scale_fill_manual(values=c("#f0999f","#46bbc0")) +
  theme_bw() +
  scale_linetype_manual(values = c("solid","dashed"))+
  labs(fill="")
Fig1b.taxa.pc1.density

Fig1b.taxa.pc2.density <-
  ggplot(pca) +
  geom_density(aes(x=PC2, group=Group, fill=Disease, linetype=Cohort),
               color="black", alpha=0.6, position = 'identity') +
  scale_fill_manual(values=c("#f0999f","#46bbc0")) +
  theme_bw() +
  scale_linetype_manual(values = c("solid","dashed"))+
  labs(fill="") + 
  coord_flip()
Fig1b.taxa.pc2.density

# Figure 1b. function #######

pca <- dat.funct
Fig1b.function.pca <- ggplot(pca,aes(PC1,PC2))+
  geom_point(size=2, aes(col=Disease,shape=Cohort))+
  scale_color_manual(values=c("#f0999f","#46bbc0"))+
  scale_shape_manual(values=c(16,15)) +
  theme_bw()+theme(axis.line = element_line(colour = "black"),
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   panel.background = element_blank())

pca$Group = paste(pca$Disease,"|", pca$Cohort,sep = "")
Fig1b.function.pca

Fig1b.function.pc1.density <-
  ggplot(pca) +
  geom_density(aes(x=PC1, group=Group, fill=Disease, linetype=Cohort),
               color="black", alpha=0.6, position = 'identity') +
  scale_fill_manual(values=c("#f0999f","#46bbc0")) +
  theme_bw() +
  scale_linetype_manual(values = c("solid","dashed"))+
  labs(fill="")
Fig1b.function.pc1.density

Fig1b.function.pc2.density <-
  ggplot(pca) +
  geom_density(aes(x=PC2, group=Group, fill=Disease, linetype=Cohort),
               color="black", alpha=0.6, position = 'identity') +
  scale_fill_manual(values=c("#f0999f","#46bbc0")) +
  theme_bw() +
  scale_linetype_manual(values = c("solid","dashed"))+
  labs(fill="") + 
  coord_flip()
Fig1b.function.pc2.density


# Figure 1c. --------------------------------------------
rm(list = ls())
library(reshape2)
dat <- read_excel("Fig 1 data.xlsx", sheet ="Fig 1c species abund")
data2<-melt(dat,id.vals=c("NAME","Disease")) 
species=factor(data2$variable,
               levels=c("Moraxella_catarrhalis","Pseudomonas_aeruginosa","Haemophilus_parahaemolyticus","Stenotrophomonas_maltophilia","Streptococcus_intermedius","Acinetobacter_johnsonii","Neisseria_subflava","Prevotella_melaninogenica","Prevotella_intermedia","Haemophilus_parainfluenzae","Mogibacterium_diversum","Fusobacterium_pseudoperio.","Fusobacterium_nucleatum","Neisseria_flavescens","Porphyromonas_gingivalis","Veillonella_parvula","Lactobacillus_oris","Prevotella_scopos","Neisseria_meningitidis","Campylobacter_concisus","Prevotella_sp_oral_taxon_299","Neisseria_elongata","Filifactor_alocis","Prevotella_denticola","Parvimonas_micra","Treponema_sp_OMZ_838","Tannerella_sp_oral_taxon_286","Prevotella_fusca","Treponema_denticola","Tannerella_forsythia","Capnocytophaga_gingivalis"))
Fig1c <- ggplot(data2,aes(x=species,y=data2$value),colour=factor(data2$Disease))+
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

Fig1c


# Figure 1d. --------------------------------------------
rm(list = ls())
library(dplyr)
dat <-read_excel("Fig 1 data.xlsx", sheet ="Fig 1d KEGG module")
meta <- dat %>% tibble::column_to_rownames("...1" ) %>% 
  t() %>% as.data.frame() %>%
  tibble::rownames_to_column("sample") %>% dplyr::select(sample, Group)
 
plotDat <- dat %>% filter(`...1` != "Group")
plotDat[,-1] <- sapply(plotDat[,-1], as.numeric)
plotDat <- plotDat %>% tibble::column_to_rownames("...1")

plotDat.Sz <- plotDat %>% dplyr::select(dplyr::starts_with("Z"))
plotDat.Gz <- plotDat %>% dplyr::select(!dplyr::starts_with("Z"))


library(gplots)
colors<-c(seq(-4,0,length=51), seq(0.1,4,length=50))
my_palette<-colorRampPalette(c("#EE9494","#F6E572","#94EE9D","#94E1EE","#94B3EE","white","white","white","white","white","white"))(n=100)
#map<-heatmap.2(as.matrix(data), breaks=colors, scale="none", trace="none", margins=c(5,7), Rowv=TRUE, dendrogram='both', col=rev(my_palette), distfun = function(x) dist(x, method="euclidean"), hclustfun = function(x) hclust(x, method="ward.D2"))
map.Gz <-heatmap.2(as.matrix(plotDat.Gz), 
                 breaks=colors, scale="row", trace="none", margins=c(5,7),
                 dendrogram='column', 
                 col=rev(my_palette),
                 distfun = function(x) dist(x, method="euclidean"), 
                 hclustfun = function(x) hclust(x, method="ward.D2"))

rev(rownames(plotDat.Gz)[map.Gz$rowInd])

plotDat.Sz<-plotDat.Sz[match(rev(rownames(plotDat.Gz)[map.Gz$rowInd]), rownames(plotDat.Sz)),]
map.Sz<-heatmap.2(as.matrix(plotDat.Sz), 
               breaks=colors, scale="row", trace="none", margins=c(5,7),
               Rowv = F,dendrogram='column', 
               col=rev(my_palette),
               distfun = function(x) dist(x, method="euclidean"), 
               hclustfun = function(x) hclust(x, method="ward.D2"))
map.Sz

# Figure 1e. -----------------------------------------------------
rm(list = ls())
# Fig1e.stats #####
data<-read_excel("Fig 1 data.xlsx" ,sheet =  "Fig 1e bin stats")
plotDat<-melt(data,id.vals="NAME")
group<-factor(plotDat$Group,levels=c("High","Medium","Low"))
Fig1e.stats <- 
  ggplot(plotDat,aes(x=group,y=value))+
  geom_boxplot(outlier.shape=NA,aes(fill=group))+
  geom_jitter(alpha=0.1,shape=20,size=3,position=position_jitter(0.2))+
  facet_wrap(~variable,scale="free")+
  scale_fill_manual(values=c("#3081bc","#9fc8e2","#dfebf5")) +
  theme_bw()+theme(axis.line = element_line(colour = "black"),
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   panel.background = element_blank())
Fig1e.stats

# Fig1e.taxa #####
data <- read_excel("Fig 1 data.xlsx", sheet = "Fig 1e bin taxa")
plotDat <- data %>% reshape2::melt() %>% 
  group_by(variable, Genus) %>% 
  summarise(sumed = sum(value), Phylum=unique(Phylum)) %>%
  mutate(freq=sumed/sum(sumed))

plotDat$Phylum <- factor(plotDat$Phylum, 
                         levels = c("Bacteroidetes","Proteobacteria","Actinobacteria",
                                    "Firmicutes","Spirochaetes","TM7"))
genus.level <- c("Prevotella","Porphyromonadaceae","Bacteroidales",
                 "Neisseria","Haemophilus","Moraxella","Campylobacter","Ralstonia","Pseudomonas","Lautropia",
                 "Rothia","Schaalia","Actinomyces",
                 "Streptococcus","Veillonella","Lachnospiraceae","Mogibacterium",
                 "Megasphaera","Selenomonas","Lactobacillus",
                 "Treponema","Saccharibacteria","Others","Unclassified")
plotDat$Genus[!plotDat$Genus  %in% genus.level] <- "Others"

plotDat$Genus <- factor(plotDat$Genus,levels = rev(genus.level))

mypalette <- c("#E25554","#E88D8D","#F2C3C3",
               "#EEA020","#F0AD4D","#F4BF76","#F7CF9E","#F6D9B9","#FBE9D8","#FCF2EA",
               "#49A9CE","#80CCE0","#BFE4EE",
               "#33A863","#4EB676","#77C392","#96CFA9",
               "#B4DCC1","#C5E3D0","#E2F1E6",
               "#D0BB8A","#968CBC","#B7B6B6","#E5E3E2")

Fig1e.taxa <- ggplot(plotDat) +
  geom_col(aes(x=variable, y=freq, fill=Genus))+
  scale_fill_manual(values = rev(mypalette))+
  theme_classic()
Fig1e.taxa