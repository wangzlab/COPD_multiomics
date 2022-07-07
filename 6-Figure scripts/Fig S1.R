library(readxl)
library(dplyr)



#Figure S1a --------------------------------------
# Figure S1a. genes ####
dat.genes <- read_excel("Fig S1 Source Data.xlsx", sheet = "Fig S1a rarefaction genes")
dat.KO <- read_excel("Fig S1 Source Data.xlsx", sheet = "Fig S1a rarefaction KOs")


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
FigS1a.genes <- ggdraw(aligned_plots[[1]]) + draw_plot(aligned_plots[[2]])
FigS1a.genes


# Figure S1a. KO ####

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
FigS1a.KO <- ggdraw(aligned_plots[[1]]) + draw_plot(aligned_plots[[2]])
FigS1a.KO



# Figure S1b. -----------------------------------------------------
rm(list = ls())
# Fig S1b.stats #####
data<-read_excel("Fig S1 Source Data.xlsx" ,sheet =  "Fig S1b bin stats")
plotDat<-melt(data,id.vals="NAME")
group<-factor(plotDat$Group,levels=c("High","Medium","Low"))
FigS1b.stats <- 
  ggplot(plotDat,aes(x=group,y=value))+
  geom_boxplot(outlier.shape=NA,aes(fill=group))+
  geom_jitter(alpha=0.1,shape=20,size=3,position=position_jitter(0.2))+
  facet_wrap(~variable,scale="free")+
  scale_fill_manual(values=c("#3081bc","#9fc8e2","#dfebf5")) +
  theme_bw()+theme(axis.line = element_line(colour = "black"),
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   panel.background = element_blank())
FigS1b.stats

# Fig S1b.taxa #####
data <- read_excel("Fig S1 Source Data.xlsx", sheet = "Fig S1b bin taxa")
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

FigS1b.taxa <- ggplot(plotDat) +
  geom_col(aes(x=variable, y=freq, fill=Genus))+
  scale_fill_manual(values = rev(mypalette))+
  theme_classic()
FigS1b.taxa