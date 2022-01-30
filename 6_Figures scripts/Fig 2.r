list.files()

library(readxl)
excel_sheets("Fig 2 data.xlsx")

library(dplyr)

# Figure 2a ------------------------------------------

# Figure 2a. metabolome #######
rm(list = ls())
dat.metab <- read_excel("Fig 2 data.xlsx", sheet = "Fig 2a metabolome PCA")

library(ggplot2)
pca <- dat.metab
Fig2a.metab.pca <- ggplot(pca,aes(PC1,PC2))+
  geom_point(size=2,aes(col=Group,shape=Cohort))+
  scale_color_manual(values=c("#f0999f","#46bbc0"))+
  scale_shape_manual(values=c(16,15)) +
  theme_bw()+theme(axis.line = element_line(colour = "black"),
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   panel.background = element_blank())
Fig2a.metab.pca

pca$GROUP = paste(pca$Group,"|", pca$Cohort,sep = "")

Fig2a.metab.pc1.density <-
  ggplot(pca) +
  geom_density(aes(x=PC1, group=GROUP, fill=Group, linetype=Cohort),
               color="black", alpha=0.6, position = 'identity') +
  scale_fill_manual(values=c("#f0999f","#46bbc0")) +
  theme_bw() +
  scale_linetype_manual(values = c("solid","dashed"))+
  labs(fill="")
Fig2a.metab.pc1.density

Fig2a.metab.pc2.density <-
  ggplot(pca) +
  geom_density(aes(x=PC2, group=GROUP, fill=Group, linetype=Cohort),
               color="black", alpha=0.6, position = 'identity') +
  scale_fill_manual(values=c("#f0999f","#46bbc0")) +
  theme_bw() +
  scale_linetype_manual(values = c("solid","dashed"))+
  labs(fill="") + 
  coord_flip()
Fig2a.metab.pc2.density



# Figure 2a. metabolome #######
excel_sheets("Fig 2 data.xlsx")
rm(list = ls())
dat.trans <- read_excel("Fig 2 data.xlsx", sheet = "Fig 2a transcriptome PCA")

#library(ggplot2)
pca <- dat.trans
Fig2a.trans.pca <- ggplot(pca,aes(PC1,PC2))+
  geom_point(size=2,aes(col=Group,shape=Cohort))+ 
  scale_color_manual(values=c("#f0999f","#46bbc0"))+
  scale_shape_manual(values=c(16,15)) +
  theme_bw()+theme(axis.line = element_line(colour = "black"),
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   panel.background = element_blank())
Fig2a.trans.pca

pca$GROUP = paste(pca$Group,"|", pca$Cohort,sep = "")

Fig2a.trans.pc1.density <-
  ggplot(pca) +
  geom_density(aes(x=PC1, group=GROUP, fill=Group, linetype=Cohort),
               color="black", alpha=0.6, position = 'identity') +
  scale_fill_manual(values=c("#f0999f","#46bbc0")) +
  theme_bw() +
  scale_linetype_manual(values = c("solid","dashed"))+
  labs(fill="")
Fig2a.trans.pc1.density

Fig2a.trans.pc2.density <-
  ggplot(pca) +
  geom_density(aes(x=PC2, group=GROUP, fill=Group, linetype=Cohort),
               color="black", alpha=0.6, position = 'identity') +
  scale_fill_manual(values=c("#f0999f","#46bbc0")) +
  theme_bw() +
  scale_linetype_manual(values = c("solid","dashed"))+
  labs(fill="") + 
  coord_flip()
Fig2a.trans.pc2.density


# Figure 2a. sputum proteome #######
excel_sheets("Fig 2 data.xlsx")
rm(list = ls())
dat.sputum <- read_excel("Fig 2 data.xlsx", sheet = "Fig 2a sputum proteome PCA")

#library(ggplot2)
pca <- dat.sputum
Fig2a.sputum.pca <- ggplot(pca,aes(PC1,PC2))+
  geom_point(size=2,aes(col=Group))+
  scale_color_manual(values=c("#f0999f","#46bbc0"))+
  scale_shape_manual(values=c(16,15)) +
  theme_bw()+theme(axis.line = element_line(colour = "black"),
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   panel.background = element_blank())
Fig2a.sputum.pca

Fig2a.sputum.pc1.density <-
  ggplot(pca) +
  geom_density(aes(x=PC1, group=Group, fill=Group),
               color="black", alpha=0.6, position = 'identity') +
  scale_fill_manual(values=c("#f0999f","#46bbc0")) +
  theme_bw() +
  scale_linetype_manual(values = c("solid","dashed"))+
  labs(fill="")
Fig2a.sputum.pc1.density

Fig2a.sputum.pc2.density <-
  ggplot(pca) +
  geom_density(aes(x=PC2, group=Group, fill=Group),
               color="black", alpha=0.6, position = 'identity') +
  scale_fill_manual(values=c("#f0999f","#46bbc0")) +
  theme_bw() +
  scale_linetype_manual(values = c("solid","dashed"))+
  labs(fill="") + 
  coord_flip()
Fig2a.sputum.pc2.density


# Figure 2a. serum proteome #######
excel_sheets("Fig 2 data.xlsx")
rm(list = ls())
dat.serum <- read_excel("Fig 2 data.xlsx", sheet = "Fig 2a serum proteome PCA")

#library(ggplot2)
pca <- dat.serum
Fig2a.serum.pca <- ggplot(pca,aes(PC1,PC2))+
  geom_point(size=2,aes(col=Group))+  
  scale_color_manual(values=c("#f0999f","#46bbc0"))+
  scale_shape_manual(values=c(16,15)) +
  theme_bw()+theme(axis.line = element_line(colour = "black"),
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   panel.background = element_blank())
Fig2a.serum.pca

Fig2a.serum.pc1.density <-
  ggplot(pca) +
  geom_density(aes(x=PC1, group=Group, fill=Group),
               color="black", alpha=0.6, position = 'identity') +
  scale_fill_manual(values=c("#f0999f","#46bbc0")) +
  theme_bw() +
  scale_linetype_manual(values = c("solid","dashed"))+
  labs(fill="")
Fig2a.serum.pc1.density

Fig2a.serum.pc2.density <-
  ggplot(pca) +
  geom_density(aes(x=PC2, group=Group, fill=Group),
               color="black", alpha=0.6, position = 'identity') +
  scale_fill_manual(values=c("#f0999f","#46bbc0")) +
  theme_bw() +
  scale_linetype_manual(values = c("solid","dashed"))+
  labs(fill="") + 
  coord_flip()
Fig2a.serum.pc2.density


# Figure 2b ------------------------------------------
excel_sheets("Fig 2 data.xlsx")
rm(list = ls())
metab.modules <- read_excel("Fig 2 data.xlsx", sheet = "Fig 2b metab modules")
metab.features <- read_excel("Fig 2 data.xlsx", sheet = "Fig 2b metab features")

# Figure 2b.density of modules ####
data <- metab.modules
library(reshape2)
data2<-melt(data,id.vals=c("NAME","CLASS","GROUP") )
md.levels <- c("MEbrown4","MEmediumpurple1","MEskyblue3","MElavenderblush3","MEslateblue",
               "MEsalmon2","MEbisque4","MEplum","MEmidnightblue","MEindianred4")

data2$variable=factor(data2$variable,levels=md.levels)
data2$GROUP4 <- paste(data2$GROUP,data2$CLASS, sep = "|")
Fig2b.density <- ggplot(data2,aes(data2$value))+
  geom_density(aes(group=data2$GROUP4,fill=data2$CLASS,alpha=0.4,linetype=data2$GROUP))+
  scale_color_manual(values=c("#f0999f","#46bbc0"))+
  facet_wrap(~data2$variable,scale="free",ncol=1)+
  scale_linetype_manual(values=c("solid","dashed"))+
  scale_x_continuous(limits=c(-0.2,0.2)) +
  theme_bw()+theme(axis.line = element_line(colour = "black"),
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   panel.background = element_blank())
Fig2b.density
ggsave(Fig2b.density, filename = "Fig2b.density.pdf",device = "pdf", width = 4, height = 10)

# Figure 2b.heatmap ####
data <- metab.features 
data.Sz <- metab.features %>% select(starts_with("Z")) 
data.Gz <- metab.features %>% select(!starts_with("Z")) %>% select(-Metabolite, -Group)

colors<-c(seq(-4,0,length=51),seq(0.1,4,length=50))
my_palette<- colorRampPalette(c("#EE9494","#F6E572","#94EE9D","#94E1EE","#94B3EE","white","white","white","white","white","white","white"))(n=100)
library(gplots)
pdf( "Fig2b.metab.hm.Sz.pdf") #save figure2b.metab.heamap.Shenzhen in pdf
heatmap.2(as.matrix(data.Sz),
          breaks=colors, scale="row", trace="none", 
          Rowv=FALSE, dendrogram='column',hclustfun = function(x) hclust(x, method="ward.D2"),
          col=rev(my_palette))
dev.off()

pdf( "Fig2b.metab.hm.Gz.pdf") #save figure2b.metab.heamap.Guangzhou in pdf
heatmap.2(as.matrix(data.Gz),
          breaks=colors, scale="row", trace="none", 
          Rowv=FALSE, dendrogram='column',hclustfun = function(x) hclust(x, method="ward.D2"),
          col=rev(my_palette))  
dev.off()


## Figure 2c ------------------------------------------
excel_sheets("Fig 2 data.xlsx")
rm(list = ls())
trans.modules <- read_excel("Fig 2 data.xlsx", sheet = "Fig 2c trans modules")
trans.features <- read_excel("Fig 2 data.xlsx", sheet = "Fig 2c trans features")

# Figure 2c.density of modules ####
data <- trans.modules
#library(reshape2)
data2<-melt(data,id.vals=c("NAME","CLASS","GROUP") )
md.levels <- c("MEcyan","MEdarkmagenta.1",
               "MEorangered4","MEskyblue2",
               "MEblueviolet","MEthistle4",
			   "MElightyellow","MEchocolate4",
			   "MEbisque4.1","MEsnow4")

data2$variable=factor(data2$variable,levels=md.levels)
data2$GROUP4 <- paste(data2$GROUP,data2$CLASS, sep = "|")
Fig2c.density <- ggplot(data2,aes(value))+
  geom_density(aes(group=GROUP4,fill=CLASS,alpha=0.4,linetype=GROUP))+
  scale_color_manual(values=c("#f0999f","#46bbc0"))+
  facet_wrap(~variable,scale="free",ncol=1)+
  scale_linetype_manual(values=c("solid","dashed"))+
  scale_x_continuous(limits=c(-0.2,0.2)) +
  theme_bw()+theme(axis.line = element_line(colour = "black"),
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   panel.background = element_blank())
Fig2c.density
ggsave(Fig2c.density, filename = "Fig2c.density.pdf",device = "pdf", width = 4, height = 10)


# Figure 2c.heatmap ####
data <- trans.features 
data.Sz <- trans.features %>% select(starts_with("Z")) 
data.Gz <- trans.features %>% select(!starts_with("Z")) %>% select(-Gene, -Module)

colors <- c(seq(-4, 0, length=51), seq(0.1, 4, length=50))
my_palette<- colorRampPalette(c("#EE9494","#F6E572","#94EE9D","#94E1EE","#94B3EE","white","white","white","white","white","white","white"))(n=100)
#library(gplots)
#pdf( "Fig2c.trans.hm.Sz.pdf") #save figure2c.trans.heamap.Shenzhen in pdf
heatmap.2(as.matrix(data.Sz),
          breaks=colors, scale="row", trace="none", 
          Rowv=FALSE, dendrogram='column',hclustfun = function(x) hclust(x, method="ward.D2"),
          col=rev(my_palette))
#dev.off()

pdf( "Fig2c.trans.hm.Gz.pdf") #save figure2c.trans.heamap.Guangzhou in pdf
heatmap.2(as.matrix(data.Gz),
          breaks=colors, scale="row", trace="none", hclustfun = function(x) hclust(x, method="ward.D2"),
          Rowv=FALSE, dendrogram='column',
          col=rev(my_palette))
dev.off()


## Figure 2d ------------------------------------------
excel_sheets("Fig 2 data.xlsx")
rm(list = ls())
prot.features <- read_excel("Fig 2 data.xlsx", sheet = "Fig 2d proteome features")
sputum.p <- prot.features %>% select(SampleID, Group, starts_with("Sputum"))
serum.p <- prot.features %>% select(SampleID, Group, starts_with("serum"))

# Figure2d. sputum proteome ####
data <- sputum.p

data.l <- reshape2::melt(data, id.vars=c("SampleID", "Group"))
sapply(data.l,class)
data.l_positive <- data.l %>% filter(value > 0)

H_median_df <- data.l_positive %>% filter(Group == "H") %>% 
  group_by(variable) %>% summarise(H_median =median(value)) %>%
  as.data.frame(stringsAsFactors=F)

plotDat <- 
  merge(data.l_positive,H_median_df,by="variable") %>%
  mutate(normalized_value = value/H_median) %>% 
  mutate(normalized.log = log10(normalized_value))


Fig2d.sputum <- ggplot(plotDat,aes(normalized.log))+
  geom_density(aes(group=Group,fill=Group,alpha=0.4))+
  scale_color_manual(values=c("#f0999f","#46bbc0"))+
  facet_wrap(~variable,scale="free",ncol = 1) +
  # scale_linetype_manual(values=c("solid","dashed"))+
  #scale_x_continuous(limits=c(-0.2,0.2))+
  theme_bw()+theme(axis.line = element_line(colour = "black"),
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   panel.background = element_blank())

Fig2d.sputum 
ggsave(Fig2d.sputum, filename = 'Fig2d.sputum.pdf', device = "pdf",width = 4,height = 10)



# Figure2d. serum proteome ####
data <- serum.p

data.l <- reshape2::melt(data, id.vars=c("SampleID", "Group"))
sapply(data.l,class)
data.l_positive <- data.l %>% filter(value > 0)

H_median_df <- data.l_positive %>% filter(Group == "H") %>% 
  group_by(variable) %>% summarise(H_median =median(value)) %>%
  as.data.frame(stringsAsFactors=F)

plotDat <- 
  merge(data.l_positive,H_median_df,by="variable") %>%
  mutate(normalized_value = value/H_median) %>% 
  mutate(normalized.log = log10(normalized_value))


Fig2d.serum <- ggplot(plotDat,aes(normalized.log))+
  geom_density(aes(group=Group,fill=Group,alpha=0.4))+
  scale_color_manual(values=c("#f0999f","#46bbc0"))+
  facet_wrap(~variable,scale="free",ncol = 1) +
  # scale_linetype_manual(values=c("solid","dashed"))+
  scale_x_continuous(limits=c(-0.2,0.2))+
  theme_bw()+theme(axis.line = element_line(colour = "black"),
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   panel.background = element_blank())

Fig2d.serum 
ggsave(Fig2d.serum, filename = 'Fig2d.serum.pdf', device = "pdf",width = 4,height = 4)
