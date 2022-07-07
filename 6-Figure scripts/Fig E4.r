library(data.table)
library(dplyr)
library(reshape2)
library(ggplot2)
library(ggdendro)
library(readxl)
library(ggpubr)


# Figure E4b ----------------------------------
dat  <- read_excel("Fig E4 Source Data.xlsx", sheet = "E4b Prop of Mediation All")

# normalize
dat.dt <- dat 


dat.dt$InflammationType <- sapply(strsplit(dat.dt$Group,"_",fixed = T),"[[",1 )
dat.dt$Path <- sapply(dat.dt$Group,
                      function(x){
                        parts <- strsplit(x, "_", fixed = T)[[1]]
                        paste(parts[2:length(parts)], collapse = "-")
                      })

dat.dt$Path <- factor(dat.dt$Path, levels = c("MetaG-MetaB","MetaB-Trans","Trans-Spuprot","Trans-Serprot"))

dat.dt$InflammationType <- factor(dat.dt$InflammationType,levels = c("NEU","EOS"))
dat.dt <- dat.dt %>% arrange(Path) %>% arrange(InflammationType)

dat.dt$Group <- factor(dat.dt$Group, levels = rev(unique(dat.dt$Group)))

FigE4b <- ggplot(dat.dt ,
                 aes(x=Group, y=Mediation_Proportion))+
  geom_violin(aes(fill=Path),trim = T, scale = "width", alpha=0.5)+
  geom_boxplot(width=0.1, outlier.shape = NA) + 
  scale_fill_manual(values = c("#F8766D","#CD9600","#7CAE00","#00BE67")) +
  theme_bw()+ theme(panel.grid = element_blank(),
                    axis.text.x = element_blank())+
  coord_flip()
FigE4b



# Figure E4c ----------------------------------

dat <- read_excel("Fig E4 Source Data.xlsx", sheet = "E4c Reverse Mediation")

plotDat <- dat %>%
  mutate(Forward = ABC_Prop, Reverse = ACB_Prop) %>%
  reshape2::melt(id.vars=c("Comparison","Type")) %>% 
  dplyr::filter(variable %in% c("Forward", "Reverse")) 

plotDat$Type <- factor(plotDat$Type, 
                       levels = c("MetaG-MetaB-NEU","MetaB-HostT-NEU","MetaG-MetaB-EOS","MetaB-HostT-EOS"))

FigE4c <- ggpaired(plotDat, x = "variable", y = "value",
                   color = "variable", line.color = "gray", line.size = 0.4)+
  stat_compare_means(paired = TRUE) +
  facet_wrap(vars(Type), scales = "free", ncol = 4)+
  theme_bw()+theme(panel.grid = element_blank())
FigE4c



# Figure E4e  ----------------------------------
# Fig E4e. LOSO density plot #######
rm(list = ls())
dat.hist <- read_excel("Fig E4 Source Data.xlsx4", sheet = "E4e LOSO")
colnames(dat.hist) <- c("Figs_map", "Figs_color","species","X1")
dat.hist <- dat.hist %>% mutate(X2=X1, Y1=-0.002, Y2=-0.009)

taxonomy_df <- read_excel("Fig E4 Source Data.xlsx4", sheet = "E4e taxonomy")

genus <- sapply(strsplit(dat.hist$species,"_", fixed = T),"[[", 1)
genus[!genus %in% taxonomy_df$Genus]

phylum <- sapply(genus, 
                 function(x){
                   if(x %in% taxonomy_df$Genus) {
                     taxonomy_df$Phylum[which(taxonomy_df$Genus == x)[1]]
                   }else "Unclassified"
                 })
phylum_color_df <- cbind.data.frame(phylum=c("Bacteroidetes","Actinobacteria","TM7","Proteobacteria","Firmicutes","other"),
                                    colors=c("#EF5656","#47B3DA","#9A8FC3","#F7A415","#2BB065","#BABABA"),
                                    stringsAsFactors=F)
dat.hist$phylum <- phylum
dat.hist$phylum_other <- sapply(dat.hist$phylum, function(x)if(x %in% phylum_color_df$phylum) x else "other")

Fig.map_ymax_df <- cbind.data.frame(Fig.map = unique(dat.hist$Figs_map),
                                    # ymax=c(0.04,0.16,0.03,0.015,0.035,0.038),
                                    ymax = c(0.06,0.19, 0.06, 0.05),
                                    stringsAsFactors=F)


for(Fig in unique(dat.hist$Figs_map)){
  
  #Fig=unique(dat.hist$Figs_map)[1]
  
  ymax = Fig.map_ymax_df$ymax[which(Fig.map_ymax_df$Fig.map == Fig)]
  
  dat_sub <-  dat.hist %>%  filter(Figs_map %in% Fig)
  
  phylum_ABC <- unique(dat_sub$phylum_other)[order(unique(dat_sub$phylum_other))]
  phylum_rank <- c(phylum_ABC[phylum_ABC != "other"],"other")
  
  dat_sub$phylum_other <- factor(dat_sub$phylum_other, 
                                 levels = phylum_rank)
  Colors <- sapply(phylum_rank, function(x) phylum_color_df$colors[which(phylum_color_df$phylum ==x)])
  
  histP <- ggplot(dat_sub) +
    geom_density(aes(x=X1,y=(..count..)/sum(..count..)),fill = "gray") + ylab("Density") + 
    xlab("") +
    theme_bw()+theme(panel.grid = element_blank()) +
    ylim(-0.01,ymax) +
    geom_segment(aes(x = X1, y = Y1, xend = X2, yend = Y2, color=phylum_other )) +
    scale_color_manual(values = Colors) + 
    theme(legend.position = "none") + ggtitle(Fig) 
  
  
  histP
  
  assign(paste("histP_", Fig, sep = ""), histP, envir = .GlobalEnv)
}

histP_MetaG_M00044
histP_MetaG_P00250
histP_MetaG_P00380
histP_MetaG_P00564


library(ggpubr)
FigE4e.LOSO <- ggarrange(histP_MetaG_P00380, histP_MetaG_M00044, histP_MetaG_P00250, histP_MetaG_P00564, ncol = 1)
#ggsave(FigE4e.LOSO,device = "pdf", filename = "FigE4e.KOcontri.pdf", width = 3, height = 10)


# Fig E4e. KO contribution #######

zscore_df <- read_excel("Fig E4 Source Data.xlsx", sheet = "E4e KO contrib" )
colnames(zscore_df) <- c("K","sp","score")

tmp <- zscore.top3_df <- zscore_df %>%
  group_by(K) %>%
  top_n(n = 3, wt = score) %>% as.data.frame()
tmp <- tmp %>% arrange(desc(score)) %>% arrange(K)

zscore.top3_df <- tmp %>% #keep the top 3 with highest score
  mutate(X = rep(c("X1","X2","X3"), length(unique(zscore_df$K)))) %>%
  mutate(genus = sapply(strsplit(sp," ", fixed = T),"[[", 1)) %>%
  mutate(phylum = sapply(genus,
                         function(x){
                           if(x %in% taxonomy_df$Genus){
                             taxonomy_df$Phylum[which(taxonomy_df$Genus == x)[1]]
                           }else "Unclassified"
                         } ) ) %>%
  mutate(phylum_other = sapply(phylum, function(x) if(x %in% phylum_color_df$phylum) x else "other"))

zscore.top3_df$phylum_other <- factor(zscore.top3_df$phylum_other, levels = phylum_rank)


FigE4e.KOcontri <- ggplot(zscore.top3_df) +
  geom_point(aes(x=X, y=K, size=score, fill=phylum_other), shape=21) +
  scale_fill_manual(values = Colors)  + 
  scale_size (range = c (3, 6)) +
  theme(panel.grid = element_blank(), axis.title = element_blank(), panel.background = element_blank(),
        axis.text.x = element_blank(), legend.position = "none", axis.ticks = element_blank())

FigE4e.KOcontri
