# Fig 4b ------------------------------------------
# Fig 4b. LOSO density plot #######
rm(list = ls())
dat.hist <- read_excel("Fig 4 data.xlsx", sheet = "Fig 4b LOSO")
colnames(dat.hist) <- c("Figs_map", "Figs_color","species","X1")
dat.hist <- dat.hist %>% mutate(X2=X1, Y1=-0.002, Y2=-0.009)

taxonomy_df <- read_excel("Fig 4 data.xlsx", sheet = "Fig 4b taxonomy")

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
Fig4b.LOSO <- ggarrange(histP_MetaG_P00380, histP_MetaG_M00044, histP_MetaG_P00250, histP_MetaG_P00564, ncol = 1)
#ggsave(Fig4b.LOSO,device = "pdf", filename = "Fig4b.KOcontri.pdf", width = 3, height = 10)


# Fig 4b. KO contribution #######

zscore_df <- read_excel("Fig 4 data.xlsx", sheet = "Fig 4b KO contrib" )
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


Fig4b.KOcontri <- ggplot(zscore.top3_df) +
  geom_point(aes(x=X, y=K, size=score, fill=phylum_other), shape=21) +
  scale_fill_manual(values = Colors)  + 
  scale_size (range = c (3, 6)) +
  theme(panel.grid = element_blank(), axis.title = element_blank(), panel.background = element_blank(),
        axis.text.x = element_blank(), legend.position = "none", axis.ticks = element_blank())

Fig4b.KOcontri
