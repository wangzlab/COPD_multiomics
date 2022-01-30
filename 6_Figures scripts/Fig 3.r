library(data.table)
library(dplyr)
library(reshape2)
library(ggplot2)
library(ggdendro)
library(readxl)
library(ggpubr)

list.files()

excel_sheets("Fig 3 data.xlsx")

#Figure 3a ---------------------------------
CorrDf_metab.metag_neu <- read_excel("Fig 3 data.xlsx", sheet = "Fig 3a metag_metab_neu")
CorrDf_metab.metag_eos <- read_excel("Fig 3 data.xlsx", sheet = "Fig 3a metag_metab_eos")

CorrDf_metab.trans_eos <- read_excel("Fig 3 data.xlsx", sheet = "Fig 3a metab_trans_eos")
CorrDf_metab.trans_neu <- read_excel("Fig 3 data.xlsx", sheet = "Fig 3a metab_trans_neu")

CorrDf_spupro.trans_eos <- read_excel("Fig 3 data.xlsx", sheet = "Fig 3a trans_spuprot_eos")
CorrDf_spupro.trans_neu <- read_excel("Fig 3 data.xlsx", sheet = "Fig 3a trans_spuprot_neu")

CorrDf_serpro.trans_eos <- read_excel("Fig 3 data.xlsx", sheet = "Fig 3a trans_serprot_eos")
CorrDf_serpro.trans_neu <- read_excel("Fig 3 data.xlsx", sheet = "Fig 3a trans_serprot_neu")

# only keep overlapped modules in each file
# NEU 
metaB.common <- intersect(unique(CorrDf_metab.metag_neu$Node1), unique(CorrDf_metab.trans_neu$Node2))
trans.common <- intersect(intersect(unique(CorrDf_metab.trans_neu$Node1), 
                                    unique(CorrDf_serpro.trans_neu$Node2)),
                          unique(CorrDf_spupro.trans_neu$Node2))

CorrDf_metab.metag_neu <- CorrDf_metab.metag_neu %>% filter(Node1 %in% metaB.common)
CorrDf_metab.trans_neu <- CorrDf_metab.trans_neu %>% filter(Node2 %in% metaB.common)

CorrDf_metab.trans_neu <- CorrDf_metab.trans_neu %>% filter(Node1 %in% trans.common)
CorrDf_serpro.trans_neu <- CorrDf_serpro.trans_neu %>% filter(Node2 %in% trans.common)
CorrDf_spupro.trans_neu <- CorrDf_spupro.trans_neu %>% filter(Node2 %in% trans.common)

# EOS
metaB.common <- intersect(unique(CorrDf_metab.metag_eos$Node1), unique(CorrDf_metab.trans_eos$Node2))
trans.common <- intersect(intersect(unique(CorrDf_metab.trans_eos$Node1), 
                                    unique(CorrDf_serpro.trans_eos$Node2)),
                          unique(CorrDf_spupro.trans_eos$Node2))

CorrDf_metab.metag_eos <- CorrDf_metab.metag_eos %>% filter(Node1 %in% metaB.common)
CorrDf_metab.trans_eos <- CorrDf_metab.trans_eos %>% filter(Node2 %in% metaB.common)

CorrDf_metab.trans_eos <- CorrDf_metab.trans_eos %>% filter(Node1 %in% trans.common)
CorrDf_serpro.trans_eos <- CorrDf_serpro.trans_eos %>% filter(Node2 %in% trans.common)
CorrDf_spupro.trans_eos <- CorrDf_spupro.trans_eos %>% filter(Node2 %in% trans.common)

keyStr_df <- cbind.data.frame(dfnameStr = c("metag","metab","trans","serpro","spupro"),
                              colStr = c("MetaG","MetaB","Trans","Cyto","Cyto"),
                              stringsAsFactors=F)
# nodes and types
nodeXs <- c("metab","serpro","spupro")  #do not change the sequences of elements in nodeXs 
nodeYs <- c("metag","trans")  #do not change the sequences of elements in nodeYs 
ENtypes <- c("eos","neu")


# create heatmap from combinations of figures :


for(enType in ENtypes){
  #  enType = ENtypes[2]
  
  for(nodex in nodeXs){
    # nodex = nodeXs[1]
    
    for(nodey in nodeYs){
      # nodey = nodeYs[2]
      
      
      corrDfName <- paste("CorrDf_",nodex,".",nodey,"_",enType,sep = "")
      if(!exists(corrDfName)) next
      
      CorrDf <- eval(parse(text = corrDfName))
      
      NodeXCol <- which(sapply(CorrDf, function(x) all(grepl(keyStr_df$colStr[which(keyStr_df$dfnameStr == nodex)], x) )) )
      colnames(CorrDf)[NodeXCol] <- "NodeX"
      
      NodeYCol <- which(sapply(CorrDf, function(x) all(grepl(keyStr_df$colStr[which(keyStr_df$dfnameStr == nodey)], x) )) )
      colnames(CorrDf)[NodeYCol] <- "NodeY"
      
      
      # organize the orders of nodex and nodey
      dat_r.w <- CorrDf %>% reshape2::dcast(NodeY ~ NodeX, value.var = "Correlation")
      rownames(dat_r.w) <- dat_r.w$NodeY; dat_r.w <- dat_r.w[-1]
      
      if(T){
        df <- t(dat_r.w) 
        x <- as.matrix(scale(df))
        dd.col <- as.dendrogram(hclust(dist(x)))
        col.ord <- order.dendrogram(dd.col)
        
        dd.row <- as.dendrogram(hclust(dist(t(x))))
        row.ord <- order.dendrogram(dd.row)
        
        xx <- scale(df)[col.ord, row.ord] 
        xx_names <- attr(xx, "dimnames") 
        #df <- as.data.frame(xx)
        ddata_x <- dendro_data(dd.row) 
        ddata_y <- dendro_data(dd.col) 
      } 
      if(!exists(paste("order_",nodey,sep = ""))) assign(paste("order_",nodey,sep = ""), xx_names[[2]],envir = .GlobalEnv)
      if(!exists(paste("order_",nodex,sep = ""))) assign(paste("order_",nodex,sep = ""), xx_names[[1]],envir = .GlobalEnv)
      
      
      order_nodex <- eval(parse(text = paste("order_",nodex,sep = "")))
      order_nodey <- eval(parse(text = paste("order_",nodey,sep = "")))
      
      if(!all(unique(CorrDf$NodeX) %in% order_nodex) ) {print(paste("not all nodes of ", nodex," were in predefined order so stop",sep = ""));break}
      if(!all(unique(CorrDf$NodeY) %in% order_nodey) ) {print(paste("not all nodes of ", nodey," were in predefined order so stop",sep = ""));break}
      
      
      CorrDf$NodeX <- factor(CorrDf$NodeX, levels = order_nodex)
      CorrDf$NodeY <- factor(CorrDf$NodeY, levels = order_nodey)
      
      CorrDf$ColorType <- sapply(c(1:nrow(CorrDf)), 
                                 function(i) {
                                   if(CorrDf$Linked[i] == "Y") return("Y") else if(CorrDf$`P-value`[i] >= 0.05) return("N_ns") else if(CorrDf$Correlation[i]>0) return("N_sig_posCorr") else return("N_sig_negCorr")
                                 })
      
      CorrDf$absCorr <- abs(CorrDf$Correlation)
      
      
      CorrDf <- CorrDf %>% filter( NodeX %in% order_nodex) %>% filter(NodeY %in% order_nodey)
      CorrDf$NodeX <- factor(CorrDf$NodeX, levels = order_nodex)
      CorrDf$NodeY <- factor(CorrDf$NodeY, levels = order_nodey)
      
      
      p <- ggplot(data = CorrDf, aes(x=NodeX,y=NodeY))+
        geom_tile(aes(fill=ColorType, alpha=absCorr),color="white") +
        theme(axis.text.x = element_text(angle = 90))+
        #scale_fill_manual(values=c("white","#e2e2e2","#cc0202"))+
        scale_fill_manual(values = c("white","#c1c1ff","#ffb6b6","#cc0202")) +
        #scale_alpha(limits = c(0.0,1.0), range = c(0,0.6))+
        theme(  panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.background = element_blank())
               # axis.title.x = element_text(colour=NA),
               # axis.title.y = element_blank())
      
      assign(paste("HeatP_",nodex,".",nodey,"_",enType,sep = ""), p, envir = .GlobalEnv)
      
    }
    
    
  }
  
  
  assign(paste("order_metab_",enType,sep = ""), order_metab, envir = .GlobalEnv)
  assign(paste("order_metag_",enType,sep = ""), order_metag, envir = .GlobalEnv)
  assign(paste("order_trans_",enType,sep = ""), order_trans, envir = .GlobalEnv)
  assign(paste("order_serpro_",enType,sep = ""), order_serpro, envir = .GlobalEnv)
  assign(paste("order_spupro_",enType,sep = ""), order_spupro, envir = .GlobalEnv)
  
  remove(order_metab,order_metag,order_trans, order_serpro, order_spupro)
}




# Fig 3a. NEU integrated plot 
ggarrange(ggarrange(HeatP_metab.trans_neu + 
                      theme(legend.position = "none", axis.text.x = element_blank(),axis.text.y = element_blank()) +
                      xlab("MetaB") + ylab("Trans"), 
                    HeatP_spupro.trans_neu + 
                      theme(legend.position = "none", axis.text.x = element_blank(),axis.text.y = element_blank()) +
                      xlab("Sputum") + ylab("Trans"),
                    HeatP_serpro.trans_neu + 
                      theme(legend.position = "none", axis.text.x  = element_blank(),axis.text.y = element_blank()) +
                      xlab("Serum") + ylab("Trans"),
                    nrow = 1, widths = c(0.5,0.3,0.2)),
          ggarrange(HeatP_metab.metag_neu + 
                      theme(legend.position = "none", axis.text.x = element_blank(),axis.text.y = element_blank()) +
                      xlab("MetaB") + ylab("MetaG"), 
                    ggplot() + geom_text(aes(x=0,y=0),label="NEU") + theme_dendro(), 
                    nrow = 1, widths = c(0.5,0.5)),
          nrow = 2)


# Fig 3a. EOS integrated plot
ggarrange(ggarrange(HeatP_spupro.trans_eos + 
                      theme(legend.position = "none", axis.text.x = element_blank(),axis.text.y = element_blank()) +
                      xlab("Sputum") + ylab("Trans"),
                    HeatP_serpro.trans_eos + 
                      theme(legend.position = "none", axis.text.x  = element_blank(),axis.text.y = element_blank()) +
                      xlab("Serum") + ylab("Trans"),
                    HeatP_metab.trans_eos + 
                      theme(legend.position = "none", axis.text.x = element_blank(),axis.text.y = element_blank()) +
                      xlab("MetaB") + ylab("Trans"), 
                    nrow = 1, widths = c(0.45,0.2,0.35)),
          ggarrange(ggplot() + geom_text(aes(x=0,y=0),label="EOS") + theme_dendro(),
                    HeatP_metab.metag_eos + 
                      theme(legend.position = "none", axis.text.x = element_blank(),axis.text.y = element_blank()) +
                      xlab("MetaB") + ylab("MetaG"), 
                    nrow = 1, widths = c(0.65,0.35)),
          nrow = 2)

# Fig 3c ------------------------------------------
# Fig 3c. LOSO density plot #######
rm(list = ls())
dat.hist <- read_excel("Fig 3 data.xlsx", sheet = "Fig 3c LOSO")
colnames(dat.hist) <- c("Figs_map", "Figs_color","species","X1")
dat.hist <- dat.hist %>% mutate(X2=X1, Y1=-0.002, Y2=-0.009)

taxonomy_df <- read_excel("Fig 3 data.xlsx", sheet = "Fig 3c taxonomy")

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
Fig3c.LOSO <- ggarrange(histP_MetaG_P00380, histP_MetaG_M00044, histP_MetaG_P00250, histP_MetaG_P00564, ncol = 1)
#ggsave(Fig3c.LOSO,device = "pdf", filename = "Fig3c.KOcontri.pdf", width = 3, height = 10)


# Fig 3c. KO contribution #######

zscore_df <- read_excel("Fig 3 data.xlsx", sheet = "Fig 3c KO contrib" )
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


Fig3c.KOcontri <- ggplot(zscore.top3_df) +
  geom_point(aes(x=X, y=K, size=score, fill=phylum_other), shape=21) +
  scale_fill_manual(values = Colors)  + 
  scale_size (range = c (3, 6)) +
  theme(panel.grid = element_blank(), axis.title = element_blank(), panel.background = element_blank(),
        axis.text.x = element_blank(), legend.position = "none", axis.ticks = element_blank())

Fig3c.KOcontri
