library(data.table)
library(dplyr)
library(reshape2)
library(ggplot2)
library(ggdendro)
library(readxl)
library(ggpubr)

excel_sheets("Fig E4 data.xlsx")

# Figure E4a ---------------------------------
CorrDf_metab.metag_neu <- read_excel("Fig E4 data.xlsx", sheet = "MetaG-MetaB Full")
CorrDf_metab.trans_neu <- read_excel("Fig E4 data.xlsx", sheet = "Trans-MetaB")
CorrDf_spupro.trans_neu <- read_excel("Fig E4 data.xlsx", sheet = "Trans-Sputpro")
CorrDf_serpro.trans_neu <- read_excel("Fig E4 data.xlsx", sheet = "Trans-Serpro")

keyStr_df <- cbind.data.frame(dfnameStr = c("metag","metab","trans","serpro","spupro"),
                              colStr = c("MetaG","MetaB","Trans","Cyto","Cyto"),
                              stringsAsFactors=F)

# nodes and types
nodeXs <- c("metab","serpro","spupro")  #do not change the sequences of elements in nodeXs 
nodeYs <- c("metag","trans")  #do not change the sequences of elements in nodeYs 
ENtypes <- c("neu")

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




# Fig E4a. NEU integrated plot
ggarrange(ggarrange(HeatP_metab.trans_neu + 
                      theme(legend.position = "none", axis.text.x = element_blank(),axis.text.y = element_blank()) +
                      xlab("MetaB") + ylab("Trans"), 
                    HeatP_spupro.trans_neu + 
                      theme(legend.position = "none", axis.text.x = element_blank(),axis.text.y = element_blank()) +
                      xlab("Sputum") + ylab("Trans"),
                    HeatP_serpro.trans_neu + 
                      theme(legend.position = "none", axis.text.x  = element_blank(),axis.text.y = element_blank()) +
                      xlab("Serum") + ylab("Trans"),
                    nrow = 1, widths = c(0.55,0.3,0.15)),
          ggarrange(HeatP_metab.metag_neu + 
                      theme(legend.position = "none", axis.text.x = element_blank(),axis.text.y = element_blank()) +
                      xlab("MetaB") + ylab("MetaG"), 
                    ggplot() + geom_text(aes(x=0,y=0),label="NEU") + theme_dendro(), 
                    nrow = 1, widths = c(0.55,0.45)),
          nrow = 2,heights = c(0.6,0.4))




excel_sheets("Fig E4 data.xlsx")


# Figure E4c ----------------------------------
dat  <- read_excel("Fig E4 data.xlsx", sheet = "Prop of Mediation All")

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

FigE4c <- ggplot(dat.dt ,
       aes(x=Group, y=Mediation_Proportion))+
  geom_violin(aes(fill=Path),trim = T, scale = "width", alpha=0.5)+
  geom_boxplot(width=0.1, outlier.shape = NA) + 
  scale_fill_manual(values = c("#F8766D","#CD9600","#7CAE00","#00BE67")) +
  theme_bw()+ theme(panel.grid = element_blank(),
                    axis.text.x = element_blank())+
  coord_flip()
FigE4c


# Figure E4d ----------------------------------

dat <- read_excel("Fig E4 data.xlsx", sheet = "Reverse Mediation")

plotDat <- dat %>%
  mutate(Forward = ABC_Prop, Reverse = ACB_Prop) %>%
  reshape2::melt(id.vars=c("Comparison","Type")) %>% 
  dplyr::filter(variable %in% c("Forward", "Reverse")) 

plotDat$Type <- factor(plotDat$Type, 
                       levels = c("MetaG-MetaB-NEU","MetaB-HostT-NEU","MetaG-MetaB-EOS","MetaB-HostT-EOS"))

FigE4d <- ggpaired(plotDat, x = "variable", y = "value",
         color = "variable", line.color = "gray", line.size = 0.4)+
  stat_compare_means(paired = TRUE) +
  facet_wrap(vars(Type), scales = "free", ncol = 4)+
  theme_bw()+theme(panel.grid = element_blank())
FigE4d