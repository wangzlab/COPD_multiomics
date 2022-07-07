library(data.table)
library(dplyr)
library(reshape2)
library(ggplot2)
library(ggdendro)
library(readxl)
library(ggpubr)

list.files()

excel_sheets("Fig 3 Source Data.xlsx")

#Figure 3a ---------------------------------
CorrDf_all <- read_excel("Fig 3 Source Data.xlsx", sheet = "Fig 3 Source Data")
CorrDf_all


# only keep overlapped modules in each file
# NEU 
metaB.common <- intersect((CorrDf_all %>% filter(Mediation_type == 'MetaG-MetaB-NEU'))$Node2,
                          (CorrDf_all %>% filter(Mediation_type == 'MetaB-HostT-NEU'))$Node1)
trans.common <- intersect((CorrDf_all %>% filter(Mediation_type == 'MetaB-HostT-NEU'))$Node2,
                          (CorrDf_all %>% filter(Mediation_type == 'HostT-Spuprot-NEU'))$Node1)

CorrDf_NEU <-bind_rows(
  CorrDf_all %>% 
    filter(Mediation_type == 'MetaG-MetaB-NEU') %>%
    filter(Node2 %in% metaB.common),
  CorrDf_all %>% 
    filter(Mediation_type == 'MetaB-HostT-NEU') %>%
    filter(Node1 %in% metaB.common) %>%
    filter(Node2 %in% trans.common),
  CorrDf_all %>% 
    filter(Mediation_type == 'HostT-Spuprot-NEU') %>%
    filter(Node1 %in% trans.common)
)


# EOS
metaB.common <- intersect((CorrDf_all %>% filter(Mediation_type == 'MetaG-MetaB-EOS'))$Node2,
                          (CorrDf_all %>% filter(Mediation_type == 'MetaB-HostT-EOS'))$Node1)
trans.common <- intersect((CorrDf_all %>% filter(Mediation_type == 'MetaB-HostT-EOS'))$Node2,
                          (CorrDf_all %>% filter(Mediation_type == 'HostT-Spuprot-EOS'))$Node1)

CorrDf_EOS <- bind_rows(
  CorrDf_all %>% 
    filter(Mediation_type == 'MetaG-MetaB-EOS') %>%
    filter(Node2 %in% metaB.common),
  CorrDf_all %>% 
    filter(Mediation_type == 'MetaB-HostT-EOS') %>%
    filter(Node1 %in% metaB.common) %>%
    filter(Node2 %in% trans.common),
  CorrDf_all %>% 
    filter(Mediation_type == 'HostT-Spuprot-EOS') %>%
    filter(Node1 %in% trans.common))


# nodes and types
nodeXs <- c("MetaB","Spuprot")  #do not change the sequences of elements in nodeXs 
nodeYs <- c("MetaG","HostT")  #do not change the sequences of elements in nodeYs 
ENtypes <- c("EOS","NEU")


# create heatmap from combinations of figures :


for(enType in ENtypes){
  #  enType = ENtypes[2]
  
  CorrDf_tmp <- eval(parse(text = paste0("CorrDf_",enType)))
  
  for(nodex in nodeXs){
    # nodex = nodeXs[1]
    
    for(nodey in nodeYs){
      # nodey = nodeYs[2]
      
      CorrDf <- CorrDf_tmp %>% filter(grepl(nodex, Mediation_type)) %>% filter(grepl(nodey, Mediation_type))
      if(nrow(CorrDf) == 0) next
      
      mediatp <- CorrDf$Mediation_type %>% unique()
      
      NodeXCol <- which(strsplit(mediatp, "-", fixed = T)[[1]] == nodex)
      colnames(CorrDf)[NodeXCol] <- "NodeX"
      
      NodeYCol <- which(strsplit(mediatp, "-", fixed = T)[[1]] == nodey)
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
                                   if(CorrDf$Linked[i] == "Y") return("Y") else if(CorrDf$`Mediation P-val`[i] >= 0.05) return("N_ns") else if(CorrDf$Correlation[i]>0) return("N_sig_posCorr") else return("N_sig_negCorr")
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
  
  
  assign(paste("order_MetaB_",enType,sep = ""), order_MetaB, envir = .GlobalEnv)
  assign(paste("order_MetaG_",enType,sep = ""), order_MetaG, envir = .GlobalEnv)
  assign(paste("order_HostT_",enType,sep = ""), order_HostT, envir = .GlobalEnv)
  assign(paste("order_Spuprot_",enType,sep = ""), order_Spuprot, envir = .GlobalEnv)
  
  remove(order_MetaB,order_MetaG,order_HostT, order_Spuprot)
}




# Fig 3a. NEU integrated plot 
ggarrange(ggarrange(HeatP_MetaB.HostT_NEU + 
                      theme(legend.position = "none", axis.text.x = element_blank(),axis.text.y = element_blank()) +
                      xlab("MetaB") + ylab("Trans"), 
                    HeatP_Spuprot.HostT_NEU + 
                      theme(legend.position = "none", axis.text.x = element_blank(),axis.text.y = element_blank()) +
                      xlab("Sputum") + ylab("Trans"),
                    nrow = 1, widths = c(0.65,0.35)),
          ggarrange(HeatP_MetaB.MetaG_NEU + 
                      theme(legend.position = "none", axis.text.x = element_blank(),axis.text.y = element_blank()) +
                      xlab("MetaB") + ylab("MetaG"), 
                    ggplot() + geom_text(aes(x=0,y=0),label="NEU") + theme_dendro(), 
                    nrow = 1, widths = c(0.65,0.35)),
          nrow = 2)


# Fig 3a. EOS integrated plot
ggarrange(ggarrange(HeatP_Spuprot.HostT_EOS + 
                      theme(legend.position = "none", axis.text.x = element_blank(),axis.text.y = element_blank()) +
                      xlab("Sputum") + ylab("Trans"),
                    HeatP_MetaB.HostT_EOS + 
                      theme(legend.position = "none", axis.text.x = element_blank(),axis.text.y = element_blank()) +
                      xlab("MetaB") + ylab("Trans"), 
                    nrow = 1, widths = c(0.45,0.35)),
          ggarrange(ggplot() + geom_text(aes(x=0,y=0),label="EOS") + theme_dendro(),
                    HeatP_MetaB.MetaG_EOS + 
                      theme(legend.position = "none", axis.text.x = element_blank(),axis.text.y = element_blank()) +
                      xlab("MetaB") + ylab("MetaG"), 
                    nrow = 1, widths = c(0.45,0.35)),
          nrow = 2, heights = c(0.6,0.4))
