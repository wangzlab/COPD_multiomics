library(readxl)
library(stringr)

sheets <- excel_sheets("Fig S3 Source Data.xlsx")

diseaseState <- c("COPD","Health")
omics <- c("taxonomy","metagenome","metabolome","transcriptome","spu proteome","ser proteome")


for(ds in diseaseState){
  # ds=diseaseState[2]
  
  for(o in omics){
    # o = omics[2]
    s = sheets[grepl(ds, sheets) & grepl(o, sheets)]
    
    dat <- read_excel("Fig S3 Source Data.xlsx", sheet = s)
    
    for(Grp in c("Smoking", "ICS")){
      # Grp = "Smoking"
      dat.tmp <- dat 
      if(!Grp %in% colnames(dat.tmp)) next
      
      if(Grp == "Smoking") Colors <- c("#E6194B","#3CB44B") else Colors <-c("#1ACDE5","#D8BA64")
      
      colnames(dat.tmp)[colnames(dat.tmp) == Grp] <- "Grp"
      
      dat.tmp$Grp <- as.factor(dat.tmp$Grp)
      p<-ggplot(dat.tmp,aes(x=PC1,y=PC2))+
        geom_point(aes(color = Grp)) +
        scale_color_manual(values = Colors)+
        scale_fill_manual(values = Colors)+
        theme_bw()+ theme(panel.grid = element_blank(),
                          legend.position = c(0.9, 0.85),
                          legend.title = element_blank(),
                          axis.text.x = element_blank(),
                          axis.text.y = element_blank(),
                          axis.ticks = element_blank())+
        ggtitle(paste(str_to_title(o), ds, Grp,sep = " "))
      p
      assign(paste("P_",paste(c(o,ds,Grp),collapse = "."),sep = ""), p, envir = .GlobalEnv)
    }
  }
}


library(ggpubr)
FigS3 <- 
  ggarrange(
    ggarrange(P_taxonomy.COPD.Smoking,P_taxonomy.Health.Smoking, P_taxonomy.COPD.ICS, ncol = 3),
    ggarrange(P_metagenome.COPD.Smoking,P_metabolome.Health.Smoking, P_metagenome.COPD.ICS, ncol = 3),
    ggarrange(P_metabolome.COPD.Smoking,P_metabolome.Health.Smoking,P_metabolome.COPD.ICS,ncol = 3),
    ggarrange(P_transcriptome.COPD.Smoking,P_transcriptome.Health.Smoking,P_transcriptome.COPD.ICS,ncol = 3),
    ggarrange(`P_spu proteome.COPD.Smoking`,`P_spu proteome.Health.Smoking`,`P_spu proteome.COPD.ICS`,ncol = 3),
    ggarrange(`P_ser proteome.COPD.Smoking`,`P_ser proteome.Health.Smoking`,`P_ser proteome.COPD.ICS`,ncol = 3),
    nrow = 6)
FigS3