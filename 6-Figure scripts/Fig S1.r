list.files()

library(readxl)

sheets <- excel_sheets("Fig S1 data.xlsx")

for(s in sheets){
  dat <-
    read_excel("Fig S1 data.xlsx", sheet = s) %>%
    mutate(Color = cut(P, breaks=c(0,0.05,1), labels=c("red3","gray")))
  colnames(dat)[1] <- "Var"
  
  library(ggplot2)
  
  dat$Var <- factor(dat$Var, levels = rev(dat$Var))
  dat$Color <- as.character(dat$Color)
  p <- ggplot(dat) +
    geom_col(aes(x=R2, y=Var, fill=Color), width = 0.6) +
    scale_fill_manual(values = unique(dat$Color)[order(unique(dat$Color))]) +
    theme_classic() +  scale_x_continuous(position = "top") +
    xlab(s) + ylab("") + 
    theme(legend.position = "none")
  
  assign(paste("P_", s,sep = ""),p, envir = .GlobalEnv)
  
}

library(ggpubr)

ggarrange(ggarrange(P_Taxonomy, P_Function, P_Metabolome, ncol = 3),
          ggarrange(P_Transcriptome,`P_Sputum Proteome`, `P_Serum Proteome`,ncol = 3),
          nrow = 2)
