list.files()

library(readxl)
# Fig S3 --------------------------------
# NEU 
dat <- read_excel("Fig S3 data.xlsx", sheet = "NEU")

library(ggplot2)
library(ggrepel)

head(dat)
labeledDots <- c("MetaG_P00380-MetaB_M056-HostT_M318","MetaG_M00044-MetaB_M100-HostT_M490","MetaG_P00250-MetaB_M006-HostT_M490")
plotDat <- 
  dat %>%
  mutate(Label = sapply(`...1`,function(x) if(x %in% labeledDots) x else NA) ) %>%
  mutate(Color = sapply(Label, function(x) if(is.na(x)) "black" else "red"))
p1.NEU <- ggplot(plotDat, aes(x=RSQ, y=RMSE)) +
  geom_point(aes(color = Color, fill=Color), shape=21) +
  scale_color_manual(values = c("black","red"))+
  scale_fill_manual(values = c("white","red"))+
  geom_text_repel(aes(label=Label))+
  theme_bw()+ theme(panel.grid = element_blank(),
                    legend.position = "none") 


p2.NEU <- ggplot(plotDat, aes(x=RSQ, y=MAE)) +
  geom_point(aes(color = Color, fill=Color), shape=21) +
  scale_color_manual(values = c("black","red"))+
  scale_fill_manual(values = c("white","red"))+
  geom_text_repel(aes(label=Label))+
  theme_bw()+ theme(panel.grid = element_blank(),
                    legend.position = "none") 


p3.NEU <- ggplot(plotDat, aes(x=RMSE, y=MAE)) +
  geom_point(aes(color = Color, fill=Color), shape=21) +
  scale_color_manual(values = c("black","red"))+
  scale_fill_manual(values = c("white","red"))+
  geom_text_repel(aes(label=Label))+
  theme_bw()+ theme(panel.grid = element_blank(),
                    legend.position = "none") 
library(ggpubr)

P.NEU <- ggarrange(p1.NEU, p2.NEU, p3.NEU, ncol = 3)



# EOS 
dat <- read_excel("Fig S3 data.xlsx", sheet = "EOS")

library(ggplot2)
library(ggrepel)

head(dat)
labeledDots <- c("MetaG_P00564-MetaB_M123-HostT_M391")
plotDat <- 
  dat %>%
  mutate(Label = sapply(`...1`,function(x) if(x %in% labeledDots) x else NA) ) %>%
  mutate(Color = sapply(Label, function(x) if(is.na(x)) "black" else "red"))
p1.EOS <- ggplot(plotDat, aes(x=RSQ, y=RMSE)) +
  geom_point(aes(color = Color, fill=Color), shape=21) +
  scale_color_manual(values = c("black","red"))+
  scale_fill_manual(values = c("white","red"))+
  geom_text_repel(aes(label=Label))+
  theme_bw()+ theme(panel.grid = element_blank(),
                    legend.position = "none") 


p2.EOS <- ggplot(plotDat, aes(x=RSQ, y=MAE)) +
  geom_point(aes(color = Color, fill=Color), shape=21) +
  scale_color_manual(values = c("black","red"))+
  scale_fill_manual(values = c("white","red"))+
  geom_text_repel(aes(label=Label))+
  theme_bw()+ theme(panel.grid = element_blank(),
                    legend.position = "none") 


p3.EOS <- ggplot(plotDat, aes(x=RMSE, y=MAE)) +
  geom_point(aes(color = Color, fill=Color), shape=21) +
  scale_color_manual(values = c("black","red"))+
  scale_fill_manual(values = c("white","red"))+
  geom_text_repel(aes(label=Label))+
  theme_bw()+ theme(panel.grid = element_blank(),
                    legend.position = "none") 
library(ggpubr)

P.EOS <-ggarrange(p1.EOS, p2.EOS, p3.EOS, ncol = 3)
