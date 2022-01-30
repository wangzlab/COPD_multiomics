library(readxl)
library(ggplot2)
library(ggrepel)

excel_sheets("Fig E3 data.xlsx")

# Figure E3a ------------------
dat <- read_excel("Fig E3 data.xlsx", sheet = "Metabolome")

FigE3a.metab <- ggplot(dat) +
  geom_point(aes(x=NEU, y=EOS), size=2)+
  theme_bw()+ theme(panel.grid = element_blank()) +
  geom_hline(yintercept = 0, linetype="twodash") +
  geom_hline(yintercept = 1, linetype="dashed") +
  geom_hline(yintercept = -1, linetype="dashed") +
  geom_vline(xintercept = -1, linetype="dashed") +
  geom_vline(xintercept = 0, linetype="twodash") +
  geom_vline(xintercept = 1, linetype="dashed") +
  xlab("Correlation NEU: Directionality x - log(P)")+
  ylab("Correlation EOS: Directionality x - log(P)")+
  ggtitle("Differential MetaB modules")
FigE3a.metab


# Figure E3b ------------------
dat <- read_excel("Fig E3 data.xlsx", sheet = "Transcriptome")

FigE3b.hostT <- ggplot(dat) +
  geom_point(aes(x=NEU, y=EOS), size=2)+
  theme_bw()+ theme(panel.grid = element_blank()) +
  geom_hline(yintercept = 0, linetype="twodash") +
  geom_hline(yintercept = 1, linetype="dashed") +
  geom_hline(yintercept = -1, linetype="dashed") +
  geom_vline(xintercept = -1, linetype="dashed") +
  geom_vline(xintercept = 0, linetype="twodash") +
  geom_vline(xintercept = 1, linetype="dashed") +
  xlab("Correlation NEU: Directionality x - log(P)")+
  ylab("Correlation EOS: Directionality x - log(P)")+
  ggtitle("Differential HostT modules")
FigE3b.hostT


# Figure E3c ------------------
dat <- read_excel("Fig E3 data.xlsx", sheet = "Sputum Proteome")

FigE3c.sputum <- ggplot(dat) +
  geom_point(aes(x=NEU, y=EOS), size=2)+
  geom_text_repel(aes(x=NEU,y=EOS, label = Label),size=3) +
  theme_bw()+ theme(panel.grid = element_blank()) +
  geom_hline(yintercept = 0, linetype="twodash") +
  geom_hline(yintercept = 1, linetype="dashed") +
  geom_hline(yintercept = -1, linetype="dashed") +
  geom_vline(xintercept = -1, linetype="dashed") +
  geom_vline(xintercept = 0, linetype="twodash") +
  geom_vline(xintercept = 1, linetype="dashed") +
  xlab("Correlation NEU: Directionality x - log(P)")+
  ylab("Correlation EOS: Directionality x - log(P)")+
  ggtitle("Differential sputum proteins")
FigE3c.sputum


# Figure E3d ------------------
dat <- read_excel("Fig E3 data.xlsx", sheet = "Serum Proteome")
colnames(dat)[1] <- "Label"

FigE3d.serum <- ggplot(dat) +
  geom_point(aes(x=NEU2, y=EOS), size=2)+
  geom_text_repel(aes(x=NEU2,y=EOS, label = Label),size=3) +
  theme_bw()+ theme(panel.grid = element_blank()) +
 # geom_hline(yintercept = 0, linetype="twodash") +
  geom_hline(yintercept = 1, linetype="dashed") +
  geom_hline(yintercept = -1, linetype="dashed") +
  geom_vline(xintercept = -1, linetype="dashed") +
 # geom_vline(xintercept = 0, linetype="twodash") +
  geom_vline(xintercept = 1, linetype="dashed") +
  xlab("Correlation NEU: Directionality x - log(P)")+
  ylab("Correlation EOS: Directionality x - log(P)")+
  ggtitle("Differential serum proteins")
FigE3d.serum

# integrate Figures
library(ggpubr)

FigE3 <- ggarrange(FigE3a.metab, FigE3b.hostT, FigE3c.sputum, FigE3d.serum)
  
