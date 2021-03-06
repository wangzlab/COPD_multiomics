# Figure E2a. Total bacterial load ####
dat <- read_excel("Fig E2 Source Data.xlsx", sheet = "Bacterial load")

library(ggpubr)
FigE2a <- ggplot(dat) +
  geom_boxplot(aes(x=Group_Site, y=log10cpn, fill=Group),outlier.shape = NA, alpha=0.5) +
  geom_jitter(aes(x=Group_Site, y=log10cpn, fill=Group),shape=21,size=2, width = 0.2) +
  theme_bw() + theme(panel.grid = element_blank()) +
  ylab("log10 (16S copy number)")

FigE2a

wilcox.test(log10cpn~Group, data = dat %>% as.data.frame() %>% dplyr::filter(Site=="Guangzhou")) 
wilcox.test(log10cpn~Group, data = dat %>% as.data.frame() %>% dplyr::filter(Site=="Shenzhen")) 


# Figure E2.b --------------------
# Figure E2b. taxonomy #####
dat <- read_excel("Fig E2 Source Data.xlsx", sheet = "Taxonomy")

FigE2b.taxonomy <- ggplot(dat) +
  geom_point(aes(x=NUE,y=EOS), size=2)+
  geom_text_repel(aes(x=NUE,y=EOS, label = Label),size=3) +
  theme_bw()+ theme(panel.grid = element_blank()) +
  geom_hline(yintercept = 0, linetype="twodash") +
  geom_hline(yintercept = 1, linetype="dashed") +
  geom_vline(xintercept = -1, linetype="dashed") +
  geom_vline(xintercept = 0, linetype="twodash") +
  geom_vline(xintercept = 1, linetype="dashed") +
  xlab("Correlation NEU: Directionality x - log(P)")+
  ylab("Correlation EOS: Directionality x - log(P)")+
  ggtitle("Differential species-level taxa")
FigE2b.taxonomy

# Figure E2b. KEGG modules #####
dat <- read_excel("Fig E2 Source Data.xlsx", sheet = "Metagenome")

FigE2b.KEGG <- ggplot(dat) +
  geom_point(aes(x=NEU,y=EOS), size=2) +
  theme_bw()+ theme(panel.grid = element_blank()) +
  geom_hline(yintercept = 0, linetype="twodash") +
  geom_hline(yintercept = 1, linetype="dashed") +
  geom_vline(xintercept = -1, linetype="dashed") +
  geom_vline(xintercept = 0, linetype="twodash") +
  geom_vline(xintercept = 1, linetype="dashed") +
  xlab("Correlation NEU: Directionality x - log(P)")+
  ylab("Correlation EOS: Directionality x - log(P)")+
  ggtitle("Differential KEGG modules")
FigE2b.KEGG



# Figure E2c. correlation #####
dat <- read_excel("Fig E2 Source Data.xlsx", sheet = "Correlation")

FigE2c <- ggplot(dat,aes(x=bin_based_FC, y=reads_based_FC)) +
  geom_point( size=2 )+
  theme_classic() + theme(panel.grid = element_blank()) +
  geom_hline(yintercept = 0, linetype="dashed") +
  geom_vline(xintercept = 0, linetype="dashed") +
  geom_smooth(method = "lm", se = F, linetype="dashed", color="black", size=0.5) +
  xlim(c(-20, 20)) + ylim(c(-3, 10)) +
  scale_y_continuous(breaks = seq(-2, 10, 2))+
  xlab("FC (Reads-based)")+
  ylab("FC (Bins-based)")
FigE2c

