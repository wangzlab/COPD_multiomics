library(readxl)

excel_sheets("Fig S5 data.xlsx")


# Figure S5 a ----------------------------
dat <- read_excel("Fig S5 data.xlsx", sheet = "Shannon")
head(dat)
library(ggplot2)

FigS5a <- ggplot(dat) +
  geom_boxplot(aes(x=Group,y=Shannon,fill=Group),alpha=0.5,outlier.shape = NA) +
  geom_jitter(aes(x=Group,y=Shannon,fill=Group),shape=21, width = 0.2,size=2) +
  theme_bw() + theme(panel.grid = element_blank())+
  xlab("") + ggtitle("Shannon index for \n K01426 contribution")

# wilcox
w=wilcox.test(Shannon~Group,data = dat)
w$p.value  #0.006


# Figure S5 b ----------------------------
dat <- read_excel("Fig S5 data.xlsx", sheet = "Contribution")
colnames(dat)[1] <- "Species"

library(RColorBrewer)
nb.cols = 17
mycolors <- colorRampPalette(brewer.pal(12, "Paired"))(nb.cols)
library(scales)
show_col(mycolors)


# rank
species_rank <- c(dat$Species[dat$Species != "Others"],"Others")


dat <- reshape2::melt(dat, id.var="Species")
dat$Species <- factor(dat$Species, levels = rev(species_rank))

FigS5b <- ggplot(dat) +
  geom_col(aes(x=variable,y=value,fill=Species))+
  scale_fill_manual(values = c("gray",rev(mycolors)))+
  theme_bw() + theme(panel.grid = element_blank())+
  xlab("") + ylab("RPKM")+ ggtitle("Contribution to K01426")



# Figure S5c ------------------------------------
dat <- read_excel("Fig S5 data.xlsx", sheet = "Differential contribution")
colnames(dat)[1] <- "Species"

library(ggrepel)
FigS5c <- ggplot(dat) +
  geom_point(aes(x=FC, y=log10P), size=2, alpha=0.5, fill="gray", shape=21)+
  theme_bw() + theme(panel.grid = element_blank())+
  geom_text_repel(aes(x=FC, y=log10P, label=Label)) +
  geom_hline(yintercept = 1.3, color="red", size=0.3, linetype="dashed")+
  geom_vline(xintercept = 0, color="darkgray", size=0.3, linetype="dashed")+
  xlab("Fold-change of contribution to K01426") + ylab("minus log10 P-value") + 
  ggtitle("Differential contribution to K01426")


library(ggpubr)
FigS5 <- ggarrange(FigS5a, FigS5b, FigS5c,ncol = 3)
FigS5
