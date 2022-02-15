# the required input file: 
# 1) Metabolomic feature abundance (metaB.module_eigengene.txt):An abundance table with compounds in rows and samples in columns.
# 2) Meta data prepared in txt file (meta.txt), containing a column named SampleID and a column indicating disease state.
# 3) "meta.mediation.NEU.txt": meta data containing NEU 
# 3) "meta.mediation.EOS.txt": meta data containing EOS

# output files generated:
# "3_significant_MetaB_modules.RData" 

library(data.table)
library(dplyr)

# load MetaB module abundance  ----------------------------

MetaB.Mod.dat <- fread("1_DimReduction/metaB.module_eigengene.txt",data.table = F) 
MetaB.Mod.dat[-1] <- sapply(MetaB.Mod.dat[-1], as.numeric)
MetaB.Mod.dat <- MetaB.Mod.dat %>% tibble::column_to_rownames("V1")

m = MetaB.Mod.dat
feature.abb_df <- cbind.data.frame(feature = rownames(m),
                                   abb = paste("feature",seq(1,nrow(m),1),sep = ""),
                                   stringsAsFactors = F)
rownames(m) <- sapply(rownames(m), function(x) feature.abb_df$abb[which(feature.abb_df$feature == x)])
m <- t(m) %>% as.data.frame(stringsAsFactors=F)


diseaseState = 'COPD'
meta_df <- fread("meta.txt",data.table = F)
colnames(meta_df)[colnames(meta_df) == diseaseState ] <- "Y"
#meta_df$Y <- as.factor(meta_df$Y)

#all(meta_df$SampleID %in% rownames(m.t))

#m.t <- t(m) %>% as.data.frame()
sig.modules <- vector("character")
for(md in colnames(m)){
  m.t_sub <- m %>% dplyr::select(all_of(md))
  # if(all(m.t_sub == 0)) next
  dat <- merge(meta_df, m.t_sub, by.x="SampleID", by.y=0,all = T) %>% dplyr::select(-SampleID)
  colnames(dat)[which(colnames(dat) == md)] <- "mod"
  #print(md)
  
  colnames(dat)
  glm.md <- glm(as.formula(paste("Y ~ ", paste(colnames(dat)[colnames(dat) != "Y"], collapse = " + "  ), sep = "" )),
                data = dat, family = "binomial")
  
  tmp <- summary(glm.md)$coefficients %>% as.data.frame()
  if( !("mod" %in% rownames(tmp)) ) next
  if(tmp$`Pr(>|z|)`[which(rownames(tmp) == "mod")] <= 0.1)  sig.modules <- append(sig.modules, md)
}

if(all(grepl("^feature", colnames(m)))) sig.modules <- sapply(sig.modules, function(x) feature.abb_df$feature[which(feature.abb_df$abb == x)])

MetaB.sigMods <- sig.modules

## ################################################################################################
# calculate correlation between modules and NEU or EOS,
# find significantly correlated modules and organize into MetaG.sigMods.NEU and MetaG.sigMods.EOS
## ################################################################################################



# load meta data --------------------------------------
meta.NEU <- fread("source.data/meta.mediation.NEU.txt")
meta.EOS <- fread("source.data/meta.mediation.EOS.txt")

# merge data and calculate correlation between module and inflammatory feature ------------------------
meta <- merge(meta.NEU, meta.EOS, by="SampleID",all = T) %>% select(SampleID, NEU, EOS)
# sampleID 不匹配，meta数据里面sampleID加上X
meta$SampleID[grepl("^\\d", meta$SampleID)] <- paste("X",meta$SampleID[grepl("^\\d", meta$SampleID)],sep = "")

MetaB.Mod.dat <- m

rp.res <- NULL
for(df.name in c("MetaB.Mod.dat")){
  # df.name = "MetaG.Mod.dat"
  
  omic.dat <- eval(parse(text = df.name))
  dat<-merge(omic.dat, meta, by.x=0, by.y="SampleID") %>% tibble::column_to_rownames("Row.names")
  
  for(clin.var in c("NEU","EOS")){
    # clin.var = "NEU"
    
    i.clin <- which(colnames(dat) == clin.var)
    
    for(i in c(1:ncol(omic.dat))){
      Mod.name = colnames(dat)[i]
      dat.sub <- dat[,c(i,i.clin)]
      dat.sub <- dat.sub[complete.cases(dat.sub),]
      Cor = cor.test(dat.sub[,1], dat.sub[,2], method = "spearman")
      r = Cor$estimate
      p = Cor$p.value
      
      vec <- c(df.name, Mod.name, clin.var, r, p)
      names(vec) <- c("Omic","ModuleName","ClinicVar","r","p")
      
      rp.res <- bind_rows(rp.res, vec)
    }# loop through Modules
  }# loop through NEU and EOS
}# loop through omic data


MetaB.sigMods.NEU <- 
  (rp.res %>%
     filter(grepl("MetaB", Omic)) %>% filter(ClinicVar == "NEU") %>%
     mutate(padj=p.adjust(p)) %>%
     filter(padj <= 0.1))$ModuleName
MetaB.sigMods.NEU <- intersect(MetaB.sigMods.NEU, MetaB.sigMods)

MetaB.sigMods.EOS <- 
  (rp.res %>% 
     filter(grepl("MetaB", Omic)) %>% filter(ClinicVar == "EOS") %>% 
     mutate(padj=p.adjust(p)) %>%
     filter(padj <= 0.1))$ModuleName
MetaB.sigMods.EOS <- intersect(MetaB.sigMods.EOS, MetaB.sigMods)



save(MetaB.sigMods.NEU, MetaB.sigMods.EOS, file = "3_significant_MetaB_modules.RData")
