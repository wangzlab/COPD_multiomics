# the required input file: 
# 1) Metabolomic feature abundance (1_metaB.module_eigengene.txt): An abundance table with compounds in rows and samples in columns.
# 2) Meta data prepared in txt file (metadata.txt), containing a column named SampleID and a column indicating disease state.

# output files generated:
# "3_significant_MetaB_modules.RData" 

library(data.table)
library(dplyr)

MetaB.Mod.dat <- fread("1_metaB.module_eigengene.txt",data.table = F) %>%
  dplyr::filter(!grepl("#",`#NAME`,fixed=T)) 
MetaB.Mod.dat[-1] <- sapply(MetaB.Mod.dat[-1], as.numeric)
MetaB.Mod.dat <- MetaB.Mod.dat %>% tibble::column_to_rownames("#NAME")
m = MetaB.Mod.dat
feature.abb_df <- cbind.data.frame(feature = rownames(m),
                                   abb = paste("feature",seq(1,nrow(m),1),sep = ""),
                                   stringsAsFactors = F)
rownames(m) <- sapply(rownames(m), function(x) feature.abb_df$abb[which(feature.abb_df$feature == x)])
m <- t(m) %>% as.data.frame(stringsAsFactors=F)


diseaseState = 'COPD'
meta_df <- fread("metadata.txt",data.table = F)
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

save(MetaB.sigMods, file = "2_significant_metaB_modules.RData")