# the required input file: 
# 1) Proteomic feature abundance (sputum_proteome.txt):An abundance table with proteins in rows and samples in columns.
# 2) Meta data prepared in txt file (metadata.txt), containing a column named SampleID and a column indicating disease state.

# output files generated:
# "3_significant_HostP_features.RData" 

library(data.table)
library(dplyr)

HostP.data <- data.frame(fread("sputum_proteome.txt"), row.names = 1 )

m = HostP.data
feature.abb_df <- cbind.data.frame(feature = rownames(m),
                                   abb = paste("feature",seq(1,nrow(m),1),sep = ""),
                                   stringsAsFactors = F)
rownames(m) <- sapply(rownames(m), function(x) feature.abb_df$abb[which(feature.abb_df$feature == x)])
m <- t(m) %>% as.data.frame(stringsAsFactors=F)


diseaseState = 'COPD'
meta_df <- fread("metadata.txt",data.table = F)
colnames(meta_df)[colnames(meta_df) == diseaseState ] <- "Y"


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
  if(tmp$`Pr(>|z|)`[which(rownames(tmp) == "mod")] <= 0.05)  sig.modules <- append(sig.modules, md)
}

if(all(grepl("^feature", colnames(m)))) sig.modules <- sapply(sig.modules, function(x) feature.abb_df$feature[which(feature.abb_df$abb == x)])

HostP.sigFeatures <- sig.modules

save(HostP.sigFeatures, file = "2_significant_hostP_features.RData")