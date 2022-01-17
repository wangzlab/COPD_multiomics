# the required input file: 
# 1) Metagenomic feature abundance (metagenome.gct): An abundance table in gctx format with KEGG KO ids in rows and samples in columns.
# 2) Metadata prepared in txt file, containing a column named SampleID and variables of confounders and a column indicating disease state.

# output files generated:
# "2_significant_metaG_modules.RData" 


library(data.table)
library(dplyr)

source("functions/io.R")


writeLines("1.Loading MetaG data")
MetaG.feature.dat <- parse.gctx(fname="metagenome.gct")@mat %>% as.data.frame()
m = MetaG.feature.dat
feature.abb_df <- cbind.data.frame(feature = rownames(m),
                                   abb = paste("feature",seq(1,nrow(m),1),sep = ""),
                                   stringsAsFactors = F)
rownames(m) <- sapply(rownames(m), function(x) feature.abb_df$abb[which(feature.abb_df$feature == x)])
m <- t(m) %>% as.data.frame(stringsAsFactors=F)

meta_df <- fread("source.data/metadata.txt",data.table = F)
disease.state = "COPD"

writeLines("2.Calculating residuals of glm models by MetaG features and co-founders")
Residuals_res <- NULL
FC.res_df <- NULL
for(md in colnames(m)){
  # md = colnames(m)[1]
  i_report = which(colnames(m) == md)
  if(i_report%%1000 == 0) writeLines(paste("Progress:",i_report," out of ",ncol(m),' metaG features', sep = ""))
  
  m.t_sub <- m %>% dplyr::select(all_of(md))
  # if(all(m.t_sub == 0)) next
  dat <- merge(meta_df[complete.cases(meta_df),], 
               m.t_sub,
               by.x="SampleID", by.y=0) %>% # remove all=T
    tibble::column_to_rownames("SampleID")
  #colnames(dat)[which(colnames(dat) == md)] <- "mod"
  #print(md)
  
  colnames(dat)
  glm.md <- glm(as.formula(paste(md, " ~ ", 
                                 paste(colnames(dat)[!colnames(dat) %in% c(md, disease.state)], collapse = " + "  ),
                                 sep = "" )),
                data = dat, family = gaussian)
  
  Residuals_res_c <- as.data.frame(resid(glm.md))    
  colnames(Residuals_res_c)[ncol(Residuals_res_c)] <- md
  Residuals_res <- bind_cols(Residuals_res, Residuals_res_c)
  
  if(length(unique(dat[,disease.state])) == 2){
    # calculate fold change of residuals
    fc_dat <- merge(Residuals_res_c, meta_df %>% select(SampleID, all_of(disease.state)),
                    by.x=0, by.y="SampleID")
    colnames(fc_dat)[colnames(fc_dat) == md] <- "feature"
    colnames(fc_dat)[colnames(fc_dat) == disease.state] <- "disease"
    
    fc_res <- fc_dat %>% group_by(disease) %>% summarise(avg = mean(feature)) 
    FC <- fc_res$avg[2]-fc_res$avg[1]
    from = fc_res$disease[1]
    to = fc_res$disease[2]
    FC.res_df_c <- cbind.data.frame("feature" = md,
                                    "disease.from" = from, 
                                    "disease.to" = to,
                                    "fold.change" = FC)
    FC.res_df <- bind_rows(FC.res_df, FC.res_df_c)
  }else{
    FC.res_df <- NULL
    writeLines("Note: can not calculate fold change because disease state has more than 2 levels. ")
  }
  
} # loop through modules (or features)

if(all(grepl("^feature", colnames(m)))){
  colnames(Residuals_res) <-
    sapply(colnames(Residuals_res), function(x) feature.abb_df$feature[which(feature.abb_df$abb == x)])
}  
if(!is.null(FC.res_df)){
  if(all(grepl("^feature", FC.res_df$feature))){
    FC.res_df$feature <-
      sapply(FC.res_df$feature, function(x) feature.abb_df$feature[which(feature.abb_df$abb == x)] )
  }
}

writeLines("3.Calculating wilcoxon rank test p for each MetaG module")
# calculate wilcoxon rank test p for each MetaG module  ----------
## import gene set databases
gene.set.databases = "source.data/KEGG_modules.gmt"
GSDB <- vector('list', length(gene.set.databases))
names(GSDB) <- gene.set.databases

for (gsdb in gene.set.databases)
  GSDB[[gsdb]] <- Read.GeneSets.db2(gsdb, thres.min = 2, thres.max = 2000)

for(i in 1:length(GSDB)){
  if(i == 1){
    gs <- GSDB[[i]]$gs
    N.gs <- GSDB[[i]]$N.gs
    gs.names <- GSDB[[i]]$gs.names
    gs.descs <- GSDB[[i]]$gs.desc
    size.G <-  GSDB[[i]]$size.G
  } else {
    gs <- append(gs, GSDB[[i]]$gs)
    gs.names <- append(gs.names, GSDB[[i]]$gs.names)
    gs.descs <- append(gs.descs, GSDB[[i]]$gs.desc)
    size.G <- append(size.G, GSDB[[i]]$size.G)
    N.gs <- N.gs + GSDB[[i]]$N.gs
  }
}

## calculate wiocoxon rank test p in modules
Module_pvalue <- NULL
for(i in c(1:length(gs))){
  #i=1
  if(i%%100 == 0) writeLines(paste("progress: ",i ," out of ", length(gs)," modules" ))
  md = gs.names[i]
  ftrs = gs[[i]]
  othr.ftrs = colnames(m)[!colnames(m) %in% ftrs]
  
  if(any(FC.res_df$feature %in% ftrs)){
    FC.res_df$inModule <-
      sapply(FC.res_df$feature, function(x) if(x %in% ftrs) "in" else "other")
  }else next  
  
  w = wilcox.test(fold.change~inModule, data = FC.res_df)
  vec_c = c("module" = unname(md), "Wilcox.p"=w$p.value) 
  Module_pvalue <- bind_rows(Module_pvalue, vec_c)
}

MetaG.sigMods <-  Module_pvalue$module[Module_pvalue$Wilcox.p < 0.1,]
save(MetaG.sigMods, file = "2_significant_metaG_modules.RData")

# quantification of MetaG modules were obtained from ssGSEA analysis,
# see https://github.com/broadinstitute/ssGSEA2.0 for more information