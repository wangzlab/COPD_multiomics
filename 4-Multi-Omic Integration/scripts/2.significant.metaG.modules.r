# the required input file: 
# 1) Metagenomic feature abundance (metagenome_rel_asinsqrt_zscore.gct):An abundance table in gctx format with KEGG KO ids in rows and samples in columns.
# 2) Meta data prepared in txt file, containing a column named SampleID and variables of co-founders and a column indicating disease state.
# 3) Metagenomic module abundance ("metagenome.gct")
# 4) "meta.mediation.NEU.txt": meta data containing NEU 
# 5) "meta.mediation.EOS.txt": meta data containing EOS

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

## ################################################################################################
# calculate correlation between modules and NEU or EOS,
# find significantly correlated modules and organize into MetaG.sigMods.NEU and MetaG.sigMods.EOS
## ################################################################################################

# load MetaG module abundance  ----------------------------
input.ds <- "1_DimReduction/metagenome.gct"
gct.unique <- NULL
dataset <- try(parse.gctx(input.ds), silent = T)
if(class(dataset) != 'try-error' ){
  
  m <- dataset@mat
  gene.names <- dataset@rid
  gene.descs <- dataset@rdesc
  sample.names <- dataset@cid
  sample.descs <- dataset@cdesc
  
} else {
  
  ## - cmapR functions stop if ids are not unique
  ## - import gct using readLines and make ids unique
  if(length(grep('rid must be unique', dataset) ) > 0) {
    gct.tmp <- readLines(input.ds)
    #first column
    rid <- gct.tmp %>% sub('\t.*','', .)
    #gct version
    ver <- rid[1]
    #data and meta data columns
    meta <- strsplit(gct.tmp[2], '\t') %>% unlist() %>% as.numeric()
    if(ver=='#1.3')
      rid.idx <- (meta[4]+3) : length(rid)
    else
      rid.idx <- 4:length(rid)
    
    #check whether ids are unique
    if(length(rid[rid.idx]) > length(unique(rid[rid.idx]))){
      warning('rids not unique! Making ids unique and exporting new GCT file...\n\n')
      #make unique
      rid[rid.idx] <- make.unique(rid[rid.idx], sep='_')
      #other columns
      rest <- gct.tmp %>% sub('.*?\t','', .)
      rest[1] <- ''
      gct.tmp2 <- paste(rid, rest, sep='\t') 
      gct.tmp2[1] <-  sub('\t.*','',gct.tmp2[1])
      
      #export
      gct.unique <- sub('\\.gct', '_unique.gct', input.ds)
      writeLines(gct.tmp2, con=gct.unique)
      
      #import using cmapR functions
      dataset <- parse.gctx(fname = gct.unique)
      
      ## extract data 
      m <- dataset@mat
      gene.names <- sub('_[0-9]{1,5}$', '', dataset@rid)
      gene.descs <- dataset@rdesc
      sample.names <- dataset@cid
      sample.descs <- dataset@cdesc
    }
    
  } else { #end if 'rid not unique'
    
    ########################################################
    ## display a more detailed error message if the import 
    ## failed due to other reasons than redundant 'rid'
    stop("\n\nError importing GCT file using 'cmapR::parse.gctx()'. Possible reasons:\n\n1) Please check whether you have the latest version of the 'cmapR' installed. Due to submission to Bioconductor the cmap team changed some naming conventions, e.g 'parse.gctx()' has been renamed to 'parse.gctx()'.\n2) The GCT file doesn't seem to be in the correct format! Please see take a look at https://clue.io/connectopedia/gct_format for details about GCT format.
\nError message returned by 'cmapR::parse.gctx()':\n\n", dataset, '\n\n')
  } 
}
MetaG.Mod.dat <- m %>% t() %>% as.data.frame()


# load meta data --------------------------------------
meta.NEU <- fread("meta.mediation.NEU.txt")
meta.EOS <- fread("meta.mediation.EOS.txt")


# merge data and calculate correlation between module and inflammatory feature ------------------------
meta <- merge(meta.NEU, meta.EOS, by="SampleID",all = T) %>% select(SampleID, NEU, EOS)
# sampleID 不匹配，meta数据里面sampleID加上X
meta$SampleID[grepl("^\\d", meta$SampleID)] <- paste("X",meta$SampleID[grepl("^\\d", meta$SampleID)],sep = "")

rp.res <- NULL
for(df.name in c("MetaG.Mod.dat")){
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


MetaG.sigMods.NEU <- 
  (rp.res %>% filter(grepl("MetaG", Omic)) %>% filter(ClinicVar == "NEU") %>% filter(p <= 0.1))$ModuleName
MetaG.sigMods.NEU <- intersect(MetaG.sigMods.NEU, MetaG.sigMods)

MetaG.sigMods.EOS <- 
  (rp.res %>% filter(grepl("MetaG", Omic)) %>% filter(ClinicVar == "EOS") %>% filter(p <= 0.1))$ModuleName
MetaG.sigMods.EOS <- intersect(MetaG.sigMods.EOS, MetaG.sigMods)



save(MetaG.sigMods.NEU, MetaG.sigMods.EOS, file = "2_significant_metaG_modules.RData")

# quantification of MetaG modules were obtained from ssGSEA analysis,
# see https://github.com/broadinstitute/ssGSEA2.0 for more information