# files required:
# 1) "1_MetaG_affects_NEU_through_MetaB.txt": resulting data frame generated from script "3.mediation_metaG.2.metaB.r"
# 2) "KEGG_modules.tab" storing information on module IDs and corresponding KOs
# 3) "KO2METABO.lists.RData" generated from script "6.KO2metabo.r"
# 4) "1_metaB.module_assign.txt" module assignment file of the metabolomic featrues, generated from script "1.WGCNA_metabolomics.r"
# 5) "metabolome.txt" metabolomic abundance table at the feature level
# 6) "metagenome.gct" metagenomic abundance table at the feature level
# 7) "Metabo.KEGGModule.match.txt" storing information on KEGG module with corresponding compound ID
#

# output:
# "4_MetaG.MetaB.modules.linked.txt" stored in the "Output" directory



mediation.res <- fread("1_MetaG_affects_NEU_through_MetaB.txt")
 
output.dir = "Output"
log.file = "MetaG.MetaB.link.log"

cat(paste("\n\n",as.character(Sys.time()), '\n'),  file=log.file, append=T)
cat("Importing data : \n", file=log.file, append=T)
## ###############################################
##
##  import data 
##
## ###############################################
library(data.table)
library(dplyr)


MetaG_module.feature <- fread("database/KEGG_modules.tab", data.table = F, header = F, col.names = c("Module","Description","Features")) 

load("database/KO2METABO.lists.RData")

MetaB_module.feature <- fread("1_metaB.module_assign.txt", data.table = F,col.names =c("Feature","Module"))

MetaB.dat <- data.frame(fread("metabolome.txt", data.table = F), row.names = 1) %>% t() %>% as.data.frame()

metabo.KEGGmodule.match <- fread("database/Metabo.KEGGModule.match.txt",data.table = F)




## import metagenomic dataset

input.ds <- "metagenome.gct"
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
MetaG.dat <- m %>% t() %>% as.data.frame()


MetaG.B.dat <- merge(MetaB.dat, MetaG.dat, by=0)

## ###############################################
##
##  perform link analysis 
##
## ###############################################

cat("Performing link analysis : \n", file=log.file, append=T)

ACME.p.co = 0.25
MetaG.MetaB.modPairs <- (mediation.res %>% dplyr::filter(ACME.p <= ACME.p.co))$Treat_Mediator_Y


links <- NULL
for(mdp in MetaG.MetaB.modPairs){  # MetaG.MetaB.modPairs
  # mdp = MetaG.MetaB.modPairs[1]
  # stoped at mdp = MetaG.MetaB.modPairs[[352]]
  
  
  
  ConfirmedLink = FALSE
  
  # metag.md = sub("MetaG\\.(.*)$","\\1",mdp[1])
  # metab.md = sub("ME(.*)$","\\1",sub("MetaB\\.(.*)$","\\1",mdp[2]) )
  parts = strsplit(mdp,"_", fixed = T)[[1]]
  metag.md = sub("MetaG\\.(.*)$","\\1",parts[1])
  metab.md = sub("ME(.*)$","\\1",sub("MetaB\\.(.*)$","\\1",parts[2]) )
  
  
  #cat(paste("\n\nAnalyzing module pairs: ", parts[1], "_",parts[2],"\n", sep = ""),file=log.file, append=T)
  i_mdp <- which(MetaG.MetaB.modPairs == mdp)
  if(i_mdp %% 1000 == 0) cat(paste(" -------------------------- progress: ", i_mdp, " out of ", length(MetaG.MetaB.modPairs)," module pairs -------------------------- \n", sep = ""), file=log.file, append=T)
  
  metag.ftr = strsplit(MetaG_module.feature$Features[which(MetaG_module.feature$Module == metag.md)], ";", fixed = T)[[1]]
  metab.ftr = MetaB_module.feature$Feature[MetaB_module.feature$Module == metab.md]
  
  
  # scenario1: if any metaG feature (KO number) has a product/substrate compound which belong to the metaB module  ------------
  
  #cat("Performing scenario 1 analysis.\n", file=log.file, append=T)
  
  # writeLines("performing scenario 1 analysis.\n")
  for(gf in metag.ftr){
    #  gf = metag.ftr[1]
    if(ConfirmedLink) next
    
    # print(gf)
    substrates = KO.Substrates.metab_list[[gf]]
    products = KO.Products.metab_list[[gf]]
    subsprod =  KO.SubsProd.metab_list[[gf]]
    
    tmp <- unique(c(substrates,products, subsprod))
    
    if(is.null(tmp) ) next
    
    
    if(!any(tmp %in% metab.ftr)) next
    
    if(T){
      bf = tmp[1]
      link = c(paste(parts[1],parts[2],sep = "_"), gf,  bf )
      names(link) <- c("MetaG.MetaB_modulePair", "MetaG.feature","MetaB.feature")
      ConfirmedLink <- T
      links <- bind_rows(links,link)
      remove(link)
    } # if not consider spearman r 
    
    
    
    if(F){
      bf_check <- tmp[tmp %in% metab.ftr]
      for(bf in bf_check){
        if(ConfirmedLink) next # to save time, we don't need to perform multiple checks if a biological link has been confirmed 
        
        if(bf %in% subsprod) { 
          link = c(paste(parts[1],parts[2],sep = "_"), gf,  bf )
          names(link) <- c("MetaG.MetaB_modulePair", "MetaG.feature", "MetaB.feature")
          ConfirmedLink <- T
          
        }else{
          
          if(all(c(gf,bf) %in% colnames(MetaG.B.dat))) cor.dat <- MetaG.B.dat %>% dplyr::select(all_of(c(gf,bf))) else next
          test.r = cor(cor.dat[,1], cor.dat[,2], method = "spearman")
          
          if(bf %in% products & test.r > 0){
            link = c(paste(parts[1],parts[2],sep = "_"), gf, bf )
            names(link) <- c("MetaG.MetaB_modulePair", "MetaG.feature","MetaB.feature")
            ConfirmedLink <- T
          }else if(bf %in% substrates & test.r < 0){
            link = c(paste(parts[1],parts[2],sep = "_"), gf, bf )
            names(link) <- c("MetaG.MetaB_modulePair", "MetaG.feature", "MetaB.feature")
            ConfirmedLink <- T
          }else next # all other scenarios don't count as biological link
          
          
        }
        
        links <- bind_rows(links,link)
        remove(link)
      }# loop through metaB features
    }# if consider spearman r 
    
    
    
  }# scenario 1： loop through metaG features 
  
  if(ConfirmedLink) next
  
  # scenario2: if any metaB feature (C number) is directly connected to the metaG module through the metabo.KEGGmodule.match file ---------------
  
  #cat("Performing scenario 2 analysis.\n", file=log.file, append=T)
  
  # writeLines("performing scenario 2 analysis.\n")
  for(bf in metab.ftr){
    #bf = metab.ftr[1]
    if(ConfirmedLink) break  # to save time, we don't need to perform multiple checks if a biological link has been confirmed 
    gms <- metabo.KEGGmodule.match$KEGGmodule[which(metabo.KEGGmodule.match$metabo == bf)]
    
    if(any(gms == metag.md)){
      link = c(paste(parts[1],parts[2],sep = "_"),NA, bf )
      names(link) <- c("MetaG.MetaB_modulePair", "MetaG.feature", "MetaB.feature")
      ConfirmedLink <- T
      links <- bind_rows(links,link)
      remove(link)
      
    }
  } # loop through metab features
  
  if(ConfirmedLink) next
  
  
  # scenario 3: suppose this metaG module(Gm1) contains a metaG feature (gf1) which has a substrate/product (bf1),  --------------------------
  #             bf1 is not a member of this metaB moduel (Bm1), 
  #             however, bf1 together with another metab feature (bf2 which belong to Bm1) are connected by co-existing in another KEGG module (Gm2), 
  #             we consider Gm1 and Bm1 biologically linked ------------------------
  
  #cat("Performing scenario 3 analysis.\n", file=log.file, append=T)
  
  #writeLines("performing scenario 3 analysis.\n")
  
  for(gf1 in metag.ftr){
    # gf = metag.ftr[1]
    if(ConfirmedLink) break
    
    substrates = KO.Substrates.metab_list[[gf]]
    products = KO.Products.metab_list[[gf]]
    subsprod =  KO.SubsProd.metab_list[[gf]]
    
    tmp <- unique(c(substrates,products, subsprod))
    
    for(bf1 in tmp){
      # bf1 = tmp[2]
      if(ConfirmedLink) break
      bf1.keggM = metabo.KEGGmodule.match$KEGGmodule[which(metabo.KEGGmodule.match$metabo == bf1)]
      if(length(bf1.keggM) == 0) next
      
      for(bf2 in metab.ftr){
        # bf2 = metab.ftr[1]
        if(ConfirmedLink) break
        
        bf2.keggM = metabo.KEGGmodule.match$KEGGmodule[which(metabo.KEGGmodule.match$metabo == bf2)]
        if(length(intersect(bf1.keggM, bf2.keggM)) > 0) {
          link = c(paste(parts[1],parts[2],sep = "_"), gf1,bf2 )
          names(link) <- c("MetaG.MetaB_modulePair", "MetaG.feature", "MetaB.feature")
          ConfirmedLink <- T
          links <- bind_rows(links,link)
          remove(link)
          
        } # if link exist
        
      }# loop through bf2 
      
    }# loop through bf1
  }# loop through gf1
  if(ConfirmedLink) next
  
  # scenario 4: suppose this metaB module (Bm1) contains a metabo feature (bf1) which is a substrate/product of a metaG feature (gf1), ----------------
  #             gf1 is not a member of this metaG module (Gm1),
  #             however, if gf1, together with another metaG feature (gf2 which belongs to Gm1) co-exist in another metaG module (Gm2),
  #             we consider Bm1 and Gm1 biologically linked ----------------
  
  #cat("Performing scenario 4 analysis.\n", file=log.file, append=T)
  
  # writeLines("performing scenario 4 analysis.\n")
  
  for(bf1 in metab.ftr){
    # bf1 = metab.ftr[1]
    
    
    tmp <- unique(c(names(KO.Products.metab_list)[grep(bf1,  KO.Products.metab_list)],
                    names(KO.Substrates.metab_list)[grep(bf1, KO.Substrates.metab_list)],
                    names(KO.SubsProd.metab_list)[grep(bf1, KO.SubsProd.metab_list)]))
    
    for(gf1 in tmp){
      # gf1=tmp[1]
      if(ConfirmedLink) break
      
      gf1.Gms <- MetaG_module.feature$Module[which(MetaG_module.feature$Features == gf1)]
      if(length(gf1.Gms) == 0)  next
      
      for(gf2 in metag.ftr){
        # gf2 = metag.ftr[1]
        if(ConfirmedLink) break
        
        gf2.Gms <- MetaG_module.feature$Module[which(MetaG_module.feature$Features == gf2)]
        if(length(intersect(gf1.Gms, gf2.Gms)) > 0){
          link = c(paste(parts[1],parts[2],sep = "_"), gf2,  bf1 )
          names(link) <- c("MetaG.MetaB_modulePair", "MetaG.feature","MetaB.feature")
          ConfirmedLink <- T
          links <- bind_rows(links,link)
          remove(link)
        } # if link exists
        
      }# loop through gf2
    }# loop through gf1
  }# loop through bf1
  
}# loop through module pairs

if(!dir.exists(output.dir)) dir.create(output.dir)
write.table(links, file = paste(output.dir,"/4_MetaG.MetaB.modules.linked.txt",sep = ""), sep = "\t", quote = F, row.names = F)

