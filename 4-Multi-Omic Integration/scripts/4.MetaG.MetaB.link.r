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

metabo.KEGGmodule.match <- fread("database/Metabo.KEGGModule.match.txt",data.table = F)






## ###############################################
##
##  perform link analysis 
##
## ###############################################

cat("Performing link analysis : \n", file=log.file, append=T)


MetaG.MetaB.modPairs <- (mediation.res %>% dplyr::filter(ACME.p <= ACME.p.co))$Treat_Mediator_Y

links_df <- NULL # a data frame to store links
for(i_mdp in c(1:length(MetaG.MetaB.modPairs))){  # MetaG.MetaB.modPairs
  mdp = MetaG.MetaB.modPairs[i_mdp]
  
  ConfirmedLink = FALSE
  
  ftr_connections <- NULL  # define a vector for a new module pair
  
  parts = strsplit(mdp,"_", fixed = T)[[1]]
  metag.md = sub("MetaG\\.(.*)$","\\1",parts[1])
  metab.md = sub("ME(.*)$","\\1",sub("MetaB\\.(.*)$","\\1",parts[2]) )
  
  if(i_mdp %% 100 == 0) cat(paste(" -------------------------- progress: ", i_mdp, " out of ", length(MetaG.MetaB.modPairs)," module pairs -------------------------- \n", sep = ""), file=log.file, append=T)
  
  metag.ftr = strsplit(MetaG_module.feature$Features[which(MetaG_module.feature$Module == metag.md)], ";", fixed = T)[[1]]
  metab.ftr = MetaB_module.feature$Feature[MetaB_module.feature$Module == metab.md]
  
  # scenario1: if any metaG feature (KO number) has a product/substrate compound which belong to the metaB module  ------------
  
  #cat("Performing scenario 1 analysis.\n", file=log.file, append=T)
  
  # writeLines("performing scenario 1 analysis.\n")
  for(gf in metag.ftr){
    #  gf = metag.ftr[1]
    
    # print(gf)
    substrates = KO.Substrates.metab_list[[gf]]
    products = KO.Products.metab_list[[gf]]
    subsprod =  KO.SubsProd.metab_list[[gf]]
    
    tmp <- unique(c(substrates,products, subsprod))
    
    if(is.null(tmp) ) next
    
    
    if(!any(tmp %in% metab.ftr)) next
    
    bf = tmp[tmp %in% metab.ftr] 
    
    link = paste(gf,bf,sep = "_")
    ConfirmedLink <- T
    ftr_connections <- append(ftr_connections, link)
    remove(link)
    
  }# scenario 1： loop through metaG features 
  
  
  if(ConfirmedLink){
    # store as data frame 
    link.vec <- c(paste(parts[1],parts[2],sep = "_"), paste(ftr_connections,collapse = ";"))
    names(link.vec) <- c("module_pair", "feature_connections")
    links_df <- bind_rows(links_df, link.vec)
    next
  } 
  
  if(F){
    # scenario2: if any metaB feature (C number) is directly connected to the metaG module through the metabo.KEGGmodule.match file ---------------
    
    #cat("Performing scenario 2 analysis.\n", file=log.file, append=T)
    
    # writeLines("performing scenario 2 analysis.\n")
    for(bf in metab.ftr){
      #bf = metab.ftr[1]
      gms <- metabo.KEGGmodule.match$KEGGmodule[which(metabo.KEGGmodule.match$metabo == bf)]
      
      if(any(gms == metag.md)){
        link = paste(metag.ftr, bf, sep = "_") # contain all KO_C link
        ftr_connections <- append(ftr_connections, link)
        ConfirmedLink <- T
        remove(link)
        
      }
    } #scenario2： loop through metab features
    
  }
  
  
  # scenario2-extended: suppose the metaG module (Gm1) contains a feature (gf1), which does not linke to any features in the metaB module (Bm1) ------------------------
  #                     however, the gf1 could be assigned to multiple KEGG modules (gm.extend1, gm.extend2, ...),
  #                     if any feature (e.g. bf1) in Bm1 is found to be connected to any gm.extend in the metabo.KEGGmodule.match table,
  #                     we consider  consider Gm1 and Bm1 biologically linked  ---------------------------------------
  
  #cat("Performing scenario 2-extended analysis.\n", file=log.file, append=T)
  
  #writeLines("performing scenario 2-extended analysis.\n")
  
  for(gf in metag.ftr){
    #gf = metag.ftr[1]
    gm.extend <- MetaG_module.feature$Module[grep(gf, MetaG_module.feature$Features)]
    
    
    for(bf in metab.ftr){
      # bf = metab.ftr[1]
      gm.byBf <- metabo.KEGGmodule.match$KEGGmodule[which(metabo.KEGGmodule.match$metabo == bf)]
      
      if(any(gm.extend %in% gm.byBf)){
        link = paste(gf, bf, sep = "_") # contain all KO_C link
        ftr_connections <- append(ftr_connections, link)
        ConfirmedLink <- T
        remove(link)
      }
      
    } # loop through metab.ftr
    
  }  # scenario 2-extended: loop through metag.ftr
  
  if(ConfirmedLink){
    # store as data frame 
    link.vec <- c(paste(parts[1],parts[2],sep = "_"), paste(ftr_connections,collapse = ";"))
    names(link.vec) <- c("module_pair", "feature_connections")
    links_df <- bind_rows(links_df, link.vec)
    next
  }  
  
  # scenario 3: suppose this metaG module(Gm1) contains a metaG feature (gf1) which has a substrate/product (bf1),  --------------------------
  #             bf1 is not a member of this metaB moduel (Bm1), 
  #             however, bf1 together with another metab feature (bf2 which belong to Bm1) are connected by co-existing in another KEGG module (Gm2), 
  #             we consider Gm1 and Bm1 biologically linked ------------------------
  
  #cat("Performing scenario 3 analysis.\n", file=log.file, append=T)
  
  #writeLines("performing scenario 3 analysis.\n")
  
  for(gf1 in metag.ftr){
    # gf = metag.ftr[1]
    
    substrates = KO.Substrates.metab_list[[gf]]
    products = KO.Products.metab_list[[gf]]
    subsprod =  KO.SubsProd.metab_list[[gf]]
    
    tmp <- unique(c(substrates,products, subsprod))
    
    for(bf1 in tmp){
      # bf1 = tmp[2]
      
      bf1.keggM = metabo.KEGGmodule.match$KEGGmodule[which(metabo.KEGGmodule.match$metabo == bf1)]
      if(length(bf1.keggM) == 0) next
      
      for(bf2 in metab.ftr){
        # bf2 = metab.ftr[1]
        
        bf2.keggM = metabo.KEGGmodule.match$KEGGmodule[which(metabo.KEGGmodule.match$metabo == bf2)]
        if(length(intersect(bf1.keggM, bf2.keggM)) > 0) {
          
          link = paste(gf1, bf2, sep = "_") # contain all KO_C link
          ftr_connections <- append(ftr_connections, link)
          ConfirmedLink <- T
          remove(link)
          
        } # if link exist
        
      }# loop through bf2： metab.ftr 
      
    }# loop through bf1：substrates and products by gf1
  }# scenario3: loop through gf1: metag.ftr
  
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
      
      gf1.Gms <-  MetaG_module.feature$Module[grep(gf1, MetaG_module.feature$Features)] 
      if(length(gf1.Gms) == 0)  next
      
      for(gf2 in metag.ftr){
        # gf2 = metag.ftr[1]
        
        gf2.Gms <- MetaG_module.feature$Module[grep(gf2, MetaG_module.feature$Features)]
        if(length(intersect(gf1.Gms, gf2.Gms)) > 0){
          
          link = paste(gf2, bf1, sep = "_") # contain all KO_C link
          ftr_connections <- append(ftr_connections, link)
          ConfirmedLink <- T
          remove(link)
          
        } # if link exists
        
      }# loop through gf2
    }# loop through gf1 (gfs connected to bf1 by products/substrate relationship)
  }# scenario 4: loop through bf1 (metab.ftr)
  
  
  if(ConfirmedLink){
    # store as data frame 
    link.vec <- c(paste(parts[1],parts[2],sep = "_"), paste(ftr_connections,collapse = ";"))
    names(link.vec) <- c("module_pair", "feature_connections")
    links_df <- bind_rows(links_df, link.vec)
  }  
  
}  # loop through MetaG.MetaB.modPairs

#save results
if(!dir.exists(output.dir)) dir.create(output.dir)
# save as table
write.table(links_df, file = paste(output.dir,"/",output.prefix,"MetaG.MetaB.modules.linked.txt",sep = ""), sep = "\t", quote = F, row.names = F)
# save as RData
# save(links_df, file = paste(output.dir,"/",output.prefix,"MetaG.MetaB.modules.linked.RData",sep = ""))