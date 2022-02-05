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
# "4_MetaG.MetaB.modules.linked.txt"



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


MetaG_module.feature <- fread("KEGG_modules.tab", data.table = F, header = F, col.names = c("Module","Description","Features")) 

load("KO2METABO.lists.RData")

MetaB_module.feature <- fread("1_metaB.module_assign.txt", data.table = F,col.names =c("Feature","Module"))

metabo.KEGGmodule.match <- fread("Metabo.KEGGModule.match.txt",data.table = F)






## ###############################################
##
##  perform link analysis 
##
## ###############################################

cat("Performing link analysis : \n", file=log.file, append=T)

ACME.p.co = 0.10
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
  
  # scenario 1: if any metaG feature has a product/substrate that belongs to the metaB module 
 
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
    # scenario 2: if any features in metaB are present in the metaG module (through the metabo.KEGGmodule.match file)
   
    for(bf in metab.ftr){
      #bf = metab.ftr[1]
      gms <- metabo.KEGGmodule.match$KEGGmodule[which(metabo.KEGGmodule.match$metabo == bf)]
      
      if(any(gms == metag.md)){
        link = paste(metag.ftr, bf, sep = "_") # contain all KO_C link
        ftr_connections <- append(ftr_connections, link)
        ConfirmedLink <- T
        remove(link)
        
      }
    } #scenario 2： loop through metaB features
    
  }
  
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
    
  }
  
  if(ConfirmedLink){
    # store as data frame 
    link.vec <- c(paste(parts[1],parts[2],sep = "_"), paste(ftr_connections,collapse = ";"))
    names(link.vec) <- c("module_pair", "feature_connections")
    links_df <- bind_rows(links_df, link.vec)
    next
  }  
  
  # scenario 3: if features in metaG have substrates/products that are linked to features in metaB by presenting in the same KEGG module
  
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
      }# loop through bf2    
    }# loop through bf1
  }# scenario3: loop through gf1: metag.ftr
  
  # scenario 4: if features in metaB are substrates/products for genes that are linked to features in metaG by presenting in the same KEGG module
  
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
    }# loop through gf1
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
