# files required:
# 1) "ko01000.keg" downloaded from KEGG
# 2) "metacyc_reactions.txt" storing information on susbstrates, products of EC numbers

# output:
# "KO2EC.list.RData", "EC2CMPD.lists.RData", "KO2CMPD.lists.RData" stored in the "database" directory

library(data.table)
library(dplyr)

## generate KO2EC list #############
dbDir = "database"
fileName <- paste(dbDir,"/ko01000.keg",sep = "")
conn <- file(fileName,open="r")
linn <-readLines(conn)
#for (i in 1:length(linn)){
#  print(linn[i])
#}
close(conn)



KO2EC_list <- vector("list",length=length(linn))
for(i in c(1:length(linn))){
  l = linn[i]
  
  if(!grepl("K\\d{5}",l,perl = T)) next
  ko <- sub(".*(K\\d{5}).*","\\1",l)
  
  if(grepl("\\[EC\\:(.*)\\]",l, perl = T)) ecs <- strsplit( sub(".*\\[EC\\:(.*)\\]$", "\\1", l, perl = T), " ", fixed = T)[[1]] else ecs <- ""
  
  if(ko %in% names(KO2EC_list)) {
    i_exist = which(names(KO2EC_list) == ko) 
    updated_ecs <-  unique(c( KO2EC_list[[i_exist]], ecs))
    updated_ecs <- updated_ecs[updated_ecs != ""]
    KO2EC_list[[i_exist]] <- updated_ecs
    next
  }
  
  KO2EC_list[[i]] <- ecs
  names(KO2EC_list)[i] <- ko
  
}


KO2EC_list <- KO2EC_list[!sapply(KO2EC_list, is.null) ] 
save(KO2EC_list, file = paste(dbDir,"/KO2EC.list.RData",sep = "")) 


## generate EC2CMPD lists ##################
metacyc_rxns <- fread(paste(dbDir,"/metacyc_reactions.txt",sep = "")) %>% filter(`EC-NUMBER` != "")
EC.Substrates_list <- vector("list", length = nrow(metacyc_rxns))
EC.Products_list <- vector("list", length = nrow(metacyc_rxns))
EC.SubsProd_list <- vector("list", length = nrow(metacyc_rxns))

for(i in c(1:nrow(metacyc_rxns))){
  
  ecs <- metacyc_rxns$`EC-NUMBER`[[i]]
  ecs <- sub("\\|?EC\\-(.*)\\|?", "\\1", ecs)
  rxn_drct = metacyc_rxns$`REACTION-DIRECTION`[[i]]
  
  
  if(grepl('LEFT-TO-RIGHT', rxn_drct) | rxn_drct == "") {
    substrates <- strsplit(metacyc_rxns$LEFT[i],";", fixed = T)[[1]]
    products <- strsplit(metacyc_rxns$RIGHT[i],";",fixed = T)[[1]]
    subs.prod <- NA
  }else if(grepl("RIGHT-TO-LEFT", rxn_drct)){
    substrates <- strsplit(metacyc_rxns$RIGHT[i], ";", fixed = T)[[1]]
    products <- strsplit(metacyc_rxns$LEFT[i],";",fixed = T)[[1]]
    subs.prod <- NA
  }else if(rxn_drct == "REVERSIBLE"){
    substrates <- NA
    products <- NA
    subs.prod <- unique(c(strsplit(metacyc_rxns$RIGHT[i], ";", fixed = T)[[1]], strsplit(metacyc_rxns$LEFT[i],";",fixed = T)[[1]] ))
  }
  
  EC.Substrates_list[[i]] <- substrates; names(EC.Substrates_list)[i] <- ecs
  EC.Products_list[[i]] <- products; names(EC.Products_list)[i] <- ecs
  EC.SubsProd_list[[i]] <- subs.prod; names(EC.SubsProd_list)[i] <- ecs
  
}

save(EC.Substrates_list, EC.Products_list, EC.SubsProd_list, 
     file = paste(dbDir,"/EC2CMPD.lists.RData",sep = ""))


## generate KO.Substrates, KO.Products and KO.SubsProd lists #########
KO.Substrates_list <- vector("list", length = length(KO2EC_list))
KO.Products_list <- vector("list", length = length(KO2EC_list))
KO.SubsProd_list <- vector("list", length = length(KO2EC_list))

for(i in c(1:length(KO2EC_list))){
  
  ko <- names(KO2EC_list)[i]
  ecs <- KO2EC_list[[ko]]
  
  products <- vector("character")
  substrates <- vector("character")
  subs.prod <- vector("character")
  
  for(ec in ecs){
    list.indexes <- which(names(EC.Products_list) == ec)
    
    for(l_ind in list.indexes){
      products <- unique(append(products, EC.Products_list[[l_ind]]) )
      substrates <- unique(append(substrates, EC.Substrates_list[[l_ind]]) )
      subs.prod <- unique(append(subs.prod, EC.SubsProd_list[[l_ind]]) )
    }
  }
  
  s.p <- intersect(products, substrates) 
  subs.prod <- unique(c(subs.prod, s.p)[!is.na(c(subs.prod, s.p))] )  # remove na and then unique
  
  substrates <- substrates[!is.na(substrates) & !(substrates %in% subs.prod)] # remove na and remove subs.prod
  products <- products[!is.na(products) & !(products %in% subs.prod)]
  
  KO.Substrates_list[[i]] <- substrates; names(KO.Substrates_list)[i] <- ko
  KO.Products_list[[i]] <- products; names(KO.Products_list)[i] <- ko
  KO.SubsProd_list[[i]] <- subs.prod; names(KO.SubsProd_list)[i] <- ko
  
}

save(KO.SubsProd_list, KO.Products_list, KO.Substrates_list,
     file = paste(dbDir,"/KO2CMPD.lists.RData",sep = ""))
