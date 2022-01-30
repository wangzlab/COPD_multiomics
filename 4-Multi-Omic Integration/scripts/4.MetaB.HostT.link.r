# files required:
# 1) "3_MetaB_affects_NEU_through_HostT.txt": resulting data frame generated from script "4.mediation_metaB.2.hostT.r"
# 2) "metabolome.txt" metabolomic abundance table on the feature level
# 3) "1_metaB.module_assign.txt" module assignment file of the metabolomic featrues, generated from script "1.WGCNA_metabolomics.r"
# 4) "transcriptome.txt" transcriptomic abundance table on the feature level
# 5) "1_hostT.module_assign.txt" module assignment file of the transcriptomic featrues, generated from script "1.WGCNA_transcriptomics.r"
# 6) "metabo2CIDm.txt" storing information on compound IDs and corresponding CIDm IDs
# 7) "all_CIDm_targets.txt" link file between CIDm ID and host targets

# output:
# "MetaB.HostT.modules.linked.txt" stored in the "Output" directory



output.dir = "Output"
log.file = "MetaB.HostT.link.log"


## ###############################################
##
##  import data 
##
## ###############################################
library(data.table)
library(dplyr)

cat(paste("\n\n",as.character(Sys.time()), '\n'),  file=log.file, append=T)
cat("Importing data : \n", file=log.file, append=T)

mediation.res = fread("3_MetaB_affects_NEU_through_HostT.txt")
MetaB.dat <- data.frame(fread("metabolome.txt"), row.names = 1) %>% t() %>% as.data.frame()
MetaB_module.feature <- fread("1_metaB.module_assign.txt", data.table = F, col.names =c("Feature","Module"))
HostT.dat <- data.frame(fread("transcriptome.txt"), row.names = 1) %>% t() %>% as.data.frame()
HostT_module.feature <- fread("1_hostT.module_assign.txt", data.table = F, col.names =c("Feature","Module"))
metabo2CIDm <- fread("database/metabo2CIDm.txt", data.table = F, header = F, col.names = c("Metabo","CIDm"))

# load and organize CIDm.receptor data 
conn <- file("database/all_cidm_receptor.txt", open="r")
linn <-readLines(conn)
close(conn)


numElements <- sapply(linn, function(x) length(strsplit(x, "\t", fixed = T)[[1]]))
table(numElements) # 4 levels of lengths

tmp <- linn[numElements == 4]
subdf1 <- unname(sapply(tmp,function(x) strsplit(x,"\t")[[1]] )) %>% t() %>% data.frame(stringsAsFactors = F)
colnames(subdf1) <- c("CIDm","ENSP","geneName","linkType")


tmp <- linn[numElements == 5]
subdf2 <- unname(sapply(tmp,function(x) strsplit(x,"\t")[[1]] )) %>% t() %>% data.frame(stringsAsFactors = F)
colnames(subdf2) <- c("CIDm","ENSP","geneName","description", "linkType")

tmp <- linn[numElements == 9]
subdf3 <- unname(sapply(tmp,function(x) strsplit(x,"\t")[[1]] )) %>% t() %>% data.frame(stringsAsFactors = F)
colnames(subdf3) <- c("CIDm","ENSP","geneName","linkType","X1","X2","X3","X4", "score") 

tmp <- linn[numElements == 10]
subdf4 <- unname(sapply(tmp,function(x) strsplit(x,"\t")[[1]] )) %>% t() %>% data.frame(stringsAsFactors = F)
colnames(subdf4) <- c("CIDm","ENSP","geneName","description", "linkType","X1","X2","X3","X4","score") 

CIDm.receptor.score.co = 800
CIDm.receptor <- bind_rows(subdf1, subdf2, subdf3, subdf4) %>% filter(score >= CIDm.receptor.score.co)


# merge data
MetaB.HostT.dat <- merge(MetaB.dat, HostT.dat, by=0)


## ###############################################
##
##  perform linke analysis 
##
## ###############################################
cat("Performing link analysis : \n", file=log.file, append=T)
# identify MetaB-HostT module pairs 
#MetaB.HostT.modPairs <- strsplit((mediation.res %>% dplyr::filter(ACME.p <= ACME.p.co))$Treat_Mediator_Y,"_",fixed = T)
ACME.p.co = 0.25
MetaB.HostT.modPairs <- (mediation.res %>% dplyr::filter(ACME.p <= ACME.p.co))$Treat_Mediator_Y


# identify links
links <- NULL
for(mdp in MetaB.HostT.modPairs){
  # mdp = MetaB.HostT.modPairs[[1]]
  
  ConfirmedLink = FALSE
  #metab.md = sub("ME(.*)$","\\1",sub("MetaB\\.(.*)$","\\1",mdp[1]) )
  #hostt.md = sub("ME(.*)$","\\1",sub("HostT\\.(.*)$","\\1",mdp[2]) )
  
  parts = strsplit(mdp,"_", fixed = T)[[1]]
  metab.md = sub("ME(.*)$","\\1",sub("MetaB\\.(.*)$","\\1",parts[1]) )
  hostt.md = sub("ME(.*)$","\\1",sub("HostT\\.(.*)$","\\1",parts[2]) )
  
  #cat(paste("\n\nAnalyzing module pairs: ", parts[1], "_",parts[2],"\n", sep = ""), file=log.file, append=T)
  i_mdp <- which(MetaB.HostT.modPairs == mdp)
  if(i_mdp %% 1000 == 0) cat(paste(" -------------------------- progress: ", i_mdp, " out of ", length(MetaB.HostT.modPairs)," module pairs -------------------------- \n", sep = ""), file=log.file, append=T)
  
  
  
  
  metab.ftr = MetaB_module.feature$Feature[MetaB_module.feature$Module == metab.md]
  hostt.ftr = HostT_module.feature$Feature[HostT_module.feature$Module == hostt.md]
  
  
  for(bf in metab.ftr){
    # bf=metab.ftr[2]
    # writeLines(paste("the ", which(metab.ftr == bf), " bf", sep = "") )
    
    if(ConfirmedLink ) break
    cidm = metabo2CIDm$CIDm[which(metabo2CIDm$Metabo == bf)]
    
    if(length(cidm) == 0) {
      writeLines(paste("MetaB feature ",bf, " (belonged to the ", parts[1]," module) doesn't have matching CIDm", sep=""))
      next
    }
    
    if(length(cidm) > 1) cidm <- as.character((data.frame(table(cidm)) %>% dplyr::arrange(desc(Freq)))$cidm[1])
    
    receptors_df = CIDm.receptor %>% filter(CIDm == cidm)
    receptors <- receptors_df$geneName
    
    receptors <- receptors[receptors != ""]
    
    if( length(receptors) == 0 ) next 
    if(!any(receptors %in% hostt.ftr)) next
    
    # verify the links
    recptr_check <- receptors[receptors %in% hostt.ftr]
    receptors_df <- receptors_df %>% filter(geneName %in%  recptr_check)
    
    for(i in c(1:nrow(receptors_df))){
      
      rcpt <- receptors_df$geneName[i]
      lkType <- receptors_df$linkType[i]
      
      if(lkType %in% c("catalysis","reaction", "expression")) next
      
      link = c(paste(parts[1], parts[2], sep = "_" ), bf, rcpt, lkType )
      names(link) <- c("MetaB.HostT_modulePair", "MetaB.feature","HostT.feature","Link.by")
      
      links <- bind_rows(links,link)
      remove(link)
      ConfirmedLink = T
      
    } # loop through all the receptors to be checked 
      

    
  }# loop through metab.features
  
  
  
} # loop through module pairs

if(!dir.exists(output.dir)) dir.create(output.dir)
write.table(links, file = paste(output.dir,"/4_MetaB.HostT.modules.NEU.linked.txt",sep = ""), sep = "\t", quote = F, row.names = F)

