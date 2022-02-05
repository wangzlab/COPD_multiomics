# files required:
# 1) "3_HostT_affects_NEU_through_HostP.txt": resulting data frame generated from script "4.mediation_hostT.2.hostP.r"
# 2) "transcriptome.txt" transcriptomic abundance table on the feature level
# 3) "sputum_proteome.txt" proteomic abundance table on the feature level
# 4) "1_hostT.module_assign.txt" module assignment file of the transcriptomic featrues, generated from script "1.WGCNA_transcriptomics.r"
# 5) "protein_info.txt" protein names and corresponding gene names
# 6) "human_pathway.gmt"  Human pathway information file (.gmt file) 

# output:
# "4_HostT.module_HostP.protein.linked.txt" stored in the "Output" directory


output.dir = "Output"
log.file = "HostT.HostP.link.log"


## ######################################################################
##
## import data 
##
## #####################################################################
library(data.table)
library(dplyr)

cat(paste("\n\n",as.character(Sys.time()), '\n'),  file=log.file, append=T)
cat("Importing data : \n", file=log.file, append=T)

mediation.res <- fread("3_HosT_affects_NEU_through_HosP.txt")


# impport host t data
m1 <- data.frame(fread("transcriptome.txt"), row.names=1)
feature.abb_df1 <- cbind.data.frame(feature = rownames(m1),
                                    abb = paste("feature",seq(1,nrow(m1),1),sep = ""),
                                    stringsAsFactors = F)
rownames(m1) <- sapply(rownames(m1), function(x) feature.abb_df1$abb[which(feature.abb_df1$feature == x)])
m1 <- t(m1) %>% as.data.frame(stringsAsFactors=F)
colnames(m1) <- sapply(colnames(m1), function(x) feature.abb_df1$feature[which(feature.abb_df1$abb == x)])

HostT.dat <- m1



# import host p data
m2 =  data.frame(fread("sputum_proteome.txt"), row.names = 1 )
feature.abb_df2 <- cbind.data.frame(feature = rownames(m2),
                                    abb = paste("feature",seq(1,nrow(m2),1),sep = ""),
                                    stringsAsFactors = F)
rownames(m2) <- sapply(rownames(m2), function(x) feature.abb_df2$abb[which(feature.abb_df2$feature == x)])
m2 <- t(m2) %>% as.data.frame(stringsAsFactors=F)
colnames(m2) <- sapply(colnames(m2), function(x) feature.abb_df2$feature[which(feature.abb_df2$abb == x)])

HostP.dat <- m2

# hostT module - feature assignment

HostT_module.feature <- fread("1_hostT.module_assign.txt", data.table = F, col.names =c("Feature","Module"))

# protein - gene information

temp <- fread("database/protein_info.txt", data.table = T) %>% dplyr::select(Target, GeneName) %>% unique()
Protein.Gene_list <- vector("list", length = nrow(temp))
names(Protein.Gene_list) <- temp$Target
for(i in c(1:nrow(temp))){
  #i=1
  prot =  temp$Target[i]
  genenames = strsplit(temp$GeneName[i], split = '[\\s\\/]', perl = T )[[1]]
  Protein.Gene_list[[prot]] <- genenames
  
}

# pathway gene information 

PTWY_list <- Read.GeneSets.db2("database/human_pathway.gmt")
Pthwy.genes_list <- PTWY_list$gs  # the gene names
PTWY_list$gs.desc[1:3] # the pathway id


# merge data
HostT.HostP.dat <- merge(HostT.dat, HostP.dat, by=0)



## ###############################################
##
##  perform linke analysis 
##
## ###############################################
cat("Performing link analysis : \n", file=log.file, append=T)
# identify HostT-HostP module pairs 
# HostT.HostP.modPairs <- strsplit((mediation.res %>% dplyr::filter(ACME.p <= ACME.p.co))$Treat_Mediator_Y,"_",fixed = T)
ACME.p.co = 0.25
HostT.HostP.modPairs <-(mediation.res %>% dplyr::filter(ACME.p <= ACME.p.co))$Treat_Mediator_Y

# identify links
links <- NULL
for(mdp in HostT.HostP.modPairs){
  # mdp = HostT.HostP.modPairs[1]
  
  
  parts = strsplit(mdp,"_", fixed = T)[[1]]
  
  
  
  #cat(paste("\n\nAnalyzing module pairs: ", parts[1], "_",parts[2],"\n", sep = ""), file=log.file, append=T)
  i_mdp <- which(HostT.HostP.modPairs == mdp)
  if(i_mdp %% 100 == 0) cat(paste(" -------------------- progress: ", i_mdp, " out of ", length(HostT.HostP.modPairs)," module pairs -------------------------- \n", sep = ""), file=log.file, append=T)
  
  
  
  
  
  
  # identify hostp associated gene
  hostp.pro = sub("ME(.*)$","\\1",sub("HostP\\.(.*)$","\\1",parts[2]) )
  hostp.gene = Protein.Gene_list[[hostp.pro]]
  
  # identify hostt associated gene 
  # identify hostt associated gene by module assignment
  hostt.md = sub("ME(.*)$","\\1",sub("HostT\\.(.*)$","\\1",parts[1]) )
  hostt.ftr = HostT_module.feature$Feature[HostT_module.feature$Module == hostt.md]
  
  # identify hostt associated gene through pathway enrichment analysis
  pathwy.pvalue <- NULL
  for(pthy in names(Pthwy.genes_list)[!grepl("^GENEGO", names(Pthwy.genes_list))]){
    #pthy = names(Pthwy.genes_list)[!grepl("^GENEGO", names(Pthwy.genes_list))][1]
    
    pthy.genes = Pthwy.genes_list[[pthy]]
    overlaped = intersect(pthy.genes, hostt.ftr)
    totalgene.no = ncol(HostT.dat)
    
    pvalue <-dhyper(x=length(overlaped), m=length(pthy.genes), n= totalgene.no - length(pthy.genes), k = length(hostt.ftr))
    p.vec_temp <-  c(pthy, pvalue)
    names(p.vec_temp) <- c("Pathway","dhyper.pvalue")
    pathwy.pvalue <- bind_rows(pathwy.pvalue, p.vec_temp)
  }
  pathwy.pvalue <- pathwy.pvalue %>% 
    mutate(fdr = p.adjust(dhyper.pvalue, method = "fdr"))  %>% 
    dplyr::arrange(fdr) %>% filter(fdr <= 0.25)
  
  if(nrow(pathwy.pvalue) == 0) {
    hostt.ftr.integrated <- hostt.ftr
  } else {
    hostt.ftr.pthwy = Pthwy.genes_list[[pathwy.pvalue$Pathway[1]]]
    hostt.ftr.integrated <- unique(c(hostt.ftr,hostt.ftr.pthwy ))
  } 
  
  
  if(length( intersect(hostt.ftr.integrated, hostp.gene)) > 0) {
    link = c(paste(parts[1], parts[2], sep = "_" ), paste(intersect(hostt.ftr.integrated, hostp.gene), collapse = ";") )
    names(link) <- c("HostT.module_HostP.protein_pair", "Genes.association")
    links <- bind_rows(links,link)
    remove(link)
  }  
  
} # loop through module pairs  


if(!dir.exists(output.dir)) dir.create(output.dir)
write.table(links, file = paste(output.dir,"/HostT.module_HostP.protein.NEU.linked.txt",sep = ""), sep = "\t", quote = F, row.names = F)
