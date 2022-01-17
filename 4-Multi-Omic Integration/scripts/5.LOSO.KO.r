# files required:
# 1) "geneDepth.txt":  gene abundance table, with the first column being gene ids which are formatted as "scaffold_number"
# 2) "ko.txt" file giving informaiton on gene ids and KO number. 
# 3) "bin_membership.txt" file giving information on bins and scaffold id
# 4) "bin_species.txt" file giving information on bins and Species annotaion

# output:
# A directory "5_LOSO_ko.abund" containing many "ko.abund_rm.species.gct" files, one for each species removed.
# Each KO abundance files will go through ssGSEA to generate abundance table for KEGG modules. 
# For information of ssGSEA, see https://github.com/broadinstitute/ssGSEA2.0 



## #####################################################
##
## load data and delete unused data to save memory
##
## #####################################################
log.file='LOSO.ko.log'

cat(paste(as.character(Sys.time()), '\n'),  file=log.file, append=T)
cat('Performing LOSO.ko analysis: \n',  file=log.file, append=T)
cat('Importing data: \n',  file=log.file, append=T)

library(data.table)
library(tibble)
library(reshape2)
library(dplyr)

gene.KO_df <- fread("ko.txt",  header = T) %>% dplyr::select(Gene, KOnumber)

geneDepth_df <- fread("F:/temp.data/geneDepth.txt", data.table = F) 
colnames(geneDepth_df)[1] <- "GeneID"

geneDepth_df <- geneDepth_df %>%   # ,nrows = 300 to test
  dplyr::filter(GeneID %in%  gene.KO_df$Gene) %>% 
  tibble::column_to_rownames("GeneID") 

gc()

geneDepth_df <- merge(geneDepth_df,
                      gene.KO_df, by.x = 0, by.y="Gene") 
colnames(geneDepth_df)[which(colnames(geneDepth_df) == "Row.names")] <- "Gene"
remove(gene.KO_df)
gc() 



## #####################################################
##
## create species - scaffold file 
##
## #####################################################
bin.scaf_df <- fread("bin_membership.txt") %>% dplyr::select(Scaffold, Bin)
bin.species_df <- fread("bin_species.txt") %>% dplyr::select(Bin, Species)

species.scaffold_df <- merge(bin.species_df, bin.scaf_df, by="Bin")
remove(bin.scaf_df, bin.species_df)
gc() 

if(!dir.exists("LOSO_ko.abund") ) dir.create("LOSO_ko.abund")

uniqSpecies = unique(species.scaffold_df$Species)
rows.scaf = sub("(_\\d+$)","", geneDepth_df$Gene)

cat(paste('Excluding', length(uniqSpecies),' species one by one: \n'),  file=log.file, append=T)
for(species in uniqSpecies){
  i=which(uniqSpecies == species)
  cat(paste(i,".  ", species, ' \n', sep = ""),  file=log.file, append=T)
  
  #  species = uniqSpecies[3]
  
  scafs = species.scaffold_df$Scaffold[species.scaffold_df$Species == species]
  
  
  geneDepth_df_sub <- geneDepth_df[!(rows.scaf %in% scafs),]
  
  res <- geneDepth_df_sub %>%
    dplyr::select(-Gene) %>%
    # reshape2::melt(id.var = "KO", )
    dplyr::group_by(KOnumber) %>%
    dplyr::summarise_if(is.numeric, sum, na.rm = TRUE) %>%
    tibble::column_to_rownames("KOnumber")
  
  
  ## export gct
  gct.tmp <- new('GCT')
  gct.tmp@mat <- data.matrix( res )
  gct.tmp@rid <- rownames(res)
  gct.tmp@cid <- colnames(res)
  gct.tmp@cdesc <- data.frame(id=colnames(res))
  gct.tmp@rdesc <- data.frame(row.names = rownames(res), Description = rownames(res))
  fnn = paste("5_LOSO_ko.abund/ko.abund_rm.", species, ".gct", sep = "")
  gct.tmp@src <- fnn
  
  write.gct(gct.tmp, ofile=fnn, appenddim = F)  
}


