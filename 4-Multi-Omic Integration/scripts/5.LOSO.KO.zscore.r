# files required:
# 1) "5_LOSO.ko.abund.dir":  the directory containing all ko abundance gct files from the LOSO.ko function. One gct file for one species removed.
# 2) "complete.ko.abund": original metagenomic feature abundane, without any species removed. gct file or txt file (rows are different features and columns are samples ).
# 3) "meta.disease" file or a data frame containig one column names "SampleID" and one column indicating Disease state, based on which the fold change of KO abundance before and after disease shall be calculated 

# output:
# "LOSO_ko.fc/LOSO.ko.zscore.txt" file containg zscore for the ko for each species removed. 




LOSO.ko.abund.dir = "./5_LOSO_ko.abund"
complete.ko.abund = "metagenome.gct"
meta.disease = fread("metadata.txt",data.table = F)
Disease.column = "Disease"
Healthy.state = "0"
Disease.state = "1"
output.dir = "."

library(data.table)

# read meta data -----------------------

if(class(meta.disease) == "character"){
  meta <- fread(meta.disease, data.table = F)
}else{
  meta <- meta.disease
}

colnames(meta)[colnames(meta) == Disease.column] <- "Y"
meta$Y <- as.character(meta$Y)

# read ko abundance in loop and calculate fold change ---------------------------
combined_FC.res <- NULL

ko.files <- list.files(LOSO.ko.abund.dir,full.names = T)
for(f in ko.files){
  # f = ko.files[1]
  
  # read ko abundance file 
  input.ds <- f
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
  ko.dat <- m %>% t() %>% as.data.frame()
  tmp <- ko.dat %>% tibble::rownames_to_column("SampleID") %>% reshape2::melt(variable.name="ko")
  ko.dat.rel <- 
    merge(tmp ,
          tmp %>% group_by(SampleID) %>% summarise(sum.sp = sum(value)),
          by = "SampleID") %>% 
    mutate(freq = value/sum.sp) 
  
  # ko.dat.rel %>% group_by(SampleID) %>% summarise(test=sum(freq))
  
  ko.dat.rel <- ko.dat.rel %>%
    reshape2::dcast(SampleID ~ ko, value.var = "freq") %>% 
    tibble::column_to_rownames("SampleID")
  
  #minRelVal = min(ko.dat.rel[ko.dat.rel != 0])
  ko.dat.rel <- ko.dat.rel + 0.00001 # to avoid dividing by 0
  dat <- merge(ko.dat.rel, meta %>% select(SampleID, Y), by.x=0, by.y="SampleID") %>%
    tibble::column_to_rownames("Row.names")
  
  
  # calculate fold change for each ko ----------------
  
  FC.dat <-dat %>% tibble::rownames_to_column("SampleID") %>%
    reshape2::melt(id.vars=c("SampleID","Y"), variable.name = "KO") %>%
    group_by(KO) %>%
    summarise(log2FC = log2(mean(value[Y==Disease.state], na.rm=T)/mean(value[Y==Healthy.state], na.rm=T) ) )
  
  #FC.dat$FC[is.na(FC.dat$FC)] <- 0 
  colnames(FC.dat)[colnames(FC.dat) == "log2FC"] <- paste( "log2FC.", sub("\\.gct", "", basename(f)),sep = "")
  
  combined_FC.res <- merge(combined_FC.res, FC.dat, by.x=0, by.y="KO", all=T) %>% tibble::column_to_rownames("Row.names")
  
  
}# loop through ko files for each species 
combined_FC.res[is.na(combined_FC.res)] <- 0

# read original metagenomic feature abundance ---------------------------
if(grepl("gct$",complete.ko.abund, perl = T)){
  input.ds <- complete.ko.abund
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
  ko.dat.orig <- m %>% t() %>% as.data.frame()
  
}else if(grepl("txt$", complete.ko.abund, perl = T)){
  
  ko.dat.orig <- fread(complete.ko.abund, data.table = F) 
  rownames(ko.dat.orig) <- ko.dat.orig[,1]
  ko.dat.orig <- ko.dat.orig[,-1]
  ko.dat.orig <- ko.dat.orig %>% t() %>% data.frame()
  
}

tmp <- ko.dat.orig %>% 
  tibble::rownames_to_column("SampleID") %>%
  reshape2::melt(variable.name="ko")

ko.dat.orig.rel <- merge(tmp,
                         tmp %>% group_by(SampleID) %>% summarise(sum.sp = sum(value)), by="SampleID" ) %>%
  mutate(freq=value/sum.sp)

ko.dat.orig.rel %>% group_by(SampleID) %>% summarise(test=sum(freq))


ko.dat.orig.rel <- ko.dat.orig.rel %>%
  reshape2::dcast(SampleID ~ ko, value.var = "freq") %>% 
  tibble::column_to_rownames("SampleID")

#minRelVal = min(ko.dat.orig.rel[ko.dat.orig.rel != 0])

ko.dat.orig.rel <- ko.dat.orig.rel + 0.00001 # to avoid dividing by 0
dat <- merge(ko.dat.orig.rel, meta %>% select(SampleID, Y), by.x=0, by.y="SampleID") %>%
  tibble::column_to_rownames("Row.names")

# calculate null fold change  for each ko ----------------

FCnull.dat <- dat %>% tibble::rownames_to_column("SampleID") %>%
  reshape2::melt(id.vars=c("SampleID","Y"), variable.name = "KO") %>%
  group_by(KO) %>%
  summarise(log2FC.null = log2(mean(value[Y==Disease.state], na.rm=T)/mean(value[Y==Healthy.state], na.rm=T) ) )

# calculate zscore for each ko --------------------
combined_FC.res.t <- t(combined_FC.res) %>% as.data.frame()

all.kos <- unique(c(as.character(FCnull.dat$KO), rownames(combined_FC.res)))

Zscore_df <- NULL
for(ko in all.kos){  # all.kos
  
  FCnull = if(ko %in% FCnull.dat$KO) FCnull.dat$log2FC.null[FCnull.dat$KO == ko] else next
  FCs_vec = combined_FC.res.t[,ko]
  FCmean = mean(FCs_vec)
  delta = FCs_vec - FCmean
  #denominator <- sqrt(sum(abs(FCs_vec - FCmean))/length(FCs_vec))
  denominator = sd(FCs_vec)
  
  if(denominator == 0) next
  
  Zscore_ko <- vector("numeric", length = ncol(combined_FC.res))
  names(Zscore_ko) <- colnames(combined_FC.res)
  for(si in colnames(combined_FC.res)){
    # si = colnames(combined_FC.res)[1]
    FCsi = combined_FC.res[ko, si]
    
    Zsi = abs(FCnull - FCsi)/denominator
    
    Zscore_ko[si] <- Zsi
  }
  
  Zscore_df <- bind_rows(Zscore_df, Zscore_ko)
  rownames(Zscore_df)[nrow(Zscore_df)] <- ko
  
} # loop through all.kos
colnames(Zscore_df) <- sub("log2FC\\.","",colnames(Zscore_df))
if(!dir.exists(output.dir)) dir.create(output.dir)
write.table(Zscore_df, file = paste(output.dir,"/5_LOSO.KO.zscore.txt", sep = ""), sep = '\t', quote = F, row.names = T)
