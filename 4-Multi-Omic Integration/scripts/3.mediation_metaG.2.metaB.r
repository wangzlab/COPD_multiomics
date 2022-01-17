# files required:
# 1) MetaG module abundance table ("1_metaG-combined.gct") generated by ssGSEA2
# 2) MetaB module abundance table ("1_metaB.module_eigengene.txt") generated by WGCNA
# 3) metadata ("meta.mediation.NUE.txt") containing clinical variable NEU
# 4) significant MetaG modules ("2_significant_metaG_modules.RData") generated from script '2.significant.metaG.modules.r'
# 5) significant MetaB modules ("2_significant_MetaB_modules.RData") generated from script '3.significant.metaB.modules.r'

# output:
# File "MetaG_affects_NEU_through_MetaB.txt" recording ACME.p,	ADE.p, and	prop.mediated of each mediation route (MetaG module -- MetaB module -- clinical variable NEU)

Treat.omic = "MetaG"
Mediator.omic = "MetaB"
Y = "NEU"

load("2_significant_metaG_modules.RData")
Treat.omic.sigModules <- MetaG.sigMods
load("2_significant_MetaB_modules.RData")
Mediator.omic.sigModules <- MetaB.sigMods
outputDir = "Output"
log.file = "mediation.MetaG2MetaB.log"

library(mediation)
m1 <- parse.gctx("1_metaG-combined.gct")@mat %>% t() %>% data.frame() 

MetaB.Mod.dat <- fread("1_metaB.module_eigengene.txt",data.table = F) %>%
  dplyr::filter(!grepl("#",`#NAME`,fixed=T)) 
MetaB.Mod.dat[-1] <- sapply(MetaB.Mod.dat[-1], as.numeric)
MetaB.Mod.dat <- MetaB.Mod.dat %>% tibble::column_to_rownames("#NAME")
m2 = MetaB.Mod.dat
feature.abb_df2 <- cbind.data.frame(feature = rownames(m2),
                                    abb = paste("feature",seq(1,nrow(m2),1),sep = ""),
                                    stringsAsFactors = F)
rownames(m2) <- sapply(rownames(m2), function(x) feature.abb_df2$abb[which(feature.abb_df2$feature == x)])
m2 <- t(m2) %>% as.data.frame(stringsAsFactors=F)
colnames(m2) <- sapply(colnames(m2), function(x) feature.abb_df2$feature[which(feature.abb_df2$abb == x)])



meta_df <- fread("meta.mediation.NUE.txt", data.table = F)

if(!dir.exists(outputDir)) dir.create(outputDir)

cat(paste(as.character(Sys.time()), '\n'),  file=log.file, append=T)
cat('Performing mediation analysis: \n',  file=log.file, append=T)

Mediation.results <- NULL
for(md1 in Treat.omic.sigModules){
  
  for(md2 in Mediator.omic.sigModules){
    
    id = paste(Treat.omic, ".", md1, "_", Mediator.omic, ".", md2, "_", Y,sep = "" )
    cat(paste(id, "\n", sep = ""),  file=log.file, append=T)
    
    dat <- base::merge(base::merge(meta_df, m1 %>% dplyr::select(all_of(md1)), by.x = "SampleID", by.y = 0),
                       m2 %>% dplyr::select(all_of(md2)), by.x = "SampleID", by.y=0)
    
    colnames(dat)[(ncol(dat)-1):ncol(dat)] <- c("Treat", "Mediator")
    
    # write.table(dat, file = paste(outputDir, "/",Treat.omic, "_", md1, "_", Mediator.omic, "_", md2,".txt",sep = ""),
    #             quote = F, row.names = F, sep = "\t")
    
    colnames(dat)[which(colnames(dat) == Y)] <- "Y"
    
    # remove na
    
    dat <- dat %>% filter(!is.na(Y)) %>% filter(!is.na(Treat)) %>% filter(!is.na(Mediator))
    
    #model.m=lm(Mediator ~ Treat+Age+Gender+CurrentSmoking+ICS,dat)
    model.m = lm(as.formula( paste("Mediator ~ Treat + ", 
                                   paste(colnames(dat)[!(colnames(dat) %in% c("SampleID","Y","Mediator","Treat"))], collapse = " + "),
                                   sep = "") ), data = dat)
    
    #model.y=lm(Y~Treat+Mediator+Age+Gender+CurrentSmoking+ICS,dat)
    model.y = lm(as.formula( paste("Y ~ Treat + Mediator + ", 
                                   paste(colnames(dat)[!(colnames(dat) %in% c("SampleID","Y","Mediator","Treat"))], collapse = " + "),
                                   sep = "") ), data = dat)
    
    summary = summary(mediate(model.m,model.y,treat="Treat",mediator="Mediator",boot=F,sims=1000))
    #capture.output(summary,file="mediator_out.txt",append=FALSE)
    res <- capture.output(summary,append=FALSE)
    
    #sub( "^()\\s", "\\1", res[7])
    tmp <-  base::strsplit(res[grep("ACME",res)],"\\s")[[1]]
    tmp <- tmp[tmp != "" & tmp!="."]
    tmp <- tmp[!grepl("*",tmp,fixed = T)] # remove stars in case the last element is star
    ACME.p <- tmp[length(tmp)]
    
    tmp <- base::strsplit(res[grep("ADE",res)],"\\s")[[1]]
    tmp <- tmp[tmp != "" & tmp!="."]
    tmp <- tmp[!grepl("*",tmp,fixed = T) ]
    ADE.p <- tmp[length(tmp)]
    
    tmp <- base::strsplit(res[grep("Prop. Mediated", res)],"\\s")[[1]]
    tmp <- tmp[tmp != "" ]
    i_str = which(grepl("Mediated", tmp))
    prop.mediated <- tmp[(i_str + 1)]
    
    spearman.r = cor(dat$Treat, dat$Mediator, method = "spearman")
    
    
    
    vec = c(id, ACME.p, ADE.p, prop.mediated, spearman.r)
    names(vec) <- c("Treat_Mediator_Y", "ACME.p", "ADE.p", "prop.mediated", "spearman.r")
    
    Mediation.results <- bind_rows(Mediation.results, vec)
    
  }
}

cat('Mediation analysis generates a result file named Treat_affects_Y_through_Mediator.txt. \n ',  file=log.file, append=T)
write.table(Mediation.results, file = paste(outputDir, "/3_", Treat.omic,"_affects_", Y, "_through_", Mediator.omic,".txt",sep = ""),
            quote = F, row.names = F, sep = "\t")

