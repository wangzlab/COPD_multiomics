## 1. Data preparation

Install packages, set the working directory, and load functions. Each analytical step is written as a R function.

Under the directory there needs to be a "function" directory with all the R function files, a "source.data" directory with all the omic quantification data and a "database" directory which contains all the database-related files.

Place these omic data tables into source.data directory: metagenome.gct, metabolome.txt, transcriptome.txt, sputum_proteome.txt and serum_proteome.txt

Place the meta.txt (all subjects with Disease and confounder columns), meta.mediation.txt (COPD subjects with dependent variable [NEU or EOS] and confounder columns) into source data directory.

The function and database directories can be directly downloaded from this github.

```R
install.packages("pacman")
install.packages("gtools")
install.packages("verification")
install.packages("doParallel")
install.packages("foreach")
install.packages("magrittr")
install.packages("tibble")
install.packages("data.table")
install.packages("dplyr")
install.packages("reshape2")
install.packages("WGCNA")
install.packages("mediation")
install.packages("tidyverse")
install.packages("ranger")

setwd("C:/Users/Wang Zhang/Desktop/COPD_multiomics/combine_analysis/pipeline/pipeline") ## reset it to your workingdir

source("functions/ssGSEA2.0.R")
source("functions/io.R")
source("functions/utils.R")
source("functions/SupportingFunc.R")
source("functions/wgcna.R")
source("functions/glm.sigModules.R")
source("functions/MediationAnalysis-parallel.R")
source("functions/ko2cmpd.R")
source("functions/ko2metabo.R")
source("functions/MetaG.MetaB.link.R")
source("functions/MetaB.HostT.link.R")
source("functions/HostT.HostP.link.R")
source("functions/LOSO.ko.R")
source("functions/LOSO.delta.r.R")
source("functions/rf_MetaG.MetaB.HostT.links.R")

if(!suppressPackageStartupMessages(require("pacman"))){
  install.packages("pacman")
}
library(pacman)
p_load(gtools)
p_load(verification)
p_load(doParallel)
p_load(foreach)
p_load(magrittr)
library(tibble)
library(data.table)
library(dplyr)
library(reshape2)
library(WGCNA)
library(mediation)
library(tidyverse)
library(ranger)
```

## 2. Dimensionality reduction

Dimensionality reduction for metagenome data

```R
ssGSEA2(input.ds = "source.data/metagenome.gct", ## the metagenome.gct is the relative abundance of KOs after arcsin z-score normalization and z-score standardization 
        gene.set.databases = "source.data/KEGG_modules.gmt",
        output.prefix = "metaG", outputDir = "1_DimReduction",
        min.overlap = 2, weight = 0, statistic = 'area.under.RES', output.score.type = "NES", nperm = 1000,
        export.signat.gct = F, 
        par = T, spare.cores = 5)
```

The enrichment scores for all KEGG modules are saved as: 1_DimReduction/metaG-combined.gct 

Dimensionality reduction for host transcriptome and metabolome data

```R
wgcna(input.ds = "source.data/transcriptome.txt", output.prefix = "hostT", outputDir = "1_DimReduction")
wgcna(input.ds = "source.data/metabolome.txt", output.prefix = "metaB", outputDir = "1_DimReduction")
```

The metabolome module eigengenes are saved as: 1_DimReduction/metaB.module_eigengene.txt

The metabolome module assignment are saved as: 1_DimReduction/metaB.module_assign.txt

The transcriptome module eigengenes are saved as: 1_DimReduction/hostT.module_eigengene.txt

The transcriptome module assignment are saved as: 1_DimReduction/hostT.module_assign.txt

## 3. Module association with COPD

```R
# For MetaG modules, run 'MetaG.Module.affectDisease' function, to 
# 1) get effect size of each KOs in association with disease in a linear model, 
# 2) rank the features by this effect size, and 
# 3) compare the ranks of features within or outside each module in a Wilcoxon rank-sum test
MetaG.mod.onDiseaseState <- MetaG.Module.affectDisease(input.ds = "source.data/metagenome.gct", 
                                                       meta.file = "source.data/meta.txt",
                                                       disease.state = "Disease",
                                                       gene.set.databases = "source.data/KEGG_modules.gmt")

MetaG.pval<-MetaG.mod.onDiseaseState$Module.pvalue
MetaG.sigMods <-  MetaG.pval$module[MetaG.pval$Wilcox.p < 0.1,]

# organize MetaB data into a dataframe with rownames being modules and colnames being samples
MetaB.Mod.dat <- fread("1_DimReduction/metaB.module_eigengene.txt",data.table = F) 
#  dplyr::filter(!grepl("#",`#NAME`,fixed=T)) 
#MetaB.Mod.dat[-1] <- sapply(MetaB.Mod.dat[-1], as.numeric)
MetaB.Mod.dat <- MetaB.Mod.dat %>% tibble::column_to_rownames("V1") %>% t() %>% data.frame()
# then identify metaB significant modules ------
MetaB.sigMods <- glm.sigModules(input.ds = MetaB.Mod.dat,
                                meta.file="source.data/meta.txt", 
                                glm.family = "binomial",
                                glm.p = 0.1)
# organize HostT data into a dataframe with rownames being modules and colnames being samples
HostT.Mod.dat <- fread("1_DimReduction/hostT.module_eigengene.txt",data.table = F)
#  dplyr::filter(!grepl("#",`#NAME`,fixed=T)) 
#HostT.Mod.dat[-1] <- sapply(HostT.Mod.dat[-1], as.numeric)
HostT.Mod.dat <- HostT.Mod.dat %>% tibble::column_to_rownames("V1") %>% t() %>% data.frame()
# then identify HostT significant modules ------
HostT.sigMods <- glm.sigModules(input.ds = HostT.Mod.dat,
                                meta.file="source.data/meta.txt",
                                glm.family = "binomial",
                                glm.p = 0.1)
# identify HostP significant modules -------
HostP.data <- data.frame(fread("source.data/sputum_cyto.txt"), row.names = 1)
HostP.sigFeatures <- glm.sigModules(input.ds = HostP.data,
                                    meta.file="source.data/meta.txt",
                                    glm.family = "binomial",
                                    glm.p = 0.1)
```

The significant module/feature IDs are stored in MetaG.sigMods, MetaB.sigMods, HostT.sigMods and HostP.sigFeatures.

## 4. Mediation analysis

Take mediation analysis for COPD neutrophilic inflammation as an example.

```R
# MetaG modules affect NEU through MetaB modules  -----------

MediationAnalysis_parallel(Treat.omic = "MetaG", Mediator.omic = "MetaB", 
                           Treat.omic.input = "1_DimReduction/metaG-combined.gct",
                           Mediator.omic.input = MetaB.Mod.dat,
                           meta.mediate = "source.data/meta.mediation.txt", Y = "NEU",
                           Treat.omic.sigModules = MetaG.sigMods,
                           Mediator.omic.sigModules = MetaB.sigMods,
                           log.file = "mediation.parallel.log",
                           outputDir = "2_Mediation",
                           threads = 25)
MetaG.MetaB.NEU_medres <- fread("2_Mediation/MetaG_affects_NEU_through_MetaB_parallel.txt")


# MetaB modules affect NEU through HostT modules -----------

MediationAnalysis_parallel(Treat.omic = "MetaB", Mediator.omic = "HostT", 
                           Treat.omic.input = MetaB.Mod.dat,
                           Mediator.omic.input = HostT.Mod.dat,
                           meta.mediate = "source.data/meta.mediation.txt", Y = "NEU",
                           Treat.omic.sigModules = MetaB.sigMods,
                           Mediator.omic.sigModules = HostT.sigMods,
                           log.file = "mediation.parallel.log",
                           outputDir = "2_Mediation",
                           threads = 25)
MetaB.HostT.NEU_medres <- fread("2_Mediation/MetaB_affects_NEU_through_HostT_parallel.txt")


# HostT modules affect NEU through HostP features ----------- 

MediationAnalysis_parallel(Treat.omic = "HostT", Mediator.omic = "HostP", 
                           Treat.omic.input = HostT.Mod.dat,
                           Mediator.omic.input = HostP.data,
                           meta.mediate = "source.data/meta.mediation.NEU.txt", Y = "NEU",
                           Treat.omic.sigModules = HostT.sigMods,
                           Mediator.omic.sigModules = HostP.sigFeatures,
                           outputDir = "2_Mediation",
                           threads = 25)
HostT.HostP.NEU_medres <- fread("2_Mediation/HostT_affects_NEU_through_HostP_parallel.txt")
```

Three files should be saved in the "2_Mediation" folder: MetaG_affects_NEU_through_MetaB_parallel.txt, MetaB_affects_NEU_through_HostT_parallel.txt, HostT_affects_NEU_through_HostP_parallel.txt

## 5. Biological links identification

```R
ko2cmpd(dbDir = "database")

#' The ko2cmpd function generates information on substrates, products and substrate&products (in the format of metacyc compound names) for each KO number. 
#' Users should provide the path to a database directory, which must contain the belowing files:
#' ko01000.keg; metacyc_reactions.txt (uploaded in this github)
#' This function should be run when the user goes through the pipeline for the first time and whenever the user updated the "ko01000.keg" and "metacyc_reactions.txt" files. 
#' Intermediate and resulting files in the format of RData will be stored in the same database directory.


ko2metabo(dbDir = "database")
#' The ko2metabo function generates information on substrates, products and substrate&products (in the format of metabolomic ids) for each KO number. 
#' Users should provide the path to a database directory, which must contain the cmpd2metabo.txt file and KO2CMPD.lists.RData file (uploaded in this github)
#' This function should be run after the ko2cmpd() function
#' Resulting files in the format of RData will be stored in the same database directory.

# The MetaG.MetaB.link function allows you to identify MetaG - MetaB module pairs with biological links. 
# user should provide a metabo.KEGGmodule.match_file

MetaG.MetaB.links <- MetaG.MetaB.link( MetaG.MetaB.NEU_medres, 
                                       MetaG_module.feature_file =  "database/KEGG_modules.tab", 
                                       KO2METABO_file = "database/KO2METABO.lists.RData", ## this is an output
                                       MetaB_module.feature_file = "1_DimReduction/metaB.module_assign.txt",
                                       MetaB_quantity_file = "source.data/metabolome.txt",
                                       MetaG_quantity_file = "source.data/metagenome.gct",
                                       metabo.KEGGmodule.match_file = "database/Metabo.KEGGModule.match.txt",
                                       output.dir = "3_Biological_Links",
                                       ACME.p.co = 0.1)

# The  MetaB.HostT.links  function allows you to identify MetaB - HostT module pairs with biological links.  
# User prepares a  metabo2CIDm file using the perl scripts provided on github.

MetaB.HostT.links <- MetaB.HostT.link(mediation.res = MetaB.HostT.NEU_medres,
                                      MetaB_quantity_file = "source.data/metabolome.txt",  
                                      MetaB_module.feature_file = "1_DimReduction/metaB.module_assign.txt", 
                                      HostT_quantity_file = "source.data/transcriptome.txt", 
                                      HostT_module.feature_file = "1_DimReduction/hostT.module_assign.txt", 
                                      METABO2CIDm_file =  "database/metabo2CIDm.txt",  
                                      CIDm.receptor_file = "database/all_cidm_receptor.txt",  
                                      output.dir = "3_Biological_Links",
                                      ACME.p.co = 0.1) 
# hostT-hostP links ------------
# The HostT.HostP.link function allows you to identify biological linkes between HostT modules and HostP features.
# User should provide a protein.gene_file and a Pthway2Gene_file. 

HostT.HostP.links <- HostT.HostP.link(mediation.res = HostT.HostP.NEU_medres,
                                      HostT_quantity_input = "source.data/transcriptome.txt",
                                      HostP_quantity_input = HostP.data,
                                      HostT_module.feature_file = "1_DimReduction/hostT.module_assign.txt",
                                      HostP_protein.gene_file = "database/protein_info.txt",
                                      Pthway2Gene_file = "database/pathway.gmt",
                                      output.dir = "3_Biological_Links",
                                      ACME.p.co = 0.1)
```


Three files should be saved in the 3_Biological_Links folder: MetaG.MetaB.modules.linked.txt, MetaB.HostT.modules.linked.txt, HostT.HostP.modules.linked.txt.

```R
# first obtain the concatenated MetaG-MetaB-HostT links ---------

MetaG.MetaB.links <- fread("3_Biological_Links/MetaG.MetaB.modules.linked.txt",select = 1, col.names = "V1") %>% unique() %>%
  mutate(MetaG.module = sapply(strsplit(V1, "_", fixed = T), "[[", 1),
         MetaB.module = sapply(strsplit(V1, "_", fixed = T), "[[", 2))


MetaB.HostT.links <- fread("3_Biological_Links/MetaB.HostT.modules.linked.txt",select = 1, col.names = "V1")  %>% 
  unique()  %>% # not unique !! need to check link script
  mutate(MetaB.module = sapply(strsplit(V1, "_", fixed = T), "[[", 1),
         HostT.module = sapply(strsplit(V1, "_", fixed = T), "[[", 2))

MetaG.MetaB.HostT.links <- NULL
for(i in c(1:nrow(MetaG.MetaB.links)) ){
  
  gb.pair = MetaG.MetaB.links$V1[i]
  
  bm = strsplit(gb.pair, "_", fixed = T)[[1]][2]
  
  bt.pairs <- MetaB.HostT.links$V1[which(MetaB.HostT.links$MetaB.module == bm)] 
  
  tmp <- expand.grid(gb.pair, bt.pairs)
  
  MetaG.MetaB.HostT.links <- bind_rows(MetaG.MetaB.HostT.links, tmp)
  
} 

MetaG.MetaB.HostT.links <- MetaG.MetaB.HostT.links %>% mutate(Var1 = as.character(Var1), Var2 = as.character(Var2)) %>%
  mutate(MetaG = sapply( strsplit(Var1,"_",fixed = T), "[[", 1) ,
         MetaB = sapply( strsplit(Var1,"_",fixed = T), "[[", 2) , 
         HostT = sapply( strsplit(Var2,"_",fixed = T), "[[", 2))


# then create a predicted variable data frame ----- 
meta <- fread("source.data/meta.mediation.txt") %>% select(SampleID, NEU) 
#meta <- fread("source.data/meta.txt") %>% select(SampleID, Disease) #%>% mutate(Y = as.factor(as.character(Disease)))

# Get Guangzhou sample
GZ.sp <-meta$SampleID[!grepl("^Z", meta$SampleID)] 

# perform random forest analysis -----------
#' The rf_MetaG.MetaB.HostT.Links function allows you to evaluate the accuracy of predicting e.g. a clinical variable by the MetaG.MetaB.HostT links through random forest modeling.  
#' This function generates a performance table recording prediction performance of each MetaG.MetaB.HostT link, as rsme and rsq for "regression", 
#' and as accuracy and roc_auc for "classification".   

# example of usage: ####

rf.Performance_cls <- rf_MetaG.MetaB.HostT.Links(MetaG.Mod.input = "1_DimReduction/metaG-combined.gct",
                                                 MetaB.Mod.input = MetaB.Mod.dat,  
                                                 HostT.Mod.input = HostT.Mod.dat,  
                                                 MetaG.MetaB.HostT.link.df = MetaG.MetaB.HostT.links,  
                                                 PredictedVar.input = meta,   
                                                 PredictionType = "regression",  
                                                 output.dir = "4_RandomForest",
                                                 Training.samples = GZ.sp) 
```

The output of the random forest analysis should be saved in 4_RandomForest/prediction.performance_byLinks.rf.txt

## 6. Driver taxa analysis

```R
# 1) generates KO abundance files with each species iteratively excluded

#' The LOSO.ko function allows you to perform LOSO (leave one species out) analysis. 
#' Users should provide 1) gene quantification data, 2) information about Species annotaion of bins, 3) connection between bins and genes through scaffold ids, and 4) KO annotaion of genes. 
#' This function generates one KO abundance file each time after excluding one Species from the gene quantification data.  

LOSO.ko(geneDepth_file = "source.data/geneDepth.txt",  
        gene.ko_file = "source.data/ko_noeuk.txt",
        bin.scaffold_file =  "source.data/all_membership.txt",
        bin.species_file = "source.data/bin_species.txt") # memory issue for large geneDepth file, can run on linux server

# 2) perform ssGSEA for each ko abundance file
if(!dir.exists("5_LOSO")) dir.create("5_LOSO")

koFiles <- list.files("5_LOSO/", full.names = T)
for(kof in koFiles){
  specs <-sub("\\.gct$", "", sub("^ko\\.abund_rm\\.", "", basename(kof)) ) 
  
  ssGSEA2(input.ds = kof, 
          gene.set.databases = "source.data/KEGG_modules.gmt",
          output.prefix = specs,
          min.overlap = 2, weight = 0, statistic = 'area.under.RES', output.score.type = "NES", nperm = 100,
          par = T,export.signat.gct = F,
          outputDir = "5_LOSO")
  
}

# 3) calculate delta spearman.r by each species
#' The LOSO.delta.r function calculate delta spearman r in the LOSO analysis. 
#' Users provide an original quantification input for MetaG modules without leaving any species out (MetaG.mod.before);   
#' a directory path containing all module quantification files after removing one species each time (MetaG.mod.after.dir);  
#' a MetaB module quantification input (MetaB.mod.input); 
#' pairs of MetaG-MetaB modules on which the impacted by the species are of interest (module.pairs).    
#' The script calculates the spearman r values between the interested MetaG-MetaB module pairs before and after removing one Species from the metagenomic data each time. 
#' A delta spearman r is calculated for each species. The delta spearman r values are z-score standardized across each MetaG-MetaB module pair. 

# create a test module pair data frame:
test.ModulePairs = fread("biological.links/MetaG.MetaB.modules.linked.txt", select = 1, data.table = F) %>% 
  mutate(MetaG.module = sapply(strsplit(MetaG.MetaB_modulePair,"_", fixed = T),"[[", 1)) %>%
  mutate(MetaB.module = sapply(strsplit(MetaG.MetaB_modulePair,"_", fixed = T),"[[", 2))


LOSO.results <- LOSO.delta.r(module.pairs = test.ModulePairs,
                             MetaG.mod.before = "1_DimReduction/metaG-combined.gct",
                             MetaG.mod.after.dir = "5_LOSO",
                             MetaB.mod.input = MetaB.Mod.dat,
                             output.dir = "5_LOSO"
                             output.prefix = "metaG_metaB")  
```

The KO-level metagenome gct files with each species iteratively removed are saved in 5_LOSO/ko.abund_rm.$species.gct

The module-level metagenome gct files with each species iteratively removed are saved in 5_LOSO/$species-combined.gct

The delta value of spearman's rho with each species iteratively removed and z-score contribution of each species to the MetaG-MetaB correlation are saved in 5_LOSO/metaG_metaB_LOSO_delta.spearman.r.txt
