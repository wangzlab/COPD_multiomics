## 1. Data preparation

Install packages, set the working directory, and load functions. 

Under the directory there needs to be a "function" directory with all the R function files, a "source.data" directory with all the omic quantification data and a "database" directory which contains all the database-related files.

The omic data file include: metagenome.gct, metabolome.txt, transcriptome.txt, sputum_proteome.txt and serum_proteome.txt

The function and database directories should be directly downloaded from this github.

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

setwd("C:/Users/Wang Zhang/Desktop/COPD_multiomics/combine_analysis/pipeline/pipeline") ## reset it according to your workingdir

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
ssGSEA2(input.ds = "source.data/metagenome.gct", 
        gene.set.databases = "source.data/KEGG_modules.gmt",
        output.prefix = "metaG", outputDir = "1_DimReduction",
        min.overlap = 2, weight = 0, statistic = 'area.under.RES', output.score.type = "NES", nperm = 50,
        export.signat.gct = F, 
        par = T, spare.cores = 5)
```

Dimensionality reduction for host transcriptome and metabolome data

```R
wgcna(input.ds = "source.data/transcriptome.txt", output.prefix = "hostT", outputDir = "1_DimReduction")
wgcna(input.ds = "source.data/metabolome.txt", output.prefix = "metaB", outputDir = "1_DimReduction")
```

