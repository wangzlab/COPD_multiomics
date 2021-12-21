## 1. Metabolome

The raw data from mass spectrometer was imported into commercial software Progenesis QI (version 2.2, hereinafter referred to as QI) for peak picking (https://www.nonlinear.com/progenesis/qi/), to obtain information of metabolites such as mass over charge, retention time and ion area. The QI workflow consists of the following steps: peak alignment, peak picking, and peak identification. 

Pre-processing of raw peak data was performed using metaX (https://www.bioconductor.org/packages/3.2/bioc/html/metaX.html), the steps include: 

- Filtering out low quality ions (first removed ions in QC sample that contain over 50% missing value, then removed ions in actual samples that contain over 80% missing value)
- Using k-nearest neighbor (KNN) method for filling the missing values
- Using probabilistic quotient normalization (PQN) method for data normalization
- Using QC-RSC (Quality control-based robust LOESS signal correction) method to alleviate the effects of peak area attenuation
- Filtering out ions in all QC samples which are RSD > 30% (the ions with RSD > 30% are fluctuate greatly in the experiment and will not be included in downstream statistical analysis)

Taken the analysis of positive ion mode for example:

```R
library(metaX)
para <- new("metaXpara")
pfile <- "m_pos.csv" ## Output from QI
sfile <- "s_pos.list" ## Output from QI
idres <- "i_pos.csv" ## Output from QI
para@outdir <- "metaX_result_pos"
para@prefix <- "pos"
para@sampleListFile <- sfile
para@ratioPairs <- "COPD:Healthy"
para <- importDataFromQI(para, file=pfile)
plsdaPara <- new("plsDAPara")
plsdaPara@scale = "pareto"
plsdaPara@cpu = 4
plsdaPara@kfold = 3
#plsdaPara@do = FALSE
res <- doQCRLSC(para, cpu=1)
missValueImputeMethod(para)<-"KNN"
p <- metaXpipe(para, plsdaPara=plsdaPara, missValueRatioQC=0.5, missValueRatioSample=0.8, cvFilter=0.3, idres=idres, qcsc=0, scale="pareto", remveOutlier=FALSE, nor.method="pqn", t=1, nor.order = 1, pclean = FALSE, doROC=FALSE)
save(p, file="pos.rda")
sessionInfo()
```

The metabolite identification was performed using the function metaboliteAnnotation() by searching HMDB (v5.0), METLIN (v3.7.1) and KEGG (v96.0) databases. 

The resultant metabolome data are uploaded as metabolome.txt