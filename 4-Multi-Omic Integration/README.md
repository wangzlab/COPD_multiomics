## 1. Dimensionality reduction

Dimensionality reduction for metagenomic data was performed using ssGSEA program. Two options are:

- Online platform (https://cloud.genepattern.org/gp/pages/index.jsf)
- local R package (https://github.com/broadinstitute/ssGSEA2.0/)

```
Input file: metagenome.gct (a z-score normalized abundance table in gctx format with KEGG KO ids in rows and samples in columns), KEGG_modules.gmt (KEGG gene-module mapping file in gmt format)

Script: Rscript ssgsea-cli.R -i metagenome.gct -o metagenome -d KEGG_modules.gmt --minoverlap 2

Output file: metagenome-combined.gct, metagenome-fdr-pvalues.gct, metagenome-parameters.gct, metagenome-scores.gct, in which the metagenome-scores.gct contains the module-level metagenomic data table.
```

Dimensionality reduction for metabolomic data was performed using WGCNA. 

```
Metabolome input: metabolome.txt (a z-score normalized abundance table with compound ID in rows and samples in column)

Script: Rscript 1.WGCNA_metabolomics.r

Output file: 1_metaB.module_assign.txt, 1_metaB.module_eigengene.txt
```

Dimensionality reduction for host transcriptomic data was performed using WGCNA.

```
Metabolome input: transcriptome.txt (a vst normalized gene expression table with gene symbol in rows and samples in column)

Script: Rscript 1.WGCNA_transcriptomics.r

Output file: 1_hostT.module_assign.txt, 1_hostT.module_eigengene.txt
```



## 2. Obtain COPD-associated modules

Differential metagenomic modules was obtained by: 1) obtaining effect size of each KOs in association with disease in a linear model adjusting demographic confounders, 2) ranking the features by their effect sizes, and 3) comparing the ranks of features within or outside each module in a Wilcoxon rank-sum test.

```
Input file: metagenome.gct, KEGG_modules.gmt, metadata.txt (containing a column named SampleID and variables of confounders and a column indicating disease state)

Script: Rscript 2.significant.metaG.modules.r

Output: 2_significant_metaG_modules.RData (a R data file containing the list of modules associated with COPD)
```

