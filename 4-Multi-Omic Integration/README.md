## 1. Dimensionality reduction

Dimensionality reduction for metagenomic data was performed using ssGSEA program. Two options are:

- Online platform (https://cloud.genepattern.org/gp/pages/index.jsf)
- local R package (https://github.com/broadinstitute/ssGSEA2.0/)

```
Input file: metagenome.gct, KEGG_modules.gmt

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

Script: Rscript 1.WGCNA_metabolomics.r

Output file: 1_metaB.module_assign.txt, 1_metaB.module_eigengene.txt
```

