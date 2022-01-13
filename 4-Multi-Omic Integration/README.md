## 1. Dimensionality reduction

Dimensionality reduction for metagenomic data was performed using ssGSEA program. Two options are:

- Online platform (https://cloud.genepattern.org/gp/pages/index.jsf)
- local R package (https://github.com/broadinstitute/ssGSEA2.0/)

```
Input: metagenome.gct (a z-score normalized abundance table in gctx format with KEGG KO ids in rows and samples in columns), KEGG_modules.gmt (KEGG gene-module mapping file in gmt format)

Script: Rscript ssgsea-cli.R -i metagenome.gct -o 1_metagenome -d KEGG_modules.gmt --minoverlap 2

Output file: 1_metagenome-combined.gct, 1_metagenome-fdr-pvalues.gct, 1_metagenome-parameters.gct, 1_metagenome-scores.gct, in which both 1_metagenome-scores.gct and 1_metagenome-combined.gct contain the module-level metagenome profile.
```

Dimensionality reduction for metabolomic data was performed using WGCNA. 

```
Input: metabolome.txt (a z-score normalized abundance table with compound ID in rows and samples in column)

Script: Rscript 1.WGCNA_metabolomics.r

Output file: 1_metaB.module_assign.txt, 1_metaB.module_eigengene.txt
```

Dimensionality reduction for host transcriptomic data was performed using WGCNA.

```
Input: transcriptome.txt (a vst normalized gene expression table with gene symbol in rows and samples in column)

Script: Rscript 1.WGCNA_transcriptomics.r

Output file: 1_hostT.module_assign.txt, 1_hostT.module_eigengene.txt
```

## 2. Obtain COPD-associated modules

Differential metagenomic modules were obtained by: 1) obtaining effect size of each KOs in association with disease in a general linear model adjusting demographic confounders, 2) ranking the features by their effect sizes, and 3) comparing the ranks of features within or outside each module in a Wilcoxon rank-sum test.

```
Input: metagenome.gct, KEGG_modules.gmt, metadata.txt (containing a column named SampleID and variables of confounders and a column indicating disease state)

Script: Rscript 2.significant.metaG.modules.r

Output: 2_significant_metaG_modules.RData (a R data file containing the list of metagenomic modules associated with COPD)
```

Differential metaB, hostT and hostP modules/features were obtained by associating the modules/features with COPD in a general linear model adjusting demographic confounders.

```
Input: 1_metaB.module_eigengene.txt (module-level metabolome profile output from step 1)

Script: Rscript 2.significant.metaB.modules.r

Output: 2_significant_metaB_modules.RData (a R data file containing the list of metabolomic modules associated with COPD)
```

```
Input: 1_hostT.module_eigengene.txt (module-level metabolome profile output from step 1)

Script: Rscript 2.significant.hostT.modules.r

Output: 2_significant_hostT_modules.RData (a R data file containing the list of metabolomic modules associated with COPD)
```

```
Input: sputum_proteome.txt (feature-level proteome profile)

Script: Rscript 2.significant.hostP.modules.r

Output: 2_significant_hostP_modules.RData (a R data file containing the list of protein features associated with COPD)
```

## 3. Mediation analysis

Mediation analysis was performed sequentially from metaG-metaB, metaB-hostT, hostT-hostP with sputum neutrophil or eosinophil percentage, respectively. Take the analysis for neutrophil as an example:

```
Input: 1_metaG-combined.gct (module-level metagenome profile output from step 1), 1_metaB.module_eigengene.txt, meta.mediation.NEU.txt (containing metadata for the endpoint for the mediation analysis [NEU] and demographic confounders), 2_significant_metaG_modules.RData and 2_significant_metaB_modules.RData (significant metaG and metaB modules output from step 2)

Script: Rscript 3.mediation_metaG.2.metaB.r

Output: MetaG_affects_NEU_through_MetaB.txt (metaG-metaB-NEU mediation analysis results, including P-value and proportion of mediation effects)
```

```
Input: 1_metaB.module_eigengene.txt, 1_hostT.module_eigengene.txt, meta.mediation.NEU.txt, , 2_significant_metaB_modules.RData and 2_significant_hostT_modules.RData (significant metaB and hostT modules output from step 2)

Script: Rscript 3.mediation_metaB.2.hostT.r

Output: MetaB_affects_NEU_through_HostT.txt (metaB-hostT-NEU mediation analysis results, including P-value and proportion of mediation effects)
```

```
Input: 1_hostT.module_eigengene.txt, sputum_proteome.txt, meta.mediation.NEU.txt, , 2_significant_hostT_modules.RData and 2_significant_hostP_modules.RData (significant hostT and hostP modules output from step 2)

Script: Rscript 3.mediation_hostT.2.hostP.r

Output: HostT_affects_NEU_through_HostP.txt (hostT-hostP-NEU mediation analysis results, including P-value and proportion of mediation effects)
```

## 4. Biological links identification

MetaG-MetaB: Biological links were identified if genes in the metaG module were involved in metabolic reaction for the metabolites in the metaB module, or if the metabolites in the metaB module were present in the metaG module, based on KEGG and MetaCyc databases.

MetaB-HostT: Biological links were identified if metabolites in the metaB modules interact with host genes in the hostT modules (activation/inhibition/binding), based on STITCH database.

HostT-HostP: Biological links were identified if the genes coding for the proteins in HostP were present in the hostT modules or its most enriched pathways, based on human pathway databases (i.e. KEGG, Reactome, MetaBase).

