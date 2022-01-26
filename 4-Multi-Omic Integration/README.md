## 1. Dimensionality reduction

Dimensionality reduction for metagenomic data was performed using ssGSEA algorithm. Two options are:

- Online platform (https://cloud.genepattern.org/gp/pages/index.jsf)
- local R package (https://github.com/broadinstitute/ssGSEA2.0/)

```
Input: metagenome.gct (a z-score normalized abundance table in gctx format with KEGG KO ids in rows and samples in columns), KEGG_modules.gmt (KEGG gene-module mapping file in gmt format)

Script: Rscript ssgsea-cli.R -i metagenome.gct -o 1_metagenome -d KEGG_modules.gmt --minoverlap 2

Output: 1_metagenome-combined.gct, 1_metagenome-fdr-pvalues.gct, 1_metagenome-parameters.gct, 1_metagenome-scores.gct, in which both 1_metagenome-scores.gct and 1_metagenome-combined.gct contain the module-level metagenome profile.
```

Dimensionality reduction for metabolomic data was performed using WGCNA. 

```
Input: metabolome.txt (a z-score normalized abundance table with compound ID in rows and samples in column)

Script: Rscript 1.WGCNA_metabolomics.r

Output: 1_metaB.module_assign.txt, 1_metaB.module_eigengene.txt
```

Dimensionality reduction for host transcriptomic data was performed using WGCNA.

```
Input: transcriptome.txt (a vst normalized gene expression table with gene symbol in rows and samples in column)

Script: Rscript 1.WGCNA_transcriptomics.r

Output: 1_hostT.module_assign.txt, 1_hostT.module_eigengene.txt
```

## 2. Obtain COPD-associated modules

Differential metagenomic modules were obtained by: 1) obtaining effect size of each KOs in association with disease in a general linear model adjusting demographic confounders, 2) ranking the features by their effect sizes, and 3) comparing the ranks of features within or outside each module in a Wilcoxon rank-sum test. The differential modules were then associated with sputum neutrophil or eosinophil percentages, and assigned as 'NEU' or 'EOS' if specifically significantly correlated with sputum neutrophil or eosinophil percentages.

```
Input: metagenome.gct, KEGG_modules.gmt, metadata.txt (containing a column named SampleID and variables of confounders and a column indicating disease state)

Script: Rscript 2.significant.metaG.modules.r

Output: 2_significant_metaG_modules.RData (a R data file containing the lists of metagenomic modules associated with COPD and specifically correlated with NEU or EOS)
```

Differential MetaB, HostT and HostP modules/features were obtained by associating the modules/features with COPD in a general linear model adjusting demographic confounders. The differential modules were then associated with sputum neutrophil or eosinophil percentages, and assigned as 'NEU' or 'EOS' if specifically significantly correlated with sputum neutrophil or eosinophil percentages.

```
Input: 1_metaB.module_eigengene.txt (module-level metabolome profile output from step 1)

Script: Rscript 2.significant.metaB.modules.r

Output: 2_significant_metaB_modules.RData (a R data file containing the list of metabolomic modules associated with COPD and specifically correlated with NEU or EOS)
```

```
Input: 1_hostT.module_eigengene.txt (module-level metabolome profile output from step 1)

Script: Rscript 2.significant.hostT.modules.r

Output: 2_significant_hostT_modules.RData (a R data file containing the list of metabolomic modules associated with COPD and specifically correlated with NEU or EOS)
```

```
Input: sputum_proteome.txt (feature-level proteome profile)

Script: Rscript 2.significant.hostP.modules.r

Output: 2_significant_hostP_modules.RData (a R data file containing the list of protein features associated with COPD and specifically correlated with NEU or EOS)
```

## 3. Mediation analysis

Mediation analysis was performed sequentially from MetaG-MetaB, MetaB-HostT, HostT-HostP with sputum neutrophil or eosinophil percentage, respectively. Take the analysis for neutrophil as an example:

```
Input: 1_metaG-combined.gct (module-level metagenome profile output from step 1), 1_metaB.module_eigengene.txt, meta.mediation.NEU.txt (containing metadata for the endpoint for the mediation analysis [NEU] and demographic confounders), 2_significant_metaG_modules.RData and 2_significant_metaB_modules.RData (significant metaG and metaB modules output from step 2)

Script: Rscript 3.mediation_metaG.2.metaB.r

Output: 3_MetaG_affects_NEU_through_MetaB.txt (metaG-metaB-NEU mediation analysis results, including P-value and proportion of mediation effects)
```

```
Input: 1_metaB.module_eigengene.txt, 1_hostT.module_eigengene.txt, meta.mediation.NEU.txt, 2_significant_metaB_modules.RData and 2_significant_hostT_modules.RData (significant metaB and hostT modules output from step 2)

Script: Rscript 3.mediation_metaB.2.hostT.r

Output: 3_MetaB_affects_NEU_through_HostT.txt (metaB-hostT-NEU mediation analysis results, including P-value and proportion of mediation effects)
```

```
Input: 1_hostT.module_eigengene.txt, sputum_proteome.txt, meta.mediation.NEU.txt, 2_significant_hostT_modules.RData and 2_significant_hostP_modules.RData (significant hostT and hostP modules output from step 2)

Script: Rscript 3.mediation_hostT.2.hostP.r

Output: 3_HostT_affects_NEU_through_HostP.txt (hostT-hostP-NEU mediation analysis results, including P-value and proportion of mediation effects)
```

## 4. Biological links identification

MetaG-MetaB: Biological links were identified if genes in the metaG module were involved in metabolic reaction for the metabolites in the metaB module, or if the metabolites in the metaB module were present in the metaG module, based on KEGG and MetaCyc databases.

- Generate link information between KEGG gene ortholog (KO) IDs and MetaCyc compound IDs.

```
Input: ko01000.keg (downloaded from KEGG), metacyc_reactions.txt (downloaded and reformatted from MetaCyc)

Script: Rscript 4.KO2cmpd.txt

Output: KO2EC.list.RData (R data file linking KO IDs and EC), EC2CMPD.lists.RData (R data file linking EC IDs and compounds), KO2CMPD.lists.RData (R data file linking KO IDs and compounds)
```

- Generate link information between KOs and the metabolite IDs in the metabolome.txt file.

```
Input: cmpd2metabo.txt (curated matching information between metabolite IDs and MetaCyc compound IDs), KO2CMPD.lists.RData (generated from above step)

Script: Rscript 4.KO2metabo.r

Output: KO2METABO.lists.RData (R data file linking KO IDs to metabolite IDs)
```

- Generate MetaG-MetaB links.

```
Input: 3_MetaG_affects_NEU_through_MetaB.txt (significant MetaG-MetaB links obtained in step 3), KEGG_modules.tab (KEGG Module-KO mapping file), KO2METABO.lists.RData (R data file linking KO IDs to metabolite IDs generated from above step), 1_metaB_module_assign.txt, metabolome.txt, metagenome.gct, Metabo.KEGGModule.match.txt (manually curated file storing KEGG compound ID-module mapping information).

Script: Rscript 4.MetaG.MetaB.link.r

Output: 4_MetaG.MetaB.modules.NEU.linked.txt (linked MetaG-MetaB modules and the MetaG and MetaB features that link the two modules)
```

MetaB-HostT: Biological links were identified if metabolites in the metaB modules interact with host genes in the hostT modules (activation/inhibition/binding), based on STITCH database.

```
Input: 3_MetaB_affects_NEU_through_HostT.txt (significant MetaB-HostT links obtained in step 3), metabolome.txt, 1_metaB.module_assign.txt, transcriptome.txt, 1_hostT.module_assign.txt, metabo2CIDm.txt (compound-CIDm ID mapping file manually curated based on STITCH database) , all_CIDm_targets.txt (CIDm ID-host target gene mapping file manually curated based on STITCH database)

Script: Rscript 4.MetaB.HostT.link.r

Output: 4_MetaB.HostT.modules.NEU.linked.txt (linked MetaB-HostT modules and the MetaB and HostT features that link the two modules)
```

HostT-HostP: Biological links were identified if the genes coding for the proteins in HostP were present in the hostT modules or its most enriched pathways, based on human pathway databases (i.e. KEGG, Reactome, MetaBase).

```
Input: 3_HostT_affects_NEU_through_HostP.txt (significant HostT-HostP links obtained in step 3), transcriptome.txt, sputum_proteome.txt, 1_hostT.module_assign.txt, protein_info.txt, human_pathway.gmt (KEGG/MetaBase pathway-gene ID mapping file)

Script: Rscript 4.HostT.HostP.link.r

Output: 4_HostT.module_HostP.protein.NEU.linked.txt (linked MetaB-HostT modules and the gene/protein ID that links the two modules)
```

## 5. Leave-one-species-out analysis

Leave-one-species-out analysis was performed to identify driver taxa for the MetaG-MetaB associations by recalculating module-level associations with each species (using bin-based or gene-based taxonomy) iteratively excluded one at a time.

- Calculate contribution of each species to MetaG-MetaB correlation.

Prepare LOSO data by 1) aggragating gene-level metagenomic profile to KO-level with genes from each species removed one at a time, 2) repeating step 1 dimensionality reduction for the KO-level profile, and 3) generate files of speciesX-combined.gct for the dimensionality reduced MetaG profile with speciesX left out.

```
Input: 1) 4_MetaG.MetaB.modules.NEU.linked.txt, 2) 1_metaG-combined.gct, 3) 1_metaB.module_eigengene.txt, 4) LOSO_metaG_DR containing all dimensionality reduced MetaG profiles named as speciesX-combined.gct, one for each species excluded

Script: Rscript 5.LOSO.delta.r

Output: 5_LOSO_NEU.delta.spearman.r.txt (which contains the delta spearman's r and the z-score for the MetaG-MetaB pair when each species was iteratively removed)
```

- Calculate contribution of each species to the turnover of KOs between COPD and control.

  Generate metagenomic KO-level profiles with each species removed one at a time.

```
Input: 1) geneDepth.txt (gene-level abundance file with the first column being gene IDs), 2) ko.txt (KO-gene mapping file), 3) bin_membership.txt (gene-bin mapping file), 2) bin_species (species-bin mapping file)

Script: Rscript 5.LOSO.KO.r

Output: 5_LOSO_ko.abund (directory containing all ko.abund_rm.speciesX.gct files, one for each species removed
```

   Calculate z-score for the contribution of each species to the turnover of each KO in COPD versus control.

```
Input: 1) 5_LOSO_ko.abund (output in last step), 2) metagenome.gct (null KO-level profile), 3) metadata.txt (metadata.txt)

Script: Rscript 5.KOSO.KO.zscore.r

Output: 5_LOSO.KO.zscore.txt (the KO by species matrix table with z-scores)
```

## 6. Random forest analysis

Perform random forest analysis using each linked MetaG-MetaB-HostT set to predict sputum neutrophil or eosinophil percentage. Take neutrophil as an example:

```
Input: 1) 4_MetaG.MetaB.modules.linked.txt, 2) 4_MetaB.HostT.modules.linked.txt, 3) meta.mediation.NEU.txt, 4) 1_metaG-combined.gct, 5) 1_metaB.module_eigengene.txt, 6) 1_hostT.module_eigengene.txt

Script: Rscript 6.random_forest.r

Output: 6_NEU_prediction.performance_byLinks.rf.txt (containing RMSE, RSQ and MAE scores for each linked MetaG-MetaB-HostT set)
```
