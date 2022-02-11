## 1. Datasets and scripts

- All required datasets are under 'datasets' or '5-Datasets' directories. For large data files that cannot be uploaded to GitHub, they are under FigShare doi:10.6084/m9.figshare.19126814.

- All scripts are under 'scripts' directory. For each step, follow the corresponding R (or perl) script (see below).

- Install required R packages:

  ```
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
  ```

  

## 2. Dimensionality reduction

Dimensionality reduction for metagenomic data was performed using ssGSEA algorithm. Two options are:

- Online platform (https://cloud.genepattern.org/gp/pages/index.jsf)
- local R package (https://github.com/broadinstitute/ssGSEA2.0/)

```
Input: metagenome.gct (a z-score normalized abundance table in gctx format with KEGG KO ids in rows and samples in columns), KEGG_modules.gmt (KEGG gene-module mapping file in gmt format)

Script: Rscript ssgsea-cli.R -i metagenome.gct -o 1_metagenome -d KEGG_modules.gmt --minoverlap 2

Output: 1_metagenome-combined.gct, 1_metagenome-fdr-pvalues.gct, 1_metagenome-parameters.gct, 1_metagenome-scores.gct (Both 1_metagenome-scores.gct and 1_metagenome-combined.gct contain the module-level metagenome profile)
```

Dimensionality reduction for metabolomic data was performed using WGCNA. 

```
Input: metabolome.txt (a z-score normalized abundance table with compound ID in rows and samples in column)

Script: 1.WGCNA_metabolomics.r

Output: 1_metaB.module_assign.txt, 1_metaB.module_eigengene.txt
```

Dimensionality reduction for host transcriptomic data was performed using WGCNA.

```
Input: transcriptome.txt (a vst normalized gene expression table with gene symbol in rows and samples in column)

Script: 1.WGCNA_transcriptomics.r

Output: 1_hostT.module_assign.txt, 1_hostT.module_eigengene.txt
```

The dimensionality reduced multi-omic data profiles (1_metagenome-combined.gct, 1_metaB.module_eigengene.txt, 1_hostT.module_eigengene.txt) was generated to be used for downstream analyses.

The number of modules generated:

|       | Features | Modules     |
| ----- | -------- | ----------- |
| MetaG | 6678     | 461         |
| MetaB | 1671     | 128 (n>=5)  |
| HostT | 19142    | 497 (n>=10) |



## 3. Obtain neutrophil or eosinophil-associated modules

Differential metagenomic modules were obtained by: 1) obtaining effect size of each KOs in association with disease in a general linear model adjusting demographic confounders, 2) ranking the features by their effect sizes, and 3) comparing the ranks of features within or outside each module in a Wilcoxon rank-sum test. The differential modules were then associated with sputum neutrophil or eosinophil percentages, and assigned as 'NEU' or 'EOS' if specifically significantly correlated with sputum neutrophil or eosinophil percentages.

```
Input: metagenome.gct, KEGG_modules.gmt, metadata.txt (containing a column named SampleID, variables of confounders, and a column indicating disease state), meta.mediation.NEU.txt or meta.mediation.EOS.txt (containing a column named SampleID, variables of confounders, a column indicating NEU or EOS)

Script: 2.significant.metaG.modules.r

Output: 2_significant_metaG_modules.RData (a R data file containing the lists of metagenomic modules associated with COPD and specifically correlated with NEU or EOS)
```

Differential MetaB, HostT and HostP modules/features were obtained by associating the modules/features with COPD in a general linear model adjusting demographic confounders. The differential modules were then associated with sputum neutrophil or eosinophil percentages, and assigned as 'NEU' or 'EOS' if specifically significantly correlated with sputum neutrophil or eosinophil percentages.

```
Input: 1_metaB.module_eigengene.txt (module-level metabolome profile output from step 1), metadata.txt, meta.mediation.NEU.txt or meta.mediation.EOS.txt

Script: 2.significant.metaB.modules.r

Output: 2_significant_metaB_modules.RData (a R data file containing the list of metabolomic modules associated with COPD and specifically correlated with NEU or EOS)
```

```
Input: 1_hostT.module_eigengene.txt (module-level transcriptome profile output from step 1), metadata.txt, meta.mediation.NEU.txt or meta.mediation.EOS.txt
    
Script: 2.significant.hostT.modules.r

Output: 2_significant_hostT_modules.RData (a R data file containing the list of metabolomic modules associated with COPD and specifically correlated with NEU or EOS)
```

```
Input: sputum_proteome.txt (feature-level proteome profile), metadata.txt, meta.mediation.NEU.txt or meta.mediation.EOS.txt

Script: 2.significant.hostP.modules.r

Output: 2_significant_hostP_modules.RData (a R data file containing the list of protein features associated with COPD and specifically correlated with NEU or EOS)
```

This step generates MetaG, MetaB, HostT modules and HostP features differentially abundant between COPD and healthy controls, and specifically associated with sputum neutrophilic or eosinophilic inflammation. 

The number of modules/features retained:

| Module/Feature | NEU                 | EOS  |
| :------------: | ------------------- | ---- |
|     MetaG      | 50 (24 with P<0.01) | 19   |
|     MetaB      | 38                  | 16   |
|     HostT      | 97                  | 29   |
|     HostP      | 21                  | 23   |



## 4. Mediation analysis

Mediation analysis was performed sequentially from MetaG-MetaB, MetaB-HostT, HostT-HostP with sputum neutrophil or eosinophil percentage, respectively. Take the analysis for neutrophil as an example:

- MetaG-MetaB-NEU: This step identifies MetaG-MetaB module pairs in which MetaG affects NEU through the mediation of MetaB.

```
Input: 1_metaG-combined.gct (module-level metagenome profile output from step 1), 1_metaB.module_eigengene.txt, meta.mediation.NEU.txt, 2_significant_metaG_modules.RData (significant metaG modules output from step 2), 2_significant_metaB_modules.RData (significant metaB modules output from step 2)

Script: 3.mediation_metaG.2.metaB.r

Output: 3_MetaG_affects_NEU_through_MetaB.txt (metaG-metaB-NEU mediation analysis results, including P-value and proportion of mediation effects)
```

- MetaB-HostT-NEU: This step identifies MetaB-HostT module pairs in which MetaB affects NEU through the mediation of HostT.

```
Input: 1_metaB.module_eigengene.txt, 1_hostT.module_eigengene.txt, meta.mediation.NEU.txt, 2_significant_metaB_modules.RData (significant metaB modules output from step 2), 2_significant_hostT_modules.RData (significant hostT modules output from step 2)

Script: 3.mediation_metaB.2.hostT.r

Output: 3_MetaB_affects_NEU_through_HostT.txt (metaB-hostT-NEU mediation analysis results, including P-value and proportion of mediation effects)
```

- HostT-HostP-NEU: This step identifies HostT-HostP module/feature pairs in which HostT affects NEU through the mediation of HostP.

```
Input: 1_hostT.module_eigengene.txt, sputum_proteome.txt, meta.mediation.NEU.txt, 2_significant_hostT_modules.RData (significant hostT modules output from step 2), 2_significant_hostP_features.RData (significant hostP features output from step 2)

Script: 3.mediation_hostT.2.hostP.r

Output: 3_HostT_affects_NEU_through_HostP.txt (hostT-hostP-NEU mediation analysis results, including P-value and proportion of mediation effects)
```

The number of links generated:

| Links       | NEU  | EOS  |
| ----------- | ---- | ---- |
| MetaG-MetaB | 428  | 304  |
| MetaB-HostT | 3396 | 464  |
| HostT-HostP | 1148 | 155  |

## 5. Biological links identification

MetaG-MetaB: Biological links were identified if genes in the MetaG module were involved in metabolic reaction for the metabolites in the MetaB module, or if the metabolites in the MetaB module were present in the MetaG module, based on KEGG and MetaCyc databases.

- Generate link information between KEGG gene ortholog (KO) IDs and MetaCyc compound IDs. 
- Before running this, we downloaded the MetaCyc database (https://metacyc.org/), which includes a file named 'reactions.dat' containing detailed metabolic reaction information (i.e. substrate, product, enzyme etc). We converted this file to a tab delimited one for future use (metacyc_reactions.txt).  
- Then we run this step, which first links KO to EC based on ko01000.keg from KEGG database, and then links EC/KO to metabolites, based on metabolic reaction information derived from MetaCyc database (that a compound is a substrate or a product of a metabolic reaction catalyzed by the protein with the corresponding EC/KO).

```
Input: ko01000.keg (downloaded from KEGG), metacyc_reactions.txt (downloaded and reformatted from MetaCyc)

Script: 4.KO2cmpd.r

Output: KO2EC.list.RData (R data file linking KO IDs and EC), EC2CMPD.lists.RData (R data file linking EC IDs and compounds), KO2CMPD.lists.RData (R data file linking KO IDs and compounds)
```

- Generate link information between KOs and the metabolite IDs in the metabolome.txt file. 
- Before running this, we first created a file matching the IDs from metabolome data and MetaCyc compounds, this was done through 1) ID conversion (i.e. use the Compound ID Conversion utility in metaboanalyst https://www.metaboanalyst.ca/), and 2) using 100% compound structure match by open babel (https://openbabel.org/wiki/Main_Page). Then we run this step to link KO IDs to the metabolite IDs in the metabolomic data.
- The MetaCyc file of interest in this step was 'compounds.dat' which we similarly converted to a tab delimited file. Then we extracted IDs (PubChem, ChEBI, KEGG, HMDB) for each compound and mapped against the metabolite IDs in metabolomic data. We used a perl script to do these which was uploaded as '4_cmpd2metabo.pl'.
- For the remaining metabolites in the metabolomic data that cannot be mapped to MetaCyc by ID mapping, we extracted their SMILE structure and searched against the compounds in MetaCyc (100% structural match), to maximize the mapping information between metabolites and MetaCyc compounds.

```
1. Generate metabolite sdf file from SMILES string list for MetaCyc compounds:
   obabel metacyc.smi –O metacyc.sdf --gen3D
2. Generate fastsearch index for the metacyc.sdf file:
   obabel metacyc.sdf -ofs 
3. For each metabolite in metabolomic data generate its own sdf file:
   obabel metabo_$i.smi -O metabo_$i.sdf --gen3D
4. Search metabolites against metacyc.fs index to get top 10 hits:
   obabel metacyc.fs –otxt –s metabo_$i.sdf –at 10 –aa > metabo_$i.results
```

```
Input: cmpd2metabo.txt (curated matching information between metabolite IDs and MetaCyc compound IDs), KO2CMPD.lists.RData (generated from above step)

Script: 4.KO2metabo.r

Output: KO2METABO.lists.RData (R data file linking KO IDs to metabolite IDs)
```

- Generate MetaG-MetaB links. 
- Before running this, we curated a file called 'Metabo.KEGGModule.match.txt' which stores information for the modules to which the KEGG compounds (C number) belong to. Then we run this step, which considers two main situations to link MetaG and MetaB modules: 1) the KO from MetaG and metabolite from MetaB are directly linked based on the above linkage data; 2) the KO and the metabolite belong to the same KEGG module (based on the 'Metabo.KEGGModule.match.txt' file).

```
Input: 3_MetaG_affects_NEU_through_MetaB.txt (significant MetaG-MetaB links obtained in step 3), KEGG_modules.tab (KEGG Module-KO mapping file), KO2METABO.lists.RData (R data file linking KO IDs to metabolite IDs generated from above step), 1_metaB_module_assign.txt, metabolome.txt, metagenome.gct, Metabo.KEGGModule.match.txt (file storing KEGG compound ID-module mapping information)

Script: 4.MetaG.MetaB.link.r

Output: 4_MetaG.MetaB.modules.NEU.linked.txt (linked MetaG-MetaB modules and the MetaG and MetaB features that link the two modules)
```

- MetaB-HostT: Biological links were identified if metabolites in the metaB modules interact with host genes in the hostT modules (activation/inhibition/binding), based on STITCH database. 
- Before running this, we first downloaded the STITCH database (http://stitch.embl.de/cgi/download.pl?UserId=LXoSgRCY5mtL&sessionId=zCuovnNcZ0QK). We mainly looked for three files (choose organism human): 
  - 9606.protein_chemical.links.detailed.v5.0.tsv which contains interactions between compounds (CIDm, CIDs initial) and human genes (Ensembl ID) and interaction scores;
  - 9606.actions.v5.0.tsv which contains modes of interaction (activate, inhibit, binding, catalysis, reaction);
  - chemical.sources.v5.0.tsv which contains compound ID mapping rules (# ChEBI, ChEMBL, KEGG, PC (PubChem Compound) and PS (PubChem Substance)).
  - Note: In STITCH, CIDm IDs are the merged compound IDs for different CIDs, so only CIDm needs to be considered when mapping to host targets.
- Then we mapped metabolite IDs from metabolomic data to CIDm IDs based on above mapping rules, which generated a file named 'metabo2CIDm.txt'.
- Then we extracted compound-target interaction information from '9606.protein_chemical.links.detailed.v5.0.tsv' and '9606.actions.v5.0.tsv' and stored it in 'all_CIDm_targets.txt' (available at doi:10.6084/m9.figshare.19126814).
- We used a perl script to do these which was uploaded as '4_metabolome2stitch.pl'.

```
Input: 3_MetaB_affects_NEU_through_HostT.txt (significant MetaB-HostT links obtained in step 3), metabolome.txt, 1_metaB.module_assign.txt, transcriptome.txt, 1_hostT.module_assign.txt, metabo2CIDm.txt (compound-CIDm ID mapping file manually curated based on STITCH database), all_CIDm_targets.txt (CIDm ID-host target gene mapping file manually curated based on STITCH database)

Script: 4.MetaB.HostT.link.r

Output: 4_MetaB.HostT.modules.NEU.linked.txt (linked MetaB-HostT modules and the MetaB and HostT features that link the two modules)
```

- HostT-HostP: Biological links were identified if the genes coding for the proteins in HostP were present in the hostT modules or its most enriched pathways, based on human pathway databases (i.e. KEGG, Reactome, MetaBase).
- Protein-Gene ID match information are stored in 'protein_info.txt'. Gene-Pathway match information are stored in 'human_pathways.gmt'.

```
Input: 3_HostT_affects_NEU_through_HostP.txt (significant HostT-HostP links obtained in step 3), transcriptome.txt, sputum_proteome.txt, 1_hostT.module_assign.txt, protein_info.txt, human_pathway.gmt (KEGG/MetaBase pathway-gene ID mapping file)

Script: 4.HostT.HostP.link.r

Output: 4_HostT.module_HostP.protein.NEU.linked.txt (linked MetaB-HostT modules and the gene/protein ID that links the two modules)
```

The number of links that were biologically associated:

| Links       | NEU  | EOS  |
| ----------- | ---- | ---- |
| MetaG-MetaB | 109  | 26   |
| MetaB-HostT | 335  | 58   |
| HostT-HostP | 135  | 22   |

## 6. Leave-one-species-out analysis

Leave-one-species-out analysis was performed to identify driver taxa for the MetaG-MetaB associations by recalculating module-level associations with each species (using bin-based or gene-based taxonomy) iteratively excluded one at a time.

- Prepare LOSO data by 
  - aggragating gene-level metagenomic profiles to KO-level with genes from each species (a total of 68 species based on binning results) removed one at a time,
  - repeating step 1 dimensionality reduction for the KO-level profiles, and,
  - generating files of speciesX-combined.gct for the dimensionality reduced MetaG profiles with speciesX left out.

- Generate metagenomic KO-level profiles with each species removed one at a time.

```
Input: geneDepth.txt (gene-level abundance file with the first column being gene IDs), ko.txt (KO-gene mapping file), bin_membership.txt (gene-bin mapping file), bin_species.txt (species-bin mapping file)

Script: 5.LOSO.KO.r

Output: 5_LOSO_ko.abund (directory containing all ko.abund_rm.speciesX.gct files, one for each species removed)
```

- Calculate contribution of each species (△Spearman's rho by excluding that species) to the MetaG-MetaB correlation of interest.

```
Input: 4_MetaG.MetaB.modules.NEU.linked.txt (the linked MetaG-MetaB modules obtained in step 4), 1_metaG-combined.gct, 1_metaB.module_eigengene.txt, LOSO_metaG_DR (containing all dimensionality reduced MetaG profiles named as speciesX-combined.gct, one for each species excluded)

Script: 5.LOSO.delta.r

Output: 5_LOSO_NEU.delta.spearman.r.txt (containing the delta spearman's rho when each species was iteratively removed)
```

- Calculate z-score for the contribution of each species to the turnover of each KO in COPD versus control.
- The z-score was calculated as the deviation of the fold-change of that KO in COPD versus controls from its original fold-change when a species was excluded, normalized by standard deviation of the fold-changes of that KO when all species were iteratively excluded.

```
Input: 5_LOSO_ko.abund (output in last step), metagenome.gct (null KO-level profile), metadata.txt

Script: 5.KOSO.KO.zscore.r

Output: 5_LOSO.KO.zscore.txt (the KO by species matrix table with z-scores)
```

## 7. Random forest analysis

Perform random forest analysis using each linked MetaG-MetaB-HostT set to predict sputum neutrophil or eosinophil percentage. Take neutrophil as an example:

- We first aggregated MetaG-MetaB and MetaB-HostT links to the full-path of MetaG-MetaB-HostT, by linking up 'KO-metabolite-host gene' feature-level information.

The number of links constituting the full paths:

|      | MetaG-MetaB | MetaB-HostT | MetaG-MetaB-HostT |
| ---- | ----------- | ----------- | ----------------- |
| NEU  | 66          | 136         | 620               |
| EOS  | 17          | 38          | 134               |



- Then we performed a random forest regression between each linked set of MetaG-MetaB-HostT modules and NEU or EOS, and ranked the module sets by model performance scores.

```
Input: 4_MetaG.MetaB.modules.NEU.linked.txt, 4_MetaB.HostT.modules.NEU.linked.txt, meta.mediation.NEU.txt or meta.mediation.EOS.txt, 1_metaG-combined.gct, 1_metaB.module_eigengene.txt, 1_hostT.module_eigengene.txt

Script: 6.random_forest.r

Output: 6_NEU_prediction.performance_byLinks.rf.txt (containing RMSE, RSQ and MAE scores for each linked MetaG-MetaB-HostT set)
```
