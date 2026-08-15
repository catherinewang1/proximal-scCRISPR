
(TODO: need to update this for changed analysis procedure. there is a new readme inside analysis/. either put the analysis in this new file or copy new file contents here and delete the new file.)


# Analysis of `papalexi-2021` using proximal causal inference


# Folders

```
papalexi-2021/
  analysis/
    1_processData/
    2_setupAY/
    3_estimateEffects/
    4_compareEstimators/
    5_compareNCs/
  saves/
    AY/
  README.md
```

# Analysis Folder



## `1_processData`


Script to run: `script_papalexi_1_processData_laptop.sh`


**Part 1**: find important genes genes, normalize count data to continuous, gather chromosome info 

Files in the folder `1_processData/`

 + `1.1_papalexi_findimportantgenes.R`
 + `1.2_papalexi_normalize.R`
 + `1.3_papalexi_getchromosome.R`




saved files include:

 + `<save_dir>/gene_deviance.csv`
 + `<save_dir>/gene_deviance_topnoTFonly.csv`
 + `<save_dir>/gene.h5` (names `gene` and `gene_norm`) normalized gene expression 
 + `<save_dir>/chromosome/gene_chromosome.csv`
 + `<save_dir>/chromosome/grna_chromosome.csv`


**Part 2**: Construct potential Negative Control variables (NCE (Z) and NCO (W))) using PCA, sparse PCA, and WGCNA.


Files in the folder `1_processData/`

 + `1.4_papalexi_sparsePCA.R` - gathers both PCA and sparse PCA
 + `1.5_papalexi_wgcna.R` 

saved files include:

 + `<save_dir>/pca/NCloadings.rds`
 + `<save_dir>/spca/NCloadings_sumabs=34.5_K=60_N=5000.rds`
 + `<save_dir>/spca/NCloadings_sumabs=8.0_K=60_N=5000.rds`
 + `<save_dir>/WGCNA/NC_wgcna_ModuleEigengene.rds`
 + `<save_dir>/WGCNA/autoWGCNA/`- other WGCNA results (e.g. visualizations and different parameters)


Also includes `1.999_EDA_UMAP.R`  

## `2_setupAY`

Script to run: `2_setupAY_script.sh`

Setup a selection of perturbations (A) and genes (Y) to test the causal effects for.
And setup a selection of NCE (Z) and NCO (W) individual genes to use as negative controls.
Specify the settings for choosing tests in `<save_dir>/AY/AY_setting.R`.

Files in the folder `2_setupAY/`

 + `2.1_choose_AY.R`

saved files include:

 + `<save_dir>/AY/<setting_name>/`
 + `<save_dir>/AY/<setting_name>/AY.csv`
 + `<save_dir>/AY/<setting_name>/AYZW_setting.rds`
 + `<save_dir>/AY/<setting_name>/AYZW_setting.txt`   



## `3_estimateEffects`

Script to run: `2_setupAY_script.sh`


Estimate effects from a variety of methods.  


Files in the folder `3_estimateEffects/`

 + `3.1_papalexi_countinuous.R`
 + `3.2_papalexi_countGLM.R`
 + `3.3_papalexi_sceptre.R`

saved files include:

For continuous results:

 + `<save_dir>/AY/<setting_name>/PCA-SPCA8.0-SPCA34.5-WGCNA-singlegene/` - run many NC choices at once
 + `<save_dir>/AY/<setting_name>/PCA-SPCA8.0-SPCA34.5-WGCNA-singlegene/effects_continuous.csv` 
 + `<save_dir>/AY/<setting_name>/PCA-SPCA8.0-SPCA34.5-WGCNA-singlegene/intermediateATEs/` 
 + `<save_dir>/AY/<setting_name>/<proximal_setting_name>/`- also save individually per NC choice
 + `<save_dir>/AY/<setting_name>/<proximal_setting_name>/effects_continuous.csv`
 + `<save_dir>/AY/<setting_name>/<proximal_setting_name>/proximal_setting.rds`
 + `<save_dir>/AY/<setting_name>/<proximal_setting_name>/proximal_setting.txt`
 + `<save_dir>/AY/<setting_name>/<proximal_setting_name>/intermediateATEs/`

For count results:

 + `<save_dir>/AY/<setting_name>/simpleCount/`
 + `<save_dir>/AY/<setting_name>/simpleCount/effects_countGLM.csv`
 + `<save_dir>/AY/<setting_name>/simpleCount/intermediateATEsCountGLM/`



## `4_compareEstimators`



## `5_compareNCs`

Files in the folder `5_compareNCs/`

 + `5.1_visualizeProxByNC.qmd`




# Extra:


Some notes and details for running code:


 + `ondisc` package updated, and the code doesn't work with the new package version (no longer has the function `ondisc::read_odm`). Download the tar.gz file of v1.1.0 from https://github.com/timothy-barry/ondisc/releases and install. 
   (I did not modify the code to fit the new package because the `papalexi-2021` data is pulled from https://github.com/Katsevich-Lab/import-papalexi-2021, which uses the older `ondisc` version).
   I saved the tar.gz locally at genData/package_ondisc_old/


