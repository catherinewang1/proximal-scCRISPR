

# Analysis of `papalexi-2021` using proximal causal inference


# Folders



## `1_processData`


Script to run all: `script_papalexi_1_processData_laptop.sh`

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


## `2_setupAY`


Setup a selection of perturbations (A) and genes (Y) to test the causal effects for.
Specify the settings for choosing tests in `<save_dir>/AY/AY_setting.R`.

Files in the folder `2_setupAY`

 + `2.1_choose_AY.R`

saved files include:

 + `<save_dir>/AY/<setting_name>/`
 + `<save_dir>/AY/<setting_name>/AY.csv`
 + `<save_dir>/AY/<setting_name>/AYZW_setting.rds`
 + `<save_dir>/AY/<setting_name>/AYZW_setting.txt`   



## `3_estimateEffects`

Estimate effects from a variety of methods.  


Files in the folder `3_estimateEffects`

 + `3.1_papalexi_countinuous.R`
 + `3.2_papalexi_countGLM.R`
 + `3.3_papalexi_sceptre.R`

saved files include:

 + `<save_dir>/AY/<setting_name>/<proximal_setting_name>/`
 + `<save_dir>/AY/<setting_name>/<proximal_setting_name>/effects_continuous.csv`
 + `<save_dir>/AY/<setting_name>/<proximal_setting_name>/proximal_setting.rds`
 + `<save_dir>/AY/<setting_name>/<proximal_setting_name>/intermediateATEs/`

## `4_compareEstimators`










# Extra:


Some notes and details for running code:


 + `ondisc` package updated, and the code doesn't work with the new package version (no longer has the function `ondisc::read_odm`). Download the tar.gz file of v1.1.0 from https://github.com/timothy-barry/ondisc/releases and install. 
   (I did not modify the code to fit the new package because the `papalexi-2021` data is pulled from https://github.com/Katsevich-Lab/import-papalexi-2021, which uses the older `ondisc` version).
   I saved the tar.gz locally at genData/package_ondisc_old/


