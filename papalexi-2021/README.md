

# Analysis of `papalexi-2021`


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



## `2_otherEstimates`



## `3_proximalEstimate`










# Extra:


Some notes and details for running code:


 + `ondisc` package updated, and the code doesn't work with the new package version (no longer has the function `ondisc::read_odm`). Download the tar.gz file of v1.1.0 from https://github.com/timothy-barry/ondisc/releases and install. 
   (I did not modify the code to fit the new package because the `papalexi-2021` data is pulled from https://github.com/Katsevich-Lab/import-papalexi-2021, which uses the older `ondisc` version).
   I saved the tar.gz locally at genData/package_ondisc_old/


