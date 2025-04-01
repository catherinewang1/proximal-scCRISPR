# proximal-scCRISPR
Proximal causal analysis on the scCRISPR papalexi-2021 dataset [(Papalexi et al., 2021)](https://www.nature.com/articles/s41588-021-00778-2), which is kindly processed/cleaned by Tim Barry at https://github.com/Katsevich-Lab/import-papalexi-2021.



## Testing for perturbation-gene causal effects
To test for causal effects beteen a selections of perturbation (A) and gene (Y) pairs, run the .sh script `papalexi-2021/analysis/run_papalexi-2021_analysis.sh`. The tests currently included are the naive linear regression of Y on A method (lmYA) and 2SLS proximal inference method implemented by Liu et al., 2024 (pci2s).



## Active p-values 


To prepare the data and to setup the AY pairs for the active p-value procedure, run the script 
`papalexi-2021/analysis/run_papalexi-2021_preprocess.sh`.
(This is basically the same script as the previous script that runs the analysis, but it excludes the calculation of the p-values. TODO: Instead, the p-values are calculated in the active method.)





