



Analyze scCRISPR screen experiments using Proximal Causal Inference




# Datasets
+ papalexi-2025: more controlled experimental setting
+ (mouse in vivo dataset): more varying/noisy setting


# Analysis Overview

Steps

+ Pre-process data
    - filter genes and perturbations (high variance/signal, non TF (TODO: where did the list of TF come from? slack attachment from Kathryn...))
    - convert from counts to continuous! 
    - choose some perturbation-gene pairs for testing 
+ Other Estimation Methods
    - estimate ATE/effects using other methods (some may be p-values and not ATE, e.g. nb est using ra counts)
    - linear regression, Negative Binomial, SCEPTRE [^barry2021], 
+ Proximal Setup and Estimation
    - construct NCE/NCOs (sparse PCA [^witten2007])
    - estimate ATE using pci2s
+ Evaluation of Estimators  
    - plots shoing performance



# Repository Structure

```
.
`-- Proximal/
    |-- README.md
    |-- papalexi-2021/
    |   |-- analysis/
    |   |   |-- 1_processData
    |   |   |-- 2_proximalSetup
    |   |   |-- 3_proximalEstimation
    |   |   `-- 4_proximalEvaluation
    |   `-- saves/
    |       `-- ...
    |-- mouse-xxxx/
    |   |   |-- 1_processData
    |   |   |-- 2_proximalSetup
    |   |   |-- 3_proximalEstimation
    |   |   `-- 4_proximalEvaluation
    |   `-- saves/
    |       `-- ...
    `-- utils/
        |-- code1.R
        |-- code2.R
        `-- ...
```





# References


[^barry2021]: Barry, Timothy et al. (2021). “SCEPTRE improves calibration and sensitivity in single-cell CRISPR screen analysis”. *Genome Biology* 22.1, p. 344. issn: 1474-760X.


[^witten2007]: Witten, Daniela M., Robert Tibshirani, and Trevor Hastie (July 2009). “A penalized matrix decomposition, with applications to sparse principal components and canonical correlation analysis”. *Biostatistics (Oxford, England)* 10.3, pp. 15–534. issn: 1468-4357.




TODO: add papalexi dataset, and mouse dataset, pcis paper/package, proximal inf papers (a collection).

