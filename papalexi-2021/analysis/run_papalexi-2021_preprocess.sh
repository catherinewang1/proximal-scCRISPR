#!/bin/sh
cd 1_processData/
Rscript 1.1_findimportantgenes_papalexi.R laptop 4000
Rscript 1.2_normalize_papalexi.R laptop 4000
Rscript 1.3_getchromosome_papalexi.R laptop
cd ../
cd 2_proximalSetup/
Rscript 2.1_sparsePCA.R laptop 4000
Rscript 2.2_choose_AYZW.R laptop A
cd ../
cd 3_proximalEstimation/
Rscript 3.X_estimate_SPCA_active.R laptop A
