#!/bin/sh
# cd 1_processData/
# Rscript 1.1_findimportantgenes_papalexi.R ubergenno 4000
# Rscript 1.2_normalize_papalexi.R ubergenno 4000
# Rscript 1.3_getchromosome_papalexi.R ubergenno
# cd ../
cd 2_proximalSetup/
# Rscript 2.1_sparsePCA.R ubergenno 4000
# Rscript 2.2_choose_AYZW.R ubergenno A
cd ../
cd 3_proximalEstimation/
Rscript 3.X_estimate_SPCA_active.R ubergenno A
