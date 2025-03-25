#!/bin/sh
Rscript 1.1_findimportantgenes_papalexi.R laptop 4000
Rscript 1.2_normalize_papalexi.R laptop 4000
Rscript 1.3_getchromosome_papalexi.R 
# Rscript 1.4_pca_umap_papalexi.R laptop 4000
Rscript 1.5_dimred_plots_papalexi.R laptop 4000
# Rscript 1.4_pca_umap_papalexi.R laptop 2000
# Rscript 1.5_dimred_plots_papalexi.R laptop 2000
# Rscript 1.4_pca_umap_papalexi.R laptop 1000
# Rscript 1.5_dimred_plots_papalexi.R laptop 1000