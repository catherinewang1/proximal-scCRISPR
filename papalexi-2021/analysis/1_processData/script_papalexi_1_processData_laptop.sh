#!/bin/sh
Rscript 1.1_papalexi_normalize.R laptop 4000
Rscript 1.2_papalexi_normalize.R laptop 4000
Rscript 1.3_papalexi_getchromosome.R laptop
Rscript 1.4_pca_umap_papalexi.R laptop 4000

# plots
# Rscript 1.5_dimred_plots_papalexi.R laptop 4000
# Rscript 1.4_pca_umap_papalexi.R laptop 2000
# Rscript 1.5_dimred_plots_papalexi.R laptop 2000
# Rscript 1.4_pca_umap_papalexi.R laptop 1000
# Rscript 1.5_dimred_plots_papalexi.R laptop 1000