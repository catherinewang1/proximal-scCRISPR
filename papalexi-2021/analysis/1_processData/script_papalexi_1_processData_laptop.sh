#!/bin/sh
Rscript 1.1_papalexi_findimportantgenes.R macbook 4000
Rscript 1.2_papalexi_normalize.R macbook 4000
Rscript 1.3_papalexi_getchromosome.R macbook
Rscript 1.4_papalexi_sparsePCA.R macbook 4000
Rscript 1.5_papalexi_wgcna.R macbook 4000

# plots
# Rscript 1.5_dimred_plots_papalexi.R laptop 4000
# Rscript 1.4_pca_umap_papalexi.R laptop 2000
# Rscript 1.5_dimred_plots_papalexi.R laptop 2000
# Rscript 1.4_pca_umap_papalexi.R laptop 1000
# Rscript 1.5_dimred_plots_papalexi.R laptop 1000