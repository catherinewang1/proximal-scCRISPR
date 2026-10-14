
# concatenate data
# python 3.1.1_jin2025_concat_data.py


# find important genes, normalize, calculate PCA loadings for Neg Controls
# Rscript 3.1.2_jin2025_processData.R macbook

# some estimates for proximal
nice -10 Rscript 3.2.1_jin2025_estimate_glm_prox.R ubergenno small       pca_both 28
nice -10 Rscript 3.2.1_jin2025_estimate_glm_prox.R ubergenno denseSmall  pca_both 28
nice -10 Rscript 3.2.1_jin2025_estimate_glm_prox.R ubergenno denseMedium pca_both 28

