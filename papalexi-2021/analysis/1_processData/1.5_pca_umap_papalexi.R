# ---------------------------------------------------------------------------- #
#                   Perform PCA then PCA/UMAP on papalexi dataset              #                            #
# Requires: prev saved normalized gene expression (HDF5)                       #
# Ouputs: (nothing) but saves dim red and plots at                             #
#         '<save_dir>/pca/'                                                    #
#         '<save_dir>/pca_umap/'                                               #
# ---------------------------------------------------------------------------- #
args = commandArgs(trailingOnly = TRUE)
    
require(assertthat) # for some assert statements
require(rhdf5)      # read/write HDF5 format
library(umap)       # perform umap
library(ggplot2)

assertthat::assert_that(length(args) > 0, msg="must give arg for specifying device eg 'Rscript <filename>.R ubergenno'")
DEVICE = args[1]
source('../PATHS.R') # load in data_dir and save_dir, depending on DEVICE value

assertthat::assert_that(length(args) >= 2, msg="provide device and NUM_IMPORTANT_GENES (currently int <=4000)")
NUM_IMPORTANT_GENES = as.integer(args[2])


# =================== START ====================================================
print(sprintf("[%s] START: PCA and PCA/UMAP papalexi (top %d genes)", Sys.time(), NUM_IMPORTANT_GENES))
h5file      = paste0(save_dir, "/gene.h5")


# =================== loading gene_norm ========================================
print(sprintf("[%s]    - loading normalized gene exp", Sys.time()))
reading_hd5file  = rhdf5::H5Fopen(name = h5file)
readin_gene_norm = reading_hd5file&'gene_norm'
gene_norm = readin_gene_norm[1:NUM_IMPORTANT_GENES, ] # dim = 4000 x 20729 = #important x #cells
rhdf5::h5closeAll()
invisible(gc(verbose=FALSE))


# =================== scaling gene_norm ========================================
print(sprintf("[%s]    - scaling normalized gene exp", Sys.time()))
gene_norm_scaled = apply(gene_norm, MARGIN = 1, FUN = scale)
invisible(gc(verbose=FALSE))


# =================== performing PCA ===========================================
print(sprintf("[%s]    - performing PCA (to 50 PCs)", Sys.time()))
t0 = Sys.time()
pca_fit = gene_norm_scaled[,] |> prcomp(rank. = 50)
t1 = Sys.time(); print(sprintf('[%s]         %.2f mins', Sys.time(), difftime(t1, t0, units ="mins")))
invisible(gc(verbose=FALSE))


# =================== saving PCA result ========================================
print(sprintf("[%s]    - saving PCA", Sys.time()))
dir.create(sprintf("%s/pca/objects/", save_dir), recursive=TRUE, showWarnings=FALSE)
saveRDS(pca_fit,
        sprintf("%s/pca/objects/pca_%s.rds", save_dir, NUM_IMPORTANT_GENES))


# =================== saving PCA plot ==========================================
print(sprintf("[%s]    - saving PCA plot", Sys.time()))

p = ggplot(NULL,
           aes(x = pca_fit$x[,'PC1'], y = pca_fit$x[,'PC2'])) +
  geom_point(alpha=.05, size = .5) +
  labs(x = 'PC1', y = 'PC2', title = NUM_IMPORTANT_GENES) +
  theme_bw() +
  theme(plot.title = element_text(hjust = .5))

ggsave(plot=p, 
       filename=sprintf('%s/pca/pca_%s.png', 
                        save_dir, NUM_IMPORTANT_GENES),
       height = 5, width = 6, dpi = 320)


# =================== performing UMAP ==========================================
print(sprintf("[%s]    - performing UMAP", Sys.time()))
t0 = Sys.time()
umap_fit = pca_fit$x |> umap(n_components = 2)
t1 = Sys.time(); print(sprintf('[%s]         %.2f mins', Sys.time(), difftime(t1, t0, units ="mins")))
invisible(gc(verbose=FALSE))


# =================== saving UMAP result =======================================
print(sprintf("[%s]    - saving UMAP", Sys.time()))
umap_fit$knn  = 'removed for space'

dir.create(sprintf("%s/pca_umap/objects/", save_dir), recursive=TRUE, showWarnings=FALSE)
saveRDS(umap_fit,
        sprintf("%s/pca_umap/objects/pca_umap_%s.rds", save_dir, NUM_IMPORTANT_GENES))


# =================== saving UMAP plot =========================================
print(sprintf("[%s]    - saving UMAP plot", Sys.time()))

p = ggplot(NULL,
           aes(x = umap_fit$layout[,1], y = umap_fit$layout[,2])) +
  geom_point(alpha=.05, size = .5) +
  labs(x = 'pca-umap1', y = 'pca-umap2', title = NUM_IMPORTANT_GENES) +
  theme_bw() +
  theme(plot.title = element_text(hjust = .5))

ggsave(plot=p, 
       filename=sprintf('%s/pca_umap/pca_umap_%s.png', 
                        save_dir, NUM_IMPORTANT_GENES),
       height = 5, width = 6, dpi = 320)


# =================== END ======================================================
print(sprintf("[%s] END", Sys.time()))
