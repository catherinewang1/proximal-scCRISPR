# ---------------------------------------------------------------------------- #
#                   desc                                                       #
# desc more                                                                    #
# ---------------------------------------------------------------------------- #
# ==================================================================================================
# =================== Sparse PCA (normalized genes) ================================================
# ==================================================================================================


# setup script parameters (e.g. paths) 
# DEVICE, NUM_IMPORTANT_GENES, K_FACTORS
args = commandArgs(trailingOnly = TRUE)
# args = c('laptop', 100, 20)



suppressPackageStartupMessages(library(crayon))
cat(crayon::blue('\n')) # for colors, sets some setting in terminal to display ANSI colors correctly

suppressPackageStartupMessages(require(assertthat)) # for some assert statements
suppressPackageStartupMessages(require(rhdf5))      # read/write HDF5 format
# suppressPackageStartupMessages(library(umap))       # perform umap
# suppressPackageStartupMessages(library(ggplot2))



assertthat::assert_that(length(args) > 0, msg="must give arg for specifying device eg 'Rscript <filename>.R ubergenno'")
DEVICE = args[1]
source('../PATHS.R') # load in data_dir and save_dir, depending on DEVICE value

assertthat::assert_that(length(args) >= 2, msg="provide device and NUM_IMPORTANT_GENES (currently int <=4000)")
NUM_IMPORTANT_GENES = as.integer(args[2])


assertthat::assert_that(length(args) > 0, msg="must give arg for specifying number of K_FACTORS for Sparse PCA")
K_FACTORS = as.integer(args[3])


myPrintColor <- function(txt, col=36) {
  # change number from 29:47 # for(col in 29:47){ cat(paste0("\033[0;", col, "m",txt,"\033[0m","\n"))}
  # cat(paste0("\033[0;", col, "m",txt,"\033[0m","\n"))
  cat(crayon::cyan(paste0(txt, '\n'))) 
}



# ==================================================================================================
# =================== Setup ========================================================================
# ==================================================================================================

# =================== START ====================================================
myPrintColor(sprintf("[%s] START: sparse PCA on papalexi (top %d genes)", Sys.time(), NUM_IMPORTANT_GENES))
h5file   = paste0(save_dir, "/gene_nonTF.h5")

# =================== loading gene_norm_nonTF ========================================
myPrintColor(sprintf("[%s]    - loading normalized gene exp", Sys.time()))
reading_hd5file  = rhdf5::H5Fopen(name = h5file)
readin_gene_norm = reading_hd5file&'gene_norm_nonTF'
gene_norm = readin_gene_norm[1:NUM_IMPORTANT_GENES, ] # dim = 4000 x 20729 = #important x #cells
rhdf5::h5closeAll()
invisible(gc(verbose=FALSE))


# =================== scaling gene_norm ========================================
myPrintColor(sprintf("[%s]    - scaling normalized gene exp", Sys.time()))
gene_norm_scaled = apply(gene_norm, MARGIN = 1, FUN = scale)
invisible(gc(verbose=FALSE))

# filling NAs in gene_norm_scaled
# anyNA(gene_norm_scaled) # but after feat sel, shouldn't have any?? some still 0sd
# gene_norm_scaled[is.na(gene_norm_scaled)] = 0 
# invisible(gc(verbose=FALSE))

# ==================================================================================================
# =================== performing Sparse PCA ================================================
# ==================================================================================================

# =================== performing sparse PCA ===========================================
myPrintColor(sprintf("[%s]    - performing sparse PCA (top ?? PCs)", Sys.time()))
t0 = Sys.time()
# is this okay? Choosing the penalty parameter based on a subsample, (for computational reasons)
cv.out <- PMA::SPC.cv(gene_norm_scaled[sample(1:nrow(gene_norm_scaled), 10000), ], nfolds = 3, vpos = TRUE) # choose penalty
myPrintColor(sprintf("[%s]       cv penalty cv.out$bestsumabsv1se = %.2f", Sys.time(), cv.out$bestsumabsv1se))
spca_fit = PMA::SPC(      x = gene_norm_scaled, 
                          sumabsv = cv.out$bestsumabsv1se,
                          # sumabsv = cv.out$bestsumabsv,
                          K = K_FACTORS,
                          orth = TRUE,
                          vpos = TRUE, # TRUE? we want `similar` genes
                          compute.pve = FALSE)
# pca_fit = gene_norm_scaled[,] |> prcomp(rank. = 50)
t1 = Sys.time(); print(sprintf('[%s]         %.2f mins', Sys.time(), difftime(t1, t0, units ="mins")))
invisible(gc(verbose=FALSE))

# `supergene` membership
smembership = apply(spca_fit$v, MARGIN = 1, FUN = which.max)
# gene_norm_scaled |> dim()
# colnames(gene_norm_scaled)
# 
# smembership
# smembership |> table()


# =================== saving sparse PCA result ========================================
myPrintColor(sprintf("[%s]    - saving sparse PCA", Sys.time()))
important_genes_nonTF = read.csv(sprintf('%s/supergene/important_genes_nonTF.csv', save_dir))
important_genes_nonTF$supergene_membership = NA
important_genes_nonTF$supergene_membership[1:length(smembership)] = smembership
write.csv(important_genes_nonTF,
          sprintf("%s/supergene/supergene_membership.csv", save_dir),
          row.names = FALSE)

# make a plot?? idk what kind of plot
# # =================== saving sparse PCA plot ==========================================
# print(sprintf("[%s]    - saving sparse PCA plot", Sys.time()))
# 
# p = ggplot(NULL,
#            aes(x = spca_fit$u[,1], y = spca_fit$u[,2])) +
#   geom_point(alpha=.05, size = .5) +
#   labs(x = 'SPC1', y = 'SPC2', title = NUM_IMPORTANT_GENES) +
#   theme_bw() +
#   theme(plot.title = element_text(hjust = .5))
# 
# ggsave(plot=p, 
#        filename=sprintf('%s/supergene/pca_%s.png', 
#                         save_dir, NUM_IMPORTANT_GENES),
#        height = 5, width = 6, dpi = 320)



# =================== construct supergene counts (sum members) from gene_nonTF =
myPrintColor(sprintf("[%s]    - construct supergene counts (sum members) from gene_nonTF", Sys.time()))
reading_hd5file  = rhdf5::H5Fopen(name = h5file)
readin_gene = reading_hd5file&'gene_nonTF'
gene = readin_gene[1:NUM_IMPORTANT_GENES, ] # dim = 4000 x 20729 = #important x #cells



# iteratively do this to save space 
# (could use dplyr, but idk how much space this would use by loading all into memory and converting to a df, maybe it's not so bad)


sgene = matrix(NA, nrow = max(smembership), ncol = ncol(gene_norm))
for(i in 1:max(smembership)) {
  sgeneiidx  = which(smembership == i)
  # test = readin_gene[idx, 1:20] ; dim(test)
  # test = readin_gene[which(smembership == i), 1:20]  ; dim(test)
  # test = readin_gene[c(1, 2), ]
  # dim(test)
  # test
  sgene[i, ] = readin_gene[sgeneiidx, ] |>
               colSums()
  
}
               
rhdf5::h5closeAll()
invisible(gc(verbose=FALSE))               
         
# =================== saving sgene expression =======================================
myPrintColor(sprintf("[%s]        - saving super gene expression", Sys.time()))

# delete existing gene_norm_nonTF if already exists
h5f         = rhdf5::H5Fopen(h5file); 
name_exists = rhdf5::H5Lexists(h5f, 'sgene_nonTF')
rhdf5::h5closeAll()
if (name_exists) {  rhdf5::h5delete(file = h5file, name = 'sgene_nonTF') }
rm(h5f, name_exists)

# save normalized gene expression of top genes in 'gene.h5' under 'gene_norm_nonTF'
if(DEVICE == 'ubergenno') { #if on ubergenno, need to realize to memory first?????
  HDF5Array::writeHDF5Array(DelayedArray::realize(sgene), filepath=h5file, name='sgene_nonTF')
} else {
  HDF5Array::writeHDF5Array(sgene, filepath=h5file, name='sgene_nonTF')
}

invisible(gc(verbose=FALSE))


      
# =================== END ======================================================
myPrintColor(sprintf("[%s] END", Sys.time()))









# # test the sparse pca function 
# install.packages('PMA')
# 
# library(PMA)
# 
# cv.out <- SPC.cv(gene_norm_scaled, nfolds = 3, vpos = TRUE)
# test = PMA::SPC(          x = gene_norm_scaled, 
#                     # sumabsv = cv.out$bestsumabsv1se,
#                     sumabsv = cv.out$bestsumabsv,
#                           K = 3,
#                        orth = TRUE,
#                        vpos = TRUE, # TRUE? we want `similarly` genes
#                 compute.pve = TRUE)
# test
# (test$u[, 1] * test$u[, 3]) |> sum()
# # hmmm genes can have nonzero coefs for each loading
# # each loading should be standardized to variance 1? So take the first loading with the maximum magnitude coef?
# test$v[, 1]
# test$v
# plot(test$prop.var.explained)
# 
# # `supergene` membership
# smembership = apply(test$v, MARGIN = 1, FUN = which.max)
# 
# table(smembership)
# 
# which.max(c(1, 2, 43, 2))






