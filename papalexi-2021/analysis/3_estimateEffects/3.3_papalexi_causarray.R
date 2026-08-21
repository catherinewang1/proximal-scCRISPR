# ------------------------------------------------------------------------------------------------ #
#     Estimate Effect (ATE) using causarray 
# e.g.     nice -10 Rscript 3.3_papalexi_causarray.R ubergenno A
# ------------------------------------------------------------------------------------------------ #
args = commandArgs(trailingOnly = TRUE)
# args = c('laptop',  'A1')
# args = c('laptop',  'A')
# args = c('macbook', 'A')
# args = c('macbook', 'A1') 




# =================== Load Libraries and Packages =========================================================
suppressPackageStartupMessages(library(assertthat)) # for some assert statements
# suppressPackageStartupMessages(library(dplyr))



# suppressPackageStartupMessages(library(future.apply))
# # options(future.globals.maxSize= 850*1024^2) #1st num is MB
# options(future.globals.maxSize= 2500*1024^2) #1st num is MB
# plan(multisession, workers = 24)
# # plan(sequential)



# =================== Read input args =========================================================
# args[1] = device e.g. macbook or laptop or ubergenno
assertthat::assert_that(length(args) > 0, msg="must give arg for specifying device eg 'Rscript <filename>.R ubergenno'")
DEVICE = args[1]
source('../PATHS.R') # load in data_dir and save_dir and CODE_DIR, depending on DEVICE value
assertthat::assert_that(!is.null(data_dir), msg='first arg must be: laptop, desktop, or ubergenno')
# 
# # args[2] = setting name e.g. A1 or A or B
# assertthat::assert_that(length(args) > 1, msg="must give arg for specifying chosen AYZW name 'Rscript <filename>.R ubergenno C'")
# AYZW_setting_name = args[2]
# 
# 
# 
# 
# 
# # =================== Start ========================================================================
# print(sprintf("[%s] START: Estimate Papalexi-2021 as counts", Sys.time()))
# print(sprintf("[%s]        w/ script: %s", Sys.time(), paste0(args, collapse = ', ')))
# 
# 
# print(sprintf("[%s]    - Load files and functions", Sys.time()))
# 
# # source(sprintf('%s/estimate_effects_pci2s_count.R', util_dir)) # for functions to estimate here
# 
# 
# # source(sprintf('%s/CBEstAll.R', util_dir)) # for functions to estimate here
# # source(sprintf('%s/CBEstAllSPCA.R', util_dir)) # add'l specifically for SPCA
# 
# 
# # load chosen AYZW names
# AY   = read.csv(sprintf('%s/AY/%s/AY.csv', save_dir, AYZW_setting_name))
# # AYZW = readRDS(sprintf('%s/cbgenes/%s/AYZW.rds', save_dir, AYZW_setting_name))
# 
# # # Only for PCA: copy to pca folder too (bc plots script reads in AY, AYZW)
# # write.csv(AY, sprintf('%s/spca/cbgenes/%s/AY.csv', save_dir, AYZW_setting_name))
# # saveRDS(AYZW, sprintf('%s/spca/cbgenes/%s/AYZW.rds', save_dir, AYZW_setting_name))
# 
# # load gene importance info: gene_name, gene_idx (idx for raw), importance_rank (rank + idx for normalized)
# # gene_importance = read.csv(sprintf('%s/gene_deviance_topnoTFonly.csv', save_dir)) |> 
# #                   dplyr::select(gene_name, gene_idx, importance_rank)
# 
# gene_importance = read.csv(sprintf('%s/gene_deviance_gene_norm.csv', save_dir)) |> 
#   dplyr::select(gene_name, gene_idx, importance_rank, gene_norm_idx)
# 
# # Load choosing AY (ZW) settings
# # AYZW_setting = readRDS(sprintf('%s/AY/%s/AYZW_setting.rds', save_dir, AYZW_setting_name))
# 
# 
# # create gene and grna ondisc managers
# gene_odm <- ondisc::read_odm(odm_fp      = paste0(data_dir, "/papalexi-2021/processed/gene/expression_matrix.odm"),
#                              metadata_fp = paste0(data_dir, "/papalexi-2021/processed/gene/metadata.rds"))
# grna_odm <- ondisc::read_odm(odm_fp      = paste0(data_dir, "/papalexi-2021/processed/grna_assignment/assignment_matrix.odm"),
#                              metadata_fp = paste0(data_dir, "/papalexi-2021/processed/grna_assignment/metadata.rds"))
# 
# # load grna assignments (load all into memory)
# # grna = grna_odm[[,1:ncol(grna_odm)]] # |> as.matrix() # ~110 x 20729 = #grnas x #cells
# # grna_rownames = grna_odm |> ondisc::get_feature_covariates() |> rownames()
# # row.names(grna) = grna_rownames
# 
# # load measured covariates (=U, pretend unmeasured for proximal estimation)
# cell_covariates = gene_odm |> ondisc::get_cell_covariates()
# # change batch info, bc lane determines rep_1, so redundant info
# # table(extra_covariates |> dplyr::select(lane, bio_rep))
# library_size = cell_covariates$n_umis
# cell_covariates = cell_covariates |> 
#   dplyr::mutate(lane_bio_rep = paste0(lane, '_', bio_rep)) |>
#   dplyr::select(-lane, -bio_rep) |>
#   dplyr::select(-n_nonzero, -n_umis)
# # Does not include offset! library size for cell i = \sum_gene count_{i, gene}
# # > head(cell_covariates, 1+1)
# #                         phase     p_mito lane_bio_rep
# #     l1_AAACCTGAGCCAGAAC    G1 0.02295577  Lane1_rep_1
# #     l1_AAACCTGAGTGGACGT    G1 0.04512939  Lane1_rep_1
# 
# 


# =================== causarray tutorial by Jinhong ========================================================================
library(Seurat)


# data_dir = "/Users/catherinewang/Documents/School/genData/papalexi"

# jincounts = read.delim(sprintf('%s/../jin/other/Counts.PBC.txt', data_dir))
# jincounts$gene |> unique() |> length()
# summary(jincounts)

sc.seurat <- readRDS(sprintf('%s/../jin/causarrayexampledata/perturbseq-exneu.rds', data_dir))

# Access the counts through Seurat when its installed version recognizes the
# serialized Assay5 object; otherwise recover the same layer and dimnames
# directly. The fallback keeps this tutorial runnable with older Seurat builds.
if ("RNA" %in% Assays(sc.seurat)) {
  counts <- GetAssayData(sc.seurat, assay = "RNA", layer = "counts")
} else {
  rna.assay <- sc.seurat@assays[["RNA"]]
  counts <- rna.assay@layers[["counts"]]
  rownames(counts) <- rownames(rna.assay@features)
  colnames(counts) <- rownames(rna.assay@cells)
}
Y <- data.frame(t(as.matrix(counts)), check.names = FALSE) # cell-by-gene matrix
metadata <- sc.seurat@meta.data

perturb <- metadata
colnames(perturb) <- gsub("Perturbation", "trt_", colnames(perturb))
perturb$trt_ <- relevel(as.factor(perturb$trt_), ref = "GFP")
A <- data.frame(
  model.matrix(~ trt_ - 1, data = perturb)[, -1, drop = FALSE],
  check.names = FALSE
) # cell-by-trt matrix; remove the first (GFP control) column
colnames(A) <- sub("^trt_", "", colnames(A))


require(reticulate)

Sys.setenv(PYTHONUNBUFFERED = TRUE)
use_condaenv('causarray')
causarray <- import("causarray")
cat(causarray$`__version__`)

# (Y, A) should be either data.frame or matrix
# optional covariates can be provided as matrices
dat <- causarray$prep_causarray_data(Y, A)
names(dat) <- c("Y", "A", "X", "X_A")
list2env(dat, .GlobalEnv)


r <- 10
# Same early-stopping settings as the Python tutorial, so both produce
# identical latent factors on this dataset.
res_gate <- causarray$fit_gcate(
  Y, X, A, r, verbose = TRUE,
  kwargs_es_1 = list(rel_tol = 2e-4, max_iters = 30L),
  kwargs_es_2 = list(rel_tol = 2e-4, max_iters = 30L)
) # a list of results from 2 stages optimization

U <- res_gate[[2]]$U
# reticulate returns `hist` as an R list, so unlist() before min().
cat(sprintf("Step 1 -- epochs: %d, best NLL: %.6f\n",
            as.integer(res_gate[[1]]$n_iter), min(unlist(res_gate[[1]]$hist))))


cat(sprintf("Step 2 -- epochs: %d, best NLL: %.6f\n",
            as.integer(res_gate[[2]]$n_iter), min(unlist(res_gate[[2]]$hist))))


offsets <- log(res_gate[[2]][['kwargs_glm']][['size_factor']]) # use the precomputed size factors
res <- causarray$LFC(Y, cbind(X, U), A, cbind(X_A, U), offset=offsets,
                     usevar="pooled", verbose=TRUE)

names(res) <- c("df_res", "estimation")
list2env(res, .GlobalEnv)




# Diagnose treatment associations and overlap --------






