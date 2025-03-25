# ---------------------------------------------------------------------------- #
#                   Do CB using genes as NCE/NCO                               #
# 3.1 - add chromosome number information                                      #
# 3.2 - choose A,Y,Z,W combinations                                            #
# 3.2 - CB Effect Estimate                                                     #
# Requires: prev saved normalized gene expression (HDF5)                       #
# ---------------------------------------------------------------------------- #
args = commandArgs(trailingOnly = TRUE)
args = c('laptop', '2000', 'Atest')

require(assertthat) # for some assert statements
library(Matrix)
library(dplyr)
library(tibble)
library(ggplot2)    # plotting
library(cowplot)


library(future.apply)
options(future.globals.maxSize= 850*1024^2) #1st num is MB
plan(multisession, workers = 8)
# plan(sequential)

# library(furrr)
# plan(multisession, workers = 2)

theme_set(theme_cowplot() +
            theme(plot.title = element_text(hjust = .5),
                  plot.subtitle = element_text(hjust = .5)))

# === Parameter Settings for CB Methods ===
ALPHA = .1
GMM_steps = 10
K_folds = 10
lambdas = (10)^(seq(from = -10, to = 3, length.out = 27))
TRIM = .1
# num_NC_pairs = c(1:8, 10, 15, 20, 25, 30) # number of NC pairs to try (<= total NC pairs, will be cut off)
num_NC_pairs = c(1:3, 5, 10, 25, 50)
# num_NC_pairs = c(1:3)
basisParams = list(basis1 = list(b_degree=1, h_degree=1, # linear
                                 b_phases=NULL, b_periods=NULL,
                                 h_phases=NULL, h_periods=NULL)# ,
                   # basis2 = list(b_degree=4, h_degree=2, # larger b deg
                   #               b_phases=NULL, b_periods=NULL,
                   #               h_phases=NULL, h_periods=NULL)
) #,
# basis3 = list(b_degree=4, h_degree=2, # add reasonable cos
#               b_phases=c(0), b_periods=c(4, 8, 16),
#               h_phases=c(0), h_periods=c(4, 8, 16)))
CB_setting_name = 'simple'
run_OCBLinOSEst = FALSE # whether to run OCBLinOne-Step estimator or not
save_intermediateATEs = 'yes' # 'yes'/'no' whether to save intermedate ATEs as they are estimated
# === === === === === === === === === ===

assertthat::assert_that(length(args) > 0, msg="must give arg for specifying device eg 'Rscript <filename>.R ubergenno'")



DEVICE = args[1]
source('../PATHS.R') # load in data_dir and save_dir and CODE_DIR, depending on DEVICE value

# DEVICE = args[1]
# # location of papalexi-2021 folder
# data_dir = switch(DEVICE,
#                   'laptop'='C:/Users/Cathe/Documents/School/Genetic Hypothesis Testing using Negative Controls/genData/papalexi',
#                   'desktop'='C:/Users/Catherine W/Documents/Research/genData/papalexi',
#                   'ubergenno'='/raid6/Catherine/papalexi')
# # location of intermediate save files/plots are written
# save_dir = switch(DEVICE,
#                   'laptop'='C:/Users/Cathe/Documents/School/Genetic Hypothesis Testing using Negative Controls/DoubleBridge/saves/papalexi_saves',
#                   'desktop'='C:/Users/Catherine W/Documents/Research/DoubleBridge/saves/papalexi_saves',
#                   'ubergenno'='/raid6/Catherine/papalexi/papalexi_saves')
# # location of utils code
# util_dir = switch(DEVICE,
#                   'laptop'='C:/Users/Cathe/Documents/School/Genetic Hypothesis Testing using Negative Controls/DoubleBridge/code/utils',
#                   'desktop'='C:/Users/Catherine W/Documents/Research/DoubleBridge/code/utils',
#                   'ubergenno'='/home/catheri2/DoubleBridge/code/utils')

assertthat::assert_that(!is.null(data_dir), msg='first arg must be: laptop, desktop, or ubergenno')

assertthat::assert_that(length(args) > 1, msg="must give arg for specifying num imp genes 'Rscript <filename>.R ubergenno 4000'")
NUM_IMPORTANTGENES = as.integer(args[2]) # should be max importance from setting... TODO: change to remove this input. extract this value from saved setting

assertthat::assert_that(length(args) > 2, msg="must give arg for specifying chosen AYZW name 'Rscript <filename>.R ubergenno C'")
AYZW_setting_name = args[3]


# save parameter setting for basis + alpha
CB_setting = list() # new copy to add ALPHA to
CB_setting$ALPHA     = ALPHA # R does not save as pointers, so ok
CB_setting$GMM_steps = GMM_steps
CB_setting$K_folds   = K_folds
CB_setting$lambdas   = lambdas
CB_setting$TRIM      = TRIM
CB_setting$basisParams  = basisParams
CB_setting$num_NC_pairs = num_NC_pairs
CB_setting$run_OCBLinOSEst = run_OCBLinOSEst
dir.create(sprintf('%s/supergene/%s/%s', save_dir, AYZW_setting_name, CB_setting_name), recursive = TRUE)
capture.output(print(CB_setting),
               file = sprintf('%s/supergene/%s/%s/CB_setting.txt',
                              save_dir, AYZW_setting_name, CB_setting_name))
saveRDS(CB_setting,
        sprintf('%s/supergene/%s/%s/CB_setting.rds',
                save_dir, AYZW_setting_name, CB_setting_name))

# =================== Start ====================================================
print(sprintf("[%s] START: CB Estimate", Sys.time()))

# # source CB utility functions (mainly CB fn for estimating)
# source(sprintf('%s/CB_utils.R', util_dir))
# Instead source new estimation methods' files
source(sprintf('%s/OCB2SLS.R', util_dir))      # 2-Stage Least Squares
source(sprintf('%s/OCBGMM.R', util_dir))       # Gen. Method of Moments
source(sprintf('%s/OCBLinBridge.R', util_dir)) # Lin Bridge Plug-In and One-Step

source(sprintf('%s/getdfPCA.R', util_dir))     # for performing PCA

source(sprintf('%s/CBEstAll.R', util_dir)) # for functions to estimate here


# load chosen AYZW names
AY   = read.csv(sprintf('%s/supergene/%s/AY.csv', save_dir, AYZW_setting_name))
AYZW = readRDS(sprintf('%s/supergene/%s/AYZW.rds', save_dir, AYZW_setting_name))

# Only for PCA: copy to pca folder too (bc plots script reads in AY, AYZW)
write.csv(AY, sprintf('%s/supergene/%s/AY.csv', save_dir, AYZW_setting_name))
saveRDS(AYZW, sprintf('%s/supergene/%s/AYZW.rds', save_dir, AYZW_setting_name))

# load gene importance info
imp_gene_names = readRDS(sprintf('%s/important_genes_name.rds', save_dir))
imp_gene_idx   = readRDS(sprintf('%s/important_genes_idx.rds',  save_dir))
imp_gene = data.frame(gene     = imp_gene_names,
                      gene_idx = imp_gene_idx,
                      gene_imp_rank = 1:length(imp_gene_names))
# load in prev saved df (safer than constructing df again?)
# imp_gene = read.csv(sprintf('%s/gene_deviance.csv', save_dir)) |>
#               arrange(desc(deviance)) |>
#               rename(gene=gene_name, gene_idx=idx) |>
#               select(gene, gene_idx, deviance) |>
#               mutate(gene_imp_rank = row_number())


# create gene and grna ondisc managers
gene_odm <- ondisc::read_odm(odm_fp      = paste0(data_dir, "/papalexi-2021/processed/gene/expression_matrix.odm"),
                             metadata_fp = paste0(data_dir, "/papalexi-2021/processed/gene/metadata.rds"))
grna_odm <- ondisc::read_odm(odm_fp      = paste0(data_dir, "/papalexi-2021/processed/grna_assignment/assignment_matrix.odm"),
                             metadata_fp = paste0(data_dir, "/papalexi-2021/processed/grna_assignment/metadata.rds"))

# load grna assignments (load all into memory)
grna = grna_odm[[,1:ncol(grna_odm)]] # |> as.matrix() # ~110 x 20729 = #grnas x #cells
grna_rownames = grna_odm |> ondisc::get_feature_covariates() |> rownames()

# load normalized gene exp
h5file      = paste0(save_dir, "/gene.h5"); print(h5file)
reading_hd5file  = rhdf5::H5Fopen(name = h5file)
readin_gene_norm = reading_hd5file&'gene_norm'
gene_norm = readin_gene_norm[1:NUM_IMPORTANTGENES, ] # dim = 4000 x 20729 = #important x #cells
rhdf5::h5closeAll()
invisible(gc(verbose=FALSE))



# =================== Define File Specific Functions ======================================
print(sprintf("[%s]    - Define File Specific Functions", Sys.time()))




# All the NT grna idx
NT_names = grna_odm |> ondisc::get_feature_covariates() |>
  filter(target_type == 'non-targeting') |> rownames()
NT_idx = which(apply(X = grna_odm[[NT_names, ]], MARGIN = 2, FUN = sum) > 0)







# ======================================
# === === === === === === === === === === === === === === === === === 
# get_importance_rank = get_importance_rank_make(imp_gene_names)
# AY_idx = 1; ZW_idx = 1; pooled=FALSE
# AY_row = AY[AY_idx, ]
# 
# 
# # gather the AYZW names chosen
# AYZW_cur = AYZW[[AY_row$A]][[AY_row$Y]][[ZW_idx]]
# 
# 
# # AY_row$A
# # AY_row$Y
# # AYZW_cur$Z_names
# # AYZW_cur$W_names
# 
# # Get A
# # -------------------------------------------
# # # idx of all 'treated' cells
# # A_idx = which(as.logical(grna_odm[[AY_row$A, ]]))
# # # subset cells of 'treated' (w A grna) and 'control' (NT gran)
# # A = grna_odm[[AY_row$A, c(A_idx, NT_idx)]]
# 
# # using grna all loaded into memory
# # idx of all 'treated' cells
# 
# A_grna_idx = which(grna_rownames == AY_row$A)
# A_idx = which(as.logical(grna[A_grna_idx, 1:ncol(grna)]))
# # subset cells of 'treated' (w A grna) and 'control' (NT gran)
# A = grna[A_grna_idx, c(A_idx, NT_idx)]
# 
# 
# # # Get Y and Z and W
# # # using h5 load in (use if not enough memory)
# # # -------------------------------------------
# # h5file      = paste0(save_dir, "/gene.h5")
# # reading_hd5file  = rhdf5::H5Fopen(name = h5file)
# # readin_gene_norm = reading_hd5file&'gene_norm'
# # Y = readin_gene_norm[get_importance_rank(AY_row$Y),     c(A_idx, NT_idx)]
# # 
# # Z = readin_gene_norm[sapply(A = AYZW_cur$Z_names,
# #                             FUN = get_importance_rank), c(A_idx, NT_idx)]
# # 
# # W = readin_gene_norm[sapply(A = AYZW_cur$W_names,
# #                             FUN = get_importance_rank), c(A_idx, NT_idx)]
# # rhdf5::h5closeAll()
# # invisible(gc(verbose=FALSE))
# 
# # Get Y and Z and W 
# # using already loaded in gene_norm (faster)
# # -------------------------------------------    
# Y = gene_norm[get_importance_rank(AY_row$Y),     c(A_idx, NT_idx)]
# 
# if(!pooled) { # individual genes' expressions (requires different AYZW format)
#   Z = gene_norm[sapply(X = AYZW_cur$Z_names, 
#                        FUN = get_importance_rank), c(A_idx, NT_idx)]
#   W = gene_norm[sapply(X = AYZW_cur$W_names, 
#                        FUN = get_importance_rank), c(A_idx, NT_idx)]
#   dfZ = data.frame(t(Z))
#   dfW = data.frame(t(W))
# } else { # pool together genes (requires different AYZW format)
#   mapZ = 
#     # dataframe of Z gene names, chrs, NCidx
#     data.frame(Y = AYZW_cur$Z_names, chr = AYZW_cur$Z_chrs, NCidx = AYZW_cur$Z_NCidx) |> 
#     # list of pooled avg gene expr
#     group_by(chr, NCidx) |> 
#     group_map(~ {data.frame(pooled = gene_norm[sapply(X = .x$Y, 
#                                                       FUN = get_importance_rank), 
#                                                c(A_idx, NT_idx)] 
#                             |> colMeans())                             })
#   # combine pooled to wide df format
#   dfZ = dplyr::bind_cols(lapply(mapZ, as.data.frame.list), .name_repair = 'unique_quiet')
#   
#   mapW = 
#     # dataframe of W gene names, chrs, NCidx
#     data.frame(Y = AYZW_cur$W_names, chr = AYZW_cur$W_chrs, NCidx = AYZW_cur$W_NCidx) |> 
#     # list of pooled avg gene expr
#     group_by(chr, NCidx) |> 
#     group_map(~ {data.frame(pooled = gene_norm[sapply(X = .x$Y, 
#                                                       FUN = get_importance_rank), 
#                                                c(A_idx, NT_idx)] 
#                             |> colMeans())                             })
#   # combine pooled to wide df format
#   dfW = dplyr::bind_cols(lapply(mapW, as.data.frame.list), .name_repair = 'unique_quiet')
# }
# 
# if(pca) {
#   dfZ = getPCArotations(dfZ, 5)
#   dfW = getPCArotations(dfW, 5)
#   
#   # colnames(dfZ) = paste0('Z', 1:ncol(dfZ))
#   # colnames(dfW) = paste0('W', 1:ncol(dfW))
#   # head(W_pca); head(Z_pca)
# }
# 
# # Assemble together
# # -------------------------------------------
# dfAY = data.frame(A = as.vector(A), 
#                   Y = Y)
# colnames(dfZ) = paste0('Z', 1:ncol(dfZ))
# colnames(dfW) = paste0('W', 1:ncol(dfW))
# 
# 
# final = cbind(dfAY, dfZ, dfW)
# 
# # Z_pca = prcomp(Z, center = TRUE, scale. = TRUE)
# # W_pca = prcomp(W, center = TRUE, scale. = TRUE)
# # 
# # Z_ = Z_pca$x[,1:r_z, drop=FALSE]
# # W_ = W_pca$x[,1:r_w, drop=FALSE]
# # colnames(Z_) = paste0('Z',1:ncol(Z_))
# # colnames(W_) = paste0('W',1:ncol(W_))
# 
# # dfZW = cbind(dfZ, dfW)
# # rm(dfZ, dfW, Z, W) # clean up
# # gc(verbose = FALSE)


# === === === === === === === === === === === === === === === === === 
# ---------------
# get_ATE_est0_nopca = get_ATE_est_make(AY=AY, AYZW=AYZW, gene_norm=gene_norm, grna_rownames=grna_rownames, 
#                                       NT_idx=NT_idx, imp_gene_names=imp_gene_names,
#                                       CB_setting=CB_setting, run_OCBLinOSEst=run_OCBLinOSEst,
#                                       save_path = switch('no',
#                                                          'yes' = sprintf('%s/supergene/%s/%s', save_dir, AYZW_setting_name, CB_setting_name),
#                                                          'no'  = NULL),
#                                       pooled=FALSE,
#                                       pca=FALSE)
# 
# get_ATE_est0_pca = get_ATE_est_make(AY=AY, AYZW=AYZW, gene_norm=gene_norm, grna_rownames=grna_rownames, 
#                                     NT_idx=NT_idx, imp_gene_names=imp_gene_names,
#                                     CB_setting=CB_setting, run_OCBLinOSEst=run_OCBLinOSEst,
#                                     save_path = switch('no',
#                                                        'yes' = sprintf('%s/supergene/%s/%s', save_dir, AYZW_setting_name, CB_setting_name),
#                                                        'no'  = NULL),
#                                     pooled=FALSE,
#                                     pca=TRUE)
# 
# ATE_nopca = get_ATE_est0_nopca(AY_idx=1, ZW_idx=1)
# ATE_pca = get_ATE_est0_pca(AY_idx=1, ZW_idx=1)
# 
# head(ATE_nopca, 10); head(ATE_pca, 10)
# ATE_nopca$ATE - ATE_pca$ATE
# ======================================


get_ATE_est0 = get_ATE_est_make(AY=AY, AYZW=AYZW, gene_norm=gene_norm, grna_rownames=grna_rownames, grna=grna, 
                                NT_idx=NT_idx, imp_gene_names=imp_gene_names,
                                CB_setting=CB_setting, run_OCBLinOSEst=run_OCBLinOSEst,
                                save_path = switch(save_intermediateATEs,
                                                   'yes' = sprintf('%s/supergene/%s/%s', save_dir, AYZW_setting_name, CB_setting_name),
                                                   'no'  = NULL),
                                pooled=FALSE,
                                pca=TRUE)


get_ATE_est <- function(AY_idx, ZW_idx) {
  res_df = tryCatch({get_ATE_est0(AY_idx=AY_idx, ZW_idx=ZW_idx)},
                    error = function(cond) {
                      # message(sprintf('Error est CondMomentOCBOS with %s', 
                      #                 gammaSetting)) 
                      return(NULL)
                    })
  # if errored, return NULL
  if(is.null(res_df)) {
    return(NULL)
  } else {
    return(res_df)
  }
}

NUM_NCENCO_per_AY = length(AYZW[[1]][[1]]) # prev defined, will be length of this sublist


ATEargs = expand.grid(AY_idx = 1:nrow(AY), 
                      ZW_idx  = 1:NUM_NCENCO_per_AY) |> 
  arrange(AY_idx, ZW_idx) # rearrange, AY pair first
# lessen the amount...
# ATEargs = expand.grid(AY_idx = c(1:30, 120:134), ZW_idx  = 1:NUM_NCENCO_per_AY)


# NUMROWS = 600 # 16 threads ? mins, 8 threads ? mins
# NUMROWS = 10 # 1 time takes 37 sec, 10 runs seq for ? mins, parallel 8 thrds ? mins
# NUMROWS = nrow(ATEargs)
# whichROWS = 300:nrow(ATEargs)
# whichROWS = 1:16
# whichROWS = 1:NUMROWS
whichROWS = 601:nrow(ATEargs)

# # =================== Get ATEs (parallel) ====================================
print(sprintf("[%s]    - Get ATEs (parallel)", Sys.time()))
t0 = Sys.time()
ATE_par = future.apply::future_mapply(get_ATE_est,
                                      AY_idx = ATEargs[whichROWS, 1], # ATEargs[1:NUMROWS, 1],
                                      ZW_idx = ATEargs[whichROWS, 2], # ATEargs[1:NUMROWS, 2],
                                      future.globals = TRUE,
                                      future.seed = 56789)
# # manually state globals
# future.globals = c('AY', 'AYZW', 'grna_rownames', 'grna', 
#                    'NT_idx', 'get_importance_rank', 'gene_norm', 
#                    'format_AYZW', 'get_CB_colnames', 'get_CB_est', 
#                    'CB', 'get_lmYA_est', 'get_ATE_est'))
t1 = Sys.time()
print(sprintf("[%s]        - %2.2f", Sys.time(), (t1 - t0)))































