# ---------------------------------------------------------------------------- #
#                   Active p-values (using sPCA of genes as NCE/NCO)           #
# Requires: prev saved normalized gene expression (HDF5)                       #
#           prev saved specified AY pairs
#           prev saved sPCA PCs (for NCE/NCO)
# Ouputs: (nothing) but saves                                                  #
# ---------------------------------------------------------------------------- #
args = commandArgs(trailingOnly = TRUE)
# args = c('laptop', 'allPerturbations')
# args = c('laptop', 'A')

require(assertthat) # for some assert statements
library(Matrix)
library(dplyr)
library(tibble)
library(ggplot2)    # plotting
library(cowplot)


library(future.apply)
options(future.globals.maxSize= 850*1024^2) #1st num is MB
# plan(multisession, workers = 8)
plan(sequential)

# library(furrr)
# plan(multisession, workers = 2)

theme_set(theme_cowplot() +
            theme(plot.title = element_text(hjust = .5),
                  plot.subtitle = element_text(hjust = .5)))

# === Parameter Settings for CB Methods ===
ALPHA = .1
GMM_steps = 10
K_folds = 5
lambdas = (10)^(seq(from = -10, to = 3, length.out = 27))
TRIM = .1
# num_NC_pairs = c(1:8, 10, 15, 20, 25, 30) # number of NC pairs to try (<= total NC pairs, will be cut off)
# num_NC_pairs = c(1, 3, 5, 10, 20)
num_NC_pairs = c(20)
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
# run_OCBLinOSEst = FALSE # whether to run OCBLinOne-Step estimator or not
save_intermediateATEs = 'yes' # 'yes'/'no' whether to save intermedate ATEs as they are estimated

# === Parameter Settings from SPCA
my_sumabsv = 5
my_K = 60
N_subsample = 'all' # subsample size, or 'all' if using all cells
N_subsample = 2000

# === Parameter Settings for which estimators to perform
which_estimators = list(lm_YA        = TRUE,
                        OCB_2SLS     = FALSE,
                        OCB_2SLS_pci2s=TRUE,
                        OCB_2SLSReg  = FALSE,
                        OCB_GMM      = FALSE,
                        OCB_GMMRw    = FALSE,
                        OCB_GMMRwReg = FALSE,
                        OCB_LinOSPI  = FALSE,
                        OCB_LinOS    = FALSE,
                        OCB_LinOStrim= FALSE)


# === Param setting for active pval
GAMMA = 0

# === === === === === === === === === ===

assertthat::assert_that(length(args) > 0, msg="must give arg for specifying device eg 'Rscript <filename>.R ubergenno'")

DEVICE = args[1]
source('../PATHS.R') # load in data_dir and save_dir and CODE_DIR, depending on DEVICE value
assertthat::assert_that(!is.null(data_dir), msg='first arg must be: laptop, desktop, or ubergenno')

assertthat::assert_that(length(args) > 1, msg="must give arg for specifying chosen AYZW name 'Rscript <filename>.R ubergenno C'")
AYZW_setting_name = args[2]


# save parameter setting for basis + alpha
CB_setting = list() # new copy to add ALPHA to
CB_setting$ALPHA     = ALPHA # R does not save as pointers, so ok
CB_setting$GMM_steps = GMM_steps
CB_setting$K_folds   = K_folds
CB_setting$lambdas   = lambdas
CB_setting$TRIM      = TRIM
CB_setting$basisParams  = basisParams
CB_setting$num_NC_pairs = num_NC_pairs
# CB_setting$run_OCBLinOSEst = run_OCBLinOSEst
dir.create(sprintf('%s/spca/cbgenes/%s/%s', save_dir, AYZW_setting_name, CB_setting_name), recursive = TRUE, showWarnings = FALSE)
capture.output(print(CB_setting),
               file = sprintf('%s/spca/cbgenes/%s/%s/CB_setting.txt',
                              save_dir, AYZW_setting_name, CB_setting_name))
saveRDS(CB_setting,
        sprintf('%s/spca/cbgenes/%s/%s/CB_setting.rds',
                save_dir, AYZW_setting_name, CB_setting_name))

# =================== Start ====================================================
print(sprintf("[%s] START: CB Estimate", Sys.time()))

# # source CB utility functions (mainly CB fn for estimating)
# source(sprintf('%s/CB_utils.R', util_dir))
# Instead source new estimation methods' files
# source(sprintf('%s/OCB2SLS.R', util_dir))      # 2-Stage Least Squares
# source(sprintf('%s/OCBGMM.R', util_dir))       # Gen. Method of Moments
# source(sprintf('%s/OCBLinBridge.R', util_dir)) # Lin Bridge Plug-In and One-Step

# source(sprintf('%s/getdfPCA.R', util_dir))     # for performing PCA

source(sprintf('%s/constructDataList.R', util_dir)) # for formatting data for other fns
source(sprintf('%s/CBEstAll.R', util_dir)) # for functions to estimate here
source(sprintf('%s/CBEstAllSPCA.R', util_dir)) # add'l specifically for SPCA

source(sprintf('%s/CBEstSPCAActive.R', util_dir)) # for active arbdep with SPCA

# load chosen AYZW names
AY   = read.csv(sprintf('%s/spca/cbgenes/%s/AY.csv', save_dir, AYZW_setting_name))
# AYZW = readRDS(sprintf('%s/cbgenes/%s/AYZW.rds', save_dir, AYZW_setting_name))

# # Only for PCA: copy to pca folder too (bc plots script reads in AY, AYZW)
# write.csv(AY, sprintf('%s/spca/cbgenes/%s/AY.csv', save_dir, AYZW_setting_name))
# saveRDS(AYZW, sprintf('%s/spca/cbgenes/%s/AYZW.rds', save_dir, AYZW_setting_name))

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

# Load choosing AY (ZW) settings
AYZW_setting = readRDS(sprintf('%s/spca/cbgenes/%s/AYZW_setting.rds', save_dir, AYZW_setting_name))
NUM_IMPORTANTGENES = AYZW_setting$MAX_Y_IMPORTANCE

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
gene_norm = readin_gene_norm[1:NUM_IMPORTANTGENES, ] # eg dim = 4000 x 20729 = #important x #cells
rhdf5::h5closeAll()
invisible(gc(verbose=FALSE))




# use Sparse PCA Loadings
if(N_subsample == 'all') {  N_subsample = ncol(gene_odm) }
NCs = readRDS(sprintf('%s/spca/NCloadings_sumabs=%.1f_K=%d_N=%d.rds', save_dir, my_sumabsv, my_K, N_subsample)) # save 
# Or the constructed clusters (averages)
# NCs = readRDS(sprintf('%s/spca/NCavg_sumabs=%.1f_K=%d_N=%d.rds', save_dir, my_sumabsv, my_K, N_subsample)) # save
# NCs = data.frame(NCs)
# colnames(NCs) = paste0('NC', 1:ncol(NCs))


# =================== Define File Specific Functions ======================================
print(sprintf("[%s]    - Define File Specific Functions", Sys.time()))




# All the NT grna idx
NT_names = grna_odm |> ondisc::get_feature_covariates() |>
  filter(target_type == 'non-targeting') |> rownames()
NT_idx = which(apply(X = grna_odm[[NT_names, ]], MARGIN = 2, FUN = sum) > 0)






# =================== Make Make Function ======================================

get_proxy_pval0 =  get_proxy_pval_make(AY=AY, gene_norm=gene_norm, NCs=NCs, grna_rownames=grna_rownames, grna=grna, 
                                    NT_idx=NT_idx, imp_gene_names=imp_gene_names,
                                    # which_estimators=which_estimators, 
                                    CB_setting=CB_setting, 
                                    save_path=NULL)

get_true_pval0  =  get_true_pval_make(AY=AY, gene_norm=gene_norm, NCs=NCs, grna_rownames=grna_rownames, grna=grna, 
                                    NT_idx=NT_idx, imp_gene_names=imp_gene_names,
                                    # which_estimators=which_estimators, 
                                    CB_setting=CB_setting, 
                                    save_path=NULL)
# get_proxy_pval0(1)
# get_true_pval0(1)

get_active_pval0 = get_active_arbdep_pval_make(get_proxy_pval=get_proxy_pval0,
                                                get_true_pval=get_true_pval0,
                                               # gamma = .5)
                                               gamma = GAMMA) # tune the querying for true pval 

# get_active_pval0(1)
# get_active_pval(1)


# another fn to ignore errors (so that the others may continue running)
get_active_pval <- function(AY_idx) {
  res = tryCatch({get_active_pval0(AY_idx)},
                    error = function(cond) {
                      # message(sprintf('Error est CondMomentOCBOS with %s', 
                      #                 gammaSetting)) 
                      return(NULL)
                    })
  # if errored, return NULL
  if(is.null(res)) {
    return(NULL)
  } else {
    return(res)
  }
}


# get_active_pval(1)

# get_ATE_est0 = get_ATE_est_NCs_make(AY=AY, gene_norm=gene_norm, NCs=NCs, grna_rownames=grna_rownames, grna=grna, 
#                                     NT_idx=NT_idx, imp_gene_names=imp_gene_names,
#                                     which_estimators=which_estimators, CB_setting=CB_setting, 
#                                     save_path=switch(save_intermediateATEs,
#                                                      'yes' = sprintf('%s/spca/cbgenes/%s/%s', save_dir, AYZW_setting_name, CB_setting_name),
#                                                      'no'  = NULL))
# test = get_ATE_est0(AY_idx = 3)



# get_ATE_est <- function(AY_idx) {
#   res_df = tryCatch({get_ATE_est0(AY_idx=AY_idx)},
#                     error = function(cond) {
#                       # message(sprintf('Error est CondMomentOCBOS with %s', 
#                       #                 gammaSetting)) 
#                       return(NULL)
#                     })
#   # if errored, return NULL
#   if(is.null(res_df)) {
#     return(NULL)
#   } else {
#     return(res_df)
#   }
# }

# NUM_NCENCO_per_AY = length(AYZW[[1]][[1]]) # prev defined, will be length of this sublist
# ATEargs = expand.grid(AY_idx = 1:nrow(AY), 
#                       ZW_idx  = 1:NUM_NCENCO_per_AY) |> 
#   arrange(AY_idx, ZW_idx) # rearrange, AY pair first

ATEargs = data.frame(AY_idx = 1:nrow(AY))
# lessen the amount...
# ATEargs = expand.grid(AY_idx = c(1:30, 120:134), ZW_idx  = 1:NUM_NCENCO_per_AY)



# NUMROWS = 10
NUMROWS = nrow(ATEargs)
# whichROWS = 1:1000
whichROWS = 1:3
# whichROWS = 1:NUMROWS
# whichROWS = 1165:NUMROWS


# !!! somehow, the parallel version is getting messed up (all the true pvals = proxy pvals exactly)
# # =================== Get ATEs (parallel) ====================================
if(F) {
  print(sprintf("[%s]    - Get ATEs (parallel)", Sys.time()))
  t0 = Sys.time()
  pvals_par = future.apply::future_mapply(get_active_pval,
                                          AY_idx = ATEargs[whichROWS, 1], # ATEargs[1:NUMROWS, 1],
                                          future.globals = TRUE,
                                          future.seed = 56789)
  # # manually state globals
  # future.globals = c('AY', 'AYZW', 'grna_rownames', 'grna', 
  #                    'NT_idx', 'get_importance_rank', 'gene_norm', 
  #                    'format_AYZW', 'get_CB_colnames', 'get_CB_est', 
  #                    'CB', 'get_lmYA_est', 'get_ATE_est'))
  t1 = Sys.time()
  print(sprintf("[%s]        - %2.2f", Sys.time(), (t1 - t0)))
  
  
  t(pvals_par)
  # columns are lists...
  pvals_df = apply(pvals_par, FUN = unlist, MARGIN=1) |> as.data.frame()
  pvals_df$true_prob = 1 - pvals_df$gamma * pvals_df$pval_proxy # add col of prob of querying true
  
  
  saveRDS(pvals_df,
          sprintf('%s/spca/cbgenes/%s/%s/ATE_activearbdep_gamma=%.2f.rds', save_dir, AYZW_setting_name, CB_setting_name, GAMMA))
  write.csv(x = pvals_df,
            row.names = FALSE, 
            file = sprintf('%s/spca/cbgenes/%s/%s/ATE_activearbdep_gamma=%.2f.csv', save_dir, AYZW_setting_name, CB_setting_name, GAMMA))
  
  
  
  # ATE_par[, 1]
  # ATE_par[[1, ]]
  # length(ATE_par)
  # unlist(ATE_par)
  # ATE_par[1] |> as.data.frame()
}


# # =================== Get ATEs (sequential) ==================================
print(sprintf("[%s]    - Get ATEs (sequential)", Sys.time()))
t0 = Sys.time()

pvals_seq = NULL
for(AY_idx in ATEargs[whichROWS, 1]) {
  pvals_seq = rbind(pvals_seq, 
                    get_active_pval(AY_idx) |> data.frame())
}


t1 = Sys.time()
print(sprintf("[%s]        - %2.2f", Sys.time(), (t1 - t0)))

saveRDS(pvals_seq,
        sprintf('%s/spca/cbgenes/%s/%s/ATE_activearbdep_gamma=%.2f.rds', save_dir, AYZW_setting_name, CB_setting_name, GAMMA))
write.csv(x = pvals_seq,
          row.names = FALSE, 
          file = sprintf('%s/spca/cbgenes/%s/%s/ATE_activearbdep_gamma=%.2f.csv', save_dir, AYZW_setting_name, CB_setting_name, GAMMA))



