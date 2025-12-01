# ---------------------------------------------------------------------------- #
#                   Do CB using genes as NCE/NCO                               #
# 3.1 - add chromosome number information                                      #
# 3.2 - choose A,Y,Z,W combinations                                            #
# 3.2 - CB Effect Estimate                                                     #
# Requires: prev saved normalized gene expression (HDF5)                       #
# Ouputs: (nothing) but saves                                                  #
#         CBGENE_AYZW in the rds file                                          #
#                         "<save_dir>/CBGENE_AYZW.rds"                         #
# ---------------------------------------------------------------------------- #
args = commandArgs(trailingOnly = TRUE)
# args = c('laptop', '2000', 'C')

require(assertthat) # for some assert statements
library(Matrix)
library(dplyr)
library(tibble)
library(ggplot2)    # plotting
library(cowplot)


library(future.apply)
options(future.globals.maxSize= 850*1024^2) #1st num is MB
plan(multisession, workers = 8)

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
# num_NC_pairs = c(1:3, 5, 10, 25, 50, 75, 100)
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
# CB_setting_name = 'default'
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
dir.create(sprintf('%s/cbgenes/%s/%s', save_dir, AYZW_setting_name, CB_setting_name))
capture.output(print(CB_setting), 
               file = sprintf('%s/cbgenes/%s/%s/CB_setting.txt', 
                              save_dir, AYZW_setting_name, CB_setting_name))
saveRDS(CB_setting,
        sprintf('%s/cbgenes/%s/%s/CB_setting.rds', 
                save_dir, AYZW_setting_name, CB_setting_name))



# =================== Start ====================================================
print(sprintf("[%s] START: CB Estimate", Sys.time()))

# # source CB utility functions (mainly CB fn for estimating)
# source(sprintf('%s/CB_utils.R', util_dir))
# Instead source new estimation methods' files
source(sprintf('%s/OCB2SLS.R', util_dir))      # 2-Stage Least Squares
source(sprintf('%s/OCBGMM.R', util_dir))       # Gen. Method of Moments
source(sprintf('%s/OCBLinBridge.R', util_dir)) # Lin Bridge Plug-In and One-Step



source(sprintf('%s/CBEstAll.R', util_dir)) # for functions to estimate here



# load chosen AYZW names
AY   = read.csv(sprintf('%s/cbgenes/%s/AY.csv', save_dir, AYZW_setting_name))
AYZW = readRDS(sprintf('%s/cbgenes/%s/AYZW.rds', save_dir, AYZW_setting_name))


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



# =================== TODO: CONVERT THESE TO UTIL FUNCTIONS! NOT FILE SPECIFIC
# ==================== e.g. does not use loaded in values
# =================== Define File Specific Functions ======================================
print(sprintf("[%s]    - Define File Specific Functions", Sys.time()))




# All the NT grna idx
NT_names = grna_odm |> ondisc::get_feature_covariates() |> 
          filter(target_type == 'non-targeting') |> rownames()
NT_idx = which(apply(X = grna_odm[[NT_names, ]], MARGIN = 2, FUN = sum) > 0)


get_ATE_est0 = get_ATE_est_make(AY=AY, AYZW=AYZW, gene_norm=gene_norm, grna_rownames=grna_rownames, grna=grna,
                               NT_idx=NT_idx, imp_gene_names=imp_gene_names,
                               CB_setting=CB_setting, run_OCBLinOSEst=run_OCBLinOSEst,
                               save_path = switch(save_intermediateATEs,
                                                  'yes' = sprintf('%s/cbgenes/%s/%s', save_dir, AYZW_setting_name, CB_setting_name),
                                                  'no'  = NULL))


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

 
# testdf = get_ATE_est(AY_idx=3, ZW_idx=2)
# res_df = testdf
# res_df = ATE_par[,1] |> data.frame()
# head(res_df); dim(res_df)
# naive_ate = res_df |> filter(method == 'lmYA') |> pull(ATE)
# # 2SLS Estimates
# ggplot(res_df |> filter(method_type == '2SLS'),
#        aes(x = numNC, y = ATE, color = method, group = method)) +
#   geom_line() +
#   geom_hline(aes(yintercept = naive_ate), color = 'red', alpha = .7, linetype = 'dashed')
# 
# # GMM Estimates
# ggplot(res_df |> filter(method_type == 'GMM'),
#        aes(x = numNC, y = ATE, color = method, group = method)) +
#   geom_line() +
#   geom_hline(aes(yintercept = naive_ate), color = 'red', alpha = .7, linetype = 'dashed') +
#   facet_grid(rows = vars(basis))
# 
# # OCBLin Estimates
# ggplot(res_df |> filter(method_type == 'OCBLin'),
#        aes(x = numNC, y = ATE, color = method, group = method)) +
#   geom_line() +
#   geom_hline(aes(yintercept = naive_ate), color = 'red', alpha = .7, linetype = 'dashed') +
#   facet_grid(rows = vars(basis))


NUM_NCENCO_per_AY = length(AYZW[[1]][[1]]) # prev defined, will be length of this sublist


ATEargs = expand.grid(AY_idx = 1:nrow(AY), 
                      ZW_idx  = 1:NUM_NCENCO_per_AY) |> 
          arrange(AY_idx, ZW_idx) # rearrange, AY pair first
# lessen the amount...
# ATEargs = expand.grid(AY_idx = c(1:30, 120:134), ZW_idx  = 1:NUM_NCENCO_per_AY)

# NUMROWS = 16 # 16 runs on 8 threads took 7 mins
# NUMROWS = 8 # 8 runs with fewer large dim took 1 min?
NUMROWS = nrow(ATEargs)
# whichROWS = 200:300 # testing which gave error
# whichROWS = 200:250 # testing which gave error
# whichROWS = 201:207 # testing which gave error
# whichROWS = 1:200 # just do a few for now
# whichROWS = 300:nrow(ATEargs)
# whichROWS = 1:16
whichROWS = 1:NUMROWS

# =================== Get ATEs (not parallel) ======================================
# print(sprintf("[%s]    - Get ATEs (not parallel)", Sys.time()))
# 
# t0 = Sys.time() 
# ATE_seq = mapply(get_ATE_est, 
#              AY_idx = ATEargs[1:NUMROWS, 1], 
#              ZW_idx = ATEargs[1:NUMROWS, 2])
# t1 = Sys.time()
# print(sprintf("[%s]        - %2.2f", Sys.time(), (t1 - t0)))
# 
# 
# write.csv(cbind(AY_idx = ATEargs[1:NUMROWS, 1], 
#                 ZW_idx = ATEargs[1:NUMROWS, 2],
#                 t(ATE_seq)), 
#           sprintf('%s/cbgenes/%s/ATE.csv', save_dir, AYZW_setting_name), row.names = FALSE)


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


# combine into one dataframe
# formatting...
# if no errors, then returns a list indexing like ATE_par[, i].
# if errors, then returns a list indexing like ATE_par[[i]]
ATE_df = NULL
if(is.null(ncol(ATE_par))) { # had errors
  # ATE_par_nulls = ATE_par
  
  for(i in 1:length(ATE_par)) {
    if(is.null(ATE_par[[i]])) {
      ATE_df = dplyr::bind_rows(ATE_df, 
                                data.frame(AY_idx = ATEargs[whichROWS[i], 1], 
                                           ZW_idx = ATEargs[whichROWS[i], 2]))
    } else {
      ATE_df = rbind(ATE_df, 
                     cbind(data.frame(AY_idx = ATEargs[whichROWS[i], 1], 
                                      ZW_idx = ATEargs[whichROWS[i], 2]),
                           ATE_par[[i]] |> data.frame()))
    }
  }
} else { # no errors
  for(i in 1:ncol(ATE_par)) {
    ATE_df = rbind(ATE_df, 
                   cbind(data.frame(AY_idx = ATEargs[whichROWS[i], 1], 
                                    ZW_idx = ATEargs[whichROWS[i], 2]),
                         ATE_par[,i] |> data.frame()))
    
    
  }
}


head(ATE_df); dim(ATE_df)


write.csv(ATE_df, 
          sprintf('%s/cbgenes/%s/%s/ATE.csv', save_dir, AYZW_setting_name, CB_setting_name), row.names = FALSE)


# =================== END ======================================================
print(sprintf("[%s] END", Sys.time()))
















# =================== Trash ======================================================
# #' important_genes_name.rds should be in order of most to least important
# #' the rank will be the index of the gene on this list
# get_importance_rank <- function(gene_name) {
#   which(imp_gene_names == gene_name)
# }

#' #' important_genes_name.rds should be in order of most to least important
#' #' the rank will be the index of the gene on this list
#' #' Maybe dont need this...
#' #' Make function that gets the importance ranking of specified gene_name using
#' #' the imput imp_gene_names parameter. The imp_gene_names should be ordered by 
#' #' importance, from most important to least important.
#' #' @param imp_gene_names (vector) of strings of gene names ordered by importance
#' #' @param gene_names (string) gene name
#' #' @example 
#' #' # suppose you want a function that returns the importance ranking according 
#' #' to myImpGeneNames 
#' #' myGetImpRank = get_importance_rank_make(myImpGeneNames)
#' #' myGetImpRank(myGeneName) # get the imp ranking of myGeneName
#' get_importance_rank_make <- function(imp_gene_names) {
#'   get_importance_rank <- function(gene_name) {
#'     which(imp_gene_names == gene_name)
#'   }
#'   return(get_importance_rank)
#' }


# get_importance_rank = get_importance_rank_make(imp_gene_names)   





# #' prep data/formatting
# #' creates a dataframe with
# #' A (boolean) of getting specified gRNA or any NT perturbation
# #' Y, Z1, ..., W1, ... (numeric) normalized gene expressions 
# format_AYZW <- function(AY_idx, ZW_idx) {
#   # AY_idx = 3; ZW_idx = 2
#   AY_row = AY[AY_idx, ]

#   # gather the AYZW names chosen
#   AYZW_cur = AYZW[[AY_row$A]][[AY_row$Y]][[ZW_idx]]


#   # AY_row$A
#   # AY_row$Y
#   # AYZW_cur$Z_names
#   # AYZW_cur$W_names

#   # Get A
#   # -------------------------------------------
#   # # idx of all 'treated' cells
#   # A_idx = which(as.logical(grna_odm[[AY_row$A, ]]))
#   # # subset cells of 'treated' (w A grna) and 'control' (NT gran)
#   # A = grna_odm[[AY_row$A, c(A_idx, NT_idx)]]

#   # using grna all loaded into memory
#   # idx of all 'treated' cells

#   A_grna_idx = which(grna_rownames == AY_row$A)
#   A_idx = which(as.logical(grna[A_grna_idx, 1:ncol(grna)]))
#   # subset cells of 'treated' (w A grna) and 'control' (NT gran)
#   A = grna[A_grna_idx, c(A_idx, NT_idx)]





#   # # Get Y and Z and W
#   # # using h5 load in (use if not enough memory)
#   # # -------------------------------------------
#   # h5file      = paste0(save_dir, "/gene.h5")
#   # reading_hd5file  = rhdf5::H5Fopen(name = h5file)
#   # readin_gene_norm = reading_hd5file&'gene_norm'
#   # Y = readin_gene_norm[get_importance_rank(AY_row$Y),     c(A_idx, NT_idx)]
#   # 
#   # Z = readin_gene_norm[sapply(A = AYZW_cur$Z_names,
#   #                             FUN = get_importance_rank), c(A_idx, NT_idx)]
#   # 
#   # W = readin_gene_norm[sapply(A = AYZW_cur$W_names,
#   #                             FUN = get_importance_rank), c(A_idx, NT_idx)]
#   # rhdf5::h5closeAll()
#   # invisible(gc(verbose=FALSE))

#   # Get Y and Z and W 
#   # using already loaded in gene_norm (faster)
#   # -------------------------------------------
#   Y = gene_norm[get_importance_rank(AY_row$Y),     c(A_idx, NT_idx)]

#   Z = gene_norm[sapply(A = AYZW_cur$Z_names, 
#                        FUN = get_importance_rank), c(A_idx, NT_idx)]
#   W = gene_norm[sapply(A = AYZW_cur$W_names, 
#                        FUN = get_importance_rank), c(A_idx, NT_idx)]


#   # Assemble together
#   # -------------------------------------------
#   dfAY = data.frame(A = as.vector(A), 
#                     Y = Y)
#   dfZ = data.frame(t(Z)); colnames(dfZ) = paste0('Z', 1:ncol(dfZ))
#   dfW = data.frame(t(W)); colnames(dfW) = paste0('W', 1:ncol(dfW))

#   return(cbind(dfAY, dfZ, dfW))

# }


#' 
#' #' prep data/formatting
#' #' creates a dataframe with
#' #' A (boolean) of getting specified gRNA or any NT perturbation
#' #' Y, Z1, ..., W1, ... (numeric) normalized gene expressions 
#' #' @param get_importance_rank (function) that takes in arg gene name (str) and 
#' #'            returns the importance ranking (which is used as the idx to get 
#' #'            values from gene_norm)
#' format_AYZW <- function(AY_idx, ZW_idx, 
#'                         AY, AYZW, 
#'                         gene_norm,
#'                         grna_rownames, #imp_gene_names
#'                         NT_idx,
#'                         get_importance_rank) {
#'   # get_importance_rank = get_importance_rank_make(imp_gene_names)   
#'   # AY_idx = 3; ZW_idx = 2
#'   AY_row = AY[AY_idx, ]
#'   
#'   
#'   # gather the AYZW names chosen
#'   AYZW_cur = AYZW[[AY_row$A]][[AY_row$Y]][[ZW_idx]]
#'   
#'   
#'   # AY_row$A
#'   # AY_row$Y
#'   # AYZW_cur$Z_names
#'   # AYZW_cur$W_names
#'   
#'   # Get A
#'   # -------------------------------------------
#'   # # idx of all 'treated' cells
#'   # A_idx = which(as.logical(grna_odm[[AY_row$A, ]]))
#'   # # subset cells of 'treated' (w A grna) and 'control' (NT gran)
#'   # A = grna_odm[[AY_row$A, c(A_idx, NT_idx)]]
#' 
#'   # using grna all loaded into memory
#'   # idx of all 'treated' cells
#'   
#'   A_grna_idx = which(grna_rownames == AY_row$A)
#'   A_idx = which(as.logical(grna[A_grna_idx, 1:ncol(grna)]))
#'   # subset cells of 'treated' (w A grna) and 'control' (NT gran)
#'   A = grna[A_grna_idx, c(A_idx, NT_idx)]
#'   
#'   
#'   
#'   
#'   
#'   # # Get Y and Z and W
#'   # # using h5 load in (use if not enough memory)
#'   # # -------------------------------------------
#'   # h5file      = paste0(save_dir, "/gene.h5")
#'   # reading_hd5file  = rhdf5::H5Fopen(name = h5file)
#'   # readin_gene_norm = reading_hd5file&'gene_norm'
#'   # Y = readin_gene_norm[get_importance_rank(AY_row$Y),     c(A_idx, NT_idx)]
#'   # 
#'   # Z = readin_gene_norm[sapply(A = AYZW_cur$Z_names,
#'   #                             FUN = get_importance_rank), c(A_idx, NT_idx)]
#'   # 
#'   # W = readin_gene_norm[sapply(A = AYZW_cur$W_names,
#'   #                             FUN = get_importance_rank), c(A_idx, NT_idx)]
#'   # rhdf5::h5closeAll()
#'   # invisible(gc(verbose=FALSE))
#'   
#'   # Get Y and Z and W 
#'   # using already loaded in gene_norm (faster)
#'   # -------------------------------------------
#'   Y = gene_norm[get_importance_rank(AY_row$Y),     c(A_idx, NT_idx)]
#'   
#'   Z = gene_norm[sapply(X = AYZW_cur$Z_names, 
#'                        FUN = get_importance_rank), c(A_idx, NT_idx)]
#'   W = gene_norm[sapply(X = AYZW_cur$W_names, 
#'                        FUN = get_importance_rank), c(A_idx, NT_idx)]
#'   
#'   
#'   # Assemble together
#'   # -------------------------------------------
#'   dfAY = data.frame(A = as.vector(A), 
#'                     Y = Y)
#'   dfZ = data.frame(t(Z)); colnames(dfZ) = paste0('Z', 1:ncol(dfZ))
#'   dfW = data.frame(t(W)); colnames(dfW) = paste0('W', 1:ncol(dfW))
#'   
#'   return(cbind(dfAY, dfZ, dfW))
#'   
#' }



#' 
#' #' makes a function that takes in AY_idx, ZW_idx to calculate ATE from a variety
#' #' of estimation methods, using the input 
#' #' @param CB_setting (list) of parameter setttings needed for CB e.g. ALPHA, GMM_steps, ...
#' get_ATE_est_make <- function(AY, AYZW, gene_norm, grna_rownames, NT_idx=NT_idx, imp_gene_names,
#'                              CB_setting) {
#'   get_importance_rank = get_importance_rank_make(imp_gene_names)  
#' 
#'   ALPHA     = CB_setting$ALPHA # R does not save as pointers, so ok
#'   GMM_steps = CB_setting$GMM_steps
#'   K_folds   = CB_setting$K_folds
#'   lambdas   = CB_setting$lambdas
#'   TRIM      = CB_setting$TRIM
#'   basisParams  = CB_setting$basisParams
#'   num_NC_pairs = CB_setting$num_NC_pairs
#'   
#'   #' returns ATE estimates from different estimation methods
#'   get_ATE_est <- function(AY_idx, ZW_idx) {
#'     # AY_idx = 3; ZW_idx = 2
#'     # === Construct df ===
#'     df_all   = format_AYZW(AY_idx=AY_idx, ZW_idx=ZW_idx, # inputs to this inner fn
#'                            AY=AY, AYZW=AYZW,             # inputs saved from outer fn
#'                            gene_norm=gene_norm,
#'                            grna_rownames=grna_rownames,
#'                            NT_idx=NT_idx,
#'                            get_importance_rank=get_importance_rank)
#'     df_all$A = as.numeric(df_all$A) # convert trtmt A to numeric 0/1 if not already (alt T/F)
#'     
#'     # === Get different estimates === 
#'     
#'     
#'     # === === === === === === === === ===  ===
#'     # === Oracle (none) and Naive (no adj) ===
#'     # === === === === === === === === ===  ===
#'     # # Oracle linear model Y ~ A + Us (NOT HERE IN REAL DATA? SHOULD WE TRY?)
#'     # U_colnames = grep('U', colnames(df), value = T)
#'     # lm_YAU  = lm(paste0('Y ~ A + ', paste(U_colnames, collapse = ' + ')), df)
#'     
#'     # === linear model Y ~ A (no confounder adj)
#'     lm_YA = lm(Y ~ A, df_all)
#'     
#'     res = list(#Oracle     = coef(lm_YAU)[['A']],
#'                 lmYA       = coef(lm_YA )[['A']])
#'     
#'     # t0 = Sys.time()
#'     for(num_NCs in num_NC_pairs) {
#'       print(num_NCs)
#'       df = df_all[, get_CB_colnames(num_NCs)] # subset only these cols
#'       
#'       # === === === === === === === === ===  ===
#'       # === Estimators using original values ===
#'       # === === === === === === === === ===  ===
#'       # === Outcome Confounding Bridge, 2 Stage LS
#'       OCB_2SLS = OCB2SLS(df)$ATE
#'       
#'       # === Outcome Confounding Bridge, 2 Stage LS w/ ElNet Reg
#'       OCB_2SLSReg_res = OCB2SLSReg(df, alpha=ALPHA, returnMiddleDFs=FALSE)
#'       OCB_2SLSReg     = OCB_2SLSReg_res[[length(OCB_2SLSReg_res)]][['ATE']]
#'       
#'       # # Bad but,... select the latest ATE estimate that is not NA (increases by #W's)
#'       # OCB_2SLSReg_res_ATEs = c()
#'       # for(i in 1:length(OCB_2SLSReg_res)) {
#'       #   OCB_2SLSReg_res_ATEs = c(OCB_2SLSReg_res_ATEs,
#'       #                            OCB_2SLSReg_res[[i]][['ATE']])
#'       #   # OCB_2SLSReg_res_ATEsOCB_2SLSReg_res[[i]][['ATE']] |> print()
#'       # }
#'       # plot(OCB_2SLSReg_res_ATEs)
#'       # OCB_2SLSReg = OCB_2SLSReg_res_ATEs[max(which(!is.na(OCB_2SLSReg_res_ATEs)))]
#'       
#'       # === save to res = list of resulting ATEs for this num_NCs in list
#'       res_numNCs = list(OCB2SLS    = OCB_2SLS,
#'                         OCB2SLSReg = OCB_2SLSReg)
#'       
#'       
#'       # === === === === === === === === === === === === ===  ===
#'       # === Estimators with Transformations/Changes of Basis ===
#'       # === === === === === === === === === === === === ===  ===
#'       for(basis in names(basisParams)) {
#'         res_numNCs_basis = list()
#'         # print(basis)
#'         dl = constructDataListv2(df,
#'                                  b_degree =basisParams[[basis]]$b_degree, 
#'                                  h_degree =basisParams[[basis]]$h_degree, type = 'simple',
#'                                  b_phases =basisParams[[basis]]$b_phases, 
#'                                  b_periods=basisParams[[basis]]$b_periods,
#'                                  h_phases =basisParams[[basis]]$h_phases, 
#'                                  h_periods=basisParams[[basis]]$h_periods)
#'         
#'         
#'         # === OCBGMM
#'         res_numNCs_basis[['OCBGMM']] = OCBGMM(dl, returnAlphas = F)$ATE
#'         # res_numNCs[[paste0('OCBGMM_', basis)]] = OCBGMM(dl, returnAlphas = F)$ATE
#'         
#'         # === OCBGMM ReWeight
#'         res_numNCs_basis[['OCBGMMRW']] = OCBGMMRw(dl, GMM_steps, returnWeights = F)$ATE
#'         # res_numNCs[[paste0('OCBGMMRW_', basis)]] = OCBGMMRw(dl, GMM_steps, returnWeights = F)$ATE
#'         
#'         # === OCBGMM ReWeight Regularization
#'         res_numNCs_basis[['OCBGMMRWReg']] = 
#'                            OCBGMMRwReg(dl, GMM_steps, K_folds, lambdas, 
#'                                        returnWeights = F, returnLambdas = F)$ATE
#'         # res_numNCs[[paste0('OCBGMMRWReg_', basis)]] = 
#'         #                    OCBGMMRwReg(dl, GMM_steps, K_folds, lambdas, 
#'         #                                returnWeights = F, returnLambdas = F)$ATE
#'         
#'         # === OCB Linear CB Plug-In AND One-Step 
#'         #      Z0,Z1 NULL --> use transformations in B0,B1 as covariates 
#'         #                     when estimating nuisance params
#'         
#'         # No trimming
#'         OCBLinOS_res = tryCatch({estCondMomentOCBOS(dataList=dl, gmmSteps=GMM_steps, Z0=NULL, Z1=NULL)},
#'                                 error = function(cond) {
#'                                   # message(sprintf('Error est CondMomentOCBOS with %s', 
#'                                   #                 gammaSetting)) 
#'                                   return(NULL)
#'                                 }) 
#'         # OCBLinOS_res = estCondMomentOCBOS(dataList=dl, gmmSteps=GMM_steps, Z0=NULL, Z1=NULL)
#'         # errors sometimes... NA for now... debug later
#'         if(is.null(OCBLinOS_res)) {
#'           res_numNCs_basis[['OCBLinPI']] = NA
#'           res_numNCs_basis[['OCBLinOS']] = NA
#'           # res[[paste0('OCBLinPI_', basis)]] = NA
#'           # res[[paste0('OCBLinOS_', basis)]] = NA
#'         } else {
#'           res_numNCs_basis[['OCBLinPI']] = OCBLinOS_res$ATE_pi
#'           res_numNCs_basis[['OCBLinOS']] = OCBLinOS_res$ATE_os
#'           # res[[paste0('OCBLinPI_', basis)]] = OCBLinOS_res$ATE_pi
#'           # res[[paste0('OCBLinOS_', basis)]] = OCBLinOS_res$ATE_os
#'         }
#'         
#'         # With Trimming (at .1)
#'         OCBLinOS_res = tryCatch({estCondMomentOCBOS(dataList=dl, gmmSteps=GMM_steps, 
#'                                                     Z0=NULL, Z1=NULL, trim=TRIM)},
#'                                 error = function(cond) {
#'                                   # message(sprintf('Error est CondMomentOCBOS with %s', 
#'                                   #                 gammaSetting)) 
#'                                   return(NULL)
#'                                 }) 
#'         if(is.null(OCBLinOS_res)) {
#'           res_numNCs_basis[['OCBLinPI']] = NA
#'           # res[[paste0('OCBLinOStrim_', basis)]] = NA
#'         } else {
#'           res_numNCs_basis[['OCBLinPI']] = OCBLinOS_res$ATE_os
#'           # res[[paste0('OCBLinOStrim_', basis)]] = OCBLinOS_res$ATE_os
#'         }
#'         
#'         res_numNCs[[basis]] = res_numNCs_basis
#'       }
#'       
#'       res[[paste0('NC', num_NCs)]] = res_numNCs # save xx NC pair res to overall res
#'       
#'     }
#'     
#'     # t1 = Sys.time(); print(t1 - t0)
#'     # 1 run took 24.38264 mins...
#'     # saveRDS(res, 
#'     #         file = sprintf('%s/cbgenes/%s/%s/onetry.rds', 
#'     #                        save_dir, AYZW_setting_name, CB_setting_name))
#'     
#'     
#'     # maybe shouldve just saved as a df...
#'     res_df = data.frame(method = 'lmYA', method_type = 'naive', numNC = NA, basis = NA, ATE = res$lmYA)
#'     
#'     for(num_NC in names(res)[-1]) {
#'       num_NC_int = substr(num_NC, start = 3, stop = nchar(num_NC)) |> as.integer()
#'       # num_NC_int = substr(num_NC, start = 1, stop = nchar(num_NC)-2) |> as.integer()
#'       res_df = rbind(res_df,
#'                      data.frame(method = 'OCB2SLS', method_type = '2SLS', 
#'                                 numNC = num_NC_int, basis = NA, ATE = res[[num_NC]][['OCB2SLS']] ))
#'       res_df = rbind(res_df,
#'                      data.frame(method = 'OCB2SLSReg', method_type = '2SLS', 
#'                                 numNC = num_NC_int, basis = NA, ATE = res[[num_NC]][['OCB2SLSReg']] ))
#'       for(basis in names(res[[num_NC]])) {
#'         for(method in names(res[[num_NC]][[basis]])) {
#'           if(method %in% c("OCBGMM", "OCBGMMRW", "OCBGMMRWReg")) {
#'             method_type = 'GMM'
#'           } else if(method %in% c("OCBLinPI", "OCBLinOS")) {
#'             method_type = 'OCBLin'
#'           }
#'           res_df = rbind(res_df,
#'                          data.frame(method = method, method_type = method_type, 
#'                                     numNC = num_NC_int, basis = basis, ATE = res[[num_NC]][[basis]][[method]] ))
#'           
#'         }
#'         
#'       }
#'     }
#'     
#'     
#'     
#'     return(res_df)
#'   }
#'   
#'   
#'   return(get_ATE_est)
#' }

#' #' NEW: TODO: GET ATE from all methods at once
#' #' COPY THE STRUCTURE FROM SIMULATION CODE FILE
#' #' (pass in everything basically... or turn to make function)
#' #' @param AY_idx (integer) which AY test (as index)
#' #' @param ZW_idx (integer) which ZW set of NCE/NCOs  (as index)
#' get_ATE_est <- function(AY_idx, ZW_idx, AY, AYZW, grna_rownames, imp_gene_names) {
#'   # AY_idx = 3; ZW_idx = 2
#'   # AY, AYZW, AY_idx, ZW_idx,
#'   # grna_rownames, imp_gene_names
#'   # df_all = format_AYZW(AY_idx=AY_idx, ZW_idx=ZW_idx,
#'   #                  AY = AY, AYZW = AYZW, 
#'   #                  grna_rownames, imp_gene_names)
#'   
#'   # === Construct df ===
#'   df_all   = format_AYZW(AY_idx=AY_idx, ZW_idx=ZW_idx)
#'   df_all$A = as.numeric(df_all$A) # convert trtmt A to numeric 0/1 if not already (alt T/F)
#'   
#'   # === Get different estimates === 
#'   
#'   
#'   # === === === === === === === === ===  ===
#'   # === Oracle (none) and Naive (no adj) ===
#'   # === === === === === === === === ===  ===
#'   # # Oracle linear model Y ~ A + Us (NOT HERE IN REAL DATA? SHOULD WE TRY?)
#'   # U_colnames = grep('U', colnames(df), value = T)
#'   # lm_YAU  = lm(paste0('Y ~ A + ', paste(U_colnames, collapse = ' + ')), df)
#'   
#'   # === linear model Y ~ A (no confounder adj)
#'   lm_YA = lm(Y ~ A, df_all)
#'   
#'   res = list(#Oracle     = coef(lm_YAU)[['A']],
#'               lmYA       = coef(lm_YA )[['A']])
#'   
#'   t0 = Sys.time()
#'   for(num_NCs in num_NC_pairs) {
#'     print(num_NCs)
#'     df = df_all[, get_CB_colnames(num_NCs)] # subset only these cols
#'     
#'     # === === === === === === === === ===  ===
#'     # === Estimators using original values ===
#'     # === === === === === === === === ===  ===
#'     # === Outcome Confounding Bridge, 2 Stage LS
#'     OCB_2SLS = OCB2SLS(df)$ATE
#'     
#'     # === Outcome Confounding Bridge, 2 Stage LS w/ ElNet Reg
#'     OCB_2SLSReg_res = OCB2SLSReg(df, alpha=ALPHA, returnMiddleDFs=FALSE)
#'     OCB_2SLSReg     = OCB_2SLSReg_res[[length(OCB_2SLSReg_res)]][['ATE']]
#'     
#'     # # Bad but,... select the latest ATE estimate that is not NA (increases by #W's)
#'     # OCB_2SLSReg_res_ATEs = c()
#'     # for(i in 1:length(OCB_2SLSReg_res)) {
#'     #   OCB_2SLSReg_res_ATEs = c(OCB_2SLSReg_res_ATEs,
#'     #                            OCB_2SLSReg_res[[i]][['ATE']])
#'     #   # OCB_2SLSReg_res_ATEsOCB_2SLSReg_res[[i]][['ATE']] |> print()
#'     # }
#'     # plot(OCB_2SLSReg_res_ATEs)
#'     # OCB_2SLSReg = OCB_2SLSReg_res_ATEs[max(which(!is.na(OCB_2SLSReg_res_ATEs)))]
#'     
#'     # === save to res = list of resulting ATEs for this num_NCs in list
#'     res_numNCs = list(OCB2SLS    = OCB_2SLS,
#'                       OCB2SLSReg = OCB_2SLSReg)
#'     
#'     
#'     # === === === === === === === === === === === === ===  ===
#'     # === Estimators with Transformations/Changes of Basis ===
#'     # === === === === === === === === === === === === ===  ===
#'     for(basis in names(basisParams)) {
#'       res_numNCs_basis = list()
#'       # print(basis)
#'       dl = constructDataListv2(df,
#'                                b_degree =basisParams[[basis]]$b_degree, 
#'                                h_degree =basisParams[[basis]]$h_degree, type = 'simple',
#'                                b_phases =basisParams[[basis]]$b_phases, 
#'                                b_periods=basisParams[[basis]]$b_periods,
#'                                h_phases =basisParams[[basis]]$h_phases, 
#'                                h_periods=basisParams[[basis]]$h_periods)
#'       
#'       
#'       # === OCBGMM
#'       res_numNCs_basis[['OCBGMM']] = OCBGMM(dl, returnAlphas = F)$ATE
#'       # res_numNCs[[paste0('OCBGMM_', basis)]] = OCBGMM(dl, returnAlphas = F)$ATE
#'       
#'       # === OCBGMM ReWeight
#'       res_numNCs_basis[['OCBGMMRW']] = OCBGMMRw(dl, GMM_steps, returnWeights = F)$ATE
#'       # res_numNCs[[paste0('OCBGMMRW_', basis)]] = OCBGMMRw(dl, GMM_steps, returnWeights = F)$ATE
#'       
#'       # === OCBGMM ReWeight Regularization
#'       res_numNCs_basis[['OCBGMMRWReg']] = 
#'                          OCBGMMRwReg(dl, GMM_steps, K_folds, lambdas, 
#'                                      returnWeights = F, returnLambdas = F)$ATE
#'       # res_numNCs[[paste0('OCBGMMRWReg_', basis)]] = 
#'       #                    OCBGMMRwReg(dl, GMM_steps, K_folds, lambdas, 
#'       #                                returnWeights = F, returnLambdas = F)$ATE
#'       
#'       # === OCB Linear CB Plug-In AND One-Step 
#'       #      Z0,Z1 NULL --> use transformations in B0,B1 as covariates 
#'       #                     when estimating nuisance params
#'       
#'       # No trimming
#'       OCBLinOS_res = tryCatch({estCondMomentOCBOS(dataList=dl, gmmSteps=GMM_steps, Z0=NULL, Z1=NULL)},
#'                               error = function(cond) {
#'                                 # message(sprintf('Error est CondMomentOCBOS with %s', 
#'                                 #                 gammaSetting)) 
#'                                 return(NULL)
#'                               }) 
#'       # OCBLinOS_res = estCondMomentOCBOS(dataList=dl, gmmSteps=GMM_steps, Z0=NULL, Z1=NULL)
#'       # errors sometimes... NA for now... debug later
#'       if(is.null(OCBLinOS_res)) {
#'         res_numNCs_basis[['OCBLinPI']] = NA
#'         res_numNCs_basis[['OCBLinOS']] = NA
#'         # res[[paste0('OCBLinPI_', basis)]] = NA
#'         # res[[paste0('OCBLinOS_', basis)]] = NA
#'       } else {
#'         res_numNCs_basis[['OCBLinPI']] = OCBLinOS_res$ATE_pi
#'         res_numNCs_basis[['OCBLinOS']] = OCBLinOS_res$ATE_os
#'         # res[[paste0('OCBLinPI_', basis)]] = OCBLinOS_res$ATE_pi
#'         # res[[paste0('OCBLinOS_', basis)]] = OCBLinOS_res$ATE_os
#'       }
#'       
#'       # With Trimming (at .1)
#'       OCBLinOS_res = tryCatch({estCondMomentOCBOS(dataList=dl, gmmSteps=GMM_steps, 
#'                                                   Z0=NULL, Z1=NULL, trim=TRIM)},
#'                               error = function(cond) {
#'                                 # message(sprintf('Error est CondMomentOCBOS with %s', 
#'                                 #                 gammaSetting)) 
#'                                 return(NULL)
#'                               }) 
#'       if(is.null(OCBLinOS_res)) {
#'         res_numNCs_basis[['OCBLinPI']] = NA
#'         # res[[paste0('OCBLinOStrim_', basis)]] = NA
#'       } else {
#'         res_numNCs_basis[['OCBLinPI']] = OCBLinOS_res$ATE_os
#'         # res[[paste0('OCBLinOStrim_', basis)]] = OCBLinOS_res$ATE_os
#'       }
#'       
#'       res_numNCs[[basis]] = res_numNCs_basis
#'     }
#'     
#'     res[[paste0('NC', num_NCs)]] = res_numNCs # save xx NC pair res to overall res
#'     
#'   }
#'   
#'   t1 = Sys.time(); print(t1 - t0)
#'   # 1 run w the 3 bases took 24.38264 mins...
#'   # saveRDS(res, 
#'   #         file = sprintf('%s/cbgenes/%s/%s/onetry.rds', 
#'   #                        save_dir, AYZW_setting_name, CB_setting_name))
#'   # 1 run w the 1-2 bases took 4 mins...
#'    
#'   
#'   res_df = data.frame(method = 'lmYA', method_type = 'naive', numNC = NA, basis = NA, ATE = res$lmYA)
#'   
#'   for(num_NC in names(res)[-1]) {
#'     # num_NC_int = substr(num_NC, start = 3, stop = nchar(num_NC))
#'     num_NC_int = substr(num_NC, start = 1, stop = nchar(num_NC)-2) |> as.integer()
#'     res_df = rbind(res_df,
#'                    data.frame(method = 'OCB2SLS', method_type = '2SLS', 
#'                               numNC = num_NC_int, basis = NA, ATE = res[[num_NC]][['OCB2SLS']] ))
#'     res_df = rbind(res_df,
#'                    data.frame(method = 'OCB2SLSReg', method_type = '2SLS', 
#'                               numNC = num_NC_int, basis = NA, ATE = res[[num_NC]][['OCB2SLSReg']] ))
#'     for(basis in names(res[[num_NC]])) {
#'       for(method in names(res[[num_NC]][[basis]])) {
#'         if(method %in% c("OCBGMM", "OCBGMMRW", "OCBGMMRWReg")) {
#'           method_type = 'GMM'
#'         } else if(method %in% c("OCBLinPI", "OCBLinOS")) {
#'           method_type = 'OCBLin'
#'         }
#'         res_df = rbind(res_df,
#'                        data.frame(method = method, method_type = method_type, 
#'                                   numNC = num_NC_int, basis = basis, ATE = res[[num_NC]][[basis]][[method]] ))
#'         
#'       }
#'       
#'     }
#'   }
#'   
#'   
#'   head(res_df); dim(res_df)
#'   naive_ate = res_df |> filter(method == 'lmYA') |> pull(ATE)
#'   # 2SLS Estimates
#'   ggplot(res_df |> filter(method_type == '2SLS'),
#'          aes(x = numNC, y = ATE, color = method, group = method)) +
#'     geom_line() +
#'     geom_hline(aes(yintercept = naive_ate), color = 'red', alpha = .7, linetype = 'dashed')
#'   
#'   # GMM Estimates
#'   ggplot(res_df |> filter(method_type == 'GMM'),
#'          aes(x = numNC, y = ATE, color = method, group = method)) +
#'     geom_line() +
#'     geom_hline(aes(yintercept = naive_ate), color = 'red', alpha = .7, linetype = 'dashed') +
#'     facet_grid(rows = vars(basis))
#'   
#'   # OCBLin Estimates
#'   ggplot(res_df |> filter(method_type == 'OCBLin'),
#'          aes(x = numNC, y = ATE, color = method, group = method)) +
#'     geom_line() +
#'     geom_hline(aes(yintercept = naive_ate), color = 'red', alpha = .7, linetype = 'dashed') +
#'     facet_grid(rows = vars(basis))
#'   
#'   res_df |> filter(method_type == 'OCBLin' & basis == 'basis3')
#'   
#' }
#' 

#' #' returns a vector of the colnames for CB when there are p NCE/NCO pairs
#' #' @return vector of chars
#' #' @example 
#' #' # p = 2 will return c('A', 'Y', 'Z1', 'Z2', 'W1', 'W2')
#' get_CB_colnames <- function(p) {
#'   c('A', 'Y', paste0('Z', 1:p), paste0('W', 1:p))
#' }


#' 
#' #' Perform CB estimate using 1:available NCE/NCO pairshttp://127.0.0.1:17257/graphics/7375ae0e-9d67-4160-a6b6-8a93b2474c44.png
#' #' @return vector of CB estimates (CB1, ..., CB(#pairs))
#' get_CB_est <- function(CBdf) {
#'   # CBdf = format_AYZW(AY_idx = AY_idx, ZW_idx = ZW_idx)
#'   
#'   num_NCENCO_pairs = (ncol(CBdf) - 2)/2
#'   CB_est = rep(NA, num_NCENCO_pairs)
#'   for(p in 1:num_NCENCO_pairs) {
#'     CB_res = CB(CBdf[, get_CB_colnames(p)]) # subset only these cols
#'     CB_est[p] = CB_res$ATE
#'   }
#'   return(CB_est)
#' }
#' 
#' 
#' #' Perform unadjusted linear regression of Y on A 
#' #' @return vector of estimate and pval
#' get_lmYA_est <- function(CBdf) {
#'   # CBdf = format_AYZW(AY_idx = AY_idx, ZW_idx = ZW_idx)
#'   lmYA = lm(Y ~ A, CBdf)
#'   lmYA$coefficients[2]
#' }
#' 
#' #' Get multiple ATE estimates
#' #' lmYA
#' #' CB1, ..., CB#pairs
#' get_ATE_est <- function(AY_idx, ZW_idx) {
#'   CBdf = format_AYZW(AY_idx = AY_idx, ZW_idx = ZW_idx)
#'   ATE_est = c(get_lmYA_est(CBdf), get_CB_est(CBdf))
#'   names(ATE_est) = c('lmYA', paste0('CB', 1:(length(ATE_est) - 1)))
#'   
#'   return(ATE_est)
#' }
#' 

# there should be #As x #(Ys per A) many AY_idx
# there should be NUM_NCENCO_per_AY number of ZW_idx

# testdf = format_AYZW(AY_idx = 3, ZW_idx = 2)
# 
# dim(testdf); head(testdf)
# dim(AY)
# 
# ATE = list() # YA_idx; ZW_idx
# NUM_NCENCO_per_AY = length(AYZW[[1]][[1]]) # prev defined, will be length of this sublist
# for(AY_idx in 1:nrow(AY)) {
#   lapply()
#   for(ZW_idx in 1:NUM_NCENCO_per_AY) {
#     # get_ATE_est(AY_idx = AY_idx, ZW_idx = ZW_idx)
#     ATE[[AY_idx]][[ZW_idx]] = get_ATE_est(AY_idx = AY_idx, ZW_idx = ZW_idx)
#   }
# }
# 
# lapply()
# 
# 
# mapply(sum, arg1=1:4, arg2=1:5)
