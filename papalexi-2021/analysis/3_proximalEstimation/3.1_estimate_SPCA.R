# ---------------------------------------------------------------------------- #
#                   Do CB using sPCA of genes as NCE/NCO                       #
# Requires: prev saved normalized gene expression (HDF5)                       #
#           prev saved specified AY pairs
#           prev saved sPCA PCs (for NCE/NCO)
# Ouputs: (nothing) but saves                                                  #
#         CBGENE_AYZW in the rds file                                          #
#                         "<save_dir>/CBGENE_AYZW.rds"                         #
# ---------------------------------------------------------------------------- #
args = commandArgs(trailingOnly = TRUE)
# args = c('laptop', 'allPerturbations')
# args = c('laptop', 'A1')

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
K_folds = 5
lambdas = (10)^(seq(from = -10, to = 3, length.out = 27))
TRIM = .1
# num_NC_pairs = c(1:8, 10, 15, 20, 25, 30) # number of NC pairs to try (<= total NC pairs, will be cut off)
num_NC_pairs = c(1, 3, 5, 10, 20)
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
source(sprintf('%s/OCB2SLS.R', util_dir))      # 2-Stage Least Squares
source(sprintf('%s/OCBGMM.R', util_dir))       # Gen. Method of Moments
source(sprintf('%s/OCBLinBridge.R', util_dir)) # Lin Bridge Plug-In and One-Step

source(sprintf('%s/getdfPCA.R', util_dir))     # for performing PCA

source(sprintf('%s/constructDataList.R', util_dir)) # for formatting data for other fns
source(sprintf('%s/CBEstAll.R', util_dir)) # for functions to estimate here
source(sprintf('%s/CBEstAllSPCA.R', util_dir)) # add'l specifically for SPCA

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










# =============================================================================

#' #' Get SE and pval of estimate when asymptotically
#' #' \sqrt{n} (betahat - beta)/\sqrt(betavar) --> N(0, 1)
#' #' @param betahat estimate of beta
#' #' @param betavar asymptotic variance of estimate
#' #' @param n sample size
#' #' @param beta true beta value or null hypothesis test for beta (typically 0)
#' get_se_pval <- function(betahat, betavar, n, beta = 0) {
#'   se = sqrt(betavar) / sqrt(n)
#'   zscore = (betahat - beta) / se
#'   pval = pnorm(-abs(zscore)) * 2
#'   return(list(se=se, pval=pval))
#' }
#' 
#' #' Make a function that takes in AY_idx (for easy looping/parallelization)
#' #' and returns ATE estimates of the specified estimators
#' #' and outer parallel function call will loop through AY rows
#' #' 
#' #' Use of an outer make function is to predefine objects used across many calls
#' #' @param AY (dataframe) of AY pair info
#' #' @param gene_norm (matrix) of gene expression (large)
#' #' @param NCs (matrix) of Negative Control values
#' #' @param grna_rownames (vector) of grna rownames from grna_odm
#' #' @param grna (matrix) of grna assignments (0 or 1) indicating if 
#' #'                      column cell received row perturbation
#' #' @param NT_idx (vector) of idx of Non-Targeting cells
#' #' @param imp_gene_names (vector) ordered vector of most important gene names
#' #' @param which_estimators (list) indicating which estimation methods to perf
#' #' @param CB_setting (list) of CB method settings
#' #' @param save_path (string) of the save path for saving intermediate ATE est
#' #'                           NULL if not to be saved
#' get_ATE_est_NCs_make <- function(AY, gene_norm, NCs,
#'                                  grna_rownames, grna, 
#'                                  NT_idx, imp_gene_names,
#'                                  which_estimators, CB_setting, save_path=NULL) {
#'   
#'   # function that returns the importance ranking of given gene name
#'   get_importance_rank = get_importance_rank_make(imp_gene_names)  
#'   
#'   
#'   
#'   
#'   # Make CB_settings accessible 
#'   ALPHA     = CB_setting$ALPHA # R does not save as pointers, so ok?
#'   GMM_steps = CB_setting$GMM_steps
#'   K_folds   = CB_setting$K_folds
#'   lambdas   = CB_setting$lambdas
#'   TRIM      = CB_setting$TRIM
#'   basisParams  = CB_setting$basisParams
#'   num_NC_pairs = CB_setting$num_NC_pairs
#'   
#'   if(!is.null(save_path)) {    dir.create(sprintf('%s/intermediateATEs/', save_path), showWarnings = FALSE)    }
#'   
#'   #' construct a dataframe of all the data 
#'   #' (to be defined inside a outer function to have access to objects)
#'   #' Requires outside named objects:
#'   #' AY, NT_idx,
#'   #' grna_rownames, grna
#'   #' gene_norm
#'   #' get_importance_rank,
#'   #' dfZ, dfW
#'   #' @param AY_idx (integer) from 1 to nrow(AY) 
#'   #' @example a = format_AYZW_inner(AY_idx = 2); head(a)
#'   format_AYZW_inner <- function(AY_idx) {
#'     AY_row = AY[AY_idx, ]
#'     # Get A
#'     # -------------------------------------------
#'     # # idx of all 'treated' cells
#'     # A_idx = which(as.logical(grna_odm[[AY_row$A, ]]))
#'     # # subset cells of 'treated' (w A grna) and 'control' (NT gran)
#'     # A = grna_odm[[AY_row$A, c(A_idx, NT_idx)]]
#'     
#'     # using grna all loaded into memory
#'     # idx of all 'treated' cells
#'     
#'     A_grna_idx = which(grna_rownames == AY_row$A)
#'     A_idx = which(as.logical(grna[A_grna_idx, 1:ncol(grna)]))
#'     # subset cells of 'treated' (w A grna) and 'control' (NT gran)
#'     A = grna[A_grna_idx, c(A_idx, NT_idx)]
#'     
#'     # Z and W 
#'     # -------------------------------------------  
#'     # already loaded in NC values outside. should not be different for every AY pair...
#'     # e.g. do this outside of loop
#'     dfZ = NCs[c(A_idx, NT_idx), seq(from = 2, to = ncol(NCs), by = 2)] # evens which((1:ncol(NCs)) %% 2 == 0)
#'     dfW = NCs[c(A_idx, NT_idx), seq(from = 1, to = ncol(NCs), by = 2)] # odds
#'     colnames(dfZ) = paste0('Z', 1:ncol(dfZ))
#'     colnames(dfW) = paste0('W', 1:ncol(dfW))
#'     
#'     # Get Y
#'     # using already loaded in gene_norm (faster)
#'     # -------------------------------------------    
#'     Y = gene_norm[get_importance_rank(AY_row$Y),     c(A_idx, NT_idx)]
#'     
#'     # Assemble together
#'     # -------------------------------------------
#'     dfAY = data.frame(A = as.vector(A), 
#'                       Y = Y)
#'     
#'     return(cbind(dfAY, dfZ, dfW))
#'   }
#'   
#'   #' Estimate ATE using variety of CB/Proximal Estimators when specifying the
#'   #' A: Exposure/Treatment
#'   #' Y: Outcome
#'   #' Negative Controls: NCE and NCOs (ZWs)
#'   get_ATE_est <- function(AY_idx) {
#'     # === Construct df ===
#'     df_all   = format_AYZW_inner(AY_idx=AY_idx)
#'     df_all$A = as.numeric(df_all$A) # convert trtmt A to numeric 0/1 if not already (alt T/F)
#'     
#'     
#'     res = data.frame(     method = character(0),
#'                           method_type = character(0),
#'                           numNC =   numeric(0),
#'                           basis = character(0),
#'                           ATE =   numeric(0),
#'                           se  =   numeric(0),
#'                           pval =   numeric(0))
#'     
#'     # === === === === === === === === ===  ===
#'     # === Oracle (none) and Naive (no adj) ===
#'     # === === === === === === === === ===  ===
#'     # # Oracle linear model Y ~ A + Us (NOT HERE IN REAL DATA? SHOULD WE TRY?)
#'     # U_colnames = grep('U', colnames(df), value = T)
#'     # lm_YAU  = lm(paste0('Y ~ A + ', paste(U_colnames, collapse = ' + ')), df)
#'     
#'     # === linear model Y ~ A (no confounder adj)
#'     
#'     if(which_estimators$lm_YA) {
#'       lm_YA = lm(Y ~ A, df_all)
#'       res = bind_rows(res, 
#'                       data.frame(
#'                         method = 'lmYA',
#'                         method_type = 'naive',
#'                         numNC = NA,
#'                         basis = NA,
#'                         ATE =    coef(lm_YA )[['A']],
#'                         se = NA,
#'                         pval= summary(lm_YA)$coefficients['A', 'Pr(>|t|)']))
#'     }
#'     
#'     
#'     # t0 = Sys.time()
#'     for(num_NCs in num_NC_pairs) {
#'       # print(num_NCs)
#'       df = df_all[, get_CB_colnames(num_NCs)] # subset only these cols
#'       
#'       
#'       # === === === === === === === === ===  ===
#'       # === Estimators using original values ===
#'       # === === === === === === === === ===  ===
#'       
#'       ## === Outcome Confounding Bridge, 2 Stage LS
#'       if(which_estimators$OCB_2SLS) {
#'         OCB_2SLS = OCB2SLS(df)$ATE
#'         
#'         res = bind_rows(res, 
#'                         data.frame(
#'                           method = 'OCB2SLS',
#'                           method_type = '2SLS',
#'                           numNC = num_NCs,
#'                           basis = NA,
#'                           ATE   = OCB_2SLS,
#'                           se    = NA,
#'                           pval  = NA))
#'       }
#'       
#'       
#'       # === Outcome Confounding Bridge, 2 Stage LS w/ pci2s package
#'       if(which_estimators$OCB_2SLS_pci2s) {
#'         OCB_2SLS_pci2s = pci2s::p2sls.lm(
#'           Y = df$Y, 
#'           A = df$A, 
#'           W = df[,grepl('W', colnames(df))], 
#'           Z = df[,grepl('Z', colnames(df))], 
#'           variance = TRUE)
#'         
#'         res = bind_rows(res, 
#'                         data.frame(
#'                           method = 'OCB2SLSpci2s',
#'                           method_type = '2SLS',
#'                           numNC = num_NCs,
#'                           basis = NA,
#'                           ATE = OCB_2SLS_pci2s$summary_second_stage['A', 'Estimate'],
#'                           se  = OCB_2SLS_pci2s$summary_second_stage['A', 'Std. Error'],
#'                           pval= OCB_2SLS_pci2s$summary_second_stage['A', 'Pr(>|z|)']))
#'       }
#'       
#'       
#'       # === Outcome Confounding Bridge, 2 Stage LS w/ ElNet Reg
#'       if(which_estimators$OCB_2SLSReg) {
#'         OCB_2SLSReg_res = OCB2SLSReg(df, alpha=ALPHA, returnMiddleDFs=FALSE)
#'         OCB_2SLSReg     = OCB_2SLSReg_res[[length(OCB_2SLSReg_res)]][['ATE']]
#'         
#'         res = bind_rows(res, 
#'                         data.frame(
#'                           method      = 'OCB2SLSReg',
#'                           method_type = '2SLS',
#'                           numNC = num_NCs,
#'                           basis = NA,
#'                           ATE   = OCB_2SLSReg,
#'                           se    = NA,
#'                           pval  = NA))
#'       }
#'       
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
#'         if(which_estimators$OCB_GMM) {
#'           OCBGMM_res = OCBGMM(dl, returnAlphas = T)
#'           
#'           GMM_coverage = checkCoverage(dl=dl, 
#'                                        alpha1 = OCBGMM_res$alpha1,
#'                                        alpha0 = OCBGMM_res$alpha0, 
#'                                        Omega  = diag(rep(1, ncol(dl$B))),
#'                                        coverage= NA, 
#'                                        ATE_est = OCBGMM_res$ATE, 
#'                                        ATE_true= 0, 
#'                                        returnVar=TRUE)
#'           
#'           OCBGMM_se_pval = get_se_pval(betahat = OCBGMM_res$ATE, 
#'                                        betavar = GMM_coverage$VarATE, 
#'                                        n = nrow(df), beta = 0)
#'           res = bind_rows(res, 
#'                           data.frame(
#'                             method      = 'OCBMEst',
#'                             method_type = 'MEst',
#'                             numNC       = num_NCs,
#'                             basis       = basis,
#'                             ATE         = OCBGMM_res$ATE,
#'                             se          = OCBGMM_se_pval$se,
#'                             pval        = OCBGMM_se_pval$pval))
#'         }
#'         
#'         
#'         
#'         
#'         
#'         # === OCBGMM ReWeight
#'         if(which_estimators$OCB_GMM) {
#'           OCBGMMRw_res = OCBGMMRw(dl, GMM_steps, returnWeights = T)
#'           
#'           GMMRw_coverage = checkCoverage(    dl=dl, 
#'                                              alpha1  = t(OCBGMMRw_res$alpha1s[GMM_steps, , drop=FALSE]),
#'                                              alpha0  = t(OCBGMMRw_res$alpha0s[GMM_steps, , drop=FALSE]), 
#'                                              Omega   = OCBGMMRw_res$Weights[[length(OCBGMMRw_res$Weights)]],
#'                                              coverage= NA, 
#'                                              ATE_est = OCBGMMRw_res$ATE, 
#'                                              ATE_true= 0, 
#'                                              returnVar=TRUE)
#'           # SE?? and pval?
#'           # OCBGMMRw_se   = sqrt(GMMRw_coverage$VarATE) / sqrt(nrow(df))
#'           # OCBGMMRw_pval = pnorm(-abs(OCBGMMRw_res$ATE - 0)/OCBGMMRw_se) * 2
#'           OCBGMMRw_se_pval = get_se_pval(betahat = OCBGMMRw_res$ATE, 
#'                                          betavar = GMMRw_coverage$VarATE, 
#'                                          n = nrow(df), beta = 0)
#'           
#'           
#'           res = bind_rows(res, 
#'                           data.frame(
#'                             method      = 'OCBMEstRW',
#'                             method_type = 'MEst',
#'                             numNC       = num_NCs,
#'                             basis       = basis,
#'                             ATE         = OCBGMMRw_res$ATE,
#'                             se          = OCBGMMRw_se_pval$se,
#'                             pval        = OCBGMMRw_se_pval$pval))
#'         }
#'         
#'         # === OCBGMM ReWeight Regularization
#'         if(which_estimators$OCB_GMM) {
#'           OCBGMMRWReg_res = OCBGMMRwReg(dl, GMM_steps, K_folds, lambdas, 
#'                                         returnWeights = F, returnLambdas = F)
#'           
#'           res = bind_rows(res, 
#'                           data.frame(
#'                             method      = 'OCBMEstRWReg',
#'                             method_type = 'MEst',
#'                             numNC       = num_NCs,
#'                             basis       = basis,
#'                             ATE         = OCBGMMRWReg_res$ATE,
#'                             se          = NA,
#'                             pval        = NA))
#'         }
#'         
#'         # NOT IMPLEMENTED: LIN OS (PI AND OS AND TRIM)
#'         # === OCB Linear CB Plug-In AND One-Step 
#'         #      Z0,Z1 NULL --> use transformations in B0,B1 as covariates 
#'         #                     when estimating nuisance params
#'         if(which_estimators$OCB_LinOS | which_estimators$OCB_LinOStrim) {
#'           # # print('Running OCB Lin OS Estimation')
#'           # # No trimming
#'           # OCBLinOS_res = tryCatch({estCondMomentOCBOS(dataList=dl, gmmSteps=GMM_steps, Z0=NULL, Z1=NULL)},
#'           #                         error = function(cond) {
#'           #                           # message(sprintf('Error est CondMomentOCBOS with %s', 
#'           #                           #                 gammaSetting)) 
#'           #                           return(NULL)
#'           #                         }) 
#'           # # OCBLinOS_res = estCondMomentOCBOS(dataList=dl, gmmSteps=GMM_steps, Z0=NULL, Z1=NULL)
#'           # # errors sometimes... NA for now... debug later
#'           # if(is.null(OCBLinOS_res)) {
#'           #   res_numNCs_basis[['OCBLinPI']] = NA
#'           #   res_numNCs_basis[['OCBLinOS']] = NA
#'           #   # res[[paste0('OCBLinPI_', basis)]] = NA
#'           #   # res[[paste0('OCBLinOS_', basis)]] = NA
#'           # } else {
#'           #   res_numNCs_basis[['OCBLinPI']] = OCBLinOS_res$ATE_pi
#'           #   res_numNCs_basis[['OCBLinOS']] = OCBLinOS_res$ATE_os
#'           #   # res[[paste0('OCBLinPI_', basis)]] = OCBLinOS_res$ATE_pi
#'           #   # res[[paste0('OCBLinOS_', basis)]] = OCBLinOS_res$ATE_os
#'           # }
#'           # 
#'           # # With Trimming (at .1)
#'           # OCBLinOS_res = tryCatch({estCondMomentOCBOS(dataList=dl, gmmSteps=GMM_steps, 
#'           #                                             Z0=NULL, Z1=NULL, trim=TRIM)},
#'           #                         error = function(cond) {
#'           #                           # message(sprintf('Error est CondMomentOCBOS with %s', 
#'           #                           #                 gammaSetting)) 
#'           #                           return(NULL)
#'           #                         }) 
#'           # if(is.null(OCBLinOS_res)) {
#'           #   res_numNCs_basis[['OCBLinOStrim']] = NA
#'           #   # res[[paste0('OCBLinOStrim_', basis)]] = NA
#'           # } else {
#'           #   res_numNCs_basis[['OCBLinOStrim']] = OCBLinOS_res$ATE_os
#'           #   # res[[paste0('OCBLinOStrim_', basis)]] = OCBLinOS_res$ATE_os
#'           # }
#'           
#'           
#'         } else if(which_estimators$OCB_LinOSPI) { # only Plug-In Estimator
#'           # # print('Running OCB Lin PI (only) Estimation')
#'           # OCBLinPI_res = tryCatch({estCondMomentOCBPI(dataList=dl, gmmSteps=GMM_steps, Z0=NULL, Z1=NULL)},
#'           #                         error = function(cond) {
#'           #                           # message(sprintf('Error est CondMomentOCBOS with %s', 
#'           #                           #                 gammaSetting)) 
#'           #                           return(NULL)
#'           #                         })
#'           # # errors sometimes... NA for now... debug later
#'           # if(is.null(OCBLinPI_res)) {
#'           #   res_numNCs_basis[['OCBLinPI']] = NA
#'           # } else {
#'           #   res_numNCs_basis[['OCBLinPI']] = OCBLinPI_res$ATE
#'           # }
#'         }
#'         
#'       }
#'       
#'       
#'       
#'     }
#'     
#'     
#'     # # 1 run: 39 sec w/  pci2s (no LinOS, all others)
#'     # # 1 run:  4 sec w/o pci2s (no LinOS, all others)
#'     # t0 = Sys.time()
#'     # get_ATE_est_NCs(AY_idx=2)
#'     # t1 = Sys.time(); print(t1 - t0)
#'     
#'     
#'     
#'     # save ATEs one by one (AY_idx, ZW_idx) if wanted
#'     if(!is.null(save_path)) {
#'       write.csv(res,
#'                 sprintf('%s/intermediateATEs/ATE_%d.csv', save_path, AY_idx),
#'                 row.names = FALSE)
#'     }
#'     
#'     return(res)
#'   }
#'   
#'   
#'   return(get_ATE_est)
#' }





  
  

  




  







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
#                                                          'yes' = sprintf('%s/spca/cbgenes/%s/%s', save_dir, AYZW_setting_name, CB_setting_name),
#                                                          'no'  = NULL),
#                                       pooled=FALSE,
#                                       pca=FALSE)
# 
# get_ATE_est0_pca = get_ATE_est_make(AY=AY, AYZW=AYZW, gene_norm=gene_norm, grna_rownames=grna_rownames, 
#                                     NT_idx=NT_idx, imp_gene_names=imp_gene_names,
#                                     CB_setting=CB_setting, run_OCBLinOSEst=run_OCBLinOSEst,
#                                     save_path = switch('no',
#                                                        'yes' = sprintf('%s/spca/cbgenes/%s/%s', save_dir, AYZW_setting_name, CB_setting_name),
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

get_ATE_est0 = get_ATE_est_NCs_make(AY=AY, gene_norm=gene_norm, NCs=NCs, grna_rownames=grna_rownames, grna=grna, 
                                    NT_idx=NT_idx, imp_gene_names=imp_gene_names,
                                    which_estimators=which_estimators, CB_setting=CB_setting, 
                                    save_path=switch(save_intermediateATEs,
                                                     'yes' = sprintf('%s/spca/cbgenes/%s/%s', save_dir, AYZW_setting_name, CB_setting_name),
                                                     'no'  = NULL))
# test = get_ATE_est0(AY_idx = 3)



get_ATE_est <- function(AY_idx) {
  res_df = tryCatch({get_ATE_est0(AY_idx=AY_idx)},
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

# NUM_NCENCO_per_AY = length(AYZW[[1]][[1]]) # prev defined, will be length of this sublist
# ATEargs = expand.grid(AY_idx = 1:nrow(AY), 
#                       ZW_idx  = 1:NUM_NCENCO_per_AY) |> 
#   arrange(AY_idx, ZW_idx) # rearrange, AY pair first

ATEargs = data.frame(AY_idx = 1:nrow(AY))
# lessen the amount...
# ATEargs = expand.grid(AY_idx = c(1:30, 120:134), ZW_idx  = 1:NUM_NCENCO_per_AY)



# NUMROWS = 10
NUMROWS = nrow(ATEargs)
# whichROWS = 300:nrow(ATEargs)
# whichROWS = 1:3
whichROWS = 1:NUMROWS
# whichROWS = 1165:NUMROWS

# # =================== Get ATEs (parallel) ====================================
print(sprintf("[%s]    - Get ATEs (parallel)", Sys.time()))
t0 = Sys.time()
ATE_par = future.apply::future_mapply(get_ATE_est,
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


saveRDS(ATE_par, sprintf('%s/spca/cbgenes/%s/%s/ATE.csv', save_dir, AYZW_setting_name, CB_setting_name))

# ATE_par[, 1]
# ATE_par[[1, ]]
# length(ATE_par)
# unlist(ATE_par)
# ATE_par[1] |> as.data.frame()
