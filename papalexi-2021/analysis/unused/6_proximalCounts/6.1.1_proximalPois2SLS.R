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
args = c('laptop', '2000', 'E_PCA')
# args = c('laptop', '2000', 'A')

suppressPackageStartupMessages(require(assertthat)) # for some assert statements
# suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(dplyr))
# suppressPackageStartupMessages(library(tibble))
# suppressPackageStartupMessages(library(ggplot2))    # plotting
# suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(ondisc))



library(future.apply)
options(future.globals.maxSize= 850*1024^2) #1st num is MB
plan(multisession, workers = 8)
# plan(sequential)

# library(furrr)
# plan(multisession, workers = 2)
# 
# theme_set(theme_cowplot() +
#             theme(plot.title = element_text(hjust = .5),
#                   plot.subtitle = element_text(hjust = .5)))

# === Parameter Settings for CB Methods ===
# ALPHA = .1
# GMM_steps = 10
# K_folds = 10
# lambdas = (10)^(seq(from = -10, to = 3, length.out = 27))
# TRIM = .1
# # num_NC_pairs = c(1:8, 10, 15, 20, 25, 30) # number of NC pairs to try (<= total NC pairs, will be cut off)
num_NCs = c(1:3, 5, 10, 15, 20, 25, 50, 100)
# num_NC_pairs = c(1:3)
# basisParams = list(basis1 = list(b_degree=1, h_degree=1, # linear
#                                  b_phases=NULL, b_periods=NULL,
#                                  h_phases=NULL, h_periods=NULL)# ,
#                    # basis2 = list(b_degree=4, h_degree=2, # larger b deg
#                    #               b_phases=NULL, b_periods=NULL,
#                    #               h_phases=NULL, h_periods=NULL)
#                    ) #,
#                    # basis3 = list(b_degree=4, h_degree=2, # add reasonable cos
#                    #               b_phases=c(0), b_periods=c(4, 8, 16),
#                    #               h_phases=c(0), h_periods=c(4, 8, 16)))
# CB_setting_name = 'simple'
# run_OCBLinOSEst = FALSE # whether to run OCBLinOne-Step estimator or not
# save_intermediateATEs = 'yes' # 'yes'/'no' whether to save intermedate ATEs as they are estimated
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
# CB_setting = list() # new copy to add ALPHA to
# CB_setting$ALPHA     = ALPHA # R does not save as pointers, so ok
# CB_setting$GMM_steps = GMM_steps
# CB_setting$K_folds   = K_folds
# CB_setting$lambdas   = lambdas
# CB_setting$TRIM      = TRIM
# CB_setting$basisParams  = basisParams
# CB_setting$num_NC_pairs = num_NC_pairs
# CB_setting$run_OCBLinOSEst = run_OCBLinOSEst
dir.create(sprintf('%s/cbgenes_count/%s', save_dir, AYZW_setting_name), showWarnings = FALSE, recursive = TRUE)
# capture.output(print(CB_setting),
#                file = sprintf('%s/cbgenes_pca/%s/%s/CB_setting.txt',
#                               save_dir, AYZW_setting_name, CB_setting_name))
# saveRDS(CB_setting,
#         sprintf('%s/cbgenes_pca/%s/%s/CB_setting.rds',
#                 save_dir, AYZW_setting_name, CB_setting_name))

# =================== Start ====================================================
print(sprintf("[%s] START: CB Estimate", Sys.time()))

# # source CB utility functions (mainly CB fn for estimating)
# source(sprintf('%s/CB_utils.R', util_dir))
# Instead source new estimation methods' files
# source(sprintf('%s/OCB2SLS.R', util_dir))      # 2-Stage Least Squares
# source(sprintf('%s/OCBGMM.R', util_dir))       # Gen. Method of Moments
# source(sprintf('%s/OCBLinBridge.R', util_dir)) # Lin Bridge Plug-In and One-Step
# 
# source(sprintf('%s/getdfPCA.R', util_dir))     # for performing PCA
# 
source(sprintf('%s/CBEstAll.R', util_dir)) # for functions gather data


source(sprintf('%s/count/OCB2SLS_count.R', util_dir)) # 2-Stage Least Squares (count)

# load chosen AYZW names
AY   = read.csv(sprintf('%s/cbgenes/%s/AY.csv', save_dir, AYZW_setting_name))
AYZW = readRDS(sprintf('%s/cbgenes/%s/AYZW.rds', save_dir, AYZW_setting_name))

# if not there: copy to pca folder too (bc plots script reads in AY, AYZW)
if(!'AY.csv' %in% list.files(sprintf('%s/cbgenes_count/%s/', save_dir, AYZW_setting_name))) {
  write.csv(AY, sprintf('%s/cbgenes_count/%s/AY.csv', save_dir, AYZW_setting_name))
}

if(!'AYZW.rds' %in% list.files(sprintf('%s/cbgenes_count/%s/', save_dir, AYZW_setting_name))) {
  saveRDS(AYZW, sprintf('%s/cbgenes_count/%s/AYZW.rds', save_dir, AYZW_setting_name))
}




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
# h5file      = paste0(save_dir, "/gene.h5"); print(h5file)
# reading_hd5file  = rhdf5::H5Fopen(name = h5file)
# readin_gene_norm = reading_hd5file&'gene_norm'
# gene_norm = readin_gene_norm[1:NUM_IMPORTANTGENES, ] # dim = 4000 x 20729 = #important x #cells
# rhdf5::h5closeAll()
# invisible(gc(verbose=FALSE))



# =================== Define File Specific Functions ======================================
print(sprintf("[%s]    - Define File Specific Functions", Sys.time()))




# All the NT grna idx
NT_names = grna_odm |> ondisc::get_feature_covariates() |>
  filter(target_type == 'non-targeting') |> rownames()
NT_idx = which(apply(X = grna_odm[[NT_names, ]], MARGIN = 2, FUN = sum) > 0)


 
# # response matrix
# # response_matrix <- gene_norm
# response_matrix <- gene_odm[[,1:ncol(gene_odm)]]
# colnames(response_matrix) <- gene_odm |> ondisc::get_cell_barcodes()
# rownames(response_matrix) <- gene_odm |> ondisc::get_feature_ids()










# ======================================
# === === === === === === === === === === === === === === === === ===

#' 
#' #' Get the gene expression count data (and grna assignment) from ondisc objects 
#' #' 
#' #' @param AY_idx (integer) index of AY pair to test
#' #' @param ZW_idx (integer) index of ZW selection 
#' #' @param AY (dataframe) of AY pair info
#' #' @param AYZW (list) of ZW info
#' #' @param gene_odm (ondisc) manager for gene expression
#' #' @param grna_odm (ondisc) manager for grna expression
#' getCountDat <- function(AY_idx, ZW_idx,
#'                         AY, AYZW,
#'                         gene_odm, grna_odm,
#'                         max_NCs = NA) {
#'   
#'   # gather the AYZW names chosen
#'   AY_row = AY[AY_idx, ]
#'   AYZW_cur = AYZW[[AY_row$A]][[AY_row$Y]][[ZW_idx]]
#'   
#'   
#'   # Get A
#'   # -------------------------------------------
#'   A_all = grna_odm[[AY_row$A,                 ]]
#'   A_idx = which(as.logical(A_all))               # idx of all 'treated' cells
#'   A     = grna_odm[[AY_row$A, c(A_idx, NT_idx)]]
#'   
#'   # Get Y and Z and W
#'   # -------------------------------------------
#'   Y = gene_odm[[AY_row$Y,  c(A_idx, NT_idx)]]
#'   
#'   
#'   if(is.na(max_NCs)) { # load all
#'     Z = t(gene_odm[[AYZW_cur$Z_names, c(A_idx, NT_idx)]]) # double [[]] pulls data
#'     W = t(gene_odm[[AYZW_cur$W_names, c(A_idx, NT_idx)]]) # into memory, slower
#'   } else { # or load up to specified max_NCs
#'     Z = t(gene_odm[[AYZW_cur$Z_names[1:max_NCs], c(A_idx, NT_idx)]]) 
#'     W = t(gene_odm[[AYZW_cur$W_names[1:max_NCs], c(A_idx, NT_idx)]]) 
#'   }
#'   
#'   colnames(Z) = paste0('Z', 1:ncol(Z))
#'   colnames(W) = paste0('W', 1:ncol(W))
#'   
#'   X = gene_odm[,c(A_idx, NT_idx)] |> ondisc::get_cell_covariates() |>
#'       dplyr::select(-n_nonzero)
#'   
#'   # put together
#'   df = cbind(A=t(A) |> as.matrix() |> as.numeric(),
#'              Y=t(Y) |> as.matrix() |> as.numeric(),
#'              Z |> as.matrix(),
#'              W |> as.matrix(),
#'              X |> as.data.frame())
#'   
#'   return(list(df = df,
#'               Z_names = colnames(Z),
#'               W_names = colnames(W),
#'               X_names = colnames(X)))
#' }

t0 = Sys.time()
# mydat = getCountDat(AY_idx = 1, ZW_idx = 2,
#                     AY = AY, AYZW = AYZW,
#                     gene_odm = gene_odm, grna_odm = grna_odm,
#                     max_NCs = max(num_NCs))
mydat = getCountDat(AY_idx = 285, ZW_idx = 2, # pass in inner args
                   # pass in outer args
                   AY = AY, AYZW = AYZW, NT_idx = NT_idx,
                   gene_odm = gene_odm, grna_odm = grna_odm, 
                   max_NCs = max(num_NCs))
print(Sys.time() - t0)

Y = mydat$df$Y
A = mydat$df$A

W = mydat$df[, mydat$W_names[1:5]] # as counts
Z = mydat$df[, mydat$Z_names[1:5]] # as counts (perhaps these should be normalized?)

t0 = Sys.time()
p2sls_result <- pci2s::p2sls.negbin(Y = Y, A = A, X = NULL,
                                    W = W, Z = Z,
                                    # nco_type = rep('negbin', ncol(W)) 
                                    nco_type = rep('poisson', ncol(W)) 
                                    )
                                    # no offset for now
                                    # nco_args = list(list(offset = rep(0, N)),
                                    #                 list(offset = rep(0, N))))
t1 = Sys.time(); print(t1 - t0)


p2sls_result$summary_first_stage
p2sls_result$summary_second_stage




# mydat$df[,mydat$Z_names] |> dim()




####################### Estimating ATEs ######################################


#' wrapper for pci2s using specific number'
#' @param nco_type (character) 'poisson' or 'negbin', or a vector of these
#' @param num_NCs
runpci2s <- function(Y, A, X, W, Z, nco_type, num_NCs=NULL) {
  
  # nco type input as single char or vec of chars
  if(length(nco_type) == 1) {
    nco_types = rep(nco_type, ncol(W))
  } else {
    nco_type = nco_type
  }
  
  if(is.null(num_NCs)) {
    # use all
    num_NCs = ncol(W)
  } 
  
  paramEstimates = c() # of the parameter coef of A
  se             = c() # of the parameter coef of A
  pval           = c() # of the parameter coef of A
  for(num_NC in num_NCs) {
    p2sls_result <- pci2s::p2sls.negbin(Y = Y, 
                                        A = A, 
                                        X = X,
                                        W = W[, 1:num_NC], 
                                        # Z = Z[, 1:num_NC],
                                        Z = Z,
                                        nco_type = nco_types[1:num_NC])
    paramEstimates = c(paramEstimates, p2sls_result$summary_second_stage['A', 'Estimate'])
    se             = c(se,             p2sls_result$summary_second_stage['A', 'Std. Error'])
    pval           = c(pval,           p2sls_result$summary_second_stage['A', 'Pr(>|z|)'])
  }
  
  res = data.frame(idx      = 1:length(num_NCs),
                   numNCO   = num_NC,
                   Estimate = paramEstimates,
                   se       = se,
                   pval     = pval)
  
  return(res)    
}


#' make function for EstATECount (oracle, naive, prox estimators for count)
#' @param AY_idx (integer) index of AY pair to test
#' @param ZW_idx (integer) index of ZW selection
#' @param AY (dataframe) of AY pair info
#' @param AYZW (list) of ZW info
#' @param NT_idx (vector) of idx of Non-Targeting cells
#' @param gene_odm (ondisc) manager for gene expression
#' @param grna_odm (ondisc) manager for grna expression#' @param include_pci2sPois (boolean) should pci2s Pois method be performed
#' @param pci2sPois_maxNC (integer) maximum number of NCs to use for pci2s Pois (to decrease time)
#' @param include_pci2sNB (boolean) should pci2s Neg Bin method be performed
#' @param pci2sNB_maxNC (integer) maximum number of NCs to use for pci2s NB (to decrease time)
makeEstATECount <- function(AY, AYZW, NT_idx,
                            gene_odm, grna_odm,
                            num_NCs,
                            max_NCs = NA,
                            include_pci2sPois=FALSE, pci2sPois_maxNC=NULL,
                            include_pci2sNB  =FALSE, pci2sNB_maxNC=NULL) {

  #' returns ATE estimates from different estimation methods
  EstATECount <- function(AY_idx, ZW_idx) {
    mydat = getCountDat(AY_idx = AY_idx, ZW_idx = ZW_idx, # pass in inner args
                        # pass in outer args
                        AY = AY, AYZW = AYZW, NT_idx = NT_idx,
                        gene_odm = gene_odm, grna_odm = grna_odm,
                        max_NCs = max_NCs)




    # Oracle glm Y ~ A + Us
    glm_YAU  = glm(paste0('Y ~ A + ', paste(mydat$X_names, collapse = ' + ')),
                   data = mydat$df, family = poisson())

    # glm Y ~ A (no confounder adj)
    glm_YA   = glm('Y ~ A',
                   data = mydat$df, family = poisson())

    # X_design = model.matrix(glm_YAU)
    # coefs = glm_YAU$coefficients; coefs[is.na(coefs)] = 0
    # (X_design %*% coefs)[1:5] |> exp()
    # predict.glm(glm_YAU, newdata = mydat$df, type = 'response')[1:5]
    
    
    glm_res = matrix(
      c(NA,
        NA,
        predict.glm(glm_YAU, newdata = mydat$df |> mutate(A=1), type = 'response') |> mean(),
        predict.glm(glm_YAU, newdata = mydat$df |> mutate(A=0), type = 'response') |> mean(),
        NA,
        NA,
        predict.glm(glm_YA , newdata = mydat$df |> mutate(A=1), type = 'response') |> mean(),
        predict.glm(glm_YA , newdata = mydat$df |> mutate(A=0), type = 'response') |> mean()),
      nrow=2, byrow=TRUE) |> data.frame()
    glm_res[, 1] = c('Oracle', 'glmYA')
    colnames(glm_res) = c('method', 'numNCO', 'EY1', 'EY0')


    modelW = makeModelPois(logTransX = FALSE)
    modelY = makeModelPois(logTransX = TRUE)

    myPoisOCB2SLS = makeOCB2SLS(modelW = modelW, modelY = modelY)

    OCB2SLS_res = myPoisOCB2SLS(Y = mydat$df$Y,
                                A = mydat$df$A,
                                W = mydat$df[, mydat$W_names],
                                Z = mydat$df[, mydat$Z_names],
                                returnIntermediateATEs = num_NCs)

    OCB2SLS_res$intermediateATEs |> mutate(method = 'Pois2SLS') |> select(method, numNCO, EY1, EY0)

    res = rbind(glm_res,
                OCB2SLS_res$intermediateATEs |>
                  mutate(method = 'Pois2SLS') |>
                  select(method, numNCO, EY1, EY0))
    
    
    if(include_pci2sPois) {
      # TODO: 
      # perform estimation
      if(!is.null(pci2sPois_maxNC)) {
        pci2sPois_num_NCs = num_NCs(  num_NCs <= pci2sPois_maxNC)
      } else {
        pci2sPois_num_NCs = num_NCs
      }
      
      runpci2s(Y = mydat$df$Y, 
               A = mydat$df$A, 
               X = NULL,
               W = mydat$df[, mydat$W_names], 
               Z = mydat$df[, mydat$Z_names],
               nco_type = 'poisson', 
               num_NCs = c(1, 2))
      
      p2sls_result <- pci2s::p2sls.negbin(Y = mydat$df$Y, 
                                          A = mydat$df$A, 
                                          X = NULL,
                                          W = mydat$df[, mydat$W_names], 
                                          Z = mydat$df[, mydat$Z_names],
                                          nco_type = rep('poisson', ncol(W)) 
      )
      
      if(returnIntermediateATEs) {
        
      }
      # add results
      
    }
    if(include_pci2sNB) {
      # TODO:  nco_type = rep('negbin', ncol(W)) 
      # perform estimation
      p2sls_result <- pci2s::p2sls.negbin(Y = mydat$df$Y, 
                                          A = mydat$df$A, 
                                          X = NULL,
                                          W = mydat$df[, mydat$W_names[1:3]], 
                                          Z = mydat$df[, mydat$Z_names[1:3]],
                                          nco_type = rep('negbin', length(mydat$W_names[1:3]))) 
      p2slsPois_result <- pci2s::p2sls.negbin(Y = mydat$df$Y, 
                                          A = mydat$df$A, 
                                          X = NULL,
                                          W = mydat$df[, mydat$W_names[1:3]], 
                                          Z = mydat$df[, mydat$Z_names[1:3]],
                                          nco_type = rep('poisson', length(mydat$W_names[1:3]))) 
      if(length(nco_type) == 1) {
        nco_types = rep(nco_type, ncol(W))
      } else {
        nco_type = nco_type
      }
      # add results
      paramEstimates = c() # of the parameter coef of A
      se             = c() # of the parameter coef of A
      pval           = c() # of the parameter coef of A
      for(num_NC in num_NCs) {
        p2sls_result <- pci2s::p2sls.negbin(Y = Y, 
                                            A = A, 
                                            X = X,
                                            W = W[, 1:num_NC], 
                                            # Z = Z[, 1:num_NC],
                                            Z = Z,
                                            nco_type = nco_types[1:num_NC])
        paramEstimates = c(paramEstimates, p2sls_result$summary_second_stage['A', 'Estimate'])
        se             = c(se,             p2sls_result$summary_second_stage['A', 'Std. Error'])
        pval           = c(pval,           p2sls_result$summary_second_stage['A', 'Pr(>|z|)'])
      }
      
      res = data.frame(idx      = 1:length(num_NCs),
                       numNCO   = num_NC,
                       Estimate = paramEstimates,
                       se       = se,
                       pval     = pval)
    }
    
    # calculate ATE and Fold Change
    res$ATE = res$EY1 - res$EY0
    res$FC  = res$EY1 / res$EY0

    return(res)

  }


}








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
#                                                          'yes' = sprintf('%s/cbgenes_pca/%s/%s', save_dir, AYZW_setting_name, CB_setting_name),
#                                                          'no'  = NULL),
#                                       pooled=FALSE,
#                                       pca=FALSE)
# 
# get_ATE_est0_pca = get_ATE_est_make(AY=AY, AYZW=AYZW, gene_norm=gene_norm, grna_rownames=grna_rownames, 
#                                     NT_idx=NT_idx, imp_gene_names=imp_gene_names,
#                                     CB_setting=CB_setting, run_OCBLinOSEst=run_OCBLinOSEst,
#                                     save_path = switch('no',
#                                                        'yes' = sprintf('%s/cbgenes_pca/%s/%s', save_dir, AYZW_setting_name, CB_setting_name),
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


myEstATECount0 = makeEstATECount(AY = AY, AYZW = AYZW, NT_idx = NT_idx,
                                gene_odm = gene_odm, grna_odm = grna_odm, 
                                num_NCs = num_NCs,
                                max_NCs = max(num_NCs))

# myEstATECount0(1, 1)
# test = getCountDat(AY_idx = 1, ZW_idx = 1,
#                    AY=AY, AYZW=AYZW, NT_idx=NT_idx,
#                    gene_odm = gene_odm, grna_odm, max_NCs = 5)
# test$df |> head()
# glm(offset = )

# get_ATE_est0 = get_ATE_est_make(AY=AY, AYZW=AYZW, gene_norm=gene_norm, grna_rownames=grna_rownames, grna=grna, 
#                                 NT_idx=NT_idx, imp_gene_names=imp_gene_names,
#                                 CB_setting=CB_setting, run_OCBLinOSEst=run_OCBLinOSEst,
#                                 save_path = switch(save_intermediateATEs,
#                                                    'yes' = sprintf('%s/cbgenes_pca/%s/%s', save_dir, AYZW_setting_name, CB_setting_name),
#                                                    'no'  = NULL),
#                                 pooled=FALSE,
#                                 pca=TRUE)


myEstATECount <- function(AY_idx, ZW_idx) {
  res_df = tryCatch({myEstATECount0(AY_idx=AY_idx, ZW_idx=ZW_idx)},
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
# whichROWS = 1:32
# whichROWS = 1:NUMROWS
whichROWS = 1:nrow(ATEargs)

# # testing w sequential
# ATE_seq = list()
# for(r in whichROWS) {
#   print(r)
#   ATE_seq[[r]] = myEstATECount0(AY_idx = ATEargs[r, 1],
#                      ZW_idx = ATEargs[r, 2])
# }


# =================== Get ATEs (parallel) ====================================
print(sprintf("[%s]    - Get ATEs (parallel)", Sys.time()))
t0 = Sys.time()
ATE_par = future.apply::future_mapply(myEstATECount,
                                      AY_idx = ATEargs[whichROWS, 1], # ATEargs[1:NUMROWS, 1],
                                      ZW_idx = ATEargs[whichROWS, 2], # ATEargs[1:NUMROWS, 2],
                                      future.globals = TRUE,
                                      future.seed = 56789, 
                                      SIMPLIFY = FALSE)
t1 = Sys.time()

print(sprintf("[%s]        - assembling together, parallel took %2.2f", Sys.time(), (t1 - t0)))
# assemble together
ATEs = NULL
for(i in 1:length(ATE_par)) {
  if(!is.null(ATE_par[[i]])) {
    ATEs = rbind(ATEs,
                 cbind(data.frame(AY_idx = ATEargs[i, 1],
                                  ZW_idx = ATEargs[i, 2]), 
                       ATE_par[[i]]))
  } else {
    print(sprintf('%d (AY_idx: %d, ZW_idx) %d was NULL', i, ATEargs[i, 1], ATEargs[i, 1]))
  }
}


# print(sprintf("[%s]        - %2.2f", Sys.time(), (t1 - t0)))


# =================== Get ATEs (parallel) ====================================
print(sprintf("[%s]    - Saving ATEs", Sys.time()))

write.csv(ATEs, sprintf('%s/cbgenes_count/%s/ATE.csv', save_dir, AYZW_setting_name), row.names = FALSE)


# =================== Start ====================================================
print(sprintf("[%s] END", Sys.time()))


