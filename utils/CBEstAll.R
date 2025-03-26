# Helper functions for fitting CB estimates on papalexi (and gasperini) data
# requires loading in CB estimation methods in other util files such as OCB2SLS, OCBGMM, OCBLinBridge



#' important_genes_name.rds should be in order of most to least important
#' the rank will be the index of the gene on this list
#' Maybe dont need this...
#' Make function that gets the importance ranking of specified gene_name using
#' the imput imp_gene_names parameter. The imp_gene_names should be ordered by 
#' importance, from most important to least important.
#' @param imp_gene_names (vector) of strings of gene names ordered by importance
#' @param gene_names (string) gene name
#' @example 
#' # suppose you want a function that returns the importance ranking according 
#' to myImpGeneNames 
#' myGetImpRank = get_importance_rank_make(myImpGeneNames)
#' myGetImpRank(myGeneName) # get the imp ranking of myGeneName
get_importance_rank_make <- function(imp_gene_names) {
  get_importance_rank <- function(gene_name) {
    which(imp_gene_names == gene_name)
  }
  return(get_importance_rank)
}



#' prep data/formatting
#' creates a dataframe with
#' A (boolean) of getting specified gRNA or any NT perturbation
#' Y, Z1, ..., W1, ... (numeric) normalized gene expressions 
#' @param get_importance_rank (function) that takes in arg gene name (str) and 
#'            returns the importance ranking (which is used as the idx to get 
#'            values from gene_norm)
#' @param pooled (boolean) pooled setting? allow for pooled input version for ZWs
#' @param pca (boolean) pca setting
#' @param pca_rank (integer) rank of PCA 
#'        (automatic input from get_ATE_est_make should be max(num_NC_pairs))
#' @param cca (boolean) cca setting
#' @param cca_rank (integer) rank of CCA 
#'        (automatic input from get_ATE_est_make should be max(num_NC_pairs))
format_AYZW <- function(AY_idx, ZW_idx, 
                        AY, AYZW, 
                        gene_norm,
                        grna_rownames, #imp_gene_names
                        grna,
                        NT_idx,
                        get_importance_rank,
                        pooled=FALSE,
                        pca=FALSE,
                        pca_rank=5,
                        cca=FALSE,
                        cca_rank=5) {
  # get_importance_rank = get_importance_rank_make(imp_gene_names)   
  # AY_idx = 3; ZW_idx = 2
  AY_row = AY[AY_idx, ]
  
  
  # gather the AYZW names chosen
  AYZW_cur = AYZW[[AY_row$A]][[AY_row$Y]][[ZW_idx]]
  
  
  # AY_row$A
  # AY_row$Y
  # AYZW_cur$Z_names
  # AYZW_cur$W_names
  
  # Get A
  # -------------------------------------------
  # # idx of all 'treated' cells
  # A_idx = which(as.logical(grna_odm[[AY_row$A, ]]))
  # # subset cells of 'treated' (w A grna) and 'control' (NT gran)
  # A = grna_odm[[AY_row$A, c(A_idx, NT_idx)]]

  # using grna all loaded into memory
  # idx of all 'treated' cells
  
  A_grna_idx = which(grna_rownames == AY_row$A)
  A_idx = which(as.logical(grna[A_grna_idx, 1:ncol(grna)]))
  # subset cells of 'treated' (w A grna) and 'control' (NT gran)
  A = grna[A_grna_idx, c(A_idx, NT_idx)]
    
  
  # # Get Y and Z and W
  # # using h5 load in (use if not enough memory)
  # # -------------------------------------------
  # h5file      = paste0(save_dir, "/gene.h5")
  # reading_hd5file  = rhdf5::H5Fopen(name = h5file)
  # readin_gene_norm = reading_hd5file&'gene_norm'
  # Y = readin_gene_norm[get_importance_rank(AY_row$Y),     c(A_idx, NT_idx)]
  # 
  # Z = readin_gene_norm[sapply(A = AYZW_cur$Z_names,
  #                             FUN = get_importance_rank), c(A_idx, NT_idx)]
  # 
  # W = readin_gene_norm[sapply(A = AYZW_cur$W_names,
  #                             FUN = get_importance_rank), c(A_idx, NT_idx)]
  # rhdf5::h5closeAll()
  # invisible(gc(verbose=FALSE))
  
  # Get Y and Z and W 
  # using already loaded in gene_norm (faster)
  # -------------------------------------------    
  Y = gene_norm[get_importance_rank(AY_row$Y),     c(A_idx, NT_idx)]

  if(!pooled) { # individual genes' expressions (requires different AYZW format)
    Z = gene_norm[sapply(X = AYZW_cur$Z_names, 
                       FUN = get_importance_rank), c(A_idx, NT_idx)]
    W = gene_norm[sapply(X = AYZW_cur$W_names, 
                       FUN = get_importance_rank), c(A_idx, NT_idx)]
    dfZ = data.frame(t(Z))
    dfW = data.frame(t(W))
  } else { # pool together genes (requires different AYZW format)
    mapZ = 
      # dataframe of Z gene names, chrs, NCidx
      data.frame(Y = AYZW_cur$Z_names, chr = AYZW_cur$Z_chrs, NCidx = AYZW_cur$Z_NCidx) |> 
      # list of pooled avg gene expr
      group_by(chr, NCidx) |> 
      group_map(~ {data.frame(pooled = gene_norm[sapply(X = .x$Y, 
                                                  FUN = get_importance_rank), 
                                           c(A_idx, NT_idx)] 
                        |> colMeans())                             })
    # combine pooled to wide df format
    dfZ = dplyr::bind_cols(lapply(mapZ, as.data.frame.list), .name_repair = 'unique_quiet')
    
    mapW = 
      # dataframe of W gene names, chrs, NCidx
      data.frame(Y = AYZW_cur$W_names, chr = AYZW_cur$W_chrs, NCidx = AYZW_cur$W_NCidx) |> 
      # list of pooled avg gene expr
      group_by(chr, NCidx) |> 
      group_map(~ {data.frame(pooled = gene_norm[sapply(X = .x$Y, 
                                                  FUN = get_importance_rank), 
                                           c(A_idx, NT_idx)] 
                        |> colMeans())                             })
    # combine pooled to wide df format
    dfW = dplyr::bind_cols(lapply(mapW, as.data.frame.list), .name_repair = 'unique_quiet')
  }

  # PCA on already gathered Z,Ws
  if(pca) {
    dfZ = getPCArotations(dfZ, pca_rank) # requires getdfPCA.R
    dfW = getPCArotations(dfW, pca_rank)
  }

  # CCA on already gathered Z,Ws
  if(cca) {
    # print('performing cca')
    cca_res = cancor(x = dfZ, y = dfW,
                        xcenter=TRUE, ycenter=TRUE) # base R fn
    dfZ = as.matrix(dfZ) %*% cca_res$xcoef[,1:min(cca_rank, ncol(cca_res$xcoef)), drop=FALSE]
    dfW = as.matrix(dfW) %*% cca_res$ycoef[,1:min(cca_rank, ncol(cca_res$ycoef)), drop=FALSE]
  }

  # Assemble together
  # -------------------------------------------
  dfAY = data.frame(A = as.vector(A), 
                    Y = Y)
  colnames(dfZ) = paste0('Z', 1:ncol(dfZ))
  colnames(dfW) = paste0('W', 1:ncol(dfW))
  
  
  return(cbind(dfAY, dfZ, dfW))
  
}





#' makes a function that takes in AY_idx, ZW_idx to calculate ATE from a variety
#' of estimation methods, using the input 
#' @param CB_setting (list) of parameter setttings needed for CB e.g. ALPHA, GMM_steps, ...
#' @param run_OCBLinOSEst (boolean) whether or not to run the OCB Lin One-Step Estimator (performs badly)
#' @param est_lmYA_only (boolean) only estimate the naive lmYA + return pvalues
#' @param pooled (boolean) pooled setting?
#' @param pca (boolean) pca setting?
get_ATE_est_make <- function(AY, AYZW, gene_norm, grna_rownames, grna, NT_idx, imp_gene_names,
                             CB_setting, run_OCBLinOSEst, save_path=NULL,
                             est_lmYA_only=FALSE,
                             pooled=FALSE,
                             pca=FALSE,
                             cca=FALSE) {
  get_importance_rank = get_importance_rank_make(imp_gene_names)  
  
  if(!est_lmYA_only) { # don't need this input if not doing CB estimates
    ALPHA     = CB_setting$ALPHA # R does not save as pointers, so ok?
    GMM_steps = CB_setting$GMM_steps
    K_folds   = CB_setting$K_folds
    lambdas   = CB_setting$lambdas
    TRIM      = CB_setting$TRIM
    basisParams  = CB_setting$basisParams
    num_NC_pairs = CB_setting$num_NC_pairs
  }
  
  
  if(!is.null(save_path)) {    dir.create(sprintf('%s/intermediateATEs/', save_path))    }
  
  #' returns ATE estimates from different estimation methods
  get_ATE_est <- function(AY_idx, ZW_idx) {
    # AY_idx = 3; ZW_idx = 2
    # === Construct df ===
    df_all   = format_AYZW(AY_idx=AY_idx, ZW_idx=ZW_idx, # inputs to this inner fn
                           AY=AY, AYZW=AYZW,             # inputs saved from outer fn
                           gene_norm=gene_norm,
                           grna_rownames=grna_rownames,
                           grna=grna,
                           NT_idx=NT_idx,
                           get_importance_rank=get_importance_rank,
                           pooled=pooled,
                           pca=pca, 
                           pca_rank=max(num_NC_pairs),
                           cca=cca,
                           cca_rank=max(num_NC_pairs))
    df_all$A = as.numeric(df_all$A) # convert trtmt A to numeric 0/1 if not already (alt T/F)
    
    # === Get different estimates === 
    
    
    # === === === === === === === === ===  ===
    # === Oracle (none) and Naive (no adj) ===
    # === === === === === === === === ===  ===
    # # Oracle linear model Y ~ A + Us (NOT HERE IN REAL DATA? SHOULD WE TRY?)
    # U_colnames = grep('U', colnames(df), value = T)
    # lm_YAU  = lm(paste0('Y ~ A + ', paste(U_colnames, collapse = ' + ')), df)
    
    # === linear model Y ~ A (no confounder adj)
    lm_YA = lm(Y ~ A, df_all)
    if(est_lmYA_only) {
        return(data.frame(
                   ATE =    coef(lm_YA )[['A']],
                   pval= summary(lm_YA)$coefficients['A', 'Pr(>|t|)']))
    }
    
    res = list(#Oracle     = coef(lm_YAU)[['A']],
                lmYA       = coef(lm_YA )[['A']])
    
    # t0 = Sys.time()
    for(num_NCs in num_NC_pairs) {
      # print(num_NCs)
      df = df_all[, get_CB_colnames(num_NCs)] # subset only these cols
      
      # === === === === === === === === ===  ===
      # === Estimators using original values ===
      # === === === === === === === === ===  ===
      # === Outcome Confounding Bridge, 2 Stage LS
      OCB_2SLS = OCB2SLS(df)$ATE
      
      # === Outcome Confounding Bridge, 2 Stage LS w/ ElNet Reg
      OCB_2SLSReg_res = OCB2SLSReg(df, alpha=ALPHA, returnMiddleDFs=FALSE)
      OCB_2SLSReg     = OCB_2SLSReg_res[[length(OCB_2SLSReg_res)]][['ATE']]
      
      # # Bad but,... select the latest ATE estimate that is not NA (increases by #W's)
      # OCB_2SLSReg_res_ATEs = c()
      # for(i in 1:length(OCB_2SLSReg_res)) {
      #   OCB_2SLSReg_res_ATEs = c(OCB_2SLSReg_res_ATEs,
      #                            OCB_2SLSReg_res[[i]][['ATE']])
      #   # OCB_2SLSReg_res_ATEsOCB_2SLSReg_res[[i]][['ATE']] |> print()
      # }
      # plot(OCB_2SLSReg_res_ATEs)
      # OCB_2SLSReg = OCB_2SLSReg_res_ATEs[max(which(!is.na(OCB_2SLSReg_res_ATEs)))]
      
      # === save to res = list of resulting ATEs for this num_NCs in list
      res_numNCs = list(OCB2SLS    = OCB_2SLS,
                        OCB2SLSReg = OCB_2SLSReg)
      
      
      # === === === === === === === === === === === === ===  ===
      # === Estimators with Transformations/Changes of Basis ===
      # === === === === === === === === === === === === ===  ===
      for(basis in names(basisParams)) {
        res_numNCs_basis = list()
        # print(basis)
        dl = constructDataListv2(df,
                                 b_degree =basisParams[[basis]]$b_degree, 
                                 h_degree =basisParams[[basis]]$h_degree, type = 'simple',
                                 b_phases =basisParams[[basis]]$b_phases, 
                                 b_periods=basisParams[[basis]]$b_periods,
                                 h_phases =basisParams[[basis]]$h_phases, 
                                 h_periods=basisParams[[basis]]$h_periods)
        
        
        # === OCBGMM
        res_numNCs_basis[['OCBGMM']] = OCBGMM(dl, returnAlphas = F)$ATE
        # res_numNCs[[paste0('OCBGMM_', basis)]] = OCBGMM(dl, returnAlphas = F)$ATE
        
        # === OCBGMM ReWeight
        res_numNCs_basis[['OCBGMMRW']] = OCBGMMRw(dl, GMM_steps, returnWeights = F)$ATE
        # res_numNCs[[paste0('OCBGMMRW_', basis)]] = OCBGMMRw(dl, GMM_steps, returnWeights = F)$ATE
        
        # === OCBGMM ReWeight Regularization
        res_numNCs_basis[['OCBGMMRWReg']] = 
                           OCBGMMRwReg(dl, GMM_steps, K_folds, lambdas, 
                                       returnWeights = F, returnLambdas = F)$ATE
        # res_numNCs[[paste0('OCBGMMRWReg_', basis)]] = 
        #                    OCBGMMRwReg(dl, GMM_steps, K_folds, lambdas, 
        #                                returnWeights = F, returnLambdas = F)$ATE
        
        # === OCB Linear CB Plug-In AND One-Step 
        #      Z0,Z1 NULL --> use transformations in B0,B1 as covariates 
        #                     when estimating nuisance params
        if(run_OCBLinOSEst) {
            # print('Running OCB Lin OS Estimation')
            # No trimming
            OCBLinOS_res = tryCatch({estCondMomentOCBOS(dataList=dl, gmmSteps=GMM_steps, Z0=NULL, Z1=NULL)},
                                    error = function(cond) {
                                      # message(sprintf('Error est CondMomentOCBOS with %s', 
                                      #                 gammaSetting)) 
                                      return(NULL)
                                    }) 
            # OCBLinOS_res = estCondMomentOCBOS(dataList=dl, gmmSteps=GMM_steps, Z0=NULL, Z1=NULL)
            # errors sometimes... NA for now... debug later
            if(is.null(OCBLinOS_res)) {
              res_numNCs_basis[['OCBLinPI']] = NA
              res_numNCs_basis[['OCBLinOS']] = NA
              # res[[paste0('OCBLinPI_', basis)]] = NA
              # res[[paste0('OCBLinOS_', basis)]] = NA
            } else {
              res_numNCs_basis[['OCBLinPI']] = OCBLinOS_res$ATE_pi
              res_numNCs_basis[['OCBLinOS']] = OCBLinOS_res$ATE_os
              # res[[paste0('OCBLinPI_', basis)]] = OCBLinOS_res$ATE_pi
              # res[[paste0('OCBLinOS_', basis)]] = OCBLinOS_res$ATE_os
            }
            
            # With Trimming (at .1)
            OCBLinOS_res = tryCatch({estCondMomentOCBOS(dataList=dl, gmmSteps=GMM_steps, 
                                                        Z0=NULL, Z1=NULL, trim=TRIM)},
                                    error = function(cond) {
                                      # message(sprintf('Error est CondMomentOCBOS with %s', 
                                      #                 gammaSetting)) 
                                      return(NULL)
                                    }) 
            if(is.null(OCBLinOS_res)) {
              res_numNCs_basis[['OCBLinOStrim']] = NA
              # res[[paste0('OCBLinOStrim_', basis)]] = NA
            } else {
              res_numNCs_basis[['OCBLinOStrim']] = OCBLinOS_res$ATE_os
              # res[[paste0('OCBLinOStrim_', basis)]] = OCBLinOS_res$ATE_os
            }
            
            
          } else { # only Plug-In Estimator
            # print('Running OCB Lin PI (only) Estimation')
            OCBLinPI_res = tryCatch({estCondMomentOCBPI(dataList=dl, gmmSteps=GMM_steps, Z0=NULL, Z1=NULL)},
                                    error = function(cond) {
                                      # message(sprintf('Error est CondMomentOCBOS with %s', 
                                      #                 gammaSetting)) 
                                      return(NULL)
                                    })
            # errors sometimes... NA for now... debug later
            if(is.null(OCBLinPI_res)) {
              res_numNCs_basis[['OCBLinPI']] = NA
            } else {
              res_numNCs_basis[['OCBLinPI']] = OCBLinPI_res$ATE
            }
        }

        res_numNCs[[basis]] = res_numNCs_basis
      }
        
      
      res[[paste0('NC', num_NCs)]] = res_numNCs # save xx NC pair res to overall res
      
    }
    
    # t1 = Sys.time(); print(t1 - t0)
    # 1 run took 24.38264 mins...
    # saveRDS(res, 
    #         file = sprintf('%s/cbgenes/%s/%s/onetry.rds', 
    #                        save_dir, AYZW_setting_name, CB_setting_name))
    
    
    # maybe shouldve just saved as a df...
    res_df = data.frame(method = 'lmYA', method_type = 'naive', numNC = NA, basis = NA, ATE = res$lmYA)
    
    for(num_NC in names(res)[-1]) {
      num_NC_int = substr(num_NC, start = 3, stop = nchar(num_NC)) |> as.integer()
      # num_NC_int = substr(num_NC, start = 1, stop = nchar(num_NC)-2) |> as.integer()
      res_df = rbind(res_df,
                     data.frame(method = 'OCB2SLS', method_type = '2SLS', 
                                numNC = num_NC_int, basis = NA, ATE = res[[num_NC]][['OCB2SLS']] ))
      res_df = rbind(res_df,
                     data.frame(method = 'OCB2SLSReg', method_type = '2SLS', 
                                numNC = num_NC_int, basis = NA, ATE = res[[num_NC]][['OCB2SLSReg']] ))
      for(basis in names(res[[num_NC]])) {
        for(method in names(res[[num_NC]][[basis]])) {
          if(method %in% c("OCBGMM", "OCBGMMRW", "OCBGMMRWReg")) {
            method_type = 'GMM'
          } else if(method %in% c("OCBLinPI", "OCBLinOS", "OCBLinOStrim")) {
            method_type = 'OCBLin'
          }
          res_df = rbind(res_df,
                         data.frame(method = method, method_type = method_type, 
                                    numNC = num_NC_int, basis = basis, ATE = res[[num_NC]][[basis]][[method]] ))
          
        }
        
      }
    }
    
    # save ATEs one by one (AY_idx, ZW_idx) if wanted
    if(!is.null(save_path)) {
    write.csv(res_df,
              sprintf('%s/intermediateATEs/ATE_%d_%d.csv', save_path, AY_idx, ZW_idx),
              row.names = FALSE)
    }
    
    return(res_df)
  }
  
  

  return(get_ATE_est)
}



#' returns a vector of the colnames for CB when there are p NCE/NCO pairs
#' @return vector of chars
#' @example 
#' # p = 2 will return c('A', 'Y', 'Z1', 'Z2', 'W1', 'W2')
get_CB_colnames <- function(p) {
  c('A', 'Y', paste0('Z', 1:p), paste0('W', 1:p))
}



#' make function for EstATECount (oracle, naive, prox estimators for count) 
#' @param AY_idx (integer) index of AY pair to test
#' @param ZW_idx (integer) index of ZW selection 
#' @param AY (dataframe) of AY pair info
#' @param AYZW (list) of ZW info
#' @param NT_idx (vector) of idx of Non-Targeting cells
#' @param gene_odm (ondisc) manager for gene expression
#' @param grna_odm (ondisc) manager for grna expression
#' @param include_pci2sPois (boolean) should pci2s Pois method be performed
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
    
    
    glm_res = matrix(
      c(# Oracle glmYA
        NA, 
        NA,
        predict.glm(glm_YAU, newdata = mydat$df |> mutate(A=1), type = 'response') |> mean(),
        predict.glm(glm_YAU, newdata = mydat$df |> mutate(A=0), type = 'response') |> mean(),
        summary(glm_YAU)$coefficients['A', 'Pr(>|z|)'],
        # glmYA
        NA,
        NA,
        predict.glm(glm_YA , newdata = mydat$df |> mutate(A=1), type = 'response') |> mean(),
        predict.glm(glm_YA , newdata = mydat$df |> mutate(A=0), type = 'response') |> mean(),
        summary( glm_YA)$coefficients['A', 'Pr(>|z|)']),
      nrow=2, byrow=TRUE) |> data.frame()
    glm_res[, 1] = c('glmYAX', 'glmYA')
    colnames(glm_res) = c('method', 'numNCO', 'EY1', 'EY0', 'pval')
    
    
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
                  mutate(method = 'Pois2SLS',
                         pval   = NA) |> 
                  select(method, numNCO, EY1, EY0, pval))

    if(include_pci2sPois) {
      # TODO: 
      # perform estimation
      p2sls_result <- pci2s::p2sls.negbin(Y = mydat$df$Y, 
                                          A = mydat$df$A, 
                                          X = NULL,
                                          W = mydat$df[, mydat$W_names], 
                                          Z = mydat$df[, mydat$Z_names],
                                          nco_type = rep('poisson', ncol(W)) 
                                    )
      # add results
    }
    if(include_pci2sNB) {
      # TODO:  nco_type = rep('negbin', ncol(W)) 
      # perform estimation
      p2sls_result <- pci2s::p2sls.negbin(Y = mydat$df$Y, 
                                          A = mydat$df$A, 
                                          X = NULL,
                                          W = mydat$df[, mydat$W_names], 
                                          Z = mydat$df[, mydat$Z_names],
                                          nco_type = rep('negbin', ncol(W)) 
                                    )
      # add results

    }
    
    # calculate ATE and Fold Change
    res$ATE = res$EY1 - res$EY0
    res$FC  = res$EY1 / res$EY0
    
    return(res)
     
  }
  
  
} 


#' Get the gene expression count data (and grna assignment) from ondisc objects 
#' 
#' @param AY_idx (integer) index of AY pair to test
#' @param ZW_idx (integer) index of ZW selection 
#' @param AY (dataframe) of AY pair info
#' @param AYZW (list) of ZW info
#' @param NT_idx (vector) of idx of Non-Targeting cells
#' @param gene_odm (ondisc) manager for gene expression
#' @param grna_odm (ondisc) manager for grna expression
getCountDat <- function(AY_idx, ZW_idx,
                        AY, AYZW, NT_idx,
                        gene_odm, grna_odm,
                        max_NCs = NA) {
  
  # gather the AYZW names chosen
  AY_row = AY[AY_idx, ]
  AYZW_cur = AYZW[[AY_row$A]][[AY_row$Y]][[ZW_idx]]

  # print(paste0(AY_idx, '-', ZW_idx, ': ', AY_row$A))
  
  # Get A
  # -------------------------------------------
  A_all = grna_odm[[AY_row$A, ]] # 1:ncol(grna_odm)]]
  A_idx = which(as.logical(A_all))               # idx of all 'treated' cells
  A     = grna_odm[[AY_row$A, c(A_idx, NT_idx)]]
  
  # Get Y and Z and W
  # -------------------------------------------
  Y = gene_odm[[AY_row$Y,  c(A_idx, NT_idx)]]
  
  
  if(is.na(max_NCs)) { # load all
    Z = Matrix::t(gene_odm[[AYZW_cur$Z_names, c(A_idx, NT_idx)]]) # double [[]] pulls data
    W = Matrix::t(gene_odm[[AYZW_cur$W_names, c(A_idx, NT_idx)]]) # into memory, slower
  } else { # or load up to specified max_NCs
    Z = Matrix::t(gene_odm[[AYZW_cur$Z_names[1:max_NCs], c(A_idx, NT_idx)]]) 
    W = Matrix::t(gene_odm[[AYZW_cur$W_names[1:max_NCs], c(A_idx, NT_idx)]]) 
  }
  
  colnames(Z) = paste0('Z', 1:ncol(Z))
  colnames(W) = paste0('W', 1:ncol(W))
  
  X = gene_odm[,c(A_idx, NT_idx)] |> ondisc::get_cell_covariates() |>
      dplyr::mutate(log_n_umis = log(n_umis)) |> 
      dplyr::select(-n_umis, 
                    -n_nonzero)
  
  # put together
  df = cbind(A=Matrix::t(A) |> as.matrix() |> as.numeric(),
             Y=Matrix::t(Y) |> as.matrix() |> as.numeric(),
             Z |> as.matrix(),
             W |> as.matrix(),
             X |> as.data.frame())
  
  return(list(df = df,
              Z_names = colnames(Z),
              W_names = colnames(W),
              X_names = colnames(X)))
}


