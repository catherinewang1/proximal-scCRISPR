
# specific functions for CBEstAll pertaining to Sparse PCA implementation 
# (should work for other runs with linear continuous values and prespecified NCs, 
#  but will look specifically for spca type name saves)



require(MASS) # negative binomial fit






# REMOVE THIS: don't use this
# #' important_genes_name.rds should be in order of most to least important
# #' the rank will be the index of the gene on this list
# #' Maybe dont need this...
# #' Make function that gets the importance ranking of specified gene_name using
# #' the imput imp_gene_names parameter. The imp_gene_names should be ordered by 
# #' importance, from most important to least important.
# #' @param imp_gene_names (vector) of strings of gene names ordered by importance
# #' @example 
# #' # suppose you want a function that returns the importance ranking according 
# #' to myImpGeneNames 
# #' myGetImpRank = get_importance_rank_make(myImpGeneNames)
# #' myGetImpRank(myGeneName) # get the imp ranking of myGeneName
# get_importance_rank_make <- function(imp_gene_names) {
#   #' @param gene_name (string) gene name
#   get_importance_rank <- function(gene_name) {
#     which(imp_gene_names == gene_name)
#   }
#   return(get_importance_rank)
# }

#' Make an function that estimates ATE using the specified methods
#' Previously: get_ATE_est_NCs_make
#' 
#' @param AY (dataframe) of AY pair info
#' @param gene_norm (matrix) of gene expression (large)
#' @param NCs (matrix) of Negative Control values
#' @param grna_rownames (vector) of grna rownames from grna_odm
#' @param grna (matrix) of grna assignments (0 or 1) indicating if 
#'                      column cell received row perturbation
#' @param NT_idx (vector) of idx of Non-Targeting cells
#' @param gene_importance (dataframe) of gene names and their ranks (rank=idx in gene_norm)
#'                         must have cols: gene_name, gene_norm_idx
#' @param which_estimators (list) indicating which estimation methods to perf
#' @param save_path (string) of the save path for saving intermediate ATE est
#'                           NULL if not to be saved
#' @param U_confounders (matrix/dataframe) of unmeasured confounders 
estimate_ATE_make <- function(AY, gene_norm, NCs,
                                 grna_rownames, grna, 
                                 NT_idx, gene_importance,
                                 which_estimators, save_path=NULL, U_confounders=NULL) {

  # # function that returns the importance ranking of given gene name
  # get_importance_rank = get_importance_rank_make(imp_gene_names)  

  if(!is.null(save_path)) {    dir.create(sprintf('%s/intermediateATEs/', save_path), showWarnings = FALSE)    }
  

  #' construct a dataframe of all the data 
  #' (to be defined inside a outer function to have access to objects)
  #' Requires outside named objects:
  #' AY, NT_idx,
  #' grna_rownames, grna
  #' gene_norm
  #' get_importance_rank,
  #' @param AY_idx (integer) from 1 to nrow(AY) 
  #' @example a = format_AYZW_inner(AY_idx = 2); head(a)
  get_AYZW_df <- function(AY_idx) {
    AY_row = AY[AY_idx, ]
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
    # A = c(rep(0, length(A_idx)), rep(1, length(NT_idx))) # these should be the same
    
    # Z and W 
    # -------------------------------------------  
    # already loaded in NC values outside. should not be different for every AY pair...
    # e.g. do this outside of loop
    dfZ = NCs[c(A_idx, NT_idx), seq(from = 2, to = ncol(NCs), by = 2)] # evens which((1:ncol(NCs)) %% 2 == 0)
    dfW = NCs[c(A_idx, NT_idx), seq(from = 1, to = ncol(NCs), by = 2)] # odds
    colnames(dfZ) = paste0('Z', 1:ncol(dfZ))
    colnames(dfW) = paste0('W', 1:ncol(dfW))
    
    # Get Y
    # using already loaded in gene_norm (faster)
    # -------------------------------------------    
    # Y = gene_norm[get_importance_rank(AY_row$Y),     c(A_idx, NT_idx)]
    Y_generank = gene_importance |> dplyr::filter(gene_name == AY_row$Y) |> pull(gene_norm_idx)
    if(length(Y_generank) != 1) {}
    Y = gene_norm[Y_generank,     c(A_idx, NT_idx)]

    
    # Assemble together
    # -------------------------------------------
    dfAYZW = cbind(data.frame(A = as.vector(A), 
                              Y = Y),
                   dfZ, dfW)
    

    
    if(is.null(U_confounders)) { 
      return(dfAYZW)
    } else {
      # Add U (if given)
      # ------------------------------------------- 
      dfU = U_confounders[c(A_idx, NT_idx), ] 
      colnames(dfU) = paste0('U', 1:ncol(dfU)) # rename cols so no mix up 
      return(cbind(dfAYZW, dfU))
    }
    


  }
  

  #' Estimate ATE using variety of Estimators when specifying the
  #' A: Exposure/Treatment
  #' Y: Outcome
  #' U: 'Unmeasured' confounders (measured for comparison)
  #' Negative Controls: NCE and NCOs (ZWs)
  estimate_ATE <- function(AY_idx) {
    # === Construct df ===
    df_all   = get_AYZW_df(AY_idx=AY_idx)
    df_all$A = as.numeric(df_all$A) # convert trtmt A to numeric 0/1 if not already (alt T/F)
    
    
    res = data.frame(     method = character(0),
                          method_type = character(0),
                          numNC =   numeric(0),
                          basis = character(0),
                          ATE =   numeric(0),
                          se  =   numeric(0),
                          pval =   numeric(0),
                          time_sec = numeric(0))
    
    # === === === === === === === === ===  ===
    # === Oracle (none) and Naive (no adj) ===
    # === === === === === === === === ===  ===
    # Oracle models using unmeas conf: Y ~ A + U1 + ...
    U_colnames = grep('U', colnames(df_all), value = T)
    unmeas_conf_formula = paste0('Y ~ A + ', paste(U_colnames, collapse = ' + '))
    # lm_YAU  = lm(paste0('Y ~ A + ', paste(U_colnames, collapse = ' + ')), df)
    
    # === linear model Y ~ A (no confounder adj)
    
    if(which_estimators$lm_YA) {
      t0 = Sys.time()
      lm_YA = lm(Y ~ A, df_all)
      t1 = Sys.time()

      res = bind_rows(res, 
                      data.frame(
                        method = 'lmYA',
                        method_type = 'naive',
                        numNC = NA,
                        basis = NA,
                        ATE =    coef(lm_YA)[['A']],
                        se =  summary(lm_YA)$coefficients['A', 'Std. Error'],
                        pval= summary(lm_YA)$coefficients['A', 'Pr(>|t|)'],
                        time_sec = difftime(t1, t0, units = 'secs') |> as.numeric()))
      rm(lm_YA, t0, t1)
    }

    if(which_estimators$lm_YAU) {
      # Cannot perf lmYAU bc U are not given!! (or U df has nothing)
      if(is.null(U_confounders) | length(U_colnames) == 0) {
        # print('skip lmYAU')
        res = bind_rows(res, 
                        data.frame(
                          method = 'lmYAU',
                          method_type = 'measuredconfounders',
                          numNC = NA,
                          basis = NA,
                          ATE   = NA,
                          se    = NA,
                          pval  = NA,
                          time_sec = NA))
      } else {
        # print('perform lmYAU')
        t0 = Sys.time()
        lm_YAU = lm(unmeas_conf_formula, df_all)
        t1 = Sys.time()

        res = bind_rows(res, 
                        data.frame(
                          method = 'lmYAU',
                          method_type = 'measuredconfounders',
                          numNC = NA,
                          basis = NA,
                          ATE =    coef(lm_YAU )[['A']],
                          se =  summary(lm_YAU)$coefficients['A', 'Std. Error'],
                          pval= summary(lm_YAU)$coefficients['A', 'Pr(>|t|)'],
                          time_sec = difftime(t1, t0, units = 'secs') |> as.numeric()))
        rm(lm_YAU, t0, t1)
      }
      
    }
    
    # Separately, compare to count analysis (pois, nb)
    # # pois_YAU: NO LONGER 'ATE'! (return the coef, compare the p-value)
    # if(which_estimators$pois_YAU) {
    #   pois_YAU = glm(Y ~ A + U, df_all, family = 'poisson')
    #   res = bind_rows(res, 
    #                   data.frame(
    #                     method = 'poisYAU',
    #                     method_type = 'measuredconfounders',
    #                     numNC = NA,
    #                     basis = NA,
    #                     ATE =    coef(pois_YAU )[['A']],
    #                     se = NA,
    #                     pval= summary(pois_YAU)$coefficients['A', 'Pr(>|t|)']))
    #   rm(pois_YAU)
    # }

    # # nb_YAU: NO LONGER 'ATE'! (return the coef, compare the p-value)
    # if(which_estimators$nb_YAU) {
    #   nb_YAU = MASS::glm.nb(Y ~ A + U, df_all)
    #   res = bind_rows(res, 
    #                   data.frame(
    #                     method = 'nbYAU',
    #                     method_type = 'measuredconfounders',
    #                     numNC = NA,
    #                     basis = NA,
    #                     ATE =    coef(nb_YAU )[['A']],
    #                     se = NA,
    #                     pval= summary(nb_YAU)$coefficients['A', 'Pr(>|t|)']))
    #   rm(nb_YAU)
    # }

    
    # t0 = Sys.time()
    for(num_NCs in num_NC_pairs) {
      # print(num_NCs)
      df = df_all[, c('A', 'Y', paste0('Z', 1:num_NCs), paste0('W', 1:num_NCs))] # subset only these cols
      
      
      # === === === === === === === === ===  ===
      # === Estimators using original values ===
      # === === === === === === === === ===  ===
      
      
      
      # === Outcome Confounding Bridge, 2 Stage LS w/ pci2s package
      if(which_estimators$OCB_2SLS_pci2s) {
        t0 = Sys.time()
        pci2s_res = pci2s::p2sls.lm(
          Y = df$Y, 
          A = df$A, 
          W = df[,grepl('W', colnames(df))], 
          Z = df[,grepl('Z', colnames(df))], 
          variance = TRUE)
        t1 = Sys.time()

        res = bind_rows(res, 
                        data.frame(
                          method = '2SLSpci2s',
                          method_type = 'proximal',
                          numNC = num_NCs,
                          basis = NA,
                          ATE = pci2s_res$summary_second_stage['A', 'Estimate'],
                          se  = pci2s_res$summary_second_stage['A', 'Std. Error'],
                          pval= pci2s_res$summary_second_stage['A', 'Pr(>|z|)'],
                         time_sec = difftime(t1, t0, units = 'secs') |> as.numeric()))
        rm(pci2s_res, t0, t1)
      }
      
     
      
      
    
      
      
      
    }
    
    
    # # 1 run: 39 sec w/  pci2s (no LinOS, all others)
    # # 1 run:  4 sec w/o pci2s (no LinOS, all others)
    # t0 = Sys.time()
    # get_ATE_est_NCs(AY_idx=2)
    # t1 = Sys.time(); print(t1 - t0)
    
    
    
    # save ATEs one by one (AY_idx, ZW_idx) if wanted
    if(!is.null(save_path)) {
      write.csv(res,
                sprintf('%s/intermediateATEs/ATE_%d.csv', save_path, AY_idx),
                row.names = FALSE)
    }
    
    return(res)
  }

  
  return(estimate_ATE)
}








#' Creates a function that estimates effects (return coefficient of treatment A, not the 'ATE'!) 
#' for count methods (poisson, negative binomial)
#' This function takes AY_idx as input, and returns a df of estimates w/ other estimation info
#' @param AY (dataframe) of AY pair info
#' @param which_estimators (list) indicating which estimation methods to perf
#' @param gene_odm (on disc manager) of gene data (<- mught change to memory loaded version for speed)
#' @param grna_odm (on disc manager) of grna data (<- mught change to memory loaded version for speed)
#' @param NT_idx (vector) of idx of Non-Targeting cells
#' @param U_confounders (matrix/dataframe) of unmeasured confounders 
#' @param save_path (string) of the save path for saving intermediate ATE est
#'                           NULL if not to be saved
estimate_effect_count_make <- function(AY, 
                                       which_estimators,
                                       gene_odm, # <- maybe change to loaded in memory (for speed)
                                       grna_odm,
                                       NT_idx,
                                       U_confounders,
                                       save_path=NULL
                                       ) {


  if(!is.null(save_path)) {    dir.create(sprintf('%s/intermediateATEsCountGLM/', save_path), showWarnings = FALSE)    }
  
  #' construct a dataframe of all the data 
  #' (to be defined inside a outer function to have access to objects)
  #' Requires outside named objects:
  #' AY, NT_idx,
  #' gene_odm, grna_odm,
  #' U_confounders,
  #' @param AY_idx (integer) from 1 to nrow(AY) 
  #' @example a = format_AYZW_inner(AY_idx = 2); head(a)
  get_df <- function(AY_idx) {
    # print('get AY row')
    AY_row = AY[AY_idx, ]
    # Get A
    # -------------------------------------------
    # print('get A')
    # idx of all 'treated' cells
    A_idx = which(as.logical(grna_odm[[AY_row$A, ]]))
    # subset cells of 'treated' (w A grna) and 'control' (NT gran)
    A = grna_odm[[AY_row$A, c(A_idx, NT_idx)]]
    
    # # using grna all loaded into memory
    # # idx of all 'treated' cells
    
    # A_grna_idx = which(grna_rownames == AY_row$A)
    # A_idx = which(as.logical(grna[A_grna_idx, 1:ncol(grna)]))
    # # subset cells of 'treated' (w A grna) and 'control' (NT gran)
    # A = grna[A_grna_idx, c(A_idx, NT_idx)]
    # A = c(rep(0, length(A_idx)), rep(1, length(NT_idx))) # should they not just be 0/1??
    
    # # Z and W (save for later for possible pci2s pois and nb version)
    # # -------------------------------------------  
    # # already loaded in NC values outside. should not be different for every AY pair...
    # # e.g. do this outside of loop
    # dfZ = NCs[c(A_idx, NT_idx), seq(from = 2, to = ncol(NCs), by = 2)] # evens which((1:ncol(NCs)) %% 2 == 0)
    # dfW = NCs[c(A_idx, NT_idx), seq(from = 1, to = ncol(NCs), by = 2)] # odds
    # colnames(dfZ) = paste0('Z', 1:ncol(dfZ))
    # colnames(dfW) = paste0('W', 1:ncol(dfW))
    
    # Get Y
    # using already loaded in gene_norm (faster)
    # -------------------------------------------
    # print('get Y')    
    Y = gene_odm[[AY_row$Y,     c(A_idx, NT_idx)]]
    

    # Assemble together
    # -------------------------------------------   
    # print('assemble df') 
    dfAY = data.frame(A = as.vector(A), 
                      Y = as.vector(Y))
    
    # print('add u confounders')
    if(is.null(U_confounders)) { 
      return(dfAY)
    } else {
      # Add U (if given)
      # ------------------------------------------- 
      dfU = U_confounders[c(A_idx, NT_idx), ] 
      colnames(dfU) = paste0('U', 1:ncol(dfU)) # rename cols so no mix up 
      return(cbind(dfAY, dfU))
    }
    
  }

  estimate_effect_count <- function(AY_idx) {

    # === Construct df ===
    # print('get df')
    df_all   = get_df(AY_idx=AY_idx)
    df_all$A = as.numeric(df_all$A) # convert trtmt A to numeric 0/1 if not already (alt T/F)
    
    
    res = data.frame(     method = character(0),
                          method_type = character(0),
                          numNC =   numeric(0),
                          basis = character(0),
                          ATE =   numeric(0),
                          se  =   numeric(0),
                          pval =   numeric(0))
    
    # print('start estimating')

    # Oracle models using (un)meas conf: Y ~ A + U1 + U2 +...
    U_colnames = grep('U', colnames(df_all), value = T)
    unmeas_conf_formula = paste0('Y ~ A + ', paste(U_colnames, collapse = ' + '))
    # print(sprintf('unmeas conf formula: %s', unmeas_conf_formula))
    # === Poisson Y ~ A      (no confounder adj)
    if(which_estimators$pois_YA) {
      # print('pois_YA')
      pois_YA = glm('Y ~ A', df_all, family = 'poisson')
      res = bind_rows(res, 
                      data.frame(
                        method = 'poisYA',
                        method_type = 'naive',
                        numNC = NA,
                        basis = NA,
                        ATE =    coef(pois_YA )[['A']],
                        se  = summary(pois_YA)$coefficients['A', 'Std. Error'],
                        pval= summary(pois_YA)$coefficients['A', 'Pr(>|z|)']))
      rm(pois_YA)
    }

    # === Poisson Y ~ A + Us (   confounder adj)
    if(which_estimators$pois_YAU) {
      # print('pois_YAU')
      pois_YAU = glm(unmeas_conf_formula, df_all, family = 'poisson')
      res = bind_rows(res, 
                      data.frame(
                        method = 'poisYAU',
                        method_type = 'measuredconfounders',
                        numNC = NA,
                        basis = NA,
                        ATE =    coef(pois_YAU )[['A']],
                        se  = summary(pois_YAU)$coefficients['A', 'Std. Error'],
                        pval= summary(pois_YAU)$coefficients['A', 'Pr(>|z|)']))
      rm(pois_YAU)
    }


    # === Negative Binomial Y ~ A      (no confounder adj)
    if(which_estimators$nb_YA) {
      # print('nb_YA')
      nb_YA = glm('Y ~ A', df_all, family = 'poisson')
      # print(summary(nb_YA)$coefficients)
      res = bind_rows(res, 
                      data.frame(
                        method = 'nbYA',
                        method_type = 'measuredconfounders',
                        numNC = NA,
                        basis = NA,
                        ATE =    coef(nb_YA )[['A']],
                        se  = summary(nb_YA)$coefficients['A', 'Std. Error'],
                        pval= summary(nb_YA)$coefficients['A', 'Pr(>|z|)']))
      rm(nb_YA)
    }

    # === Negative Binomial Y ~ A + Us (   confounder adj)
    if(which_estimators$nb_YAU) {
      # print('nb_YAU')
      nb_YAU = glm(unmeas_conf_formula, df_all, family = 'poisson')
      res = bind_rows(res, 
                      data.frame(
                        method = 'nbYAU',
                        method_type = 'measuredconfounders',
                        numNC = NA,
                        basis = NA,
                        ATE =    coef(nb_YAU )[['A']],
                        se  = summary(nb_YAU)$coefficients['A', 'Std. Error'],
                        pval= summary(nb_YAU)$coefficients['A', 'Pr(>|z|)']))
      rm(nb_YAU)
    }
  

    # save ATEs one by one (AY_idx, ZW_idx) if wanted
    if(!is.null(save_path)) {
      write.csv(res,
                sprintf('%s/intermediateATEsCountGLM/ATE_%d.csv', save_path, AY_idx),
                row.names = FALSE)
    }
     
    return(res)

  }

  return(estimate_effect_count)
}
































# =============================================================================
#                  TRASH
# =============================================================================

#' Get SE and pval of estimate when asymptotically
#' \sqrt{n} (betahat - beta)/\sqrt(betavar) --> N(0, 1)
#' @param betahat estimate of beta
#' @param betavar asymptotic variance of estimate
#' @param n sample size
#' @param beta true beta value or null hypothesis test for beta (typically 0)
get_se_pval <- function(betahat, betavar, n, beta = 0) {
  se = sqrt(betavar) / sqrt(n)
  zscore = (betahat - beta) / se
  pval = pnorm(-abs(zscore)) * 2
  return(list(se=se, pval=pval))
}

#' Make a function that takes in AY_idx (for easy looping/parallelization)
#' and returns ATE estimates of the specified estimators
#' and outer parallel function call will loop through AY rows
#' 
#' Use of an outer make function is to predefine objects used across many calls
#' @param AY (dataframe) of AY pair info
#' @param gene_norm (matrix) of gene expression (large)
#' @param NCs (matrix) of Negative Control values
#' @param grna_rownames (vector) of grna rownames from grna_odm
#' @param grna (matrix) of grna assignments (0 or 1) indicating if 
#'                      column cell received row perturbation
#' @param NT_idx (vector) of idx of Non-Targeting cells
#' @param imp_gene_names (vector) ordered vector of most important gene names
#' @param which_estimators (list) indicating which estimation methods to perf
#' @param CB_setting (list) of CB method settings
#' @param save_path (string) of the save path for saving intermediate ATE est
#'                           NULL if not to be saved
get_ATE_est_NCs_make_OLD <- function(AY, gene_norm, NCs,
                                 grna_rownames, grna, 
                                 NT_idx, imp_gene_names,
                                 which_estimators, CB_setting, save_path=NULL) {
  
  # function that returns the importance ranking of given gene name
  get_importance_rank = get_importance_rank_make(imp_gene_names)  
  
  
  
  
  # Make CB_settings accessible 
  ALPHA     = CB_setting$ALPHA # R does not save as pointers, so ok?
  GMM_steps = CB_setting$GMM_steps
  K_folds   = CB_setting$K_folds
  lambdas   = CB_setting$lambdas
  TRIM      = CB_setting$TRIM
  basisParams  = CB_setting$basisParams
  num_NC_pairs = CB_setting$num_NC_pairs
  
  if(!is.null(save_path)) {    dir.create(sprintf('%s/intermediateATEs/', save_path), showWarnings = FALSE)    }
  
  #' construct a dataframe of all the data 
  #' (to be defined inside a outer function to have access to objects)
  #' Requires outside named objects:
  #' AY, NT_idx,
  #' grna_rownames, grna
  #' gene_norm
  #' get_importance_rank,
  #' dfZ, dfW
  #' @param AY_idx (integer) from 1 to nrow(AY) 
  #' @example a = format_AYZW_inner(AY_idx = 2); head(a)
  format_AYZW_inner <- function(AY_idx) {
    AY_row = AY[AY_idx, ]
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
    
    # Z and W 
    # -------------------------------------------  
    # already loaded in NC values outside. should not be different for every AY pair...
    # e.g. do this outside of loop
    dfZ = NCs[c(A_idx, NT_idx), seq(from = 2, to = ncol(NCs), by = 2)] # evens which((1:ncol(NCs)) %% 2 == 0)
    dfW = NCs[c(A_idx, NT_idx), seq(from = 1, to = ncol(NCs), by = 2)] # odds
    colnames(dfZ) = paste0('Z', 1:ncol(dfZ))
    colnames(dfW) = paste0('W', 1:ncol(dfW))
    
    # Get Y
    # using already loaded in gene_norm (faster)
    # -------------------------------------------    
    Y = gene_norm[get_importance_rank(AY_row$Y),     c(A_idx, NT_idx)]
    
    # Assemble together
    # -------------------------------------------
    dfAY = data.frame(A = as.vector(A), 
                      Y = Y)
    
    return(cbind(dfAY, dfZ, dfW))
  }
  
  #' Estimate ATE using variety of CB/Proximal Estimators when specifying the
  #' A: Exposure/Treatment
  #' Y: Outcome
  #' Negative Controls: NCE and NCOs (ZWs)
  get_ATE_est <- function(AY_idx) {
    # === Construct df ===
    df_all   = format_AYZW_inner(AY_idx=AY_idx)
    df_all$A = as.numeric(df_all$A) # convert trtmt A to numeric 0/1 if not already (alt T/F)
    
    
    res = data.frame(     method = character(0),
                          method_type = character(0),
                          numNC =   numeric(0),
                          basis = character(0),
                          ATE =   numeric(0),
                          se  =   numeric(0),
                          pval =   numeric(0))
    
    # === === === === === === === === ===  ===
    # === Oracle (none) and Naive (no adj) ===
    # === === === === === === === === ===  ===
    # # Oracle linear model Y ~ A + Us (NOT HERE IN REAL DATA? SHOULD WE TRY?)
    # U_colnames = grep('U', colnames(df), value = T)
    # lm_YAU  = lm(paste0('Y ~ A + ', paste(U_colnames, collapse = ' + ')), df)
    
    # === linear model Y ~ A (no confounder adj)
    
    if(which_estimators$lm_YA) {
      lm_YA = lm(Y ~ A, df_all)
      res = bind_rows(res, 
                      data.frame(
                        method = 'lmYA',
                        method_type = 'naive',
                        numNC = NA,
                        basis = NA,
                        ATE =    coef(lm_YA )[['A']],
                        se = NA,
                        pval= summary(lm_YA)$coefficients['A', 'Pr(>|t|)']))
    }
    
    
    # t0 = Sys.time()
    for(num_NCs in num_NC_pairs) {
      # print(num_NCs)
      df = df_all[, get_CB_colnames(num_NCs)] # subset only these cols
      
      
      # === === === === === === === === ===  ===
      # === Estimators using original values ===
      # === === === === === === === === ===  ===
      
      ## === Outcome Confounding Bridge, 2 Stage LS
      if(which_estimators$OCB_2SLS) {
        OCB_2SLS = OCB2SLS(df)$ATE
        
        res = bind_rows(res, 
                        data.frame(
                          method = 'OCB2SLS',
                          method_type = '2SLS',
                          numNC = num_NCs,
                          basis = NA,
                          ATE   = OCB_2SLS,
                          se    = NA,
                          pval  = NA))
      }
      
      
      # === Outcome Confounding Bridge, 2 Stage LS w/ pci2s package
      if(which_estimators$OCB_2SLS_pci2s) {
        OCB_2SLS_pci2s = pci2s::p2sls.lm(
          Y = df$Y, 
          A = df$A, 
          W = df[,grepl('W', colnames(df))], 
          Z = df[,grepl('Z', colnames(df))], 
          variance = TRUE)
        
        res = bind_rows(res, 
                        data.frame(
                          method = 'OCB2SLSpci2s',
                          method_type = '2SLS',
                          numNC = num_NCs,
                          basis = NA,
                          ATE = OCB_2SLS_pci2s$summary_second_stage['A', 'Estimate'],
                          se  = OCB_2SLS_pci2s$summary_second_stage['A', 'Std. Error'],
                          pval= OCB_2SLS_pci2s$summary_second_stage['A', 'Pr(>|z|)']))
      }
      
      
      # === Outcome Confounding Bridge, 2 Stage LS w/ ElNet Reg
      if(which_estimators$OCB_2SLSReg) {
        OCB_2SLSReg_res = OCB2SLSReg(df, alpha=ALPHA, returnMiddleDFs=FALSE)
        OCB_2SLSReg     = OCB_2SLSReg_res[[length(OCB_2SLSReg_res)]][['ATE']]
        
        res = bind_rows(res, 
                        data.frame(
                          method      = 'OCB2SLSReg',
                          method_type = '2SLS',
                          numNC = num_NCs,
                          basis = NA,
                          ATE   = OCB_2SLSReg,
                          se    = NA,
                          pval  = NA))
      }
      
      
      # # Bad but,... select the latest ATE estimate that is not NA (increases by #W's)
      # OCB_2SLSReg_res_ATEs = c()
      # for(i in 1:length(OCB_2SLSReg_res)) {
      #   OCB_2SLSReg_res_ATEs = c(OCB_2SLSReg_res_ATEs,
      #                            OCB_2SLSReg_res[[i]][['ATE']])
      #   # OCB_2SLSReg_res_ATEsOCB_2SLSReg_res[[i]][['ATE']] |> print()
      # }
      # plot(OCB_2SLSReg_res_ATEs)
      # OCB_2SLSReg = OCB_2SLSReg_res_ATEs[max(which(!is.na(OCB_2SLSReg_res_ATEs)))]
      
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
        if(which_estimators$OCB_GMM) {
          OCBGMM_res = OCBGMM(dl, returnAlphas = T)
          
          GMM_coverage = checkCoverage(dl=dl, 
                                       alpha1 = OCBGMM_res$alpha1,
                                       alpha0 = OCBGMM_res$alpha0, 
                                       Omega  = diag(rep(1, ncol(dl$B))),
                                       coverage= NA, 
                                       ATE_est = OCBGMM_res$ATE, 
                                       ATE_true= 0, 
                                       returnVar=TRUE)
          
          OCBGMM_se_pval = get_se_pval(betahat = OCBGMM_res$ATE, 
                                       betavar = GMM_coverage$VarATE, 
                                       n = nrow(df), beta = 0)
          res = bind_rows(res, 
                          data.frame(
                            method      = 'OCBMEst',
                            method_type = 'MEst',
                            numNC       = num_NCs,
                            basis       = basis,
                            ATE         = OCBGMM_res$ATE,
                            se          = OCBGMM_se_pval$se,
                            pval        = OCBGMM_se_pval$pval))
        }
        
        
        
        
        
        # === OCBGMM ReWeight
        if(which_estimators$OCB_GMM) {
          OCBGMMRw_res = OCBGMMRw(dl, GMM_steps, returnWeights = T)
          
          GMMRw_coverage = checkCoverage(    dl=dl, 
                                             alpha1  = t(OCBGMMRw_res$alpha1s[GMM_steps, , drop=FALSE]),
                                             alpha0  = t(OCBGMMRw_res$alpha0s[GMM_steps, , drop=FALSE]), 
                                             Omega   = OCBGMMRw_res$Weights[[length(OCBGMMRw_res$Weights)]],
                                             coverage= NA, 
                                             ATE_est = OCBGMMRw_res$ATE, 
                                             ATE_true= 0, 
                                             returnVar=TRUE)
          # SE?? and pval?
          # OCBGMMRw_se   = sqrt(GMMRw_coverage$VarATE) / sqrt(nrow(df))
          # OCBGMMRw_pval = pnorm(-abs(OCBGMMRw_res$ATE - 0)/OCBGMMRw_se) * 2
          OCBGMMRw_se_pval = get_se_pval(betahat = OCBGMMRw_res$ATE, 
                                         betavar = GMMRw_coverage$VarATE, 
                                         n = nrow(df), beta = 0)
          
          
          res = bind_rows(res, 
                          data.frame(
                            method      = 'OCBMEstRW',
                            method_type = 'MEst',
                            numNC       = num_NCs,
                            basis       = basis,
                            ATE         = OCBGMMRw_res$ATE,
                            se          = OCBGMMRw_se_pval$se,
                            pval        = OCBGMMRw_se_pval$pval))
        }
        
        # === OCBGMM ReWeight Regularization
        if(which_estimators$OCB_GMM) {
          OCBGMMRWReg_res = OCBGMMRwReg(dl, GMM_steps, K_folds, lambdas, 
                                        returnWeights = F, returnLambdas = F)
          
          res = bind_rows(res, 
                          data.frame(
                            method      = 'OCBMEstRWReg',
                            method_type = 'MEst',
                            numNC       = num_NCs,
                            basis       = basis,
                            ATE         = OCBGMMRWReg_res$ATE,
                            se          = NA,
                            pval        = NA))
        }
        
        # NOT IMPLEMENTED: LIN OS (PI AND OS AND TRIM)
        # === OCB Linear CB Plug-In AND One-Step 
        #      Z0,Z1 NULL --> use transformations in B0,B1 as covariates 
        #                     when estimating nuisance params
        if(which_estimators$OCB_LinOS | which_estimators$OCB_LinOStrim) {
          # # print('Running OCB Lin OS Estimation')
          # # No trimming
          # OCBLinOS_res = tryCatch({estCondMomentOCBOS(dataList=dl, gmmSteps=GMM_steps, Z0=NULL, Z1=NULL)},
          #                         error = function(cond) {
          #                           # message(sprintf('Error est CondMomentOCBOS with %s', 
          #                           #                 gammaSetting)) 
          #                           return(NULL)
          #                         }) 
          # # OCBLinOS_res = estCondMomentOCBOS(dataList=dl, gmmSteps=GMM_steps, Z0=NULL, Z1=NULL)
          # # errors sometimes... NA for now... debug later
          # if(is.null(OCBLinOS_res)) {
          #   res_numNCs_basis[['OCBLinPI']] = NA
          #   res_numNCs_basis[['OCBLinOS']] = NA
          #   # res[[paste0('OCBLinPI_', basis)]] = NA
          #   # res[[paste0('OCBLinOS_', basis)]] = NA
          # } else {
          #   res_numNCs_basis[['OCBLinPI']] = OCBLinOS_res$ATE_pi
          #   res_numNCs_basis[['OCBLinOS']] = OCBLinOS_res$ATE_os
          #   # res[[paste0('OCBLinPI_', basis)]] = OCBLinOS_res$ATE_pi
          #   # res[[paste0('OCBLinOS_', basis)]] = OCBLinOS_res$ATE_os
          # }
          # 
          # # With Trimming (at .1)
          # OCBLinOS_res = tryCatch({estCondMomentOCBOS(dataList=dl, gmmSteps=GMM_steps, 
          #                                             Z0=NULL, Z1=NULL, trim=TRIM)},
          #                         error = function(cond) {
          #                           # message(sprintf('Error est CondMomentOCBOS with %s', 
          #                           #                 gammaSetting)) 
          #                           return(NULL)
          #                         }) 
          # if(is.null(OCBLinOS_res)) {
          #   res_numNCs_basis[['OCBLinOStrim']] = NA
          #   # res[[paste0('OCBLinOStrim_', basis)]] = NA
          # } else {
          #   res_numNCs_basis[['OCBLinOStrim']] = OCBLinOS_res$ATE_os
          #   # res[[paste0('OCBLinOStrim_', basis)]] = OCBLinOS_res$ATE_os
          # }
          
          
        } else if(which_estimators$OCB_LinOSPI) { # only Plug-In Estimator
          # # print('Running OCB Lin PI (only) Estimation')
          # OCBLinPI_res = tryCatch({estCondMomentOCBPI(dataList=dl, gmmSteps=GMM_steps, Z0=NULL, Z1=NULL)},
          #                         error = function(cond) {
          #                           # message(sprintf('Error est CondMomentOCBOS with %s', 
          #                           #                 gammaSetting)) 
          #                           return(NULL)
          #                         })
          # # errors sometimes... NA for now... debug later
          # if(is.null(OCBLinPI_res)) {
          #   res_numNCs_basis[['OCBLinPI']] = NA
          # } else {
          #   res_numNCs_basis[['OCBLinPI']] = OCBLinPI_res$ATE
          # }
        }
        
      }
      
      
      
    }
    
    
    # # 1 run: 39 sec w/  pci2s (no LinOS, all others)
    # # 1 run:  4 sec w/o pci2s (no LinOS, all others)
    # t0 = Sys.time()
    # get_ATE_est_NCs(AY_idx=2)
    # t1 = Sys.time(); print(t1 - t0)
    
    
    
    # save ATEs one by one (AY_idx, ZW_idx) if wanted
    if(!is.null(save_path)) {
      write.csv(res,
                sprintf('%s/intermediateATEs/ATE_%d.csv', save_path, AY_idx),
                row.names = FALSE)
    }
    
    return(res)
  }
  
  
  return(get_ATE_est)
}
