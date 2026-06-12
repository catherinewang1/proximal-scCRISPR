
  #' 
  #' construct a dataframe of all the data 
  #' (outside general function)
  #' 
  #' @param A_name (character) name of grna/perturbation treatment
  #' @param Y_name (character) name of gene outcome
  #' @param NT_idx (vector) of indices for cells receiving non-targeting perturbations
  #' @param grna (matrix) of 0/1 of size #grnas x #cells indicating grna assignment
  #'    make sure the rownames of grna are the names of the grna (e.g. A_name is one of the rownames)
  #' @param gene_norm (matrix) of gene expr outcomes (should be continuous for this pci2s)
  #' @param NCs (matrix or list)
  #'    matrix if using PCA, Sparse PCA, or WGCNA (or etc...)
  #'    list if using individualgene: input here as a list with the structure:
  #'          NCs$Z_names   (vector) of characters of gene names
  #'          NCs$W_names   (vector) of characters of gene names
  #' @example a = format_AYZW_inner(AY_idx = 2); head(a)
  get_AYZW_df_pci2sbyNCs <- function(A_name, Y_name, NT_idx,
                                     grna, 
                                     gene_norm=NULL, 
                                     NCs=NULL, 
                                     U_confounders=NULL) {
    
    # Get A
    # -------------------------------------------
    # idx of all 'treated' cells    
    # A_grna_idx = which(row.names(grna) == A_name)       # idx of A grna
    # A_idx      = which(as.logical(grna[A_grna_idx, ]))  # idx of cells receiving this A grna
    A_idx      = which(as.logical(grna[A_name, ]))        # idx of cells receiving this A grna
    # subset cells of 'treated' (w A grna) and 'control' (NT gran)
    A = grna[A_name, c(A_idx, NT_idx)]
    # A = c(rep(0, length(A_idx)), rep(1, length(NT_idx))) # these should be the same
    
    # Get Z and W 
    # -------------------------------------------  
    if(is.matrix(NCs)) {
      dfZ = NCs[c(A_idx, NT_idx), seq(from = 2, to = ncol(NCs), by = 2)] # evens which((1:ncol(NCs)) %% 2 == 0)
      dfW = NCs[c(A_idx, NT_idx), seq(from = 1, to = ncol(NCs), by = 2)] # odds
    } else if(is.list(NCs) && ('Z_names' %in% names(NCs)) && ('W_names' %in% names(NCs)) ) {
      dfZ = gene_norm[NCs[['Z_names']], c(A_idx, NT_idx)] |> t() |> data.frame()
      dfW = gene_norm[NCs[['W_names']], c(A_idx, NT_idx)] |> t() |> data.frame()
    } else {
      print(sprint('(A:%s,Y:%s) Bad NCs input in get_AYZW_df_pci2sbyNCs', A_name, Y_name))
      return()
    }
    colnames(dfZ) = paste0('Z', 1:ncol(dfZ))
    colnames(dfW) = paste0('W', 1:ncol(dfW))
    
    # Get Y
    # using already loaded in gene_norm (faster)
    # -------------------------------------------    
    # Y = gene_norm[get_importance_rank(AY_row$Y),     c(A_idx, NT_idx)]
    # Y_generank = gene_importance |> dplyr::filter(gene_name == AY_row$Y) |> pull(gene_norm_idx)
    # if(length(Y_generank) != 1) {}
    # Y = gene_norm[Y_generank,     c(A_idx, NT_idx)]
    Y = gene_norm[Y_name,     c(A_idx, NT_idx)]
    
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

#' Make an function that estimates ATE using the various Negative Controls constructed
#' Previously: get_ATE_est_NCs_make
#' 
#' @param AY (dataframe) of AY pair info
#' @param gene_norm (matrix) of gene expression (large)
#' @param grna (matrix) of grna assignments (0 or 1) indicating if 
#'                      column cell received row perturbation
#' @param NT_idx (vector) of idx of Non-Targeting cells
#' @param NCs_list (list) of Negative Control matrices
#' @param NCs_settings (list) of Negative Control settings, e.g.
#'       $num_NC_pairs: how many pairs of NCs to use in proximal est
#'       $which_estimators: indicating which estimation methods to perf
#' @param save_path (string) of the save path for saving intermediate ATE est
#'                           NULL if not to be saved
#' @param U_confounders (matrix/dataframe) of unmeasured confounders 
estimate_ATE_pci2sbyNC_make <- function(AY, gene_norm, grna, NT_idx,
                                        NCs_list, NCs_settings,
                                        save_path=NULL, U_confounders=NULL) {

  # # function that returns the importance ranking of given gene name
  # get_importance_rank = get_importance_rank_make(imp_gene_names)  

  if(!is.null(save_path)) {    dir.create(sprintf('%s/intermediateATEs/', save_path), showWarnings = FALSE)    }
  

  
  

  #' Estimate ATE using variety of Estimators when specifying the
  #' A: Exposure/Treatment
  #' Y: Outcome
  #' U: 'Unmeasured' confounders (measured for comparison)
  #' Negative Controls: NCE and NCOs (ZWs)
  #' @param (integer) an integer to specify which row of AY to test/estimate for
  estimate_ATE <- function(AY_idx) {
    res = data.frame(       NC_type=character(0),
                            method = character(0),
                            method_type = character(0),
                            numNC =   numeric(0),
                            basis = character(0),
                            ATE =   numeric(0),
                            se  =   numeric(0),
                            pval =   numeric(0),
                            time_sec = numeric(0))
    
    performed_lmYA  = FALSE # only need to estimate using lm one time but can be done with any of the NCs_list loops
    performed_lmYAU = FALSE # set these to TRUE when these methods are done
    for(NC_name in names(NCs_list)) {
      # print(sprintf('[%s] %s', Sys.time(), NC_name))
      NCs = NCs_list[[NC_name]]
      
      # === Construct df ===
      if(!is.matrix(NCs)) { # edit NCs into smaller list with just Z_names and W_names, and limit num_NC_pairs actually used
        # print(sprintf('%s Making NC as list of Z and Ws!', NC_name))
        NCs_new = list(Z_names = NCs[[AY[AY_idx, 'A']]][[AY[AY_idx, 'Y']]][[1]]$Z_names[1:max(NCs_settings[[NC_name]]$num_NC_pairs)],
                       W_names = NCs[[AY[AY_idx, 'A']]][[AY[AY_idx, 'Y']]][[1]]$W_names[1:max(NCs_settings[[NC_name]]$num_NC_pairs)])
        NCs = NCs_new
        rm(NCs_new)
      }
      
      df_all = get_AYZW_df_pci2sbyNCs(
            A_name       =AY[AY_idx, 'A'], 
            Y_name       =AY[AY_idx, 'Y'], 
            NT_idx       =NT_idx,
            grna         =grna, 
            gene_norm    =gene_norm, 
            NCs          =NCs, 
            U_confounders=U_confounders) 
      
      df_all$A = as.numeric(df_all$A) # convert trtmt A to numeric 0/1 if not already (alt T/F)
      
      
      
      which_estimators = NCs_settings[[NC_name]]$which_estimators
      # === === === === === === === === ===  ===
      # === Oracle (none) and Naive (no adj) ===
      # === === === === === === === === ===  ===
      # Oracle models using unmeas conf: Y ~ A + U1 + ...
      U_colnames = grep('U', colnames(df_all), value = T)
      unmeas_conf_formula = paste0('Y ~ A + ', paste(U_colnames, collapse = ' + '))
      
      # === linear model Y ~ A (no confounder adj)
      if(which_estimators$lm_YA && !performed_lmYA) { # if specified and not already done
        t0 = Sys.time()
        lm_YA = lm(Y ~ A, df_all)
        t1 = Sys.time()

        res = bind_rows(res, 
                        data.frame(
                          NC_type=NA,
                          method = 'lmYA',
                          method_type = 'naive',
                          numNC = NA,
                          basis = NA,
                          ATE =    coef(lm_YA)[['A']],
                          se =  summary(lm_YA)$coefficients['A', 'Std. Error'],
                          pval= summary(lm_YA)$coefficients['A', 'Pr(>|t|)'],
                          time_sec = difftime(t1, t0, units = 'secs') |> as.numeric()))
        rm(lm_YA, t0, t1)
        performed_lmYA = TRUE
      }
      # === linear model Y ~ A + U (with confounder adj)
      if(which_estimators$lm_YAU && !performed_lmYAU) { # if specified and not already done
        # Cannot perf lmYAU bc U are not given!! (or U df has nothing)
        if(is.null(U_confounders) | length(U_colnames) == 0) {
          # print('skip lmYAU')
          res = bind_rows(res, 
                          data.frame(
                            NC_type=NA,
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

        performed_lmYAU = TRUE 
        
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
      for(num_NCs in NCs_settings[[NC_name]]$num_NC_pairs) {
        # print(num_NCs)
        chosen_cols =  c('A', 'Y', paste0('Z', 1:num_NCs), paste0('W', 1:num_NCs))
        if(!all(chosen_cols %in% colnames(df_all))) { # if not all NCs are available (e.g. not enough NCs for the specified num_NCs)
          next
        }
        df = df_all[, chosen_cols] # subset only these cols
        
        
        # === === === === === === === === ===  === === ===
        # === Proximal Estimators (only pci2s here)    ===
        # === === === === === === === === ===  === === ===
            
        # === Proximal/Outcome Confounding Bridge, 2 Stage LS w/ pci2s package
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
                            NC_type     = NC_name,
                            method      = '2SLSpci2s',
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





