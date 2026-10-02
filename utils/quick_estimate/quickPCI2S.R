







#' estimate ATEs/logFoldChanges/Coefficients using loaded in vars in environment
#' Specifically: 
#'    grna, gene, gene_norm, AY, validsinglegenes
#'    libsize_log, control_cell_idx, NCs_pca, numNCs
#'    
#' res is dataframe with cols:
#' - grna, gene, type 
#' - ... etc ...
#' - tstat: test statistic either the t value (lm) or z value (glm, proximal)
#' @param AY_idx (integer) index of the AY test 
#' @param AY (dataframe) of AY tests, must have columnss 'grna' and 'gene'
#' @param numNCs (vector) of number of Negative Control pairs to be used for proximal estimates
#' @param control_cell_idx (vector) of integer indices for the control cells
#'       (columns of gene_norm, gene, grna) and (rows of covariates)
#' @param gene_norm (matrix) of gene x cell continuous gene expression
#' @param gene (matrix) of gene x cell count gene expression
#' @param grna (matrix) of perturbation/grna x cell assignment
#' @param libsize_log (vector) of log library size (log nUMIs)
#' @param NCs_pca (matrix/dataframe) of PCA loadings for NC variables
#' @param validsinglegenes (vector) of characters 
#' @param covariates (matrix/dataframe) of cell x features of potential 
#'                   covariates to add to the models
#' @param save_intermediateATEs (boolean) whether to save this current AY test in intermediateATEs/ subfolder
#' @param skip_if_saved_intermediateATEs (boolean) whether to skip computation if this AY test is already saved in intermediateATEs/ subfolder
#' @param run_proximal_linear_pca (boolean) whether to run linear proximal with PCA NCs
#' @param run_proximal_negbin_pca (boolean) whether to run neg bin proximal with PCA NCs
#' @param run_proximal_linear_singlegenes (boolean) whether to run linear proximal with singlegenes NCs
#' @param run_proximal_negbin_singlegenes (boolean) whether to run neg bin proximal with singlegenes NCs
#' @param subsample_control_size (integer) when the control sample size is very large, the methods take a very long time, and it would be better
#' to just subsample cells. e.g. 200k control cells available vs <1k treatment cells. you can subsample like 5000
#' set to NA to not subsample
#'
estimate_effects <- function(
                # test information
                AY_idx, AY, 
                # proximal number of NCs
                numNCs, 
                # data
                control_cell_idx, gene_norm, gene, grna, 
                libsize_log, 
                # info for NCs or covariates
                NCs_pca=NULL,
                validsinglegenes=NULL,
                covariates=NULL,
                # saving intermediateATEs
                save_intermediateATEs=TRUE, 
                skip_if_saved_intermediateATEs=FALSE,
                intermediateATEs_folder=NULL,
                # which proximal estimates
                run_proximal_linear_pca         = TRUE,
                run_proximal_negbin_pca         = TRUE,
                run_proximal_linear_singlegenes = TRUE,
                run_proximal_negbin_singlegenes = TRUE,
                subsample_control_size = NA) {
    # print(sprintf('ATE%s.csv', AY_idx))
    # print(list.files(intermediateATEs_folder))
    # skip this run if this test is already saved
    if(skip_if_saved_intermediateATEs & 
       (sprintf('ATE%s.csv', AY_idx) %in% list.files(intermediateATEs_folder))
       ) {

        res = read.csv( sprintf('%s/ATE%s.csv', intermediateATEs_folder, AY_idx) ) |> dplyr::select(-X) # removed saved row idx

        return(res)
    }
    res = NULL
    A_name = AY[AY_idx, 'grna']
    Y_name = AY[AY_idx, 'gene']
  
  
    
  
    # ===== assemble df =====
    # get idx of relevant grna
    if(grepl(pattern='Non_target', x = A_name)) {
      # 
      NT_vector = strsplit(A_name, '&')[[1]]#; print(NT_vector)
      NT_matrix = grna[NT_vector, ]
      # print(dim(grna[NT_vector, ]))
      NT_01 = colSums(grna[NT_vector, ]) |> as.vector()
      NT_idx = which(NT_01 == 1)
      
      control_cell_idx_noNTtrt = setdiff(selected_control_cells, NT_idx)
      if(is.na(subsample_control_size) | (subsample_control_size >= length(control_cell_idx))) {
        selected_control_cells = control_cell_idx_noNTtrt
      } else {
        selected_control_cells = sample(x=control_cell_idx_noNTtrt, size=subsample_control_size, replace=FALSE)
      }
      AY_data_idx = c(NT_idx, setdiff(selected_control_cells, NT_idx))

      A_01 = NT_01[AY_data_idx]

      rm(control_cell_idx_noNTtrt, selected_control_cells)
    } else {
      # use grna
      if(is.na(subsample_control_size) | (subsample_control_size >= length(control_cell_idx))) {
        selected_control_cells = control_cell_idx
      } else {
        selected_control_cells = sample(x=control_cell_idx, size=subsample_control_size, replace=FALSE)
      }
      
      AY_data_idx = c(which(grna[A_name, ]==1), selected_control_cells)

      A_01 = grna[A_name, AY_data_idx]

      rm(selected_control_cells)
    }

  # return(A_01)
    df = data.frame(
      A           =        A_01,
      Y           =        gene[Y_name, AY_data_idx], 
      Y_norm      =   gene_norm[Y_name, AY_data_idx],
      libsize_log = libsize_log[      AY_data_idx]
    )

  # return(df)
  
  
  lm_formula     = 'Y_norm ~ A'
  negbin_formula = 'Y_count ~ A + offset(libsize_log)'  
  if(!is.null(covariates)) {
    cur_covariates = covariates[cur_idx, ]
    colnames(cur_covariates) = paste0('X', 1:ncol(cur_covariates)) # rename columns to avoid potential naming mixups/collisions

    df = cbind(df, cur_covariates)

    cov_suffix = paste0(paste0('+', colnames(cur_covariates)), collapse = '') 
    lm_formula     = paste0(    lm_formula, cov_suffix)
    negbin_formula = paste0(negbin_formula, cov_suffix)
    rm(cov_suffix)
  } else {
    cur_covariates = NULL # set NULL for input for pci2s X
  }
  
  # return(df)
  
  # ===== estimate =====
  # ===== \__Linear Model Fit =====
  t0 = Sys.time()
  lmfit = lm(formula='Y_norm ~ A', data = df)
  lmfit_summary = summary(lmfit)
  t1 = Sys.time()
  
  res = bind_rows(res, 
                  data.frame(
                    NC_type     = NA,
                    method      = 'lm',
                    method_type = 'glm',
                    numNC       = NA,
                    ATE = lmfit_summary$coefficients['A', 'Estimate'],
                    se  = lmfit_summary$coefficients['A', 'Std. Error'],
                    tstat= lmfit_summary$coefficients['A', 't value'],
                    pval= lmfit_summary$coefficients['A', 'Pr(>|t|)'],
                    time_sec = difftime(t1, t0, units = 'secs') |> as.numeric()))
  rm(lmfit, lmfit_summary, t1, t0)
  # return(res)
  
  # ===== \__ Negative Binomial Fit  ====
  t0 = Sys.time()
  # nbfit = MASS::glm.nb(formula='Y ~ A + offset(libsize_log)', data = df)
  # nbfit_summary = summary(nbfit)
  nbfit = tryCatch({  # give negbin results if works
          MASS::glm.nb(formula='Y ~ A + offset(libsize_log)', data = df)
        },
        error = function(cond) { # give NA results if errored
          return(NULL)
        })
  t1 = Sys.time()
  if(is.null(nbfit)) {
    res = bind_rows(res, 
                  data.frame(
                    NC_type     = NA,
                    method      = 'negbin',
                    method_type = 'glm',
                    numNC       = NA,
                    ATE = NA,
                    se  = NA,
                    tstat= NA,
                    pval= NA,
                    time_sec = difftime(t1, t0, units = 'secs') |> as.numeric()))
  } else {
    nbfit_summary = summary(nbfit)
    res = bind_rows(res, 
                  data.frame(
                    NC_type     = NA,
                    method      = 'negbin',
                    method_type = 'glm',
                    numNC       = NA,
                    ATE = nbfit_summary$coefficients['A', 'Estimate'],
                    se  = nbfit_summary$coefficients['A', 'Std. Error'],
                    tstat= nbfit_summary$coefficients['A', 'z value'],
                    pval= nbfit_summary$coefficients['A', 'Pr(>|z|)'],
                    time_sec = difftime(t1, t0, units = 'secs') |> as.numeric()))
  }
  
  
  # res = bind_rows(res, 
  #                 data.frame(
  #                   NC_type     = NA,
  #                   method      = 'negbin',
  #                   method_type = 'glm',
  #                   numNC       = NA,
  #                   ATE = nbfit_summary$coefficients['A', 'Estimate'],
  #                   se  = nbfit_summary$coefficients['A', 'Std. Error'],
  #                   tstat= nbfit_summary$coefficients['A', 'z value'],
  #                   pval= nbfit_summary$coefficients['A', 'Pr(>|z|)'],
  #                   time_sec = difftime(t1, t0, units = 'secs') |> as.numeric()))
  rm(nbfit, nbfit_summary, t1, t0)
  # return(res)

  # ===== \__ Proximal PCA ====
  # Proximal with PCA loadings
  dfZ = NCs_pca[AY_data_idx, seq(from = 2, to = min(2*max(numNCs), ncol(NCs_pca)), by = 2)] # evens which((1:ncol(NCs)) %% 2 == 0)
  dfW = NCs_pca[AY_data_idx, seq(from = 1, to = min(2*max(numNCs), ncol(NCs_pca)), by = 2)] # odds
  colnames(dfZ) = paste0('Z', 1:ncol(dfZ))
  colnames(dfW) = paste0('W', 1:ncol(dfW))
  for(numNCs_ in numNCs)  {
    if(run_proximal_linear_pca) {
        # ===== \__\__Proximal Linear ====
        t0 = Sys.time()
        pci2s_res = tryCatch({  # give pci2s results if works
          pci2s::p2sls.lm(
            Y = df$Y_norm, 
            A = df$A, 
            W = dfZ[,1:numNCs_] |> as.data.frame(), 
            Z = dfW[,1:numNCs_] |> as.data.frame(), 
            X = cur_covariates,
            variance = TRUE)
        },
        error = function(cond) { # give NA results if errored
          return(list(summary_second_stage=
                        matrix(c(NA, NA, NA, NA), nrow=1, ncol=4, 
                               dimnames = list(c('A'), c('Estimate', 'Std. Error', 'z value', 'Pr(>|z|)')))))
        })
        t1 = Sys.time()
        
        res = bind_rows(res, 
                        data.frame(
                          NC_type     = 'PCA',
                          method      = 'proximallinear',
                          method_type = 'proximal',
                          numNC       = numNCs_,
                          ATE  = pci2s_res$summary_second_stage['A', 'Estimate'],
                          se   = pci2s_res$summary_second_stage['A', 'Std. Error'],
                          tstat= pci2s_res$summary_second_stage['A', 'z value'],
                          pval = pci2s_res$summary_second_stage['A', 'Pr(>|z|)'],
                          time_sec = difftime(t1, t0, units = 'secs') |> as.numeric()))
        
        rm(pci2s_res, t0, t1)
    }
    
    # ===== \__\__Proximal Negative Binomial ====
    if(run_proximal_negbin_pca){   
        t0 = Sys.time()
        pci2s_res = tryCatch({  # give pci2s results if works
          pci2s::p2sls.negbin(
            Y = df$Y, 
            A = df$A, 
            W = dfZ[,1:numNCs_] |> as.data.frame(), 
            Z = dfW[,1:numNCs_] |> as.data.frame(), 
            X = cur_covariates,
            offset = df$libsize_log,
            nco_type = rep("linear", numNCs_),
            # nco_args = lapply(X = 1:num_NCs, FUN = function(x){list(init=NA, offset=log(df_all$library_size))}), # No offsets for continuous versions
            variance = TRUE,
            verbose = FALSE)
        },
        error = function(cond) { # give NA results if errored
          return(list(summary_second_stage=
                        matrix(c(NA, NA, NA, NA), nrow=1, ncol=4, 
                               dimnames = list(c('A'), c('Estimate', 'Std. Error', 'z value', 'Pr(>|z|)')))))
        })
        t1 = Sys.time()
        
        res = bind_rows(res, 
                        data.frame(
                          NC_type     = 'PCA',
                          method      = 'proximalnegbin',
                          method_type = 'proximal',
                          numNC       = numNCs_,
                          ATE  = pci2s_res$summary_second_stage['A', 'Estimate'],
                          se   = pci2s_res$summary_second_stage['A', 'Std. Error'],
                          tstat= pci2s_res$summary_second_stage['A', 'z value'],
                          pval = pci2s_res$summary_second_stage['A', 'Pr(>|z|)'],
                          time_sec = difftime(t1, t0, units = 'secs') |> as.numeric()))
        
        rm(pci2s_res, t0, t1)
    }
    
  }
  rm(dfZ, dfW, numNCs_)
  
  
  
  # ===== \__ Proximal singlegenes ====
  if(run_proximal_linear_singlegenes | run_proximal_negbin_singlegenes) {
    if(is.null(validsinglegenes)) {
        print('Error with prox singlegenes: did not give validsinglegenes input')
    } else {
        # ===== \__\__Proximal Linear ====
        # Proximal with singlegenes
        NC_names = setdiff(validsinglegenes, c(A_name, Y_name)) # set of genes to use as NCs (not outcome and not grna target [assume grna name is the targeted gene name])
        NC_names = sample(NC_names, 2*max(numNCs))
        # Z_names = NC_names[1:max(numNCs)]
        # W_names = NC_names[(max(numNCs)+1):(2*max(numNCs))]
        
        NCs_singlegene = gene_norm[NC_names, AY_data_idx] |> t() |> as.data.frame()
        dfZ = NCs_singlegene[ , seq(from = 2, to = min(2*max(numNCs), ncol(NCs_singlegene)), by = 2)] # evens which((1:ncol(NCs)) %% 2 == 0)
        dfW = NCs_singlegene[ , seq(from = 1, to = min(2*max(numNCs), ncol(NCs_singlegene)), by = 2)] # odds
        colnames(dfZ) = paste0('Z', 1:ncol(dfZ))
        colnames(dfW) = paste0('W', 1:ncol(dfW))
        for(numNCs_ in numNCs)  {
            if(run_proximal_linear_singlegenes) {
                t0 = Sys.time()
                pci2s_res = tryCatch({  # give pci2s results if works
                pci2s::p2sls.lm(
                  Y = df$Y_norm, 
                  A = df$A, 
                  W = dfZ[,1:numNCs_], 
                  Z = dfW[,1:numNCs_], 
                  X = cur_covariates,
                  variance = TRUE)
                },
                error = function(cond) { # give NA results if errored
                return(list(summary_second_stage=
                              matrix(c(NA, NA, NA, NA), nrow=1, ncol=4, 
                                     dimnames = list(c('A'), c('Estimate', 'Std. Error', 'z value', 'Pr(>|z|)')))))
                })
                t1 = Sys.time()

                res = bind_rows(res, 
                              data.frame(
                                NC_type     = 'singlegene',
                                method      = 'proximallinear',
                                method_type = 'proximal',
                                numNC       = numNCs_,
                                ATE  = pci2s_res$summary_second_stage['A', 'Estimate'],
                                se   = pci2s_res$summary_second_stage['A', 'Std. Error'],
                                tstat= pci2s_res$summary_second_stage['A', 'z value'],
                                pval = pci2s_res$summary_second_stage['A', 'Pr(>|z|)'],
                                time_sec = difftime(t1, t0, units = 'secs') |> as.numeric()))

                rm(pci2s_res, t0, t1)
            }
              
          
            # ===== \__\__Proximal Negative Binomial ====
            if(run_proximal_negbin_singlegenes){ 
              t0 = Sys.time()
              pci2s_res = tryCatch({  # give pci2s results if works
                pci2s::p2sls.negbin(
                  Y = df$Y, 
                  A = df$A, 
                  W = dfZ[,1:numNCs_] |> as.data.frame(), 
                  Z = dfW[,1:numNCs_] |> as.data.frame(), 
                  X = cur_covariates,
                  offset = df$libsize_log,
                  nco_type = rep("linear", numNCs_),
                  # nco_args = lapply(X = 1:num_NCs, FUN = function(x){list(init=NA, offset=log(df_all$library_size))}), # No offsets for continuous versions
                  variance = TRUE,
                  verbose = FALSE)
              },
              error = function(cond) { # give NA results if errored
                return(list(summary_second_stage=
                              matrix(c(NA, NA, NA, NA), nrow=1, ncol=4, 
                                     dimnames = list(c('A'), c('Estimate', 'Std. Error', 'z value', 'Pr(>|z|)')))))
              })
              t1 = Sys.time()
              
              res = bind_rows(res, 
                              data.frame(
                                NC_type     = 'singlegene',
                                method      = 'proximalnegbin',
                                method_type = 'proximal',
                                numNC       = numNCs_,
                                ATE  = pci2s_res$summary_second_stage['A', 'Estimate'],
                                se   = pci2s_res$summary_second_stage['A', 'Std. Error'],
                                tstat= pci2s_res$summary_second_stage['A', 'z value'],
                                pval = pci2s_res$summary_second_stage['A', 'Pr(>|z|)'],
                                time_sec = difftime(t1, t0, units = 'secs') |> as.numeric()))
              
              rm(pci2s_res, t0, t1)
           }
        }
        rm(dfZ, dfW, numNCs_)
    } # end numNCs_
  }   # end singlegenes
  
  
  # ===== return ==== 
  res = res |> mutate(AY_idx = AY_idx, grna = A_name, gene = Y_name, type = AY[AY_idx, 'type'], 
                      n_trtmnt = sum(df$A==1), n_cntrl=sum(df$A==0), .before = 1)
  
  if(save_intermediateATEs) {
    write.csv(res, sprintf('%s/ATE%s.csv', intermediateATEs_folder, AY_idx))
  }
  return(res)
}






