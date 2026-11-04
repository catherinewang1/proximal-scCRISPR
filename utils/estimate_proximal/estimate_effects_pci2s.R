
#' Gather the used NC types from the specified estimator settings
#' used when subsetting NCs and useful for removing unnecessary NCs
gather_used_nctypes <- function(estimator_settings) {
  estimator_settings_nctypes = c() # all the NC_types used in the specified runs, will skip input NCs if none are used
  for(estimator_name in names(estimator_settings)) {
    if(estimator_settings[[estimator_name]][['RUN']] && !is.null(estimator_settings[[estimator_name]][['NC_type']])) {
      estimator_settings_nctypes = c(estimator_settings_nctypes, estimator_settings[[estimator_name]][['NC_type']])
    }
  }
  unique(estimator_settings_nctypes)
}



#' Estimate the effect specifying the 
#' - estimation method (model, model name, numNCs, etc...)
#' - and the relevant data: Y, A, X, Negative Controls (NCs)
#' 
#' Y: Outcome
#' A: Exposure/Treatment
#' X: Covariates (measured) 
#' Negative Controls: NCE and NCOs (ZWs)
#' @param estimator_name (character) of method name
#' GLMs will have special names 'lmYA', 'lmYAU', 'nbYA', 'nbYAU'.
#' Otherwise, will perform proximal and use this as the estimator name and use NC
#' @param Y (vector) of response  
#' @param A (vector) of treatment
#' @param X (matrix/dataframe) of confounders 
#' @param NC (list) of matrix/dataframe of negative controls. The names should match with the names of which_estimators
estimate_effect_one_estimator <- function(estimator_name, estimator_setting, Y, A,  X=NULL, NC=NULL, libsize_log=NULL) {
  
  # === === === === === === === === === ===
  # === Linear/GLM Models   === === === ===
  # === === === === === === === === === ===

  # print(sprintf('start estimate_effect_one_estimator: %s', estimator_name))

  # === \__Linear Model: Y ~ A (no confounder adj) ===
  if(estimator_name == 'lmYA') { # if specified and not already done
    t0 = Sys.time()
    lmfit = lm('Y ~ A', data.frame(Y=Y, A=A))
    lmfit_summary = summary(lmfit)
    t1 = Sys.time()
    
    
    return(data.frame(n_trtmnt = sum(A==1),
                      n_cntrl  = sum(A==0),
                      method      = 'lmYA',
                      method_type = 'lm',
                      NC_type     = NA,
                      numNC       = NA,
                      estimate    = lmfit_summary$coefficients['A', 'Estimate'],
                      se          = lmfit_summary$coefficients['A', 'Std. Error'],
                      tstat       = lmfit_summary$coefficients['A', 't value'],
                      pval        = lmfit_summary$coefficients['A', 'Pr(>|t|)'],
                      time_sec = difftime(t1, t0, units = 'secs') |> as.numeric()))
  } 
  
  # === \__Linear Model: Y ~ A + X1 + ... (w/ confounder adj) ===
  if(estimator_name == 'lmYAX') {
    # must give X to perf lmYAX 
    if(is.null(X)) {
      # print('skip lmYAU')
      return(data.frame(method      = 'lmYAX',
                        method_type = 'lm'))
    } else {
      # print('perform lmYAU')
      t0 = Sys.time()
      df = cbind(data.frame(Y=Y, A=A), X |> `colnames<-`(paste0('X', seq_len(ncol(X)))))
      lmfit = lm(paste0('Y~A+', paste0('X', seq_len(ncol(X)), collapse='+')), df)
      lmfit_summary = summary(lmfit)
      t1 = Sys.time()

      return(data.frame(n_trtmnt = sum(A==1),
                        n_cntrl  = sum(A==0),
                        method     = 'lmYAU',
                        method_type = 'lm',
                        numNC = NA,
                        estimate    = lmfit_summary$coefficients['A', 'Estimate'],
                        se          = lmfit_summary$coefficients['A', 'Std. Error'],
                        tstat       = lmfit_summary$coefficients['A', 't value'],
                        pval        = lmfit_summary$coefficients['A', 'Pr(>|t|)'],
                        time_sec = difftime(t1, t0, units = 'secs') |> as.numeric()))
    }
  }

  # === \__Negative Binomial: Y ~ A (no confounder adj) ===
  if(estimator_name == 'nbYA') {
    t0 = Sys.time()
    # nbfit = MASS::glm.nb(formula='Y ~ A + offset(libsize_log)', data = df)
    # nbfit_summary = summary(nbfit)
    nbfit = tryCatch({  # give negbin results if works
      MASS::glm.nb(formula='Y ~ A + offset(libsize_log)', data = data.frame(Y=Y, A=A, libsize_log=libsize_log))
    },
    error = function(cond) { # give NA results if errored
      return(NULL)
    })
    t1 = Sys.time()
    if(is.null(nbfit)) {
      return(data.frame(method      = 'nbYA',
                        method_type = 'nb',
                        time_sec = difftime(t1, t0, units = 'secs') |> as.numeric()))
    } else {
      nbfit_summary = summary(nbfit)
      return(data.frame(n_trtmnt = sum(A==1),
                        n_cntrl  = sum(A==0),
                        method      = 'nbYAU',
                        method_type = 'nb',
                        NC_type     = NA,
                        numNC       = NA,
                        estimate = nbfit_summary$coefficients['A', 'Estimate'],
                        se       = nbfit_summary$coefficients['A', 'Std. Error'],
                        tstat    = nbfit_summary$coefficients['A', 'z value'],
                        pval     = nbfit_summary$coefficients['A', 'Pr(>|z|)'],
                        time_sec = difftime(t1, t0, units = 'secs') |> as.numeric()))
    }
  }

  # === \__Negative Binomial: Y ~ A + log(libsize) + X1 + ... (with confounder adj) ===
  if(estimator_name == 'nbYAX') {
    # must give X to perf nbYAX 
    if(is.null(X)) {
      return(data.frame(method      = 'nbYAX',
                        method_type = 'nb'))
    } 
    t0 = Sys.time()
    
    nbfit = tryCatch({  # give negbin results if works
      MASS::glm.nb(formula=paste0('Y ~ A + offset(libsize_log) + ', paste0('X', seq_len(ncol(X)), collapse='+')), 
                   data  = cbind(data.frame(Y=Y, A=A, libsize_log=libsize_log), X |> `colnames<-`(paste0('X', seq_len(ncol(X))))))
    },
    error = function(cond) { # give NA results if errored
      return(NULL)
    })
    t1 = Sys.time()
    if(is.null(nbfit)) {
      return(data.frame(
              method      = 'nbYAX',
              method_type = 'nb',
              time_sec = difftime(t1, t0, units = 'secs') |> as.numeric())
              )
    } else {
      nbfit_summary = summary(nbfit)
      return(data.frame(
              n_trtmnt = sum(A==1),
              n_cntrl  = sum(A==0),
              method      = 'nbYAX',
              method_type = 'nb',
              NC_type     = NA,
              numNC       = NA,
              estimate = nbfit_summary$coefficients['A', 'Estimate'],
              se       = nbfit_summary$coefficients['A', 'Std. Error'],
              tstat    = nbfit_summary$coefficients['A', 'z value'],
              pval     = nbfit_summary$coefficients['A', 'Pr(>|z|)'],
              time_sec = difftime(t1, t0, units = 'secs') |> as.numeric())
            )
    }
  }
  
  # === === === === === === === === === ===
  # === Proximal Estimators === === === ===
  # === === === === === === === === === ===
  # stop early if settings/inputs don't align with proximal estimates
  if(is.null(estimator_setting$NUM_NCs) || is.null(NC) || !estimator_setting$IS_PROX) {return()} 
  

  res = NULL  
  # ===== \__for each num NC ====
  for(num_NC in estimator_setting$NUM_NCs) {
    # print(num_NCs)

    # ===== \__\__Proximal Linear ====
    if(estimator_setting$PROXIMAL_METHOD == 'proximallinear') {
      t0 = Sys.time()
      pci2s_res = tryCatch({  # give pci2s results if works
        pci2s::p2sls.lm(Y = Y, 
                        A = A, 
                        W = NC[, 1:num_NC*2 - 1], 
                        Z = NC[, 1:num_NC*2], 
                        variance = TRUE)
        },
        error = function(cond) {  NULL   }
      )
      t1 = Sys.time()
      if(is.null(pci2s_res)) {
        res = dplyr::bind_rows(res, 
                               data.frame(method     = 'proximallinear',
                                          method_type = 'proximal',
                                          time_sec    = difftime(t1, t0, units = 'secs') |> as.numeric()))
      } else {
        res = dplyr::bind_rows(res, 
                data.frame(
                  n_trtmnt = sum(A==1),
                  n_cntrl  = sum(A==0),
                  method      = 'proximallinear',
                  method_type = 'proximal',
                  NC_type     = estimator_setting$NC_type,
                  numNC       = num_NC,
                  estimate    = pci2s_res$summary_second_stage['A', 'Estimate'],
                  se          = pci2s_res$summary_second_stage['A', 'Std. Error'],
                  tstat       = pci2s_res$summary_second_stage['A', 'z value'],
                  pval        = pci2s_res$summary_second_stage['A', 'Pr(>|z|)'],
                  time_sec    = difftime(t1, t0, units = 'secs') |> as.numeric())
              )
      }
    }

    # ===== \__\__Proximal NegBin ====
    if(estimator_setting$PROXIMAL_METHOD == 'proximalnegbin') {
      t0 = Sys.time()
      pci2s_res = tryCatch({  # give pci2s results if works
        pci2s::p2sls.negbin(
              Y = Y, 
              A = A, 
              W = NC[, 1:num_NC*2 - 1], 
              Z = NC[, 1:num_NC*2    ], 
              X = X,
              offset = libsize_log,
              nco_type = rep("linear", num_NC),
              # nco_args = lapply(X = 1:num_NCs, FUN = function(x){list(init=NA, offset=log(df_all$library_size))}), # No offsets for continuous versions
              variance = TRUE,
              verbose = FALSE)
        },
        error = function(cond) {  NULL   }
      )
      t1 = Sys.time()
      if(is.null(pci2s_res)) {
        res = dplyr::bind_rows(res, 
                                data.frame(method     = 'proximalnegbin',
                                          method_type = 'proximal',
                                          time_sec    = difftime(t1, t0, units = 'secs') |> as.numeric()))
      } else {
        res = dplyr::bind_rows(res, 
                data.frame(
                  n_trtmnt = sum(A==1),
                  n_cntrl  = sum(A==0),
                  method      = 'proximalnegbin',
                  method_type = 'proximal',
                  NC_type     =  estimator_setting$NC_type,
                  numNC       = num_NC,
                  estimate    = pci2s_res$summary_second_stage['A', 'Estimate'],
                  se          = pci2s_res$summary_second_stage['A', 'Std. Error'],
                  tstat       = pci2s_res$summary_second_stage['A', 'z value'],
                  pval        = pci2s_res$summary_second_stage['A', 'Pr(>|z|)'],
                  time_sec    = difftime(t1, t0, units = 'secs') |> as.numeric())
              )
      }
    }
  }
  return(res)
    
}








#' Estimate the effect specifying the 
#' - estimation method (model, model name, numNCs, etc...)
#' - and the relevant data: Y, A, X, Negative Controls (NCs)
#' 
#' Y: Outcome
#' A: Exposure/Treatment
#' X: Covariates (measured) 
#' Negative Controls: NCE and NCOs (ZWs)
#' @param method (character) of method name ('lmYA', 'lmYAU', )
#' @param Y_norm (vector) of continuous response  
#' @param Y_count (vector) of count response
#' @param A (vector) of treatment
#' @param X (matrix/dataframe) of confounders 
#' @param NCs (list) of matrix/dataframe of negative controls. The names should match with the names of which_estimators
estimate_effect_one_test <- function(estimator_settings, Y_norm, Y_count, A, 
                                     X=NULL, NCs=NULL, libsize_log=NULL) { 
  res = data.frame( n_trtmnt = numeric(0),
                    n_cntrl  = numeric(0),
                    method      = character(0),
                    method_type = character(0),
                    NC_type  = character(0),
                    numNC    =   numeric(0),
                    estimate =   numeric(0),
                    se       =   numeric(0),
                    pval     =   numeric(0),
                    time_sec = numeric(0))
    
  for(estimator_name in names(estimator_settings)) {
    # print(paste0('estimate_effect_one_test: ', estimator_name))
    estimator_setting = estimator_settings[[estimator_name]]

    # skip if estimator is specified to not run
    if(is.null(estimator_setting$RUN) || !estimator_setting$RUN) { # element can exist but can also indicate false
      next
    } 
    
    # check that if method is a proximal method, there must be a NCs list value (to specify how NCs are gathered)
    if(!is.null(estimator_setting$IS_PROX) && estimator_setting$IS_PROX && !(estimator_setting$NC_type %in% names(NCs))) {
      print(sprintf('%s method must have a given NC', estimator_name))
      next
    }
    # print('set Y')
    if(!is.null(estimator_setting$Y_TYPE) && estimator_setting$Y_TYPE == 'continuous') {
      # print('set Y as continuous')
      Y = Y_norm
    } else if(!is.null(estimator_setting$Y_TYPE) && estimator_setting$Y_TYPE == 'count') {
      # print('set Y as count')
      Y = Y_count
    } else {
      print('in estimator_settings, specify the type of Y_TYPE as "count" or "continuous". default to continuous')
      Y = Y_norm
    }
    
    # print('estimate for one estimator')
    # estimate and bind to results
    res = rbind(res, 
                tryCatch(
                  expr = {    
                    estimate_effect_one_estimator(
                        estimator_name=estimator_name, 
                        estimator_setting=estimator_setting, 
                        Y=Y, 
                        A=A,
                        X=X, 
                        NC=NCs[[estimator_setting$NC_type]],
                        libsize_log=libsize_log)    
                  },
                  error = function(cond) {    NULL   }
                )
    )

  }
  
  return(res)

}

#' gather the relevant idx of cells for a certain A->() test.
#' will get the idx of cells recieivng A_name if present in grna. Or else with gather from the list of A_names.
#' Returns the list of relevant treated and then control idx
#' @param A_name (character)
#' @param grna (matrix)
#' @param control_cell_idx (vector) of control cells
#' @output list(  trt_idx = c(...),
#'              untrt_idx = c(...))
gather_AYidx <- function(A_name, grna, control_cell_idx, subsample_control_size=NA) {
  control_cell_idx_noNTtrt = control_cell_idx
  if(A_name %in% row.names(grna)) {
    A_idx = which(grna[A_name, ] == 1) |> unname()
  } else {
    # needs to be w pattern <name>&<name>... where <name> is a member of grna row.names
    # parse NT combo name (e.g. 'NT1&NT2')
    grnacombo_vector = strsplit(A_name, '&')[[1]]
    
    grnacombo_01 = colSums(grna[grnacombo_vector, , drop=FALSE]) |> as.vector()
    A_idx = which(grnacombo_01 == 1) |> unname()
    
    control_cell_idx_noNTtrt = setdiff(control_cell_idx_noNTtrt, A_idx) # remove these from control cell set
    
    rm(grnacombo_vector, grnacombo_01)
  }
  
  # subsample from controls if desired
  if(!is.na(subsample_control_size) && is.integer(subsample_control_size) && subsample_control_size < length(control_cell_idx_noNTtrt)) {
    control_cell_idx_noNTtrt = sample(control_cell_idx_noNTtrt, size=subsample_control_size)
  }
  return(list(  trt_idx = A_idx,
              untrt_idx = control_cell_idx_noNTtrt))
}

#' Estimate effects gene-by-gene (parallel). Good when there are many A's per Y response test 
#' 
#' When called from estimate_effects(...), pass in most arguments
#' 
estimate_effects_bygene <- function(estimator_settings, AY, NCs,
                                    grna, gene_norm, gene_counts, control_cell_idx, 
                                    libsize_log=NULL,
                                    covariates=NULL, 
                                    subsample_control_size=NA,
                                    save_intermediateATEs=TRUE, 
                                    skip_if_saved_intermediateATEs=FALSE,
                                    intermediateATEs_folder=NULL) {
  # all Y response names, and put positive gene responses first
  Y_names = unique(AY$gene)
  if("type" %in% colnames(AY)) {
    Y_pos = AY |> dplyr::filter(type == 'positive') |> pull(gene) |> unique()
    Y_names = c(Y_pos, setdiff(Y_names, Y_pos))
  }
  # for each Y response, estimate (A1->Y) (A2->Y) ... in parallel
  effects_df = NULL
  for(Y_name in Y_names) {
    # print(Y_name)
    intermediateATE_fn = sprintf('ATE_%s.csv', Y_name)
    if(skip_if_saved_intermediateATEs && !is.null(intermediateATEs_folder) && (intermediateATE_fn %in% list.files(intermediateATEs_folder))) {
      res_df = read.csv( sprintf('%s/%s', intermediateATEs_folder, intermediateATE_fn)) # removed saved row idx
    } else {
      A_names = AY |> filter(gene == Y_name) |> pull(grna)
      

      used_nctypes = gather_used_nctypes(estimator_settings=estimator_settings)

      # update the NCs to matrices if not already (e.g. randomly choose a subset of singelgenes to use as NCs)
      NCs_mat = list()
      for(estimator_name in names(estimator_settings)) {
        estimator_setting = estimator_settings[[estimator_name]]

        if(!estimator_setting$RUN) {   next   }
        if(is.null(estimator_setting$NC_type)) {   next   }
        if(!estimator_setting$NC_type %in% names(NCs)) {   next  } # should be specified, will err
        if(estimator_setting$NC_type %in% names(NCs_mat)) {   next   } # already subsetted

        if(estimator_setting$PROVIDES_NC_MATRIX) {
          # NCs gives already a (cell x #NC) matrix
          NCs_mat[[estimator_setting$NC_type]] = NCs[[estimator_setting$NC_type]]
        } else {
          # NCs gives a vector of potential singlegenes for these (A1, A2, ...) --> Y tests
          singlegenes_Y = setdiff(NCs[[estimator_setting$NC_type]], Y_name) # randomly sample some singlegenes to be used for these tests (smaller subset in order to reduce the size of objects distributed)
          NCs_mat[[estimator_setting$NC_type]] = gene_norm[sample(x = singlegenes_Y, 
                                                        size = min(length(singlegenes_Y), round(2*2*max(estimator_setting$NUM_NCs))), 
                                                        replace = FALSE),                                  ] |> t()
        }
      }      

      # estimate_bygene_Y_name <- function(A_name) {
      #   estimate_bygene(A_name,
      #                   grna    = grna, 
      #                   Y_count = gene_counts[Y_name, ],
      #                   Y_norm  =   gene_norm[Y_name], 
      #                   estimator_settings = estimator_settings,
      #                   NCs_mat            = NCs_mat,
      #                   control_cell_idx   = control_cell_idx,
      #                   libsize_log        = libsize_log)
      # }
      # # run in parallel (hopefully distributing fewer globals?)
      # res_par = future.apply::future_mapply(
      #     FUN      = estimate_bygene_Y_name,
      #     A_name   = A_names, 
      #     # future.globals = TRUE,
      #     future.seed = 123456, 
      #     SIMPLIFY = FALSE
      #   )
      # print('test one run')

      
      # return(estimate_effects_bygene_helper(
      #   A_name = A_names[1],
      #   grna    = grna,
      #   Y_count = gene_counts[Y_name, ],
      #   Y_norm  =   gene_norm[Y_name, ],
      #   estimator_settings = estimator_settings,
      #   NCs_mat            = NCs_mat,
      #   covariates         = covariates,
      #   control_cell_idx   = control_cell_idx,
      #   libsize_log        = libsize_log
      # ))

      # print('running parallel')
      # run in parallel (hopefully distributing fewer globals?)
      res_par = future.apply::future_mapply(
        FUN      = estimate_effects_bygene_helper,
        A_name   = A_names,
        MoreArgs = list(grna    = grna,
                        Y_count = gene_counts[Y_name, ],
                        Y_norm  =   gene_norm[Y_name, ],
                        estimator_settings = estimator_settings,
                        NCs_mat            = NCs_mat,
                        covariates         = covariates,
                        control_cell_idx   = control_cell_idx,
                        libsize_log        = libsize_log),
        # future.globals = TRUE,
        future.seed = 123456,
        SIMPLIFY = FALSE
      )
      # print('gather par results')
      # gather parallel results into a dataframe
      res_df = NULL
      for(i in 1:length(res_par)) { 
        cur_effects_par = res_par[[i]]
        if(!is.null(cur_effects_par)) {
          res_df = dplyr::bind_rows(res_df, cur_effects_par |> dplyr::mutate(gene=Y_name, grna=A_names[i], .before=1))
        }
        rm(cur_effects_par)
      }
      
      if(save_intermediateATEs && !is.null(intermediateATEs_folder)) {
        write.csv(res_df, sprintf('%s/%s', intermediateATEs_folder, intermediateATE_fn), row.names=FALSE)
      } 
    }
    effects_df = rbind(effects_df, res_df)
  }
  
  return(effects_df)
}


#' @param A_name (character) of treatment name. Should be a rowname of grna or a combination of rownames separated by '&'
#' @param grna (matrix) of perturbation/grna x cell assignment
#' @param Y_norm (vector) of continuous response  
#' @param Y_count (vector) of count response
#' @param estimator_settings (list) named list of settings specifying estimators 
#' @param NCs_mat (list) of matrices. All should be matrices now. The previous function call for singlegenes methods should
#'    have subset for the valid single genes to pass through, but also may provide a larger set of genes s.t. each estimate
#'    will be less likely to be affected by a 
estimate_effects_bygene_helper <- function(A_name, grna, Y_norm, Y_count,
                                           estimator_settings, NCs_mat, control_cell_idx, 
                                           covariates=NULL, libsize_log=NULL, subsample_control_size=NA) {
  tryCatch(
    expr = {
      # gather the relevant samples' (cells') idx: 
      AY_data_idx_list = gather_AYidx(A_name=A_name, grna=grna, control_cell_idx=control_cell_idx, subsample_control_size=subsample_control_size)

      AY_data_idx = c(AY_data_idx_list[['trt_idx']], AY_data_idx_list[['untrt_idx']])
      
      A10 = c(rep(1, times = length(AY_data_idx_list[['trt_idx'  ]])),
              rep(0, times = length(AY_data_idx_list[['untrt_idx']])))

      # subset within the NCs: e.g. subset for the relevant idx
      # current used cell idx for this test
      NCs_mat_AY = list()
      for(estimator_name in names(estimator_settings)) {
        estimator_setting = estimator_settings[[estimator_name]]

        if(!estimator_setting$RUN) {   next   }
        if(is.null(estimator_setting$NC_type)) {   next   }
        if(!estimator_setting$NC_type %in% names(NCs)) {   next  } # should be specified and given in NCs if a prox method, will err
        if(estimator_setting$NC_type %in% names(NCs_mat_AY)) {   next   } # already subsetted
        
        # all should be matrices now. The outside call should have sampled for a smaller subset of valid single genes 
        # (reduce size of objects spread to different threads)
        if(!is.null(estimator_setting$NC_RANDOMIZEORDER) && estimator_setting$NC_RANDOMIZEORDER) { # randomize the given NC order (good for singlegenes)
          NCs_mat_AY[[estimator_setting$NC_type]] = NCs_mat[[estimator_setting$NC_type]][AY_data_idx, sample(ncol(NCs_mat[[estimator_setting$NC_type]]), size=2*max(estimator_setting$NUM_NCs))]
        } else {
          NCs_mat_AY[[estimator_setting$NC_type]] = NCs_mat[[estimator_setting$NC_type]][AY_data_idx, ]
        }        
      } 
      
      Y_norm_AY  = if(is.null( Y_norm)) NULL else  Y_norm[AY_data_idx]
      Y_count_AY = if(is.null(Y_count)) NULL else Y_count[AY_data_idx]


      # print(Y_norm[AY_data_idx[1:4]])
      res = estimate_effect_one_test(
            estimator_settings = estimator_settings,
            Y_norm  = Y_norm_AY, 
            Y_count = Y_count_AY,
            A       = A10,
            X       = if(is.null(covariates)) NULL else X[AY_data_idx, ],
            NCs     = if(length(NCs_mat_AY) == 0) NULL else NCs_mat_AY, 
            libsize_log=if(is.null(libsize_log)) NULL else libsize_log[AY_data_idx]
          )
      
      # res = tryCatch(
      #   expr = {
      #     estimate_effect_one_test(
      #       estimator_settings = estimator_settings,
      #       Y_norm  = if(is.null( Y_norm)) NULL else  Y_norm[, AY_data_idx], 
      #       Y_count = if(is.null(Y_count)) NULL else Y_count[, AY_data_idx],
      #       A       = A10, 
      #       X       = if(is.null(covariates)) NULL else X[AY_data_idx, ],
      #       NCs     = if(length(NCs_mat_AY) == 0) NULL else NCs_mat_AY, 
      #       libsize_log=if(is.null(libsize_log)) NULL else libsize_log[AY_data_idx]
      #     )
      #   },
      #   error = function(cond) {  return(NULL)  }
      # )
      return(res)
    },
    error = function(cond) { return(NULL) }
  )
}






#' Estimate effects by AY (parallel). Good when there are few A's per Y response test
#' When called from estimate_effects(...), pass in most arguments
estimate_effects_byAY <- function(estimator_settings, AY, NCs,
                                    grna, gene_norm, gene_counts, control_cell_idx, 
                                    libsize_log=NULL,
                                    covariates=NULL, 
                                    subsample_control_size=NA,
                                    save_intermediateATEs=TRUE,
                                    skip_if_saved_intermediateATEs=FALSE,
                                    intermediateATEs_folder=NULL) {

  # create inner function to parallel call over
  inner <- function(AY_idx) {
    tryCatch(
      expr = {
        estimate_effects_byAY_helper(
            AY_idx=AY_idx, 
            estimator_settings=estimator_settings, AY=AY, NCs=NCs,
            grna=grna, 
            gene_norm=gene_norm, gene_counts=gene_counts, control_cell_idx=control_cell_idx, 
            libsize_log=libsize_log,
            covariates=covariates, 
            subsample_control_size=subsample_control_size,
            save_intermediateATEs=save_intermediateATEs, 
            skip_if_saved_intermediateATEs=skip_if_saved_intermediateATEs,
            intermediateATEs_folder=intermediateATEs_folder
          )
      },
      error = function(cond) {   return(NULL)   }
    )
  }

  # call over inner function   
  effects_par = future.apply::future_mapply(FUN = inner,
                                          AY_idx = 1:nrow(AY), 
                                          future.globals = TRUE,
                                          future.seed = 123456, 
                                          SIMPLIFY = FALSE)
  
  # assemble results
  effects_df = NULL
  for(i in 1:length(effects_par)) { 
    effects_df = dplyr::bind_rows(effects_df, effects_par[[i]]) # ok if NULL result
  }
  
  # return results
  return(effects_df)
}

#' Helper for estimate_effects_byAY which tests for each AY_idx test 
estimate_effects_byAY_helper <- function(AY_idx, 
                                    estimator_settings, AY, NCs,
                                    grna, gene_norm, gene_counts, control_cell_idx, 
                                    libsize_log=NULL,
                                    covariates=NULL, 
                                    subsample_control_size=NA,
                                    save_intermediateATEs=TRUE, 
                                    skip_if_saved_intermediateATEs=FALSE,
                                    intermediateATEs_folder=NULL) {
  
  intermediateATE_fn = sprintf('ATE_%s.csv', AY_idx)
  if(skip_if_saved_intermediateATEs && (intermediateATE_fn %in% list.files(intermediateATEs_folder))) {
    res_df = read.csv( sprintf('%s/%s', intermediateATEs_folder, intermediateATE_fn)) # removed saved row idx
    return(res_df)
  }

  A_name = AY[AY_idx, 'grna']
  Y_name = AY[AY_idx, 'gene']
  # gather the relevant samples' (cells') idx: 
  AY_data_idx_list = gather_AYidx(A_name=A_name, grna=grna, control_cell_idx=control_cell_idx, subsample_control_size=subsample_control_size)
  AY_data_idx = c(AY_data_idx_list[['trt_idx']], AY_data_idx_list[['untrt_idx']])
  A10 = c(rep(1, times = length(AY_data_idx_list[['trt_idx'  ]])),
          rep(0, times = length(AY_data_idx_list[['untrt_idx']])))
  
  
  # subset within the NCs: e.g. subset for the relevant idx
  # current used cell idx for this test
  NCs_mat_AY = list()
  for(estimator_name in names(estimator_settings)) {
    estimator_setting = estimator_settings[[estimator_name]]
    
    if(!estimator_setting$RUN) {   next   }
    if(is.null(estimator_setting$NC_type)) {   next   }
    if(!estimator_setting$NC_type %in% names(NCs)) {   next  } # should be specified and given in NCs if a prox method, will err
    if(estimator_setting$NC_type %in% names(NCs_mat_AY)) {   next   } # already subsetted


    if(!is.null(estimator_setting$PROVIDES_NC_MATRIX) && estimator_setting$PROVIDES_NC_MATRIX) {
      cur_NC_mat = NCs[[estimator_setting$NC_type]][AY_data_idx, ]
      if(!is.null(estimator_setting$NC_RANDOMIZEORDER) && estimator_setting$NC_RANDOMIZEORDER) {
        cur_NC_mat = cur_NC_mat[, sample(ncol(cur_NC_mat), size=min(ncol(cur_NC_mat), 2*max(estimator_setting$NUM_NCs)))]
      }   
      NCs_mat_AY[[estimator_setting$NC_type]] = cur_NC_mat[, 1:min(ncol(cur_NC_mat), 2*max(estimator_setting$NUM_NCs))]
      rm(cur_NC_mat)
    } else {
      # NCs gives a vector of potential singlegenes for these (A1, A2, ...) --> Y tests
      chosen_NCgenes = setdiff(NCs[[estimator_setting$NC_type]], Y_name) # randomly sample some singlegenes to be used for these tests (smaller subset in order to reduce the size of objects distributed)
      if(!is.null(estimator_setting$NC_RANDOMIZEORDER) && estimator_setting$NC_RANDOMIZEORDER) { # randomly take necessary NCs
        chosen_NCgenes = sample(x = chosen_NCgenes, size = min(length(chosen_NCgenes), 2*max(estimator_setting$NUM_NCs)))
      } else {  # choose first few NCs
        chosen_NCgenes = chosen_NCgenes[1:min(length(chosen_NCgenes), 2*max(estimator_setting$NUM_NCs))]
      }      
      NCs_mat_AY[[estimator_setting$NC_type]] = gene_norm[chosen_NCgenes, AY_data_idx] |> t() 
      rm(chosen_NCgenes)
    }      
  }
  
  
  res = tryCatch(
    expr = {
      estimate_effect_one_test(
          estimator_settings = estimator_settings, 
          Y_norm  = if(is.null(gene_norm  )) NULL else   gene_norm[Y_name, AY_data_idx], 
          Y_count = if(is.null(gene_counts)) NULL else gene_counts[Y_name, AY_data_idx],
          A       = A10, 
          X       = if(is.null(covariates)) NULL else covariates[AY_data_idx, ],
          NCs     = if(length(NCs_mat_AY) == 0) NULL else NCs_mat_AY, 
          libsize_log=if(is.null(libsize_log)) NULL else libsize_log[AY_data_idx]
        )
    },
    error = function(cond) { return(NULL) }
  )
   
  res = res |> dplyr::mutate(grna = A_name, gene = Y_name, .before=1)

  if(save_intermediateATEs && !is.null(intermediateATEs_folder) && dir.exists(intermediateATEs_folder)) {
    write.csv(res, sprintf('%s/%s', intermediateATEs_folder, intermediateATE_fn), row.names=FALSE)
  } 
  return(res)
}


#' Main function to call to estimate effects using the specified estimators with estimator settings and the specified tests
#' 
#' @param estimator_settings (list) named list of settings specifying estimators desired to be run and the relevant estimators' parameters
#'      such as NUM_NC, ...
#' Each element must have the following named elements:
#'      RUN=bool whether to run this method or not
#'      Y_TYPE \in {'continuous', 'count'}
#' and proximal methods should have more, such as IS_PROX, NUM_NC, ...
#' @param AY (dataframe) of AY tests, must have columnss 'grna' and 'gene'
#' @param NCs (list) of matrices/dataframes of NC variables (e.g. PCA loadings)
#'            or of gene names (for singlegene NCs)
#' @param grna (matrix) of perturbation/grna x cell assignment
#' @param gene_norm (matrix) of gene x cell continuous gene expression
#' @param gene_counts (matrix) of gene x cell count gene expression
#' @param control_cell_idx (vector) of integer indices for the control cells
#'       (columns of gene_norm, gene, grna) and (rows of covariates)
#' @param parallel_by_gene (boolena) whether to run parallel by gene (rather than by each AY row). This will be helpful when there are many tests per Y response
#' @param libsize_log (vector) of log library size (log nUMIs)
#' @param covariates (matrix/dataframe) of cell x features of potential 
#'                   covariates to add to the models
#' @param subsample_control_size (integer) when the control sample size is very large, the methods take a very long time, and it would be better
#' @param save_intermediateATEs (boolean) whether to save this current AY test in intermediateATEs/ subfolder
#' @param skip_if_saved_intermediateATEs (boolean) whether to skip computation if this AY test is already saved in intermediateATEs/ subfolder
estimate_effects <- function(estimator_settings, AY, NCs,
                                     grna, gene_norm, gene_counts, control_cell_idx, 
                                     parallel_by_gene = FALSE,
                                     libsize_log=NULL,
                                     covariates=NULL, 
                                     subsample_control_size=NA,
                                     save_intermediateATEs=TRUE, 
                                     skip_if_saved_intermediateATEs=FALSE,
                                     intermediateATEs_folder=NULL) {

  if(save_intermediateATEs && !is.null(intermediateATEs_folder)) {    dir.create(intermediateATEs_folder, showWarnings = FALSE)    }
  
  # detect if the tests are dense or not (more dense = more worth it to perf parallel by gene)
  if(parallel_by_gene == 'auto') {
    median_tests_per_Y = AY |> dplyr::group_by(gene) |> dplyr::summarize(count = n()) |> pull(count) |> median()
    parallel_by_gene = (median_tests_per_Y >= 10)
    rm(median_tests_per_Y)
  }

  
  if(parallel_by_gene) { # by each response Y
    return(estimate_effects_bygene(...))
  } else { # by each AY row
    return(estimate_effects_byAY(...))
  }
}


      


      