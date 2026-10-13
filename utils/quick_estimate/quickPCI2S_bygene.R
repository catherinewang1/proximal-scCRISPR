




#' estimate ATEs/logFoldChanges/Coefficients by gene
#' By gene means: input the specified gene (response), and all the relevant perturbations (treatments) are
#' then gathered from the df of AY tests. These tests are all performed at once. This is hopefully faster
#' because I think the computation is slow because of the query/distribution of the large gene_norm/gene_counts 
#' datasets.
#' 
#' 
#' res is dataframe with cols:
#' - grna, gene, type 
#' - ... etc ...
#' - tstat: test statistic either the t value (lm) or z value (glm, proximal)
#' 
#' For inputs, changed to 
#' Y_info = list(Y_name=character, 
#'               Y_norm=numeric vector,
#'               Y_count=integer vector)
#' @param Y_info (integer) index of the AY test 
#' @param AY (dataframe) of AY tests, must have columnss 'grna' and 'gene'
#' @param numNCs (vector) of number of Negative Control pairs to be used for proximal estimates
#' @param control_cell_idx (vector) of integer indices for the control cells
#'       (columns of gene_norm, gene, grna) and (rows of covariates)
#' @param gene_norm (matrix) of gene x cell continuous gene expression
#' @param gene_counts (matrix) of gene x cell count gene expression
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
estimate_effects_bygene <- function(
                # test information
                Y_info, 
                AY, 
                # proximal number of NCs
                numNCs, 
                # data
                control_cell_idx, 
                # gene_norm, gene_counts, 
                grna, 
                libsize_log, 
                # info for NCs or covariates
                NCs_pca=NULL,
                Y_validsinglegenes=NULL,
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

  intermediateATE_fn = sprintf('ATE_%s.csv', Y_info$Y_name)
  if(skip_if_saved_intermediateATEs & (intermediateATE_fn %in% list.files(intermediateATEs_folder))) {
    res = read.csv( sprintf('%s/%s', intermediateATEs_folder, intermediateATE_fn)) # removed saved row idx
    return(res)
  }

  A_names = AY |> filter(gene == Y_info$Y_name) |> pull(grna)
  # map of A_name to the cell indicies that received this perturbation (or combo of perturbations)
  A_map_idx = list()
  control_cell_idx_noNTtrt = control_cell_idx
  for(A_name in A_names) {
    if(A_name %in% row.names(grna)) {
      A_map_idx[[A_name]]= which(grna[A_name, ] == 1) |> unname()
    } else {
      # needs to be w pattern <name>&<name>... where <name> is a member of grna row.names
      # parse NT combo name (e.g. 'NT1&NT2')
      grnacombo_vector = strsplit(A_name, '&')[[1]]
      
      grnacombo_01 = colSums(grna[grnacombo_vector, , drop=FALSE]) |> as.vector()
      grnacombo_idx = which(grnacombo_01 == 1) |> unname()
      
      A_map_idx[[A_name]] = grnacombo_idx
      
      control_cell_idx_noNTtrt = setdiff(control_cell_idx_noNTtrt, grnacombo_idx) # remove these from control cell set
      
      rm(grnacombo_vector, grnacombo_01, grnacombo_idx)
    }
  }

  

  # create a function to estimate per specified treatment 
  inner_estimate_effect <- function(A_name) {
    AY_idx = which((AY$grna == A_name) & (AY$gene == Y_name))    
    res = NULL
    
    # current used cell idx for this test
    cur_idx = c(control_cell_idx_noNTtrt, A_map_idx[[A_name]])
    
    df = data.frame(
      A = c(rep(0, times = length(control_cell_idx_noNTtrt)), 
            rep(1, times = length(A_map_idx[[A_name]]))),
      Y_count     = Y_info$Y_count[cur_idx], 
      Y_norm      = Y_info$Y_norm[cur_idx],
      libsize_log = libsize_log[cur_idx], 
      check.names=FALSE
    )
    
    # add covariates to 
    lm_formula     = 'Y_norm ~ A'
    negbin_formula = 'Y_count ~ A + offset(libsize_log)'
    if(!is.null(covariates)) {
      cur_covariates = covariates[cur_idx, ]
      colnames(cur_covariates) = paste0('X', 1:ncol(cur_covariates)) # rename columns to avoid potential naming mixups/collisions
      
      # append to df
      df = cbind(df, cur_covariates)
      # append to glm formulas
      cov_suffix = paste0(paste0('+', colnames(cur_covariates)), collapse = '') 
      lm_formula     = paste0(    lm_formula, cov_suffix)
      negbin_formula = paste0(negbin_formula, cov_suffix)
      rm(cov_suffix)
    } else {
      cur_covariates = NULL # set NULL for input for pci2s X
    }
    
    
    # ===== estimate =====
    # ===== \__Linear Model Fit =====
    t0 = Sys.time()
    lmfit = lm(formula=lm_formula, data = df)
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
      MASS::glm.nb(formula=negbin_formula, data = df)
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
      rm(nbfit_summary)
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
    rm(nbfit, t1, t0)
    
    # ===== \__ Proximal PCA ====
    # Proximal with PCA loadings
    dfZ = NCs_pca[cur_idx, seq(from = 2, to = min(2*max(numNCs), ncol(NCs_pca)), by = 2)] # evens which((1:ncol(NCs)) %% 2 == 0)
    dfW = NCs_pca[cur_idx, seq(from = 1, to = min(2*max(numNCs), ncol(NCs_pca)), by = 2)] # odds
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
        # return(df)
        t0 = Sys.time()
        pci2s_res = tryCatch({  # give pci2s results if works
          pci2s::p2sls.negbin(
            Y = df$Y_count, 
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
      if(is.null(Y_validsinglegenes)) {
        print('Error with prox singlegenes: did not give Y_validsinglegenes input')
      } else {
        # print('run_proximal_linear_singlegenes')
        # ===== \__\__Proximal Linear ====
        NCs_singlegene = Y_validsinglegenes[sample(nrow(Y_validsinglegenes), 2*max(numNCs)), cur_idx] |> t() |> as.data.frame()
        dfZ = NCs_singlegene[ , seq(from = 2, to = min(2*max(numNCs), ncol(NCs_singlegene)), by = 2)] # evens which((1:ncol(NCs)) %% 2 == 0)
        dfW = NCs_singlegene[ , seq(from = 1, to = min(2*max(numNCs), ncol(NCs_singlegene)), by = 2)] # odds
        colnames(dfZ) = paste0('Z', 1:ncol(dfZ))
        colnames(dfW) = paste0('W', 1:ncol(dfW))
        rm(NCs_singlegene)
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
                Y = df$Y_count, 
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
    res = res |> mutate(grna = A_name, gene = Y_name, 
                        # type = AY |> filter(grna == A_name & gene == Y_name) |> head(1) |> pull(type), 
                        type = AY[AY_idx, 'type'],
                        n_trtmnt = sum(df$A==1), n_cntrl=sum(df$A==0), .before = 1)
    
    return(res)
  }

  # return(inner_estimate_effect)

  # run in parallel (hopefully distributing fewer globals?)
  res_par = future.apply::future_mapply(inner_estimate_effect,
                                        A_name = A_names, 
                                        # future.globals = TRUE,
                                        future.seed = 123456, 
                                        SIMPLIFY = FALSE)

  # gather parallel results into a dataframe
  res_df = NULL
  for(i in 1:length(res_par)) { 
    cur_effects_par = res_par[[i]]
    if(!is.null(cur_effects_par)) {
      res_df = dplyr::bind_rows(res_df, cur_effects_par)
    }
    rm(cur_effects_par)
  }
  
  if(save_intermediateATEs) {
    write.csv(res, sprintf('%s/%s', intermediateATEs_folder, intermediateATE_fn), row.names=FALSE)
  } 

  return(res_df)

}






