


#' Quickly fit Linear Regression and a Negative Binomial fit for the
#' currently specified test (AY_idx) row of (AY_df). 
#' 
#' 
#' (modified from 2.2_jin_alldiscovery_estimates)
#'   control_cell_idx, gene_norm, gene, grna <-- move to arguments
#' @param AY_idx (integer) idx of the test in the df of all discovery pairs
#' @param AY_df (dataframe) of AY tests, must have columnss 'grna' and 'gene'
#' @param control_cell_idx (vector) of integer indices for the control cells
#'       (columns of gene_norm, gene, grna) and (rows of covariates)
#' @param gene_norm (matrix) of gene x cell continuous gene expression
#' @param gene (matrix) of gene x cell count gene expression
#' @param grna (matrix) of perturbation/grna x cell assignment
#' @param libsize_log (vector) of log library size (log nUMIs)
#' @param covariates (matrix/dataframe) of cell x features of potential 
#'                   covariates to add to the models
fit_glms <- function(AY_idx, AY_df, control_cell_idx, gene_norm, gene, grna, libsize_log, covariates=NULL) {
  # AY_idx = 2
  
  # Set up data: Y_count, Y_norm, A, libsize_log, (covariates (optional))
  cur_grna = AY_df[AY_idx, 'grna']
  cur_gene = AY_df[AY_idx, 'gene']
  
  cur_idx = c(which(grna[cur_grna, ]==1), control_cell_idx)
  df = data.frame(Y_count =      gene[cur_gene, cur_idx],
                  Y_norm  = gene_norm[cur_gene, cur_idx],
                  A       =      grna[cur_grna, cur_idx], 
                  libsize_log =     libsize_log[cur_idx])

  lm_formula     = 'Y_norm ~ A'
  negbin_formula = 'Y_count ~ A + offset(libsize_log)'  
  if(!is.null(covariates)) {
    cur_covariates = covariates[cur_idx, ]
    colnames(cur_covariates) = paste0('X', 1:ncol(cur_covariates)) # rename columns to avoid potential naming mixups/collisions

    df = cbind(df, cur_covariates)

    cov_suffix = paste0(paste0('+', colnames(cur_covariates)), collapse = '') 
    lm_formula     = paste0(    lm_formula, cov_suffix)
    negbin_formula = paste0(negbin_formula, cov_suffix)
    rm(cur_covariates, cov_suffix)
  } 
  
  # Lin Reg
  linreg_fit = lm(lm_formula, data = df)
  linreg_fit_summ = summary(linreg_fit)
  linreg_res = (data.frame(linreg_fit_summ$coefficients))[-1, ]
  colnames(linreg_res) = c('estimate', 'se', 'tstat', 'pval')
  row.names(linreg_res) = NULL
  
  # Neg Bin
  nbfit = MASS::glm.nb(formula=negbin_formula, data = df)
  nbfit_summary = summary(nbfit)
  nb_res = (data.frame(nbfit_summary$coefficients))[-1, ]
  colnames(nb_res) = c('estimate', 'se', 'tstat', 'pval')
  row.names(nb_res) = NULL
  

  # Return  
  return(rbind(linreg_res |> mutate(method =     'lm', grna=cur_grna, gene=cur_gene, n_trtmnt=sum(df$A==1), n_cntrl=sum(df$A==0), .before='estimate'), 
               nb_res     |> mutate(method = 'negbin', grna=cur_grna, gene=cur_gene, n_trtmnt=sum(df$A==1), n_cntrl=sum(df$A==0), .before='estimate')))
}




