# Functions to call for the active procedure 
# Probably the easiest way without rewriting new code, is to call the previous fns

# requires loading in functions in file `CBEstAllSPCA.R
# ie. get_ATE_est_NCs_make





#' Make fn that gets the PROXY / BIASED p-value 
#' In this case, this is the lmYA p-value
#' from simple linear regression of Y on A
#' @param (all the params for get_ATE_est_NCs_make, except for which_estimators)
#' @returns function- that takes in AY idx, and returns the lmYA pval 
get_proxy_pval_make <- function(AY, gene_norm, NCs,
                                 grna_rownames, grna, 
                                 NT_idx, imp_gene_names,
                                 # which_estimators, 
                                 CB_setting, save_path=NULL) {

    # create a list where all methods are FALSE, except for lmYA
    which_estimators_onlylmYA = list(lm_YA        = TRUE,
                                     OCB_2SLS     = FALSE,
                                     OCB_2SLS_pci2s=FALSE,
                                     OCB_2SLSReg  = FALSE,
                                     OCB_GMM      = FALSE,
                                     OCB_GMMRw    = FALSE,
                                     OCB_GMMRwReg = FALSE,
                                     OCB_LinOSPI  = FALSE,
                                     OCB_LinOS    = FALSE,
                                     OCB_LinOStrim= FALSE)

    # make the get_ATE_est_NCs_make function using the new which_estimators, 
    # and inputting all the other provided params
    get_ATE_est_NCs_lmYA = get_ATE_est_NCs_make(AY=AY, 
                                            gene_norm=gene_norm,
                                            NCs=NCs,
                                            grna_rownames=grna_rownames, 
                                            grna=grna, 
                                            NT_idx=NT_idx, 
                                            imp_gene_names=imp_gene_names, 
                                            which_estimators=which_estimators_onlylmYA, # <-- changed
                                            CB_setting=CB_setting, 
                                            save_path=save_path)
    
    #' takes in AY idx and returns lmYA p-value
    #' @param idx (integer) index of the AY pair in the list
    #'                      of AY pairs
    get_proxy_pval <- function(idx) {
        # call fn to analyze
        res = get_ATE_est_NCs_lmYA(idx)
        # extract p-value (should be a df with 1 row?)
        return(res[1, "pval"])
    }
    
    # return the inner function (which takes in AY idx, and returns the lmYA p-value)
    return(get_proxy_pval)
} 


#' Make fn that gets the TRUE / UNBIASED p-value 
#' In this case, this is the proximal p-value
#'
#' @param (all the params for get_ATE_est_NCs_make, except for which_estimators)
#' @returns function- that takes in AY idx, and returns the OCB_2SLS_pci2s pval 
get_true_pval_make <- function(AY, gene_norm, NCs,
                                 grna_rownames, grna, 
                                 NT_idx, imp_gene_names,
                                 # which_estimators, 
                                 CB_setting, save_path=NULL) {

    # create a list where all methods are FALSE, except for OCB_2SLS_pci2s
    which_estimators_onlypci2s = list(lm_YA        = FALSE,
                                     OCB_2SLS     = FALSE,
                                     OCB_2SLS_pci2s=TRUE,
                                     OCB_2SLSReg  = FALSE,
                                     OCB_GMM      = FALSE,
                                     OCB_GMMRw    = FALSE,
                                     OCB_GMMRwReg = FALSE,
                                     OCB_LinOSPI  = FALSE,
                                     OCB_LinOS    = FALSE,
                                     OCB_LinOStrim= FALSE)

    # make the get_ATE_est_NCs_make function using the new which_estimators, 
    # and inputting all the other provided params
    get_ATE_est_NCs_pci2s = get_ATE_est_NCs_make(AY=AY, 
                                            gene_norm=gene_norm,
                                            NCs=NCs,
                                            grna_rownames=grna_rownames, 
                                            grna=grna, 
                                            NT_idx=NT_idx, 
                                            imp_gene_names=imp_gene_names, 
                                            which_estimators=which_estimators_onlypci2s, # <-- changed
                                            CB_setting=CB_setting, 
                                            save_path=save_path)
    
    #' takes in AY idx and returns OCB_2SLS_pci2s p-value
    #' @param idx (integer) index of the AY pair in the list
    #'                      of AY pairs
    get_true_pval <- function(idx) {
        # call fn to analyze
        res = get_ATE_est_NCs_pci2s(idx)
        # extract p-value (should be a df with 1 row?)
        return(res[1, "pval"])
    }
    
    # return the inner function (which takes in AY idx, and returns the OCB_2SLS_pci2s p-value)
    return(get_true_pval)
} 




#' Calculate the active pvalue in the arbitrarily dependent case
#' For proxy Q, true P, and random T|Q \sim Bern(1-gamma Q) for gamma \in (0,1]:
#'  (1-T) Q + T (1-gamma)^{-1} P
#'
#' @param get_proxy_pval (function) that takes input AY idx, and returns proxy pvalue
#' @param get_true_pval  (function) that takes input AY idx, and returns true pvalue
#' @param gamma (numeric) \in (0, 1]
#'
#' (prob also return some other values, e.g. if true pval was called, computation time, etc...)
#' @return active pval
get_active_arbdep_pval_make <- function(get_proxy_pval, get_true_pval, gamma = .5) {

    get_active_arbdep_pval <- function(AY_idx) {
        t0 = Sys.time()

        # first call proxy
        pval_proxy = get_proxy_pval(AY_idx)
        # then maybe call true
        T_ = rbinom(n=1, size=1, prob= 1 - gamma * pval_proxy)
        # (1-T) Q + T (1-gamma)^{-1} P
        if(T_ > .5) {
            pval_true = get_true_pval(AY_idx)
            active_pval = ( 1/(1-gamma)) * pval_true
        } else {  
            pval_true = NA
            active_pval = pval_proxy
        }

        t1 = Sys.time()

        res = list(AY_idx = AY_idx,
                   pval = active_pval,
                   queried_true = T_,
                   pval_proxy = pval_proxy, # redundant info but ok
                   pval_true  = pval_true,
                   gamma = gamma, 
                   time = difftime(t1, t0, units = 'secs') |> as.numeric())
        return(res)
    }

    # return the inner function (which takes in AY idx, and returns the active p-value)
    return(get_active_arbdep_pval)
}


