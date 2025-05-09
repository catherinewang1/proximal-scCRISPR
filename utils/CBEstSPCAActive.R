# Functions to call for the active procedure 
# Probably the easiest way without rewriting new code, is to call the previous fns

# requires loading in functions in file `CBEstAllSPCA.R
# ie. get_ATE_est_NCs_make


require(lpdensity) # estimation of density for indep p-value active

# =================== Independent p-values ============================================

#' @param lpdensEstimate (matrix) result from lpdens call with (result)$Estimate
#' @examples 
#' est1 = lpdensity::lpdensity(data = dfsim_null0$x, bwselect = "imse-dpi", grid = seq(from = 0, to = 1, length.out = 1000))
#' myDensEstimate <- makeGetDensityEstimate(est1$Estimate)
makeGetDensityEstimate <- function(lpdensEstimate) {
  #' to get estimates, we just take the closest grid point (close enough hopefully)
  #' @param x (numeric) between 0,1 to get the estimated density at f(x)
  getDensityEstimate <- function(x) {
    lpdensEstimate[base::findInterval(x, lpdensEstimate[,'grid']), 'f_p']
  }
  
  return(getDensityEstimate)
}



#' Set c parameter. Then creates and 
#' outputs a function that takes in density qx value
#' and returns the probability of auditing, ax
#' (From eqn, choose c = -L where L is the lower bound for qx)
#' (e.g. qx \in [.05, \infty) --> use c = -.05)
#' @example 
#' ax01 = create_ax(-.01) # L = .01
#' qx = seq(from = .001, to = 5, length.out = 200)
#' plot(qx, ax01(qx), type = 'l')
create_ax <- function(c) {
  ax <- function(qx) {   c / qx + 1   }
  ax
}


#' Calculate the active pvalue in the arbitrarily dependent case
#' using saved proxy and true p-values (e.g. when both are already run, randomly
#' decide to use true pval with some prob, dep on parameter gamma)
#' @param pvals_df (dataframe) with at least the columns
#'      pval_proxy (numeric) proxy pvalue
#'      pval_true  (numeric)  true pvalue
#' @param pval_proxy_null (vector) of pvalues using proxy method on null hypothesis tests
#'                                 density of the null distn will be estimated on this
get_active_indep_pval_using_saved <- function(pvals_df, pval_proxy_null) {
        
    # estimate density of null proxy
    lpdensResult = lpdensity::lpdensity(data = pval_proxy_null, bwselect = "imse-dpi", grid = seq(from = 0, to = 1, length.out = 1000))
    myDensEstimate = makeGetDensityEstimate(lpdensResult$Estimate)
    
    # estimated densities
    estqx = myDensEstimate(pvals_df$pval_proxy)
    # audit/query probabilities
    get_query_prob = create_ax(-min(estqx, na.rm = TRUE)) # L = min of estimated dens
    query_prob = get_query_prob(estqx)

    T_ = rbinom(n = nrow(pvals_df), size = 1, prob = query_prob)
    pval_active_indep  = (1 - T_) * pvals_df$pval_proxy + (T_) * pvals_df$pval_true
    
    # add active time info if proxy and true time info is provided
    time_info_provided = all(c('time_proxy', 'time_true') %in% colnames(pvals_df))
    if(time_info_provided) {
        # actual time   (using if actually queried)
        time_active         = pvals_df$time_proxy +         T_ * pvals_df$time_true
        # expected time (usimg probability of query) 
        time_active_expect  = pvals_df$time_proxy + query_prob * pvals_df$time_true
    } else {
        time_active        = NA
        time_active_expect = NA
    }
    
    return(data.frame(pval_active_indep = pval_active_indep,
                      time_active_indep = time_active,
                      time_active_indep_expect = time_active_expect))
}







# =================== Arbitrary Dependence p-values ===================================


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
        # print(sprintf("[%s] get_proxy_pval called", Sys.time()))
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
        # print(sprintf("[%s] get_true_pval called", Sys.time()))
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
#' @return (named list) active pval + other info
get_active_arbdep_pval_make <- function(get_proxy_pval, get_true_pval, gamma = .5) {

    get_active_arbdep_pval <- function(AY_idx) {
        t0 = Sys.time()
        
        # first call proxy
        t0_proxy   = Sys.time()
        pval_proxy = get_proxy_pval(AY_idx)
        t_proxy    = difftime(Sys.time(), t0_proxy, units = 'secs') |> as.numeric()

        # then maybe call true
        T_ = rbinom(n=1, size=1, prob= 1 - gamma * pval_proxy)
        # (1-T) Q + T (1-gamma)^{-1} P
        if(T_ > .5) {
            # print(sprintf("[%s] calling true pvalue fn", Sys.time()))
            t0_true   = Sys.time()
            pval_true = get_true_pval(AY_idx)
            t_true    = difftime(Sys.time(), t0_true, units = 'secs') |> as.numeric()
            # print(sprintf("[%s]     true pval = %.2f", Sys.time(), pval_true))

            pval_active = ( 1/(1-gamma)) * pval_true
        } else {  
            pval_true   = NA
            t_true      = NA

            pval_active = pval_proxy
        }
        
        t1 = Sys.time()
        
        res = list(AY_idx       = AY_idx,
                   queried_true = T_,
                   pval_active  = pval_active,
                   pval_proxy   = pval_proxy, # redundant info but ok
                   pval_true    = pval_true,
                   gamma        = gamma, 
                   time_active  = difftime(t1, t0, units = 'secs') |> as.numeric(),
                   time_proxy   = t_proxy,
                   time_true    = t_true)
        return(res)
    }

    # return the inner function (which takes in AY idx, and returns the active p-value)
    return(get_active_arbdep_pval)
}




#' Calculate the active pvalue in the arbitrarily dependent case
#' using saved proxy and true p-values (e.g. when both are already run, randomly
#' decide to use true pval with some prob, dep on parameter gamma)
#'
#' For proxy Q, true P, and random T|Q \sim Bern(1-gamma Q) for gamma \in (0,1]:
#'  (1-T) Q + T (1-gamma)^{-1} P
#' @param pvals_df (dataframe) with at least the columns
#'      AY_idx     (numeric) AY test index
#'      pval_proxy (numeric) proxy pvalue
#'      pval_true  (numeric)  true pvalue
#' @param gamma (numeric) \in (0, 1]
#'
#' @return (named list) active pval + other info
get_active_arbdep_pval_using_saved_make <- function(pvals_df, gamma = .5) {
    
    # if proxy and true time info is provided
    time_info_provided = all(c('time_proxy', 'time_true') %in% colnames(pvals_df))

    get_active_arbdep_pval_using_saved <- function(AY_idx) {
        
        AY_row     = pvals_df[(pvals_df[, 'AY_idx'] == AY_idx), ]
        pval_proxy = AY_row[1, 'pval_proxy']
        pval_true  = AY_row[1, 'pval_true']

        # maybe call true
        queried_true_prob = 1 - gamma * pval_proxy
        T_ = rbinom(n=1, size=1, prob= queried_true_prob)


        # (1-T) Q + T (1-gamma)^{-1} P
        # if(T_ > .5) {
        #     pval_true = get_true_pval(AY_idx)
        #     pval_active = ( 1/(1-gamma)) * pval_true
        # } else {  
        #     pval_true   = NA
        #     pval_active = pval_proxy
        # }
        pval_active = (1-T_) * pval_proxy + T_ * (1-gamma)^{-1} * pval_true

        # add active time info if proxy and true time info is provided
        if(time_info_provided) {
            # actual time (using if queried)
            time_active        = AY_row[1, 'time_proxy'] +                T_ * AY_row[1, 'time_true']
            # expected time (usimg probability of query)
            time_active_expect = AY_row[1, 'time_proxy'] + queried_true_prob * AY_row[1, 'time_true']
        } else {
            time_active        = NA
            time_active_expect = NA
        }
        
        res = list(AY_idx       = AY_idx,
                   queried_true = T_,
                   queried_true_prob = queried_true_prob,
                   pval_active  = pval_active,
                   pval_proxy   = pval_proxy, # redundant info but ok
                   pval_true    = pval_true,
                   gamma        = gamma, 
                   time_active  = time_active,
                   time_active_expect = time_active_expect)
        return(res)
    }
    
    return(get_active_arbdep_pval_using_saved)
}


