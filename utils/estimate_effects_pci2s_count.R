


if(F) {
NC_name_raw = 'simpleCount-proximalCountCountCount'
NC_names = strsplit(x = NC_name_raw, split = '-') |> unlist()
proximal_settings[NC_names]
NCs_settings            = proximal_settings[NC_names]
NC_name = NC_names[2]

# === only extract cols and idx of gene and grna we are interested in (to hopefully save time)

# create count NCE and count NCO
NCE_type='count'; NCO_type='count' # NCE_type='continuous'; NCO_type='continuous'
if(NCE_type=='count' & NCO_type=='count') {
  AYZW = readRDS(sprintf('%s/AY/%s/AYZW.rds', save_dir, AYZW_setting_name))
  max_num_NC_pairs = proximal_settings[['proximalCountCountCount']]$num_NC_pairs |> max() # max number of ZW used 
  for(A_name in names(AYZW)) {
    for(Y_name in names(AYZW[[A_name]])) {
      all_Ys = c(all_Ys, 
                 AYZW[[A_name]][[Y_name]][[1]]$Z_names[1:max_num_NC_pairs] ,
                 AYZW[[A_name]][[Y_name]][[1]]$W_names[1:max_num_NC_pairs] ) |> unique()
      
    }
  }
  rm(A_name, Y_name, AYZW, max_num_NC_pairs)
}







# === TESTING PCI2S Negative Binomial

N <- 2000
expit  <-  function(x) exp(x)/(1 + exp(x))
U <- runif(N); X <- runif(N)
A <- rbinom(N, 1, expit(-3 + 5 * U + 1 * X))
Y <- rpois(N, exp(0.2 + 1.5 * U + 0.2 * X + 0.2 * A))
Z <- rnorm(N, 2 * U + 0.5 * X)
W2 <- rexp(N, 0.1 + 1 * U + 0.1 * X)
D2 <- as.numeric(W2 < 5)
W2[D2 == 0] <- 5
W <- cbind(W1 = rnbinom(N, size = 25,
                        mu = exp(2.5 * U + 0.2 * X)),
           W2)
p2sls_result <- pci2s::p2sls.negbin(Y = Y, A = A, X = X,
                                    W = W, Z = Z,
                                    nco_type = c("negbin", "ah"),
                                    nco_args = list(list(offset = rep(0, N)),
                                                    list(offset = rep(0, N),
                                                         event = D2)))
p2sls_result$summary_first_stage
p2sls_result$summary_second_stage
p2sls_result$ESTIMATE
p2sls_result$SE


p2sls_result$summary_second_stage['A', ] # <- what we are interested in (ATE + SE + pval)


# 3 options:
# - continuous NCE + continuous NCO
# - continuous NCE + count NCO
# - count NCE + count NCO
# with the continuous ones, use PCA 
# with count, use individual genes

# try with one AY test



AY_idx = 1
A_name = AY[AY_idx, 'A']
Y_name = AY[AY_idx, 'Y']

# - continuous NCE + continuous NCO
# - continuous NCE + count NCO


# === count NCE + count NCO ===
# idx of cells to use for this A-Y test
A_idx = which(as.logical(grna_odm[[A_name, ]]))



# Get A
# ------------------------------------------- #
A_idx = which(as.logical(grna_odm[[A_name, ]])) # idx of all 'treated' cells   
control_idx = setdiff(NT_idx, A_idx)                  # idx of control cells (NT without A)
AY_data_idx  = c(A_idx, control_idx)                  # idx of data for this AY test
# A = grna_odm[[A_name, AY_data_idx]] # subset cells of 'treated' (w A grna) and 'control' (NT grna)
A = c(rep(0, length(A_idx)), rep(1, length(control_idx))) # these should be the same


# Get Z and W 
# -------------------------------------------  #
Z_names = AYZW[[A_name]][[Y_name]][[1]]$Z_names[1:max_num_NC_pairs]
W_names = AYZW[[A_name]][[Y_name]][[1]]$W_names[1:max_num_NC_pairs]
dfZ = gene_odm[[Z_names, AY_data_idx]] |> t() |> as.matrix()
dfW = gene_odm[[W_names, AY_data_idx]] |> t() |> as.matrix()
colnames(dfZ) = paste0('Z', 1:ncol(dfZ))
colnames(dfW) = paste0('W', 1:ncol(dfW))


# prob good for future, but keep simple for now
# NCs = readRDS(sprintf('%s/AY/%s/AYZW.rds', save_dir, AYZW_setting_name))
# if(NC_name == 'proximalCountCountCount') {
#   NCs_new = list(Z_names = NCs[[AY[AY_idx, 'A']]][[AY[AY_idx, 'Y']]][[1]]$Z_names[1:max(NCs_settings[[NC_name]]$num_NC_pairs)],
#                  W_names = NCs[[AY[AY_idx, 'A']]][[AY[AY_idx, 'Y']]][[1]]$W_names[1:max(NCs_settings[[NC_name]]$num_NC_pairs)])
#   NCs = NCs_new
#   rm(NCs_new)
# }
# 
# 
# if(is.matrix(NCs)) {
#   dfZ = NCs[AY_data_idx, seq(from = 2, to = ncol(NCs), by = 2)] # evens which((1:ncol(NCs)) %% 2 == 0)
#   dfW = NCs[AY_data_idx, seq(from = 1, to = ncol(NCs), by = 2)] # odds
# } else if(is.list(NCs) && ('Z_names' %in% names(NCs)) && ('W_names' %in% names(NCs)) ) {
#   dfZ = gene_norm[NCs[['Z_names']], AY_data_idx] |> t() |> data.frame()
#   dfW = gene_norm[NCs[['W_names']], AY_data_idx] |> t() |> data.frame()
# } else {
#   print(sprint('(A:%s,Y:%s) Bad NCs input in get_AYZW_df_pci2sbyNCs', A_name, Y_name))
#   return()
# }


# Get Y
# -------------------------------------------  #
Y = gene_odm[[Y_name, AY_data_idx]] |> as.vector()

for(num_NCs in NCs_settings[[NC_name]]$num_NC_pairs) {
  print(sprintf('pci2s negbin #NCs=%02.f', num_NCs))
  pci2s_res = pci2s::p2sls.negbin(
    Y = Y |> as.vector(), 
    A = A, 
    W = dfW[, 1:num_NCs], 
    Z = dfZ[, 1:num_NCs], 
    variance = TRUE,
    nco_type = rep("negbin", num_NCs))
  
}


res = data.frame(       NC_type=character(0),
                        method = character(0),
                        method_type = character(0),
                        numNC =   numeric(0),
                        basis = character(0),
                        ATE =   numeric(0),
                        se  =   numeric(0),
                        pval =   numeric(0),
                        time_sec = numeric(0))

for(num_NCs in NCs_settings[[NC_name]]$num_NC_pairs) {
  # print(num_NCs)
  # chosen_cols =  c('A', 'Y', paste0('Z', 1:num_NCs), paste0('W', 1:num_NCs))
  # if(!all(chosen_cols %in% colnames(df_all))) { # if not all NCs are available (e.g. not enough NCs for the specified num_NCs)
  #   next
  # }
  # df = df_all[, chosen_cols] # subset only these cols
  
  
  
  
  # === === === === === === === === ===  === === ===
  # === Proximal Estimators (only pci2s here)    ===
  # === === === === === === === === ===  === === ===
  
  # === Proximal/Outcome Confounding Bridge, 2 Stage LS w/ pci2s package
  if(which_estimators$OCB_2SLS_pci2s_countcountcount) {
    print(sprintf('pci2s negbin #NCs=%02.f', num_NCs))
    t0 = Sys.time()
    pci2s_res = pci2s::p2sls.negbin(
      Y = Y, 
      A = A, 
      W = dfW[, 1:num_NCs], 
      Z = dfZ[, 1:num_NCs], 
      variance = TRUE,
      nco_type = rep("negbin", num_NCs))
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





# Assemble together
# ------------------------------------------- #
dfAYZW = cbind(data.frame(A = A, 
                          Y = as.vector(Y)),
               dfZ, dfW)


cbind(dfZ, dfW)

if(is.null(U_confounders)) { 
  return(dfAYZW)
} else {
  # Add U (if given)
  # ------------------------------------------- #
  dfU = U_confounders[AY_data_idx, ] 
  colnames(dfU) = paste0('U', 1:ncol(dfU)) # rename cols so no mix up 
  return(cbind(dfAYZW, dfU))
}



}


#' 
#' construct a dataframe of all the data 
#' (outside general function)
#' 
#' @param A_name (character) name of grna/perturbation treatment
#' @param Y_name (character) name of gene outcome
#' @param NT_idx (vector) of indices for cells receiving non-targeting perturbations
#' @param grna_odm (matrix-like) of 0/1 of size #grnas x #cells indicating grna assignment
#'    make sure the rownames of grna are the names of the grna (e.g. A_name is one of the rownames)
#' @param gene_odm (matrix-like) of gene expr outcomes (should be continuous for this pci2s)
#' @param NC_types (character) of what the NCs should be. 
#'       'CountCount' --> count NCE and count NCO using the other genes' raw expression
#'       'ContinuousContinuous' --> continuous NCE and continuous NCO using the 
#'         other input `NCs` (e.g. PC or WGCNA or normalized singlegenes)
#' @param NCs (matrix or list)
#'    matrix if using PCA, Sparse PCA, or WGCNA (or etc...)
#'    list if using individualgene: input here as a list with the structure:
#'          NCs$Z_names   (vector) of characters of gene names
#'          NCs$W_names   (vector) of characters of gene names
#'    Detected through 'if(!is.matrix(NCs))'
#' @param gene_norm
#' @param U_confounders
#' @param library_size
#' @param max_num_NC_pairs (integer) maximum number of NC pairs 
#'     (e.g. if NCs has 100 Z_names and W_names available, but the specified models
#'      only use at most 10 NC pairs, then only extract 10)
#' @example a = format_AYZW_inner(AY_idx = 2); head(a)
get_AYZW_df_pci2s_negbin <- function(A_name, Y_name, NT_idx,
                                   grna_odm, 
                                   gene_odm, 
                                   NC_types,
                                   NCs=NULL, 
                                   gene_norm=NULL,
                                   U_confounders=NULL,
                                   library_size=NULL,
                                   max_num_NC_pairs=NULL) {
  if((!is.null(NC_types)) && !(NC_types %in% c('CountCount', 'ContinuousContinuous'))) {
    print("Bad Input 'NC_types' for function 'get_AYZW_pci2s_negbin': must be 'countcount' or 'continuouscontinuous' (others not implemented yet)")
    return()
  }
  
  # if NC_types=='countcount', then we use other individual genes' counts, and so the NCs must be a list of the Z_names and W_names
  # if NC_types=='continuouscontinuous', then
  #     if PCA/SPCA/WGCNA (is.matrix(NCs)), use the given NCs which must be a matrix
  #     if singlegenes   (!is.matrix(NCs)), the input gene_norm must be given and NCs must be a list of the Z_names and W_names
  
  
  
  
  # Get A ------------------------------------- 
  # ------------------------------------------- #
  A_idx = which(as.logical(grna_odm[[A_name, ]]))       # idx of all 'treated' cells   
  control_idx = setdiff(NT_idx, A_idx)                  # idx of control cells (NT without A)
  AY_data_idx  = c(A_idx, control_idx)                  # idx of data for this AY test
  
  A = as.integer(grna_odm[[A_name, AY_data_idx]]) # subset cells of 'treated' (w A grna) and 'control' (NT grna)
  # A = c(rep(0, length(A_idx)), rep(1, length(control_idx))) # these should be the same
  
  
  # Get Y -------------------------------------
  # ------------------------------------------- #   
  Y = gene_odm[[Y_name, AY_data_idx]] |> as.vector()
  
  
  # Get Z and W  ------------------------------
  # ------------------------------------------- #
  
  # tried to remove the NC_type arg... but can't really, because there are single individual genes for both counts and continuous 
  if(is.null(NC_types)) { # e.g. for simpleCount with regular GLMs
    # do nothing
    Z = NULL
    W = NULL
  } else if(NC_types == 'CountCount') {
    # Count Z and Count W
    
    if(is.null(max_num_NC_pairs)) { # use all the Z and W given
      Z_names = NCs$Z_names
      W_names = NCs$W_names
    } else { # only take needed
      Z_names = NCs$Z_names[1:max_num_NC_pairs]
      W_names = NCs$W_names[1:max_num_NC_pairs]
    }
    Z = gene_odm[[Z_names, AY_data_idx]] |> t() |> as.matrix()
    W = gene_odm[[W_names, AY_data_idx]] |> t() |> as.matrix()
    colnames(Z) = paste0('Z', 1:ncol(Z))
    colnames(W) = paste0('W', 1:ncol(W))
  } else if(NC_types == 'ContinuousContinuous') {
    # Continuous Z and continuous W: could be
    #  - PCA (or WGCNA or SPCA) <-- NCs should be matrix
    #  - singlegenes  <-- NCs should be list of Z and W names
    if(is.matrix(NCs)) {
      Z = NCs[AY_data_idx, seq(from = 2, to = ncol(NCs), by = 2)] # evens which((1:ncol(NCs)) %% 2 == 0)
      W = NCs[AY_data_idx, seq(from = 1, to = ncol(NCs), by = 2)] # odds
    } else if(is.list(NCs) && ('Z_names' %in% names(NCs)) && ('W_names' %in% names(NCs)) ) {
      if(is.null(max_num_NC_pairs)) { # use all the Z and W given
        Z_names = NCs$Z_names
        W_names = NCs$W_names
      } else {                        # only take needed
        Z_names = NCs$Z_names[1:max_num_NC_pairs]
        W_names = NCs$W_names[1:max_num_NC_pairs]
      }
      if(is.null(gene_norm)) {        print('For (ContinuousContinuous & singlegenes) as NCs, must provide gene_norm input'); return()     }
      Z = gene_norm[Z_names, AY_data_idx] |> t() |> as.matrix()
      W = gene_norm[W_names, AY_data_idx] |> t() |> as.matrix()
    } else {
      print(sprint('(A:%s,Y:%s) Bad NCs input in get_AYZW_df_pci2s_negbin', A_name, Y_name))
      return()
    }
    colnames(Z) = paste0('Z', 1:ncol(Z))
    colnames(W) = paste0('W', 1:ncol(W))
  } else {
    print('bad NC_types input (should never get here though, check inputs first)')
    return()
  }
  
  
  df_all = cbind(A=A,
                   Y=Y,
                   Z=Z,
                   W=W) |> data.frame()
  
  if(!is.null(U_confounders)) { 
    # Add U (if given)  -------------------------
    # ------------------------------------------- #
    dfU = U_confounders[AY_data_idx, ] 
    colnames(dfU) = paste0('U', 1:ncol(dfU)) # rename cols so no mix up 
    df_all = cbind(df_all, dfU)
  } 
  
  if(!is.null(library_size)) {
    # Add library size (if given)  --------------
    # ------------------------------------------- #
    df_all$library_size = library_size[AY_data_idx] 
  }
  
  return(df_all)
  
  
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
estimate_ATE_pci2snegbin_make <- function(AY, gene_odm, grna_odm, gene_norm, NT_idx,
                                        NCs_list, NCs_settings, 
                                        library_size, U_confounders=NULL,
                                        save_path=NULL) {
  
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
    # Setup -----------------------------------------------------------------------------------------------------------------------
    res = data.frame(       NC_type=character(0),
                            method = character(0),
                            method_type = character(0),
                            numNC =   numeric(0),
                            basis = character(0),
                            ATE =   numeric(0),
                            se  =   numeric(0),
                            pval =   numeric(0),
                            time_sec = numeric(0))
    
    performed_pois_YA  = FALSE # only need to estimate using lm one time but can be done with any of the NCs_list loops
    performed_pois_YAU = FALSE # set these to TRUE when these methods are done
    performed_nb_YA    = FALSE # (for each of the proximalNegBinXXXX settings, the pois_YA/.../nb_YAU flags may be set s.t. these 
    performed_nb_YAU   = FALSE #  glms are performed within each, but we only need to run it 1x if running a lot at once          )
    
    # Calc for each NC_name (e.g. simpleCount, proximalNegBinCountCount, etc ...) -------------------------------------------------
    for(NC_name in names(NCs_list)) { # e.g. for this pci2s negbin- proximalNegBinCountCount or proximalNegBinContinuousContinuous
      # print(sprintf('[%s] %s', Sys.time(), NC_name))
      NCs = NCs_list[[NC_name]]
      
      # === Construct df ------------------------------------------------------------------------------------------------------
      if(!is.matrix(NCs) && (is.list(NCs))) { # edit NCs into smaller list with just Z_names and W_names, and limit num_NC_pairs actually used
        # print(sprintf('%s Making NC as list of Z and Ws!', NC_name))
        NCs_new = list(Z_names = NCs[[AY[AY_idx, 'A']]][[AY[AY_idx, 'Y']]][[1]]$Z_names[1:max(NCs_settings[[NC_name]]$num_NC_pairs)],
                       W_names = NCs[[AY[AY_idx, 'A']]][[AY[AY_idx, 'Y']]][[1]]$W_names[1:max(NCs_settings[[NC_name]]$num_NC_pairs)])
        NCs = NCs_new
        rm(NCs_new)
      }
      
      df_all = get_AYZW_df_pci2s_negbin(
        A_name       =AY[AY_idx, 'A'], 
        Y_name       =AY[AY_idx, 'Y'], 
        NT_idx       =NT_idx,
        grna_odm     =grna_odm,
        gene_odm     =gene_odm,
        NC_types     =switch(NC_name, 
                             "simpleCount"                       = NULL, 
                             "proximalNegBinCountCount"          ='CountCount',
                             "proximalNegBinPCAPCA"              ="ContinuousContinuous",
                             "proximalNegBinSinglegeneSinglegene"="ContinuousContinuous",
                             "proximalNegBinSPCASPCA"            ="ContinuousContinuous", # hopefully, these will not 
                             "proximalNegBinWGCNAWGNCA"          ="ContinuousContinuous", # be used 
                             NULL),
        NCs          =NCs, 
        gene_norm    =gene_norm, 
        U_confounders=U_confounders,
        library_size =library_size,
        max_num_NC_pairs=NULL) # numNCs limited before called
      
      df_all$A = as.numeric(df_all$A) # convert trtmt A to numeric 0/1 if not already (alt T/F)

      
      
      
      # ATE Estimate -----------------------------------------------------------------------------------------------
      which_estimators = NCs_settings[[NC_name]]$which_estimators
      
      # === === === === === === === === ===  ===
      # === * GLM (Pois and NB) models =========
      # === === === === === === === === ===  ===
      
      # Oracle models using (un)meas conf: Y ~ A + offset(log(library_size)) U1 + U2 +...
      U_colnames = grep('U', colnames(df_all), value = T)
      unmeas_conf_formula = paste0('Y ~ A + offset(log(library_size)) + ', paste(U_colnames, collapse = ' + '))
      # print(sprintf('unmeas conf formula: %s', unmeas_conf_formula))
      
      
      
      # === Poisson Y ~ A      (no confounder adj)
      if((!is.null(which_estimators$pois_YA )) && which_estimators$pois_YA && (!performed_pois_YA)) { # if specified and not already done
        # print('pois_YA')
        t0 = Sys.time()
        pois_YA = glm('Y ~ A + offset(log(library_size))', df_all, family = 'poisson')
        t1 = Sys.time()
        
        res = bind_rows(res, 
                        data.frame(
                          method = 'poisYA',
                          method_type = 'naive',
                          numNC = NA,
                          basis = NA,
                          ATE =    coef(pois_YA )[['A']],
                          se  = summary(pois_YA)$coefficients['A', 'Std. Error'],
                          pval= summary(pois_YA)$coefficients['A', 'Pr(>|z|)'],
                          time_sec = difftime(t1, t0, units = 'secs') |> as.numeric()))
        
        rm(pois_YA, t0, t1)
        performed_pois_YA = TRUE
      }
      
      # === Poisson Y ~ A + Us (   confounder adj)
      if( !is.null(which_estimators$pois_YAU) && which_estimators$pois_YAU && !performed_pois_YAU) { # if specified and not already done
        # print('pois_YAU')
        t0 = Sys.time()
        pois_YAU = glm(unmeas_conf_formula, df_all, family = 'poisson')
        t1 = Sys.time()
        
        res = bind_rows(res, 
                        data.frame(
                          method = 'poisYAU',
                          method_type = 'measuredconfounders',
                          numNC = NA,
                          basis = NA,
                          ATE =    coef(pois_YAU )[['A']],
                          se  = summary(pois_YAU)$coefficients['A', 'Std. Error'],
                          pval= summary(pois_YAU)$coefficients['A', 'Pr(>|z|)'],
                          time_sec = difftime(t1, t0, units = 'secs') |> as.numeric()))
        rm(pois_YAU, t0, t1)
        performed_pois_YAU = TRUE
      }
      
      
      # === Negative Binomial Y ~ A      (no confounder adj)
      if(!is.null(which_estimators$nb_YA) && which_estimators$nb_YA && !performed_nb_YA) { # if specified and not already done
        # print('nb_YA')
        t0 = Sys.time()
        nb_YA = MASS::glm.nb('Y ~ A + offset(log(library_size))', df_all)
        t1 = Sys.time()
        
        # print(summary(nb_YA)$coefficients)
        res = bind_rows(res, 
                        data.frame(
                          method = 'nbYA',
                          method_type = 'naive',
                          numNC = NA,
                          basis = NA,
                          ATE =    coef(nb_YA )[['A']],
                          se  = summary(nb_YA)$coefficients['A', 'Std. Error'],
                          pval= summary(nb_YA)$coefficients['A', 'Pr(>|z|)'],
                          time_sec = difftime(t1, t0, units = 'secs') |> as.numeric()))
        rm(nb_YA, t0, t1)
        performed_nb_YA = TRUE
      }
      
      # === Negative Binomial Y ~ A + Us (   confounder adj)
      if(!is.null(which_estimators$nb_YAU) && which_estimators$nb_YAU & !performed_nb_YAU) { # if specified and not already done
        # print('nb_YAU')
        t0 = Sys.time()
        nb_YAU = MASS::glm.nb(unmeas_conf_formula, df_all)
        t1 = Sys.time()
        
        res = bind_rows(res, 
                        data.frame(
                          method = 'nbYAU',
                          method_type = 'measuredconfounders',
                          numNC = NA,
                          basis = NA,
                          ATE =    coef(nb_YAU )[['A']],
                          se  = summary(nb_YAU)$coefficients['A', 'Std. Error'],
                          pval= summary(nb_YAU)$coefficients['A', 'Pr(>|z|)'],
                          time_sec = difftime(t1, t0, units = 'secs') |> as.numeric()))
        rm(nb_YAU, t0, t1)
        performed_nb_YAU = TRUE
      }
      
      
      # === === === === === === === === ===  === === ===
      # === * Proximal Neg Bin Estimators ======================
      # ===       (pci2s::p2sls.negbin)      === === ===
      # === === === === === === === === ===  === === ===
      if(!is.null(NCs_settings[[NC_name]]$num_NC_pairs)) {
        for(num_NCs in NCs_settings[[NC_name]]$num_NC_pairs) {
          # print(num_NCs)
          chosen_cols =  c('A', 'Y', paste0('Z', 1:num_NCs), paste0('W', 1:num_NCs))
          if(!all(chosen_cols %in% colnames(df_all))) { # if not all NCs are available (e.g. not enough NCs for the specified num_NCs)
            next
          }
          # df = df_all[, chosen_cols] # subset only these cols, no need to do this
          
          
          # === pci2s NegBin Y and Count NCs (need to specify nco types as negbin too) ===
          if(!is.null(which_estimators$OCB_2SLS_pci2s_NegBinCountCount) && 
                which_estimators$OCB_2SLS_pci2s_NegBinCountCount          
             ) {
            # print(sprintf('pci2s negbin (countcount) #NCs=%02.f', num_NCs))
            t0 = Sys.time()
            pci2s_res = pci2s::p2sls.negbin(
              Y = df_all$Y, 
              A = df_all$A, 
              W = df_all[,(grepl('W', colnames(df_all))) & (colnames(df_all) %in% chosen_cols)], 
              Z = df_all[,(grepl('Z', colnames(df_all))) & (colnames(df_all) %in% chosen_cols)], 
              offset = log(df_all$library_size),
              nco_type = rep("negbin", num_NCs),
              nco_args = lapply(X = 1:num_NCs, FUN = function(x){list(init=NA, offset=log(df_all$library_size))}),
              variance = TRUE,
              verbose = FALSE)
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
          
          # === pci2s NegBin Y and continuous NCs (need to specify nco types as continuous) ===
          if((!is.null(which_estimators$OCB_2SLS_pci2s_NegBinPCAPCA) 
                && which_estimators$OCB_2SLS_pci2s_NegBinPCAPCA              ) ||
               (!is.null(which_estimators$OCB_2SLS_pci2s_NegBinSinglegeneSinglegene) 
                && which_estimators$OCB_2SLS_pci2s_NegBinSinglegeneSinglegene)
             ) {
            # print(sprintf('pci2s negbin (continuouscontinuous) #NCs=%02.f', num_NCs))
            t0 = Sys.time()
            pci2s_res = pci2s::p2sls.negbin(
              Y = df_all$Y, 
              A = df_all$A, 
              W = df_all[,(grepl('W', colnames(df_all))) & (colnames(df_all) %in% chosen_cols)], 
              Z = df_all[,(grepl('Z', colnames(df_all))) & (colnames(df_all) %in% chosen_cols)], 
              offset = log(df_all$library_size),
              nco_type = rep("linear", num_NCs),
              nco_args = lapply(X = 1:num_NCs, FUN = function(x){list(init=NA, offset=log(df_all$library_size))}),
              variance = TRUE,
              verbose = FALSE)
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







