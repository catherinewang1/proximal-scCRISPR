



# set all num_NC_pairs if wanted (or set to NULL)
# default_num_NC_pairs = c(1, 5, 10) # set small values for testing functions
default_num_NC_pairs = NULL

proximal_continuous_settings = list(
    'PCA'     = list(
                 # === NC Parameters
                 NC_type      = 'PCA',
                 NC_name      = 'PCA',
                 num_NC_pairs = if(!is.null(default_num_NC_pairs)) default_num_NC_pairs else c(1, 3, 5, 8, 10, 15, 20),
                 # extra params for sparsePCA 
                 my_sumabsv   = NA, 
                 my_K         = NA,
                 N_subsample  = NA, 
                 # === procedure parameters
                 save_intermediateATEs = 'yes',        # 'yes'/'no' whether to save intermedate ATEs as they are estimated (should be removed, specify in overall script. keep for later if it is helpful to do one by one, but probably not.)
                 # === Parameter Settings for which estimators to perform
                 which_estimators = list(
                          lm_YA        = TRUE,
                          lm_YAU       = TRUE,
                          # pois_YAU     = TRUE,
                          # nb_YAU       = TRUE,
                          # OCB_2SLS     = FALSE,
                          OCB_2SLS_pci2s=TRUE #,
                          # OCB_2SLSReg  = FALSE,
                          # OCB_GMM      = FALSE,
                          # OCB_GMMRw    = FALSE,
                          # OCB_GMMRwReg = FALSE,
                          # OCB_LinOSPI  = FALSE,
                          # OCB_LinOS    = FALSE,
                          # OCB_LinOStrim= FALSE
                        )

                 ),
    'SPCA8.0'  = list(
                 # === NC Parameters
                 NC_type      = 'SPCA',
                 NC_name      = 'SPCA8.0',
                 num_NC_pairs = if(!is.null(default_num_NC_pairs)) default_num_NC_pairs else c(1, 3, 5, 8, 10, 15, 20),
                 # extra params for sparsePCA 
                 my_sumabsv   = 8, 
                 my_K         = 60,
                 N_subsample  = 5000, # subsample size, or 'all' if using all cells
                 # === procedure parameters
                 save_intermediateATEs = 'yes',        # 'yes'/'no' whether to save intermedate ATEs as they are estimated
                 # === Parameter Settings for which estimators to perform
                 which_estimators = list(
                          lm_YA        = TRUE,
                          lm_YAU       = TRUE,
                          # pois_YAU     = TRUE,
                          # nb_YAU       = TRUE,
                          # OCB_2SLS     = FALSE,
                          OCB_2SLS_pci2s=TRUE #,
                          # OCB_2SLSReg  = FALSE,
                          # OCB_GMM      = FALSE,
                          # OCB_GMMRw    = FALSE,
                          # OCB_GMMRwReg = FALSE,
                          # OCB_LinOSPI  = FALSE,
                          # OCB_LinOS    = FALSE,
                          # OCB_LinOStrim= FALSE
                        )

                 ),    
    'SPCA34.5' = list(
                 # === NC Parameters
                 NC_type      = 'SPCA',
                 NC_name      = 'SPCA34.5',
                 num_NC_pairs = if(!is.null(default_num_NC_pairs)) default_num_NC_pairs else c(1, 3, 5, 8, 10, 15, 20),
                 # extra params for sparsePCA 
                 my_sumabsv   = 34.5, 
                 my_K         = 60,
                 N_subsample  = 5000, # subsample size, or 'all' if using all cells
                 # === procedure parameters
                 save_intermediateATEs = 'yes',        # 'yes'/'no' whether to save intermedate ATEs as they are estimated
                 # === Parameter Settings for which estimators to perform
                 which_estimators = list(
                          lm_YA        = TRUE,
                          lm_YAU       = TRUE,
                          # pois_YAU     = TRUE,
                          # nb_YAU       = TRUE,
                          # OCB_2SLS     = FALSE,
                          OCB_2SLS_pci2s=TRUE #,
                          # OCB_2SLSReg  = FALSE,
                          # OCB_GMM      = FALSE,
                          # OCB_GMMRw    = FALSE,
                          # OCB_GMMRwReg = FALSE,
                          # OCB_LinOSPI  = FALSE,
                          # OCB_LinOS    = FALSE,
                          # OCB_LinOStrim= FALSE
                        )

                 ), 
    'WGCNA'   = list(
                 # === NC Parameters
                 NC_type      = 'WGCNA',
                 NC_name      = 'WGCNA',
                 num_NC_pairs = if(!is.null(default_num_NC_pairs)) default_num_NC_pairs else c(1, 3, 5, 7, 9),
                 # extra params for sparsePCA 
                 my_sumabsv   = NA, 
                 my_K         = NA,
                 N_subsample  = NA, 
                 # === procedure parameters
                 save_intermediateATEs = 'yes',        # 'yes'/'no' whether to save intermedate ATEs as they are estimated
                 # === Parameter Settings for which estimators to perform
                 which_estimators = list(
                          lm_YA        = TRUE,
                          lm_YAU       = TRUE,
                          # pois_YAU     = TRUE,
                          # nb_YAU       = TRUE,
                          # OCB_2SLS     = FALSE,
                          OCB_2SLS_pci2s=TRUE #,
                          # OCB_2SLSReg  = FALSE,
                          # OCB_GMM      = FALSE,
                          # OCB_GMMRw    = FALSE,
                          # OCB_GMMRwReg = FALSE,
                          # OCB_LinOSPI  = FALSE,
                          # OCB_LinOS    = FALSE,
                          # OCB_LinOStrim= FALSE
                        )

                 ),
    'singlegene'   = list(
                 # === NC Parameters
                 NC_type      = 'singlegene',
                 NC_name      = 'singlegene',
                 num_NC_pairs = if(!is.null(default_num_NC_pairs)) default_num_NC_pairs else c(3, 5, 10, 25, 50),
                 # extra params for sparsePCA 
                 my_sumabsv   = NA, 
                 my_K         = NA,
                 N_subsample  = NA, 
                 # === procedure parameters
                 save_intermediateATEs = 'yes',        # 'yes'/'no' whether to save intermedate ATEs as they are estimated
                 # === Parameter Settings for which estimators to perform
                 which_estimators = list(
                          lm_YA        = TRUE,
                          lm_YAU       = TRUE,
                          # pois_YAU     = TRUE,
                          # nb_YAU       = TRUE,
                          # OCB_2SLS     = FALSE,
                          OCB_2SLS_pci2s=TRUE #,
                          # OCB_2SLSReg  = FALSE,
                          # OCB_GMM      = FALSE,
                          # OCB_GMMRw    = FALSE,
                          # OCB_GMMRwReg = FALSE,
                          # OCB_LinOSPI  = FALSE,
                          # OCB_LinOS    = FALSE,
                          # OCB_LinOStrim= FALSE
                        )

                 )

    )




rm(default_num_NC_pairs)