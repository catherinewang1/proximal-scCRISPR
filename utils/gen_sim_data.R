# Simulate Data for:
# Exposure A
# Response Y
# Unmeasured Confounding U
# Negative Control Exposure Z
# Negative Control Outcome W

#         beta
#      A ----> Y      
#      ^\     /^
#gammaX  \   / gammaY
#          U
#gammaZ   /  \ gammaW
#       \/   \/
#       Z     W


require(assertthat)
# require(MatchIt)

# helper function: inverse logit to transform R to [0,1]
#        1
#  ------------
#   1 + exp(-x)
invlogit <- function(x) {1/(1+exp(-x))}

#' helper function: softmax to transform vect of R to probs
#' softmax(z)_i = exp(z_i) / sum_j(exp(z_j))
#' @param z (vector) of numeric 
softmax <- function(z) {
  exp(z) / sum(exp(z))
}


#####################################################################################
#####            Simulating Data 
##### Includes:
#####       - sim_data_one: simulate one point from simple model
#####       - sim_data: simulate points from simple model
#####       - sim_cell_one: simulate one point from cell model
#####       - sim_cell_data: simulate points from cell model
#####       - match_cells: match 
#####################################################################################

#' Simulate ONE point of Data according to the causal dag above
#' with the following parameters for the linear effects
#' @param U_mean (vector) vector of unconfounders' mean
#' @param U_sd  (vector) vector of unonfounders' sd
#' @param beta  (numeric) linear effect of A on Y
#' @param gammaA (vector) of numeric of length=dim(U) linear effect (logistic) of U on A 
#' @param gammaY (vector) of numeric of length=dim(U) linear effect            of U on Y
#' @param gammaZ (vector) of numeric of length=dim(U) linear effect            of U on Z
#' @param gammaW (vector) of numeric of length=dim(U) linear effect            of U on W
#' @param muA (numeric) muA constant for X draw invlogit(muA + ...) 
#' @param sigmaY (numeric) > 0 sd of normal noise of Y
#' @param sigmaZ (numeric) > 0 sd of normal noise of Z
#' @param sigmaW (numeric) > 0 sd of normal noise of W
sim_data_one <- function(U_mean, U_sd, beta, gammaA, gammaY, gammaZ, gammaW, muA=0, sigmaY=1, sigmaZ=1, sigmaW=1) {
    # draw U: vector of R
    U = sapply(1:length(U_mean), FUN = function(i) {rnorm(n=1, mean=U_mean[i], sd=U_sd[i])})

    # draw W: R
    W = sum(gammaW*U) + rnorm(n=1, mean=0, sd=sigmaW)

    # draw Z: {0,1}
    # Z = rbinom(n=1, size=1, prob=invlogit(sum(gammaZ*U)))
    # draw Z: R
    Z = sum(gammaZ*U) + rnorm(n=1, mean=0, sd=sigmaZ)

    # draw A: {0,1}
    A = rbinom(n=1, size=1, prob=invlogit(muA + gammaA*U))

    # draw Y: R
    Y = sum(gammaY*U) + X*beta + rnorm(n=1, mean=0, sd=sigmaY)
     
    return(c(Y, A, Z, W, U))
    # format
    # U_df = data.frame(t(matrix(U)))
    # colnames(U_df) = paste0('U', 1:length(U_mean))

    # return(cbind(data.frame(Y=Y, W=W, A=A, Z=Z), U_df))
}


#' Simulate Data according to the causal dag above
#' with the following parameters for the linear effects
#' @param N (integer) number of data points
#' @param U_dim (vector) vector of unconfounders' mean
#' @param U_sd  (vector) vector of unonfounders' sd
#' @param beta  (numeric) linear effect of A on Y
#' @param gammaA (vector) of numeric of length=dim(U) linear effect (logistic) of U on A 
#' @param gammaY (vector) of numeric of length=dim(U) linear effect            of U on Y
#' @param gammaZ (vector) of numeric of length=dim(U) linear effect            of U on Z
#' @param gammaW (vector) of numeric of length=dim(U) linear effect            of U on W
sim_data <- function(N, U_mean, U_sd, beta, gammaA, gammaY, gammaZ, gammaW, muA=0, sigmaY=1, sigmaZ=1, sigmaW=1) {
    # check inputs
    assertthat::assert_that(length(U_mean) == length(U_sd))
    assertthat::assert_that(length(gammaX) == length(U_mean))
    assertthat::assert_that(length(gammaY) == length(U_mean))
    assertthat::assert_that(length(gammaZ) == length(U_mean))
    assertthat::assert_that(length(gammaW) == length(U_mean))
    assertthat::assert_that(length(muA)    == 1)

    # generate data
    dat = replicate(N, sim_data_one(U_mean = U_mean,
                                    U_sd   = U_sd,
                                    beta   = beta,
                                    gammaA = gammaA,
                                    gammaY = gammaY,
                                    gammaZ = gammaZ,
                                    gammaW = gammaW,
                                    muA    = muA,
                                    sigmaY = sigmaY,
                                    sigmaZ = sigmaZ,
                                    sigmaW = sigmaW))
    # format data
    dat = data.frame(t(dat))
    colnames(dat) = c('Y','A', 'Z', 'W', paste0('U', 1:length(U_mean)))
    return(dat)
}


#' Simulate ONE point of Data according to the causal dag above
#' with the following parameters for the linear effects
#' @param U_mean (vector) vector of unconfounders' mean
#' @param U_sd  (vector) vector of unonfounders' sd
#' @param beta  (numeric) linear effect of A on Y
#' @param gammaA (vector) of numeric of length=dim(U) linear effect (logistic) of U on A 
#' @param gammaY (vector) of numeric of length=dim(U) linear effect            of U on Y
#' @param gammaZ (matrix) of numeric of dim= dim(U)  x (#Zs) linear effect     of U on Z
#' @param gammaW (matrix) of numeric of dim= dim(U)  x (#Ws) linear effect     of U on W
#' @param muA    (numeric) mean intercept for A
#' @param sigmaY (numeric) sd of Y's noise
#' @param sigmaZ (vector or numeric) > 0 of sd of Z's noise
#' @param sigmaW (vector or numeric) > 0 of sd of W's noise
sim_data_multidimU_one <- function(U_mean, U_sd, beta, gammaA, gammaY, gammaZ, gammaW, 
                                   muA=0, sigmaY=1, sigmaZ=1, sigmaW=1) {
    # draw U: vector of R
    U = sapply(1:length(U_mean), FUN = function(i) {rnorm(n=1, mean=U_mean[i], sd=U_sd[i])})

    # draw W: R^ncol(gammaW)=#Ws
    W = t(U) %*% gammaW + rnorm(ncol(gammaW), mean = 0, sd = sigmaW)

    # draw Z: {0,1}
    # Z = rbinom(n=1, size=1, prob=invlogit(sum(gammaZ*U)))
    # draw Z: R
    # Z = sum(gammaZ*U) + rnorm(n=1, mean=0, sd=1)
    Z = t(U) %*% gammaZ + rnorm(ncol(gammaZ), mean = 0, sd = sigmaZ)

    # draw A: {0,1}
    A = rbinom(n=1, size=1, prob=invlogit(sum(gammaA*U)))

    # draw Y: R
    Y = sum(gammaY*U) + A*beta + rnorm(n=1, mean=0, sd=sigmaY)
     
    return(c(Y, A, Z, W, U))
    # format
    # U_df = data.frame(t(matrix(U)))
    # colnames(U_df) = paste0('U', 1:length(U_mean))

    # return(cbind(data.frame(Y=Y, W=W, A=A, Z=Z), U_df))
}


#' Simulate Data according to the causal dag above
#' with the following parameters for the linear effects
#' @param N (integer) number of data points
#' @param U_dim (vector) vector of unconfounders' mean
#' @param U_sd  (vector) vector of unonfounders' sd
#' @param beta  (numeric) linear effect of A on Y
#' @param gammaA (vector) of numeric of length=dim(U) linear effect (logistic) of U on A 
#' @param gammaY (vector) of numeric of length=dim(U) linear effect            of U on Y
#' @param gammaZ (vector) of numeric of length=dim(U) linear effect (logistic) of U on Z
#' @param gammaW (matrix) of numeric of dim= dim(U)  x (#Ws) linear effect     of U on W
#' @param muA    (numeric) mean intercept for A
#' @param sigmaY (numeric) sd of Y's noise
#' @param sigmaZ (vector or numeric) > 0 of sd of Z's noise
#' @param sigmaW (vector or numeric) > 0 of sd of W's noise
sim_data_multidimU <- function(N, U_mean, U_sd, beta, gammaA, gammaY, gammaZ, gammaW,
                               muA=0, sigmaY=1, sigmaZ=1, sigmaW=1) {
    # check inputs
    assertthat::assert_that(length(U_mean) == length(U_sd))
    assertthat::assert_that(length(gammaA) == length(U_mean))
    assertthat::assert_that(length(gammaY) == length(U_mean))
    assertthat::assert_that(  nrow(gammaZ) == length(U_mean))
    assertthat::assert_that(  nrow(gammaW) == length(U_mean))
    assertthat::assert_that(length(muA)    == 1)

    # generate data
    dat = replicate(N, sim_data_multidimU_one(U_mean = U_mean,
                                    U_sd   = U_sd,
                                    beta   = beta,
                                    gammaA = gammaA,
                                    gammaY = gammaY,
                                    gammaZ = gammaZ,
                                    gammaW = gammaW,
                                    muA    = muA,
                                    sigmaY = sigmaY,
                                    sigmaZ = sigmaZ,
                                    sigmaW = sigmaW))
    # format data
    dat = data.frame(t(dat))
    colnames(dat) = c('Y','A',
                      paste0('Z', 1:ncol(gammaZ)),
                      paste0('W', 1:ncol(gammaW)),
                      paste0('U', 1:length(U_mean)))
    return(dat)
}

# #' Simulate one cell. 
# #' gRNA1 is A/has beta effect on GENE1 (is Y)
# #' @param U_mean (vector) vector of unconfounders' mean
# #' @param U_sd   (vector) vector of unonfounders' sd
# #' @param beta   (numeric) linear effect of X on Y
# #' @param gRNA_mean (vector) vector of gRNA mean (in logistic fn)
# #' @param gammaGRNA (matrix) #gRNA x dim U for chance of each cell 
# #'                  for receiving each gRNA influenced by U
# #'                  gammaGRNA %*% U 
# #' @param gammaGENE (matrix) #genes x dim U for linear effect of U
# #'                 on gene expression
# #'
# sim_cell_one <- function(U_mean, U_sd, beta, gRNA_mean, gammaGRNA, gammaGENE) {
#   # draw U: vector of R
#   U = sapply(1:length(U_mean), FUN = function(i) {rnorm(n=1, mean=U_mean[i], sd=U_sd[i])})
  
#   # draw gRNA: vector of length(gRNA_mean)=nrow(gammaGRNA). 
#   gRNA_probs = apply(gammaGRNA %*% U + gRNA_mean, MARGIN = 2, softmax)
#   gRNA_i     = sample.int(nrow(gammaGRNA), size = 1, prob = gRNA_probs)
#   gRNA       = rep(0, nrow(gammaGRNA)); gRNA[gRNA_i] = 1
  
#   # draw G: gene expressions vector of R
#   gene_means = gammaGENE %*% U
#   GENE = sapply(gene_means, FUN = function(m) {rnorm(n = 1, mean=m, sd = 1)}) |>
#     matrix(nrow=nrow(gene_means))
  
#   GENE[1] = GENE[1] + gRNA[1]*beta # add effect of gRNA1
  
#   return(c(gRNA, GENE, U))
  
#   # sapply(matrix(1:20, nrow=5), FUN = function(m) {rnorm(n = 1, mean=m, sd = .01)}) |>
#   #   matrix(nrow=5)
  
# }


# #' Simulate cells
# #' gRNA1 (is X) has beta effect on GENE1 (is Y)
# #' @param N (integer) number of data points
# #' @param U_mean (vector) vector of unconfounders' mean
# #' @param U_sd   (vector) vector of unonfounders' sd
# #' @param beta   (numeric) linear effect of X on Y
# #' @param gRNA_mean (vector) vector of gRNA mean (in logistic fn)
# #' @param gammaGRNA (matrix) #gRNA x dim U for chance of each cell 
# #'                  for receiving each gRNA influenced by U
# #'                  gammaGRNA %*% U 
# #' @param gammaGENE (matrix) #genes x dim U for linear effect of U
# #'                 on gene expression
# #'
# sim_cell <- function(N, U_mean, U_sd, beta, gRNA_mean, gammaGRNA, gammaGENE) {
#   # check inputs
#   assertthat::assert_that((N %% 1 == 0) & (N > 0))
#   assertthat::assert_that(length(U_mean)  == length(U_sd))
#   assertthat::assert_that(ncol(gammaGRNA) == length(U_mean))
#   assertthat::assert_that(ncol(gammaGENE) == length(U_mean))
  
#   # generate data
#   dat = replicate(N, sim_cell_one(U_mean = U_mean,
#                                   U_sd   = U_sd,
#                                   beta   = beta,
#                                   gRNA_mean = gRNA_mean,
#                                   gammaGRNA = gammaGRNA,
#                                   gammaGENE = gammaGENE))
#   # format data
#   dat = data.frame(t(dat))
#   colnames(dat) = c(paste0('gRNA', 1:nrow(gammaGRNA)),
#                     paste0('GENE', 1:nrow(gammaGENE)), 
#                     paste0('U',    1:length(U_mean)))
#   return(dat)
# }


# #' Matches cells 
# #' @param df (dataframe) must have a gRNA1 and GENE1 column
# #' @param NT_names (vector of chars) colnames of NonTargeting gRNA
# #' @param Z_names (vector of chars) colnames of gRNA for Z (Z = 1 if any of 
# #'                these columns are 1)
# #' @param match_cov_names (vector of chars) colnames of covariates to match on
# #' @return dataframe with columns: i, Y, X, Z, W, mi
# #' @example 
# #' NT_names = paste0('gRNA', 2:3)
# #' Z_names = c('gRNA3') # Z = 1 if any of these cols is 1
# #' GENE_names = paste0('GENE', 2:20)
# #' cell_matched = match_cells(df=cell_df, 
# #' NT_names = paste0('gRNA', 2:3),
# #' Z_names  = c('gRNA3'),
# #' GENE_names = paste0('GENE', 2:20))
# #' calc_ACE(cell_matched )
# match_cells <- function(df, NT_names, Z_names, match_cov_names) {
  
#   # idx of Targeting and NonTargeting 
#   Tidx  = which(df$gRNA1 == 1)
#   NTidx = which(apply(df[,NT_names], MARGIN=1, FUN=sum) == 1)
#   NTidx_0 = NTidx[1:round((length(NTidx) - length(Tidx))/10)]
#   Originalidx = c(Tidx, NTidx_0)
  
#   # make new col indicating if original cell or matched cell 
#   df$OriginalCell = 0
#   df$OriginalCell[Originalidx] = 1
  
#   # make new col indicating the Z variable
#   if(length(Z_names) == 1) { df$Z = df[,Z_names] 
#   } else {                  df$Z = apply(df[,Z_names], MARGIN=1, FUN=sum)}
  
  
#   # Match Original cells to Matched cells
#   cell_matches = MatchIt::matchit(
#     as.formula(paste0('OriginalCell ~ ',
#                       paste0(match_cov_names, collapse = ' + '))),
#     data = df)
#   # TESTING FUNCTION: Matching on U1
#   # cell_matches = MatchIt::matchit(
#   #   as.formula(paste0('OriginalCell ~ U1')), 
#   #   data = df)
  
#   # construct dataframe
#   Original_df = data.frame(i = rownames(cell_matches$match.matrix))
#   Original_df = cbind(Original_df, 
#                       df[Original_df$i, c('GENE1', 'gRNA1', 'U1')] |> 
#                         rename(X = gRNA1, Y = GENE1, U_Original=U1))
  
#   Matched_df = data.frame(mi = cell_matches$match.matrix)
#   Matched_df = cbind(df[Matched_df$mi, c('Z', 'GENE1', 'U1')] |>
#                        rename(W = GENE1, U_Matched=U1),
#                      Matched_df)
  
  
  
#   return(cbind(Original_df, Matched_df))
# }



# #' Perform pca dimension reduction on input df (all the cols of df)
# #' @param df (dataframe)
# #' @param rank (integer) rank of the pca dimension reduction
# dimreduce_cells <- function(df, rank) {
#   pca_fit = df |> prcomp(rank. = rank, center = TRUE, scale = TRUE)
#   pca_df = data.frame(pca_fit$x)
#   colnames(pca_df) = paste0('PCA', 1:rank)
#   return(pca_df)
  
# }



# #####################################################################################
# #####            Estimating ACE/ATE with binary Exposure and binary NCExposure
# ##### Includes:
# #####       - calc_RD: calculate Risk Difference
# #####       - calc_ERD: calculated Expected Risk Difference
# #####       - calc_ACE: calculate Average Causal Effect (aka ATE)
# #####################################################################################

# #' Calculate Risk Difference (RD) as defined in Miao, Shi, Tchetgen Tchetgen
# #'          RD_{AB|C} := E(B|A=1, C) - E(B|A=0, C)
# #' Risk difference of A when comparing B=1 vs B=0 and conditioning on C
# #' @param A (vector) of A values (binary 0/1)
# #' @param B (vector) of B values 
# #' @param C (vector) of C values (typically binary 0/1)
# #' @param c (numeric) value of C to condition on 
# calc_RD <- function(A, B, C, c) {
#   #E(B|A=1, C=c)
#   E1 = mean(B[which(A == 1 & C == c)])
  
#   #E(B|A=0, C=c)
#   E2 = mean(B[which(A == 0 & C == c)])
  
#   return(E1 - E2)
# }


# #' Calculate Expected Risk Difference E_C(RD_{AB|C})
# #' @param A (vector) of A values (binary 0/1)
# #' @param B (vector) of B values 
# #' @param C (vector) of C values (typically binary 0/1)
# calc_ERD <- function(A, B, C) {
#   c_s = unique(C)
#   RDs = c(); c_count = c()
#   for(c in c_s) {
#     RD_c = calc_RD(A, B, C, c=c)  
#     RDs = c(RDs, RD_c)
#     c_count = c(c_count, sum(C == c))
#   }
#   return(sum(RDs*c_count)/length(C))
# }

# #' Calculate the Average Causal Effect of X on Y with binary X, Z
# #' (Example 8)
# #' ERD(X, Y, Z) - (ERD(Z, Y, X) / ERD(Z, W, X)) * ERD(X, W, Z)
# #' @param df (dataframe) dataframe with variables Y, W, X, Z
# calc_ACE_binaryZ <- function(df) {
#   X = df$X
#   Y = df$Y
#   Z = df$Z
#   W = df$W
  
#   # gamma2 =  calc_RD(Z, Y, X, 0) / calc_RD(Z, W, X, 0)
#   # gamma3 = (calc_RD(Z, Y, X, 1) / calc_RD(Z, W, X, 1)) - gamma2
#   # calc_ERD(X, Y, Z) - (calc_ERD(Z, Y, X) / calc_ERD(Z, W, X)) * calc_ERD(X, W, Z)
#   E1 = calc_ERD(X, Y, Z)
#   gamma2hat = (calc_ERD(Z, Y, X) / calc_ERD(Z, W, X))
#   E2 = calc_ERD(X, W, Z)
#   # print(sprintf('%.2f -- %.2f -- %.2f', E1, gamma2hat, E2))
#   return(E1 - gamma2hat*E2) 
# }



# #' Performs pca on Z, W separately and puts back
# #' into dataframe format w/ same names used for OCB/prox functions
# #'
# #'
# #' @param df
# #' @param r_z (integer) rank/#pcs for pca on Z
# #' @param r_w (integer) rank/#pcs for pca on W
# #' @param colnames
# #' @return dataframe with original A,Us,
# #'               and with new PCA Zs, Ws
# getdfPCA <- function(df, r_z, r_w) {
#   # r_z = 6; r_w = 6
#   Z_colnames = grep('Z', colnames(df), value = T)
#   W_colnames = grep('W', colnames(df), value = T)

#   Z = df[, Z_colnames]
#   W = df[, W_colnames]

#   Z_pca = prcomp(Z, center = TRUE, scale. = TRUE)
#   W_pca = prcomp(W, center = TRUE, scale. = TRUE)

#   Z_ = Z_pca$x[,1:r_z, drop=FALSE]
#   W_ = W_pca$x[,1:r_w, drop=FALSE]
#   colnames(Z_) = paste0('Z',1:ncol(Z_))
#   colnames(W_) = paste0('W',1:ncol(W_))

#   df_ = cbind(df[, !(colnames(df) %in% c(Z_colnames, W_colnames))],
#               Z_,
#               W_)

#   return(df_)
# }




