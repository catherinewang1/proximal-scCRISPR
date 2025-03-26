#####################################################################################
#####            Estimating ACE/ATE with binary Exposure and continuous Outcome
##### Using 2 Stage Least Squares (2 versions: LS and Regularized [Lasso,Ridge,elastic net,...])
##### Includes:
#####       - OCB2SLS
#####       - OCB2SLSreg
#####################################################################################



#' Calculate the Average Causal Effect of A on Y 
#' and multidimensional Unmeasured Confounding U
#' NOTE: Intermediate ATE calcs use ALL Z variables 
#' (e.g. 2nd stage is W3 ~ Z1 + Z2 + Z3 + Z4 + Z5, and
#'       NOT          W3 ~ Z1 + Z2 + Z3'               )
#' @param df (dataframe) dataframe with variables Y, W, A, Z
#' @return list of:
#'   EY = function that takes in a and estimates E(Y^a)
#'   EY1 = est expected value of potential outcome Y^1
#'   EY0 = est expected value of potential outcome Y^0
#'   ATE = estimate Average Causal Effect
OCB2SLS <- function(df, returnIntermediateATEs=FALSE) {
  
  Znames = grep('Z', colnames(df), value = TRUE)
  Wnames = grep('W', colnames(df), value = TRUE)
  # Estimate E(Wj | Z, A)
  W_ZA = matrix(NA, nrow = nrow(df), ncol = length(Wnames))
  for(j in 1:length(Wnames)) {
    W_ZA[, j] = lm(sprintf('W%d ~ %s + A', 
                           j, paste0(Znames, collapse = ' + ')),
                   df)$fitted.values
  }
  # format estimated E(Wj|Z,A) and observed Y, A into one dataframe 
  df2 = data.frame(W_ZA)
  colnames(df2) = Wnames
  df2 = cbind(df2, df|>select(Y, A))
  
  # coefs for: intercept, W1, ..., WJ, A
  est_coef = as.vector(lm(sprintf('Y ~ %s + A', 
                                  paste0(Wnames, collapse = ' + ')), 
                          df2)$coefficients)
  
  EY <- function(a) {
    # coefs^T [1, mean(W1), ..., mean(WJ), x]
    sum(est_coef * c(1, colMeans(df[,Wnames, drop=FALSE]), a))
  }
  EY0 = EY(0)
  EY1 = EY(1)

 
  res = list(EY=EY,
              EY0=EY0,
              EY1=EY1,
              ATE=EY1 - EY0)
  

  # ATE using 1, 2, ..., J NCOs
  if(returnIntermediateATEs) {
    # A coef for A from using up to j = 1, 2, ...
    ATEs = rep(NA, length(Wnames))
    for(j in 1:length(Wnames)) {
      ATEs[j] = lm(sprintf('Y ~ %s + A', 
                                    paste0(Wnames[1:j], collapse = ' + ')), 
                            df2)$coefficients[['A']]
    }

    res$ATEs = ATEs
  }
  
  return(res)  
}




#' Calculate the Average Causal Effect of A on Y with continuous Z
#' and multidimensional Unmeasured Confounding U
#' DOES ALL 1:numW ACE calculations (for performance)
#' NOTE: Intermediate ATE calcs use ALL Z variables 
#' (e.g. 2nd stage is W3 ~ Z1 + Z2 + Z3 + Z4 + Z5, and
#'  NOT               W3 ~ Z1 + Z2 + Z3'               )
#' @param df (dataframe) dataframe with variables Y, W, A, Z
#' @param alpha (numeric) [0,1] in glmnet (0=ridge, (0,1)=elastic net, 1=lasso)
#' @param returnMiddleDFs (boolean) whether or not to return dfs made in the middle
#' @return list of list of:
#'   EY = function that takes in x and estimates E(Y^x)
#'   EY1 = est expected value of potential outcome Y^1
#'   EY0 = est expected value of potential outcome Y^0
#'   ACE = estimate Average Causal Effect
OCB2SLSReg <- function(df, alpha = 0, returnMiddleDFs=FALSE) {
  
  Znames = grep('Z', colnames(df), value = TRUE)
  Wnames = grep('W', colnames(df), value = TRUE)

  
  W_ZA = matrix(NA, nrow = nrow(df), ncol = length(Wnames))  # Estimates E(Wj | Z, X)
  AZ   = df[,c('A', paste0('Z', 1:length(Znames)))] |> as.matrix() # covariate matrix
  Wjs_bad = c() # when est of Wj is bad (e.g. constant #), we exclude in 2nd Stage LS
  res = list()

  for(j in 1:length(Wnames)) {
    # ----------- Estimate E(Wj | Z, X) -----------
    # response vector
    Wj = df[,paste0('W', j)]
    
    # ???? ridge regression (alpha = 0) does poorly?
    # cross validation for lambda penalty selection
    cv_model  = glmnet::cv.glmnet(x=AZ, y=Wj, alpha = alpha) # plot(cv_model)
    # fit final model
    # model     = glmnet::glmnet(x = AZ, y = Wj, alpha = alpha, lambda = 0)
    model     = glmnet::glmnet(x = AZ, y = Wj, alpha = alpha, lambda = cv_model$lambda.1se) # TODO: change 1se
    W_ZA[, j] = predict(model, newx = AZ)
    
    if(sd(W_ZA[, j]) == 0) {Wjs_bad = c(Wjs_bad, j)} # bad hat Wj if constant

    # ----------- calculate ATE -----------
    # format estimated E(Wj|Z,A) and observed Y, A into one dataframe (exclude bad hat Wjs)
    df2 = data.frame(W_ZA[,setdiff(1:j, Wjs_bad)])
    colnames(df2) = Wnames[setdiff(1:j, Wjs_bad)]
    df2 = cbind(df2, df|>select(Y, A))
    
    # coefs for: intercept, W1, ..., WJ, A
    est_coef = as.vector(lm(sprintf('Y ~ %s + A', 
                                    paste0(Wnames[setdiff(1:j, Wjs_bad)], collapse = ' + ')), 
                            df2)$coefficients)
    est_coef[is.na(est_coef)] = 0 # sometimes design mat not full rank, some coef NA, replace with 0
    # print(est_coef)
    Wmeans = colMeans(df[,Wnames[setdiff(1:j, Wjs_bad)], drop=FALSE])

    EY <- function(a) {
      # coefs^T [1, mean(W1), ..., mean(WJ), x]
      sum(est_coef * c(1, Wmeans, a))
    }
    EY0 = EY(0)
    EY1 = EY(1)
    
    res[[j]] = list( EY=EY,
                     EY0=EY0,
                     EY1=EY1,
                     ATE=EY1 - EY0)
    if(returnMiddleDFs) {
      # print(dim(df2))
      res[[j]][['df2']] = df2
    }
  }
  
  return(res)
}




# #' Estimate the Treatment Confounding Bridge function
# #' and multidimensional Unmeasured Confounding U
# #' and uses elastic net regularization to estimate E(W_j | X,Z)
# #' and uses logistic regression to estimate E(W_j | X,Z)
# #' DOES ALL 1:numW ACE calculations (for performance)
# #' @param df (dataframe) dataframe with variables Y, W, X, Z
# #' @param alpha (numeric) [0,1] in glmnet (0=ridge, (0,1)=elastic net, 1=lasso)
# #' @return list of list of:
# #'   EY = function that takes in x and estimates E(Y^x)
# #'   EY1 = est expected value of potential outcome Y^1
# #'   EY0 = est expected value of potential outcome Y^0
# #'   ACE = estimate Average Causal Effect
# est_TrtCB <- function(df) {
#   Znames = grep('Z', colnames(df), value = TRUE)
#   Wnames = grep('W', colnames(df), value = TRUE)

  
#   W_ZX = matrix(NA, nrow = nrow(df), ncol = length(Wnames))  # Estimates E(Wj | Z, X)
#   XZ   = df[,c('X', paste0('Z', 1:numZ))] |> as.matrix() # covariate matrix
#   res = list()

#   for(j in 1:length(Wnames)) {
#     # ----------- Estimate E(Wj | Z, X) -----------
#     # response vector
#     Wj = df[,paste0('W', j)]
    
#     # ???? ridge regression (alpha = 0) does poorly?
#     # cross validation for lambda penalty selection
#     cv_model  = glmnet::cv.glmnet(x=XZ, y=Wj, alpha = alpha) # plot(cv_model)
#     # fit final model
#     # model     = glmnet::glmnet(x = XZ, y = Wj, alpha = alpha, lambda = 0)
#     model     = glmnet::glmnet(x = XZ, y = Wj, alpha = alpha, lambda = cv_model$lambda.min)
#     W_ZX[, j] = predict(model, newx = XZ)
  
#     # ----------- calculate ACE -----------
#     # format estimated E(Wj|Z,X) and observed Y, X into one dataframe 
#     df2 = data.frame(W_ZX[,1:j])
#     colnames(df2) = Wnames[1:j]
#     df2 = cbind(df2, df|>select(Y, X))
    
#     # coefs for: intercept, W1, ..., WJ, X
#     est_coef = as.vector(lm(sprintf('Y ~ %s + X', 
#                                     paste0(Wnames[1:j], collapse = ' + ')), 
#                             df2)$coefficients)
    
#     Wmeans = colMeans(df[,Wnames[1:j], drop=FALSE])

#     EY <- function(x) {
#       # coefs^T [1, mean(W1), ..., mean(WJ), x]
#       sum(est_coef * c(1, Wmeans, x))
#     }
#     EY0 = EY(0)
#     EY1 = EY(1)
    
#     res[[j]] = list( EY=EY,
#                      EY0=EY0,
#                      EY1=EY1,
#                      ACE=EY1 - EY0)
#   }
  
#   return(res)
# }
