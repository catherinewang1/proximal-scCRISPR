#####################################################################################
#####            Estimating ACE/ATE with binary Exposure and continuous NCExposure using GMM
##### Includes:
#####       - 
#####################################################################################


require(EQL)


# # ============= Helper Functions
# # polybasis_transform: constructs a n x degree matrix
# # poly_eval: evaluates input vals


# #' Constructs a n x degree matrix for each Hermite Polynomial 
# #' basis j function evaluated at value i  b_j(vals[i])
# #' @param vals (numeric vector) values to evaluate polynomials at 
# #' @param degree (integer) highest degree of Hermite polynomial
# #' @param type (string) polynomial basis type (hermite, simple [1, X^2, X^3, ...])
# polybasis_transform <- function(vals, degree, type = 'hermite') {
#   switch(type,
#     hermite = {
#       inner_func <- function(d) {EQL::hermite(x=vals, n=d, prob=TRUE)}
#     },
#     simple  = {
#       inner_func <- function(d) {vals**d}
#     },
#     {print('bad [type] input (hermite, simple)'); return()})

#   return(sapply(X=0:degree, FUN=inner_func)) 
# }

# #' Evaluates the input vals at the specified Hermite Polynomial 
# #' basis degree with the specified alpha coefficients
# #' result[i] = \sum_k=1^degree \alpha[k] poly_k (val[i])
# #' @param vals (numeric vector) values to evaluate polynomials at 
# #' @param degree (integer) highest degree of Hermite polynomial
# #' @param alpha (numeric vector) coefficients
# #' @return result[i] = \sum_k=1^degree \alpha[k] poly_k (val[i])
# poly_eval <- function(vals, degree, alpha) {
#   assertthat::assert_that(length(alpha) == degree + 1, msg='alpha must be a tall vector of length == degree + 1 (for constant 1)')
#   polybasis_transform(vals, degree) %*% alpha
# }


# #' construct a n x degree matrix for cos bases
# #' y = cos( B (x + phase shift))
# #' where B is the coefficient that corresponds to the specified 
# #' period = 2 pi / B <==> B = 2 pi / period
# #' @param vals (numeric vector) values to evaluate basis fns at 
# #' @param phases (vector) of phase shifts (for basis functions)
# #' @param periods (vector) of period lengths (for basis functions)
# #' @examples
# #' v = seq(from = -5, to = 5, length.out = 100)
# #'     
# #' # changing phases
# #' vcos = cosbasis_transform(v, phases = c(0, 1, 2, 3), periods = 2*base::pi )
# #' head(vcos)
# #' plot(v, vcos[,1], type = 'l', col = "orange")
# #' lines(v, vcos[,2], col = "red")
# #' lines(v, vcos[,3], col = "purple")
# #' lines(v, vcos[,4], col = "blue")
# #' 
# #' # changing periods
# #' vcos = cosbasis_transform(v, phases = c(0), periods = 2*base::pi * c(.5, 1, 2, 4) )
# #' head(vcos)
# #' plot(v, vcos[,1], type = 'l', col = "orange")
# #' lines(v, vcos[,2], col = "red")
# #' lines(v, vcos[,3], col = "purple")
# #' lines(v, vcos[,4], col = "blue")
# cosbasis_transform <- function(vals, phases, periods) {
#    X = matrix(NA, nrow = length(vals), ncol = length(periods) * length(phases))
#    cur_col = 1
#    for(phase in phases) {
#     for(period in periods) {
#       X[,cur_col] = cos((2 * base::pi / period) * (vals + phase))
#       cur_col = cur_col + 1 # ; print(cur_col)
#     }
#    }
#    return(X)
# }



# Utils: Polynomial-Polynomial Bases (for B and H classes)
# (matching on polynomial bases B AND hypothesis class is polynomial H)


#==================== Defining functions for matching bases

#' Calculate the alpha_star (for A=1 and A=0) for the 
#' Outcome Confounding Bridge function h
#' with the given subset of Z,W,Y and 
#' with the specified degrees for B (b_degree) and D (h_degree)
#' @param B (matrix) n x # basis fns Z evaluated at basis fns
#' @param D (matrix) n x # basis fns W evaluated at basis fns
#' @param Y (vector) observed values for Outcome
#' @param Weights (matrix) b_degree x b_degree of weights
#' Solution is: 
#'    (D^T B Weights B^T D)^-1   D^T B Weights B^T Y
calcAlphaStar <- function(B, D, Y, Weights=NULL) {
  BtD = t(B) %*% D
  if(is.null(Weights)) { # if no weighting (identity)
    # sometimes cannot invert matrix (close to singular), return NA w/ msg if so
    invM = tryCatch({solve(t(BtD) %*% BtD)},
             error = function(cond) {
               message(sprintf('calcAlphaStar: failed to solve inverse at #B = %d, #D = %d (no weights)', 
                               ncol(B)-1, ncol(D)-1 )) 
               return(NA)
             }) 
    if(any(is.na(invM))) {
        return(rep(NA, ncol(D)))
    }
    
    alpha_star = invM %*% t(BtD) %*% t(B) %*% Y
  } else { # if weighting
    invM = tryCatch({solve(t(BtD) %*% Weights %*% BtD)},
             error = function(cond) {
               message(sprintf('calcAlphaStar: failed to solve inverse at #B = %d, #D = %d (with weights)', 
                               ncol(B)-1, ncol(D)-1)) 
               return(NA)
             }) 
    if(any(is.na(invM))) {
        return(rep(NA, ncol(D)))
    }
    
    alpha_star = invM %*% t(BtD) %*% Weights %*% t(B) %*% Y
  }

  
  alpha_star
}




#' Calcaulate alpha* when the regularization norm is L2 (Ridge)
#' (X^TX + lambda I)^{-1} X^T Y 
#' for special X and Ys
#' @param B (matrix) n x # basis fns Z evaluated at basis fns
#' @param D (matrix) n x # basis fns W evaluated at basis fns
#' @param Y (vector) observed values for Outcome
#' @param Weights_sqrt (matrix) b_degree x b_degree of weights square rooted
#' @param lambdas 
calcAlphaStarL2Reg <- function(B, D, Y, Weights_sqrt=NULL, lambdas=NULL) {
  # B = dataList$B0[-cv_idx0[[k]],]
  # Y = dataList$Y0[-cv_idx0[[k]]]
  # D = dataList$D0[-cv_idx0[[k]],]
  # lambda = .2
  
  if(is.null(Weights_sqrt)) {Weights_sqrt = diag(1, ncol = ncol(B))}
  if(is.null(lambdas)) {lambdas = 0}
  
  WBt =  Weights_sqrt %*% t(B) 
  Y_tilde = WBt %*% Y / nrow(B) # divide by n for scaling...
  X_tilde = WBt %*% D / nrow(B) # after other calcs for numeric?
  
  
  # alphas = list()
  alphas = matrix(NA, nrow = length(lambdas), ncol = ncol(D))
  for(i in 1:length(lambdas)) {
    lambda = lambdas[i]
    invM = tryCatch({solve(t(X_tilde) %*% X_tilde + lambda * diag(1, nrow = ncol(X_tilde)))},
                    error = function(cond) {
                      message(sprintf('calcAlphaStarL2Reg: failed to solve inverse at lambda=%.2f', 
                                      lambda )) 
                      return(NA)
                    })
    if(any(is.na(invM))) {
      alpha_star = rep(NA, ncol(D))
    } else {
      alpha_star = invM %*% t(X_tilde) %*% Y_tilde
    }
    alphas[i, ] = alpha_star
  }
  
  return(alphas)
}

#' Calc (prob out of sample/CV of k fold of alpha trained on -k folds)
#' @alphas (matrix) with shape (dim of alpha (in R^K)) x (number of lambdas)
#' @param B (matrix) n x # basis fns Z evaluated at basis fns
#' @param D (matrix) n x # basis fns W evaluated at basis fns
#' @param Y (vector) observed values for Outcome
#' @param Weights_sqrt (matrix) b_degree x b_degree of weights square rooted
calcCVErr <- function(alphas, B, D, Y, Weights_sqrt=NULL) {
  
  if(is.null(Weights_sqrt)) {Weights_sqrt = diag(1, ncol = ncol(B))}
  
  WBt =  Weights_sqrt %*% t(B) 
  Y_tilde = WBt %*% Y / nrow(B) # divide by n for scaling...
  X_tilde = WBt %*% D / nrow(B) # after other calcs for numeric?
  
  sqErrs = rep(NA, nrow(alphas))
  for(i in 1:nrow(alphas)) {
    # || Y_tilde - X_tilde alpha||_2^2
    sqErrs[i] = mean((Y_tilde - X_tilde %*% alphas[i, ])**2) 
  }
  return(sqErrs)
}

#' Estimate the ATE using the Outcome Confounding Bridge (OCB) function
#' which is estimated by matching the conditional moments
#' @param df (dataframe)
#' @param b_degree (integer) degree of 
#' @param Weights (matrix) b_degree x b_degree of weights
#' @param returnPreds (boolean) return predicted values or not
#' @param returnAlphas (boolean) return estimated alphas or not
#' estCondMomentPolyOCB(df, b_degree=1, h_degree=1)
estCondMomentPolyOCB <- function(df, b_degree, h_degree, Weights=NULL,
                                 returnPreds=FALSE, returnAlphas=FALSE) {
  # A == 1
  df1 = df |> filter(A == 1)
  alpha1 = calcAlphaStar_old(Z=df1$Z1, W=df1$W1, Y=df1$Y, 
                             b_degree=b_degree, h_degree=h_degree, Weights=Weights)
  hat_h1 = poly_eval(df$W1, h_degree, alpha1) # get estimates for ALL individuals?
  
  # A == 0
  df0 = df |> filter(A == 0)
  alpha0 = calcAlphaStar_old(Z=df0$Z1, W=df0$W1, Y=df0$Y, 
                             b_degree=b_degree, h_degree=h_degree, Weights=Weights)
  hat_h0 = poly_eval(df$W1, h_degree, alpha0) # get estimates for ALL individuals?
  
  # return
  res = list(ATE=mean(hat_h1 - hat_h0))
  if(returnPreds) { # append pred vals if asked for
    res = c(res, 
            list(hat_h1=hat_h1,
                 hat_h0=hat_h0))
  } 
  if(returnAlphas) {
    res = c(res, 
            list(alpha1=alpha1,
                 alpha0=alpha0))
  }
  return(res)
}


#' Estimate the ATE using the Outcome Confounding Bridge (OCB) function
#' which is estimated by matching the conditional moments
#' (updated: rename to OCBGMM, changed input to dataList like other fns)
#' @param dataList (list) of values containing data 
#' Need matrices/vectors/numerics (named)
#' Y1, Y0: response
#' B1, B0: conditional moment basis functions evaluated at Zs
#' D1, D0: function class basis functions evaluated at Ws
#' p1, p0: proportion of samples in each treatment group
#' # @param b_degree (integer) degree of opponent bases 
#' # @param h_degree (integer) degree of h (confounding bridge)
#' @param Weights (matrix) #b x #b (#b = ncol(dataList$B1)) of weights 
#'                         (equal weights by default)
#' @param returnPreds (boolean) return predicted values or not
#' @param returnAlphas (boolean) return estimated alphas or not
#' estCondMomentPolyOCB(df, b_degree=1, h_degree=1)
OCBGMM <- function(dataList, Weights=NULL, returnPreds=FALSE, returnAlphas=FALSE) {
  alpha1 = calcAlphaStar(B=dataList$B1, D=dataList$D1, Y=dataList$Y1, Weights=Weights)
  alpha0 = calcAlphaStar(B=dataList$B0, D=dataList$D0, Y=dataList$Y0, Weights=Weights)

  hat_h1 = rbind(dataList$D1, dataList$D0) %*% alpha1 # get estimates for ALL individuals?
  hat_h0 = rbind(dataList$D1, dataList$D0) %*% alpha0 # get estimates for ALL individuals?
    
  # return
  res = list(ATE=mean(hat_h1 - hat_h0),
             EY0=mean(hat_h0),
             EY1=mean(hat_h1))
  if(returnPreds) { # append pred vals if requested
    res = c(res, 
            list(hat_h1=hat_h1,
                 hat_h0=hat_h0))
  } 
  if(returnAlphas) { # append alphas if requested
    res = c(res, 
            list(alpha1=alpha1,
                 alpha0=alpha0))
  }
  return(res)
}


#' calculates the new 'best' weight matrix (inverse of covariance) based 
#' on previous estimate of alpha at treatment A=a
calc_weight <- function(Ya, Ba, Da, alphaa) {
  vec = (Da %*% alphaa - Ya) |> as.vector() # vector to scale each 1xJ by 
  Sigma_matrix = Ba * vec   # n x J (for estimating best weights)
  Sigma = cov(Sigma_matrix) # Covariance matrix J x J
  # Weights_a = solve(Sigma)   # Weight matrix = cov inverse
  Weights_a = tryCatch({solve(Sigma)},
             error = function(cond) {
               message(sprintf('calc_weight: failed to solve inverse of Cov (%d, %d). Just return NA for now', 
                               ncol(Ba)-1, ncol(Da)-1)) 
               # print(Sigma); return(solve(Sigma + diag(rep(1, ncol(Sigma)))))
               # print(Sigma); return(solve(diag(diag(Sigma))))
               return(NA)
             }) 
  return(Weights_a)
}



#' Estimate the ATE using the Outcome Confounding Bridge (OCB) function
#' which is estimated by matching the conditional moments
#' (One dimensional!, for now... handled by changing input df --> dataList which can be customized)
#' Rename from estCondMomentPolyOCBReweight to OCBRW
#' # @param df (dataframe)
#' @param dataList (list) of values containing data 
#' Need matrices/vectors/numerics (named)
#' Y1, Y0: response
#' B, B1, B0: conditional moment basis functions evaluated at Zs
#' D, D1, D0: function class basis functions evaluated at Ws
#' p1, p0: proportion of samples in each treatment group
#' # @param b_degree (integer) degree of opponent bases 
#' # @param h_degree (integer) degree of h (confounding bridge)
#' @param gmmSteps (integer) number of times left to perform 
#'                            estimate weighting <--> estimate alpha
#' # @param type (string) polynomial basis type (hermite, simple [1, X^2, X^3, ...])
#' @param returnWeights (boolean) return list of Weight matrices or not
OCBGMMRw <- function(dataList, gmmSteps, returnWeights=FALSE) {
  # # define values/dataframes for computing
  # df1 = df |> filter(A == 1); df0 = df |> filter(A == 0) # subset df by trtmnt
  # p1 = mean(df$A == 1); p0 = mean(df$A == 0)             # prop in each trtmnt
  # Y1 = df1$Y; Y0 = df0$Y                             # response of each trtmnt
  
  # B  = polybasis_transform( df$Z1, b_degree, type=type) # match up to _th degree
  # # B1 = polybasis_transform(df1$Z1, b_degree, type=type) # match up to _th degree
  # # B0 = polybasis_transform(df0$Z1, b_degree, type=type) # match up to _th degree
  # B1 = B[which(df$A == 1), ]; B0 = B[which(df$A == 0), ]
  
  # D  = polybasis_transform( df$W1, h_degree, type=type) # allow flexibility up to _ degrees 
  # # D1 = polybasis_transform(df1$W1, h_degree, type=type) # allow flexibility up to _ degrees 
  # # D0 = polybasis_transform(df0$W1, h_degree, type=type) # allow flexibility up to _ degrees 
  # D1 = D[which(df$A == 1), ]; D0 = D[which(df$A == 0), ]
  
  #' #' Helper fn for estCondMomentPolyOCBReWeight
  #' #' calculates the ATE, alpha1, alpha0 for a specific weighting
  #' #' @param  A (0 or 1) which treatment
  #' #' @param Weights (matrix) weighting matrix
  #' estCondMomentPolyOCBReWeight_ <- function(A, Weights) {
  #'   if(A==1) {
  #'     dfa = df1
  #'   } else if (A==0) {
  #'     dfa = df0
  #'   } else{
  #'     print('bad A arg'); return()
  #'   }
  #'   
  #'   alphaa = calcAlphaStar(Z=dfa$Z1, W=dfa$W1, Y=dfa$Y, 
  #'                          b_degree=b_degree, h_degree=h_degree, Weights=Weights)
  #'   hat_ha = poly_eval(df$W1, h_degree, alphaa) # get estimates for ALL individuals?
  #'   
  #' }
  
  # save intermediate values
  hat_EY1s=rep(NA, gmmSteps)
  hat_EY0s=rep(NA, gmmSteps)
  ATEs    = rep(NA, gmmSteps)
  alpha1s = matrix(NA, gmmSteps, ncol(dataList$D1))
  alpha0s = matrix(NA, gmmSteps, ncol(dataList$D1)) 
  if(returnWeights) {Weights_list = list()}
  
  # iteratively calc alpha <--> calc weights
  Weights = diag(rep(1, ncol(dataList$B1))) # initial weight (equal weights)
  for(step in 1:gmmSteps) {
    # print(step); print(Weights)
    # ocb = estCondMomentPolyOCB(df, b_degree, h_degree, Weights, returnAlphas = TRUE)
    alpha1 = calcAlphaStar(B=dataList$B1, D=dataList$D1, Y=dataList$Y1, Weights=Weights)
    alpha0 = calcAlphaStar(B=dataList$B0, D=dataList$D0, Y=dataList$Y0, Weights=Weights)
    
    # hat_h1 = dataList$D %*% alpha1 # get estimates for ALL individuals?
    # hat_h0 = dataList$D %*% alpha0 # get estimates for ALL individuals?
    hat_h1 = rbind(dataList$D1, dataList$D0) %*% alpha1 # get estimates for ALL individuals?
    hat_h0 = rbind(dataList$D1, dataList$D0) %*% alpha0 # get estimates for ALL individuals?
    
    hat_EY1s[step] = mean(hat_h1)
    hat_EY0s[step] = mean(hat_h0)
    ATEs[step] = mean(hat_h1 - hat_h0) # ATEs[step] = hat_EY1s[step] - hat_EY0s[step] 
    alpha1s[step, ] = alpha1
    alpha0s[step, ] = alpha0
    if(returnWeights) {Weights_list[[step]] = Weights}
    
    Weights_a1 = calc_weight(Ya=dataList$Y1, Ba=dataList$B1, Da=dataList$D1, alphaa=alpha1)
    Weights_a0 = calc_weight(Ya=dataList$Y0, Ba=dataList$B0, Da=dataList$D0, alphaa=alpha0)
    if(is.na(Weights_a1) || is.na(Weights_a0)) {break}
    Weights = Weights_a1 * dataList$p1 + Weights_a0 * dataList$p0 # avg 
  }
  
  res = list(ATE=ATEs[step], # last est of ATE
             ATEs=ATEs,
             EY1 =hat_EY1s[step],
             EY1s=hat_EY1s,
             EY0 =hat_EY0s[step],
             EY0s=hat_EY0s,
             alpha1s=alpha1s,
             alpha0s=alpha0s)
  if(returnWeights) {res$Weights = Weights_list}

  return(res)
}


#' Estimate the ATE using the Outcome Confounding Bridge (OCB) function
#' which is estimated by matching the conditional moments
#' with additional reweighting and regularization
#' (One dimensional!, for now... handled by chaning input df --> dataList which can be customized)
#' Rename from estCondMomentPolyOCBReweightReg to OCBGMMRwReg'
#' @param dataList (list) of values containing data 
#' Need matrices/vectors/numerics (named)
#' Y1, Y0: response
#' B, B1, B0: conditional moment basis functions evaluated at Zs
#' D, D1, D0: function class basis functions evaluated at Ws
#' p1, p0: proportion of samples in each treatment group
#' @param gmmSteps (integer) number of times left to perform 
#'                            estimate weighting <--> estimate alpha
#' @param K_folds (integer) number of folds for Cross Validation
#' @param returnWeights (boolean) return list of Weight matrices or not
#' @param returnLambdas (boolean) return chosen lambdas or not
OCBGMMRwReg <- function(dataList, gmmSteps, K_folds, lambdas, 
                                            returnWeights=FALSE, returnLambdas=FALSE,
                                            returnCVATEs=FALSE) {
  # Setup:  
  # save intermediate values
  ATEs    = rep(NA, gmmSteps + 1)
  alpha1s = matrix(NA, gmmSteps + 1, ncol(dataList$D1))
  alpha0s = matrix(NA, gmmSteps + 1, ncol(dataList$D1)) 
  if(returnWeights) {Weights_list = list()}
  if(returnLambdas) {cvLambdas = list(lambda0s = rep(NA, gmmSteps + 1), 
                                      lambda1s = rep(NA, gmmSteps + 1));
                     cvSqErrs = list(sqErrs1 = matrix(NA, nrow=gmmSteps, ncol=length(lambdas)),
                                     sqErrs0 = matrix(NA, nrow=gmmSteps, ncol=length(lambdas)))}
  if(returnCVATEs) {cvATEs   = matrix(nrow = gmmSteps, ncol = length(lambdas))}
  # split into K folds (split here? (more dependence, less computation) or split at each lambda calc? (less dependence, more computation))
  cv_idx0 = caret::createFolds(y = dataList$Y0, k = K_folds, list = TRUE)
  cv_idx1 = caret::createFolds(y = dataList$Y1, k = K_folds, list = TRUE)
  

  # Step 0: Initialize Omega=I, lambda=0
  Weights = diag(rep(1, ncol(dataList$B1))) # initial weight (equal weights)
  lambda0 = 0; lambda1 = 0 # separate regularization params
  # lambda0 = 10^(-1); lambda1 = 10^(-1) # start with small regularization param
  
  # Step 1: Estimate alpha = alpha(Omega, lambda=0)
  alpha1 = calcAlphaStar(B=dataList$B1, D=dataList$D1, Y=dataList$Y1, Weights=Weights)
  alpha0 = calcAlphaStar(B=dataList$B0, D=dataList$D0, Y=dataList$Y0, Weights=Weights)
  
  
  # if(anyNA(alpha1)) { # if did not succeed to calc alpha with lambdas=0, do with small lambda=.1
  #   lambda1 = .1
  #   alpha1 = calcAlphaStarL2Reg(B = dataList$B1, D = dataList$D1, Y = dataList$Y1,
  #                               Weights_sqrt = Weights, lambdas = lambda1)} 
  # if(anyNA(alpha0)) {
  #   lambda0 = .1
  #   alpha0 = calcAlphaStarL2Reg(B = dataList$B0, D = dataList$D0, Y = dataList$Y0,
  #                               Weights_sqrt = Weights, lambdas = lambda0)
  #   
  # }
  
  hat_h1 = rbind(dataList$D1, dataList$D0) %*% alpha1 # get estimates for ALL individuals?
  hat_h0 = rbind(dataList$D1, dataList$D0) %*% alpha0 # get estimates for ALL individuals?
  
  # save intermediate values
  ATEs[1] = mean(hat_h1 - hat_h0)
  alpha1s[1, ] = alpha1
  alpha0s[1, ] = alpha0
  if(returnWeights) {Weights_list[[1]] = Weights}
  if(returnLambdas) {cvLambdas$lambda1s[1] = lambda1;
                     cvLambdas$lambda0s[1] = lambda1}
  
  # repeat step 2-3 
  for(step in 2:(gmmSteps + 1)) {
    # Step 2: Reweight Omega = REWEIGHT(alpha, Data)
      Weights_a1 = calc_weight(Ya=dataList$Y1, Ba=dataList$B1, Da=dataList$D1, alphaa=alpha1)
      Weights_a0 = calc_weight(Ya=dataList$Y0, Ba=dataList$B0, Da=dataList$D0, alphaa=alpha0)
      if(any(is.na(Weights_a1)) && any(is.na(Weights_a0))) {         # both NA
        step = step - 1 # (for indexing most recent ATE)
        break
      } else if(any(is.na(Weights_a1)) && !any(is.na(Weights_a0))) { # 1 is NA, 0 is OK
        Weights = Weights_a0
      } else if(!any(is.na(Weights_a1)) && any(is.na(Weights_a0))) { # 1 is OK, 0 is NA
        Weights = Weights_a1
      } else {                                                       # both OK
        # avg
        Weights = Weights_a1 * dataList$p1 + Weights_a0 * dataList$p0 
      }
      # if(any(is.na(Weights_a1)) && any(is.na(Weights_a0))) {break} # CHANGED: ok if 1 is NA
      # Weights = Weights_a1 * dataList$p1 + Weights_a0 * dataList$p0 # avg 
    
    # Step 3: CV to choose lambda (--> also gets alpha)
    Weights_sqrt = expm::sqrtm(Weights)
    sqErrs1 = rep(0, length(lambdas))
    sqErrs0 = rep(0, length(lambdas))
    for(k in 1:K_folds) {
      # For A = 1:
      # train on (-k)th folds
      alphas1 = calcAlphaStarL2Reg(B = dataList$B1[-cv_idx1[[k]],], 
                                   D = dataList$D1[-cv_idx1[[k]],],
                                   Y = dataList$Y1[-cv_idx1[[k]]],
                                   Weights_sqrt = Weights_sqrt,
                                   lambdas = lambdas)
      # calculate err on (k)th fold
      sqErrs1_k = calcCVErr(alphas = alphas1, 
                            B = dataList$B1[cv_idx1[[k]],], 
                            D = dataList$D1[cv_idx1[[k]],],
                            Y = dataList$Y1[cv_idx1[[k]]],
                            Weights_sqrt = Weights_sqrt)
      sqErrs1 = sqErrs1 + sqErrs1_k
      
      # For A = 0:
      # train on (-k)th folds
      alphas0 = calcAlphaStarL2Reg(B = dataList$B0[-cv_idx0[[k]],], 
                                   D = dataList$D0[-cv_idx0[[k]],],
                                   Y = dataList$Y0[-cv_idx0[[k]]],
                                   Weights_sqrt = Weights_sqrt,
                                   lambdas = lambdas)
      # calculate err on (k)th fold
      sqErrs0_k = calcCVErr(alphas = alphas0, 
                            B = dataList$B0[cv_idx0[[k]],], 
                            D = dataList$D0[cv_idx0[[k]],],
                            Y = dataList$Y0[cv_idx0[[k]]],
                            Weights_sqrt = Weights_sqrt)
      # print(sqErrs0_k); plot(log(lambdas), sqErrs0_k, type = 'l')
      sqErrs0 = sqErrs0 + sqErrs0_k
    }
    lambda1 = lambdas[which.min(sqErrs1)]
    lambda0 = lambdas[which.min(sqErrs0)]
    

    alpha1 = calcAlphaStarL2Reg(B = dataList$B1, 
                                D = dataList$D1,
                                Y = dataList$Y1,
                                Weights_sqrt = Weights_sqrt,
                                lambdas = lambda1) |> t() # transpose to be tall
    alpha0 = calcAlphaStarL2Reg(B = dataList$B0, 
                                D = dataList$D0,
                                Y = dataList$Y0,
                                Weights_sqrt = Weights_sqrt,
                                lambdas = lambda0) |> t() # transpose to be tall
    # alpha1 =  alphas1[which.min(sqErrs1), ]
    # alpha0 =  alphas0[which.min(sqErrs0), ]
    
    hat_h1 = rbind(dataList$D1, dataList$D0) %*% alpha1 # get estimates for ALL individuals?
    hat_h0 = rbind(dataList$D1, dataList$D0) %*% alpha0 # get estimates for ALL individuals?
    
    
    
    # save intermediate values
    ATEs[step] = mean(hat_h1 - hat_h0)
    alpha1s[step, ] = alpha1
    alpha0s[step, ] = alpha0
    if(returnWeights) {Weights_list[[step]] = Weights}
    if(returnLambdas) {cvLambdas$lambda1s[step] = lambda1;
                       cvLambdas$lambda0s[step] = lambda0;
                       cvSqErrs$sqErrs1[step - 1, ] = sqErrs1;
                       cvSqErrs$sqErrs0[step - 1, ] = sqErrs0}
    
    # if save intermediate cv ATEs (more computation)
    if(returnCVATEs) {
      alpha1CV = calcAlphaStarL2Reg(B = dataList$B1, 
                                    D = dataList$D1,
                                    Y = dataList$Y1,
                                    Weights_sqrt = Weights_sqrt,
                                    lambdas = lambdas) |> t() # transpose to be tall
      alpha0CV = calcAlphaStarL2Reg(B = dataList$B0, 
                                    D = dataList$D0,
                                    Y = dataList$Y0,
                                    Weights_sqrt = Weights_sqrt,
                                    lambdas = lambdas) |> t() # transpose to be tall
      
      hat_h1 = rbind(dataList$D1, dataList$D0) %*% alpha1CV # get estimates for ALL individuals?
      hat_h0 = rbind(dataList$D1, dataList$D0) %*% alpha0CV # get estimates for ALL individuals?
      cvATEs[step - 1, ] = colMeans(hat_h1) - colMeans(hat_h0)
    }
  }
  
  # assemble results together
  res = list(ATE=ATEs[step], # last est of ATE
             ATEs=ATEs,
             alpha1s=alpha1s,
             alpha0s=alpha0s)
  if(returnWeights) {res$Weights   = Weights_list}
  if(returnLambdas) {res$cvLambdas = cvLambdas;
                     res$cvSqErrs  = cvSqErrs}
  if(returnCVATEs)  {res$cvATEs    = cvATEs}
  return(res)
}



# #' Constructing Data List for Poly OCB variation 1
# #' Transform each Zk and Wj with each basis function b's and d's respectively
# #' with NO interactions
# #' 
# #' B and D's are centered and scaled (except for intercept)
# #' Need matrices/vectors/numerics
# #' Y1, Y0: response
# #' B, B1, B0: conditional moment basis functions evaluated at Zs
# #' D, D1, D0: function class basis functions evaluated at Ws
# #' p1, p0: proportion of samples in each treatment group
# #' @param df (dataframe)
# #' @param b_degree (integer) max degree of conditional basis functions
# #' @param h_degree (integer) max degree of h (confounding bridge) basis fns
# #' @param type (string) type of basis functions (e.g. hermite, simple, )
# constructDataListv1 <- function(df, b_degree, h_degree, type='hermite') {
#   Znames = grep('Z', colnames(df), value = TRUE); numZ = length(Znames)
#   Wnames = grep('W', colnames(df), value = TRUE); numW = length(Wnames)
  
#   # define values/dataframes for computing
#   df1 = df |> filter(A == 1); df0 = df |> filter(A == 0) # subset df by trtmnt
#   p1 = mean(df$A == 1); p0 = mean(df$A == 0)             # prop in each trtmnt
#   Y1 = df1$Y; Y0 = df0$Y                             # response of each trtmnt
  
#   # Construct B (basis b's evaluated at Zs)
#   B = matrix(NA, nrow = nrow(df), ncol = 1 + length(Znames) * b_degree) 
#   B[, 1] = 1
#   for(j in 1:length(Znames)) {
#     B_zname = polybasis_transform(df[, Znames[j]], b_degree, type=type) # match up to _th degree
#     # print(((j-1)*b_degree + 2):(j*b_degree + 1)) # this is what happens when indexing by 1 and inclusive of last point in range...
#     # B[,((j-1)*b_degree + 2):(j*b_degree + 1)] = B_zname[,-1] # add cols to B
#     B[,((j-1)*b_degree + 2):(j*b_degree + 1)] = scale(B_zname[,-1], center = TRUE, scale = TRUE) # add cols to B (centered&scaled)
#   }
#   # B[,2:ncol(B)] = scale(B[,2:ncol(B)], center = TRUE, scale = TRUE) # center and standardize each column (except 1)
#   B1 = B[which(df$A == 1), ]; B0 = B[which(df$A == 0), ]
  
#   # D  = polybasis_transform( df$W1, h_degree, type=type) # allow flexibility up to _ degrees 
#   # D1 = polybasis_transform(df1$W1, h_degree, type=type) # allow flexibility up to _ degrees 
#   # D0 = polybasis_transform(df0$W1, h_degree, type=type) # allow flexibility up to _ degrees 
  
#   D = matrix(NA, nrow = nrow(df), ncol = 1 + length(Wnames) * h_degree) 
#   D[, 1] = 1
#   for(j in 1:length(Wnames)) {
#     D_zname = polybasis_transform(df[, Wnames[j]], h_degree, type=type) # match up to _th degree
#     # print(((j-1)*h_degree + 2):(j*h_degree + 1)) # this is what happens when indexing by 1 and inclusive of last point in range...
#     # D[,((j-1)*h_degree + 2):(j*h_degree + 1)] = D_zname[,-1] # add cols to B
#     D[,((j-1)*h_degree + 2):(j*h_degree + 1)] = scale(D_zname[,-1], center = TRUE, scale = TRUE)  # add cols to B
#   }
#   # D[,2:ncol(D)] = scale(D[,2:ncol(D)], center = TRUE, scale = TRUE)
#   D1 = D[which(df$A == 1), ]; D0 = D[which(df$A == 0), ]
  
#   return(list(Y1=Y1, Y0=Y0,
#               # B=B, 
#               B1=B1, B0=B0,
#               # D=D, 
#               D1=D1, D0=D0, 
#               p1=p1, p0=p0))
# }



# #' Constructing Data List for Poly OCB variation 2 (with Cos bases)
# #' Transform each Zk and Wj with each basis function b's and d's respectively
# #' with NO interactions
# #' 
# #' B and D's are centered and scaled (except for intercept)
# #' Need matrices/vectors/numerics
# #' Y1, Y0: response
# #' B, B1, B0: conditional moment basis functions evaluated at Zs
# #' D, D1, D0: function class basis functions evaluated at Ws
# #' p1, p0: proportion of samples in each treatment group
# #' @param df (dataframe)
# #' @param b_degree (integer) max degree of conditional basis functions
# #' @param h_degree (integer) max degree of h (confounding bridge) basis fns
# #' @param b_phases  (vector) of phase shifts (numeric) for b basis functions
# #' @param b_periods (vector) of period lengths (numeric) for b basis functions
# #' @param h_phases  (vector) of phase shifts (numeric) for d basis functions
# #' @param h_periods (vector) of period lengths (numeric) for d basis functions
# #' @param type (string) type of poly basis functions (e.g. hermite, simple, )
# constructDataListv2 <- function(df, 
#                                 b_degree=0, h_degree=0, 
#                                 b_phases=NULL, b_periods=NULL,
#                                 h_phases=NULL, h_periods=NULL, 
#                                 type='hermite') {
#   Znames = grep('Z', colnames(df), value = TRUE); numZ = length(Znames)
#   Wnames = grep('W', colnames(df), value = TRUE); numW = length(Wnames)
  
#   # define values/dataframes for computing
#   df1 = df |> filter(A == 1); df0 = df |> filter(A == 0) # subset df by trtmnt
#   p1 = mean(df$A == 1); p0 = mean(df$A == 0)             # prop in each trtmnt
#   Y1 = df1$Y; Y0 = df0$Y                             # response of each trtmnt
  
#   B = matrix(NA, nrow = nrow(df), ncol = 1) 
#   B[, 1] = 1
#   # cur_col = 2
#   for(j in 1:length(Znames)) {
#     B_zname_poly = NULL
#     B_zname_cos  = NULL
    
#     if(b_degree > 0) {
#       B_zname_poly = polybasis_transform(df[, Znames[j]], degree=b_degree, type=type)[,-1]
#     }
    
#     if(!is.null(b_phases) && !is.null(b_periods)) {
#       B_zname_cos = cosbasis_transform(vals=df[, Znames[j]], phases=b_phases, periods=b_periods)
      
#     }
#     B = cbind(B, B_zname_poly, B_zname_cos)
#     # B = cbind(B, scale(cbind(B_zname_poly, B_zname_cos), center = TRUE, scale = TRUE))
#     # B[,cur_col:(ncol(B_zname) + cur_col - 2)] = scale(B_zname[,-1], center = TRUE, scale = TRUE)
#     # B[,((j-1)*b_degree + 2):(j*b_degree + 1)] = scale(B_zname[,-1], center = TRUE, scale = TRUE) # add cols to B (centered&scaled)
#     # cur_col = cur_col + ncol(B_zname)
#   }
#   B[,-1] = scale(B[,-1], center=TRUE, scale=TRUE) 
#   B1 = B[which(df$A == 1), ]; B0 = B[which(df$A == 0), ]
  
#   # Construct D (basis d's evaluated at Ws)
#   D = matrix(NA, nrow = nrow(df), ncol = 1) 
#   D[, 1] = 1
#   for(j in 1:length(Wnames)) {
#     D_zname_poly = NULL; D_zname_cos = NULL
#     if(h_degree > 0) {
#       D_zname_poly = polybasis_transform(df[, Wnames[j]], degree=h_degree, type=type)[,-1] 
#     }

#     if(!is.null(h_phases) && !is.null(h_periods)) {
#       D_zname_cos = cosbasis_transform(vals=df[, Wnames[j]], phases=h_phases, periods=h_periods)
#     }
#     D = cbind(D, D_zname_poly, D_zname_cos)
#   }
#   D[,-1] = scale(D[,-1], center=TRUE, scale=TRUE)
#   D1 = D[which(df$A == 1), ]; D0 = D[which(df$A == 0), ]
  
#   return(list(Y1=Y1, Y0=Y0,
#               B=B,
#               B1=B1, B0=B0,
#               D=D,
#               D1=D1, D0=D0, 
#               p1=p1, p0=p0))
# }


#' Constructing dataframe out of datalist format 
#' (GMM fns use datalist format, OCB2SLS use dataframe)
#' 
#' useful for getting basis transformations
#' @param dataList (list) of values containing data 
#' Need matrices/vectors/numerics (named)
#' Y1, Y0: response
#' B, B1, B0: conditional moment basis functions evaluated at Zs
#' D, D1, D0: function class basis functions evaluated at Ws
#' p1, p0: proportion of samples in each treatment group
#' # @param b_degree (integer) degree of opponent bases 
#' # @param h_degree (integer) degree of h (confounding bridge)
#' @return a dataframe with variables Y, W, A, Z
#' @examples
#' dl = constructDataListv2(df, b_degree=4, h_degree=2, type = 'simple',
#'                         b_phases=c(0), b_periods=c(4, 8, 16, 24),
#'                         h_phases=c(0), h_periods=c(4, 8, 16, 24))
#' df_bases = get_df_from_dl(dl)
#' OCB_2SLSReg_CosPoly_res = OCB2SLSReg(df_bases, alpha=.1)
get_df_from_dl <- function(dl) {
  # reformat a dataframe for OCB2SLS using new bases (poly and cos)
  B = rbind(dl$B0, dl$B1)
  B = B[, -1] # remove first col of 1s
  colnames(B) = paste0('Z', 1:ncol(B))
  D = rbind(dl$D0, dl$D1)
  D = D[, -1] # remove first col of 1s
  colnames(D) = paste0('W', 1:ncol(D))
  
  A = data.frame(A = c(rep(0, nrow(dl$B0)), rep(1, nrow(dl$B1))))
  Y = data.frame(Y = c(dl$Y0, dl$Y1))
  df_bases = cbind(Y, A, B, D)
  return(df_bases)
}





#' Get the estimated Inverse of M = E(  b(Z) d(W)^T  )
#' separate function, so that each can be estimated on each subset of data?
getMhatinv <- function(B, D) {
  # --------- calculate Mhat where M = E(  b(Z) d(W)^T)
  # ones = matrix(1, nrow=1, ncol=ncol(D))
  # mySum = matrix(0, nrow = ncol(B), ncol = ncol(D))
  # for(i in 1:nrow(B)) {
  #   mySum = mySum + 
  #            # Omega_half %*% # cancels out w g(O; ) term: (AB)^-1= B^-1 A^-1 
  #           t(B[i, , drop=FALSE]) %*% 
  #           D[i, , drop=FALSE] # (ones * Y[i] - D[i, , drop=FALSE]) # should only be b(Z)k d(W)j?
  #           
  # }
  # Mhat0 = (1/nrow(B)) * mySum 
  
  # faster way to calc
  Mhat = (t(B) %*% D) / (nrow(B))
  
  
  
  # Mhat_inv = solve(Mhat)
  # change to pseudoinverse? (Moore-Penrose)
  # Mhat_inv = corpcor::pseudoinverse(Mhat)
  Mhat_svd_res = svd(Mhat) # get diag/sing vals to get 'importance', then we remove <.05
  Mhat_inv = corpcor::pseudoinverse(Mhat, tol = sum(Mhat_svd_res$d) * .05) # perhaps try to stabilize inverse calc
  
  return(Mhat_inv)
}

#' Get the values f(\mathcal{O};\alpha) for linear OCB and m-estimator for \alpha
#' 
#' @description 
#' Get the values f(\mathcal{O};\alpha) which characterizes the variance of the 
#' final estimate of the mean potential outcome when the \alpha parameter is
#' estimated using m-estimation
#' This uses the specific linear confounding bridge form h(W;\alpha) = d(W)^T\alpha
#' and our estimating equations E(g(\mathcal{O};\alpha)) = 0
#' for g = \sqrt{Omega} b(Z)(Y - h(W))
#' 
#' These fO may be used to calculate the variance of the mean potential outcome
#' at specific treatment level $A=a$ or combined with the outputs from 
#' A=1 and A=0 for the variance of the ATE estimate (take Var_n(fO1 - fO0))
#' 
#' updated equation: no longer need the Omega matrix because terms cancel out
#' @param alpha (vector) parameterizing h(W;\alpha) ocb function
#' @param Y (vector) of Y responses
#' @param B (matrix) of basis transformations of Z where B_{ij} = b_j(Z_i)
#' @param D (matrix) of basis transformations of W where D_{ik} = d_k(W_i)
#' @param Omega (matrix) square/symmetric/PSD weight matrix of dimension JxJ
#' 
#' @examples
#'  
#' gammaSetting = 'A1'; basis = 'basis1' 
#' GMM_steps = 5
#' # gammaSetting = 'A2'
#' df = sim_data_multidimU(N = N,
#'                         U_mean = gammas[[gammaSetting]]$U_mean,
#'                         U_sd   = gammas[[gammaSetting]]$U_sd,
#'                         beta   = gammas[[gammaSetting]]$beta,
#'                         gammaA = gammas[[gammaSetting]]$gammaA,
#'                         gammaY = gammas[[gammaSetting]]$gammaY,
#'                         gammaZ = gammas[[gammaSetting]]$gammaZ,
#'                         gammaW = gammas[[gammaSetting]]$gammaW,
#'                         muA    = gammas[[gammaSetting]]$muA)
#' # df_all = data.frame(Y = dat$Y, A = dat$A, dat$Z, dat$W, dat$U)
#' 
#' dl = constructDataListv2(df,
#'                          b_degree=basisParams[[basis]]$b_degree, 
#'                          h_degree=basisParams[[basis]]$h_degree, type = 'simple',
#'                          b_phases=basisParams[[basis]]$b_phases, 
#'                          b_periods=basisParams[[basis]]$b_periods,
#'                          h_phases=basisParams[[basis]]$h_phases, 
#'                          h_periods=basisParams[[basis]]$h_periods)
#' 
#' # get inputs
#' # OCBGMM ReWeight
#' OCBMEstRes = OCBGMMRw(dl, GMM_steps, returnWeights = T)
#' f_hat_res = getFhat(alpha = t(OCBMEstRes$alpha1s[GMM_steps, , drop=FALSE]),
#'                   Y = dl$Y1, #+ rnorm(n = length(dl$Y1)),
#'                   B = dl$B1,
#'                   D = dl$D1,
#'                   Omega = OCBMEstRes$Weights[[GMM_steps]] ) 
#' 
#' V2_hat_res = mean((f_hat_res - mean(f_hat_res))**2)
#' 1.96 * (sqrt(V2_hat_res) / sqrt(N))
#' 
#' 
getFhat <- function(alpha, Y, B, D, Mhat_inv) { #, Omega=NULL) {
  # square root of Omega weight matrix
  # Omega_half = expm::sqrtm(Omega) #  Omega_half %*% Omega_half
  
  
  
  
  
  
  # --------- Calc f_hat d(W)^T alpha + d(W)^T M b(Z)(Y - d(W)^T alpha)
  # slower way to calc f_hat
  # f_hat = rep(NA, length = nrow(B))
  # for(i in 1:nrow(B)) {
  #   f_hat[i] =
  #     D[i, , drop=FALSE] %*% (
  #       alpha
  #       +
  #         Mhat_inv %*% 
  #         # Omega_half %*% 
  #         t(B[i, , drop=FALSE]) * (Y[i] - as.numeric(D[i, , drop=FALSE] %*% alpha))
  #     )
  # }
  
  # another way to calc f_hat, should be quicker
  # predictions of Y: d(W)^T alpha 
  # f(O) = Ypreds + d(W) Dhatinv b(Z) Ypreds
  Ypreds = D %*% alpha |> as.vector()
  # # plot of residuals
  # hist(Y - Ypreds) # at least they are centered...
  D_Mhatinv = D %*% Mhat_inv # %*% Omega_half
  DMB = rep(NA, length = nrow(B))
  for(i in 1:nrow(B)) {
    DMB[i] = D_Mhatinv[i, , drop=FALSE]  %*% t(B[i, , drop=FALSE])
  }
  
  f_hat2 = Ypreds + DMB * (Y - Ypreds)
  
  # check both ways are the same
  # mean((f_hat2 - mean(f_hat2))**2)
  # mean((f_hat - mean(f_hat))**2)
  # f_hat[1:5]
  # f_hat2[1:5]
  return(f_hat2)
}

#' Get the Fhat of the ATE when stacking alpha = c(alpha1, alpha0)
getFhatATE <- function(dl, alpha1, alpha0, Mhat_inv1, Mhat_inv0) {
  # debugging
  # alpha = rbind(alpha1, alpha0)
  # Y = Y10
  # B = rbind(cbind(dl$B1,  matrix(0, nrow = nrow(dl$B1), ncol = ncol(dl$B0))),
  #           cbind(matrix(0, nrow = nrow(dl$B0), ncol = ncol(dl$B1)), dl$B0))
  # D = cbind(rbind(dl$D1, dl$D0), -rbind(dl$D1, dl$D0))
  # Mhat_inv = rbind(cbind(Mhat_inv1, matrix(0, nrow = nrow(Mhat_inv1), ncol = ncol(Mhat_inv0))),
  #                  cbind(matrix(0, nrow = nrow(Mhat_inv0), ncol = ncol(Mhat_inv1)), Mhat_inv0))
  
  alpha = rbind(alpha1, alpha0)
  Mhat_inv = rbind(cbind(Mhat_inv1, 
                         matrix(0, nrow = nrow(Mhat_inv1), ncol = ncol(Mhat_inv0))),
                  cbind(matrix(0, nrow = nrow(Mhat_inv0), ncol = ncol(Mhat_inv1)), 
                        Mhat_inv0))
  Y = c(dl$Y1, dl$Y0)
  # [B A, B(1-A)] = bind together observed B with 0's for I(A=a) parts
  B = rbind(cbind(dl$B1,  
                  matrix(0, nrow = nrow(dl$B1), ncol = ncol(dl$B0))),
            cbind(matrix(0, nrow = nrow(dl$B0), ncol = ncol(dl$B1)), 
                  dl$B0))
  # A = c(rep(1, times = nrow(dl$B1)), 
  #       rep(0, times = nrow(dl$B0)))
  
  # DH = derivative of H = derivative of h1 - h0 = [d(W) -d(W)]
  DH = cbind( rbind(dl$D1, dl$D0), 
             -rbind(dl$D1, dl$D0))
  # [D A, D(1-A)] = bind together observed D with 0's for I(A=a) parts
  D_EY10 = rbind(cbind(dl$D1,  matrix(0, nrow = nrow(dl$D1), ncol = ncol(dl$D0))),
                 cbind(matrix(0, nrow = nrow(dl$D0), ncol = ncol(dl$D1)), dl$D0))
  
  
  
  # h(W;alpha1) - h(W;alpha0)
  Y1minusY0preds = DH %*% alpha
  
  
  vec = rep(NA, times = length(Y))
  for(i in 1:length(vec)) {
    vec[i] = DH[i, ] %*%
             Mhat_inv %*% 
             t(B[i, , drop=FALSE]) * as.numeric(Y[i]  -  D_EY10[i, , drop=FALSE] %*% alpha)
  }
  
  
  FhatATE = Y1minusY0preds + vec
  # # what if minus?
  # FhatATE = Y1minusY0preds - vec
  return(FhatATE)
}


#' check empirical estimating equation's solution (how close to 0)
checkEmpEstEq <- function(Y, B, D, alpha) {
  # B = dl$B0
  # Y = dl$Y0
  # D = dl$D0
  # alpha = OCBGMM_res$alpha0
  
  mySum = matrix(0, ncol=1, nrow = ncol(B))
  for(i in 1:nrow(B)) {
    mySum = mySum + 
      t(B[i, , drop=FALSE]) * as.numeric(Y[i] - D[i, ] %*% alpha)
  }
  mySum = mySum / nrow(B)
  return(mySum)  
}



#' Calculating the variance of M-estimator (+ construct CI)
#' 
#' from data list (dl), fits GMM/M-estimation based on GMM 
#' parameters such as basis and GMM_steps, constructs CI with specified coverage,
#' and checks if ATE_true lies within constructed CI if ATE_true is not NULL
#' @param dl (list) of data objects used in GMM estimation function
#' @param alpha1 (vector) of alpha estimate A=1
#' @param alpha0 (vector) of alpha estimate A=0
#' @param Omega (matrx) weight matrix PD
#' @param coverage (numeric) coverage proportion (eg 95) for the CIs
#' @param ATE_est (numeric) the estimated ATE
#' @param ATE_true (numeric) the true ATE (e.g. beta=2)
#' @param returnVar (boolean) whether or not to return the estimated variance
checkCoverage <- function(dl, alpha1, alpha0, Omega, coverage, ATE_est, ATE_true=NULL, returnVar=FALSE) {
  # df = sim_data_multidimU(N = N,
  #                         U_mean = gammas[[gammaSetting]]$U_mean,
  #                         U_sd   = gammas[[gammaSetting]]$U_sd,
  #                         beta   = gammas[[gammaSetting]]$beta,
  #                         gammaA = gammas[[gammaSetting]]$gammaA,
  #                         gammaY = gammas[[gammaSetting]]$gammaY,
  #                         gammaZ = gammas[[gammaSetting]]$gammaZ,
  #                         gammaW = gammas[[gammaSetting]]$gammaW,
  #                         muA    = gammas[[gammaSetting]]$muA)
  # # df_all = data.frame(Y = dat$Y, A = dat$A, dat$Z, dat$W, dat$U)
  # 
  # dl = constructDataListv2(df,
  #                          b_degree=basisParams[[basis]]$b_degree, 
  #                          h_degree=basisParams[[basis]]$h_degree, type = 'simple',
  #                          b_phases=basisParams[[basis]]$b_phases, 
  #                          b_periods=basisParams[[basis]]$b_periods,
  #                          h_phases=basisParams[[basis]]$h_phases, 
  #                          h_periods=basisParams[[basis]]$h_periods)
  
  
  
  # OCBGMM ReWeight
  # OCBMEstRes = OCBGMMRw(dl, GMM_steps, returnWeights=T)
  
  
  # estimate of variance, we can take over the entire dataset, just change the alpha's?
  # f_hat1 = getFO(alpha = t(OCBMEstRes$alpha1s[GMM_steps, , drop=FALSE]),
  #                Y = dl$Y1, #+ rnorm(n = length(dl$Y1)),
  #                B = dl$B1,
  #                D = dl$D1,
  #                Omega = OCBMEstRes$Weights[[GMM_steps]] )
  # 
  # f_hat0 = getFO(alpha = t(OCBMEstRes$alpha0s[GMM_steps, , drop=FALSE]),
  #                Y = dl$Y0, #+ rnorm(n = length(dl$Y1)),
  #                B = dl$B0,
  #                D = dl$D0,
  #                Omega = OCBMEstRes$Weights[[GMM_steps]] )
  
  # alpha0 = t(OCBMEstRes$alpha0s[GMM_steps, , drop=FALSE])
  # alpha1 = t(OCBMEstRes$alpha1s[GMM_steps, , drop=FALSE])
  # Omega = OCBMEstRes$Weights[[GMM_steps]]
  
  
  Y10 = c(dl$Y1, dl$Y0)
  B10 = rbind(dl$B1, dl$B0)
  D10 = rbind(dl$D1, dl$D0)
  
  Mhat_inv1 = getMhatinv(B = dl$B1, D = dl$D1)
  Mhat_inv0 = getMhatinv(B = dl$B0, D = dl$D0)
  
  f_hat1 = getFhat(alpha = alpha1,
                   Y = Y10, 
                   B = B10,
                   D = D10, 
                   Mhat_inv = Mhat_inv1
                   # Omega =  Omega
                   )
  
  f_hat0 = getFhat(alpha = alpha0,
                   Y = Y10, 
                   B = B10,
                   D = D10,
                   Mhat_inv = Mhat_inv0
                   # Omega = Omega
                   )
  
  
  
  
  f_hatATE = getFhatATE(dl=dl, 
             alpha1=alpha1, 
             alpha0=alpha0, 
             Mhat_inv1=Mhat_inv1, 
             Mhat_inv0=Mhat_inv0)
  
  # variance of hat ATE ? one of these?
  V2_ATE = mean((f_hatATE - mean(f_hatATE))**2)
  # V2_ATE = mean(((f_hat1-f_hat0) - mean(f_hat1-f_hat0))**2)
  
  
  # Construct CI
  
  # testing coverage at .95
  z_coverage = qnorm(1 - (1 - coverage)/2, mean = 0, sd = 1)
  # OCBMEstRes$ATE + 1.96 * (sqrt(V2_ATE) / sqrt(N))
  Nsamples = length(Y10) # ??? different sample sizes for each 
  CI  = ATE_est + c(-1, 1) * z_coverage * (sqrt(V2_ATE) / sqrt(Nsamples))
  res = list(CI=CI)
  
  
  # is the true ATE (e.g.=2) in the CI
  if(!is.null(ATE_true)) {
    res$covered = (CI[1] <= ATE_true) && (ATE_true <= CI[2])
  } else {
    res$covered = NULL
  }
  
  if(returnVar) {
    # variance of hat EY1 
    V2_EY1 = mean((f_hat1 - mean(f_hat1))**2)
    # variance of hat EY0 
    V2_EY0 = mean((f_hat0 - mean(f_hat0))**2)
    
    res$VarEY1 = V2_EY1
    res$VarEY0 = V2_EY0
    res$VarATE = V2_ATE
  }
  
  return(res)
  
}