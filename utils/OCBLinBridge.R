
# Assuming linear bridge function
# h(W) = d(W)^T alpha
# 


########################################################################
# Helper function:
# 
########################################################################

#' For Plug-In and One-Step functions, need 
#' matrix of predicted values dW E(dWi|Z,A)
#' and vector mu (E(Y|Z,A)).
#' This function returns the matrix of predicted E(dWi|Z,A)
#' Assumes separated A's already! e.g. input is Z |> filter(A==a)
#' @param D matrix of dWi's, we want to estimate E(dWi|Z,A)
#' @param Z conditioning on Zs (can include transforms if wanted)
#' @output thetaHatMatrix (matrix of estimated E(dWi|Z,A))
getInputTheta <- function(D, Z) {
  # D = dl$D0
  # Z = Z0
  
  # add intercept to Z if not already there
  if(!all(Z[,1] == 1)) {Z = cbind(rep(1, nrow(Z)),Z)}
  
  # matrix of predicted values of dW
  thetaHatMatrix = matrix(NA, nrow = nrow(D), ncol = ncol(D))
  
  i = 1 # current col 
  # fill constant col, if D had constant col
  if(all(D[,1] == 1)) {
    thetaHatMatrix[, 1] = 1
    i = i + 1
  } 
  
  # fill out rest of columns dWi's
  while(i <= ncol(D)) {
    fit = glm.fit(y = D[,i], x = Z, family = gaussian(), intercept = TRUE)
    thetaHatMatrix[, i] = fit$fitted.values
    i = i + 1
  }
  
  return(thetaHatMatrix)
  
}

#' For Plug-In and One-Step functions, need 
#' matrix of predicted values dW E(dWi|Z,A)
#' and vector predicted mu (E(Y|Z,A)).
#' This function returns the vector of predicted E(Y|Z,A)
#' Assumes separated A's already! e.g. input is Z |> filter(A==a)
#' @param Y vector of Y's, we want to estimate E(Y|Z,A)
#' @param Z conditioning on Zs (can include transforms if wanted)
#' @output muVector (matrix of estimated E(Y|Z,A))
getInputMu <- function(Y, Z) {
  # Y = dl$Y0
  # Z = Z0
  # add intercept to Z if not already there
  if(!all(Z[,1] == 1)) {Z = cbind(rep(1, nrow(Z)),Z)}
  
  fit = glm.fit(x = Z, y = Y, family = gaussian())
  muVector = fit$fitted.values
  return(muVector)
}


#' P(A=a_i|Z) using Logistic Regression
#' and trim estimated probs to be \in [trim, 1 - trim]
#' Note: P(A=a_i|Z) and NOT P(A=1|Z)
#' this is the probability of the observed a_i
#' @param A (vector) of length n of binary 0/1 indicating treatment
#' @param Z (matrix) of dim n x (#Z or further transf of Z)
#'                   of Z values (or further transf of Z)
#' @param trim (numeric) in [0,1] to trim the estimated probabilities
#'        p \in [0,1] --> p \in [trim, 1 - trim]
#' @example 
#' testPAZ = getInputPaz(A = c(rep(0, nrow(Z0)), 
#'                        rep(1, nrow(Z1))),
#'                  Z = rbind(Z0, Z1))
#' testPAZ0 = testPAZ[1:nrow(Z0)]
#' testPAZ1 = testPAZ[-(1:nrow(Z0))]
getInputPaz <- function(A, Z, trim=0) {
  
  # Z = rbind(Z0, Z1)
  # A = c(rep(0, nrow(Z0)), rep(1, nrow(Z1)))
  
  logfit = glm.fit(y=A, x=Z, family=binomial())
  p1z = logfit$fitted.values      # probability of P(A=1|Z)
  p1z = pmax(trim, pmin(p1z, 1 - trim)) # trim ps to be \in [trim, 1 - trim]
  paz = A * p1z + (1-A) * (1-p1z) # probability of P(A=ai|Z)
  
  return(paz)
}



#' IF of E Left Side 
#' Only put inputs for specific A=a
#' (e.g. mu, theta are predictions only for A=a)
#' Will have to adjust paz inputs to be P(A=A_i) (eg vs P(A=1))
#' (IA=a/paZ) bt Omega b 
#' (d thetat + theta dt - 2 theta thetatt)
#' @param paz (vector) probability of treatment given Z
#'            index is in order of mu and theta
#' @output list of matrices of 
#'                   nrow = #samples (that have trtmt A=a)
#'                   ncol = length of alpha
#'          Inefficient for now...
#'          later: probably one matrix of n x len(alpha)^2
#'                and reform
#' @example 
#' testIFEL = getIFEL(B = dl$B0, 
#'                  D = dl$D0,
#'                  Omega = Omega_Identity, 
#'                  # mu    = getInputMu(   Y = dl$Y0, Z = Z0), 
#'                  theta = getInputTheta(D = dl$D0, Z = Z0),
#'                  paz   = rep(.5, times = nrow(dl$B0)))
getIFEL <- function(B, D, Omega, theta, paz) {
  
  # # test w inputs
  # B = dl$B0
  # D = dl$D0
  # Omega = Omega_Identity 
  # mu    = getInputMu(   Y = dl$Y0, Z = Z0) 
  # theta = getInputTheta(D = dl$D0, Z = Z0)
  # paz   = rep(.5, times = nrow(theta)) # fake input
  
  
  n = nrow(theta)
  
  # IFEL = rep(0, times = n) 
  IFEL = list()
  
  # centering term E [ theta bt Omega b thetat]
  # est w/ empirical avg P_n {theta bt Omega b thetat}
  center = matrix(0, nrow = ncol(theta), ncol = ncol(theta))
  for(i in 1:n) {
    center = center + 
             theta[i, ] %*% t(B[i, ]) %*%
             Omega %*% 
             B[i, ] %*% t(theta[i, ])
  }
  center = center / n
  
  # First parts of IF
  # IA=a/paZ (bt Omega b) 
  #  * (d thetat + theta dT - 2 theta thetat)
  #  + theta bt Omega b thetat
  for(i in 1:n) {
    # Indicator{A=a} = 1 always bc of input
    first = (1/paz[i]) * t(B[i, ]) %*% Omega %*% B[i, ] 
    first = first |> as.numeric()
    second =  D[i, ]     %*% t(theta[i, ]) +
              theta[i, ] %*% t(D[i, ])  -
            2*theta[i, ] %*% t(theta[i, ])
    
    # theta  bt Omega b thetat
    third = theta[i, ] %*% t(B[i, ]) %*% 
            Omega %*% 
            B[i, ] %*% t(theta[i, ])
    
    IFEL[[i]] = first * second + third - center
    
  }
  
  return(IFEL)
}




#' IF of E Right Side 
#' Only put inputs for specific A=a
#' (e.g. mu, theta are predictions only for A=a)
#' Will have to adjust paz inputs to be P(A=A_i) (eg vs P(A=1))
#' each IFER_i is of dim J x 1
#' @param paz (vector) probability of treatment given Z
#'            index is in order of mu and theta
#' @output matrix of dim n x J
#'                   nrow = #samples (that have trtmt A=a)
#'                   ncol = length of alpha
#' (so can extract IF for sample i with result[i, ])
#' @example
#' testIFER = getIFER(B = dl$B0, 
#'                  D = dl$D0, 
#'                  Y = dl$Y0,
#'                  Omega = Omega_Identity, 
#'                  mu    = getInputMu(   Y = dl$Y0, Z = Z0),
#'                  theta = getInputTheta(D = dl$D0, Z = Z0),
#'                  paz   = rep(.5, times = nrow(dl$B0)))
#' dim(testIFER)
#' head(testIFER)
getIFER <- function(B, D, Y, Omega, mu, theta, paz) {
  
  # # test w inputs
  # B = dl$B0
  # D = dl$D0
  # Omega = Omega_Identity 
  # mu    = getInputMu(   Y = dl$Y0, Z = Z0) 
  # theta = getInputTheta(D = dl$D0, Z = Z0)
  # paz   = rep(.5, times = nrow(theta)) # fake input
  
  
  n = nrow(theta)
  
  # IFER = rep(0, times = n) 
  IFER = matrix(0, nrow = n, ncol = ncol(theta)) # n x length(alpha)
  
  # centering term E [ theta bt Omega b mut]
  # est w/ empirical avg P_n {theta bt Omega b mut}
  center = matrix(0, nrow = ncol(theta), ncol = 1)
  for(i in 1:n) {
    center = center + 
             theta[i, ] %*% t(B[i, ]) %*%
             Omega %*% 
             B[i, ] * mu[i]
  }
  center = center / n
  # center = t(center) # transpose to be 1 x J
  
  # First parts of IF
  # IA=a/paZ (bt Omega b) 
  #  * (d thetat + theta dT - 2 theta thetat)
  #  + theta bt Omega b thetat
  for(i in 1:n) {
    # Indicator{A=a} = 1 always bc of input
    first = (1/paz[i]) * t(B[i, ]) %*% Omega %*% B[i, ] 
    first = first |> as.numeric()
    second =  D[i, ]     * mu[i] +
              theta[i, ] * Y[i]  -
            2*theta[i, ] * mu[i]
    
    
    # theta  bt Omega b thetat
    third = theta[i, ] %*% t(B[i, ]) %*% 
            Omega %*% 
            B[i, ] %*% mu[i]
    
    IFER[i, ] = first * second + third - center
    
  }
  
  return(IFER)
}





getEL <- function(B, Omega, theta) {
   # alpha_pi = (EL)^{-1} ER
  EL = matrix(0, nrow = ncol(theta), ncol = ncol(theta)) 
  for(i in 1:nrow(B)) {
    # Z[i, ]
    EL = EL + 
         #   t x 1         1 x b        b x b     b x 1       1 x t       
         (theta[i, ]) %*% t(B[i, ]) %*% Omega %*% B[i, ] %*% t(theta[i, ])
         # theta(Z[i, ]) %*% t(B[i, ]) %*% Omega %*% B[i, ] %*% t(theta(Z))
  }
  EL = EL / (nrow(B))
  return(EL)
}





########################################################################
# Main functions
# (that implement methods)
########################################################################

#'
#' @param Y (vector) of length n
#' @param B (matrix) n x dim b transformation
#' @param D (matrix) n x dim dW
#' @param mu (vector) n x 1 of predicted values of Y given Zs,A
#' @param theta (matrix) n x #length(d(W)) = n x #transformed NCEs
#'                       of predicted values of d(W)s given Zs,A
#' @param Omega (matrix) dim b x dim b weight matrix first weight matrix
#'                default is identity
#' @param steps (integer) number of reweighting steps 
PlugIn <- function(Y,B,D, mu, theta, Omega=NULL, steps=1) {
  
  # B = dl$B 
  # Omega = diag(1, ncol = 3, nrow = 3)
  # mu = testMuVector
  # theta = testThetaMatrix
  
  if(is.null(Omega)) {
    Omega = diag(rep(1, ncol(B)))
  }
  
  # Using matrices
  dim_dW = ncol(theta)
  
  for(step in 1:steps) {
    # alpha_pi = (EL)^{-1} ER
    EL = matrix(0, nrow = dim_dW, ncol = dim_dW) 
    for(i in 1:nrow(B)) {
      # Z[i, ]
      EL = EL + 
           #   t x 1         1 x b        b x b     b x 1       1 x t       
           (theta[i, ]) %*% t(B[i, ]) %*% Omega %*% B[i, ] %*% t(theta[i, ])
           # theta(Z[i, ]) %*% t(B[i, ]) %*% Omega %*% B[i, ] %*% t(theta(Z))
    }
    EL = EL / (nrow(B))
    
    ER = matrix(0, nrow = dim_dW, ncol = 1) 
    for(i in 1:nrow(B)) {
      # Z[i, ]
      ER = ER + 
           (theta[i, ]) %*% t(B[i, ]) %*% Omega %*% B[i, ] %*% t(mu[i])
           # theta(Z[i, ]) %*% t(B[i, ]) %*% Omega %*% B[i, ] %*% t(mu(Z))
    }
    ER = ER / (nrow(B))
    
    alpha = solve(EL) %*% ER
    
    # print(Omega)
    # print(alpha)
    
    # update the weight matrix
    Omega = calc_weight(Ya=Y, 
                        Ba=B, 
                        Da=D, 
                        alphaa=alpha)
    
  }
  
  return(alpha)
}

#' One Step Estimator for alpha
#' @param alpha_pi plug in estimate, matrix #dim alpha x 1
#' @param IFEL
#' @param IFER
OneStep <- function(alpha_pi, IFEL, IFER, EL) {
  # alpha_pi = alpha0
  # IFEL = testIFEL0
  # IFER = testIFER0
  # EL = testEL0
  
  inside = rep(0, times = nrow(alpha_pi))
  for(i in 1:nrow(IFER)) {
    inside = inside + 
             (IFER[i, ] - IFEL[[i]] %*% alpha_pi)
               
  }
  inside = inside / nrow(IFER)
  return(solve(EL) %*% inside)
}

########################################################################
# Wrapper functions 
# (for ease of use)
########################################################################

#' wrapper for estimating alpha Plug In
#' @param dataList (list) of named data objects needed
#' @param gmmSteps (integer) number of steps to take of reweighting Omega matrix
#' @param Z0 (matrix) of Z (or transformed Z) values used to estimate nuisance functions
#'           e.g. (E(Y|A=0, Z0)) for treatment A=0
#'           If unspecified (NULL), then use named B0 matrix from dataList
#' @param Z1 (matrix) of Z (or transformed Z) values used to estimate nuisance functions
#'           e.g. (E(Y|A=1, Z1)) for treatment A=1
#'           If unspecified (NULL), then use named B1 matrix from dataList
estCondMomentOCBPI <- function(dataList, gmmSteps, Z0=NULL, Z1=NULL) {
    if(is.null(Z0)) {Z0 = dataList$B0}
    if(is.null(Z1)) {Z1 = dataList$B1}
    # variables used many times
    mu0    = getInputMu(   Y = dataList$Y0, Z = Z0) # hat E(Y|Z,A=0)
    mu1    = getInputMu(   Y = dataList$Y1, Z = Z1) # 
    theta0 = getInputTheta(D = dataList$D0, Z = Z0) # hat E(Wi's|Z,A=0)
    theta1 = getInputTheta(D = dataList$D1, Z = Z1) #
    Omega_Identity = diag(1, nrow=ncol(dataList$B0), ncol=ncol(dataList$B0))
    alpha0 = PlugIn(Y = dataList$Y0,
                    B = dataList$B0,
                    D = dataList$D0,
                    mu    = mu0,
                    theta = theta0,
                    Omega = Omega_Identity,
                    steps = gmmSteps)
    alpha1 = PlugIn(Y = dataList$Y1,
                    B = dataList$B1,
                    D = dataList$D1,
                    mu    = mu1,
                    theta = theta1,
                    Omega = Omega_Identity,
                    steps = gmmSteps)

    res = list(alpha0 = alpha0,
               alpha1 = alpha1,
               ATE    = mean(dataList$D %*% (alpha1 - alpha0)))
    return(res)
}



#' wrapper for estimating alpha OneStep 
#' (and can also return the resulting alpha Plug In)
#' @param dataList (list) of named data objects needed
#' @param gmmSteps (integer) number of steps to take of reweighting Omega matrix
#' @param Z0 (matrix) of Z (or transformed Z) values used to estimate nuisance functions
#'           e.g. (E(Y|A=0, Z0)) for treatment A=0
#'           If unspecified (NULL), then use named B0 matrix from dataList
#' @param Z1 (matrix) of Z (or transformed Z) values used to estimate nuisance functions
#'           e.g. (E(Y|A=1, Z1)) for treatment A=1
#'           If unspecified (NULL), then use named B1 matrix from dataList
#' @param trim (numeric) to trim estimated probs if wanted (0 indicates no trimming)
estCondMomentOCBOS <- function(dataList, gmmSteps, Z0=NULL, Z1=NULL, trim=0) {
    if(is.null(Z0)) {Z0 = dataList$B0}
    if(is.null(Z1)) {Z1 = dataList$B1}
    
    # variables used many times
    mu0    = getInputMu(   Y = dataList$Y0, Z = Z0) # hat E(Y|Z,A=0)
    mu1    = getInputMu(   Y = dataList$Y1, Z = Z1) # 
    theta0 = getInputTheta(D = dataList$D0, Z = Z0) # hat E(Wi's|Z,A=0)
    theta1 = getInputTheta(D = dataList$D1, Z = Z1) #
    paz = getInputPaz(A = c(rep(0, nrow(Z0)),      # hat p(A=ai|Z)
                            rep(1, nrow(Z1))),
                      Z = rbind(Z0, Z1),
                      trim = trim)
    paz0 = paz[  1:nrow(Z0)]
    paz1 = paz[-(1:nrow(Z0))]
    Omega_Identity = diag(1, nrow=ncol(dataList$B0), ncol=ncol(dataList$B0))
    
    # For A=0
    alpha_pi0 = PlugIn(Y = dataList$Y0,
                    B = dataList$B0,
                    D = dataList$D0,
                    mu    = mu0,
                    theta = theta0,
                    Omega = Omega_Identity,
                    steps = gmmSteps)

    IFEL0 = getIFEL(B = dataList$B0, 
                  D = dataList$D0,
                  Omega = Omega_Identity, 
                  theta = theta0,
                  paz   = paz0)

    IFER0 = getIFER(B = dataList$B0, 
                  D = dataList$D0, 
                  Y = dataList$Y0,
                  Omega = Omega_Identity, 
                  mu    = mu0,
                  theta = theta0,
                  paz   = paz0)

    EL0 = getEL(B = dataList$B0, 
                Omega = Omega_Identity, 
                theta = theta0)

    alpha_os0 = OneStep(alpha_pi=alpha_pi0,
                        IFEL = IFEL0,
                        IFER = IFER0,
                        EL = EL0)

    # For A=1
    alpha_pi1 = PlugIn(Y = dataList$Y1,
                    B = dataList$B1,
                    D = dataList$D1,
                    mu    = mu1,
                    theta = theta1,
                    Omega = Omega_Identity,
                    steps = gmmSteps)
    
    IFEL1 = getIFEL(B = dataList$B1, 
                  D = dataList$D1,
                  Omega = Omega_Identity, 
                  theta = theta1,
                  paz   = paz1)

    IFER1 = getIFER(B = dataList$B1, 
                  D = dataList$D1, 
                  Y = dataList$Y1,
                  Omega = Omega_Identity, 
                  mu    = mu1,
                  theta = theta1,
                  paz   = paz1)

    EL1 = getEL(B = dataList$B1, 
                Omega = Omega_Identity, 
                theta = theta1)

    alpha_os1 = OneStep(alpha_pi=alpha_pi1,
                        IFEL = IFEL1,
                        IFER = IFER1,
                        EL = EL1)
    
    # Return named list of desired results
    res = list(alpha_pi0 = alpha_pi0,
               alpha_pi1 = alpha_pi1,
                 ATE_pi  = mean(dataList$D %*% (alpha_pi1 - alpha_pi0)),
               alpha_os0 = alpha_os0,
               alpha_os1 = alpha_os1,
                 ATE_os  = mean(dataList$D %*% (alpha_os1 - alpha_os0)))
    return(res)
}


