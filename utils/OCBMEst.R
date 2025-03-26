

# Starting over with solving for the M-Estimator from scratch... 9/25/24





#' Estimate the ATE using the Outcome Confounding Bridge (OCB) function
#' which is estimated by matching the conditional moments
#' @param dataList (list) of values containing data 
#' Need matrices/vectors/numerics (named)
#' Y1, Y0: response
#' B1, B0: conditional moment basis functions evaluated at Zs
#' D1, D0: function class basis functions evaluated at Ws
#' p1, p0: proportion of samples in each treatment group
#' # @param b_degree (integer) degree of opponent bases 
#' # @param h_degree (integer) degree of h (confounding bridge)
# ' @param Weights (matrix) #b x #b (#b = ncol(dataList$B1)) of weights 
#'                         (equal weights by default)
#' @param returnPreds (boolean) return predicted values or not
#' @param returnAlphas (boolean) return estimated alphas or not
#' estCondMomentPolyOCB(df, b_degree=1, h_degree=1)
mest <- function(dataList, Weights=NULL, returnPreds=FALSE, returnAlphas=FALSE, returnVar=FALSE) {
    # E(b(Z) Y) = E(b(Z) d(W)^T) alpha
    # --> calculate alpha from \hat alpha = P_n(b(Z) d(W)^T)^-1  P_n(b(Z) Y)
    
    alpha1 = mestCalcAlpha(B=dataList$B1, D=dataList$D1, Y=dataList$Y1, Weights=Weights)
    alpha0 = mestCalcAlpha(B=dataList$B0, D=dataList$D0, Y=dataList$Y0, Weights=Weights)

    # NO WEIGHTS YET
    
    hat_h1 = rbind(dataList$D1, dataList$D0) %*% alpha1 
    hat_h0 = rbind(dataList$D1, dataList$D0) %*% alpha0 

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
    if(returnVar) {
        

        Mhat_inv1 = getMhatinv(B = dataList$B1, D = dataList$D1)
        Mhat_inv0 = getMhatinv(B = dataList$B0, D = dataList$D0)

        f_hatATE = getFhatATE(dl=dataList, 
             alpha1=alpha1, 
             alpha0=alpha0, 
             Mhat_inv1=Mhat_inv1, 
             Mhat_inv0=Mhat_inv0)

        res = c(res, 
                list(VarEY1=mestCalcVarA(dataList=dataList, alpha=alpha1, a=1),
                     VarEY0=mestCalcVarA(dataList=dataList, alpha=alpha0, a=0),
                     VarATE=var(f_hatATE) |> as.numeric()))
    }
    return(res)
}

#' Calculate the alpha_star (for A=1 and A=0) for the 
#' Outcome Confounding Bridge function h
#' with the given subset of Z,W,Y and 
#' with the specified degrees for B (b_degree) and D (h_degree)
#' @param B (matrix) n x # basis fns Z evaluated at basis fns
#' @param D (matrix) n x # basis fns W evaluated at basis fns
#' @param Y (vector) observed values for Outcome
#' @param Weights (matrix) b_degree x b_degree of weights
#' Solution is: 
#'    \hat alpha = P_n(b(Z) d(W)^T)^-1  P_n(b(Z) Y)
mestCalcAlpha <- function(B, D, Y, Weights=NULL) {
    # NO WEIGHTS YET

    # slow but should work
    sumbd = matrix(0,  nrow = ncol(B), ncol = ncol(D)) # maybe matrix mult t(B) %*% D 
    for(i in 1:nrow(B)) {

       sumbd = sumbd + as.matrix(B[i, ], ncol = nrow(B)) %*% D[i, ] 
    }
    invbd = solve((1/nrow(B)) * sumbd)


    sumby = matrix(0, nrow = ncol(B), ncol = 1) # maybe matrix mult diag(Y) %*% B |> colSum()
    for(i in 1:nrow(B)) {
       sumby = sumby + as.matrix(B[i, ], ncol = nrow(B)) * Y[i]
    }    

    alphahat = invbd %*% ((1/nrow(B)) * sumby)
    
    return(alphahat)
}



# For A=a ==========
mestCalcVarA <- function(dataList, alpha, a, Omega=NULL) {
    # TODO: rename Omega = Weights
    Fhat = mestCalcFA(dataList, alpha, a, Omega=Omega)
    return(var(Fhat)) # sample variance of f(O; alpha)
}

mestCalcFA <- function(dataList, alpha, a, Omega=NULL) {
    # h(W;\alpha) + (deriv h) (deriv E(g))^{-1} g
    # d(W)^T alpha + d(W)^T {E(\sqrt{Omega} b(Z) d(W)^T)}^{-1} b(Z)(Y - d(W)^T alpha)

    # setup matrices of all data
    D = rbind(dataList$D1, dataList$D0)
    B = rbind(dataList$B1, dataList$B0) 
    A = c(rep(1, nrow(dataList$D1)), 
          rep(0, nrow(dataList$D0)))
    Y = c(dataList$Y1, dataList$Y0)

    if(is.null(Omega)) {
        Omega_sqrt = diag(rep(1, ncol(B))) # identity 
    } else {
        Omega_sqrt = (Omega_sqrt)**(1/(1+1)) # solve? 
    }
    
    # D %*% alpha 
    
    # rbind(B1, B0) %*% alpha 

    

    # # {E(\sqrt{Omega} b(Z) d(W)^T)}
    # M = Omega_sqrt %*% t(B) %*% D 
    
    # solve(M)
    
    # D %*% solve(M) t(B) 

    # (Y - D %*% alpha)


    # # {E(\sqrt{Omega} b(Z) d(W)^T)}
    # M = Omega_sqrt %*% t(dataList$B) %*% D 
    # rbind(dataList$D1, dataList$D0) %*% solve(M) %*% (dataList$B1 * (dataList$Y1 - dataList$D1 %*% alpha))
    
    # h(W;\alpha) + (deriv h) (deriv E(g))^{-1} g
    

    
    
    # M = {E(\sqrt{Omega} b(Z) d(W)^T)}
    Mhat = matrix(0, nrow = ncol(B), ncol = ncol(D))
    # slo but trying to be correct
    for(i in 1:nrow(D)) { # from 1 to n
        Mhat = Mhat + as.matrix(B[i, ], ncol = ncol(B)) %*% D[i, ] * as.integer(A[i] == a)
    }
    Mhat = -(1/nrow(D)) * Omega_sqrt %*% Mhat
    
    Mhat_inv = solve(Mhat)
    
    # h(W;\alpha) 
    H = as.vector(D %*% alpha)

    # (deriv h) = deriv d()^T alpha = d()
    
    # G = (deriv h) (deriv E(g))^{-1} g
    # g = sqrt(Omega) b(Z) (Y - d()^T alpha) I{A=a}
    G = rep(NA, times = nrow(D)) 
    for(i in 1:nrow(D)) {
      G[i] = D[i, ] %*% 
        Mhat_inv %*% 
        (Omega_sqrt %*%  matrix(B[i, ], nrow = ncol(B), ncol = 1) %*%
           (Y[i] - D[i, ] %*% alpha) *
           as.integer(A[i] == a)) 
    }
    
    Fhat = H + G
    
    return(Fhat)

    
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




# For ATE ==========
#' Calculate the variance of the estimator
mestCalcVar <- function() {

    f = mestCalcF(...)
    var(f)
}


getFhatATE<- function(dl, alpha1, alpha0, Mhat_inv1, Mhat_inv0) {
  # H(W;\alpha) + [deriv H] [deriv G]^{-1} g
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



mestCalcH <- function() {

}

mestCalcDerivH <- function() {

}

mestCalcG <- function() {
    
}

mestCalcDerivG <- function() {

}


