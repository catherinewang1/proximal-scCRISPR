








#' Makes a function that performs Proximal Inference with the Outcome Confounding Bridge using 2SLS 
#' designed to be modular 
#' with the specified models/algorithms/procedures for each stage
#' modelW and modelY can be defined similarly, but with different parameters to save space/time
#' they can be defined with the same custom functions to more easily set up
#' @param modelW (function) should take in Y (response), X (design matrix), A (trtmt) and return a 
#'        list of:
#'             EY1, EY0 = mean potential outcomes 
#'             predA, predA1, predA0 = pred vals at obs A, setting A=1, setting A=0, respectively
#'        
#'        only requires predicted values at A=A (at observed trtmt A)
#' @param modelY (function) should take in Y (response), X (design matrix), A (trtmt), and
#'        returnIntermediateATEs (boolean, or a list of integers)
#'        and return a list of:
#'             EY1, EY0 = mean potential outcomes 
#'             predA, predA1, predA0 = pred vals at obs A, setting A=1, setting A=0, respectively
#' modelY should also allow parameter returnIntermediateATEs (boolean, or a list of integers), 
#' which should additionally return a data.frame of intermediate EY1, EY0, ATE values 
#' based on the specified number of NCOs if given (returnIntermediateATEs as vector) or TRUE 
#' indicates all 1, 2, ..., maxNCOs
#' (e.g. returnIntermediateATEs=c(1, 4, 6) would return 
#' list(EY1, EY0, intermediateATEs=data.frame(numNCOs=c(1, 4, 6), EY1=..., EY0=..., ATE=...))
#' @returns function that performs 2SLS with the specified modelW and modelY algs
makeOCB2SLS <- function(modelW, modelY) {

    #' function to return that actually performs 2SLS estimation using the 
    #' initially specificied modelW and modelY (the algorithms used for 2SLS)
    #' @param Y (vector) of response
    #' @param A (vector) of treatment assignment
    #' @param W (matrix/dataframe) of Negative Control Outcomes
    #' @param Z (matrix/dataframe) of Negative Control Exposures
    #' @param returnIntermediateATEs (boolean or vector of integers) 
    #'     should it additionally return a data.frame of intermediate EY1, EY0, ATE values 
    #'     based on the specified number of NCOs if given (returnIntermediateATEs as vector), 
    #'     or TRUE indicates all 1, 2, ..., maxNCOs
    #'     (e.g. returnIntermediateATEs=c(1, 4, 6) would return 
    #'     list(EY1, EY0, intermediateATEs=dataframe(numNCOs=c(1, 4, 6), EY1=..., EY0=..., ATE=...))
    #'          FALSE=using max NCO
    #'          TRUE =using 1, 2, ..., max NCO
    #'          vector=using this specified number
    #' @returns list(EY1, EY0, ATE, [intermediateATEs <data.frame>])
    innerOCB2SLS <- function(Y, A, W, Z, returnIntermediateATEs=FALSE) {
        # ==========================================================================================
        # === 1st Stage: estimate E(Wj | Z, A) using modelW algorithm
        # ==========================================================================================
        ZA = cbind(Z, A)
        W_ZA = matrix(NA, nrow = nrow(W), ncol = ncol(W))
        for(j in 1:ncol(W)) {
            W_ZA[, j] = modelW(Y = W[,j], A = A, X = ZA)$predA
        }
        
        # ==========================================================================================
        # === 2nd Stage: estimate Y ~ A + E(Wj | Z, A)
        # ==========================================================================================
        # secondStageFit = modelY(Y = Y, X = W_ZA, A = A, returnIntermediateATEs=returnIntermediateATEs)
        # return(secondStageFit)

        # --------------------- max NCOs -----------------------------------------------------------
        # return one result using max NCOs 
        if(is.logical(returnIntermediateATEs) & !returnIntermediateATEs) { # returnIntermediateATEs=FALSE
            secondStageRes = modelY(Y = Y, A = A, X = W_ZA, returnPreds=FALSE)            
            return(c(secondStageRes,
                     list(ATE = secondStageRes$EY1 - secondStageRes$EY0)))
        } 
        
        # --------------------- intermediate NCOs --------------------------------------------------
        # set idx to evaluate intermediate ATEs at
        if(is.logical(returnIntermediateATEs)) { # returnIntermediateATEs=TRUE
            # all 1:ALL NCOs (should not be used often, and only for small maxNCO)
            NCOidx = 1:ncol(W)
        } else { # returnIntermediateATEs = vector of ints (hopefully)
            NCOidx = returnIntermediateATEs
        }

        # return results using NCOidx  
        EY1s = rep(NA, length(NCOidx))
        EY0s = rep(NA, length(NCOidx))
        for(j in 1:length(NCOidx)) {
            secondStageResj = modelY(Y = Y, A = A, X = W_ZA[, 1:NCOidx[j]], returnPreds=FALSE) 
            EY1s[j] = secondStageResj$EY1
            EY0s[j] = secondStageResj$EY0
        }

        intermediateATEs = data.frame(idx    = 1:length(NCOidx),
                                      numNCO = NCOidx,
                                      EY1    = EY1s,
                                      EY0    = EY0s,
                                      ATE    = EY1s - EY0s)
        
        res = list(EY1=intermediateATEs[length(NCOidx), 'EY1'], # maxNCO
                   EY0=intermediateATEs[length(NCOidx), 'EY0'],
                   intermediateATEs=intermediateATEs)

        return(res)
    }
    
    return(innerOCB2SLS)
}



#' Estimate using Linear regression for various stages of 2SLS 
#' 
#' e.g. a modelW or modelY
#' set modelW = modelLin
#' set modelY = modelLin
#' @param Y (vector) of response
#' @param A (vector) of treatment assignment
#' @param X (matrix/dataframe) of predicted Negative Control Outcomes E(Wj | Z, A) 
#' @param returnPreds (boolean) return predicted values at each sample point OR just return EY1 and EY0's (saves space)
#' @return list(EY1, EY0, [predA, predA1, predA0])
modelLin <- function(Y, A, X, returnPreds=TRUE) {
    
    AX  = cbind(A, X)
    fit = glm.fit(Y = Y, X = AX, family = gaussian)

    predA   = predict.glm(fit, type = 'response')
    AX[, 1] = 1 # set trmt = 1
    predA1  = predict.glm(fit, newdata = AX, type = 'response')
    AX[, 1] = 0 # set trmt = 0
    predA0  = predict.glm(fit, newdata = AX, type = 'response')

    res = list(EY1 = mean(predA1), 
               EY0 = mean(predA0))
    if(returnPreds) {
        res = c(res, 
                list(predA =predA,
                     predA1=predA1,
                     predA0=predA0))
    }
    return(res)

}

#' Make function that estimates using Linear Reg w Regularization for various stages of 2SLS 
#' 
#' e.g. a modelW or modelY
#' modelLinLasso = makeModelLinReg(alpha=1)
#' modelLinRidge = makeModelLinReg(alpha=0)
#' modelLinElNet = makeModelLinReg(alpha=.3)
#' set modelW = modelLinLasso
#' set modelY = modelLinLasso
#' @param Y (vector) of response
#' @param A (vector) of treatment assignment
#' @param X (matrix/dataframe) of predicted Negative Control Outcomes E(Wj | Z, A) 
#' @param returnPreds (boolean) return predicted values at each sample point OR just return EY1 and EY0's (saves space)
#' @return list(EY1, EY0, [predA, predA1, predA0])
#' @importFrom glmnet cv.glmnet glmnet
makeModelLinRegularize <- function(alpha, nfolds=10) {
    modelLinRegularize <- function(Y, A, X, returnPreds=TRUE) {
        AX      = cbind(A, X)
        cv_fit  = glmnet::cv.glmnet(y = Y, x = AX, alpha = alpha, nfolds=nfolds) # plot(cv_fit)
        fit     = glmnet::glmnet(y = Y, x = AX, alpha = alpha, lambda = cv_fit$lambda.1se) # 1se rule 
        
        predA   = predict(fit, newx = AX)
        AX[, 1] = 1 # set trmt = 1
        predA1  = predict(fit, newdata = AX)
        AX[, 1] = 0 # set trmt = 0
        predA0  = predict(fit, newdata = AX)
    }
    
    return(modelLinRegularize)
}


#' Estimate using Poisson GLM regression for various stages of 2SLS 
#' 
#' e.g. a modelW or modelY
#' set modelW = modelPois
#' set modelY = modelPois
#' @param Y (vector) of response
#' @param A (vector) of treatment assignment
#' @param X (matrix/dataframe) of predicted Negative Control Outcomes E(Wj | Z, A) 
#' @param returnPreds (boolean) return predicted values at each sample point OR just return EY1 and EY0's (saves space)
#' @return list(EY1, EY0, [predA, predA1, predA0])
modelPois <- function(Y, A, X, returnPreds=TRUE) {
    # Annoying: predict.glm only works for return from glm and not glm.fit
    # prefer glm.fit bc of y, x input format with no formula 
    # (and no creating new dataframe that uses space)
    # can calculate the predicted values ourselves, but more risky
    # AX  = cbind(A, X) |> as.data.frame()
    # colnames(AX) = c('A', paste0('X', 1:(ncol(AX)-1)))
    # rm(A);rm(X)# ; gc(verbose = FALSE) # remove for space
    # fit = glm.fit(y = Y, x = AX, family = poisson())
    # predA   = predict.glm(fit, newdata = AX, type = 'response') # would have to
    # AX[, 1] = 1 # set trmt = 1                                  # change to a
    # predA1  = predict.glm(fit, newdata = AX, type = 'response') # custom predict
    # AX[, 1] = 0 # set trmt = 0                                  # function
    # predA0  = predict.glm(fit, newdata = AX, type = 'response') # https://stat.ethz.ch/pipermail/r-help/2004-September/058282.html


    # using the glm function w formula, then predict
    X = log(X)
    YAX = cbind(Y, A, X) |> as.data.frame()
    colnames(YAX) = c('Y', 'A', paste0('X', 1:(ncol(YAX)-2)))
    rm(Y); rm(A); rm(X)
    fit = glm('Y ~ .', data =  YAX, family = poisson()) # Y on all A + X's

    predA   = predict.glm(fit, newdata = YAX, type = 'response')
    YAX$A = 1 # set trmt = 1
    predA1  = predict.glm(fit, newdata = YAX, type = 'response')
    YAX$A = 0 # set trmt = 0
    predA0  = predict.glm(fit, newdata = YAX, type = 'response')

    res = list(EY1 = mean(predA1), 
               EY0 = mean(predA0))
    if(returnPreds) {
    res = c(res, 
            list(predA =predA,
                 predA1=predA1,
                 predA0=predA0))
    }
    return(res)

}


#' For the second stage, we want to log-transform the predicted E(W|Z,A)
makeModelPois <- function(logTransX=FALSE) {
    #' Estimate using Poisson GLM regression for various stages of 2SLS 
    #' 
    #' e.g. a modelW or modelY
    #' set modelW = modelPois
    #' set modelY = modelPois
    #' @param Y (vector) of response
    #' @param A (vector) of treatment assignment
    #' @param X (matrix/dataframe) of predicted Negative Control Outcomes E(Wj | Z, A) 
    #' @param returnPreds (boolean) return predicted values at each sample point OR just return EY1 and EY0's (saves space)
    #' @return list(EY1, EY0, [predA, predA1, predA0])
    modelPois <- function(Y, A, X, returnPreds=TRUE) {
        # Annoying: predict.glm only works for return from glm and not glm.fit
        # prefer glm.fit bc of y, x input format with no formula 
        # (and no creating new dataframe that uses space)
        # can calculate the predicted values ourselves, but more risky
        # AX  = cbind(A, X) |> as.data.frame()
        # colnames(AX) = c('A', paste0('X', 1:(ncol(AX)-1)))
        # rm(A);rm(X)# ; gc(verbose = FALSE) # remove for space
        # fit = glm.fit(y = Y, x = AX, family = poisson())
        # predA   = predict.glm(fit, newdata = AX, type = 'response') # would have to
        # AX[, 1] = 1 # set trmt = 1                                  # change to a
        # predA1  = predict.glm(fit, newdata = AX, type = 'response') # custom predict
        # AX[, 1] = 0 # set trmt = 0                                  # function
        # predA0  = predict.glm(fit, newdata = AX, type = 'response') # https://stat.ethz.ch/pipermail/r-help/2004-September/058282.html


        # using the glm function w formula, then predict
        if(logTransX) {
            X = log(X)
        }         
        YAX = cbind(Y, A, X) |> as.data.frame()
        colnames(YAX) = c('Y', 'A', paste0('X', 1:(ncol(YAX)-2)))
        rm(Y); rm(A); rm(X)
        fit = glm('Y ~ .', data =  YAX, family = poisson()) # Y on all A + X's

        predA   = predict.glm(fit, newdata = YAX, type = 'response')
        YAX$A = 1 # set trmt = 1
        predA1  = predict.glm(fit, newdata = YAX, type = 'response')
        YAX$A = 0 # set trmt = 0
        predA0  = predict.glm(fit, newdata = YAX, type = 'response')

        res = list(EY1 = mean(predA1), 
                   EY0 = mean(predA0))
        if(returnPreds) {
        res = c(res, 
                list(predA =predA,
                     predA1=predA1,
                     predA0=predA0))
        }
        return(res)

    }
    modelPois
}

#' Estimate using Negative Binomial GLM regression for various stages of 2SLS 
#' 
#' e.g. a modelW or modelY
#' set modelW = modelPois
#' set modelY = modelPois
#' @param Y (vector) of response
#' @param A (vector) of treatment assignment
#' @param X (matrix/dataframe) of predicted Negative Control Outcomes E(Wj | Z, A) 
#' @param returnPreds (boolean) return predicted values at each sample point OR just return EY1 and EY0's (saves space)
#' @importFrom MASS glm.nb
modelNegBin <- function(Y, A, X, returnPreds=TRUE) {
    
    AYX  = cbind(A, Y, X)
    colnames(AYX) = c('A', 'Y', paste0('X', 1:ncol(X)))
    fit = MASS::glm.nb(formula = sprintf('Y ~ %s',  paste0(colnames(AYX)[-2], collapse = ' + ')),
                       data    = AYX)
    
    predA   = predict.glm(fit, type = 'response')
    AYX[, 1] = 1 # set trmt = 1
    predA1  = predict.glm(fit, newdata = AYX, type = 'response')
    AYX[, 1] = 0 # set trmt = 0
    predA0  = predict.glm(fit, newdata = AYX, type = 'response')

    res = list(EY1 = mean(predA1), 
               EY0 = mean(predA0))
    if(returnPreds) {
        res = c(res, 
                list(predA =predA,
                     predA1=predA1,
                     predA0=predA0))
    }
    return(res)

}


# #' Reconfigured: More code put into make2SLS... so it is easier to define these smaller modelXXX fns
# #' Estimate using Poisson GLM regression for various stages of 2SLS 
# #' 
# #' e.g. a modelW or modelY
# #' set modelW = modelPois(Y, A, X,  returnPreds=TRUE, returnIntermediateATEs=FALSE)
# #' set modelY = modelPois(Y, A, X, returnPreds=FALSE, returnIntermediateATEs=c(1, 3, 5, 10)
# #' @param returnPreds (boolean) return predicted values at each sample point OR just return EY1 and EY0's (saves space)
# #' @param returnIntermediateATEs (boolean or vector of ints) 
# #'      FALSE=using max NCO
# #'      TRUE =using 1, 2, ..., max NCO
# #'      vector=using this specified number
# #'      specifying TRUE or vector will result in returnPreds=FALSE! only give mean potential outcomes (saves space)
# #' @return list()
# modelPois <- function(Y, A, X, returnPreds=TRUE, returnIntermediateATEs=FALSE) {
    
#     # ================================ max NCOs ================================
#     # return one result using max NCOs 
#     if(is.logical(returnIntermediateATEs) & !returnIntermediateATEs) { # returnIntermediateATEs=FALSE
#         AX  = cbind(A, X)
#         fit = glm(Y = Y, X = AX, family = 'poisson')

#         predA = predict(fit, type = 'response')
#         AX[, 1] = 1 # set trmt = 1
#         predA1 = predict(fit, newdata = AX, type = 'response')
#         AX[, 1] = 0 # set trmt = 0
#         predA0 = predict(fit, newdata = AX, type = 'response')

#         res = list(EY1 = mean(predA1), 
#                    EY0 = mean(predA0))
#         if(returnPreds) {
#            res = c(res, 
#                    list(predA=predA,
#                         predA1=predA1,
#                         predA0=predA0))
#         }
#         return(res)
#     } 
    
#     # ================================ intermediate NCOs ========================
#     # set idx to evaluate intermediate ATEs at
#     if(is.logical(returnIntermediateATEs)) { # returnIntermediateATEs=TRUE
#         # all 1:ALL NCOs (should not be used often, and only for small maxNCO)
#         NCOidx = 1:ncol(X)
#     } else { # returnIntermediateATEs = vector of ints (hopefully)
#         NCOidx = returnIntermediateATEs
#     }

#     # return results using NCOidx  
#     EY1s = rep(NA, length(NCOidx))
#     EY0s = rep(NA, length(NCOidx))
#     for(j in NCOidx) {
#         AX = cbind(A, X[, 1:j])
#         fit = glm(Y = Y, X = AX, family = 'poisson')

#         predA = predict(fit, type = 'response')
#         AX[, 1] = 1 # set trmt = 1
#         predA1 = predict(fit, newdata = AX, type = 'response')
#         AX[, 1] = 0 # set trmt = 0
#         predA0 = predict(fit, newdata = AX, type = 'response')

#         EY1s[j] = mean(predA1)
#         EY0s[j] = mean(predA0)
#     }
#     intermediateATEs = data.frame(idx    = 1:length(NCOidx),
#                                   numNCO = NCOidx,
#                                   EY1    = EY1s,
#                                   EY0    = EY0s,
#                                   ATE    = EY1s - EY0s)
    
#     res = list(EY1=intermediateATEs[length(NCOidx), 'EY1'], # maxNCO
#                EY0=intermediateATEs[length(NCOidx), 'EY0'],
#                intermediateATEs=intermediateATEs)

#     return(res)

# }


#' wrapper for pci2s using specific number'
#' @param nco_type (character) 'poisson' or 'negbin', or a vector of these
#' @param num_NCs
runpci2s <- function(Y, A, X, W, Z, nco_type, num_NCs=NULL) {

    # nco type input as single char or vec of chars
    if(length(nco_type) == 1) {
        nco_types = rep(nco_type, ncol(W))
    } else {
        nco_type = nco_type
    }
    
    if(is.null(num_NCs)) {
        # use all
        num_NCs = ncol(W)
    } 

    paramEstimates = c() # of the parameter coef of A
    se             = c() # of the parameter coef of A
    pval           = c() # of the parameter coef of A
    for(num_NC in num_NCs) {
    p2sls_result <- pci2s::p2sls.negbin(Y = Y, 
                                        A = A, 
                                        X = X,
                                        W = W[, 1:num_NC], 
                                        # Z = Z[, 1:num_NC],
                                        Z = Z,
                                        nco_type = nco_types[1:num_NC])
    paramEstimates = c(paramEstimates, p2sls_result$summary_second_stage['A', 'Estimate'])
    se             = c(se,             p2sls_result$summary_second_stage['A', 'Std. Error'])
    pval           = c(pval,           p2sls_result$summary_second_stage['A', 'Pr(>|z|)'])
    }

    res = data.frame(idx      = 1:length(num_NCs),
                   numNCO   = num_NC,
                   Estimate = paramEstimates,
                   se       = se,
                   pval     = pval)
    
    return(res)    
}


