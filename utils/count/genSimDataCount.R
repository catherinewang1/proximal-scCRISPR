# Simulate Data for:
# Exposure A
# Response Y (COUNTS)
# Unmeasured Confounding U
# Negative Control Exposure Z (COUNTS)
# Negative Control Outcome W (COUNTS)

#         beta
#      A ----> Y      
#      ^\     /^
#gammaX  \   / gammaY
#          U
#gammaZ   /  \ gammaW
#       \/   \/
#       Z     W



# helper function: inverse logit to transform R to [0,1]
#        1
#  ------------
#   1 + exp(-x)
invlogit <- function(x) {1/(1+exp(-x))}



#' Simulate Counts from Poisson
#' @param N (integer) number of data points
#' @param U_dim (vector) vector of unconfounders' mean
#' @param U_sd  (vector) vector of unonfounders' sd
#' @param beta  (numeric) linear effect of A on Y
#' @param gammaA (vector) of numeric of length=dim(U) linear effect (logistic) of U on A 
#' @param gammaY (vector) of numeric of length=dim(U) linear effect            of U on Y
#' @param gammaZ (vector) of numeric of length=dim(U) linear effect            of U on Z
#' @param gammaW (vector) of numeric of length=dim(U) linear effect            of U on W
#' @param returnAsList (boolean) return result as a List (ith named Y,A,...) or as one dataframe
simDataCount <- function(N, U_mean, U_sd, beta, gammaA, gammaY, gammaZ, gammaW,
                         muA=0, returnAsList = FALSE) {# , dispersionY=1, dispersionZ=1, dispersionW=1) {
    
    # debugging
    # N = 20
    # U_mean = c(0, -1)
    # U_sd = c(.5, 1.25)
    # beta = 1
    # gammaA = c(1, -.5)
    # gammaY = c(.2, .8)
    # gammaZ = matrix(c(1, .2,
    #                   .1, 1), nrow=2, byrow = TRUE)
    # gammaW = matrix(c(.3, .9,
    #                   .4, -1), nrow=2, byrow = TRUE)
    # muA = 0
    
	  # check inputs
    assertthat::assert_that(length(U_mean) == length(U_sd))
    assertthat::assert_that(length(gammaA) == length(U_mean))
    assertthat::assert_that(length(gammaY) == length(U_mean))
    assertthat::assert_that(  nrow(gammaZ) == length(U_mean))
    assertthat::assert_that(  nrow(gammaW) == length(U_mean))
    assertthat::assert_that(length(muA)    == 1)


    
    # draw U: matrix of N x dim(U)
    U = mapply(rnorm, n = N, mean = U_mean, sd = U_sd)
    colnames(U) = paste0('U', 1:ncol(U))
    # colMeans(U); sd(U[,1]); sd(U[,2])
    
    # draw A {0,1}: vector of length N
    Agamma = muA + U %*% gammaA
    A = apply(invlogit(Agamma), c(1, 2), function(lambda) { rbinom(n = 1, size = 1, prob=lambda) })
    # plot(invlogit(Agamma), A)
    
    
    # Drawing Poisson's
    # \mu = E(Y|X) = g^-1(X \beta)?
    # X \beta = g(\mu)
    # Poisson: g = log --> X \beta = log(\mu) = log(E(Y|X)) 
    # so for Y|X ~ Pois(\lambda), the mean is \lambda = \mu = E(Y|X) = exp(X \beta)
    
    # Draw Y (Poisson) 
    Ygamma = U %*% gammaY + A * beta
    # Ygamma = U[,1] * gammaY[1] + U[, 2] * gammaY[2] + A * beta
    # Ygamma = A * beta
    Y = apply(exp(Ygamma), c(1, 2), function(lambda) {rpois(n = 1, lambda = lambda)})
    # hist(exp(Ygamma)); hist(Y)
    
    # Draw Z (Poisson)
    Zgamma = U %*% gammaZ 
    Z = apply(Zgamma, c(1, 2), function(x) {rpois(n = 1, lambda = exp(x))})
    colnames(Z) = paste0('Z', 1:ncol(Z))
    
    # Draw W (Poisson)
    Wgamma = U %*% gammaW
    W = apply(Wgamma, c(1, 2), function(x) {rpois(n = 1, lambda = exp(x))})
    colnames(W) = paste0('W', 1:ncol(W))

    # put together
    if(returnAsList) {
        return(list(Y=Y,A=A, Z=Z, W=W, U=U))
    } else {
        return(data.frame(Y, A, Z, W, U))  
    }
    
    
}
# 
# 
# 
# # potentially for NB distn
# dispersionZ = 2
# dispersionW = 2
# 
# 
# 
# 
# 
# 
# 
# 
# head(dat)
# 
# 
# 
# 
# fit = glm('Y ~ A', data = dat, family = poisson())
# fit = glm('Y ~ A', data = dat, family = poisson(link = "log"))
# 
# 
# 
# fit
# summary(fit)
# 
# fit = glm('Y ~ A + U1 + U2', data = dat, family = poisson())
# summary(fit)
# 
# 
# # Negative Binomial? Additional Dispersion parameters for Z and W
# 
# 
# rnbinom(n = N, size = dispersionZ, mu = exp((U %*% gammaZ)[, 1]))
# rnbinom(n = N, size = dispersionZ, mu = 4)
# 
# 
# 
# 
# #sample size
# n <- 10
# #regression coefficients
# beta0 <- 1
# beta1 <- 0.2
# #generate covariate values
# x <- runif(n=n, min=0, max=1.5)
# #compute mu's
# mu <- exp(beta0 + beta1 * x)
# #generate Y-values
# y <- rpois(n=n, lambda=mu)
# #data set
# data <- data.frame(y=y, x=x)
# 
# fit = glm("y ~ x", data, family = poisson())
# 
# summary(fit)
# 
