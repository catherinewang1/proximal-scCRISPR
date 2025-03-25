
# Reproducible code

# removing previous install
remove.packages('pci2s')

# reinstalling
# devtools::install_github('https://github.com/KenLi93/pci2s') # or
remotes::install_github('https://github.com/KenLi93/pci2s')

library('pci2s')

# ------------ Example Code from function `pci.lm` doc
N <- 2000
expit  <-  function(x) exp(x)/(1 + exp(x))
U <- runif(N); X <- runif(N)
A <- rbinom(N, 1, expit(-3 + 5 * U + 1 * X))
Y <- rnorm(N, 0.2 + 1.5 * U + 0.2 * X + 0.2 * A, 0.5)
Z <- rnorm(N, 2 * U + 0.5 * X)
W2 <- rexp(N, 0.1 + 1 * U + 0.1 * X)
D2 <- as.numeric(W2 < 5)
W2[D2 == 0] <- 5
W <- cbind(W1 = rnbinom(N, size = 25,
                        mu = exp(2.5 * U + 0.2 * X)),
           W2)

# Errors with the message: 
# Error in colnames(W1) : object 'W1' not found
pci_result <- pci2s::pci.lm(Y = Y, A = A, X = X,
                     W = W, Z = Z,
                     nco_type = c("negbin", "ah"),
                     nco_args = list(list(offset = rep(0, N)),
                                     list(offset = rep(0, N),
                                          event = D2)))


# ------------ Example Code from function `pci.ah` doc works!
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
# Errors with the message: 
# Error in colnames(W1) : object 'W1' not found
pci_result <- pci.loglin(Y = Y, A = A, X = X,
                         W = W, Z = Z,
                         nco_type = c("negbin", "ah"),
                         nco_args = list(list(offset = rep(0, N)),
                                         list(offset = rep(0, N),
                                              event = D2)))

# ------------ Example Code from function `pci.negbin` doc works!
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
pci_result <- pci.negbin(Y = Y, A = A, X = X,
                         W = W, Z = Z,
                         nco_type = c("negbin", "ah"),
                         nco_args = list(list(offset = rep(0, N)),
                                         list(offset = rep(0, N),
                                              event = D2)))
pci_result$summary_first_stage
pci_result$summary_second_stage

# ------------ Example Code from function `pci.ah` doc works!
N <- 2000
U <- runif(N); X <- runif(N)
expit  <-  function(x) exp(x)/(1 + exp(x))
A <- rbinom(N, 1, expit(-3 + 5 * U + 1 * X))
Y <- rexp(N, 0.2 + 1.5 * U + 0.2 * X + 0.2 * A)
D <- as.numeric(Y < 5)
Y[D == 0] <- 5
Z <- rnorm(N, 2 * U + 0.5 * X)
W2 <- rexp(N, 0.1 + 1 * U + 0.1 * X)
D2 <- as.numeric(W2 < 5)
W2[D2 == 0] <- 5
W <- cbind(W1 = rnbinom(N, size = 25, mu = exp(2.5 * U + 0.2 * X)),
           W2)
pci_result <- pci.ah(Y = Y, D = D, A = A, X = X,
                     W = W, Z = Z, variance = TRUE,
                     nco_type = c("negbin", "ah"),
                     nco_args = list(list(offset = rep(0, N)),
                                     list(offset = rep(0, N),
                                          event = D2)))

pci_result$summary_first_stage
pci_result$summary_second_stage




# --------- update?
N <- 2000
expit  <-  function(x) exp(x)/(1 + exp(x))
U <- runif(N); X <- runif(N)
A <- rbinom(N, 1, expit(-3 + 5 * U + 1 * X))
Y <- rnorm(N, 0.2 + 1.5 * U + 0.2 * X + 0.2 * A, 0.5)
Z <- rnorm(N, 2 * U + 0.5 * X)
W2 <- rexp(N, 0.1 + 1 * U + 0.1 * X)
D2 <- as.numeric(W2 < 5)
W2[D2 == 0] <- 5
W <- cbind(W1 = rnbinom(N, size = 25,
                        mu = exp(2.5 * U + 0.2 * X)),
           W2)
p2sls_result <- pci2s::p2sls.lm(Y = Y, A = A, X = X,
                         W = W, Z = Z,
                         nco_type = c("negbin", "ah"),
                         nco_args = list(list(offset = rep(0, N)),
                                         list(offset = rep(0, N),
                                              event = D2)))
p2sls_result$summary_first_stage
p2sls_result$summary_second_stage

# --------- p2sls.negbin
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


W[, 1] = W[, 2] + rnorm(nrow(W))
p2sls_result <- pci2s::p2sls.negbin(Y = Y, A = A, X = NULL,
                                    W = W, Z = Z,
                                    nco_type = c("linear", "linear"),
                                    nco_args = list(list(offset = rep(0, N)),
                                                    list(offset = rep(0, N))))


