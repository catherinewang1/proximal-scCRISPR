

suppressPackageStartupMessages(library(crayon))
cat(crayon::blue(' '))


# cat(blue('hi'))
# crayon::reset('resetting')


myPrintColor <- function(txt, col=36) {
  # change number from 29:47 # for(col in 29:47){ cat(paste0("\033[0;", col, "m",txt,"\033[0m","\n"))}
  cat(paste0("\033[0;", col, "m",txt,"\033[0m","\n"))
}

myPrintColor('hi')
myPrintColor(sprintf("[%s] START: 7_supergenes/1.0_createsupergenes.R", Sys.time()))
myPrintColor(sprintf("[%s]        - performing normalization", Sys.time()))
myPrintColor(sprintf("[%s]        - performing normalization", Sys.time()))
message('hi')
myPrintColor(sprintf("[%s]        - performing normalization", Sys.time()))
cat('    \n \n  ')
myPrintColor(sprintf("[%s]        - performing normalization", Sys.time()))
myPrintColor(sprintf("[%s]          performing normalization", Sys.time()))
myPrintColor(sprintf("[%s]          loading list of TF genes", Sys.time()))
# read in xlsx sheet
# data_dir = "C:/Users/Cathe/Documents/School/Genetic Hypothesis Testing using Negative Controls/genData/papalexi"
# tf_raw = readxl::read_xlsx(path = paste0(data_dir, "/../extra/transcriptionfactorlist.xlsx"),  # given from Kathryn over slack (suppl of a paper?)
#                            sheet = 2) |> suppressMessages() # suppress messages on how they renamed columns
# # clean up a bit
# #      2nd col ('...2') is gene name, 4th col is TF indicator 'Is TF?'
# #      First row is not data (2 rows of colnames)
# tf = tf_raw[-1 , c('...2', 'Is TF?')]
# colnames(tf) = c('gene_name', 'TF')
# # table(tf$TF) # only levels are No, Yes
# #        No  Yes 
# #       1126 1639 
# TF_names = tf |> dplyr::filter(TF == 'Yes') |> dplyr::pull(gene_name)
myPrintColor(sprintf("[%s]        - performing normalization", Sys.time()))
myPrintColor(sprintf("[%s]        - performing normalization", Sys.time()))

message('hi')



# Reproducible code

remove.packages('pci2s')
devtools::install_github('https://github.com/KenLi93/pci2s')
remotes::install_github('https://github.com/KenLi93/pci2s')

library('pci2s')
# Example Code for function `pci.lm`
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
pci_result <- pci.lm(Y = Y, A = A, X = X,
                     W = W, Z = Z,
                     nco_type = c("negbin", "ah"),
                     nco_args = list(list(offset = rep(0, N)),
                                     list(offset = rep(0, N),
                                          event = D2)))
# Errors with the message: 
# Error in colnames(W1) : object 'W1' not found


