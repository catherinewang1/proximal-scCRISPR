#' Performs pca on Z, W separately and puts back
#' into dataframe format w/ same names used for OCB/prox functions
#'
#'
#' @param df
#' @param r_z (integer) rank/#pcs for pca on Z
#' @param r_w (integer) rank/#pcs for pca on W
#' @param colnames
#' @return dataframe with original A,Us,
#'               and with new PCA Zs, Ws
getdfPCA <- function(df, r_z, r_w) {
  # r_z = 6; r_w = 6
  Z_colnames = grep('Z', colnames(df), value = T)
  W_colnames = grep('W', colnames(df), value = T)

  Z = df[, Z_colnames]
  W = df[, W_colnames]

  Z_pca = prcomp(Z, center = TRUE, scale. = TRUE)
  W_pca = prcomp(W, center = TRUE, scale. = TRUE)

  Z_ = Z_pca$x[,1:r_z, drop=FALSE]
  W_ = W_pca$x[,1:r_w, drop=FALSE]
  colnames(Z_) = paste0('Z',1:ncol(Z_))
  colnames(W_) = paste0('W',1:ncol(W_))

  df_ = cbind(df[, !(colnames(df) %in% c(Z_colnames, W_colnames))],
              Z_,
              W_)

  return(df_)
}


#' Get top r PCA rotations of matrix M 
#' @param M (matrix)
#' @param r (rank)
getPCArotations <- function(M, r) {
  M_pca = prcomp(M, center = TRUE, scale. = TRUE)
  M_    = M_pca$x[,1:r, drop=FALSE]
  colnames(M_) = paste0('PC',1:ncol(M_))
  return(M_)
}
