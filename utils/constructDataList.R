



# ============= Helper Functions
# polybasis_transform: constructs a n x degree matrix
# poly_eval: evaluates input vals


#' Constructs a n x degree matrix for each Hermite Polynomial 
#' basis j function evaluated at value i  b_j(vals[i])
#' @param vals (numeric vector) values to evaluate polynomials at 
#' @param degree (integer) highest degree of Hermite polynomial
#' @param type (string) polynomial basis type (hermite, simple [1, X^2, X^3, ...])
polybasis_transform <- function(vals, degree, type = 'hermite') {
  switch(type,
    hermite = {
      inner_func <- function(d) {EQL::hermite(x=vals, n=d, prob=TRUE)}
    },
    simple  = {
      inner_func <- function(d) {vals**d}
    },
    {print('bad [type] input (hermite, simple)'); return()})

  return(sapply(X=0:degree, FUN=inner_func)) 
}

#' Evaluates the input vals at the specified Hermite Polynomial 
#' basis degree with the specified alpha coefficients
#' result[i] = \sum_k=1^degree \alpha[k] poly_k (val[i])
#' @param vals (numeric vector) values to evaluate polynomials at 
#' @param degree (integer) highest degree of Hermite polynomial
#' @param alpha (numeric vector) coefficients
#' @return result[i] = \sum_k=1^degree \alpha[k] poly_k (val[i])
poly_eval <- function(vals, degree, alpha) {
  assertthat::assert_that(length(alpha) == degree + 1, msg='alpha must be a tall vector of length == degree + 1 (for constant 1)')
  polybasis_transform(vals, degree) %*% alpha
}


#' construct a n x degree matrix for cos bases
#' y = cos( B (x + phase shift))
#' where B is the coefficient that corresponds to the specified 
#' period = 2 pi / B <==> B = 2 pi / period
#' @param vals (numeric vector) values to evaluate basis fns at 
#' @param phases (vector) of phase shifts (for basis functions)
#' @param periods (vector) of period lengths (for basis functions)
#' @examples
#' v = seq(from = -5, to = 5, length.out = 100)
#'     
#' # changing phases
#' vcos = cosbasis_transform(v, phases = c(0, 1, 2, 3), periods = 2*base::pi )
#' head(vcos)
#' plot(v, vcos[,1], type = 'l', col = "orange")
#' lines(v, vcos[,2], col = "red")
#' lines(v, vcos[,3], col = "purple")
#' lines(v, vcos[,4], col = "blue")
#' 
#' # changing periods
#' vcos = cosbasis_transform(v, phases = c(0), periods = 2*base::pi * c(.5, 1, 2, 4) )
#' head(vcos)
#' plot(v, vcos[,1], type = 'l', col = "orange")
#' lines(v, vcos[,2], col = "red")
#' lines(v, vcos[,3], col = "purple")
#' lines(v, vcos[,4], col = "blue")
cosbasis_transform <- function(vals, phases, periods) {
   X = matrix(NA, nrow = length(vals), ncol = length(periods) * length(phases))
   cur_col = 1
   for(phase in phases) {
    for(period in periods) {
      X[,cur_col] = cos((2 * base::pi / period) * (vals + phase))
      cur_col = cur_col + 1 # ; print(cur_col)
    }
   }
   return(X)
}




#' Constructing Data List for Poly OCB variation 1
#' Transform each Zk and Wj with each basis function b's and d's respectively
#' with NO interactions
#' 
#' B and D's are centered and scaled (except for intercept)
#' Need matrices/vectors/numerics
#' Y1, Y0: response
#' B, B1, B0: conditional moment basis functions evaluated at Zs
#' D, D1, D0: function class basis functions evaluated at Ws
#' p1, p0: proportion of samples in each treatment group
#' @param df (dataframe)
#' @param b_degree (integer) max degree of conditional basis functions
#' @param h_degree (integer) max degree of h (confounding bridge) basis fns
#' @param type (string) type of basis functions (e.g. hermite, simple, )
constructDataListv1 <- function(df, b_degree, h_degree, type='hermite') {
  Znames = grep('Z', colnames(df), value = TRUE); numZ = length(Znames)
  Wnames = grep('W', colnames(df), value = TRUE); numW = length(Wnames)
  
  # define values/dataframes for computing
  df1 = df |> filter(A == 1); df0 = df |> filter(A == 0) # subset df by trtmnt
  p1 = mean(df$A == 1); p0 = mean(df$A == 0)             # prop in each trtmnt
  Y1 = df1$Y; Y0 = df0$Y                             # response of each trtmnt
  
  # Construct B (basis b's evaluated at Zs)
  B = matrix(NA, nrow = nrow(df), ncol = 1 + length(Znames) * b_degree) 
  B[, 1] = 1
  for(j in 1:length(Znames)) {
    B_zname = polybasis_transform(df[, Znames[j]], b_degree, type=type) # match up to _th degree
    # print(((j-1)*b_degree + 2):(j*b_degree + 1)) # this is what happens when indexing by 1 and inclusive of last point in range...
    # B[,((j-1)*b_degree + 2):(j*b_degree + 1)] = B_zname[,-1] # add cols to B
    B[,((j-1)*b_degree + 2):(j*b_degree + 1)] = scale(B_zname[,-1], center = TRUE, scale = TRUE) # add cols to B (centered&scaled)
  }
  # B[,2:ncol(B)] = scale(B[,2:ncol(B)], center = TRUE, scale = TRUE) # center and standardize each column (except 1)
  B1 = B[which(df$A == 1), ]; B0 = B[which(df$A == 0), ]
  
  # D  = polybasis_transform( df$W1, h_degree, type=type) # allow flexibility up to _ degrees 
  # D1 = polybasis_transform(df1$W1, h_degree, type=type) # allow flexibility up to _ degrees 
  # D0 = polybasis_transform(df0$W1, h_degree, type=type) # allow flexibility up to _ degrees 
  
  D = matrix(NA, nrow = nrow(df), ncol = 1 + length(Wnames) * h_degree) 
  D[, 1] = 1
  for(j in 1:length(Wnames)) {
    D_zname = polybasis_transform(df[, Wnames[j]], h_degree, type=type) # match up to _th degree
    # print(((j-1)*h_degree + 2):(j*h_degree + 1)) # this is what happens when indexing by 1 and inclusive of last point in range...
    # D[,((j-1)*h_degree + 2):(j*h_degree + 1)] = D_zname[,-1] # add cols to B
    D[,((j-1)*h_degree + 2):(j*h_degree + 1)] = scale(D_zname[,-1], center = TRUE, scale = TRUE)  # add cols to B
  }
  # D[,2:ncol(D)] = scale(D[,2:ncol(D)], center = TRUE, scale = TRUE)
  D1 = D[which(df$A == 1), ]; D0 = D[which(df$A == 0), ]
  
  return(list(Y1=Y1, Y0=Y0,
              # B=B, 
              B1=B1, B0=B0,
              # D=D, 
              D1=D1, D0=D0, 
              p1=p1, p0=p0))
}



#' Constructing Data List for Poly OCB variation 2 (with Cos bases)
#' Transform each Zk and Wj with each basis function b's and d's respectively
#' with NO interactions
#' 
#' B and D's are centered and scaled (except for intercept)
#' Need matrices/vectors/numerics
#' Y1, Y0: response
#' B, B1, B0: conditional moment basis functions evaluated at Zs
#' D, D1, D0: function class basis functions evaluated at Ws
#' p1, p0: proportion of samples in each treatment group
#' @param df (dataframe)
#' @param b_degree (integer) max degree of conditional basis functions
#' @param h_degree (integer) max degree of h (confounding bridge) basis fns
#' @param b_phases  (vector) of phase shifts (numeric) for b basis functions
#' @param b_periods (vector) of period lengths (numeric) for b basis functions
#' @param h_phases  (vector) of phase shifts (numeric) for d basis functions
#' @param h_periods (vector) of period lengths (numeric) for d basis functions
#' @param type (string) type of poly basis functions (e.g. hermite, simple, )
constructDataListv2 <- function(df, 
                                b_degree=0, h_degree=0, 
                                b_phases=NULL, b_periods=NULL,
                                h_phases=NULL, h_periods=NULL, 
                                type='hermite') {
  Znames = grep('Z', colnames(df), value = TRUE); numZ = length(Znames)
  Wnames = grep('W', colnames(df), value = TRUE); numW = length(Wnames)
  
  # define values/dataframes for computing
  df1 = df |> filter(A == 1); df0 = df |> filter(A == 0) # subset df by trtmnt
  p1 = mean(df$A == 1); p0 = mean(df$A == 0)             # prop in each trtmnt
  Y1 = df1$Y; Y0 = df0$Y                             # response of each trtmnt
  
  B = matrix(NA, nrow = nrow(df), ncol = 1) 
  B[, 1] = 1
  # cur_col = 2
  for(j in 1:length(Znames)) {
    B_zname_poly = NULL
    B_zname_cos  = NULL
    
    if(b_degree > 0) {
      B_zname_poly = polybasis_transform(df[, Znames[j]], degree=b_degree, type=type)[,-1]
    }
    
    if(!is.null(b_phases) && !is.null(b_periods)) {
      B_zname_cos = cosbasis_transform(vals=df[, Znames[j]], phases=b_phases, periods=b_periods)
      
    }
    B = cbind(B, B_zname_poly, B_zname_cos)
    # B = cbind(B, scale(cbind(B_zname_poly, B_zname_cos), center = TRUE, scale = TRUE))
    # B[,cur_col:(ncol(B_zname) + cur_col - 2)] = scale(B_zname[,-1], center = TRUE, scale = TRUE)
    # B[,((j-1)*b_degree + 2):(j*b_degree + 1)] = scale(B_zname[,-1], center = TRUE, scale = TRUE) # add cols to B (centered&scaled)
    # cur_col = cur_col + ncol(B_zname)
  }
  B[,-1] = scale(B[,-1], center=TRUE, scale=TRUE) 
  B1 = B[which(df$A == 1), ]; B0 = B[which(df$A == 0), ]
  
  # Construct D (basis d's evaluated at Ws)
  D = matrix(NA, nrow = nrow(df), ncol = 1) 
  D[, 1] = 1
  for(j in 1:length(Wnames)) {
    D_zname_poly = NULL; D_zname_cos = NULL
    if(h_degree > 0) {
      D_zname_poly = polybasis_transform(df[, Wnames[j]], degree=h_degree, type=type)[,-1] 
    }

    if(!is.null(h_phases) && !is.null(h_periods)) {
      D_zname_cos = cosbasis_transform(vals=df[, Wnames[j]], phases=h_phases, periods=h_periods)
    }
    D = cbind(D, D_zname_poly, D_zname_cos)
  }
  D[,-1] = scale(D[,-1], center=TRUE, scale=TRUE)
  D1 = D[which(df$A == 1), ]; D0 = D[which(df$A == 0), ]
  
  return(list(Y1=Y1, Y0=Y0,
              B=B,
              B1=B1, B0=B0,
              D=D,
              D1=D1, D0=D0, 
              p1=p1, p0=p0))
}
