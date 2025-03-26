# -------------------------------------------------------------------------------- #
#    Prepare for Proximal/Confounding Bridge using sparse PCA of genes as NCE/NCO  #
# Perform Sparse PCA in order to construct Negative Controls by combining individual genes #
# 2.1 - get sparse PCA PCs                                                         #
# 2.2 - choose A,Y,Z,W combinations                                                #
# Requires: prev saved normalized gene expression (HDF5)                           #
#           prev saved chromosome information                                      #
#           TF and non TF information                                              #
# Ouputs: (nothing) but saves                                                      #
#         CBGENE_AYZW in the rds file                                              #
#                         "<save_dir>/cbgenes/<AYZW_setting_name>/CBGENE_AYZW.rds" #
# -------------------------------------------------------------------------------- #
args = commandArgs(trailingOnly = TRUE)

library(PMA)


RUN_SPCA_CV = FALSE # run sPCA tuning with cross-validation (ow used hard coded defaults)
RUN_SPCA    = TRUE  # run sPCA (ow load in saved)
# number of sampled cells used for creating sparse PCs (NA for all)
N_subsample = 2000
# N_subsample = NA 

assertthat::assert_that(length(args) > 0, msg="must give arg for specifying device eg 'Rscript <filename>.R ubergenno'")

DEVICE = args[1]
source('../PATHS.R') # load in data_dir and save_dir, depending on DEVICE value

assertthat::assert_that(length(args) >= 2, msg="provide device and TEST_GENENUM (currently int <=4000)")
NUM_IMPORTANT_GENES = as.integer(args[2])


dir.create(sprintf('%s/spca/', save_dir), showWarnings = FALSE)

# ==================================================================================================
# =================== Load Normalized Gene Expr ====================================================
# ==================================================================================================
print(sprintf("[%s]    - Loading Normalized Gene Expression w/ %d Genes", Sys.time(), NUM_IMPORTANT_GENES))



# =================== loading top genes ===========================================================
print(sprintf("[%s]          loading top genes", Sys.time()))

# load normalized gene exp
h5file      = paste0(save_dir, "/gene.h5"); print(h5file)
reading_hd5file  = rhdf5::H5Fopen(name = h5file)
readin_gene_norm = reading_hd5file&'gene_norm'
gene_norm = readin_gene_norm[1:NUM_IMPORTANT_GENES, ] # dim = 4000 x 20729 = #important x #cells
rhdf5::h5closeAll()
invisible(gc(verbose=FALSE))


# =================== load list of genes that are Transcription Factors ===========================
print(sprintf("[%s]          loading list of TF genes", Sys.time()))
# read in xlsx sheet
tf_raw = readxl::read_xlsx(path = paste0(data_dir, "/../extra/transcriptionfactorlist.xlsx"),  # given from Kathryn over slack (suppl of a paper?)
                           sheet = 2) |> suppressMessages() # suppress messages on how they renamed columns
# clean up a bit
#      2nd col ('...2') is gene name, 4th col is TF indicator 'Is TF?'
#      First row is not data (2 rows of colnames)
tf = tf_raw[-1 , c('...2', 'Is TF?')]
colnames(tf) = c('gene_name', 'TF')
# table(tf$TF) # only levels are No, Yes
#        No  Yes 
#       1126 1639 
TF_names = tf |> dplyr::filter(TF == 'Yes') |> dplyr::pull(gene_name)


# =================== load targeted gene names ====================================================
grna_odm <- ondisc::read_odm(odm_fp      = paste0(data_dir, "/papalexi-2021/processed/grna_assignment/assignment_matrix.odm"),
                             metadata_fp = paste0(data_dir, "/papalexi-2021/processed/grna_assignment/metadata.rds"))

grnaTarget_names = ondisc::get_feature_covariates(grna_odm)[, 'target'] |> unique()

# =================== remove these genes from the top XX genes ====================================
geneRemove_names = c(TF_names, grnaTarget_names) |> unique()

# check
# length(geneRemove_names)
# length(TF_names)
# length(grnaTarget_names)

myGenenames     = readRDS(sprintf('%s/important_genes_name.rds', save_dir))
myGenenames_idx = readRDS(sprintf('%s/important_genes_idx.rds',  save_dir))

# keep idx genes not in remove list
keep_idx = sapply(X   = 1:length(myGenenames), 
                  FUN = function(i) {!(myGenenames[i] %in% geneRemove_names)})

sum(keep_idx) # number of genes kept
length(myGenenames) - sum(keep_idx) # number of genes removed

myGenenames_df = data.frame(importance_rank = 1:length(myGenenames),
                        gene            = myGenenames,
                        idx             = myGenenames_idx,
                        TF          = sapply(X   = 1:length(myGenenames), 
                                             FUN = function(i) {(myGenenames[i] %in% TF_names)}),
                        grna_target = sapply(X   = 1:length(myGenenames), 
                                             FUN = function(i) {(myGenenames[i] %in% grnaTarget_names)}))

write.csv(x = myGenenames_df, 
          file = sprintf('%s/important_genes_name_TF_Target.csv', save_dir), 
          row.names = FALSE) # top important genes noting the TF and gRNA targets


# myGenenames_df = data.frame(myGenenames, myGenenames_idx)

# dim(gene_norm)

set.seed(12345)
gene_norm_noTFTargets = gene_norm[myGenenames_df |> dplyr::filter((!TF) & (!grna_target)) |> dplyr::pull(importance_rank), ]
# if(N_subsample == ncol(gene_norm)) {
if(is.na(N_subsample)) {  # NA --> use all cells    
  gene_norm_SPC = t(gene_norm_noTFTargets) # all 21k cells
} else {
  gene_norm_SPC = t(gene_norm_noTFTargets[, sample(1:ncol(gene_norm), N_subsample)]) # sample 5k/XX out of 21k cells
}
rm(gene_norm_noTFTargets); gc()
# gene_norm_SPC = t(gene_norm[myGenenames_df |> dplyr::filter((!TF) & (!grna_target)) |> dplyr::pull(importance_rank),
#                             sample(1:ncol(gene_norm), N_subsample)]) # sample 5k out of 21k cells
dim(gene_norm_SPC)


# ==================================================================================================
# =================== Perform Sparse PCA ===========================================================
# ==================================================================================================

# =================== Use SPC.cv to choose tuning parameters: ======================================
# Use SPC.cv to choose tuning parameters:
if(RUN_SPCA_CV) {
  cv.out <- SPC.cv(gene_norm_SPC, 
                   
                   sumabsvs = c(seq(1.2, 5, len = 5), seq(6, floor(sqrt(ncol(gene_norm_SPC))), len = 5)),
                   orth=TRUE)
  print(cv.out)
  plot(cv.out)
  print(cv.out$bestsumabsv) # testing out on top 1000 genes sampled 5000 cells: 5 and 4.58 but default values are 1.2-5. Should be between 1 and sqrt(p). which is 1000 here... 31.6
  print(cv.out$bestsumabsv1se)
  
  saveRDS(cv.out, 
          sprintf('%s/spca/cvout.rds', save_dir)) # save cv res

  my_sumabsv = cv.out$bestsumabsv1se
  my_K = 60
} else {
  my_sumabsv = 5 # 8
  my_K = 60
}


if(RUN_SPCA) {
  # run sPCA and save result
  out.orth <- SPC(gene_norm_SPC,
                  sumabsv=my_sumabsv, # tuning parameter
                  # sumabsv=cv.out$bestsumabsv1se, # 33.5 not strong enough. v's not sparse and number of nonzero coefs too large
                  K=my_K, 
                  orth=TRUE)
  saveRDS(out.orth, 
          sprintf('%s/spca/outorth_sumabs=%.1f_K=%d_N=%d.rds', save_dir, my_sumabsv, my_K, N_subsample)) # save res
  
} else {
  # load in previously saved result
  # out.orth = readRDS(sprintf('%s/spca/outorth_sumabs=%.1f_K=%d_N=%d.rds', save_dir, 5, 60))
  out.orth = readRDS(sprintf('%s/spca/outorth_sumabs=%.1f_K=%d_N=%d.rds', save_dir, my_sumabsv, my_K, N_subsample))
} 


# =================== Not used methods to construct NCs ===========================================================
if(T) {
  print(out.orth, verbose=TRUE)


  out.orth
  names(out.orth)

  plot(out.orth$prop.var.explained, ylim = c(0, 1))
  out.orth$cnames
  out.orth$meanx
  out.orth$vpos

  # how many features in each Sparse PC:
  dim(out.orth$v); dim(out.orth$v != 0)  # num cells x num Sparse PCs
  plot(colSums(out.orth$v != 0), 
       main = '#Genes per Sparse PC', xlab = 'Sparse PC', ylab = '#Genes')

  # assign gene to NC idx without overlap between NCs (e.g. gene5 has nonzero coef in SparsePC2 and 3 --> assign only to 2nd NC)
  out.orth$v |> dim()
  abs(out.orth$v)

  out.orth$v[1, ]
  which.max(out.orth$v[1, ])




  # by gene, how many PCs does this gene have nonzero coefs in
  plot(rowSums(out.orth$v != 0), 
       main = '#nonzero Sparse PCs per Gene', xlab = 'importance rank (ish)', ylab = '#Sparse PCs')


  out.orth$v == 0
  out.orth$v[1:10, 1:2]
  colSums(out.orth$v[, 1:2] * out.orth$v[, 1:2])

  # map each of the Sparse PCs to a NC. Then, remove repeats based on coef abs val
  which.absmax.nonzero <- function(vec) {
    # if the maximum abs value is 0, then return NA
    i = which.max(abs(vec))
    if(vec[i] == 0) {
      return(NA)
    } else {
      return(i)
    }
  }
  gene_NC_noupdate = apply(out.orth$v, MARGIN = 1, FUN = which.absmax.nonzero) # which.max still returns something when ties...







  # i think we just need to increase the sparsity parameter
  # inefficient but updates coef vector after assigning the gene
  cur_v = out.orth$v |> abs()
  gene_NC = rep(NA, nrow(cur_v))
  for(g_idx in 1:nrow(cur_v)) {
    # find which PC has max coef
    maxcoef = max(cur_v[g_idx, ]) 
    
    # if nonzer coef anywhere, assign to a NC
    if(maxcoef != 0) {
      whichNC = which.max(cur_v[g_idx, ])  
      gene_NC[g_idx] = whichNC
      
      # zero out other PCs coefs and normalize. hmm should also update assigned PC? (move the weight? no for now)
      othernonzero_PC = setdiff(which(cur_v[g_idx, ] != 0), whichNC)
      for(nonzero_PC in othernonzero_PC) {
        cur_v[g_idx, nonzero_PC] = 0
        cur_v[     , nonzero_PC] = cur_v[, nonzero_PC] / sum(cur_v[, nonzero_PC]**2)
      }
    }
    
    
  }
  # checking should be just 1 or 0: by gene, how many 'PC's does this gene have nonzero coefs in
  plot(rowSums(cur_v != 0), 
       main = '#nonzero Updated Sparse PCs per Gene', xlab = 'importance rank (ish)', ylab = '#Sparse PCs')


  plot(gene_NC_noupdate)
  plot(gene_NC)
  table(gene_NC); table(gene_NC_noupdate)


  gene_NC |> table() |> as.numeric() |> hist(breaks = seq(1, to = 100, by = 1))

  plot(data.frame(table(gene_NC))[,1],
       data.frame(table(gene_NC))[,2])


  out.orth$v
  (out.orth$v[, 1] != 0) |> sum()

}


# =================== Construct NCs ================================================================
gene_norm_noTFTargets = gene_norm[myGenenames_df |> dplyr::filter((!TF) & (!grna_target)) |> dplyr::pull(importance_rank), ]
# construct NCs as the actual loadings
NC_loadings = t(gene_norm_noTFTargets) %*% out.orth$v
# gene_norm_noTFTargets |> dim()
# out.orth$v|> dim()
NC_loadings |> dim()

# check orthogonality (not exactly Identity matrix bc SPC used sample)
t(NC_loadings) %*% NC_loadings / nrow(NC_loadings)



heatmap( t(NC_loadings) %*% NC_loadings / nrow(NC_loadings), 
         Rowv=NA, Colv=NA, col = heat.colors(256),  margins=c(5,10))


# OR construct NCs as the average 
NC_avg = matrix(NA, 
                nrow = ncol(gene_norm_noTFTargets), # number of cells 
                ncol = ncol(out.orth$v)) # number of PCs
for(j in 1:ncol(out.orth$v)) {
  NC_avg[ ,j] = gene_norm_noTFTargets[which(gene_NC == j), ,drop=FALSE] |> colMeans()
}

NC_avg  |> dim()


# =================== Save =========================================================================
print(sprintf("[%s]         Saving NCs", Sys.time()))
saveRDS(NC_loadings, 
        sprintf('%s/spca/NCloadings_sumabs=%.1f_K=%d_N=%d.rds', save_dir, my_sumabsv, my_K, N_subsample)) # save 
saveRDS(NC_avg, 
        sprintf(     '%s/spca/NCavg_sumabs=%.1f_K=%d_N=%d.rds', save_dir, my_sumabsv, my_K, N_subsample)) # save







# ==================================================================================================
# =================== Extra: Code Testing ==========================================================
# ==================================================================================================
if(F) {
  
  

########################## TESTING SPARSE PCA PACKAGES

install.packages('PMA')
install.packages('sparsepca')

###### Package: PMA 
###### https://cran.r-project.org/web/packages/PMA/PMA.pdf
###### authors: Witten et al.,

help(PMA::SPC)
help(PMA)
library(PMA)
help(SPC)


set.seed(1)
u <- matrix(c(rnorm(50), rep(0,150)),ncol=1)
v <- matrix(c(rnorm(75),rep(0,225)), ncol=1)
x <- u%*%t(v)+matrix(rnorm(200*300),ncol=300)
# Perform Sparse PCA - that is, decompose a matrix w/o penalty on rows
# and w/ L1 penalty on columns
# First, we perform sparse PCA and get 4 components, but we do not
# require subsequent components to be orthogonal to previous components
out <- SPC(x,sumabsv=3, K=4)
print(out,verbose=TRUE)
# We could have selected sumabsv by cross-validation, using function SPC.cv
# Now, we do sparse PCA using method in Section 3.2 of WT&H(2008) for getting
# multiple components - that is, we require components to be orthogonal
out.orth <- SPC(x,sumabsv=3, K=4, orth=TRUE)
print(out.orth,verbose=TRUE)
par(mfrow=c(1,1))
plot(out$u[,1], out.orth$u[,1], xlab="", ylab="")
# Note that the first components w/ and w/o orth option are identical,
# since the orth option only affects the way that subsequent components
# are found
print(round(t(out$u)%*%out$u,4)) # not orthogonal
print(round(t(out.orth$u)%*%out.orth$u,4)) # orthogonal

# Use SPC.cv to choose tuning parameters:
cv.out <- SPC.cv(x)
print(cv.out)
plot(cv.out)
out <- SPC(x, sumabsv=cv.out$bestsumabsv)
print(out)
# or we could do
out <- SPC(x, sumabsv=cv.out$bestsumabsv1se)
print(out)





###### Package: sparsepca
###### https://cran.r-project.org/web/packages/sparsepca/sparsepca.pdf
###### authors: Erichson et al., 
###### another implementation??

}
