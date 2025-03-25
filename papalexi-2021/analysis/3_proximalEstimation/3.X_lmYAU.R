# ---------------------------------------------------------------------------- #
#                   Do lmYAU (using known confounders)                         #
# 3.X lmYAU Fits                                                               #
# Requires: prev saved normalized gene expression (HDF5)  and list of AY pairs #
# Ouputs: (nothing) but saves                                                  #
#         ATE_lmYAU.csv                                                        #
#         "<save_dir>/spca/cbgenes/<AYZW_setting_name>/lmYAU/ATE_lmYAU.csv"    #
# ---------------------------------------------------------------------------- #
args = commandArgs(trailingOnly = TRUE)
args = c('laptop', 'allPerturbations')

require(assertthat) # for some assert statements
library(Matrix)
library(dplyr)
library(ggplot2)
library(cowplot)
library(future.apply)
options(future.globals.maxSize= 850*1024^2) #1st num is MB
plan(multisession, workers = 4)
# plan(sequential)

theme_set(theme_cowplot() +
            theme(plot.title = element_text(hjust = .5),
                  plot.subtitle = element_text(hjust = .5),
                  strip.background = element_rect(color = 'black', fill = 'white')))



assertthat::assert_that(length(args) > 0, msg="must give arg for specifying device eg 'Rscript <filename>.R ubergenno'")

DEVICE = args[1]
source('../PATHS.R') # load in data_dir and save_dir and CODE_DIR, depending on DEVICE value

assertthat::assert_that(!is.null(data_dir), msg='first arg must be: laptop, desktop, or ubergenno')

assertthat::assert_that(length(args) > 1, msg="must give arg for specifying chosen AYZW name 'Rscript <filename>.R ubergenno C'")
AYZW_setting_name = args[2]

dir.create(sprintf('%s/spca/cbgenes/%s/lmYAU', save_dir, AYZW_setting_name), recursive = TRUE, showWarnings = FALSE)


# =================== Start ====================================================
print(sprintf("[%s] START: lmYAU Estimate", Sys.time()))


# load chosen AYZW names
AY   = read.csv(sprintf('%s/spca/cbgenes/%s/AY.csv', save_dir, AYZW_setting_name))

# load gene importance info
imp_gene_names = readRDS(sprintf('%s/important_genes_name.rds', save_dir))
imp_gene_idx   = readRDS(sprintf('%s/important_genes_idx.rds',  save_dir))
imp_gene = data.frame(gene     = imp_gene_names,
                      gene_idx = imp_gene_idx,
                      gene_imp_rank = 1:length(imp_gene_names))

# Load choosing AY (ZW) settings
AYZW_setting = readRDS(sprintf('%s/spca/cbgenes/%s/AYZW_setting.rds', save_dir, AYZW_setting_name))
NUM_IMPORTANTGENES = AYZW_setting$MAX_Y_IMPORTANCE

# create gene and grna ondisc managers
gene_odm <- ondisc::read_odm(odm_fp      = paste0(data_dir, "/papalexi-2021/processed/gene/expression_matrix.odm"),
                             metadata_fp = paste0(data_dir, "/papalexi-2021/processed/gene/metadata.rds"))
grna_odm <- ondisc::read_odm(odm_fp      = paste0(data_dir, "/papalexi-2021/processed/grna_assignment/assignment_matrix.odm"),
                             metadata_fp = paste0(data_dir, "/papalexi-2021/processed/grna_assignment/metadata.rds"))

# load grna assignments (load all into memory)
grna = grna_odm[[,1:ncol(grna_odm)]] # |> as.matrix() # ~110 x 20729 = #grnas x #cells
grna_rownames = grna_odm |> ondisc::get_feature_covariates() |> rownames()

# load normalized gene exp
h5file      = paste0(save_dir, "/gene.h5"); print(h5file)
reading_hd5file  = rhdf5::H5Fopen(name = h5file)
readin_gene_norm = reading_hd5file&'gene_norm'
gene_norm = readin_gene_norm[1:NUM_IMPORTANTGENES, ] # eg dim = 4000 x 20729 = #important x #cells
rhdf5::h5closeAll()
invisible(gc(verbose=FALSE))




# =================== Define File Specific Functions ======================================
print(sprintf("[%s]    - Define File Specific Functions", Sys.time()))


# All the NT grna idx
NT_names = grna_odm |> ondisc::get_feature_covariates() |>
  filter(target_type == 'non-targeting') |> rownames()
NT_idx = which(apply(X = grna_odm[[NT_names, ]], MARGIN = 2, FUN = sum) > 0)

cell_covariates = ondisc::get_cell_covariates(gene_odm)

#' get the importance rank of specified gene_name
get_importance_rank <- function(gene_name) {
  which(imp_gene_names == gene_name)
}

#' Construct dataframe of columns A,Y,U (treatment, response, covariates)
#' Requires outside named objects:
#' AY, NT_idx,
#' grna_rownames, grna
#' gene_norm
#' get_importance_rank,
#' @param AY_idx (integer) from 1 to nrow(AY) idx of the AY perturbation pair in the AY df
#' @example a = format_AYZW_inner(AY_idx = 2); head(a)
format_AYU <- function(AY_idx) {
  AY_row = AY[AY_idx, ]
  # Get A
  # using grna all loaded into memory
  # -------------------------------------------
  # idx of all 'treated' cells
  A_grna_idx = which(grna_rownames == AY_row$A)
  A_idx = which(as.logical(grna[A_grna_idx, 1:ncol(grna)]))  
  
  # subset cells of 'treated' (w A grna) and 'control' (NT gran)
  A = grna[A_grna_idx, c(A_idx, NT_idx)]
  
  # Get Y
  # using already loaded in gene_norm (faster)
  # -------------------------------------------    
  Y = gene_norm[get_importance_rank(AY_row$Y),     c(A_idx, NT_idx)]
  
  # Get U
  # -------------------------------------------
  U = cell_covariates[c(A_idx, NT_idx), ]
  
  # Assemble together
  # -------------------------------------------
  dfAYU = data.frame(A = as.vector(A), 
                     Y = Y,
                     U)
  
  return(dfAYU)
}

#' Perform linear regression of Y on A + <covariates> 
#' @param AY_idx (integer) idx of the AY perturbation pair in the AY df
#' @param covariates (vector) vector of strings of covariates used for linear regression 
#'        if NULL (default), then use all covariates with no tranformation
#'        e.g. not including n_nonzero: c('n_umis', 'lane', 'bio_rep', 'phase', 'p_mito')
#'             and will be input into the formula like this: 
#'             paste0(c('n_umis', 'lane', 'bio_rep', 'phase', 'p_mito'), collapse = ' + ')
get_ATE_est_lmYAU <- function(AY_idx, covariates=NULL) {
  # AY_idx = 8755; covariates = c('n_umis', 'lane', 'bio_rep', 'phase', 'p_mito')
  
  # format df
  dfAYU = format_AYU(AY_idx) # head(dfAYU); dim(dfAYU)
  
  # perform linear regression without covariates (naive)
  lmYA = lm('Y ~ A', data = dfAYU)
  lmYA_sum = summary(lmYA)
  res = data.frame(method      = 'lmYA',
                   method_type = 'linearmodel',
                   numNC       = NA,
                   basis       = NA,
                   ATE  = lmYA_sum$coefficients['ATRUE',   'Estimate'],
                   se   = lmYA_sum$coefficients['ATRUE', 'Std. Error'],
                   pval = lmYA_sum$coefficients['ATRUE',   'Pr(>|t|)'])
  
  
  # perform linear regression with known covariates
  if(is.null(covariates)) {
    myFormula = 'Y ~ A + .'
  } else {
    myFormula = paste0('Y ~ A + ',  paste0(covariates, collapse = ' + '))
  }
  
  lmYAU = lm(myFormula, data = dfAYU)
  lmYAU_sum = summary(lmYAU)
  
  res = rbind(res, 
              data.frame(method      = 'lmYAU',
                         method_type = 'linearmodel',
                         numNC       = NA,
                         basis       = NA,
                         ATE  = lmYAU_sum$coefficients['ATRUE',   'Estimate'],
                         se   = lmYAU_sum$coefficients['ATRUE', 'Std. Error'],
                         pval = lmYAU_sum$coefficients['ATRUE',   'Pr(>|t|)']))
  
  return(res)
}


get_ATE_est_lmYAU(1000, covariates = NULL)



get_ATE_est_catcherror <- function(AY_idx, covariates=NULL) {
  res_df = tryCatch({get_ATE_est_lmYAU(AY_idx=AY_idx, covariates=covariates)},
                    error = function(cond) {
                      # message(sprintf('Error est CondMomentOCBOS with %s', 
                      #                 gammaSetting)) 
                      return(NULL)
                    })
  # if errored, return NULL
  if(is.null(res_df)) {
    print(sprintf('AY_idx = %d errored: returned null', AY_idx))
    return(data.frame(method      = NA,
                      method_type = NA,
                      numNC       = NA,
                      basis       = NA,
                      ATE         = NA,
                      se          = NA,
                      pval        = NA))
  } else {
    return(res_df)
  }
}

ATEargs = data.frame(AY_idx = 1:nrow(AY))
# NUMROWS = 10
NUMROWS = nrow(ATEargs)
# whichROWS = 300:nrow(ATEargs)
# whichROWS = 1000:8000
whichROWS = 1:NUMROWS
# whichROWS = 6000:NUMROWS

# 6000-END


# # =================== Get ATEs (parallel) ====================================
print(sprintf("[%s]    - Get ATEs (parallel)", Sys.time()))

myCovariateList = list(none = NULL,
                       all  = c('n_nonzero', 'n_umis', 'lane', 'bio_rep', 'phase', 'p_mito'),
                       variation1)

cell_covariates |> colnames()
# myCovariates = c('n_nonzero', 'n_umis', 'lane', 'bio_rep', 'phase', 'p_mito')
# myCovariates = c('n_umis', 'lane', 'bio_rep', 'phase', 'p_mito')
myCovariates = c('log(n_umis)', 'lane', 'bio_rep', 'phase', 'p_mito')
# myCovariates = NULL

t0 = Sys.time()
ATE_par = future.apply::future_mapply(get_ATE_est_catcherror,
                                      MoreArgs = list(covariates = myCovariates),
                                      AY_idx = ATEargs[whichROWS, 1], # ATEargs[1:NUMROWS, 1],
                                      future.globals = TRUE,
                                      future.seed = 56789)
t1 = Sys.time()
print(sprintf("[%s]        - %2.2f", Sys.time(), (t1 - t0)))


# format
ATE = NULL
for(j in 1:ncol(ATE_par)) {
  ATE = rbind(ATE,
              ATE_par[,j] |> data.frame() |> mutate(AY_idx = j))
}

write.csv(x = ATE,
          file = sprintf('%s/spca/cbgenes/%s/lmYAU/ATE_lmYAU.csv', save_dir, AYZW_setting_name),
          row.names = FALSE)


# =================== END ======================================================
print(sprintf("[%s] END", Sys.time()))



# Some plots just for lmYAU


# ATE = read.csv(sprintf('%s/spca/cbgenes/%s/lmYAU/ATE_lmYAU.csv', save_dir, AYZW_setting_name))
df = merge(AY |> mutate(AY_idx = 1:nrow(AY)), 
           ATE, 
           by = 'AY_idx')
head(df); dim(df)



# hist(df |> filter(type == 'negative' &  method == 'lmYAU') |> pull(pval), breaks = seq(0, 1, length.out = 100))
# hist(df |> filter(type == 'negative' &  method == 'lmYA') |> pull(pval), breaks = seq(0, 1, length.out = 100))
# df |> arrange(pval)


myType = 'negative'
# df |> filter(type == myType & method == 'lmYA' )  |> arrange(pval) |> head()

lmYApvals  = df |> filter(type == myType & method == 'lmYA' ) |> pull(pval) |> sort()
lmYAUpvals = df |> filter(type == myType & method == 'lmYAU') |> pull(pval) |> sort()



plot_df = rbind(data.frame(EstMethod = 'lmYA',
                           pvals = lmYApvals,
                           theoretical = ppoints(length(lmYApvals))),
                data.frame(EstMethod = 'lmYAU',
                           pvals = lmYAUpvals,
                           theoretical = ppoints(length(lmYAUpvals)))) |> 
  data.frame() |>
  mutate(EstMethod = factor(EstMethod, levels = c('lmYA', 'lmYAU'))) |>
  mutate(pvals_trans = -log(pvals, base = 10),
         theoretical_trans = -log(theoretical, base = 10))



# qqplot of pvals vs unif
p1 = ggplot(plot_df, aes(sample=pvals, group = EstMethod, color = EstMethod)) +
  geom_abline(slope = 1, intercept = 0, color = 'purple3', linewidth = 1, alpha = .8) +
  # stat_qq_line(distribution = stats::qunif, color = 'purple', line.p = c(0, 1)) +
  stat_qq(distribution = stats::qunif, size = 1, alpha = .8) +
  scale_color_brewer(palette = "Set2") +
  labs(title = 'QQ Plot of p-values against Unif(0,1)',
       x = 'theoretical Unif(0,1)', y = 'sample p-values') +
  theme(legend.position = 'inside', 
        legend.position.inside = c(.1, .75))
p1
ggsave(plot = p1,
       filename = sprintf('%s/spca/cbgenes/%s/lmYAU/pval_qqplot.pdf', save_dir, AYZW_setting_name),
       height = 5, width = 5)


# qqplot of transformed pvals vs unif
p2 = ggplot(plot_df, 
            aes(# x = theoretical, y = pvals,
              x = theoretical_trans, y = pvals_trans,
              group = EstMethod, color = EstMethod)) +
  geom_point() +
  # stat_qq_line(distribution = stats::qunif, color = 'purple', line.p = c(0, 1)) +
  # stat_qq(distribution = stats::qunif, size = 1, alpha = .8) +
  geom_abline(slope = 1, intercept = 0, color = 'purple3', linewidth = 1, alpha = .8) +
  scale_color_brewer(palette = "Set2") +
  labs(title = '-log(p-values) vs -log(Unif(0,1))',
       x = 'theoretical -log(Unif(0,1))', y = 'sample -log(pval)') +
  theme(legend.position = 'inside', 
        legend.position.inside = c(.1, .75))
p2
ggsave(plot = p2,
       filename = sprintf('%s/spca/cbgenes/%s/lmYAU/pval_qqplot_trans.pdf', save_dir, AYZW_setting_name),
       height = 5, width = 5)

