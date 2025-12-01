# ---------------------------------------------------------------------------- #
# just get the p-val from lmYA estimates
# ---------------------------------------------------------------------------- #

args = commandArgs(trailingOnly = TRUE)
# args = c('laptop', '2000', 'C_ubergenno')
# args = c('laptop', '2000', 'A')
# args = c('laptop', '2000', 'D') # lots of negative AY tests
# args = c('laptop', '2000', 'E_PCA')

require(assertthat) # for some assert statements
library(Matrix)
library(dplyr)
library(tibble)
library(ggplot2)    # plotting
library(cowplot)


library(future.apply)
# options(future.globals.maxSize= 850*1024^2) #1st num is MB
# plan(multisession, workers = 8)
plan(sequential)


theme_set(theme_cowplot() +
            theme(plot.title = element_text(hjust = .5),
                  plot.subtitle = element_text(hjust = .5)))


assertthat::assert_that(length(args) > 0, msg="must give arg for specifying device eg 'Rscript <filename>.R ubergenno'")



DEVICE = args[1]
source('../PATHS.R') # load in data_dir and save_dir and CODE_DIR, depending on DEVICE value

# DEVICE = args[1]
# # location of papalexi-2021 folder
# data_dir = switch(DEVICE,
#                   'laptop'='C:/Users/Cathe/Documents/School/Genetic Hypothesis Testing using Negative Controls/genData/papalexi', 
#                   'desktop'='C:/Users/Catherine W/Documents/Research/genData/papalexi', 
#                   'ubergenno'='/raid6/Catherine/papalexi')
# # location of intermediate save files/plots are written
# save_dir = switch(DEVICE,
#                   'laptop'='C:/Users/Cathe/Documents/School/Genetic Hypothesis Testing using Negative Controls/DoubleBridge/saves/papalexi_saves', 
#                   'desktop'='C:/Users/Catherine W/Documents/Research/DoubleBridge/saves/papalexi_saves', 
#                   'ubergenno'='/raid6/Catherine/papalexi/papalexi_saves')
# # location of utils code
# util_dir = switch(DEVICE,
#                   'laptop'='C:/Users/Cathe/Documents/School/Genetic Hypothesis Testing using Negative Controls/DoubleBridge/code/utils', 
#                   'desktop'='C:/Users/Catherine W/Documents/Research/DoubleBridge/code/utils', 
#                   'ubergenno'='/home/catheri2/DoubleBridge/code/utils')

assertthat::assert_that(!is.null(data_dir), msg='first arg must be: laptop, desktop, or ubergenno')

assertthat::assert_that(length(args) > 1, msg="must give arg for specifying num imp genes 'Rscript <filename>.R ubergenno 4000'")
NUM_IMPORTANTGENES = as.integer(args[2]) # should be max importance from setting... TODO: change to remove this input. extract this value from saved setting

assertthat::assert_that(length(args) > 2, msg="must give arg for specifying chosen AYZW name 'Rscript <filename>.R ubergenno C'")
AYZW_setting_name = args[3]



# =================== Start ====================================================
print(sprintf("[%s] START: lmYA Pval", Sys.time()))

# # source CB utility functions (mainly CB fn for estimating)
# source(sprintf('%s/CB_utils.R', util_dir))
# Instead source new estimation methods' files
# source(sprintf('%s/OCB2SLS.R', util_dir))      # 2-Stage Least Squares
# source(sprintf('%s/OCBGMM.R', util_dir))       # Gen. Method of Moments
# source(sprintf('%s/OCBLinBridge.R', util_dir)) # Lin Bridge Plug-In and One-Step



source(sprintf('%s/CBEstAll.R', util_dir)) # for functions to estimate here



# load chosen AYZW names
AY   = read.csv(sprintf('%s/cbgenes/%s/AY.csv', save_dir, AYZW_setting_name))
AYZW = readRDS(sprintf('%s/cbgenes/%s/AYZW.rds', save_dir, AYZW_setting_name))


# load gene importance info
imp_gene_names = readRDS(sprintf('%s/important_genes_name.rds', save_dir))
imp_gene_idx   = readRDS(sprintf('%s/important_genes_idx.rds',  save_dir))
imp_gene = data.frame(gene     = imp_gene_names,
                      gene_idx = imp_gene_idx,
                      gene_imp_rank = 1:length(imp_gene_names))
# load in prev saved df (safer than constructing df again?)
# imp_gene = read.csv(sprintf('%s/gene_deviance.csv', save_dir)) |> 
#               arrange(desc(deviance)) |>
#               rename(gene=gene_name, gene_idx=idx) |>
#               select(gene, gene_idx, deviance) |> 
#               mutate(gene_imp_rank = row_number())


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
gene_norm = readin_gene_norm[1:NUM_IMPORTANTGENES, ] # dim = 4000 x 20729 = #important x #cells
rhdf5::h5closeAll()
invisible(gc(verbose=FALSE))



# =================== Define File Specific Functions ======================================
print(sprintf("[%s]    - Define File Specific Functions", Sys.time()))


# All the NT grna idx
NT_names = grna_odm |> ondisc::get_feature_covariates() |> 
  filter(target_type == 'non-targeting') |> rownames()
NT_idx = which(apply(X = grna_odm[[NT_names, ]], MARGIN = 2, FUN = sum) > 0)


get_ATE_est0 = get_ATE_est_make(AY=AY, AYZW=AYZW, gene_norm=gene_norm, grna_rownames=grna_rownames, 
                                NT_idx=NT_idx, imp_gene_names=imp_gene_names,
                                CB_setting=NULL, run_OCBLinOSEst=FALSE,
                                save_path = NULL,
                                est_lmYA_only=TRUE)

# lm might error if all cells got treatment/control e.g. grna_odm[['PDL1g3', ]] |> sum()
get_ATE_est <- function(AY_idx, ZW_idx) {
  res_df = tryCatch({get_ATE_est0(AY_idx=AY_idx, ZW_idx=ZW_idx)},
                    error = function(cond) {
                      # message(sprintf('Error est CondMomentOCBOS with %s', 
                      #                 gammaSetting)) 
                      return(NULL)
                    })
  # if errored, return NULL
  if(is.null(res_df)) {
    return(data.frame(ATE=NA, pval=NA))
  } else {
    return(res_df)
  }
}


# just AY idx changes (bc the same for all ZW_idx)
ATEargs = expand.grid(AY_idx = 1:nrow(AY), 
                      ZW_idx  = 1) 
# NUMROWS = 1000
# NUMROWS = 16 # 16 runs on 8 threads took 7 mins
# NUMROWS = 8 # 8 runs with fewer large dim took 1 min?
NUMROWS = nrow(ATEargs)
whichROWS = 1:NUMROWS

# =================== Get ATEs (parallel) ====================================
print(sprintf("[%s]    - Get ATEs (parallel, not actually)", Sys.time()))
t0 = Sys.time()
ATE_par = future.apply::future_mapply(get_ATE_est,
                                      AY_idx = ATEargs[whichROWS, 1], # ATEargs[1:NUMROWS, 1],
                                      ZW_idx = ATEargs[whichROWS, 2], # ATEargs[1:NUMROWS, 2],
                                      future.globals = TRUE,
                                      future.seed = 56789)
t1 = Sys.time()
print(sprintf("[%s]        - %2.2f", Sys.time(), (t1 - t0)))


ATE_df = NULL
for(i in 1:ncol(ATE_par)) {
  ATE_df = rbind(ATE_df, 
                 cbind(data.frame(AY_idx = ATEargs[whichROWS[i], 1], 
                                  ZW_idx = ATEargs[whichROWS[i], 2]),
                       ATE_par[,i] |> data.frame()))
  
  
}
ATE_df = cbind(AY[whichROWS, ], ATE_df)

write.csv(ATE_df, 
          sprintf('%s/cbgenes/%s/lmYA.csv', save_dir, AYZW_setting_name), row.names = FALSE)



# =================== END ======================================================
print(sprintf("[%s] END", Sys.time()))


# =================== Some Plots of the pvals/ATEs =============================
# print(sprintf("[%s]    - Some Plots", Sys.time()))
# ggplot(ATE_df,
#        aes(x = pval)) +
#   geom_histogram()
# 
# 
# ggplot(ATE_df,
#        aes(sample = pval)) +
#   geom_qq(distribution = stats::qunif) +
#   geom_abline(aes(slope = 1, intercept = 0)) +
#   facet_wrap(vars(type), ncol = 1)
# 
# 
# transform_pval <- function(x){-log(x, base = 10)}
# 
# unifpointsref = runif(n = 10000, min = 0, max = 1) |> transform_pval()
# qqplot(unifpointsref,
#        ATE_df |> filter(type == 'negative') |> pull(pval) |> transform_pval());
# qqline(unifpointsref, col = 'red')
# 
# qqplot(runif(n =  1000, min = 0, max = 1) |> transform_pval(),
#        runif(n =  nrow(ATE_df), min = 0, max = 1) |> transform_pval());
# qqline(runif(n =  1000, min = 0, max = 1) |> transform_pval())
# 
# 
# 
# ggplot(ATE_df) +
#   geom_histogram(aes(x = pval)) +
#   # geom_histogram(aes(x = pval |> transform_pval())) +
#   facet_wrap(vars(type), ncol = 1, scales = 'free_y')
# 
# 
# ggplot(ATE_df |> filter(type == 'negative')) +
#   geom_histogram(aes(x = pval |> transform_pval()))
# 
# 
# ggplot(ATE_df |> filter(type == 'positive')) +
#   geom_histogram(aes(x = pval))



