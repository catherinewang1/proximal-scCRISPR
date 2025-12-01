# ------------------------------------------------------------------------------------------------ #
#     Estimate Effect (ATE) with continuous (normalized) outcome Y 
# linear regression
# linear regression with 'un'-measured confounders
# proximal (outcome) using pci2s package
# ------------------------------------------------------------------------------------------------ #
args = commandArgs(trailingOnly = TRUE)
# args = c('laptop', 'A1')
# args = c('laptop', 'A')




suppressPackageStartupMessages(library(assertthat)) # for some assert statements
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(ggplot2))    # plotting
suppressPackageStartupMessages(library(cowplot))
# also make sure you have pci2s package installed: remotes::install_github("KenLi93/pci2s")



library(future.apply)
# options(future.globals.maxSize= 850*1024^2) #1st num is MB
options(future.globals.maxSize= 1250*1024^2) #1st num is MB
plan(multisession, workers = 8)
# plan(sequential)


theme_set(theme_cowplot() +
            theme(plot.title = element_text(hjust = .5),
                  plot.subtitle = element_text(hjust = .5)))


# === Parameter Settings for Proximal Methods ===
num_NC_pairs = c(1, 3, 5, 10, 15, 20)
save_intermediateATEs = 'yes' # 'yes'/'no' whether to save intermedate ATEs as they are estimated

proximal_setting_name = 'simple'

# === Parameter Settings from SPCA
my_sumabsv = 34.5 # used for getting the name of the spca saved filename
my_K = 60
N_subsample = 5000 # subsample size, or 'all' if using all cells

# === Parameter Settings for which estimators to perform
which_estimators = list(lm_YA        = TRUE,
	                      lm_YAU       = TRUE,
	                      # pois_YAU     = TRUE,
	                      # nb_YAU       = TRUE,
                        # OCB_2SLS     = FALSE,
                        OCB_2SLS_pci2s=TRUE #,
                        # OCB_2SLSReg  = FALSE,
                        # OCB_GMM      = FALSE,
                        # OCB_GMMRw    = FALSE,
                        # OCB_GMMRwReg = FALSE,
                        # OCB_LinOSPI  = FALSE,
                        # OCB_LinOS    = FALSE,
                        # OCB_LinOStrim= FALSE
                        )
# === === === === === === === === === ===



assertthat::assert_that(length(args) > 0, msg="must give arg for specifying device eg 'Rscript <filename>.R ubergenno'")
DEVICE = args[1]
source('../PATHS.R') # load in data_dir and save_dir and CODE_DIR, depending on DEVICE value
assertthat::assert_that(!is.null(data_dir), msg='first arg must be: laptop, desktop, or ubergenno')


assertthat::assert_that(length(args) > 1, msg="must give arg for specifying chosen AYZW name 'Rscript <filename>.R ubergenno C'")
AYZW_setting_name = args[2]


# save parameter settings
CB_setting = list() 
CB_setting$num_NC_pairs = num_NC_pairs
CB_setting$proximal_setting_name = 'simple'

dir.create(sprintf('%s/AY/%s/%s', save_dir, AYZW_setting_name, proximal_setting_name), recursive = FALSE, showWarnings = FALSE)
capture.output(print(CB_setting),
               file = sprintf('%s/AY/%s/%s/proximal_setting.txt',
                              save_dir, AYZW_setting_name, proximal_setting_name))
saveRDS(CB_setting,
        sprintf('%s/AY/%s/%s/proximal_setting.rds',
                save_dir, AYZW_setting_name, proximal_setting_name))




# =================== Start ====================================================
print(sprintf("[%s] START: Estimate", Sys.time()))



source(sprintf('%s/estimate_effects.R', util_dir)) # for functions to estimate here


# source(sprintf('%s/CBEstAll.R', util_dir)) # for functions to estimate here
# source(sprintf('%s/CBEstAllSPCA.R', util_dir)) # add'l specifically for SPCA


# load chosen AYZW names
AY   = read.csv(sprintf('%s/AY/%s/AY.csv', save_dir, AYZW_setting_name))
# AYZW = readRDS(sprintf('%s/cbgenes/%s/AYZW.rds', save_dir, AYZW_setting_name))

# # Only for PCA: copy to pca folder too (bc plots script reads in AY, AYZW)
# write.csv(AY, sprintf('%s/spca/cbgenes/%s/AY.csv', save_dir, AYZW_setting_name))
# saveRDS(AYZW, sprintf('%s/spca/cbgenes/%s/AYZW.rds', save_dir, AYZW_setting_name))

# load gene importance info: gene_name, gene_idx (idx for raw), importance_rank (rank + idx for normalized)
gene_importance = read.csv(sprintf('%s/gene_deviance_topnoTFonly.csv', save_dir)) |> 
                  dplyr::select(gene_name, gene_idx, importance_rank)


# imp_gene_names = readRDS(sprintf('%s/important_genes_name.rds', save_dir))
# imp_gene_idx   = readRDS(sprintf('%s/important_genes_idx.rds',  save_dir))
# imp_gene = data.frame(gene     = imp_gene_names,
#                       gene_idx = imp_gene_idx,
#                       gene_imp_rank = 1:length(imp_gene_names))
# load in prev saved df (safer than constructing df again?)
# imp_gene = read.csv(sprintf('%s/gene_deviance.csv', save_dir)) |>
#               arrange(desc(deviance)) |>
#               rename(gene=gene_name, gene_idx=idx) |>
#               select(gene, gene_idx, deviance) |>
#               mutate(gene_imp_rank = row_number())

# Load choosing AY (ZW) settings
AYZW_setting = readRDS(sprintf('%s/AY/%s/AYZW_setting.rds', save_dir, AYZW_setting_name))
NUM_IMPORTANTGENES = AYZW_setting$MAX_Y_IMPORTANCE

# create gene and grna ondisc managers
gene_odm <- ondisc::read_odm(odm_fp      = paste0(data_dir, "/papalexi-2021/processed/gene/expression_matrix.odm"),
                             metadata_fp = paste0(data_dir, "/papalexi-2021/processed/gene/metadata.rds"))
grna_odm <- ondisc::read_odm(odm_fp      = paste0(data_dir, "/papalexi-2021/processed/grna_assignment/assignment_matrix.odm"),
                             metadata_fp = paste0(data_dir, "/papalexi-2021/processed/grna_assignment/metadata.rds"))

# load grna assignments (load all into memory)
grna = grna_odm[[,1:ncol(grna_odm)]] # |> as.matrix() # ~110 x 20729 = #grnas x #cells
grna_rownames = grna_odm |> ondisc::get_feature_covariates() |> rownames()


# load measured covariates (=U, pretend unmeasured for proximal estimation)
cell_covariates = gene_odm |> ondisc::get_cell_covariates()
# change batch info, bc lane determines rep_1, so redundant info
# table(extra_covariates |> dplyr::select(lane, bio_rep))
cell_covariates = cell_covariates |> 
                      dplyr::mutate(lane_bio_rep = paste0(lane, '_', bio_rep)) |>
                      dplyr::select(-lane, -bio_rep) |>
                      dplyr::select(-n_nonzero, -n_umis)

# load normalized gene exp
h5file      = paste0(save_dir, "/gene.h5"); print(h5file)
reading_hd5file  = rhdf5::H5Fopen(name = h5file)
readin_gene_norm = reading_hd5file&'gene_norm'
gene_norm = readin_gene_norm[1:NUM_IMPORTANTGENES, ] # eg dim = 4000 x 20729 = #important x #cells
rhdf5::h5closeAll()
invisible(gc(verbose=FALSE))





# use Sparse PCA Loadings
if(N_subsample == 'all') {  N_subsample = ncol(gene_odm) }
NCs = readRDS(sprintf('%s/spca/NCloadings_sumabs=%.1f_K=%d_N=%d.rds', save_dir, my_sumabsv, my_K, N_subsample)) # save 
# Or the constructed clusters (averages)
# NCs = readRDS(sprintf('%s/spca/NCavg_sumabs=%.1f_K=%d_N=%d.rds', save_dir, my_sumabsv, my_K, N_subsample)) # save
# NCs = data.frame(NCs)
# colnames(NCs) = paste0('NC', 1:ncol(NCs))




# All the NT grna idx
NT_names = grna_odm |> ondisc::get_feature_covariates() |>
  filter(target_type == 'non-targeting') |> rownames()
NT_idx = which(apply(X = grna_odm[[NT_names, ]], MARGIN = 2, FUN = sum) > 0)










# =============================================================================


# source(sprintf('%s/estimate_effects.R', util_dir)) # for functions to estimate here



estimate_ATE_0 = estimate_ATE_make(AY                      = AY, 
                                   gene_norm               = gene_norm, 
                                   NCs                     = NCs, 
                                   grna_rownames           = grna_rownames, 
                                   grna                    = grna, 
                                   NT_idx                  = NT_idx, 
                                   gene_importance         = gene_importance,
                                   which_estimators        = which_estimators, 
                                   save_path = switch(save_intermediateATEs,
                                                     'yes' = sprintf('%s/AY/%s/%s', save_dir, AYZW_setting_name, proximal_setting_name),
                                                     'no'  = NULL),
                                   U_confounders           = cell_covariates)
# test = estimate_ATE_0(AY_idx = 3)



estimate_ATE <- function(AY_idx) {
  res_df = tryCatch({estimate_ATE_0(AY_idx=AY_idx)},
                    error = function(cond) {
                      # message(sprintf('Error est CondMomentOCBOS with %s', 
                      #                 gammaSetting)) 
                      return(NULL)
                    })
  # if errored, return NULL
  if(is.null(res_df)) {
    return(NULL)
  } else {
    return(res_df)
  }
}


# NUM_NCENCO_per_AY = length(AYZW[[1]][[1]]) # prev defined, will be length of this sublist
# ATEargs = expand.grid(AY_idx = 1:nrow(AY), 
#                       ZW_idx  = 1:NUM_NCENCO_per_AY) |> 
#   arrange(AY_idx, ZW_idx) # rearrange, AY pair first

ATEargs = data.frame(AY_idx = 1:nrow(AY))
# lessen the amount...
# ATEargs = expand.grid(AY_idx = c(1:30, 120:134), ZW_idx  = 1:NUM_NCENCO_per_AY)



# NUMROWS = 10
NUMROWS = nrow(ATEargs)
# whichROWS = 300:nrow(ATEargs)
# whichROWS = 1:1164
whichROWS = 1:NUMROWS
# whichROWS = 1165:NUMROWS

# =================== Get ATEs (parallel) ====================================
print(sprintf("[%s]    - Get ATEs (parallel)", Sys.time()))
t0 = Sys.time()
ATE_par = future.apply::future_mapply(estimate_ATE,
                                      AY_idx = ATEargs[whichROWS, 1], # ATEargs[1:NUMROWS, 1],
                                      future.globals = TRUE,
                                      future.seed = 123456)
# # manually state globals
# future.globals = c('AY', 'AYZW', 'grna_rownames', 'grna', 
#                    'NT_idx', 'get_importance_rank', 'gene_norm', 
#                    'format_AYZW', 'get_CB_colnames', 'get_CB_est', 
#                    'CB', 'get_lmYA_est', 'get_ATE_est'))
t1 = Sys.time()
print(sprintf("[%s]        - %2.2f", Sys.time(), (t1 - t0)))

# ATE_par[, 1]
# ATE_par[[1, ]]
# length(ATE_par)
# unlist(ATE_par)
# ATE_par[1] |> as.data.frame()


ATE_df = NULL
for(AY_idx in whichROWS) {
  if(!is.null(ATE_par[[AY_idx]])) {
    ATE_df = rbind(ATE_df,
                   cbind(data.frame(AY_idx=AY_idx),
                         AY[AY_idx, ] |> `rownames<-`( NULL ),
                         ATE_par[[AY_idx]]))
  }
}

write.csv(x = ATE_df, file = sprintf('%s/AY/%s/%s/effects_continuous.csv', save_dir, AYZW_setting_name, proximal_setting_name), row.names = FALSE)






