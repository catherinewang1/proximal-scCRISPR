# ------------------------------------------------------------------------------------------------ #
#     Estimate Effect (ATE) with count (original raw counts) outcome Y 
# Poisson GLM
# Poisson GLM           regression with 'un'-measured confounders
# Negative Binomial GLM
# Negative Binoimial GLM regression with 'un'-measured confounders
# e.g. nice -10 Rscript 3.2_papalexi_countGLM.R ubergenno A
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



suppressPackageStartupMessages(library(future.apply))
# options(future.globals.maxSize= 850*1024^2) #1st num is MB
options(future.globals.maxSize= 2500*1024^2) #1st num is MB
plan(multisession, workers = 8)
# plan(sequential)


theme_set(theme_cowplot() +
            theme(plot.title = element_text(hjust = .5),
                  plot.subtitle = element_text(hjust = .5)))


# === === === === === === === === === ===
save_intermediateATEs = 'yes' # 'yes'/'no' whether to save intermediate ATEs as they are estimated
proximal_setting_name = 'simpleCount' # <-- Always this for this GLM counts script (variable: proximal_setting_name~=NC_name in other script)
# === === === === === === === === === ===

assertthat::assert_that(length(args) > 0, msg="must give arg for specifying device eg 'Rscript <filename>.R ubergenno'")
DEVICE = args[1]
source('../PATHS.R') # load in data_dir and save_dir and CODE_DIR, depending on DEVICE value
assertthat::assert_that(!is.null(data_dir), msg='first arg must be: laptop, desktop, or ubergenno')

assertthat::assert_that(length(args) > 1, msg="must give arg for specifying chosen AYZW name 'Rscript <filename>.R ubergenno C'")
AYZW_setting_name = args[2]

# TODO: prob should change to just proximal_settings... (bc includes glm counts settings)
source(sprintf('%s/AY/proximal_continuous_settings.r', save_dir)) # loads in list: proximal_continuous_settings

PROXIMAL_SETTINGS = proximal_continuous_settings[[proximal_setting_name]] # the current proximal settings
PROXIMAL_SETTINGS$DEVICE = DEVICE

# save parameter settings
dir.create(sprintf('%s/AY/%s/%s', save_dir, AYZW_setting_name, proximal_setting_name), recursive = FALSE, showWarnings = FALSE)
capture.output(print(PROXIMAL_SETTINGS),
               file = sprintf('%s/AY/%s/%s/proximal_setting.txt',
                              save_dir, AYZW_setting_name, proximal_setting_name))
saveRDS(PROXIMAL_SETTINGS,
        sprintf('%s/AY/%s/%s/proximal_setting.rds',
                save_dir, AYZW_setting_name, proximal_setting_name))

# =================== Start ========================================================================
print(sprintf("[%s] START: Estimate Papalexi-2021 as counts", Sys.time()))
print(sprintf("[%s]        w/ script: %s", Sys.time(), paste0(args, collapse = ', ')))




print(sprintf("[%s]    - Load files and functions", Sys.time()))

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
# gene_importance = read.csv(sprintf('%s/gene_deviance_topnoTFonly.csv', save_dir)) |> 
#                   dplyr::select(gene_name, gene_idx, importance_rank)

gene_importance = read.csv(sprintf('%s/gene_deviance_gene_norm.csv', save_dir)) |> 
                  dplyr::select(gene_name, gene_idx, importance_rank, gene_norm_idx)

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
row.names(grna) = grna_rownames

# load measured covariates (=U, pretend unmeasured for proximal estimation)
cell_covariates = gene_odm |> ondisc::get_cell_covariates()
# change batch info, bc lane determines rep_1, so redundant info
# table(extra_covariates |> dplyr::select(lane, bio_rep))
library_size = cell_covariates$n_umis
cell_covariates = cell_covariates |> 
                      dplyr::mutate(lane_bio_rep = paste0(lane, '_', bio_rep)) |>
                      dplyr::select(-lane, -bio_rep) |>
                      dplyr::select(-n_nonzero, -n_umis)
# Does not include offset! library size for cell i = \sum_gene count_{i, gene}
# > head(cell_covariates, 1+1)
#                         phase     p_mito lane_bio_rep
#     l1_AAACCTGAGCCAGAAC    G1 0.02295577  Lane1_rep_1
#     l1_AAACCTGAGTGGACGT    G1 0.04512939  Lane1_rep_1



# # load normalized gene exp
# h5file      = paste0(save_dir, "/gene.h5"); print(h5file)
# reading_hd5file  = rhdf5::H5Fopen(name = h5file)
# readin_gene_norm = reading_hd5file&'gene_norm'
# gene_norm = readin_gene_norm[, 1:ncol(gene_odm)] # eg dim = 4000 x 20729 = #important x #cells
# rhdf5::h5closeAll()
# invisible(gc(verbose=FALSE))





# use Sparse PCA Loadings
# if(N_subsample == 'all') {  N_subsample = ncol(gene_odm) }
# NCs = readRDS(sprintf('%s/spca/NCloadings_sumabs=%.1f_K=%d_N=%d.rds', save_dir, my_sumabsv, my_K, N_subsample)) # save 
# Or the constructed clusters (averages)
# NCs = readRDS(sprintf('%s/spca/NCavg_sumabs=%.1f_K=%d_N=%d.rds', save_dir, my_sumabsv, my_K, N_subsample)) # save
# NCs = data.frame(NCs)
# colnames(NCs) = paste0('NC', 1:ncol(NCs))




# All the NT grna idx
NT_names = grna_odm |> ondisc::get_feature_covariates() |>
  filter(target_type == 'non-targeting') |> rownames()
NT_idx = which(apply(X = grna_odm[[NT_names, ]], MARGIN = 2, FUN = sum) > 0)







# =================== Define fns for Effects (parallel) ============================================
print(sprintf("[%s]    - Define fns for Effects (parallel)", Sys.time()))

# source(sprintf('%s/estimate_effects.R', util_dir)) # for functions to estimate here



estimate_effect_0 = estimate_effect_count_make(
	                              AY               = AY, 
                                which_estimators = PROXIMAL_SETTINGS$which_estimators,
                                gene_odm         = gene_odm, # <- maybe change to loaded in memory (for speed)
                                grna_odm         = grna_odm,
                                NT_idx           = NT_idx,
	                              library_size     = library_size,
                                U_confounders    = cell_covariates,
                                save_path        = switch(save_intermediateATEs,
                                                    'yes' = sprintf('%s/AY/%s/%s', save_dir, AYZW_setting_name, proximal_setting_name),
                                                    'no'  = NULL)
                                )


# test = estimate_effect_0(AY_idx = 3)



estimate_effect <- function(AY_idx) {
  res_df = tryCatch({estimate_effect_0(AY_idx=AY_idx)},
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

# =================== Get Effects (parallel) =======================================================
print(sprintf("[%s]    - Get Effects (parallel)", Sys.time()))
t0 = Sys.time()
ATE_par = future.apply::future_mapply(estimate_effect,
                                      AY_idx = ATEargs[whichROWS, 1], # ATEargs[1:NUMROWS, 1],
                                      future.globals = TRUE,
                                      future.seed = 123456, 
                                      SIMPLIFY = FALSE)
# # manually state globals
# future.globals = c('AY', 'AYZW', 'grna_rownames', 'grna', 
#                    'NT_idx', 'get_importance_rank', 'gene_norm', 
#                    'format_AYZW', 'get_CB_colnames', 'get_CB_est', 
#                    'CB', 'get_lmYA_est', 'get_ATE_est'))
t1 = Sys.time()
print(sprintf("[%s]        - %.2f mins", Sys.time(), difftime(t1, t0, units = "mins")))

saveRDS(ATE_par, file = sprintf('%s/AY/%s/%s/effects_countGLM_ATE_par.rds', save_dir, AYZW_setting_name, proximal_setting_name)) # save this parallel res as rds


# =================== Combine ATEs (into one df) ===================================================
print(sprintf("[%s]    - Combine ATEs", Sys.time()))
# ATE_par[, 1]
# ATE_par[[1, ]]
# length(ATE_par)
# unlist(ATE_par)
# ATE_par[1] |> as.data.frame()

ATE_df = NULL
for(whichROWS_idx in 1:length(whichROWS)) { # whichROWS_idx is always 1, 2, ...
  cur_ATE_par = ATE_par[[whichROWS_idx]]
  
  if(!is.null(cur_ATE_par)) {
    # AY test is whichROWS[whichROWS_idx]  = AY_idx!!! 
    AY_idx = whichROWS[whichROWS_idx]
    
    ATE_df = rbind(ATE_df,
                   cbind(data.frame(AY_idx=AY_idx),
                         AY[AY_idx, ] |> `rownames<-`( NULL ),
                         cur_ATE_par))
    rm(AY_idx)
  }
}

write.csv(x = ATE_df, file = sprintf('%s/AY/%s/%s/effects_countGLM.csv', save_dir, AYZW_setting_name, proximal_setting_name), row.names = FALSE)




# =================== END ==========================================================================
print(sprintf("[%s] END: %s", Sys.time(), AYZW_setting_name))



