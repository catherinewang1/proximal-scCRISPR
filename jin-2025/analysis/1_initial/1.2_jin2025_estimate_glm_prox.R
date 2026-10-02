# ------------------------------------------------------------------------------------------------ #
#     Initial Estimates on Jin-2025 dataset
#     Mouse brain data from Xin Jin's lab ~ 2025 and downloaded from dropbox 
# Email: Ergan to Jin-Hong to Catherine:
#     Please check the following two links of the data available. 
#     Those data is shared in 2 batches by Seoyeon, where one is containing 10 perturbations and another one contains 48 perturbations.
#     
# 
# in 1.0 notebook
# - Pre-processing (gene importance and normalization)
# - Perform PCA 
# in 1.1 script
# - initial lin reg and neg bin glm estimates
# - some estimates
# ------------------------------------------------------------------------------------------------ #
args = commandArgs(trailingOnly = TRUE)
# args = c('macbook')


assertthat::assert_that(length(args) > 0, msg="must give arg for specifying device eg 'Rscript <filename>.R ubergenno'")
DEVICE = args[1]
source('../PATHS.R') # load in data_dir and save_dir and CODE_DIR, depending on DEVICE value
save_dir = paste0(save_dir, '/1_initial')
assertthat::assert_that(!is.null(data_dir), msg='first arg must be: laptop, desktop, or ubergenno')



library(dplyr)
# library(Seurat)
# library(SeuratObject)
# library(anndataR)
library(ggplot2)
library(h5lite)
library(DelayedMatrixStats)
library(sparseMatrixStats)

library(future.apply)
options(future.globals.maxSize= 2500*1024^2) #1st num is MB
# plan(multisession, workers = 24)
plan(multisession, workers = 8)
# plan(sequential)

AY_SETTING_NAME = 'small'

GENE_NUMNONZEROCELLS_MIN = 7500
GENE_IMPORTANCERANK_MAX = 2000

AY_NUM_NEGATIVE  = 150
AY_NUM_DISCOVERY = 150

set.seed(12345)

# =================== START ========================================================================


source(sprintf('%s/quick_estimate/quickGLMs.R', util_dir))
source(sprintf('%s/quick_estimate/quickPCI2S.R', util_dir))


dir.create(sprintf('%s/AY/', save_dir), showWarnings = FALSE)
dir.create(sprintf('%s/AY/%s/', save_dir, AY_SETTING_NAME), showWarnings = FALSE)
dir.create(sprintf('%s/AY/%s/initiallinreg/', save_dir, AY_SETTING_NAME), showWarnings = FALSE)
dir.create(sprintf('%s/AY/%s/intermediateATEs/', save_dir, AY_SETTING_NAME), showWarnings = FALSE)

# ============================================================================== #
# =================== Load data from h5ad   ====================================
# ============================================================================== #



hdf5_filename = sprintf('%s/jin2025_subsetgenes.hdf5', save_dir)
h5_filename   = sprintf('%s/jin2025.h5', save_dir)

# HDF5Array::HDF5Array()
# gene_norm = rhdf5::h5read(hdf5_filename, name = 'gene_norm/gene_norm')
# gene_norm = h5lite::h5_read(hdf5_filename, name = 'gene_norm/gene_norm')
# gene_counts = h5lite::h5_read(h5_filename, "counts/counts")

# rhdf5::h5ls(hdf5_filename)
# gene_norm  = rhdf5::h5read(hdf5_filename, name = 'gene_norm/gene_norm')
# gene_counts = rhdf5::h5read(hdf5_filename, name = 'counts/counts')
# row.names(gene_counts) = rhdf5::h5read(hdf5_filename, name = 'counts/col_genenames')
# colnames(gene_counts) = rhdf5::h5read(hdf5_filename, name = 'counts/row_cellnames')
# 
# 
# row.names(gene_norm) = rhdf5::h5read(hdf5_filename, name = 'gene_norm/.gene_norm_dimnames/1')
# colnames(gene_norm) = rhdf5::h5read(hdf5_filename, name = 'gene_norm/.gene_norm_dimnames/2')
# 




gene_norm = HDF5Array::HDF5Array(hdf5_filename, name='gene_norm/gene_norm')
row.names(gene_norm) = rhdf5::h5read(hdf5_filename, name = 'gene_norm/.gene_norm_dimnames/1') |> as.vector()
colnames( gene_norm) = rhdf5::h5read(hdf5_filename, name = 'gene_norm/.gene_norm_dimnames/2') |> as.vector()


gene_counts = HDF5Array::HDF5Array(hdf5_filename, name='gene_counts/gene_counts')
row.names(gene_counts) = rhdf5::h5read(hdf5_filename, name = 'gene_counts/.gene_counts_dimnames/1') |> as.vector()
colnames( gene_counts) = rhdf5::h5read(hdf5_filename, name = 'gene_counts/.gene_counts_dimnames/2') |> as.vector()


h5lite::h5_ls(h5_filename)

# grna = rhdf5::h5read(h5_filename, name = 'grna/grna', native = TRUE) 
# grna = apply(grna, MARGIN=c(1, 2), FUN=as.integer)
# row.names(grna) = rhdf5::h5read(h5_filename, name = 'grna/row_grnanames')
# colnames(grna) = rhdf5::h5read(h5_filename, name = 'grna/col_cellnames')


rhdf5::h5closeAll()
grna     = h5lite::h5_read(h5_filename, name = 'grna/grna')
grnaNT   = h5lite::h5_read(h5_filename, name = 'grnaNT/grnaNT')
metadata = h5lite::h5_read(h5_filename, name = 'cell_metadata/cell_metadata')

rhdf5::h5closeAll()


gene_dev_df = read.csv(file = sprintf('%s/gene_dev.csv', save_dir)) 


# ================ checks
assertthat::assert_that(all(row.names(gene_norm) == row.names(gene_counts))) # same genes and in same order
assertthat::assert_that(all(colnames(gene_norm) == colnames(gene_counts))) # same genes and in same order

# setdiff(row.names(gene_counts), row.names(gene_norm))
# 
# setdiff(row.names(gene_counts), normalized_genes_names)
# setdiff(row.names(gene_norm), normalized_genes_names)
# setdiff(row.names(gene_norm), row.names(gene_counts))

# ============================================================================== #
# =============  Construct AY tests: ===========================================
# ============================================================================== #

# ================ some prep for constructing tests

genes_nonzero = row.names(gene_counts)[(rowSums(gene_counts == 0) >= GENE_NUMNONZEROCELLS_MIN)] # genes that meet the nonzero threshold
all_Ys_df = gene_dev_df |>  
  filter(gene_name %in% genes_nonzero) |>  
  filter(is_normalized_gene) |> 
  # filter(!is_grna_target & !is_transcription_factor) |>
  arrange(importance_rank)

validsinglegenes = all_Ys_df |> 
  filter(!is_grna_target & !is_transcription_factor) |> # no grna target and no TF
  filter(importance_rank <= GENE_IMPORTANCERANK_MAX) |> pull(gene_name) # valid names for single genes


all_Ys = all_Ys_df |>filter(importance_rank <= GENE_IMPORTANCERANK_MAX | is_grna_target) |> pull(gene_name)


# =============  \__Construct Negative tests: ===============================================

# For now, just do 2 NT against others to try to make the sample sizes similar?
# ggplot(NULL, aes(x = metadata |> filter(gene_target == 'Non_target') |> pull(guide_call) )) +
#   geom_bar()
# metadata |> filter(gene_target == 'Non_target') |> pull(guide_call) |> hist()

NT_names = grnaNT |> row.names() 
NT_pairs = replicate(AY_NUM_NEGATIVE*1.5, sample(NT_names, 2, replace = FALSE)) |> t() |> data.frame() |> distinct()

NT_treatment_name = paste0(NT_pairs[, 1], '&', NT_pairs[, 2])[1:AY_NUM_NEGATIVE] # would err if the number specified is too large 

NT_outcome_name = sample(all_Ys, size = length(NT_treatment_name))

AY_negative = data.frame(grna = NT_treatment_name, 
                         gene = NT_outcome_name, 
                         type = 'negative')

rm(NT_names, NT_pairs, NT_treatment_name, NT_outcome_name)

# ============= \__Construct Discovery tests: ============================================

# dim(gene_count)
# gene_numnonzerocells = rowSums(gene_count); hist(gene_numnonzerocells)
AY_discovery = data.frame(grna = sample(setdiff(row.names(grna), 'Non_target'), size = AY_NUM_DISCOVERY, replace=TRUE),
                          gene = sample(all_Ys,                                 size = AY_NUM_DISCOVERY, replace=TRUE),
                          type = 'discovery')

# ============= \__Construct Positive tests: ============================================
pos_grna = row.names(grna)[row.names(grna) %in% row.names(gene_norm)]
AY_positive = data.frame(grna = pos_grna, gene = pos_grna,type='positive')

rm(pos_grna)


# ============= \__Combine AY tests: ============================================
AY = rbind( AY_positive, 
            AY_negative, 
           AY_discovery)

write.csv(AY, sprintf('%s/AY/%s/AY.csv', save_dir, AY_SETTING_NAME), row.names=FALSE)


rm(AY_negative, AY_discovery, AY_positive)


# ============================================================================== #
# ============= Realize Relevant Gene info in memory: ==========================
# ============================================================================== #
# required_Y = unique(AY$gene) # only these response genes' expression are required 
# 
# gene_counts_subset = gene_counts[required_Y, ] |> DelayedArray::realize()
# gene_norm_subset   =   gene_norm[required_Y, ] |> DelayedArray::realize()
# 
# 
# rm(gene_counts, gene_norm)



# ============================================================================== #
# ============= Initial LinReg and NegBin: =====================================
# ============================================================================== #

# if('initial_glms.csv' %in% list.files(sprintf('%s/AY/%s/initiallinreg/', save_dir, AY_SETTING_NAME))) {
#   print(sprintf("[%s]        - Loading glm (one fit per gene) on all pairs", Sys.time()))
#   initial_glms = read.csv(sprintf('%s/AY/%s/initiallinreg/all_glm_one.csv', save_dir, subfoldername))
# } else {
#   print(sprintf("[%s]        - Fitting glm (one fit per gene) on all pairs", Sys.time()))
#   t0 = Sys.time()
#   initial_glms_par = future.apply::future_mapply(fit_glms,
#                                             cur_gene = many_Ys, 
#                                             future.globals = TRUE,
#                                             future.seed = 123456, 
#                                             SIMPLIFY = FALSE)
#   t1 = Sys.time()
#   print(sprintf("[%s]        - %2.2f mins", Sys.time(), difftime(t1, t0, units = "mins")))
#   
#   initial_glms = NULL
#   for(idx in 1:length(initial_glms_par)) {   initial_glms = rbind(initial_glms, initial_glms_par[[idx]])   }
#   write.csv(initial_glms, sprintf('%s/AY/%s/initiallinreg/initial_glms.csv', save_dir, subfoldername))
#   rm(initial_glms_par)
# }
# 
# 
# 
# (AY_idx, AY_df, control_cell_idx, gene_norm, gene, grna, libsize_log, covariates=NULL)





# ============================================================================== #
# =================== Get ATEs (parallel) ======================================
# ============================================================================== #

numNCs = c(1,3, 5, 10)

control_cell_idx = which(grna['Non_target', ] == 1) |> as.vector() # just the overall control cell idx. use setdiff for negative tests

combined_grna = rbind(grna, grnaNT)

libsize_log = log(metadata[, 'num_rna_umi'] + 1e-5)
# libsize_log = rep(NA, ncol(grna))

NCs_pca = readRDS(sprintf('%s/pca/NCloadings.rds', save_dir)) |> as.matrix() # PCA loadings


# some setup
estimate_effects_errorhandling <- function(AY_idx) {
  res_df = tryCatch({estimate_effects(AY_idx=AY_idx, 
                                      # test information
                                      AY=AY, 
                                      # proximal number of NCs
                                      numNCs=numNCs, 
                                      # data
                                      control_cell_idx=control_cell_idx, 
                                      gene_norm=gene_norm, 
                                      gene=gene_counts, 
                                      grna=combined_grna,
                                      libsize_log=libsize_log, 
                                      # info for NCs or covariates
                                      NCs_pca=NCs_pca,
                                      validsinglegenes=validsinglegenes,
                                      covariates=NULL,
                                      # saving intermediateATEs
                                      save_intermediateATEs=TRUE, 
                                      skip_if_saved_intermediateATEs=TRUE,
                                      intermediateATEs_folder = sprintf('%s/AY/%s/intermediateATEs/', save_dir, AY_SETTING_NAME),
                                      # which proximal estimates
                                      run_proximal_linear_pca         = TRUE,
                                      run_proximal_negbin_pca         = FALSE,
                                      run_proximal_linear_singlegenes = TRUE,
                                      run_proximal_negbin_singlegenes = FALSE,
                                      subsample_control_size = 5000)
                    },
                    error = function(cond) {
                      # message(sprintf('Error est ATE w pci2s with %s', 
                      #                 AY_idx)) 
                      return(NULL)
                    })
  # if errored, return NULL
  if(is.null(res_df)) {
    return(NULL)
  } else {
    return(res_df)
  }
}




# # test a run
# a = estimate_effects(AY_idx=4,
#                  # test information
#                  AY=AY,
#                  # proximal number of NCs
#                  numNCs=c(1, 5),
#                  # data
#                  control_cell_idx=control_cell_idx,
#                  gene_norm=gene_norm,
#                  gene=gene_counts,
#                  # gene_norm=gene_norm_subset,
#                  # gene=gene_counts_subset,
#                  grna=combined_grna,
#                  libsize_log=libsize_log,
#                  # info for NCs or covariates
#                  NCs_pca=NCs_pca,
#                  validsinglegenes=validsinglegenes,
#                  covariates=NULL,
#                  # saving intermediateATEs
#                  save_intermediateATEs=FALSE,
#                  skip_if_saved_intermediateATEs=TRUE,
#                  intermediateATEs_folder = sprintf('%s/AY/%s/intermediateATEs/', save_dir, AY_SETTING_NAME),
#                  # which proximal estimates
#                  run_proximal_linear_pca         = TRUE,
#                  run_proximal_negbin_pca         = FALSE,
#                  run_proximal_linear_singlegenes = TRUE,
#                  run_proximal_negbin_singlegenes = FALSE)


# estimate_effects_errorhandling(4)
# 
# a = combined_grna[c("Non_target_399", "Non_target_313"), ]





print(sprintf("[%s]    - Get Estimates (parallel)", Sys.time()))

whichROWS = 1:60
# whichROWS = 1:nrow(AY)
print(sprintf("[%s]        + nrow(AY)=%s", Sys.time(), nrow(AY)))
print(sprintf("[%s]        + #unique grna=%s", Sys.time(), length(unique(AY$grna))))
print(sprintf("[%s]        + #unique gene=%s", Sys.time(), length(unique(AY$gene))))

t0 = Sys.time()
effects_par = future.apply::future_mapply(estimate_effects_errorhandling,
                                          AY_idx = whichROWS, 
                                          future.globals = TRUE,
                                          future.seed = 123456, 
                                          SIMPLIFY = FALSE)
t1 = Sys.time()
print(sprintf("[%s]        - %2.2f mins", Sys.time(), difftime(t1, t0, units = "mins")))

saveRDS(effects_par, file = sprintf('%s/AY/%s/effects_par.rds', save_dir, AY_SETTING_NAME)) # save this parallel res as rds



# =================== \__Combine ATEs (into one df) ===================================================
print(sprintf("[%s]    - Combine ATEs", Sys.time()))

effects_df = NULL
for(whichROWS_idx in 1:length(whichROWS)) { # whichROWS_idx is always 1, 2, ...
  cur_effects_par = effects_par[[whichROWS_idx]]
  if(!is.null(cur_effects_par)) {
    # AY test is whichROWS[whichROWS_idx]  = AY_idx!!! 
    effects_df = dplyr::bind_rows(effects_df, cur_effects_par)
  }
  rm(cur_effects_par)
}
write.csv(x = effects_df, file = sprintf('%s/AY/%s/effects.csv', save_dir, AY_SETTING_NAME), row.names = FALSE)















