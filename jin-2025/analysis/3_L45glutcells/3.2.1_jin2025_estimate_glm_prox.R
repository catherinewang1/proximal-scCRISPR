# ------------------------------------------------------------------------------------------------ #
#     Initial Estimates on Jin-2025 dataset on glut cell type
#     Mouse brain data from Xin Jin's lab ~ 2025 and downloaded from dropbox 
# Email: Ergan to Jin-Hong to Catherine:
#     Please check the following two links of the data available. 
#     Those data is shared in 2 batches by Seoyeon, where one is containing 10 perturbations and another one contains 48 perturbations.
#     
# 
# 2.1 preprocessing
# - Pre-processing (gene importance and normalization)
# - Perform PCA 
# 2.2 estimates
# - 2.2.1 glms + proximal
# - 2.2.2 deseq2
# - 2.2.3 causarray
# - 2.2.4 sceptre
# ------------------------------------------------------------------------------------------------ #
args = commandArgs(trailingOnly = TRUE)
# args = c('macbook', 'smallDense', '8')


assertthat::assert_that(length(args) > 0, msg="must give arg for specifying device eg 'Rscript <filename>.R ubergenno'")
DEVICE = args[1]
source('../PATHS.R') # load in data_dir and save_dir and CODE_DIR, depending on DEVICE value
assertthat::assert_that(!is.null(data_dir), msg='first arg must be: laptop, desktop, or ubergenno')

assertthat::assert_that(length(args) > 0, msg="must give arg for specifying AY setting eg 'Rscript <filename>.R ubergenno small'")
AY_SETTING_NAME = args[2]

assertthat::assert_that(length(args) > 0, msg="must give arg for specifying num threads eg 'Rscript <filename>.R ubergenno small 8'")
THREADS = as.integer(args[3])




# suppressPackageStartupMessages(library(Seurat))
# suppressPackageStartupMessages(library(SeuratObject))
# suppressPackageStartupMessages(library(anndataR))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
# suppressPackageStartupMessages(library(h5lite))
                                                                                                                                   
suppressPackageStartupMessages(library(DelayedMatrixStats))
suppressPackageStartupMessages(library(sparseMatrixStats))

suppressPackageStartupMessages(library(future.apply))
options(future.globals.maxSize= 5000*1024^2) #1st num is MB
# plan(multisession, workers = 24)
plan(multisession, workers = THREADS)
# plan(sequential)



PARAMS = list( 
  'TAG'       = 'L45Glut',
  'subfolder' = '3_L45glutcells',
  'seed'      = 12345,
  'A_MIN_TRTCELLS'  = 40, # filter for grnas with at least this number of cells (keep low here, set higher later))
  'Y_MAX_GENE_IMPORTANCE' = 1500
)
save_dir = paste0(save_dir, '/', PARAMS[['subfolder']])
set.seed(PARAMS[['seed']])


# AY_SETTING_NAME = 'denseSmall' # key of one of the following settings, prob change this to an arg parameter

AY_SETTINGS = list(
  'small' = list(
    AY_NUM_NEGATIVE  = 150,
    AY_NUM_DISCOVERY = 150,
    numNCs = c(1, 3, 5, 10),
    prox_methods = list(run_proximal_linear_pca         = TRUE,
                        run_proximal_negbin_pca         = FALSE,
                        run_proximal_linear_singlegenes = TRUE,
                        run_proximal_negbin_singlegenes = FALSE),
    skip_if_saved_intermediateATEs=TRUE
  ),
  
  'smallNB' = list(
    AY_NUM_NEGATIVE  = 150,
    AY_NUM_DISCOVERY = 150,
    numNCs = c(1, 3, 5, 10),
    prox_methods = list(run_proximal_linear_pca         = TRUE,
                        run_proximal_negbin_pca         = TRUE,
                        run_proximal_linear_singlegenes = TRUE,
                        run_proximal_negbin_singlegenes = TRUE),
    skip_if_saved_intermediateATEs=TRUE
  ),
  
  'medium' = list(
    AY_NUM_NEGATIVE  = 1000,
    AY_NUM_DISCOVERY = 1000,
    numNCs = c(1, 3, 5, 8, 10, 15),
    prox_methods = list(run_proximal_linear_pca         = TRUE,
                        run_proximal_negbin_pca         = FALSE,
                        run_proximal_linear_singlegenes = TRUE,
                        run_proximal_negbin_singlegenes = FALSE),
    skip_if_saved_intermediateATEs=TRUE
  ),
  
  'denseSmall' = list(
    DENSE = TRUE,
    NUM_Y = 100,
    AY_NUM_NEGATIVE  = NA,
    AY_NUM_DISCOVERY = NA,
    numNCs = c(1, 3, 5, 8, 10, 15),
    prox_methods = list(run_proximal_linear_pca         = TRUE,
                        run_proximal_negbin_pca         = FALSE,
                        run_proximal_linear_singlegenes = TRUE,
                        run_proximal_negbin_singlegenes = FALSE),
    skip_if_saved_intermediateATEs=FALSE
  ), 
  
  'denseMedium' = list(
    DENSE = TRUE,
    NUM_Y = 1000,
    AY_NUM_NEGATIVE  = NA,
    AY_NUM_DISCOVERY = NA,
    numNCs = c(1, 3, 5, 10),
    prox_methods = list(run_proximal_linear_pca         = TRUE,
                        run_proximal_negbin_pca         = FALSE,
                        run_proximal_linear_singlegenes = TRUE,
                        run_proximal_negbin_singlegenes = FALSE),
    skip_if_saved_intermediateATEs=FALSE
  )
  
)


SETTING = AY_SETTINGS[[AY_SETTING_NAME]]



# =================== START ========================================================================
print(sprintf("[%s] START: glm and proximal estimates for Jin2025", Sys.time()))

source(sprintf('%s/quick_estimate/quickGLMs.R', util_dir))
source(sprintf('%s/quick_estimate/quickPCI2S.R', util_dir))


dir.create(sprintf('%s/AY/', save_dir), showWarnings = FALSE)
dir.create(sprintf('%s/AY/%s/', save_dir, AY_SETTING_NAME), showWarnings = FALSE)
dir.create(sprintf('%s/AY/%s/initiallinreg/', save_dir, AY_SETTING_NAME), showWarnings = FALSE)
dir.create(sprintf('%s/AY/%s/intermediateATEs/', save_dir, AY_SETTING_NAME), showWarnings = FALSE)

# ============================================================================== #
# =================== Load data from h5 ========================================
# ============================================================================== #
print(sprintf("[%s]    - load h5 (prev saved from processData)", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
h5file   = sprintf("%s/jin2025_%s.h5", save_dir, PARAMS[['TAG']])


# =================== Load gene_norm
gene_norm = HDF5Array::HDF5Array(h5file, name='gene_norm/gene_norm')
row.names(gene_norm) = rhdf5::h5read(h5file, name = 'gene_norm/.gene_norm_dimnames/1') |> as.vector()
colnames( gene_norm) = rhdf5::h5read(h5file, name = 'gene_norm/.gene_norm_dimnames/2') |> as.vector()

# =================== Load gene_counts
gene_counts = HDF5Array::HDF5Array(h5file, name='gene_counts/gene_counts')
row.names(gene_counts) = rhdf5::h5read(h5file, name = 'gene_counts/.gene_counts_dimnames/1') |> as.vector()
colnames( gene_counts) = rhdf5::h5read(h5file, name = 'gene_counts/.gene_counts_dimnames/2') |> as.vector()


# h5lite::h5_ls(h5file)
# rhdf5::h5closeAll()
# rhdf5::h5ls(h5file)

# =================== Load grna, grnaNT, metadata, gene_dev_df

# using h5lite package, it is slightly easier, but the package is only available for newer versions of R
# grna     = h5lite::h5_read(h5file, name = 'grna/grna') # grna x cells
# metadata = h5lite::h5_read(h5file, name = 'cell_metadata/cell_metadata')

# use rhdf5 bc h5lite not available for R version 4.1 on server
grna            = rhdf5::h5read(h5file, name = 'grna/grna') |> t()
class(grna)     = "logical"
row.names(grna) = rhdf5::h5read(h5file, name = 'grna/grna_rownames') 
colnames( grna) = rhdf5::h5read(h5file, name = 'grna/grna_colnames') 

metadata            = rhdf5::h5read(h5file, name = 'cell_metadata/cell_metadata')
row.names(metadata) = rhdf5::h5read(h5file, name = 'cell_metadata/cell_metadata_rownames')


gene_dev_df = read.csv(file = sprintf('%s/gene_dev.csv', save_dir)) 


rhdf5::h5closeAll()

# ================ checks
assertthat::assert_that(all(row.names(gene_norm) == row.names(gene_counts))) # same genes and in same order
assertthat::assert_that(all(colnames(gene_norm) == colnames(gene_counts))) # same genes and in same order

# setdiff(row.names(gene_counts), row.names(gene_norm))
# 
# setdiff(row.names(gene_counts), normalized_genes_names)
# setdiff(row.names(gene_norm), normalized_genes_names)
# setdiff(row.names(gene_norm), row.names(gene_counts))




# check cellnames 
assertthat::assert_that(all(colnames(gene_norm) == colnames(gene_counts))) # same cells and in same order
assertthat::assert_that(all(colnames(gene_norm) == colnames(grna))) # same cells and in same order
# assertthat::assert_that(all(colnames(gene_norm) == colnames(grnaNT))) # same cells and in same order
assertthat::assert_that(all(colnames(gene_norm) == row.names(metadata))) # same cells and in same order


# ============================================================================== #
# =============  Construct AY tests: ===========================================
# ============================================================================== #

# ================ some prep for constructing tests

# genes_nonzero = row.names(gene_counts)[(rowSums(gene_counts == 0) >= GENE_NUMNONZEROCELLS_MIN)] # genes that meet the nonzero threshold, should be lready filtered
all_Ys_df = gene_dev_df |>  
  # filter(gene_name %in% genes_nonzero) |>  
  filter(is_normalized_gene) |> 
  # filter(!is_grna_target & !is_transcription_factor) |>
  arrange(importance_rank)

validsinglegenes = all_Ys_df |> 
  filter(!is_grna_target & !is_transcription_factor) |> # no grna target and no TF
  filter(importance_rank <= PARAMS$Y_MAX_GENE_IMPORTANCE) |> pull(gene_name) # valid names for single genes


all_Ys = all_Ys_df |>filter(importance_rank <= PARAMS$Y_MAX_GENE_IMPORTANCE | is_grna_target) |> pull(gene_name)


# all valid treatments (not all grnas, eg if there weren't enough cells)
# all_As = row.names(grna)[rowSums(grna) >= PARAMS$A_MIN_TRTCELLS] |> setdiff('Non_target') |> sort()

true_As = row.names(grna)[rowSums(grna) >= PARAMS$A_MIN_TRTCELLS]
true_As = true_As[!grepl('Non_target', true_As)]

fake_As = grep('Non_target', setdiff(row.names(grna), 'Non_target'), value=TRUE)


# construct AY tests base on whether dense setting is chosen or not
if(!is.null(SETTING$DENSE) && SETTING$DENSE) {
  # first select all the grna targets
  # then choose other Ys from all_Ys. all_Ys should be arranged by importance rank, so choose the top important/most deviance genes as responses
  selected_Ys = c(true_As, setdiff(all_Ys, true_As))[1:SETTING$NUM_Y]
  
  
  AY = expand.grid(grna = c(true_As, fake_As), gene = selected_Ys, stringsAsFactors = FALSE) |> mutate(type = 'discovery')
  AY[(AY$grna == AY$gene), 'type'] = 'positive'
  AY[grepl(pattern='Non_target', x = AY$grna), 'type'] = 'negative'
} else {
  # =============  \__Construct Negative tests: ===============================================
  AY_negative = data.frame(# grna = sample(setdiff(row.names(grna), 'Non_target'), size = AY_NUM_DISCOVERY, replace=TRUE),
    grna = sample(fake_As, size=SETTING$AY_NUM_NEGATIVE*1.2, replace=TRUE),
    gene = sample( all_Ys, size=SETTING$AY_NUM_NEGATIVE*1.2, replace=TRUE),
    type = 'negative') |> 
    filter(grna != gene) |>        # random choice of grnas and genes may result in pos tests, so sample more
    distinct() |> 
    slice_sample(n=SETTING$AY_NUM_NEGATIVE) # than needed, then remove and take the requested number of tests
  
  
  # ============= \__Construct Discovery tests: ============================================
  AY_discovery = data.frame(grna = sample(true_As, size=SETTING$AY_NUM_DISCOVERY*1.2, replace=TRUE),
                            gene = sample( all_Ys, size=SETTING$AY_NUM_DISCOVERY*1.2, replace=TRUE),
                            type = 'discovery') |> 
    filter(grna != gene) |>        # random choice of grnas and genes may result in pos tests, so sample more
    distinct() |> 
    slice_sample(n=SETTING$AY_NUM_DISCOVERY) # than needed, then remove and take the requested number of tests
  
  # ============= \__Construct Positive tests: ============================================
  pos_grna = row.names(grna)[row.names(grna) %in% row.names(gene_norm)]
  AY_positive = data.frame(grna = pos_grna, gene = pos_grna,type='positive')
  
  rm(pos_grna)
  
  
  # ============= \__Combine AY tests: ============================================
  AY = rbind( AY_positive, 
              AY_negative, 
              AY_discovery)
  rm(AY_negative, AY_discovery, AY_positive)
}

write.csv(AY, sprintf('%s/AY/%s/AY.csv', save_dir, AY_SETTING_NAME), row.names=FALSE)






# ============================================================================== #
# ============= Initial LinReg and NegBin:                  ====================
# ============= Compare across control cell subsample sizes ====================
# ============================================================================== #

# # positive tests across different NUM_SUBSAMPLE_CONTROL_CELLS sizes
# if('initial_glms_posneg_comparecontrolsubsamplesizes.csv' %in% list.files(sprintf('%s/AY/%s/initiallinreg/', save_dir, AY_SETTING_NAME))) {
#   # just skip if already done
#   initial_glms_pos_df = read.csv(sprintf('%s/AY/%s/initiallinreg/initial_glms_posneg_comparecontrolsubsamplesizes.csv', save_dir, AY_SETTING_NAME))
#   
# } else {
#   # realize some for this comparison
#   required_Y = unique(AY$gene) # only these response genes' expression are required
#   gene_counts_subset = gene_counts[required_Y, ] |> DelayedArray::realize()
#   gene_norm_subset   =   gene_norm[required_Y, ] |> DelayedArray::realize()
#   
#   # prep some used values
#   # combined_grna = rbind(grna, grnaNT)
#   control_cell_idx = which(grna['Non_target', ] == 1) |> as.vector() # just the overall control cell idx. use setdiff for negative tests
#   libsize_log = log(metadata[, 'num_rna_umi'] + 1e-5)
#   
#   fit_glms_errorhandling <- function(AY_idx, subsamplesize) {
#     tryCatch({fit_glms(AY_idx=AY_idx, 
#                        AY=AY, 
#                        control_cell_idx=control_cell_idx, 
#                        gene_norm=gene_norm_subset, 
#                        gene_counts=gene_counts_subset, 
#                        grna=grna, 
#                        libsize_log=libsize_log, 
#                        covariates=NULL, 
#                        subsample_control_size=subsamplesize)
#     },
#     error = function(cond) {
#       return(NULL)
#     })
#   }
# }



# ============================================================================== #
# ============= Realize Relevant Gene info in memory: ==========================
# ============================================================================== #
required_Y = unique(AY$gene) # only these response genes' expression are required
required_Y = unique(c(required_Y, validsinglegenes)) |> sort() # also for single genes


# # separate NTs for control and for null tests
# subsample_ctrl_cells = sample(which(colSums(grnaNT[NT_grnas_for_control, ])==1)  |> names(), round(NUM_SUBSAMPLE_CONTROL_CELLS*1.1), replace=FALSE)# take a few more than necessary
# subsample_cells = c(which(grna['Non_target', ]==0) |> names(), 
#                     which(colSums(grnaNT[NT_grnas_for_null, ])==1)  |> names(),
#                     subsample_ctrl_cells) # all treated cells + null pert cells + the prev sampled cells

# if NTs were mixed up
# subsample_ctrl_cells = sample(which(grna['Non_target', ]==1) |> names(), round(NUM_SUBSAMPLE_CONTROL_CELLS*2), replace=FALSE)# take a few more than necessary
# subsample_cells = c(which(grna['Non_target', ]==0) |> names(), subsample_ctrl_cells) # all treated cells + the prev sampled cells

gene_counts_subset = gene_counts[required_Y, ] |> DelayedArray::realize()
gene_norm_subset   =   gene_norm[required_Y, ] |> DelayedArray::realize()

# rm(subsample_ctrl_cells)
# rm(gene_counts, gene_norm)
# gc()


format(object.size(gene_counts_subset), units='Gb')
format(object.size(gene_norm_subset), units='Gb')

# ============================================================================== #
# ============= Prep variables: ================================================
# ============================================================================== #

# prep reused variables in estimating fns
control_cell_idx = which(grna['Non_target', ] == 1) |> as.vector() 
libsize_log      = log(metadata[, 'num_rna_umi'] + 1e-5)

# prep reused variables in estimating fns (prox)
NCs_pca = readRDS(sprintf('%s/pca/NCloadings.rds', save_dir)) |> as.matrix() # PCA loadings

rm(metadata)
gc()

# ================ checks
assertthat::assert_that(all(row.names(gene_norm_subset) == row.names(gene_counts_subset))) # same genes and in same order


# check cellnames 
assertthat::assert_that(all(colnames(gene_norm_subset) == colnames(gene_counts_subset))) # same cells and in same order
# assertthat::assert_that(all(colnames(gene_norm_subset) == colnames(combined_grna))) # same cells and in same order



# ============================================================================== #
# ============= Initial LinReg and NegBin:                    ==================
# ============= after choosing a control cell subsample size  ==================
# ============================================================================== #

# if('initial_glms.csv' %in% list.files(sprintf('%s/AY/%s/initiallinreg/', save_dir, AY_SETTING_NAME))) {
#   print(sprintf("[%s]        - Loading glm (one fit per gene) on all pairs", Sys.time()))
#   initial_glms = read.csv(sprintf('%s/AY/%s/initiallinreg/initial_glms.csv', save_dir, AY_SETTING_NAME))
# } else {
#   print(sprintf("[%s]        - Fitting glm (one fit per gene) on all pairs", Sys.time()))
#   
#   fit_glms_errorhandling <- function(AY_idx) {
#     tryCatch({fit_glms(AY_idx=AY_idx, 
#                        AY=AY, 
#                        control_cell_idx=control_cell_idx, 
#                        gene_norm=gene_norm_subset, 
#                        gene_counts=gene_counts_subset, 
#                        grna=grna, 
#                        libsize_log=libsize_log, 
#                        covariates=NULL, 
#                        subsample_control_size=NA)
#     },
#     error = function(cond) {
#       return(NULL)
#     })
#   }
#   
#   
#   t0 = Sys.time()
#   initial_glms_par = future.apply::future_mapply(
#     fit_glms_errorhandling,
#     AY_idx = 1:nrow(AY), 
#     future.globals = TRUE,
#     future.seed = 123456, 
#     SIMPLIFY = FALSE)
#   t1 = Sys.time()
#   print(sprintf("[%s]        - %2.2f mins", Sys.time(), difftime(t1, t0, units = "mins")))
#   
#   initial_glms = NULL
#   for(idx in 1:length(initial_glms_par)) {   initial_glms = rbind(initial_glms, initial_glms_par[[idx]])   }
#   write.csv(initial_glms, sprintf('%s/AY/%s/initiallinreg/initial_glms.csv', save_dir, AY_SETTING_NAME))
#   rm(initial_glms_par, t0, t1, fit_glms_errorhandling)
# }





# # 134k if use all. 50000 speeds up a lot, and lower actually loses power... for now, set at 50k cells
# fit_glms(AY_idx=4, AY=AY, control_cell_idx=control_cell_idx, gene_norm=gene_norm, gene_counts=gene_counts, grna=grna, libsize_log=libsize_log, covariates=NULL, 
#          subsample_control_size=NA)
# fit_glms(AY_idx=4, AY=AY, control_cell_idx=control_cell_idx, gene_norm=gene_norm, gene_counts=gene_counts, grna=grna, libsize_log=libsize_log, covariates=NULL, 
#          subsample_control_size=50000)
# fit_glms(AY_idx=4, AY=AY, control_cell_idx=control_cell_idx, gene_norm=gene_norm, gene_counts=gene_counts, grna=grna, libsize_log=libsize_log, covariates=NULL, 
#          subsample_control_size=10000)
# fit_glms(AY_idx=4, AY=AY, control_cell_idx=control_cell_idx, gene_norm=gene_norm, gene_counts=gene_counts, grna=grna, libsize_log=libsize_log, covariates=NULL, 
#          subsample_control_size=5000)



# ============================================================================== #
# =================== Get ATEs by gene (parallel) ==============================
# ============================================================================== #
print(sprintf("[%s]    - Get Estimates by gene (parallel)", Sys.time()))

print(sprintf("[%s]        + nrow(AY)=%s", Sys.time(), nrow(AY)))
print(sprintf("[%s]        + #unique grna=%s", Sys.time(), length(unique(AY$grna))))
print(sprintf("[%s]        + #unique gene=%s", Sys.time(), length(unique(AY$gene))))

# source(sprintf('%s/quick_estimate/quickPCI2S_bygene.R', util_dir))

effects_df = NULL
for(Y_name in unique(AY$gene)) {
  print(sprintf("[%s]    - Y: %s", Sys.time(), Y_name))
  Y_info = list(Y_name = Y_name, 
                Y_norm = gene_norm[Y_name, ], # replace w gene_norm_subset maybe
                Y_count = gene_counts[Y_name, ])
  # choose a smaller set, but still larger than max numNCs so that there is less randomness from bad selection
  Y_validsinglegenes = gene_norm[sample(setdiff(validsinglegenes, Y_name), size = round(2*2*max(SETTING$numNCs)), replace = FALSE), ]
  
  # estimate effects: parallel calls within this function
  cur_effects_df = tryCatch({estimate_effects_bygene(Y_info=Y_info, 
                                                     # test information
                                                     AY=AY, 
                                                     # proximal number of NCs
                                                     numNCs=SETTING$numNCs, 
                                                     # data
                                                     control_cell_idx=control_cell_idx, 
                                                     grna=grna,
                                                     libsize_log=libsize_log, 
                                                     # info for NCs or covariates
                                                     NCs_pca=NCs_pca,
                                                     Y_validsinglegenes=Y_validsinglegenes,
                                                     covariates=NULL,
                                                     # saving intermediateATEs
                                                     save_intermediateATEs=TRUE, 
                                                     skip_if_saved_intermediateATEs=SETTING$skip_if_saved_intermediateATEs,
                                                     intermediateATEs_folder = sprintf('%s/AY/%s/intermediateATEs/', save_dir, AY_SETTING_NAME),
                                                     # which proximal estimates
                                                     run_proximal_linear_pca         = SETTING$prox_methods$run_proximal_linear_pca,
                                                     run_proximal_negbin_pca         = SETTING$prox_methods$run_proximal_negbin_pca,
                                                     run_proximal_linear_singlegenes = SETTING$prox_methods$run_proximal_linear_singlegenes,
                                                     run_proximal_negbin_singlegenes = SETTING$prox_methods$run_proximal_negbin_singlegenes,
                                                     subsample_control_size = NA)
  },
  error = function(cond) {
    # message(sprintf('Error est ATE w pci2s with %s', 
    #                 AY_idx)) 
    return(NULL)
  })
  
  if(!is.null(cur_effects_df)) {
    effects_df = rbind(effects_df, cur_effects_df)
  } 
}


# ============================================================================== #
# =================== Get ATEs by AY_idx (parallel) ============================
# ============================================================================== #
if(F) {
# numNCs = c(1,3)
# 
# control_cell_idx = which(grna['Non_target', ] == 1) |> as.vector() # just the overall control cell idx. use setdiff for negative tests
# 
# combined_grna = rbind(grna, grnaNT)
# 
# libsize_log = log(metadata[, 'num_rna_umi'] + 1e-5) # constant? 1e-5 for now, could also just do 1
# # libsize_log = rep(NA, ncol(grna))
# 
# NCs_pca = readRDS(sprintf('%s/pca/NCloadings.rds', save_dir)) |> as.matrix() # PCA loadings


# some setup
estimate_effects_errorhandling <- function(AY_idx) {
  res_df = tryCatch({estimate_effects(AY_idx=AY_idx, 
                                      # test information
                                      AY=AY, 
                                      # proximal number of NCs
                                      numNCs=SETTING$numNCs, 
                                      # data
                                      control_cell_idx=control_cell_idx, 
                                      gene_norm=gene_norm_subset, 
                                      gene_counts=gene_counts_subset, 
                                      grna=grna,
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
                                      run_proximal_linear_pca         = SETTING$prox_methods$run_proximal_linear_pca,
                                      run_proximal_negbin_pca         = SETTING$prox_methods$run_proximal_negbin_pca,
                                      run_proximal_linear_singlegenes = SETTING$prox_methods$run_proximal_linear_singlegenes,
                                      run_proximal_negbin_singlegenes = SETTING$prox_methods$run_proximal_negbin_singlegenes,
                                      subsample_control_size = NA)
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


# source(sprintf('%s/quick_estimate/quickPCI2S.R', util_dir))
# 
# # test a run
# a = estimate_effects(AY_idx=2,
#                  # test information
#                  AY=AY,
#                  # proximal number of NCs
#                  numNCs=c(1, 3),
#                  # data
#                  control_cell_idx=control_cell_idx,
#                  gene_norm=gene_norm_subset,
#                  gene_counts=gene_counts_subset,
#                  grna=grna,
#                  libsize_log=libsize_log,
#                  # info for NCs or covariates
#                  NCs_pca=NCs_pca,
#                  validsinglegenes=validsinglegenes,
#                  covariates=NULL,
#                  # saving intermediateATEs
#                  save_intermediateATEs=TRUE,
#                  skip_if_saved_intermediateATEs=FALSE,
#                  intermediateATEs_folder = sprintf('%s/AY/%s/intermediateATEs/', save_dir, AY_SETTING_NAME),
#                  # which proximal estimates
#                  run_proximal_linear_pca         = TRUE,
#                  run_proximal_negbin_pca         = TRUE,
#                  run_proximal_linear_singlegenes = TRUE,
#                  run_proximal_negbin_singlegenes = TRUE,
#                  subsample_control_size = 25000)
# 
# estimate_effects_errorhandling(2)
# 
# a = combined_grna[c("Non_target_399", "Non_target_313"), ]





print(sprintf("[%s]    - Get Estimates (parallel)", Sys.time()))

# whichROWS = 1:60
# whichROWS = 70:80
whichROWS = 1:nrow(AY)
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
    # AY test is whichROWS[whichROWS_idx]  = AY_idx 
    effects_df = dplyr::bind_rows(effects_df, cur_effects_par)
  }
  rm(cur_effects_par)
}
write.csv(x = effects_df, file = sprintf('%s/AY/%s/effects.csv', save_dir, AY_SETTING_NAME), row.names = FALSE)



}






# ===========================================================================
print(sprintf("[%s] END", Sys.time()))



# ===========================================================================
# ===========================================================================
# ===========================================================================
# ===========================================================================
# ===========================================================================

if(F) {
# writing function


AY = rbind(AY, data.frame(grna = c("Zbtb7a", 'Uba5', "Non_target_64&Non_target_370&Non_target_4&Non_target_219&Non_target_107&Non_target_290&Non_target_246&Non_target_349&Non_target_157&Non_target_319"),
                          gene = c('Prdm10', 'Prdm10', 'Prdm10'), type = 'discovery')) |> distinct()

# writing function

Y_infos = list()
for(Y_name in unique(AY$gene)) {
  Y_infos[[Y_name]] = list(Y_name = Y_name, 
                           Y_norm = gene_norm[Y_name, ],
                           Y_count = gene_counts[Y_name, ])
  # choose a smaller set, but still larger than max numNCs so that there is less randomness from bad selection
  Y_validsinglegenes = gene_norm[sample(setdiff(validsinglegenes, Y_name), size = round(1.5*2*max(SETTING$numNCs)), replace = FALSE), ]
  
  estimate_effects_bygene()
}





  

Y_info = Y_infos[[33]]
Y_info$Y_name



# === start fn
res = NULL


A_names = AY |> filter(gene == Y_info$Y_name) |> pull(grna)
A_names = c(A_names, 'Stxbp1&Srebf2')

# construct A's and collect the used cell idx

# for the treatments already in grna
# A_01 = grna[A_name[A_name %in% row.names(grna))], ] 

# the indices of the chosen cells (cells receiving this perturbation)
AY_data_idx = which(colSums(grna[A_names[A_names %in% row.names(grna)], ]) == 1) |> unname()
# the A treatment in factor form
A_factor = row.names(grna)[apply(grna[, AY_data_idx], MARGIN=2, FUN = function(x){which(x==1)})]


control_cell_idx_noNTtrt = control_cell_idx
# for the treatments not already in grna 
# needs to handle combined perturbations that aren't pre-constructed in grna (e.g. some NT1&NT2 are already constructed in grna. But some are)
for(A_name in A_names[! A_names %in% row.names(grna)]) {
  # needs to be w pattern <name>&<name>... where <name> is a member of grna row.names
  # parse NT combo name (e.g. 'NT1&NT2')
  grnacombo_vector = strsplit(A_name, '&')[[1]]
  
  grnacombo_01 = colSums(grna[grnacombo_vector, , drop=FALSE]) |> as.vector()
  grnacombo_idx = which(grnacombo_01 == 1) |> unname()
  
  control_cell_idx_noNTtrt = setdiff(control_cell_idx_noNTtrt, grnacombo_idx) # remove these from control cell set
  
  AY_data_idx = c(AY_data_idx, grnacombo_idx)
  A_factor = c(A_factor, rep(A_name, times = length(grnacombo_idx)))
  
  rm(grnacombo_vector, grnacombo_01, grnacombo_idx)
}

# deal with subsampling if requested
if(is.na(subsample_control_size) | (subsample_control_size >= length(control_cell_idx_noNTtrt))) {
  selected_control_cells = control_cell_idx_noNTtrt
} else {
  selected_control_cells = sample(x=control_cell_idx_noNTtrt, size=subsample_control_size, replace=FALSE)
}




AY_data_idx = c(AY_data_idx, selected_control_cells)
A_factor    = c(A_factor, rep('control_cell', times = length(selected_control_cells))) |> factor() |> relevel('control_cell')

A_01 = model.matrix(~ A_ - 1, data = data.frame(A_=A_factor), check.names=FALSE)

df = data.frame(
  A_01,
  Y_count     = Y_info$Y_count[AY_data_idx], 
  Y_norm      = Y_info$Y_norm[AY_data_idx],
  libsize_log = libsize_log[AY_data_idx], 
  check.names=FALSE
)

# add covariates to 
if(!is.null(covariates)) {
  covariates = covariates[AY_data_idx, ]
  colnames(covariates) = paste0('X', 1:ncol(covariates)) # rename columns to avoid potential naming mixups/collisions
  df = cbind(df, covariates)
  cov_suffix = paste0(paste0('+', colnames(covariates)), collapse = '') 
} else {
  covariates = NULL # set NULL for input for pci2s X
}


# === start over...


# map of A_name to the cell indicies that received this perturbation (or combo of perturbations)
A_map_idx = list()
for(A_name in A_names) {
  if(A_name %in% row.names(grna)) {
    A_map_idx[[A_name]]= which(grna[A_name, ] == 1) |> unname()
  } else {
    # needs to be w pattern <name>&<name>... where <name> is a member of grna row.names
    # parse NT combo name (e.g. 'NT1&NT2')
    grnacombo_vector = strsplit(A_name, '&')[[1]]
    
    grnacombo_01 = colSums(grna[grnacombo_vector, , drop=FALSE]) |> as.vector()
    grnacombo_idx = which(grnacombo_01 == 1) |> unname()
    
    A_map_idx[[A_name]] = grnacombo_idx
    
    control_cell_idx_noNTtrt = setdiff(control_cell_idx_noNTtrt, grnacombo_idx) # remove these from control cell set
    
    rm(grnacombo_vector, grnacombo_01, grnacombo_idx)
  }
}




run_proximal_linear_pca         = TRUE
run_proximal_negbin_pca         = TRUE
run_proximal_linear_singlegenes = TRUE
run_proximal_negbin_singlegenes = TRUE


# create a function to estimate per specified treatment 
inner_estimate_effect <- function(A_name) {
  # A_name = A_names[1]
  
  res = NULL
  
  # current used cell idx for this test
  cur_idx = c(control_cell_idx_noNTtrt, A_map_idx[[A_name]])
  
  df = data.frame(
    A = c(rep(0, times = length(control_cell_idx_noNTtrt)), 
          rep(1, times = length(A_map_idx[[A_name]]))),
    Y_count     = Y_info$Y_count[cur_idx], 
    Y_norm      = Y_info$Y_norm[cur_idx],
    libsize_log = libsize_log[cur_idx], 
    check.names=FALSE
  )
  
  # add covariates to 
  lm_formula     = 'Y_norm ~ A'
  negbin_formula = 'Y_count ~ A + offset(libsize_log)'
  if(!is.null(covariates)) {
    cur_covariates = covariates[cur_idx, ]
    colnames(cur_covariates) = paste0('X', 1:ncol(cur_covariates)) # rename columns to avoid potential naming mixups/collisions
    
    # append to df
    df = cbind(df, cur_covariates)
    # append to glm formulas
    cov_suffix = paste0(paste0('+', colnames(cur_covariates)), collapse = '') 
    lm_formula     = paste0(    lm_formula, cov_suffix)
    negbin_formula = paste0(negbin_formula, cov_suffix)
    rm(cov_suffix)
  } else {
    cur_covariates = NULL # set NULL for input for pci2s X
  }
  
  
  # ===== estimate =====
  # ===== \__Linear Model Fit =====
  t0 = Sys.time()
  lmfit = lm(formula=lm_formula, data = df)
  lmfit_summary = summary(lmfit)
  t1 = Sys.time()
  
  res = bind_rows(res, 
                  data.frame(
                    NC_type     = NA,
                    method      = 'lm',
                    method_type = 'glm',
                    numNC       = NA,
                    ATE = lmfit_summary$coefficients['A', 'Estimate'],
                    se  = lmfit_summary$coefficients['A', 'Std. Error'],
                    tstat= lmfit_summary$coefficients['A', 't value'],
                    pval= lmfit_summary$coefficients['A', 'Pr(>|t|)'],
                    time_sec = difftime(t1, t0, units = 'secs') |> as.numeric()))
  rm(lmfit, lmfit_summary, t1, t0)
  # return(res)
  
  # ===== \__ Negative Binomial Fit  ====
  t0 = Sys.time()
  # nbfit = MASS::glm.nb(formula='Y ~ A + offset(libsize_log)', data = df)
  # nbfit_summary = summary(nbfit)
  nbfit = tryCatch({  # give negbin results if works
    MASS::glm.nb(formula=negbin_formula, data = df)
  },
  error = function(cond) { # give NA results if errored
    return(NULL)
  })
  t1 = Sys.time()
  if(is.null(nbfit)) {
    res = bind_rows(res, 
                    data.frame(
                      NC_type     = NA,
                      method      = 'negbin',
                      method_type = 'glm',
                      numNC       = NA,
                      ATE = NA,
                      se  = NA,
                      tstat= NA,
                      pval= NA,
                      time_sec = difftime(t1, t0, units = 'secs') |> as.numeric()))
  } else {
    nbfit_summary = summary(nbfit)
    res = bind_rows(res, 
                    data.frame(
                      NC_type     = NA,
                      method      = 'negbin',
                      method_type = 'glm',
                      numNC       = NA,
                      ATE = nbfit_summary$coefficients['A', 'Estimate'],
                      se  = nbfit_summary$coefficients['A', 'Std. Error'],
                      tstat= nbfit_summary$coefficients['A', 'z value'],
                      pval= nbfit_summary$coefficients['A', 'Pr(>|z|)'],
                      time_sec = difftime(t1, t0, units = 'secs') |> as.numeric()))
    rm(nbfit_summary)
  }
  
  
  # res = bind_rows(res, 
  #                 data.frame(
  #                   NC_type     = NA,
  #                   method      = 'negbin',
  #                   method_type = 'glm',
  #                   numNC       = NA,
  #                   ATE = nbfit_summary$coefficients['A', 'Estimate'],
  #                   se  = nbfit_summary$coefficients['A', 'Std. Error'],
  #                   tstat= nbfit_summary$coefficients['A', 'z value'],
  #                   pval= nbfit_summary$coefficients['A', 'Pr(>|z|)'],
  #                   time_sec = difftime(t1, t0, units = 'secs') |> as.numeric()))
  rm(nbfit, t1, t0)
  
  # ===== \__ Proximal PCA ====
  # Proximal with PCA loadings
  dfZ = NCs_pca[cur_idx, seq(from = 2, to = min(2*max(numNCs), ncol(NCs_pca)), by = 2)] # evens which((1:ncol(NCs)) %% 2 == 0)
  dfW = NCs_pca[cur_idx, seq(from = 1, to = min(2*max(numNCs), ncol(NCs_pca)), by = 2)] # odds
  colnames(dfZ) = paste0('Z', 1:ncol(dfZ))
  colnames(dfW) = paste0('W', 1:ncol(dfW))
  for(numNCs_ in numNCs)  {
    if(run_proximal_linear_pca) {
      # ===== \__\__Proximal Linear ====
      t0 = Sys.time()
      pci2s_res = tryCatch({  # give pci2s results if works
        pci2s::p2sls.lm(
          Y = df$Y_norm, 
          A = df$A, 
          W = dfZ[,1:numNCs_] |> as.data.frame(), 
          Z = dfW[,1:numNCs_] |> as.data.frame(), 
          X = cur_covariates,
          variance = TRUE)
      },
      error = function(cond) { # give NA results if errored
        return(list(summary_second_stage=
                      matrix(c(NA, NA, NA, NA), nrow=1, ncol=4, 
                             dimnames = list(c('A'), c('Estimate', 'Std. Error', 'z value', 'Pr(>|z|)')))))
      })
      t1 = Sys.time()
      
      res = bind_rows(res, 
                      data.frame(
                        NC_type     = 'PCA',
                        method      = 'proximallinear',
                        method_type = 'proximal',
                        numNC       = numNCs_,
                        ATE  = pci2s_res$summary_second_stage['A', 'Estimate'],
                        se   = pci2s_res$summary_second_stage['A', 'Std. Error'],
                        tstat= pci2s_res$summary_second_stage['A', 'z value'],
                        pval = pci2s_res$summary_second_stage['A', 'Pr(>|z|)'],
                        time_sec = difftime(t1, t0, units = 'secs') |> as.numeric()))
      
      rm(pci2s_res, t0, t1)
    }
    
    # ===== \__\__Proximal Negative Binomial ====
    if(run_proximal_negbin_pca){   
      # return(df)
      t0 = Sys.time()
      pci2s_res = tryCatch({  # give pci2s results if works
        pci2s::p2sls.negbin(
          Y = df$Y_count, 
          A = df$A, 
          W = dfZ[,1:numNCs_] |> as.data.frame(), 
          Z = dfW[,1:numNCs_] |> as.data.frame(), 
          X = cur_covariates,
          offset = df$libsize_log,
          nco_type = rep("linear", numNCs_),
          # nco_args = lapply(X = 1:num_NCs, FUN = function(x){list(init=NA, offset=log(df_all$library_size))}), # No offsets for continuous versions
          variance = TRUE,
          verbose = FALSE)
      },
      error = function(cond) { # give NA results if errored
        return(list(summary_second_stage=
                      matrix(c(NA, NA, NA, NA), nrow=1, ncol=4, 
                             dimnames = list(c('A'), c('Estimate', 'Std. Error', 'z value', 'Pr(>|z|)')))))
      })
      t1 = Sys.time()
      
      res = bind_rows(res, 
                      data.frame(
                        NC_type     = 'PCA',
                        method      = 'proximalnegbin',
                        method_type = 'proximal',
                        numNC       = numNCs_,
                        ATE  = pci2s_res$summary_second_stage['A', 'Estimate'],
                        se   = pci2s_res$summary_second_stage['A', 'Std. Error'],
                        tstat= pci2s_res$summary_second_stage['A', 'z value'],
                        pval = pci2s_res$summary_second_stage['A', 'Pr(>|z|)'],
                        time_sec = difftime(t1, t0, units = 'secs') |> as.numeric()))
      
      rm(pci2s_res, t0, t1)
    }
    
  }
  rm(dfZ, dfW, numNCs_)
  
  
  # ===== \__ Proximal singlegenes ====
  if(run_proximal_linear_singlegenes | run_proximal_negbin_singlegenes) {
    if(is.null(Y_validsinglegenes)) {
      print('Error with prox singlegenes: did not give Y_validsinglegenes input')
    } else {
      # print('run_proximal_linear_singlegenes')
      # ===== \__\__Proximal Linear ====
      NCs_singlegene = Y_validsinglegenes[sample(nrow(Y_validsinglegenes), 2*max(numNCs)), cur_idx] |> t() |> as.data.frame()
      dfZ = NCs_singlegene[ , seq(from = 2, to = min(2*max(numNCs), ncol(NCs_singlegene)), by = 2)] # evens which((1:ncol(NCs)) %% 2 == 0)
      dfW = NCs_singlegene[ , seq(from = 1, to = min(2*max(numNCs), ncol(NCs_singlegene)), by = 2)] # odds
      colnames(dfZ) = paste0('Z', 1:ncol(dfZ))
      colnames(dfW) = paste0('W', 1:ncol(dfW))
      rm(NCs_singlegene)
      for(numNCs_ in numNCs)  {
        if(run_proximal_linear_singlegenes) {
          t0 = Sys.time()
          pci2s_res = tryCatch({  # give pci2s results if works
            pci2s::p2sls.lm(
              Y = df$Y_norm, 
              A = df$A, 
              W = dfZ[,1:numNCs_], 
              Z = dfW[,1:numNCs_], 
              X = cur_covariates,
              variance = TRUE)
          },
          error = function(cond) { # give NA results if errored
            return(list(summary_second_stage=
                          matrix(c(NA, NA, NA, NA), nrow=1, ncol=4, 
                                 dimnames = list(c('A'), c('Estimate', 'Std. Error', 'z value', 'Pr(>|z|)')))))
          })
          t1 = Sys.time()
          
          res = bind_rows(res, 
                          data.frame(
                            NC_type     = 'singlegene',
                            method      = 'proximallinear',
                            method_type = 'proximal',
                            numNC       = numNCs_,
                            ATE  = pci2s_res$summary_second_stage['A', 'Estimate'],
                            se   = pci2s_res$summary_second_stage['A', 'Std. Error'],
                            tstat= pci2s_res$summary_second_stage['A', 'z value'],
                            pval = pci2s_res$summary_second_stage['A', 'Pr(>|z|)'],
                            time_sec = difftime(t1, t0, units = 'secs') |> as.numeric()))
          
          rm(pci2s_res, t0, t1)
        }
        
        
        # ===== \__\__Proximal Negative Binomial ====
        if(run_proximal_negbin_singlegenes){ 
          t0 = Sys.time()
          pci2s_res = tryCatch({  # give pci2s results if works
            pci2s::p2sls.negbin(
              Y = df$Y_count, 
              A = df$A, 
              W = dfZ[,1:numNCs_] |> as.data.frame(), 
              Z = dfW[,1:numNCs_] |> as.data.frame(), 
              X = cur_covariates,
              offset = df$libsize_log,
              nco_type = rep("linear", numNCs_),
              # nco_args = lapply(X = 1:num_NCs, FUN = function(x){list(init=NA, offset=log(df_all$library_size))}), # No offsets for continuous versions
              variance = TRUE,
              verbose = FALSE)
          },
          error = function(cond) { # give NA results if errored
            return(list(summary_second_stage=
                          matrix(c(NA, NA, NA, NA), nrow=1, ncol=4, 
                                 dimnames = list(c('A'), c('Estimate', 'Std. Error', 'z value', 'Pr(>|z|)')))))
          })
          t1 = Sys.time()
          
          res = bind_rows(res, 
                          data.frame(
                            NC_type     = 'singlegene',
                            method      = 'proximalnegbin',
                            method_type = 'proximal',
                            numNC       = numNCs_,
                            ATE  = pci2s_res$summary_second_stage['A', 'Estimate'],
                            se   = pci2s_res$summary_second_stage['A', 'Std. Error'],
                            tstat= pci2s_res$summary_second_stage['A', 'z value'],
                            pval = pci2s_res$summary_second_stage['A', 'Pr(>|z|)'],
                            time_sec = difftime(t1, t0, units = 'secs') |> as.numeric()))
          
          rm(pci2s_res, t0, t1)
        }
      }
      rm(dfZ, dfW, numNCs_)
    } # end numNCs_
  }   # end singlegenes
  
  
  # ===== return ==== 
  res = res |> mutate(grna = A_name, gene = Y_name, 
                      type = AY |> filter(grna == A_name & gene == Y_name) |> head(1) |> pull(type), 
                      n_trtmnt = sum(df$A==1), n_cntrl=sum(df$A==0), .before = 1)
  
  if(save_intermediateATEs) {
    write.csv(res, sprintf('%s/ATE_%s.csv', intermediateATEs_folder, which((AY$grna == A_name) & (AY$gene == Y_name))))
  }
  return(res)
}

res_all = NULL
for(A_name in A_names) {
  inner_estimate_effect(A_name)
}


res_par = future.apply::future_mapply(inner_estimate_effect,
                                          A_name = A_names, 
                                          future.globals = TRUE,
                                          future.seed = 123456, 
                                          SIMPLIFY = FALSE)


res_df = NULL
for(i in 1:length(res_par)) { 
  cur_effects_par = res_par[[i]]
  if(!is.null(cur_effects_par)) {
    res_df = dplyr::bind_rows(res_df, cur_effects_par)
  }
  rm(cur_effects_par)
}

return(res_df)
}

