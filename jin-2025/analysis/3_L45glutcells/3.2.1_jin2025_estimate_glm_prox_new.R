# ------------------------------------------------------------------------------------------------ #
#     Initial Estimates on Jin-2025 dataset on glut cell type
#     Mouse brain data from Xin Jin's lab ~ 2025 and downloaded from dropbox 
# Email: Ergan to Jin-Hong to Catherine:
#     Please check the following two links of the data available. 
#     Those data is shared in 2 batches by Seoyeon, where one is containing 10 perturbations and another one contains 48 perturbations.
#     
# 
# 3.1 preprocessing
# - Pre-processing (gene importance and normalization)
# - Perform PCA 
# 3.2 estimates
# - 2.2.1 glms + proximal
# - 2.2.2 deseq2
# - 2.2.3 causarray
# - 2.2.4 sceptre
# Rscript 
# ------------------------------------------------------------------------------------------------ #
args = commandArgs(trailingOnly = TRUE)
# args = c('macbook', 'small2',  '8')


assertthat::assert_that(length(args) > 0, msg="must give arg for specifying device eg 'Rscript <filename>.R ubergenno'")
DEVICE = args[1]
source('../PATHS.R') # load in data_dir and save_dir and CODE_DIR, depending on DEVICE value
assertthat::assert_that(!is.null(data_dir), msg='first arg must be: laptop, desktop, or ubergenno')

assertthat::assert_that(length(args) > 1, msg="must give arg for specifying AY setting eg 'Rscript <filename>.R ubergenno small'")
AY_SETTING_NAME = args[2]


assertthat::assert_that(length(args) > 2, msg="must give arg for specifying num threads eg 'Rscript <filename>.R ubergenno small 8'")
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

if(THREADS == 1) {
  plan(sequential)
} else {
  # plan(multisession, workers = 24)
  plan(multisession, workers = THREADS)
}



PARAMS = list( 
  'TAG'       = 'L45Glut',
  'subfolder' = '3_L45glutcells',
  'seed'      = 12345,
  'A_MIN_TRTCELLS'  = 40, # filter for grnas with at least this number of cells (keep low here, set higher later))
  'Y_MAX_GENE_IMPORTANCE' = 1500
)
save_dir = paste0(save_dir, '/', PARAMS[['subfolder']])
set.seed(PARAMS[['seed']])




# SKIP_IF_SAVED_INTERMEDIATEATES = TRUE
SKIP_IF_SAVED_INTERMEDIATEATES = FALSE

# Set Parameters for Estimation: Choosing AY 
# AY_SETTING_NAME = 'denseSmall' # key of one of the following settings, prob change this to an arg parameter
# for non-dense settings, probably do sequential by setting THREADS=1
AY_SETTINGS = list(
  'small' = list(
    AY_NUM_NEGATIVE  = 150,
    AY_NUM_DISCOVERY = 150,
    numNCs = c(1, 3, 5, 10)
  ),
  'small2' = list( # very small, to test all of the positive tests and many methods and parameter settings
    AY_NUM_NEGATIVE  = 100,
    AY_NUM_DISCOVERY = 100,
    numNCs = c(1, 3, 5, 10)
  ),
  'smallNB' = list(
    AY_NUM_NEGATIVE  = 150,
    AY_NUM_DISCOVERY = 150,
    numNCs = c(1, 3, 5, 10)
  ),
  
  'medium' = list(
    AY_NUM_NEGATIVE  = 1000,
    AY_NUM_DISCOVERY = 1000,
    numNCs = c(1, 3, 5, 8, 10, 15)
  ),
  
  'denseSmall' = list(
    DENSE = TRUE,
    NUM_Y = 100,
    AY_NUM_NEGATIVE  = NA,
    AY_NUM_DISCOVERY = NA,
    numNCs = c(1, 3, 5, 8, 10, 15)
  ), 
  
  'denseMedium' = list(
    DENSE = TRUE,
    NUM_Y = 1000,
    AY_NUM_NEGATIVE  = NA,
    AY_NUM_DISCOVERY = NA,
    numNCs = c(1, 3, 5, 10)
  )
  
)


assertthat::assert_that(AY_SETTING_NAME %in% names(AY_SETTINGS), msg=sprintf("AY_SETTINGS input must be a member of the available setting names. Available AY_SETTINGS are: %s", 
                                                                             paste0(names(AY_SETTINGS), collapse = ', ')))
SETTING = AY_SETTINGS[[AY_SETTING_NAME]]


DEFAULT_NUM_NC = NULL
DEFAULT_NUM_NC = c(1, 3, 5, 8, 10, 15)



# create a list of settings for which estimators to run
# named list where the names are the estimator method name
# each element is a list with at least the named elements:
#   RUN=bool whether to run this method or not
#   Y_TYPE=type of Y
# for proximal methods
ESTIMATOR_SETTINGS = list(
  'lmYA' =list(RUN=TRUE,
               Y_TYPE='continuous'), 
  'lmYAX'=list(RUN=FALSE,
               Y_TYPE='continuous'), 
  'nbYA' =list(RUN=TRUE,
               Y_TYPE='count'), 
  'nbYAX'=list(RUN=FALSE,
               Y_TYPE='count'), 
  'PCA_UntrtCells' =list(
    RUN=FALSE,
    IS_PROX=TRUE, 
    NUM_NCs=if(!is.null(DEFAULT_NUM_NC)) {DEFAULT_NUM_NC} else {c(1, 3, 5, 8, 10, 15)},
    PROVIDES_NC_MATRIX=TRUE,   # in other list, must provide matrix of NCs (false=vector of gene names for singlegenes)
    NC_RANDOMIZEORDER=FALSE,
    PROXIMAL_METHOD = 'proximallinear',
    NC_type = 'PCA_UntrtCells',
    Y_TYPE='continuous'
  ),
  'PCA_MixCells' =list(
    RUN=TRUE,
    IS_PROX=TRUE,
    NUM_NCs=if(!is.null(DEFAULT_NUM_NC)) {DEFAULT_NUM_NC} else {c(1, 3, 5, 8)},
    PROVIDES_NC_MATRIX=TRUE,   # in other list, must provide matrix of NCs (false=vector of gene names for singlegenes)
    NC_RANDOMIZEORDER=FALSE,
    PROXIMAL_METHOD = 'proximallinear',
    NC_type = 'PCA_MixCells',
    Y_TYPE='continuous'
  ),
  'PCA_MixCells_NB' =list(
    RUN=FALSE,
    IS_PROX=TRUE,
    NUM_NCs=if(!is.null(DEFAULT_NUM_NC)) {DEFAULT_NUM_NC} else {c(1, 3, 5, 8, 10, 15)},  
    PROVIDES_NC_MATRIX=TRUE,   # in other list, must provide matrix of NCs (false=vector of gene names for singlegenes)
    NC_RANDOMIZEORDER=FALSE,
    PROXIMAL_METHOD = 'proximalnegbin',
    NC_type = 'PCA_MixCells',
    Y_TYPE='count'
  ),
  'singlegene'= list(
    RUN=TRUE,
    IS_PROX=TRUE,
    NUM_NCs=if(!is.null(DEFAULT_NUM_NC)) {DEFAULT_NUM_NC} else {c(1, 3, 5, 8)},  
    PROVIDES_NC_MATRIX=FALSE, # (false=vector of gene names for singlegenes)
    NC_RANDOMIZEORDER=TRUE, # whether to randomize order of given NC (yes for singlegenes)
    PROXIMAL_METHOD = 'proximallinear',
    NC_type = 'singlegene',
    Y_TYPE='continuous'
  ),    
  'singlegene_NB'= list(
    RUN=FALSE,
    IS_PROX=TRUE,
    NUM_NCs=if(!is.null(DEFAULT_NUM_NC)) {DEFAULT_NUM_NC} else {c(1, 3, 5, 8, 10, 15)},  
    PROVIDES_NC_MATRIX=FALSE, # (false=vector of gene names for singlegenes)
    NC_RANDOMIZEORDER=TRUE, # whether to randomize order of given NC (yes for singlegenes)
    PROXIMAL_METHOD = 'proximalnegbin',
    NC_type = 'singlegene',
    Y_TYPE='count'
  ), 
  'singlegene_housekeeping'= list(
    RUN=FALSE,
    IS_PROX=TRUE,
    NUM_NCs=if(!is.null(DEFAULT_NUM_NC)) {DEFAULT_NUM_NC} else {c(1, 3, 5, 8, 10, 15)},
    PROVIDES_NC_MATRIX=FALSE,
    NC_RANDOMIZEORDER=TRUE,
    PROXIMAL_METHOD = 'proximallinear',
    NC_type = 'singlegene_housekeeping',
    Y_TYPE='continuous'
  )         
)






# =================== START ========================================================================
print(sprintf("[%s] START: glm and proximal estimates for Jin2025", Sys.time()))

source(sprintf('%s/estimate_proximal/estimate_effects_pci2s.r', util_dir))



dir.create(sprintf('%s/AY/', save_dir), showWarnings = FALSE)
dir.create(sprintf('%s/AY/%s/', save_dir, AY_SETTING_NAME), showWarnings = FALSE)
# dir.create(sprintf('%s/AY/%s/initiallinreg/', save_dir, AY_SETTING_NAME), showWarnings = FALSE)
# dir.create(sprintf('%s/AY/%s/intermediateATEs/', save_dir, AY_SETTING_NAME), showWarnings = FALSE)


intermediateATEs_folder = sprintf('%s/AY/%s/intermediateATEs/', save_dir, AY_SETTING_NAME)
dir.create(intermediateATEs_folder, showWarnings = FALSE)

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
# check gene names
assertthat::assert_that(all(row.names(gene_norm) == row.names(gene_counts))) # same genes and in same order
assertthat::assert_that(all(colnames(gene_norm) == colnames(gene_counts))) # same genes and in same order





# check cellnames 
assertthat::assert_that(all(colnames(gene_norm) == colnames(gene_counts))) # same cells and in same order
assertthat::assert_that(all(colnames(gene_norm) == colnames(grna)))        # same cells and in same order
assertthat::assert_that(all(colnames(gene_norm) == row.names(metadata)))   # same cells and in same order


# ============================================================================== #
# =============  Construct AY tests: ===========================================
# ============================================================================== #
print(sprintf("[%s]    - construct AY tests", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))

# ================ some prep for constructing tests
print(sprintf("[%s]        + some prep selecting possible As and Ys", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
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
  print(sprintf("[%s]        + select AYs under dense setting", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
  # first select all the grna targets
  # then choose other Ys from all_Ys. all_Ys should be arranged by importance rank, so choose the top important/most deviance genes as responses
  selected_Ys = c(true_As, setdiff(all_Ys, true_As))[1:SETTING$NUM_Y]
  
  
  AY = expand.grid(grna = c(true_As, fake_As), gene = selected_Ys, stringsAsFactors = FALSE) |> mutate(type = 'discovery')
  AY[(AY$grna == AY$gene), 'type'] = 'positive'
  AY[grepl(pattern='Non_target', x = AY$grna), 'type'] = 'negative'
} else {
  print(sprintf("[%s]        + select AYs under not dense setting", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
  # =============  \__Construct Negative tests: ===============================================
  AY_negative = data.frame(# grna = sample(setdiff(row.names(grna), 'Non_target'), size = AY_NUM_DISCOVERY, replace=TRUE),
    grna = sample(fake_As, size=SETTING$AY_NUM_NEGATIVE*1.2, replace=TRUE),
    gene = sample( all_Ys, size=SETTING$AY_NUM_NEGATIVE*1.2, replace=TRUE),
    type = 'negative') |> 
    filter(grna != gene) |>        # random choice of grnas and genes may result in pos tests, so sample more
    distinct() |> 
    slice_sample(n=SETTING$AY_NUM_NEGATIVE) # than needed, then remove and take the requested number of tests
  
  
  # ============= \__Construct Discovery tests: ============================================
  AY_discovery = data.frame(grna = sample(true_As, size=SETTING$AY_NUM_DISCOVERY*1.3, replace=TRUE),
                            gene = sample( all_Ys, size=SETTING$AY_NUM_DISCOVERY*1.3, replace=TRUE),
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


saveRDS(list(AY_SETTING_NAME=AY_SETTING_NAME, AY_SETTING=SETTING, ESTIMATOR_SETTINGS=ESTIMATOR_SETTINGS), 
        sprintf('%s/AY/%s/run_settings.rds', save_dir, AY_SETTING_NAME))



# ============================================================================== #
# ============= Prep variables: ================================================
# ============================================================================== #


# ============= Realize Relevant Gene info in memory: ========================== #
required_Y = unique(AY$gene) # only these response genes' expression are required
required_Y = unique(c(required_Y, validsinglegenes)) |> sort() # also for single genes

gene_counts_subset = gene_counts[required_Y, ] |> DelayedArray::realize()
gene_norm_subset   =   gene_norm[required_Y, ] |> DelayedArray::realize()

rm(required_Y)
rm(gene_counts, gene_norm)
invisible(gc())


format(object.size(gene_counts_subset), units='Gb')
format(object.size(gene_norm_subset), units='Gb')


# ============= prep reused variables in estimating fns ========================== #
# prep reused variables in estimating fns
control_cell_idx = which(grna['Non_target', ] == 1) |> as.vector() 
libsize_log      = log(metadata[, 'num_rna_umi'] + 1e-5)

housekeepinggenes = setdiff(validsinglegenes, gene_dev_df |> filter(is_housekeeping & is_normalized_gene & !is_grna_target & !is_transcription_factor))

# for all the proximal methods, this must provide either a matrix of NC values or a vector of gene names as NCs
max_num_ncs = 0
for(n_ in names(ESTIMATOR_SETTINGS)) {max_num_ncs = max(max_num_ncs, ESTIMATOR_SETTINGS[[n_]][['NUM_NCs']])}
NCs = list(
  'PCA_UntrtCells' = as.matrix(readRDS(sprintf('%s/pca/NCloadings_UntrtCells.rds',       save_dir)))[, 1:(2*max_num_ncs)],  # PCA loadings from PCA on Untreated Cells
  'PCA_MixCells'   = as.matrix(readRDS(sprintf('%s/pca/NCloadings_mixTrtUntrtCells.rds', save_dir)))[, 1:(2*max_num_ncs)],  # PCA loadings from PCA on Mix of Cells
  'singlegene'              = validsinglegenes, 
  'singlegene_housekeeping' = housekeepinggenes
)


rm(metadata)
invisible(gc())

# ================ checks
assertthat::assert_that(all(row.names(gene_norm_subset) == row.names(gene_counts_subset))) # same genes and in same order
# check cellnames 
assertthat::assert_that(all(colnames(gene_norm_subset) == colnames(gene_counts_subset))) # same cells and in same order
assertthat::assert_that(is.null(NCs[['PCA_UntrtCells']]) || all(colnames(gene_norm_subset) == row.names(NCs[['PCA_UntrtCells']])))
assertthat::assert_that(is.null(NCs[['PCA_MixCells']])   || all(colnames(gene_norm_subset) == row.names(NCs[['PCA_MixCells']])))


# ============================================================================== #
# =================== Get ATEs by gene (parallel) ==============================
# ============================================================================== #
print(sprintf("[%s]    - Get Estimates by gene (parallel)", Sys.time()))

print(sprintf("[%s]        + nrow(AY)=%s", Sys.time(), nrow(AY)))
print(sprintf("[%s]        + #unique grna=%s", Sys.time(), length(unique(AY$grna))))
print(sprintf("[%s]        + #unique gene=%s", Sys.time(), length(unique(AY$gene))))




# call function here


whichRows = c(1:10)
whichRows = 1:nrow(AY)
# source(sprintf('%s/estimate_proximal/estimate_effects_pci2s.r', util_dir))
# test the bygene version
effects_df = estimate_effects_bygene(estimator_settings=ESTIMATOR_SETTINGS,
                        AY = AY[whichRows, ],
                        NCs = NCs,
                        grna=grna,
                        gene_norm=gene_norm_subset,
                        gene_counts=gene_counts_subset,
                        control_cell_idx=control_cell_idx,
                       libsize_log=libsize_log,
                       covariates=NULL,
                       subsample_control_size=10000, # should already be subsampled
                       save_intermediateATEs=TRUE,
                       skip_if_saved_intermediateATEs=SKIP_IF_SAVED_INTERMEDIATEATES,
                       intermediateATEs_folder=intermediateATEs_folder)



# whichRows = c(1:10)
# source(sprintf('%s/estimate_proximal/estimate_effects_pci2s.r', util_dir))
# test by AY version

# whichRows = 1:nrow(AY)
# effects_df = estimate_effects_byAY(estimator_settings=ESTIMATOR_SETTINGS, 
#                           AY = AY[whichRows, ], 
#                           NCs = NCs,
#                           grna=grna, 
#                           gene_norm=gene_norm_subset, 
#                           gene_counts=gene_counts_subset, 
#                           control_cell_idx=control_cell_idx, 
#                           libsize_log=libsize_log,
#                           covariates=NULL, 
#                           subsample_control_size=10000, # should already be subsampled
#                           save_intermediateATEs=TRUE,
#                           skip_if_saved_intermediateATEs=FALSE,
#                           intermediateATEs_folder=intermediateATEs_folder)




write.csv(effects_df, file = sprintf('%s/AY/%s/effects.csv', save_dir, AY_SETTING_NAME), row.names = FALSE)




# effects_df$time_sec |> sum() / length(whichRows) / 60 # approx this many mins (3.4 mins) per test when trying many methods and method settings
# ggplot(effects_df, aes(x = sprintf('%02.f',numNC), y = pval |> log(), color = NC_type, group = NC_type)) + geom_point() + geom_line() + facet_grid(rows = vars(method), cols = vars(gene)) + coord_cartesian(ylim = c(-10, 0.1))



