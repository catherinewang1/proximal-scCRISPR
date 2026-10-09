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
# args = c('macbook')


assertthat::assert_that(length(args) > 0, msg="must give arg for specifying device eg 'Rscript <filename>.R ubergenno'")
DEVICE = args[1]
source('../PATHS.R') # load in data_dir and save_dir and CODE_DIR, depending on DEVICE value
save_dir = paste0(save_dir, '/2_glutcells')
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
options(future.globals.maxSize= 5000*1024^2) #1st num is MB
# plan(multisession, workers = 24)
plan(multisession, workers = 8)
# plan(sequential)


# AY_SETTING_NAME = 'medium2'
# AY_NUM_NEGATIVE  = 1000
# AY_NUM_DISCOVERY = 1000
# numNCs = c(1, 3, 5, 8, 10, 15)

# AY_SETTING_NAME = 'medium'
# AY_NUM_NEGATIVE  = 1000
# AY_NUM_DISCOVERY = 1000
# numNCs = c(1, 3, 5, 8, 10, 15)


# AY_SETTING_NAME = 'small'
# AY_NUM_NEGATIVE  = 150
# AY_NUM_DISCOVERY = 150
# numNCs = c(1, 3, 5, 8, 10)
# 
# 
# prox_methods = list(run_proximal_linear_pca         = TRUE,
#                     run_proximal_negbin_pca         = TRUE,
#                     run_proximal_linear_singlegenes = TRUE,
#                     run_proximal_negbin_singlegenes = TRUE)




# # Runs with NB proximal 
# AY_SETTING_NAME = 'smallNB'
# AY_NUM_NEGATIVE  = 150
# AY_NUM_DISCOVERY = 150
# numNCs = c(1, 3, 5, 10)


# Runs with NB proximal 
AY_SETTING_NAME = 'mediumNB'
AY_NUM_NEGATIVE  = 250
AY_NUM_DISCOVERY = 250
numNCs = c(1, 3, 5, 8, 10, 15)

prox_methods = list(run_proximal_linear_pca         = TRUE,
                    run_proximal_negbin_pca         = TRUE,
                    run_proximal_linear_singlegenes = TRUE,
                    run_proximal_negbin_singlegenes = TRUE)



GENE_NUMNONZEROCELLS_MIN = 7500
GENE_IMPORTANCERANK_MAX = 2000
GENE_MIN_NONZERO_TRTCELLS = 10000 # filter for genes with at least this number of TREATED cells (not Non_target)

GRNA_MIN_CELLS = 250 # filter for grnas with at least this number of cells


NUM_SUBSAMPLE_CONTROL_CELLS = 10000



set.seed(12345)

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
h5file   = sprintf('%s/jin2025_glut.h5', save_dir)

# =================== Load gene_norm
gene_norm = HDF5Array::HDF5Array(h5file, name='gene_norm/gene_norm')
row.names(gene_norm) = rhdf5::h5read(h5file, name = 'gene_norm/.gene_norm_dimnames/1') |> as.vector()
colnames( gene_norm) = rhdf5::h5read(h5file, name = 'gene_norm/.gene_norm_dimnames/2') |> as.vector()

# =================== Load gene_counts
gene_counts = HDF5Array::HDF5Array(h5file, name='gene_counts/gene_counts')
row.names(gene_counts) = rhdf5::h5read(h5file, name = 'gene_counts/.gene_counts_dimnames/1') |> as.vector()
colnames( gene_counts) = rhdf5::h5read(h5file, name = 'gene_counts/.gene_counts_dimnames/2') |> as.vector()


h5lite::h5_ls(h5file)
rhdf5::h5closeAll()

# =================== Load grna, grnaNT, metadata, gene_dev_df
grna     = h5lite::h5_read(h5file, name = 'grna/grna')
grnaNT   = h5lite::h5_read(h5file, name = 'grnaNT/grnaNT')
metadata = h5lite::h5_read(h5file, name = 'cell_metadata/cell_metadata')

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
assertthat::assert_that(all(colnames(gene_norm) == colnames(grnaNT))) # same cells and in same order
assertthat::assert_that(all(colnames(gene_norm) == row.names(metadata))) # same cells and in same order


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


# all valid treatments (not all grnas, eg if there weren't enough cells)
all_As = row.names(grna)[rowSums(grna) >= GRNA_MIN_CELLS] |> setdiff('Non_target') |> sort()



# =============  \__Construct Negative tests: ===============================================

# For now, just do 2 NT against others to try to make the sample sizes similar?
# ggplot(NULL, aes(x = metadata |> filter(gene_target == 'Non_target') |> pull(guide_call) )) +
#   geom_bar()
# metadata |> filter(gene_target == 'Non_target') |> pull(guide_call) |> hist()

NT_names = grnaNT |> row.names() 





# 'Treatment' = groups of NT perturbations (sample size more similar to actual disc and pos tests)
# choose groups of NT s.t. the n_trtmnt is in the range of 400-1000

# grnaNT_counts = rowSums(grnaNT)


# because of the way the control cells are subsampled, these lose a lot of sample 


# separate NT perturbations into truly control or as Null tests
set.seed(54321)
NT_grnas_for_control = sample(NT_names, size = round(.75*length(NT_names)) )
NT_grnas_for_null    = setdiff(NT_names, NT_grnas_for_control)


prop_pairs    = .50
prop_triplets = .25
AY_NUM_NEGATIVE_TRIPLETS = round(prop_triplets*AY_NUM_NEGATIVE)
AY_NUM_NEGATIVE_PAIRS    = round(prop_pairs*AY_NUM_NEGATIVE)
AY_NUM_NEGATIVE_SINGLE   = AY_NUM_NEGATIVE - AY_NUM_NEGATIVE_TRIPLETS - AY_NUM_NEGATIVE_PAIRS


# 'Treatment' = triplets of NT perturbations (sample size more similar to actual disc and pos tests)
NT_triplets = replicate(AY_NUM_NEGATIVE_TRIPLETS*1.5, sample(NT_grnas_for_null, 3, replace = FALSE)) |> t() |> data.frame() |> distinct()
NT_treatment_name = paste0(NT_triplets[, 1], '&', NT_triplets[, 2], '&', NT_triplets[, 3])[1:AY_NUM_NEGATIVE_TRIPLETS] # would err if the number specified is too large 
NT_outcome_name = sample(all_Ys, size = length(NT_treatment_name), replace = TRUE)
AY_negative_triplets = data.frame(
  grna = NT_treatment_name, 
  gene = NT_outcome_name, 
  type = 'negative')


# 'Treatment' = pairs of NT perturbations (sample size more similar to actual disc and pos tests)
NT_pairs = replicate(AY_NUM_NEGATIVE_PAIRS*1.5, sample(NT_grnas_for_null, 2, replace = FALSE)) |> t() |> data.frame() |> distinct()
NT_treatment_name = paste0(NT_pairs[, 1], '&', NT_pairs[, 2])[1:AY_NUM_NEGATIVE_PAIRS] # would err if the number specified is too large 
NT_outcome_name = sample(all_Ys, size = length(NT_treatment_name), replace = TRUE)
AY_negative_pairs = data.frame(
  grna = NT_treatment_name, 
  gene = NT_outcome_name, 
  type = 'negative')

# 'Treatment' = individual NT perturbations  (less mixup of NT pert)
AY_negative_single = data.frame(
  grna = sample(NT_grnas_for_null, size=AY_NUM_NEGATIVE_SINGLE*1.2, replace=TRUE), 
  gene = sample(all_Ys,            size=AY_NUM_NEGATIVE_SINGLE*1.2, replace = TRUE), 
  type = 'negative') |> distinct() |> head(AY_NUM_NEGATIVE_SINGLE)

AY_negative = rbind(AY_negative_triplets, AY_negative_pairs, AY_negative_single)





rm(NT_names, NT_pairs, NT_treatment_name, NT_outcome_name)

# ============= \__Construct Discovery tests: ============================================

# dim(gene_count)
# gene_numnonzerocells = rowSums(gene_count); hist(gene_numnonzerocells)
AY_discovery = data.frame(# grna = sample(setdiff(row.names(grna), 'Non_target'), size = AY_NUM_DISCOVERY, replace=TRUE),
                          grna = sample(all_As, size=AY_NUM_DISCOVERY*1.2, replace=TRUE),
                          gene = sample(all_Ys, size=AY_NUM_DISCOVERY*1.2, replace=TRUE),
                          type = 'discovery') |> 
               filter(grna != gene) |>        # random choice of grnas and genes may result in pos tests, so sample more
               distinct() |> 
               slice_sample(n=AY_NUM_DISCOVERY) # than needed, then remove and take the requested number of tests

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
# ============= Initial LinReg and NegBin:                  ====================
# ============= Compare across control cell subsample sizes ====================
# ============================================================================== #

# positive tests across different NUM_SUBSAMPLE_CONTROL_CELLS sizes
if('initial_glms_posneg_comparecontrolsubsamplesizes.csv' %in% list.files(sprintf('%s/AY/%s/initiallinreg/', save_dir, AY_SETTING_NAME))) {
  # just skip if already done
  initial_glms_pos_df = read.csv(sprintf('%s/AY/%s/initiallinreg/initial_glms_posneg_comparecontrolsubsamplesizes.csv', save_dir, AY_SETTING_NAME))
  
} else {
  # realize some for this comparison
  required_Y = unique(AY$gene) # only these response genes' expression are required
  gene_counts_subset = gene_counts[required_Y, ] |> DelayedArray::realize()
  gene_norm_subset   =   gene_norm[required_Y, ] |> DelayedArray::realize()
  
  # prep some used values
  combined_grna = rbind(grna, grnaNT)
  control_cell_idx = which(grna['Non_target', ] == 1) |> as.vector() # just the overall control cell idx. use setdiff for negative tests
  libsize_log = log(metadata[, 'num_rna_umi'] + 1e-5)
  
  fit_glms_errorhandling <- function(AY_idx, subsamplesize) {
    tryCatch({fit_glms(AY_idx=AY_idx, 
                       AY=AY, 
                       control_cell_idx=control_cell_idx, 
                       gene_norm=gene_norm_subset, 
                       gene_counts=gene_counts_subset, 
                       grna=combined_grna, 
                       libsize_log=libsize_log, 
                       covariates=NULL, 
                       subsample_control_size=subsamplesize)
    },
    error = function(cond) {
      return(NULL)
    })
  }
  # fit_glms_errorhandling(100, 500)
  # fit_glms_errorhandling(179, 500)
  # fit_glms(AY_idx=179, 
  #          AY=AY, 
  #          control_cell_idx=control_cell_idx, 
  #          gene_norm=gene_norm_subset, 
  #          gene_counts=gene_counts_subset, 
  #          grna=combined_grna, 
  #          libsize_log=libsize_log, 
  #          covariates=NULL, 
  #          subsample_control_size=400)
  
  initial_glms_pos_df = NULL
  for(cur_subsample_size in c(NA, 100000, 75000, 50000, 25000, 15000, 10000, 5000, 1000, 500)) {
    initial_glms_pos_par = future.apply::future_mapply(
          fit_glms_errorhandling, 
          MoreArgs = list(subsamplesize=cur_subsample_size),
          AY_idx = which(AY$type == 'positive' | AY$type == 'negative'), 
          future.globals = TRUE,
          future.seed = 123456, 
          SIMPLIFY = FALSE
        )
    
    for(idx in 1:length(initial_glms_pos_par)) {   initial_glms_pos_df = rbind(initial_glms_pos_df, initial_glms_pos_par[[idx]])   }
  }
  
  
  write.csv(initial_glms_pos_df, sprintf('%s/AY/%s/initiallinreg/initial_glms_posneg_comparecontrolsubsamplesizes.csv', save_dir, AY_SETTING_NAME))

  initial_glms_pos_df = merge(AY, initial_glms_pos_df, by = c('grna', 'gene')) # add test type
  
  
  initial_glms_pos_df[, 'n_cntrl_clean'] = initial_glms_pos_df$n_cntrl
  initial_glms_pos_df[initial_glms_pos_df$n_cntrl >=120000, 'n_cntrl_clean'] = '>130k'
  initial_glms_pos_df$n_cntrl_clean = factor(initial_glms_pos_df$n_cntrl_clean, 
                                             levels = c(initial_glms_pos_df |> filter(n_cntrl <= 120000) |> pull(n_cntrl) |> unique() |> sort() |> as.character(),  '>130k'))
  
  
  
  ggplot(initial_glms_pos_df, aes(x = n_cntrl_clean, y = pval)) + geom_violin() + facet_grid(rows = vars(method, type)) + theme_bw()
  
  
  ggplot(initial_glms_pos_df, aes(y= pval)) + 
    geom_histogram(binwidth = .01) + 
    facet_grid(rows = vars(method, type), cols = vars(n_cntrl_clean), scales = 'free_x') + 
    labs(title = 'Hist of pvals by #controls subsampled') +
    theme_bw()
  
  ggsave(sprintf('%s/AY/%s/initiallinreg/histpvals_by_numsubsample.pdf', save_dir, AY_SETTING_NAME), height=10, width=10)
  
  
  ggplot(initial_glms_pos_df, aes(y= log(pval+1e-10))) + 
    geom_histogram(binwidth = 1) + 
    facet_grid(rows = vars(method, type), cols = vars(n_cntrl_clean), scales = 'free_x') +
    labs(title = 'Hist of log(pvals) by #controls subsampled') +
    theme_bw()
  
  ggsave(sprintf('%s/AY/%s/initiallinreg/histpvals_by_numsubsample_log10.pdf', save_dir, AY_SETTING_NAME), height=10, width=10)

  ggplot(initial_glms_pos_df, aes(y= se)) + 
    geom_histogram(binwidth = .01) + 
    facet_grid(rows = vars(method, type), cols = vars(n_cntrl_clean), scales = 'free_x') +
    labs(title = 'Hist of standard errs by #controls subsampled') +
    theme_bw()
  
  ggsave(sprintf('%s/AY/%s/initiallinreg/histstandarderr_by_numsubsample.pdf', save_dir, AY_SETTING_NAME), height=10, width=10)
  
  
  ggplot(initial_glms_pos_df, aes(y= se)) + 
    geom_histogram(binwidth = .0025) + 
    facet_grid(rows = vars(method, type), cols = vars(n_cntrl_clean), scales = 'free_x') +
    coord_cartesian(ylim = c(0, .2)) +
    labs(title = 'Hist of standard errs by #controls subsampled') +
    theme_bw()
  
  ggsave(sprintf('%s/AY/%s/initiallinreg/histstandarderr_by_numsubsample_zoomed.pdf', save_dir, AY_SETTING_NAME), height=10, width=10)
  
  
  # subsampling for 10k cells seem ok (still have power, standard errors are pretty close to full sample)
}



# ============================================================================== #
# ============= Realize Relevant Gene info in memory: ==========================
# ============================================================================== #
required_Y = unique(AY$gene) # only these response genes' expression are required
required_Y = unique(c(required_Y, validsinglegenes)) |> sort() # also for single genes

set.seed(12345678)

# separate NTs for control and for null tests
subsample_ctrl_cells = sample(which(colSums(grnaNT[NT_grnas_for_control, ])==1)  |> names(), round(NUM_SUBSAMPLE_CONTROL_CELLS*1.1), replace=FALSE)# take a few more than necessary
subsample_cells = c(which(grna['Non_target', ]==0) |> names(), 
                    which(colSums(grnaNT[NT_grnas_for_null, ])==1)  |> names(),
                    subsample_ctrl_cells) # all treated cells + null pert cells + the prev sampled cells

# if NTs were mixed up
# subsample_ctrl_cells = sample(which(grna['Non_target', ]==1) |> names(), round(NUM_SUBSAMPLE_CONTROL_CELLS*2), replace=FALSE)# take a few more than necessary
# subsample_cells = c(which(grna['Non_target', ]==0) |> names(), subsample_ctrl_cells) # all treated cells + the prev sampled cells

gene_counts_subset = gene_counts[required_Y, subsample_cells] |> DelayedArray::realize()
gene_norm_subset   =   gene_norm[required_Y, subsample_cells] |> DelayedArray::realize()

rm(subsample_ctrl_cells)
rm(gene_counts, gene_norm)
gc()


format(object.size(gene_counts_subset), units='Gb')
format(object.size(gene_norm_subset), units='Gb')

# ============================================================================== #
# ============= Prep variables: ================================================
# ============================================================================== #

# prep reused variables in estimating fns
combined_grna    = rbind(grna[, subsample_cells], grnaNT[, subsample_cells])
control_cell_idx = which(combined_grna['Non_target', ] == 1) |> as.vector() 
libsize_log      = log(metadata[subsample_cells, 'num_rna_umi'] + 1e-5)

# prep reused variables in estimating fns (prox)
NCs_pca = readRDS(sprintf('%s/pca/NCloadings.rds', save_dir)) |> as.matrix() # PCA loadings

rm(grna, grnaNT, metadata)
gc()

# ================ checks
assertthat::assert_that(all(row.names(gene_norm_subset) == row.names(gene_counts_subset))) # same genes and in same order


# check cellnames 
assertthat::assert_that(all(colnames(gene_norm_subset) == colnames(gene_counts_subset))) # same cells and in same order
assertthat::assert_that(all(colnames(gene_norm_subset) == colnames(combined_grna))) # same cells and in same order



# ============================================================================== #
# ============= Initial LinReg and NegBin:                    ==================
# ============= after choosing a control cell subsample size  ==================
# ============================================================================== #

if('initial_glms.csv' %in% list.files(sprintf('%s/AY/%s/initiallinreg/', save_dir, AY_SETTING_NAME))) {
  print(sprintf("[%s]        - Loading glm (one fit per gene) on all pairs", Sys.time()))
  initial_glms = read.csv(sprintf('%s/AY/%s/initiallinreg/initial_glms.csv', save_dir, AY_SETTING_NAME))
} else {
  print(sprintf("[%s]        - Fitting glm (one fit per gene) on all pairs", Sys.time()))
  
  fit_glms_errorhandling <- function(AY_idx) {
    tryCatch({fit_glms(AY_idx=AY_idx, 
                       AY=AY, 
                       control_cell_idx=control_cell_idx, 
                       gene_norm=gene_norm_subset, 
                       gene_counts=gene_counts_subset, 
                       grna=combined_grna, 
                       libsize_log=libsize_log, 
                       covariates=NULL, 
                       subsample_control_size=NUM_SUBSAMPLE_CONTROL_CELLS)
    },
    error = function(cond) {
      return(NULL)
    })
  }
  
  
  t0 = Sys.time()
  initial_glms_par = future.apply::future_mapply(
    fit_glms_errorhandling,
    AY_idx = 1:nrow(AY), 
    future.globals = TRUE,
    future.seed = 123456, 
    SIMPLIFY = FALSE)
  t1 = Sys.time()
  print(sprintf("[%s]        - %2.2f mins", Sys.time(), difftime(t1, t0, units = "mins")))
  
  initial_glms = NULL
  for(idx in 1:length(initial_glms_par)) {   initial_glms = rbind(initial_glms, initial_glms_par[[idx]])   }
  write.csv(initial_glms, sprintf('%s/AY/%s/initiallinreg/initial_glms.csv', save_dir, AY_SETTING_NAME))
  rm(initial_glms_par, t0, t1, fit_glms_errorhandling)
}





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
# =================== Get ATEs (parallel) ======================================
# ============================================================================== #
# 
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
                                      numNCs=numNCs, 
                                      # data
                                      control_cell_idx=control_cell_idx, 
                                      gene_norm=gene_norm_subset, 
                                      gene_counts=gene_counts_subset, 
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
                                      run_proximal_linear_pca         = prox_methods$run_proximal_linear_pca,
                                      run_proximal_negbin_pca         = prox_methods$run_proximal_negbin_pca,
                                      run_proximal_linear_singlegenes = prox_methods$run_proximal_linear_singlegenes,
                                      run_proximal_negbin_singlegenes = prox_methods$run_proximal_negbin_singlegenes,
                                      subsample_control_size = NUM_SUBSAMPLE_CONTROL_CELLS)
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
#                  grna=combined_grna,
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
    # AY test is whichROWS[whichROWS_idx]  = AY_idx!!! 
    effects_df = dplyr::bind_rows(effects_df, cur_effects_par)
  }
  rm(cur_effects_par)
}
write.csv(x = effects_df, file = sprintf('%s/AY/%s/effects.csv', save_dir, AY_SETTING_NAME), row.names = FALSE)



# ===========================================================================
print(sprintf("[%s] END", Sys.time()))








