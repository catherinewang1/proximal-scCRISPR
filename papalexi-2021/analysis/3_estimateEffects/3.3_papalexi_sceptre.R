# ---------------------------------------------------------------------------- #
#           Other Estimator: SCEPTRE w measured confounders                    #
# Oracle Estimators:                                                           #
# - simple lm(Y ~ A + U)                                                       #
# - SCEPTRE (using raw counts)                                                 #
# SCEPTRE is the ideal estimator, but SCEPTRE uses raw counts and assumes a    #
# non-identity link function (for NB or Pois?). So, we should really compare   #
# these linear proximal methods                                                #
# ---------------------------------------------------------------------------- #
args = commandArgs(trailingOnly = TRUE)
args = c('ubergenno', 'notneeded', 'A')

require(assertthat) # for some assert statements
library(Matrix)
library(dplyr)
library(tibble)
library(ggplot2)    # plotting
library(cowplot)
library(sceptre)
# library(sceptredata)

# 
# library(future.apply)
# options(future.globals.maxSize= 850*1024^2) #1st num is MB
# plan(multisession, workers = 8)
# plan(sequential)

theme_set(theme_cowplot() +
            theme(plot.title = element_text(hjust = .5),
                  plot.subtitle = element_text(hjust = .5),
                  strip.background = element_rect(color = 'black', fill = 'white')))

DEVICE = args[1]
source('../PATHS.R') # load in data_dir and save_dir and CODE_DIR, depending on DEVICE value

assertthat::assert_that(!is.null(data_dir), msg='first arg must be: laptop, desktop, or ubergenno')

assertthat::assert_that(length(args) > 1, msg="must give arg for specifying num imp genes 'Rscript <filename>.R ubergenno 4000'")
NUM_IMPORTANTGENES = as.integer(args[2]) # should be max importance from setting... TODO: change to remove this input. extract this value from saved setting

assertthat::assert_that(length(args) > 2, msg="must give arg for specifying chosen AYZW name 'Rscript <filename>.R ubergenno C'")
AYZW_setting_name = args[3]

SCEPTRE_savepath = sprintf('%s/AY/%s/SCEPTRE/', save_dir, AYZW_setting_name)
dir.create(SCEPTRE_savepath, recursive = TRUE, showWarnings = FALSE)


# =================== Start ====================================================
print(sprintf("[%s] START: Oracle (SCEPTRE)", Sys.time()))


# source(sprintf('%s/CBEstAll.R', util_dir)) # for functions to estimate here

# load chosen AYZW names
AY   = read.csv(sprintf('%s/AY/%s/AY.csv', save_dir, AYZW_setting_name))
# AYZW = readRDS(sprintf('%s/AY/%s/AYZW.rds', save_dir, AYZW_setting_name))

# load gene importance info
# imp_gene_names = readRDS(sprintf('%s/important_genes_name.rds', save_dir))
# imp_gene_idx   = readRDS(sprintf('%s/important_genes_idx.rds',  save_dir))
# imp_gene = data.frame(gene     = imp_gene_names,
#                       gene_idx = imp_gene_idx,
#                       gene_imp_rank = 1:length(imp_gene_names))

# gene_importance = read.csv(sprintf('%s/gene_deviance_gene_norm.csv', save_dir)) |> 
#   dplyr::select(gene_name, gene_idx, importance_rank, gene_norm_idx)


# create gene and grna ondisc managers
gene_odm <- ondisc::read_odm(odm_fp      = paste0(data_dir, "/papalexi-2021/processed/gene/expression_matrix.odm"),
                             metadata_fp = paste0(data_dir, "/papalexi-2021/processed/gene/metadata.rds"))
grna_odm <- ondisc::read_odm(odm_fp      = paste0(data_dir, "/papalexi-2021/processed/grna_assignment/assignment_matrix.odm"),
                             metadata_fp = paste0(data_dir, "/papalexi-2021/processed/grna_assignment/metadata.rds"))

# load grna assignments (load all into memory)
grna = grna_odm[[,1:ncol(grna_odm)]] # |> as.matrix() # ~110 x 20729 = #grnas x #cells
grna_rownames = ondisc::get_feature_ids(grna_odm)
rownames(grna) = grna_rownames

# # load normalized gene exp
# h5file      = paste0(save_dir, "/gene.h5"); print(h5file)
# reading_hd5file  = rhdf5::H5Fopen(name = h5file)
# readin_gene_norm = reading_hd5file&'gene_norm'
# gene_norm = readin_gene_norm[1:NUM_IMPORTANTGENES, ] # dim = 4000 x 20729 = #important x #cells
# rownames(gene_norm) = imp_gene_names[1:1:NUM_IMPORTANTGENES]
# rhdf5::h5closeAll()
# invisible(gc(verbose=FALSE))


# original gene expr counts
# gene = gene_odm[[,1:ncol(gene_odm)]]

# =========== Create sceptre object ============================================



# ==============================================================================
# 1.2 Import data from a collection of R objects
# https://timothy-barry.github.io/sceptre-book/import-data.html#import-data-from-a-collection-of-r-objects
# ==============================================================================

# data(highmoi_example_data)
# # response matrix
# response_matrix <- highmoi_example_data$response_matrix
# # grna matrix
# grna_matrix <- highmoi_example_data$grna_matrix
# # batch information
# extra_covariates <- highmoi_example_data$extra_covariates
# # response names
# response_names <- highmoi_example_data$gene_names
# # gRNA target data frame
# grna_target_data_frame <- grna_target_data_frame_highmoi


# cell confounders
cell_covariates = gene_odm |> ondisc::get_cell_covariates()
# rownames(gene_norm) = imp_gene_names[1:nrow(gene_norm)]

# response matrix
# response_matrix <- gene_norm
response_matrix <- gene_odm[[,1:ncol(gene_odm)]]
colnames(response_matrix) <- gene_odm |> ondisc::get_cell_barcodes()
rownames(response_matrix) <- gene_odm |> ondisc::get_feature_ids()
# grna matrix
grna_matrix <- grna
# batch information
extra_covariates <- cell_covariates
# change batch info, bc lane determines rep_1, so redundant info
# table(extra_covariates |> dplyr::select(lane, bio_rep))
extra_covariates <- cell_covariates |> 
                      dplyr::mutate(lane_bio_rep = paste0(lane, '_', bio_rep)) |>
                      dplyr::select(-lane, -bio_rep) |>
                      dplyr::select(-n_nonzero, -n_umis)
         


# response names
# gRNA target data frame
grna_odm_feature_covariates = ondisc::get_feature_covariates(grna_odm)
grna_target_data_frame = data.frame(grna_id = rownames(grna_odm_feature_covariates),
                                    grna_target = grna_odm_feature_covariates$target)



sceptre_object <- import_data(
  response_matrix        = response_matrix,
  grna_matrix            = grna_matrix,
  grna_target_data_frame = grna_target_data_frame,
  moi                    = 'low',
  extra_covariates       = extra_covariates
)

sceptre_object


# ==============================================================================
# 2. Set analysis parameters
# https://timothy-barry.github.io/sceptre-book/sceptre.html#sec-whole_game_set_analysis_parameters
# ==============================================================================
positive_control_pairs <- construct_positive_control_pairs(sceptre_object)
head(positive_control_pairs)
# 
# discovery_pairs <- construct_cis_pairs(
#   sceptre_object = sceptre_object,
#   positive_control_pairs = positive_control_pairs,
#   distance_threshold = 5e6
# )


# Automatically gather possible tests, then add our selected tests (from AY)
discovery_pairs_auto <- construct_trans_pairs(
  sceptre_object = sceptre_object,
  positive_control_pairs = positive_control_pairs,
  pairs_to_exclude = "pc_pairs"
)



# Add AY (it should already be automatically added but just in case)
# format AY to discovery_pairs (grna_target vs actual grna (A). Y = response_id)

# > head(discovery_pairs)
# grna_target response_id
# <char>      <char>
#   1:        CUL3  AL627309.1
#   2:       CMTM6  AL627309.1
# AY$Y %in% discovery_pairs$response_id |> mean()




discovery_pairs_AY = 
  merge(AY |> filter(type != 'positive'), 
        ondisc::get_feature_covariates(grna_odm) |> tibble::rownames_to_column(var = "A"), 
        all.x = TRUE, all.y = FALSE, by = "A") |> 
  dplyr::mutate(grna_target = target,
                response_id = Y) |> 
  dplyr::select("grna_target", "response_id")




# discovery_pairs_auto = discovery_pairs_auto[1:1000, ] # Testing: just do some

discovery_pairs = rbind(discovery_pairs_AY, discovery_pairs_auto) |> dplyr::distinct()



# dim(discovery_pairs_auto)
# dim(discovery_pairs_AY)
# dim(discovery_pairs)
# 
# dim(discovery_pairs |> distinct())
# 
# 
# ondisc::get_feature_covariates(grna_odm) |> tibble::rownames_to_column(var = "A")  |> head()




# side <- "left" # left for testing decrease in expr
side = "both" # change to "both" bc other methods test both!


sceptre_object <- set_analysis_parameters(
  sceptre_object = sceptre_object,
  discovery_pairs = discovery_pairs,
  positive_control_pairs = positive_control_pairs,
  side = side, 
  grna_integration_strategy = 'singleton',
  formula_object = 'default'
)
print(sceptre_object) # output suppressed for brevity


sceptre_object@covariate_data_frame |> anyNA()
sceptre_object@covariate_data_frame |> head()

# ==============================================================================
# 3. Assign gRNAs to cells 
# https://timothy-barry.github.io/sceptre-book/sceptre.html#sec-sceptre_assign_grnas
# ==============================================================================

sceptre_object <- assign_grnas(sceptre_object = sceptre_object, parallel = FALSE,
                               method = 'maximum', 
                               min_grna_n_umis_threshold = 1, 
                               umi_fraction_threshold = .8)
print(sceptre_object) # output suppressed for brevity
plot(sceptre_object) 



# === === === === === === === === === === === === === === === === === === === ===
# # Debugging 
# # Histograms of the gRNA count distributions- only 0/1 bc low-MOI
# plot_grna_count_distributions(sceptre_object)
# 
# 
# # don't know how to set this properly, just setting this manually. default min thresh is 5, but impossible for lowMOI (receives just 1), so, set to just 1?
# # ?? every time assign_grnas changes the thres to 5
# sceptre_object@grna_assignment_hyperparameters$min_grna_n_umis_threshold = .9
# 
# sceptre_object <- assign_grnas(sceptre_object = sceptre_object, parallel = FALSE,
#                                method = 'maximum', 
#                                min_grna_n_umis_threshold = 1, 
#                                umi_fraction_threshold = .8)
# # sceptre_object <- assign_grnas(sceptre_object = sceptre_object, parallel = FALSE,
#                                # method = 'default')
# # sceptre_object <- assign_grnas(sceptre_object = sceptre_object, parallel = FALSE,
# #                                method = 'thresholding')
# print(sceptre_object) # output suppressed for brevity
# plot(sceptre_object) # errs
# 
# # wait, why assign grnas? if already have known assignment?
# 
# sceptre_object
# 
# 
# 
# plots = sceptre::plot_assign_grnas(sceptre_object, return_indiv_plots = TRUE)
# names(plots)
# length(plots)
# # https://github.com/Katsevich-Lab/sceptre/blob/main/R/plotting_functions.R
# # lines 165-279
# plots[[1]] # errs
# plots[[2]] # works
# plots[[3]] # blank, looking at source code, yes this is intended for 'maximum' method
# 
# 
# # ISSUE: is that the fn throws away cells that received 0 or too many mixed grnas (relevant in high MOI)
# # based on sceptre_object@grna_assignment_hyperparameters$umi_fraction_threshold (determines if too mixed)
# # and 
# # based on sceptre_object@grna_assignment_hyperparameters$min_grna_n_umis_threshold hyperparameter (determines of zero) which was set to 5
# # so removed all bc each cell received at most 1.
# 
# 
# # -=-=-=-=-=-=-=-=-=- copying code for first plot, called 'plot_a' -=-=-=-=-=-=-=-=-=-
# # https://github.com/Katsevich-Lab/sceptre/blob/8abf87aaa390dcf5188fbf5bbe7e700b05e53697/R/s4_helpers.R#L224
# load_row <- function(mat, id) {
#   if (methods::is(mat, "odm")) {
#     mat[id, ]
#   } else if (methods::is(mat, "dgRMatrix")) {
#     load_csr_row(j = mat@j, p = mat@p, x = mat@x, row_idx = which(id == rownames(mat)), n_cells = ncol(mat))
#   }
# }
# 
# 
# 
# n_grnas_to_plot = 3L; grnas_to_plot = NULL; point_size = 0.9; transparency = 0.8; return_indiv_plots = FALSE 
# n_points_to_plot_per_umi <- 1000
# n_grnas_to_plot_panel_b <- 1000
# if (!sceptre_object@functs_called["assign_grnas"]) {
#   stop("This `sceptre_object` has not yet had `assign_grnas` called on it.")
# }
# init_assignments <- sceptre_object@initial_grna_assignment_list
# grna_matrix <- get_grna_matrix(sceptre_object) |> set_matrix_accessibility(make_row_accessible = TRUE) # is 'dgRMatrix', looking at source code, nothing happens, line 17, https://github.com/Katsevich-Lab/sceptre/blob/8abf87aaa390dcf5188fbf5bbe7e700b05e53697/R/preprocessing_functions.R#L1
# grna_ids <- names(init_assignments)
# assigned <- vapply(init_assignments, length, FUN.VALUE = integer(1)) >= 1
# grna_ids <- grna_ids[assigned]
# # sample grnas to plot
# if (is.null(grnas_to_plot)) {
#   grnas_to_plot <- sample(x = grna_ids, size = min(nrow(grna_matrix), n_grnas_to_plot), replace = FALSE)
# } else {
#   if (!(all(grnas_to_plot %in% grna_ids))) stop("gRNA IDs must be a subset of the rownames of the gRNA matrix.")
# }
# to_plot_a <- lapply(X = grnas_to_plot, function(curr_grna_id) {
#   assignment <- cells_w_zero_or_twoplus_grnas <- logical(length = ncol(grna_matrix)) # logical vecs w/ one entry per cell
#   assignment[init_assignments[[curr_grna_id]]] <- TRUE # for this grna, `assignment` indicates which cells got this grna initially
#   cells_w_zero_or_twoplus_grnas[sceptre_object@cells_w_zero_or_twoplus_grnas] <- TRUE # indicates which cells have 0/2+ grnas
#   # g <- load_row(grna_matrix, curr_grna_id)
#   g <- grna_matrix[curr_grna_id, ]
#   df <- data.frame(
#     g = g,
#     assignment = ifelse(assignment, "pert", "unpert") |> factor(),
#     grna_id = curr_grna_id |> factor(),
#     cells_w_zero_or_twoplus_grnas = cells_w_zero_or_twoplus_grnas
#   )
#   # if assignment method maximum, remove cells with 0/2+ grnas
#   if (sceptre_object@grna_assignment_method == "maximum") df <- df |> dplyr::filter(!cells_w_zero_or_twoplus_grnas)
#   return(df)
# }) |> data.table::rbindlist()
# 
# # downsample the unperturbed cells
# to_plot_a <- to_plot_a |>
#   dplyr::group_by(g) |>
#   dplyr::sample_n(size = min(n_points_to_plot_per_umi, dplyr::n())) |>
#   dplyr::ungroup()
# 
# # plot a
# p_a <- ggplot2::ggplot(data = to_plot_a, mapping = ggplot2::aes(x = assignment, y = g, col = assignment)) +
#   ggplot2::geom_jitter(alpha = transparency, size = point_size) +
#   ggplot2::facet_wrap(nrow = 1, facets = grna_id ~ ., scales = "free_y") +
#   ggplot2::scale_y_continuous(trans = "log1p", breaks = c(0, 1, 5, 15, 50, 200, 1000, 5000, 10000)) +
#   ggplot2::xlab("Assignment") +
#   ggplot2::ylab("gRNA count") +
#   get_my_theme() +
#   ggplot2::theme(legend.position = "none") +
#   ggplot2::theme(axis.title.x = ggplot2::element_blank()) +
#   ggplot2::scale_color_manual(values = c("firebrick1", "darkorchid1"))
# 
# 
# # cells_w_multiple_grnas <- which(sceptre_object@import_grna_assignment_info$max_grna_frac_umis <= umi_fraction_threshold)
# # cells_w_zero_grnas <- which(grna_n_umis < min_grna_n_umis_threshold)
# cells_w_multiple_grnas <- which(sceptre_object@import_grna_assignment_info$max_grna_frac_umis <= sceptre_object@grna_assignment_hyperparameters$umi_fraction_threshold)
# cells_w_zero_grnas <- which(sceptre_object@covariate_data_frame$grna_n_nonzero < sceptre_object@grna_assignment_hyperparameters$min_grna_n_umis_threshold)
# 
# 
# 
# sceptre_object@covariate_data_frame$grna_n_nonzero
# sceptre_object@grna_assignment_hyperparameters$min_grna_n_umis_threshold
#  
# # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# === === === === === === === === === === === === === === === === === === === ===

# ==============================================================================
# 4. Run quality control 
# https://timothy-barry.github.io/sceptre-book/sceptre.html#sec-sceptre_qc
# ==============================================================================

sceptre_object <- sceptre::run_qc(sceptre_object, p_mito_threshold = 0.075)
print(sceptre_object) # output suppressed for brevity
plot(sceptre_object) 

# ==============================================================================
# 5. Run calibration check
# https://timothy-barry.github.io/sceptre-book/sceptre.html#sec-sceptre_calibration_check
# ==============================================================================
sceptre_object <- run_calibration_check(sceptre_object, parallel = FALSE)
print(sceptre_object) # output suppressed for brevity

plot(sceptre_object)

# ==============================================================================
# 6. Run power check
# https://timothy-barry.github.io/sceptre-book/sceptre.html#sec-sceptre_run_power_check
# ==============================================================================
sceptre_object <- run_power_check(sceptre_object, parallel = FALSE)
print(sceptre_object) # output suppressed for brevity

plot(sceptre_object)

# ==============================================================================
# 7. Run discovery analysis
# https://timothy-barry.github.io/sceptre-book/sceptre.html#sec-sceptre_run_discovery_analysis
# ==============================================================================
t0 = Sys.time()
sceptre_object <- run_discovery_analysis(sceptre_object, parallel = FALSE)
t1 = Sys.time()
print(sceptre_object) # output suppressed for brevity


plot(sceptre_object)

# ==============================================================================
# 8. Write outputs to directory
# https://timothy-barry.github.io/sceptre-book/sceptre.html#sec-sceptre_write_outputs_to_directory
# ==============================================================================
write_outputs_to_directory(
  sceptre_object = sceptre_object, 
  directory = SCEPTRE_savepath

)

SCEPTRE_time_benchmark = 
  list(time_sec = difftime(t1, t0, units = 'secs'), 
       num_tests = nrow(sceptre_object@discovery_pairs),
       avg_time_per_test = difftime(t1, t0, units = 'secs')/nrow(sceptre_object@discovery_pairs),
       info = 'SCEPTRE on papalexi timing for discovery pair analysis')

saveRDS(SCEPTRE_time_benchmark, file = sprintf('%s/SCEPTRE_time_benchmark.rds', SCEPTRE_savepath))

# sceptre_object@discovery_result
# test_load_results = readRDS(sprintf('%s/results_run_discovery_analysis.rds', SCEPTRE_savepath))
# test_load_results




