# ---------------------------------------------------------------------------- #
#   Check SCEPTRE w/ measured confounders 
#   by artificially inserting actual confounding effects by and performing
#   SCEPTRE with and without these fake confounders (while still including the
#   actual measured covariates)
# 
# 1. Creating 1 continuous confounder:
#      FakeU1 = rnorm(0, 1)
#   For notation, keep the previous measured covariates as U1, U2, ...
# 2. Insert confounding effects 
#   confounding effect sizes, for perturbation p and gene g
#     betaA[p]      = runif(-1, 1)
#     betaFakeU1[g] = runif(-1, 1)
#   change perturbation assignment and gene expression
#     cell i, perturbation p, gene g
#   grna_odm_Fake[[perturbation p, cell i]] = Bernoulli(grna_odm[[p, i]] + betaA[p])
#   gene_odm_Fake[[gene g, cell i]] = gene_odm_Fake[[gene g, cell i]] + NegBin(betaFakeU1[g] + offset(log(libsize?)))
# 3. Then, we perform X runs of SCEPTRE
#    1) SCEPTRE with U1 + U2 + ...
#    2) SCEPTRE with U1 + U2 + ... + FakeU1 + FakeU2
#
# 
# ---------
# 
# 1. Creating 2 continuous confounders:
#      FakeU1 = rnorm(0, 1)
#      FakeU2 = rnorm(2, 2)
#      FakeU1 \indep FakeU2
#   For notation, keep the previous measured covariates as U1, U2, ...
# 2. Insert confounding effects 
#   confounding effect sizes
#     betaA1 = runif(0, 1)
#     betaFakeU1 = runif(-1, 1)
#     betaFakeU2 = runif(0, 1)
#   change perturbation assignment and gene expression
#     cell i, perturbation p, gene g
#   grna_odm_Fake[[perturbation p, cell i]] = Bernoulli(grna_odm + betaA[p])
# 3. Then, we perform X runs of SCEPTRE
#    1) SCEPTRE with U1 + U2 + ...
#    2) SCEPTRE with U1 + U2 + ... + FakeU1 + FakeU2
# ---------
# ---------------------------------------------------------------------------- #
args = commandArgs(trailingOnly = TRUE)
# args = c('ubergenno')
# args = c('macbook', 'singleton')



RUN_PARALLEL = TRUE

require(assertthat) # for some assert statements
library(Matrix)
library(dplyr)
library(tibble)
library(ggplot2)    # plotting
library(cowplot)
library(sceptre)


theme_set(theme_cowplot() +
            theme(plot.title = element_text(hjust = .5),
                  plot.subtitle = element_text(hjust = .5),
                  strip.background = element_rect(color = 'black', fill = 'white')))

DEVICE = args[1]
source('../PATHS.R') # load in data_dir and save_dir and CODE_DIR, depending on DEVICE value

assertthat::assert_that(!is.null(data_dir), msg='first arg must be: laptop, desktop, or ubergenno')

assertthat::assert_that(length(args) >= 2, msg='second arg for grna_integration_strategy: union or singleton')
GRNA_INTEGRATION_STRATEGY = args[2]


# =================================================================================================#
# =================== FUNCTIONS =======================================================================
# =================================================================================================#



# =================================================================================================#
# =================== + Add Fake Confounding ========================================================
# =================================================================================================#

#' construct fake data with 1 new confounder by updating a given dataset
#' @param confounding_size_A (numeric) size of confounding coefficient on 
#'  treatment assignments. coefficients are drawn from
#'   unif(-1, 1)*confounding_size_A
#' @param 
#' @example
#' chosen_grnas = c('NTg1', 'NTg2', 'MARCH8g4', 'MYCg2')
#' chosen_grnas_idx = which(colSums(grna_odm[[chosen_grnas, ]])>0)
#' new_data = construct_fake_data(cur_grna= grna_odm[[chosen_grnas, chosen_grnas_idx]],
#'                                cur_gene=gene_odm[[1:100, chosen_grnas_idx]],
#'                                libsize=gene_odm@cell_covariates$n_umis[chosen_grnas_idx])
construct_fake_data <- function(cur_grna, cur_gene, libsize=NULL, 
                                confounding_size_A=1, confounding_size_Y=1) {
  # # for testing on a subset
  # chosen_grnas = c('NTg1', 'NTg2', 'MARCH8g4', 'MYCg2')
  # chosen_grnas_idx = which(colSums(grna_odm[[chosen_grnas, ]])>0)
  # cur_grna = grna_odm[[chosen_grnas, chosen_grnas_idx]]
  # cur_gene = gene_odm[[1:100, chosen_grnas_idx]]
  # 
  # libsize = gene_odm@cell_covariates$n_umis[chosen_grnas_idx]
  # construct_fake_data(cur_grna= grna_odm[[chosen_grnas, chosen_grnas_idx]],
  #                     cur_gene=gene_odm[[1:100, chosen_grnas_idx]],
  #                     libsize=gene_odm@cell_covariates$n_umis[chosen_grnas_idx])
  num_cells = ncol(cur_grna)
  num_perturbations = nrow(cur_grna)
  num_genes = nrow(cur_gene)
  
  if(is.null(libsize)) {libsize=rep(.001, num_cells)} # if libsize not specified, basically set to 0
  
  # for each cell i, the current assignment is grna_odm[[, i]].
  # if current assignment is NT, then with some random chance, choose 
  
  betaA      = runif(n=num_perturbations, min=-1, max=1) * confounding_size_A
  betaFakeU1 = runif(n=num_genes,         min=-1, max=1) * confounding_size_Y
  fakeU1 = rnorm(n=num_cells, mean=0, sd=2)
  
  new_grna = matrix(NA, nrow = num_perturbations, ncol = num_cells)
  # new_gene = matrix(NA, nrow = num_genes,         ncol = num_cells)
  for(i in 1:num_cells) {
    # if(i%%100==0) {print(i)}
    # potentially change treatment assignment
    pert_prob_mean = as.numeric(cur_grna[, i]) + betaA*fakeU1[i]    # like \mu = X \beta
    pert_prob      = exp(pert_prob_mean) / sum(exp(pert_prob_mean)) # exp(\mu_p)/sum(exp(\mu_p)) 
    new_grna[, i]  = rmultinom(n=1, size = 1, prob = pert_prob)
    
    # change gene expression 
    # new_gene[, i]  = cur_gene[, i] + MASS::rnegbin(n=num_genes, mu = exp(betaFakeU1*fakeU1[i]), theta = 4) # theta here is just set, but may not match w data
    # cur_gene[, i] + MASS::rnegbin(n=num_genes, mu = exp(betaFakeU1) + log(libsize[i]), theta = 4) # theta here is just set, but may not match w data
    # cur_gene[, i] + MASS::rnegbin(n=1, mu = exp(betaFakeU1 + log(libsize[i])), theta = 4) # theta here is just set, but may not match w data
    
  }
  row.names(new_grna) = row.names(cur_grna)
  
  # #pert x #cells
  # pert_mean_matrix = as.matrix(betaA) %*% t(as.matrix(fakeU1))
  
  
  # change gene expression here (faster)
  # #genes x #cells 
  # how to include the library size?
  # gene_expr_mean = exp(as.matrix(betaFakeU1) %*% t(as.matrix(fakeU1)))
  # gene_expr_mean = exp(as.matrix(betaFakeU1) %*% t(as.matrix(fakeU1)))+ rep(1, num_genes) %*% t(log(libsize))
  gene_expr_mean = exp(as.matrix(betaFakeU1) %*% t(as.matrix(fakeU1))) + rep(1, num_genes) %*% t(libsize)
  
  
  gene_expr_mean_vec  = as.vector(gene_expr_mean) # columnwise: x[1,1], x[2,1], x[3,1], ...
  gene_expr_to_add_vec = rnbinom(n=length(gene_expr_mean_vec), mu = gene_expr_mean_vec, size=4)
  gene_expr_to_add_mat = matrix(gene_expr_to_add_vec, nrow=nrow(gene_expr_mean), ncol=ncol(gene_expr_mean), byrow = FALSE) # same direction as as.vector (columnwise)
  new_gene = as.matrix(cur_gene + gene_expr_to_add_mat)
  
  rm(gene_expr_mean, gene_expr_mean_vec, gene_expr_to_add_vec, gene_expr_to_add_mat)
  return(list(grna=new_grna,
              gene=new_gene, 
              fakeU1=fakeU1))
}



# =================================================================================================#
# =================== + SCEPTRE in a function =======================================================
# =================================================================================================#


run_sceptre <- function(response_matrix, 
                        grna_matrix, 
                        extra_covariates,
                        grna_target_data_frame,
                        positive_control_pairs,
                        discovery_pairs,
                        write_path=NULL) {
  
  # =================================================================================================#
  # =================== * 1.2 Import data from a collection of R objects ============================
  # https://timothy-barry.github.io/sceptre-book/import-data.html#import-data-from-a-collection-of-r-objects
  print(sprintf('[%s]: SCEPTRE 1. Import data from a collection of R objects', format(Sys.time(), digits = 0)))
  
  if(is.null(extra_covariates)) {
    print(sprintf('[%s]:         add covariates', format(Sys.time(), digits = 0)))
    sceptre_object <- import_data(
      response_matrix        = response_matrix,
      grna_matrix            = grna_matrix,
      grna_target_data_frame = grna_target_data_frame,
      moi                    = 'low'
    )
  } else { # don't include extra_covariates
    print(sprintf('[%s]:         no covariates', format(Sys.time(), digits = 0)))
    sceptre_object <- import_data(
      response_matrix        = response_matrix,
      grna_matrix            = grna_matrix,
      grna_target_data_frame = grna_target_data_frame,
      moi                    = 'low',
      extra_covariates       = extra_covariates
    )
  }
  
  
  
  # =================================================================================================#
  # =================== * 2. Set analysis parameters ================================================
  # https://timothy-barry.github.io/sceptre-book/sceptre.html#sec-whole_game_set_analysis_parameters
  print(sprintf('[%s]: SCEPTRE 2. Set analysis parameters', format(Sys.time(), digits = 0)))
  # here, use prev saved AY tests from standard run. So code to run that was removed here in this file
  
  # side <- "left" # left for testing decrease in expr
  side = "both" # change to "both" bc other methods test both!
  
  
  sceptre_object <- set_analysis_parameters(
    sceptre_object = sceptre_object,
    discovery_pairs = discovery_pairs, # remove? just don't have any discovery pairs? no include
    positive_control_pairs = positive_control_pairs,
    side = side, 
    # grna_integration_strategy = 'singleton', # single grnas together?
    # grna_integration_strategy = 'union', # hopefully combines NTs?
    grna_integration_strategy = GRNA_INTEGRATION_STRATEGY,
    formula_object = 'default'
  )
  
  # =================================================================================================#
  # =================== * 3. Assign gRNAs to cells 
  # https://timothy-barry.github.io/sceptre-book/sceptre.html#sec-sceptre_assign_grnas
  print(sprintf('[%s]: SCEPTRE 3. Assign gRNAs to cells', format(Sys.time(), digits = 0)))
  sceptre_object <- assign_grnas(sceptre_object = sceptre_object, parallel = RUN_PARALLEL,
                                 method = 'maximum', 
                                 min_grna_n_umis_threshold = 1, 
                                 umi_fraction_threshold = .8)
  
  # =================================================================================================#
  # =================== * 4. Run quality control ====================================================
  # https://timothy-barry.github.io/sceptre-book/sceptre.html#sec-sceptre_qc
  print(sprintf('[%s]: SCEPTRE 4. Run quality control', format(Sys.time(), digits = 0)))
  sceptre_object <- sceptre::run_qc(sceptre_object, p_mito_threshold = 0.075)
  print(sceptre_object) # output suppressed for brevity
  plot(sceptre_object) 
  
  
  # =================================================================================================#
  # =================== * 5. Run calibration check ==================================================
  # https://timothy-barry.github.io/sceptre-book/sceptre.html#sec-sceptre_calibration_check
  print(sprintf('[%s]: SCEPTRE 5. Run calibration check', format(Sys.time(), digits = 0)))
  sceptre_object <- sceptre::run_calibration_check(sceptre_object, parallel = RUN_PARALLEL)
  
  # =================================================================================================#
  # =================== * 6. Run power check ========================================================
  # https://timothy-barry.github.io/sceptre-book/sceptre.html#sec-sceptre_run_power_check
  print(sprintf('[%s]: SCEPTRE 6. Run power check', format(Sys.time(), digits = 0)))
  sceptre_object <- run_power_check(sceptre_object, parallel = RUN_PARALLEL)
  
  # =================================================================================================#
  # =================== * 7. Run discovery analysis =================================================
  # https://timothy-barry.github.io/sceptre-book/sceptre.html#sec-sceptre_run_discovery_analysis
  print(sprintf('[%s]: SCEPTRE 7. Run discovery analysis', format(Sys.time(), digits = 0)))
  t0 = Sys.time()
  sceptre_object <- run_discovery_analysis(sceptre_object, parallel = RUN_PARALLEL)
  t1 = Sys.time()
  
  # =================================================================================================#
  # =================== * 8. Write outputs to directory =============================================
  # https://timothy-barry.github.io/sceptre-book/sceptre.html#sec-sceptre_write_outputs_to_directory
  print(sprintf('[%s]: SCEPTRE 8. Write outputs to directory', format(Sys.time(), digits = 0)))
  sceptre::write_outputs_to_directory(
    sceptre_object = sceptre_object, 
    directory = write_path
  )
  
  SCEPTRE_time_benchmark = 
    list(time_sec = difftime(t1, t0, units = 'secs'), 
         num_tests = nrow(sceptre_object@discovery_pairs),
         avg_time_per_test = difftime(t1, t0, units = 'secs')/nrow(sceptre_object@discovery_pairs),
         info = 'SCEPTRE on papalexi timing for discovery pair analysis')
  
  saveRDS(SCEPTRE_time_benchmark, file = sprintf('%s/SCEPTRE_time_benchmark.rds', write_path))
  
  
  return(sceptre_object)
}

# =================================================================================================#
# =================== START =======================================================================
# =================================================================================================#
print(sprintf("[%s] START: SCEPTRE", Sys.time()))
print(sprintf("[%s] run with no covariates:  %s  %s", Sys.time(), DEVICE, GRNA_INTEGRATION_STRATEGY))



# =================== + Set up saving dir + save setting ======================================
SCEPTRE_savepath_withcovariates= sprintf('%s/sceptre/%s/withcovariates/', save_dir, GRNA_INTEGRATION_STRATEGY) # should have already run
SCEPTRE_savepath_fake          = sprintf('%s/sceptre/%s/fake/',           save_dir, GRNA_INTEGRATION_STRATEGY)
dir.create(SCEPTRE_savepath_fake, recursive = TRUE, showWarnings = FALSE)

dir.create(sprintf('%s/nofakeU/', SCEPTRE_savepath_fake),
           sprintf('%s/withfakeU/', SCEPTRE_savepath_fake))

# create gene and grna ondisc managers
gene_odm <- ondisc::read_odm(odm_fp      = paste0(data_dir, "/papalexi-2021/processed/gene/expression_matrix.odm"),
                             metadata_fp = paste0(data_dir, "/papalexi-2021/processed/gene/metadata.rds"))
grna_odm <- ondisc::read_odm(odm_fp      = paste0(data_dir, "/papalexi-2021/processed/grna_assignment/assignment_matrix.odm"),
                             metadata_fp = paste0(data_dir, "/papalexi-2021/processed/grna_assignment/metadata.rds"))

# load grna assignments (load all into memory)
grna = grna_odm[[,1:ncol(grna_odm)]] # |> as.matrix() # ~110 x 20729 = #grnas x #cells
grna_rownames = ondisc::get_feature_ids(grna_odm)
rownames(grna) = grna_rownames


# =================== + Common sceptre setup ========================================================



# cell confounders
cell_covariates = gene_odm |> ondisc::get_cell_covariates()
extra_covariates <- cell_covariates |> 
  dplyr::mutate(lane_bio_rep = paste0(lane, '_', bio_rep)) |>
  dplyr::select(-lane, -bio_rep) |>
  dplyr::select(-n_nonzero, -n_umis)
rm(cell_covariates)

libsize = gene_odm@cell_covariates$n_umis


# grna_target_data_frame: gRNA target data frame
grna_odm_feature_covariates = ondisc::get_feature_covariates(grna_odm)
grna_target_data_frame = data.frame(grna_id     = rownames(grna_odm_feature_covariates),
                                    grna_target = grna_odm_feature_covariates$target)
rm(grna_odm_feature_covariates)

# Set discovery and positive and calibration as previous sceptre analysis's AY test pairs
#  (read in prev analysis for sceptre)
sceptre_wcov_positive    = readRDS(sprintf('%s/results_run_power_check.rds',        SCEPTRE_savepath_withcovariates))
sceptre_wcov_discovery   = readRDS(sprintf('%s/results_run_discovery_analysis.rds', SCEPTRE_savepath_withcovariates))
sceptre_wcov_calibration = readRDS(sprintf('%s/results_run_calibration_check.rds',  SCEPTRE_savepath_withcovariates))

positive_control_pairs = sceptre_wcov_positive  |> dplyr::select(grna_target, response_id) |> dplyr::distinct()
discovery_pairs        = sceptre_wcov_discovery |> dplyr::select(grna_target, response_id) |> dplyr::distinct()



# response matrix
response_matrix <- gene_odm[[,1:ncol(gene_odm)]]
colnames(response_matrix) <- gene_odm |> ondisc::get_cell_barcodes()
rownames(response_matrix) <- gene_odm |> ondisc::get_feature_ids()

# grna matrix
grna_matrix <- grna |> as.matrix() + 0 # from boolean to 0/1

# =================================================================================================#
# =================== + Test a run with the real sceptre analysis ==============
# =================================================================================================#


# =================== +--- SCEPTRE without U ======================================================
sceptre_noFakeU = run_sceptre(
  response_matrix       = response_matrix, 
  grna_matrix           = grna_matrix, 
  extra_covariates      = NULL,
  grna_target_data_frame= grna_target_data_frame,
  positive_control_pairs= positive_control_pairs,
  discovery_pairs       = discovery_pairs,
  write_path            = sprintf('%s/testreal_nocovariates/', SCEPTRE_savepath_fake)) 

# =================== +--- SCEPTRE with U ======================================================
sceptre_noFakeU = run_sceptre(
  response_matrix       = response_matrix, 
  grna_matrix           = grna_matrix, 
  extra_covariates      = extra_covariates,
  grna_target_data_frame= grna_target_data_frame,
  positive_control_pairs= positive_control_pairs,
  discovery_pairs       = discovery_pairs,
  write_path            = sprintf('%s/testreal_withcovariates/', SCEPTRE_savepath_fake)) 

# =================================================================================================#
# =================== + Add Fake and analyze for different magnitudes of confounding ==============
# =================================================================================================#

confounding_magnitudes = c(0, .05, .1, .2, .5)
for(confounding_magnitude in confounding_magnitudes) {
  # set up save dir
  SCEPTRE_savepath_fake_conf = sprintf('%s/confoundstrength=%05.2f/', SCEPTRE_savepath_fake, confounding_magnitude)
  SCEPTRE_savepath_fake_nofakeU   = sprintf('%s/nofakeU/', SCEPTRE_savepath_fake_conf)
  SCEPTRE_savepath_fake_withfakeU = sprintf('%s/withfakeU/', SCEPTRE_savepath_fake_conf)
  dir.create(SCEPTRE_savepath_fake_conf, showWarnings = FALSE)
  dir.create(SCEPTRE_savepath_fake_nofakeU, showWarnings = FALSE)
  dir.create(SCEPTRE_savepath_fake_withfakeU, showWarnings = FALSE)
  
  # =================== +--- Add fake confounding ========================================================
  data_fake = construct_fake_data(cur_grna = grna_matrix, 
                                  cur_gene = response_matrix, 
                                  # libsize  = libsize, # this is not right...
                                  confounding_size_A = confounding_magnitude, 
                                  confounding_size_Y = confounding_magnitude)
  gc()
  # format(object.size(data_fake), units = "Gb")
  
  saveRDS(object = data_fake$fakeU1, 
          file = sprintf('%s/fakeU1.rds', SCEPTRE_savepath_fake_conf))
  # =================== +--- SCEPTRE without Fake U ======================================================
  sceptre_noFakeU = run_sceptre(
    response_matrix       = data_fake$gene, 
    grna_matrix           = data_fake$grna, 
    extra_covariates       = extra_covariates,
    grna_target_data_frame= grna_target_data_frame,
    positive_control_pairs= positive_control_pairs,
    discovery_pairs       = discovery_pairs,
    write_path            = SCEPTRE_savepath_fake_nofakeU) 
  
  
  # =================== +- SCEPTRE with    Fake U ======================================================
  
  sceptre_withFakeU = run_sceptre(
    response_matrix       = data_fake$gene, 
    grna_matrix           = data_fake$grna, 
    extra_covariates       = extra_covariates |> dplyr::mutate(fakeU1=data_fake$fakeU1),
    grna_target_data_frame= grna_target_data_frame,
    positive_control_pairs= positive_control_pairs,
    discovery_pairs       = discovery_pairs,
    write_path            = SCEPTRE_savepath_fake_withfakeU) 
  
  rm(data_fake, sceptre_noFakeU, sceptre_withFakeU)
}

# 
# # =================================================================================================#
# # =================== + Add fake confounding ========================================================
# # =================================================================================================#
# data_fake = construct_fake_data(cur_grna = grna_matrix, 
#                                 cur_gene = response_matrix, 
#                                 libsize = libsize)
# gc()
# 
# # =================================================================================================#
# # =================== + SCEPTRE without Fake U ======================================================
# # =================================================================================================#
# sceptre_noFakeU = run_sceptre(
#   response_matrix       = data_fake$gene, 
#   grna_matrix           = data_fake$grna, 
#   extra_covariates       = extra_covariates,
#   grna_target_data_frame= grna_target_data_frame,
#   positive_control_pairs= positive_control_pairs,
#   discovery_pairs       = discovery_pairs,
#   write_path            = sprintf('%s/nofakeU/', SCEPTRE_savepath_fake)) 
# 
# 
# # =================================================================================================#
# # =================== + SCEPTRE with    Fake U ======================================================
# # =================================================================================================#
# 
# sceptre_withFakeU = run_sceptre(
#   response_matrix       = data_fake$gene, 
#   grna_matrix           = data_fake$grna, 
#   extra_covariates       = extra_covariates |> dplyr::mutate(fakeU1=data_fake$fakeU1),
#   grna_target_data_frame= grna_target_data_frame,
#   positive_control_pairs= positive_control_pairs,
#   discovery_pairs       = discovery_pairs,
#   write_path            = sprintf('%s/withfakeU/', SCEPTRE_savepath_fake)) 
# 
# 
# 
# # =================== END ====================================================
# print(sprintf("[%s] END", Sys.time()))
# 



N = nrow(sceptre_noFakeU@calibration_result)
plot(x = (1:N-.5)/N,y =  sort(sceptre_noFakeU@calibration_result$p_value))


N = nrow(sceptre_noFakeU@power_result)
plot(x = (1:N-.5)/N,y =  sort(sceptre_noFakeU@power_result$p_value))


plot_qqplot_quick <- function(vec) {
  N = length(vec) 
  plot(x = sort(vec), y = (1:N-.5)/N)
  abline(a = 0, b = 1)
}

plot_qqplot_quick(vec = sceptre_noFakeU@calibration_result |> filter(!is.na(p_value)) |> pull(p_value))
plot_qqplot_quick(vec = sceptre_noFakeU@discovery_result |> filter(!is.na(p_value)) |> pull(p_value))
plot_qqplot_quick(vec = sceptre_noFakeU@power_result |> filter(!is.na(p_value)) |> pull(p_value))



plot_qqplot_log_quick <- function(vec) {
  N = length(vec) 
  plot(y = -log(sort(vec) + 1e-5, base = 10), 
       x = -log((1:N-.5)/N + 1e-5, base = 10))
  abline(a = 0, b = 1)
}


plot_qqplot_log_quick(vec = sceptre_noFakeU@calibration_result |> filter(!is.na(p_value)) |> pull(p_value))
plot_qqplot_log_quick(vec = sceptre_noFakeU@discovery_result |> filter(!is.na(p_value)) |> pull(p_value))
plot_qqplot_log_quick(vec = sceptre_noFakeU@power_result |> filter(!is.na(p_value)) |> pull(p_value))




