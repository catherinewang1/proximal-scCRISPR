# ------------------------------------------------------------------------------------------------ #
#     Estimate ATE/Log FC
# 
# - 
# - 
# ------------------------------------------------------------------------------------------------ #
args = commandArgs(trailingOnly = TRUE)
# args = c('macbook')


assertthat::assert_that(length(args) > 0, msg="must give arg for specifying device eg 'Rscript <filename>.R ubergenno'")
DEVICE = args[1]
source('../PATHS.R') # load in data_dir and save_dir and CODE_DIR, depending on DEVICE value
assertthat::assert_that(!is.null(data_dir), msg='first arg must be: laptop, desktop, or ubergenno')

library(future.apply)
# options(future.globals.maxSize= 850*1024^2) #1st num is MB
options(future.globals.maxSize= 2500*1024^2) #1st num is MB
# plan(multisession, workers = 24)
plan(multisession, workers = 15)
# plan(sequential)

# ============================================================================== #
# =================== LOAD ==========================================
# ============================================================================== #


# process to #features x #cells to be consistent
# 
# grna: perturbation-by-cell 0/1 trt assignment
# gene:      gene-by-cell count gene expression 
# gene_norm: gene-by-cell normalized gene expression
# 
# but others may be #cells x #features
# metadata: #cells x 13 variables


# =================== \__Load counts: from causarray tutorial by Jinhong ========================================================================
sc.seurat <- readRDS(sprintf('%s/causarrayexampledata/perturbseq-exneu.rds', data_dir))

# Access the counts through Seurat when its installed version recognizes the
# serialized Assay5 object; otherwise recover the same layer and dimnames
# directly. The fallback keeps this tutorial runnable with older Seurat builds.
if ("RNA" %in% Assays(sc.seurat)) {
  gene <- GetAssayData(sc.seurat, assay = "RNA", layer = "counts")
} else {
  rna.assay <- sc.seurat@assays[["RNA"]]
  gene <- rna.assay@layers[["counts"]]
  rownames(gene) <- rownames(rna.assay@features)
  colnames(gene) <- rownames(rna.assay@cells)
  rm(rna.assay)
}

metadata <- sc.seurat@meta.data
perturb <- metadata
colnames(perturb) <- gsub("Perturbation", "trt_", colnames(perturb))
perturb$trt_ <- relevel(as.factor(perturb$trt_), ref = "GFP")
A <- data.frame(
  model.matrix(~ trt_ - 1, data = perturb)[, -1, drop = FALSE],
  check.names = FALSE
) # cell-by-trt matrix; remove the first (GFP control) column
colnames(A) <- sub("^trt_", "", colnames(A))

grna = as.matrix(t(A))

rm(sc.seurat, perturb, A) # only keep: gene, grna, metadata

# =================== \__Load continuous: normalized counts ==============================================

gene_norm = readRDS(sprintf('%s/gene_norm.rds', save_dir))
gene_norm = as.matrix(gene_norm)
(dim(gene_norm) == dim(gene)) |> all() # same dimensions (all genes)


# =================== \__Load gene info (dev, TF, ...) ==============================================

gene_dev_df = read.csv(file = sprintf('%s/gene_dev.csv', save_dir)) 

(gene_dev_df$gene_name == row.names(gene))      |> table() # should be all TRUE
(gene_dev_df$gene_name == row.names(gene_norm)) |> table() # 



gene_dev_df |> head()
# gene_name gene_idx  deviance importance_rank is_grna_target genename_cleaned is_transcription_factor
#      Nrgn        1 154587.25               2          FALSE             nrgn                   FALSE
#    Snhg11        2  27156.25              40          FALSE           snhg11                   FALSE
#      Actb        3 370086.62               1          FALSE             actb                   FALSE
#     Mef2c        4  56872.76              16          FALSE            mef2c                    TRUE
#     Cd24a        5  39437.06              26          FALSE            cd24a                   FALSE
#   Camk2n1        6  20781.46              61          FALSE          camk2n1                   FALSE

# ============================================================================== #
# =================== Choose AY tests ==========================================
# ============================================================================== #


# --------  discovery tests  ------ # 
# filter for perturbations (in this example dataset, they seem ok already, but we can filter here too)
rowSums(grna) |> hist(breaks = 50)
GRNA_SAMPLESIZE_MIN = 80 # just set a thresh here based on distn

all_As = row.names(grna)[rowSums(grna) >= GRNA_SAMPLESIZE_MIN]


# filter for genes: larger #cells with non-zero expr and higher rank
rowSums(gene == 0) |> hist(breaks = 50)
GENE_NUMNONZEROCELLS_MIN = 500 # num of cells w nonzero counts
GENE_IMPORTANCERANK_MAX = 1500 # maximum importance rank of gene, ranked by gene deviance

all_Ys = row.names(gene)[(rowSums(gene == 0) >= GENE_NUMNONZEROCELLS_MIN) & 
                         (gene_dev_df$importance_rank <= GENE_IMPORTANCERANK_MAX)]


NUM_DISC_TESTS = 2000
disc_tests_df = expand.grid(grna=all_As, gene=all_Ys, stringsAsFactors = FALSE) |> dplyr::filter(grna != gene)
disc_tests_df = disc_tests_df[sample(nrow(disc_tests_df), size = NUM_DISC_TESTS), ] |> mutate(type = 'discovery')


# --------  'positive' tests  ------ # 
# also require 'positive' tests, even if they did not pass the qc (some grna targeted genes are not here in the subset...)
pos_tests_vec = row.names(grna)[row.names(grna) %in% row.names(gene)]
pos_tests_df = data.frame(grna = pos_tests_vec, gene = pos_tests_vec, type = 'positive')

# --------  assemble together and save  ------ # 
AY = rbind(pos_tests_df, disc_tests_df)

dir.create(sprintf('%s/AY/', save_dir))
dir.create(sprintf('%s/AY/exampledata/', save_dir))
write.csv(x = AY, file = sprintf('%s/AY/exampledata/AY.csv', save_dir), row.names = FALSE) # csv 

rm(GRNA_SAMPLESIZE_MIN, GENE_NUMNONZEROCELLS_MIN, GENE_IMPORTANCERANK_MAX, NUM_DISC_TESTS)
rm(all_As, all_Ys, disc_tests_df, pos_tests_vec, pos_tests_df, gene_dev_df) # prob keep all_Ys to have a collection of genes to use as NCs



# ============================================================================== #
# =================== Some Estimates (parallel) ================================
# ============================================================================== #
dir.create(sprintf('%s/AY/exampledata/intermediateATEs/', save_dir))


NCs = readRDS(sprintf('%s/pca/NCloadings.rds', save_dir))  # PCA loadings
numNCs = c(1, 3, 5, 10, 15, 20)
libsize_log = log((metadata$nUMI |> as.numeric())+ 1e-20)
control_cell_idx = which(colSums(grna) == 0) # cells with GFP

# control_cell_idx2 = which(metadata$Perturbation == 'GFP') # cells with GFP  table(control_cell_idx == control_cell_idx2)


# remove all non-relevant objects in environment
# rm(setdiff(objects(), 
#            c('grna', 'gene', 'gene_norm', 'AY', 
#              'libsize_log', 'control_cell_idx', 'NCs', 'numNCs',
#              'DEVICE', 'save_dir', 'util_dir', 'data_dir')))



#' estimate ATEs/logFoldChanges/Coefficients using loaded in vars in environment
#' Specifically: 
#'    grna, gene, gene_norm, AY
#'    libsize_log, control_cell_idx, NCs, numNCs
#'    
#' res is dataframe with cols:
#' - grna, gene, type 
#' - ... etc ...
#' - tstat: test statistic either the t value (lm) or z value (glm, proximal)
#' @param AY_idx (integer) index of the AY test 
#' 
estimate_effects <- function(AY_idx, save_intermediateATEs=TRUE) {
  res = NULL
  A_name = AY[AY_idx, 'grna']
  Y_name = AY[AY_idx, 'gene']
  
  # ===== assemble df =====
  grna_cell_idx = which(grna[A_name, ] == 1)
  AY_data_idx = c(control_cell_idx, grna_cell_idx)
  
  df = data.frame(# A2 = c(rep(0, length(control_cell_idx)), rep(1, length(grna_cell_idx))), # another way to get 0/1
    A           = grna[A_name, AY_data_idx],
    Y           = gene[Y_name, AY_data_idx], 
    Y_norm      = gene_norm[Y_name, AY_data_idx],
    libsize_log = libsize_log[AY_data_idx])
  
  # Proximal with PCA loadings
  dfZ = NCs[AY_data_idx, seq(from = 2, to = min(2*max(numNCs), ncol(NCs)), by = 2)] # evens which((1:ncol(NCs)) %% 2 == 0)
  dfW = NCs[AY_data_idx, seq(from = 1, to = min(2*max(numNCs), ncol(NCs)), by = 2)] # odds
  colnames(dfZ) = paste0('Z', 1:ncol(dfZ))
  colnames(dfW) = paste0('W', 1:ncol(dfW))
  
  
  # ===== estimate =====
  # ===== \__Linear Model without any covariates =====
  t0 = Sys.time()
  lmfit = lm(formula='Y_norm ~ A', data = df)
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
  
  
  # ===== \__ Negative Binomial Fit without any covariates ====
  t0 = Sys.time()
  nbfit = MASS::glm.nb(formula='Y ~ A + offset(libsize_log)', data = df)
  nbfit_summary = summary(nbfit)
  t1 = Sys.time()
  
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
  rm(nbfit, nbfit_summary, t1, t0)
  

  # ===== \__ Proximal PCA ====
  for(numNCs_ in numNCs)  {
    t0 = Sys.time()
    pci2s_res = pci2s::p2sls.lm(
      Y = df$Y_norm, 
      A = df$A, 
      W = dfZ[,1:numNCs_], 
      Z = dfW[,1:numNCs_], 
      variance = TRUE)
    t1 = Sys.time()
    
    res = bind_rows(res, 
                    data.frame(
                      NC_type     = 'PCA',
                      method      = '2SLSpci2s',
                      method_type = 'proximal',
                      numNC       = numNCs_,
                      ATE = pci2s_res$summary_second_stage['A', 'Estimate'],
                      se  = pci2s_res$summary_second_stage['A', 'Std. Error'],
                      tstat= pci2s_res$summary_second_stage['A', 'z value'],
                      pval= pci2s_res$summary_second_stage['A', 'Pr(>|z|)'],
                      time_sec = difftime(t1, t0, units = 'secs') |> as.numeric()))
    
    rm(pci2s_res, t0, t1)
  }
  rm(dfZ, dfW, numNCs_)
  
  # ===== return ==== 
  res = res |> mutate(AY_idx = AY_idx, grna = A_name, gene = Y_name, type = AY[AY_idx, 'type'], .before = 1)
  
  if(save_intermediateATEs) {
    write.csv(res, sprintf('%s/AY/exampledata/intermediateATEs/ATE%s.csv', save_dir, AY_idx))
  }
  return(res)
}



# test a run
# estimate_effects(4)


# =================== Get ATEs (parallel) ==========================================================
print(sprintf("[%s]    - Get Estimates (parallel)", Sys.time()))

# whichROWS = 1:5
whichROWS = 1:nrow(AY)


t0 = Sys.time()
effects_par = future.apply::future_mapply(estimate_effects,
                                      AY_idx = whichROWS, 
                                      future.globals = TRUE,
                                      future.seed = 123456, 
                                      SIMPLIFY = FALSE)
t1 = Sys.time()
print(sprintf("[%s]        - %2.2f mins", Sys.time(), difftime(t1, t0, units = "mins")))


saveRDS(effects_par, file = sprintf('%s/AY/exampledata/effects_par.rds', save_dir)) # save this parallel res as rds



# =================== Combine ATEs (into one df) ===================================================
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
write.csv(x = effects_df, file = sprintf('%s/AY/exampledata/effects.csv', save_dir), row.names = FALSE)








# ==== TRASH ==== TRASH ==== TRASH ==== TRASH ==== TRASH ==== TRASH ==== TRASH 
# ============================= TRASH ========================================
# ==== TRASH ==== TRASH ==== TRASH ==== TRASH ==== TRASH ==== TRASH ==== TRASH 











# ==== TRASH ==== TRASH ==== TRASH ==== TRASH ==== TRASH ==== TRASH ==== TRASH 
# ============================= TRASH ========================================
# ==== TRASH ==== TRASH ==== TRASH ==== TRASH ==== TRASH ==== TRASH ==== TRASH 

if(F) {
  
  # ============================================================================== #
  # =================== Some Estimates (seq) =====================================
  # ============================================================================== #
  
  NCs = readRDS(sprintf('%s/pca/NCloadings.rds', save_dir))  # PCA loadings
  numNCs = c(1, 3, 5, 10, 15, 20)
  
  libsize_log = log((metadata$nUMI |> as.numeric())+ 1e-20)
  
  
  control_cell_idx = which(colSums(grna) == 0) # cells with GFP
  # control_cell_idx2 = which(metadata$Perturbation == 'GFP') # cells with GFP  table(control_cell_idx == control_cell_idx2)
  
  res_all = NULL
  # whichROWS = 1:nrow(AY)
  whichROWS = 1:5
  for(AY_idx in whichROWS) {
    
    res = NULL
    A_name = AY[AY_idx, 'grna']
    Y_name = AY[AY_idx, 'gene']
    
    # assemble df
    grna_cell_idx = which(grna[A_name, ] == 1)
    AY_data_idx = c(control_cell_idx, grna_cell_idx)
    
    df = data.frame(# A2 = c(rep(0, length(control_cell_idx)), rep(1, length(grna_cell_idx))), # another way to get 0/1
      A           = grna[A_name, AY_data_idx],
      Y           = gene[Y_name, AY_data_idx], 
      Y_norm      = gene_norm[Y_name, AY_data_idx],
      libsize_log = libsize_log[AY_data_idx])
    
    # Negative Binomial Fit without any covariates
    formula = 'Y ~ A + offset(libsize_log)'
    t0 = Sys.time()
    nbfit = MASS::glm.nb(formula=formula, data = df)
    nbfit_summary = summary(nbfit)
    t1 = Sys.time()
    
    res = bind_rows(res, 
                    data.frame(
                      NC_type     = NA,
                      method      = 'negbin',
                      method_type = 'glm',
                      numNC       = NA,
                      ATE = nbfit_summary$coefficients['A', 'Estimate'],
                      se  = nbfit_summary$coefficients['A', 'Std. Error'],
                      zval= nbfit_summary$coefficients['A', 'z value'],
                      pval= nbfit_summary$coefficients['A', 'Pr(>|z|)'],
                      time_sec = difftime(t1, t0, units = 'secs') |> as.numeric()))
    
    
    # Proximal with PCA loadings
    dfZ = NCs[AY_data_idx, seq(from = 2, to = min(2*max(numNCs), ncol(NCs)), by = 2)] # evens which((1:ncol(NCs)) %% 2 == 0)
    dfW = NCs[AY_data_idx, seq(from = 1, to = min(2*max(numNCs), ncol(NCs)), by = 2)] # odds
    colnames(dfZ) = paste0('Z', 1:ncol(dfZ))
    colnames(dfW) = paste0('W', 1:ncol(dfW))
    
    for(numNCs_ in numNCs)  {
      t0 = Sys.time()
      pci2s_res = pci2s::p2sls.lm(
        Y = df$Y_norm, 
        A = df$A, 
        W = dfZ[,1:numNCs_], 
        Z = dfW[,1:numNCs_], 
        variance = TRUE)
      t1 = Sys.time()
      
      res = bind_rows(res, 
                      data.frame(
                        NC_type     = 'PCA',
                        method      = '2SLSpci2s',
                        method_type = 'proximal',
                        numNC       = numNCs_,
                        ATE = pci2s_res$summary_second_stage['A', 'Estimate'],
                        se  = pci2s_res$summary_second_stage['A', 'Std. Error'],
                        zval= pci2s_res$summary_second_stage['A', 'z value'],
                        pval= pci2s_res$summary_second_stage['A', 'Pr(>|z|)'],
                        time_sec = difftime(t1, t0, units = 'secs') |> as.numeric()))
      
      rm(pci2s_res, t0, t1)
    }
    rm(dfZ, dfW)
    
    res_all = bind_rows(res_all, 
                        res |> mutate(AY_idx = AY_idx, grna = A_name, gene = Y_name, type = AY[AY_idx, 'type'], .before = 1))
    
  }
  
  # ============================================================================== #
  # =================== Negative Binomial GLM (for a lot) ========================
  # ============================================================================== #
  
  # sometimes the names are bad for the glm formulas (e.g. AB-C). Remove non letter and number characters 
  clean_name <- function(x) {  gsub(pattern='[^a-zA-Z0-9]', replacement='', x=x)  }
  
  
  
  # Neg Bin fits without any covariates
  libsize = data.frame(libsize_log=log(metadata$nUMI |> as.numeric() + 1e-20))
  # libsize_log = log((metadata$nUMI |> as.numeric())+ 1e-20)
  
  
  
  
  # which_genes = colnames(Y) # all genes
  N_genes = 1000   # or randomly choose some
  set.seed(12345)
  
  perturbations_names = unique(metadata$Perturbation) |> sort()
  pert_targetgene_in_subsetdata = intersect(perturbations_names, colnames(Y)) # always test these bc they are pert targets
  which_genes = union(pert_targetgene_in_subsetdata, 
                      sample(setdiff(colnames(Y), perturbations_names),  replace = F, size = N_genes - length(pert_targetgene_in_subsetdata)))
  
  # fit and save res: 1 glm per gene (vs individual pert-gene fits)
  res = NULL
  for(Y_name in which_genes) {
    
    df = cbind(Y[, Y_name, drop=FALSE] |> `colnames<-`(paste0('gene_', clean_name(Y_name))), 
               libsize,
               A
    )
    
    
    # formula = paste0('gene_', Y_name, ' ~ offset(libsize_log) + ', paste0(colnames(A), collapse = ' + '))
    formula = paste0('gene_', clean_name(Y_name), 
                     ' ~ offset(libsize_log) + ', 
                     paste0(sapply(colnames(A), clean_name), collapse = ' + '))
    nbfit = MASS::glm.nb(formula=formula, data = df)
    
    res = rbind(res, 
                summary(nbfit)$coefficients[-1, ] |> 
                  `colnames<-`(c('estimate', 'se', 'zvalue', 'pvalue')) |> 
                  data.frame() |> 
                  dplyr::mutate(grna=colnames(A), gene = Y_name, .after="estimate") )
    
    
    # summary(nbfit)$coefficients[colnames(A), ] |> 
    #   `colnames<-`(c('estimate', 'se', 'zvalue', 'pvalue')) |> 
    #   data.frame() |> 
    #   tibble::rownames_to_column("grna") |> 
    #   dplyr::mutate(gene = Y_name, .after="grna")
  }
  
  
  
  
  
}
