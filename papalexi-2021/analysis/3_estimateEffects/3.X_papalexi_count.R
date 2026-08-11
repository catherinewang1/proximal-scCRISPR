# ------------------------------------------------------------------------------------------------ #
#     Estimate Effect (ATE) with count (original raw counts) outcome Y 
# e.g.     nice -10 Rscript 3.2_papalexi_count.R ubergenno A
# GLM Methods:
#  - Poisson GLM
#  - Poisson GLM           regression with 'un'-measured confounders
#  - Negative Binomial GLM
#  - Negative Binoimial GLM regression with 'un'-measured confounders
# Proximal Methods:
#  - singlegenes
#  - PCA
#  - Sparse PCA 34.5
#  - WGCNA
# ------------------------------------------------------------------------------------------------ #
args = commandArgs(trailingOnly = TRUE)
# args = c('laptop',  'A1')
# args = c('laptop',  'A')
# args = c('macbook', 'A')



args = c('macbook', 'A1', 'simpleCount-proximalNegBinCountCount-proximalNegBinPCAPCA-proximalNegBinSinglegeneSinglegene') # do many at the same time


# =================== Load Libraries and Packages =========================================================
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

# =================== Read input args =========================================================
save_intermediateATEs = 'yes' # 'yes'/'no' whether to save intermedate ATEs as they are estimated

# args[1] = device e.g. macbook or laptop or ubergenno
assertthat::assert_that(length(args) > 0, msg="must give arg for specifying device eg 'Rscript <filename>.R ubergenno'")
DEVICE = args[1]
source('../PATHS.R') # load in data_dir and save_dir and CODE_DIR, depending on DEVICE value
assertthat::assert_that(!is.null(data_dir), msg='first arg must be: laptop, desktop, or ubergenno')

# args[2] = setting name e.g. A1 or A or B
assertthat::assert_that(length(args) > 1, msg="must give arg for specifying chosen AYZW name 'Rscript <filename>.R ubergenno C'")
AYZW_setting_name = args[2]

# args[3] = analysis methods separated by '-' e.g. simpleCount-proximalNegBinCountCount
# should be the names of the list made in saves/AY/proximal_continuous_settings.r
assertthat::assert_that(length(args) > 2, msg="must give arg for specifying chosen analysis methods 'Rscript <filename>.R ubergenno C SPCA8.0' (see options listed in AY/proximal_continuous_settings.r)")
NC_name_raw = args[3]
NC_names = strsplit(x = NC_name_raw, split = '-') |> unlist()

source(sprintf('%s/AY/proximal_continuous_settings.r', save_dir)) # loads in list: proximal_continuous_settings

assertthat::assert_that(all(NC_names %in% names(proximal_continuous_settings)), msg = "bad NC_name input, list mult like PCA-SPCA8.0-singlegene (see options listed in AY/proximal_continuous_settings.r)")


for(NC_name in NC_names) { # save this run's settings from proximal_continuous_settings + add device
  PROXIMAL_SETTINGS = proximal_continuous_settings[[NC_name]] # the current proximal settings
  PROXIMAL_SETTINGS$DEVICE = DEVICE
  # save parameter settings
  dir.create(sprintf('%s/AY/%s/%s', save_dir, AYZW_setting_name, NC_name), recursive = FALSE, showWarnings = FALSE)
  capture.output(print(PROXIMAL_SETTINGS),
                 file = sprintf('%s/AY/%s/%s/proximal_setting.txt',
                                save_dir, AYZW_setting_name, NC_name))
  saveRDS(PROXIMAL_SETTINGS,
          sprintf('%s/AY/%s/%s/proximal_setting.rds',
                  save_dir, AYZW_setting_name, NC_name))
  
  rm(NC_name, PROXIMAL_SETTINGS) # clean environment
}



intermediateATEs_folder = sprintf('%s/AY/%s/intermediateATEsCount', save_dir, AYZW_setting_name)
if(save_intermediateATEs == 'yes') {
  dir.create(intermediateATEs_folder, showWarnings = FALSE)
  cat(NC_name_raw, file = sprintf('%s/%s.txt', intermediateATEs_folder, NC_name_raw))
}


# =================== Start ========================================================================
print(sprintf("[%s] START: Estimate Papalexi-2021 as counts", Sys.time()))
print(sprintf("[%s]        w/ script: %s", Sys.time(), paste0(args, collapse = ', ')))




print(sprintf("[%s]    - Load files and functions", Sys.time()))

source(sprintf('%s/estimate_effects_pci2s_count.R', util_dir)) # for functions to estimate here


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

# Load choosing AY (ZW) settings
# AYZW_setting = readRDS(sprintf('%s/AY/%s/AYZW_setting.rds', save_dir, AYZW_setting_name))


# create gene and grna ondisc managers
gene_odm <- ondisc::read_odm(odm_fp      = paste0(data_dir, "/papalexi-2021/processed/gene/expression_matrix.odm"),
                             metadata_fp = paste0(data_dir, "/papalexi-2021/processed/gene/metadata.rds"))
grna_odm <- ondisc::read_odm(odm_fp      = paste0(data_dir, "/papalexi-2021/processed/grna_assignment/assignment_matrix.odm"),
                             metadata_fp = paste0(data_dir, "/papalexi-2021/processed/grna_assignment/metadata.rds"))

# load grna assignments (load all into memory)
# grna = grna_odm[[,1:ncol(grna_odm)]] # |> as.matrix() # ~110 x 20729 = #grnas x #cells
# grna_rownames = grna_odm |> ondisc::get_feature_covariates() |> rownames()
# row.names(grna) = grna_rownames

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







NCs_list = list()
for(NC_name in NC_names) {
  PROXIMAL_SETTINGS = proximal_continuous_settings[[NC_name]]
  
  if(NC_name == 'simpleCount') {                               # regular GLMs (very fast, can just keep)
    NCs = NA
  } else if(NC_name == 'proximalNegBinCountCount') {           # prox with NegBin Y, count NCE, count NCO (individual genes)
    # set the NCs as list of chosen Z_names and W_names for each AY test
    NCs = readRDS(sprintf('%s/AY/%s/AYZW.rds', save_dir, AYZW_setting_name))
  } else if(NC_name == 'proximalNegBinPCAPCA') {               # prox with NegBin Y, continuous NCE, continuous NCO (PCA)
    # use PCA Loadings
    NCs = readRDS(sprintf('%s/pca/NCloadings.rds', save_dir)) |> as.matrix()
  } else if(NC_name == 'proximalNegBinSinglegeneSinglegene') { # prox with NegBin Y, continuous NCE, continuous NCO (singlegenes)
    # load in gene_norm
    
    # load normalized gene exp
    h5file      = paste0(save_dir, "/gene.h5")#; print(h5file)
    reading_hd5file  = rhdf5::H5Fopen(name = h5file)
    readin_gene_norm = reading_hd5file&'gene_norm'
    gene_norm = readin_gene_norm[, 1:ncol(gene_odm)] # eg dim = 4000 x 20729 = #important x #cells
    rhdf5::h5closeAll()
    invisible(gc(verbose=FALSE))
    rm(h5file, reading_hd5file, readin_gene_norm)
    
    # set the NCs as list of chosen Z_names and W_names for each AY test
    NCs = readRDS(sprintf('%s/AY/%s/AYZW.rds', save_dir, AYZW_setting_name))
    
    # subset gene_norm to only the individual genes we need (to make the usage space smaller)
    # s.t. the environment is smaller? then it is potentially faster? 
    # e.g. only include gene_norm that are needed: (single genes will add a lot of Ys...)
    # all_As = AY$A |> unique() # grna size is small enough
    all_singlegenes = c()
    max_num_NC_pairs = proximal_continuous_settings[['proximalNegBinSinglegeneSinglegene']]$num_NC_pairs |> max() # max number of ZW used 
    for(A_name in names(NCs)) {
      for(Y_name in names(NCs[[A_name]])) {
        all_singlegenes = c(all_singlegenes, 
                            NCs[[A_name]][[Y_name]][[1]]$Z_names[1:max_num_NC_pairs] ,
                            NCs[[A_name]][[Y_name]][[1]]$W_names[1:max_num_NC_pairs] ) |> unique()
        
      }
    }
    rm(A_name, Y_name, max_num_NC_pairs)
    
    # dim(gene_norm)        # 4017 20729
    # row.names(gene_norm)  # NULL
    # head(gene_importance) # <-- look here! cols: gene_name gene_idx importance_rank gene_norm_idx
    gene_importance_subset      = gene_importance |>                                # idx of used genes in gene_norm      
      dplyr::filter(gene_name %in% all_singlegenes) |> 
      dplyr::arrange(gene_norm_idx)
    gene_norm_subset            = gene_norm[gene_importance_subset$gene_norm_idx, ] # subset from gene_norm
    row.names(gene_norm_subset) =           gene_importance_subset$gene_name        # name the rows to gene names
    
    # could also speed up by subselecting cells, but bc we will likely use all treatments --> all cells, so not worth it for now
    # it would be helpful if we do the matrix form though
    
    
    rm(gene_norm, gene_importance_subset, all_singlegenes)
    
  } else {
    print(sprintf('Bad NC_name: %s', NC_name)) # should've already been caught at the beginning. no? or i meant to and i didnt do it
  }
  
  NCs_list[[NC_name]] = NCs
  
  rm(NC_name, PROXIMAL_SETTINGS) # clean environment
}







# All the NT grna idx
NT_names = grna_odm |> ondisc::get_feature_covariates() |>
  filter(target_type == 'non-targeting') |> rownames()
NT_idx = which(apply(X = grna_odm[[NT_names, ]], MARGIN = 2, FUN = sum) > 0)
rm(NT_names)






# =================== Define fns for Effects (parallel) ============================================
print(sprintf("[%s]    - Define fns for Effects (parallel)", Sys.time()))

source(sprintf('%s/estimate_effects_pci2s_count.R', util_dir)) # for functions to estimate here


estimate_effect_0 = estimate_ATE_pci2snegbin_make(AY = AY, 
                                                  gene_odm = gene_odm, 
                                                  grna_odm = grna_odm, 
                                                  gene_norm= gene_norm_subset, 
                                                  NT_idx   = NT_idx,
                                                  NCs_list = NCs_list,                                 # list of actual NCs (or the names of ZW)
                                                  NCs_settings=proximal_continuous_settings[NC_names], # settings for chosen methods
                                                  library_size=library_size, 
                                                  U_confounders=cell_covariates,
                                                  save_path=switch(save_intermediateATEs,
                                                                   'yes' = intermediateATEs_folder,
                                                                   'no'  = NULL)
                                                  ) 
test = estimate_effect_0(AY_idx = 3)



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




