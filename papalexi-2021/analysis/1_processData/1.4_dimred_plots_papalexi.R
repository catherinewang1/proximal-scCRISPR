# ---------------------------------------------------------------------------- #
#                   Plot 2-d dim reductions (pca, pca-umap)                    #
# Requires: prev saved 2-d dim reductions                                      #
# Ouputs: (nothing) but saves plots at                                         #
#         '<data_dir>/papalexi_saves/pca/'                                     #
#         '<data_dir>/papalexi_saves/umap/'                                    #
#         '<data_dir>/papalexi_saves/pca_umap/'                                #
# ---------------------------------------------------------------------------- #
args = commandArgs(trailingOnly = TRUE)
    
require(assertthat) # for some assert statements
# require(rhdf5)      # read/write HDF5 format
require(ondisc)
# library(umap)       # perform umap
library(ggplot2)
library(dplyr)

assertthat::assert_that(length(args) > 0, msg="must give arg for specifying device eg 'Rscript <filename>.R ubergenno'")
DEVICE = args[1]
source('../PATHS.R') # load in data_dir and save_dir and CODE_DIR and util_dir, depending on DEVICE value
# # location of papalexi-2021 folder
# data_dir = switch(DEVICE,
#                'laptop'='C:/Users/Cathe/Documents/School/Genetic Hypothesis Testing using Negative Controls/genData/papalexi', 
#                'desktop'='C:/Users/Catherine W/Documents/Research/genData/papalexi', 
#                'ubergenno'='/raid6/Catherine/papalexi')
# # location of intermediate save files/plots are written
# save_dir = switch(DEVICE,
#                'laptop'='C:/Users/Cathe/Documents/School/Genetic Hypothesis Testing using Negative Controls/DoubleBridge/saves/papalexi_saves', 
#                'desktop'='C:/Users/Catherine W/Documents/Research/DoubleBridge/saves/papalexi_saves', 
#                'ubergenno'='/raid6/Catherine/papalexi/papalexi_saves')
# # location of utils code
# CODE_DIR = switch(DEVICE,
#                   'laptop'='C:/Users/Cathe/Documents/School/Genetic Hypothesis Testing using Negative Controls/DoubleBridge/code/utils', 
#                   'desktop'='C:/Users/Catherine W/Documents/Research/DoubleBridge/code/utils', 
#                   'ubergenno'='/raid6/home/catheri2/DoubleBridge/code/utils')

assertthat::assert_that(!is.null(data_dir), msg='first arg must be: laptop, desktop, or ubergenno')

assertthat::assert_that(length(args) >= 2, msg="provide device and NUM_IMPORTANT_GENES (currently int <=4000)")
NUM_IMPORTANT_GENES = as.integer(args[2])





# =================== Start ========================================================================
print(sprintf("[%s] START: Normalize papalexi", Sys.time()))

papalexi_dir = paste0(data_dir, '/papalexi-2021')
gene_odm <- read_odm(odm_fp      = paste0(papalexi_dir, "/processed/gene/expression_matrix.odm"),
                     metadata_fp = paste0(papalexi_dir, "/processed/gene/metadata.rds"))
grna_odm <- read_odm(odm_fp      = paste0(papalexi_dir, "/processed/grna_assignment/assignment_matrix.odm"),
                     metadata_fp = paste0(papalexi_dir, "/processed/grna_assignment/metadata.rds"))

chosen_genes = c('FTH1', 'EEF1A1', 'FAU') 
chosen_grna_T  = c('IFNGR2g1', 'ATF2g1', 'CD86g1')
all_grna_NT = grna_odm |> get_feature_covariates() |> filter(target_type == 'non-targeting') |> rownames()





# =================== Constructing Dataframes of Dim Reds ========================
print(sprintf("[%s]    - Constructing dataframes of dim reds", Sys.time()))

source(sprintf('%s/dataProcessing/5_dimred_plots_utils.R', util_dir))

papalexi_pca      = readRDS(sprintf('%s/pca/objects/pca_%s.rds', save_dir, NUM_IMPORTANT_GENES))
papalexi_pca_umap = readRDS(sprintf('%s/pca_umap/objects/pca_umap_%s.rds', save_dir, NUM_IMPORTANT_GENES))

# dataframe indicating gRNA each cell received
grna_assignment = get_grna_assignment(grna_odm = grna_odm)

# function for getting dim red for same (non-dimred method) params 
my_get_dimred = get_dimred_bydimredmethod(grna_assignment = grna_assignment, 
                                          cell_covariates = gene_odm |> get_cell_covariates(),
                                          gRNA_names = chosen_grna_T,
                                          NT_names = all_grna_NT,
                                          NT_exposed_size=NULL,
                                          seed = 12345) 

pca_cells      = my_get_dimred(df_dimred = as.data.frame(papalexi_pca$x[,1:2]))
# umap_cells     = my_get_dimred(df_dimred = as.data.frame(papalexi_umap$layout)) # no? only umap after pca
pca_umap_cells = my_get_dimred(df_dimred = as.data.frame(papalexi_pca_umap$layout))


# =================== Plotting Dim Reductions ========================
print(sprintf("[%s]    - Plotting dim red", Sys.time()))
# function for plotting dim red for same (non-dimred method) params
my_plot_dimred = plot_dimred_bydimredmethod(save_dir=save_dir, 
                                            num_imp_genes = NUM_IMPORTANT_GENES, 
                                            gRNA_names = chosen_grna_T,
                                            NT_names = all_grna_NT)

my_plot_dimred(dimred_type='pca',      gRNAcells=pca_cells$gRNAcells,      allcells=pca_cells$allcells)
# my_plot_dimred(dimred_type='umap',     gRNAcells=umap_cells$gRNAcells,     allcells=umap_cells$allcells)
my_plot_dimred(dimred_type='pca_umap', gRNAcells=pca_umap_cells$gRNAcells, allcells=pca_umap_cells$allcells)



# =================== END ==========================================================================
print(sprintf("[%s] END", Sys.time()))

