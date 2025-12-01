# ------------------------------------------------------------------------------------------------ #
#                        Find Important Genes from papalexi-2021 dataset
#     Find top XXX genes in the papalexi dataset (processed by Tim Barry) with Will Townes'
#     method by fitting a multinomial model and using the glm deviance as the
#     normalized values (..., fam  = "binomial", type = "deviance")
# Requires: papalexi dataset saved in ondisc format at <data_dir>/papalexi-2021/
# Outputs: (nothing) but saves genes' deviance and important genes' 
#       <save_dir>/papalexi_saves/gene_deviance.csv   
#       <save_dir>/papalexi_saves/important_genes_idx.rds
#       <save_dir>/papalexi_saves/important_genes_name.rds
# ------------------------------------------------------------------------------------------------ #
args = commandArgs(trailingOnly = TRUE)
    
require(assertthat) # for some assert statements
require(ondisc)     # loading in data
require(rhdf5)      # read/save in HDF5 format
require(HDF5Array)  # (^same)
require(scry)       # normalizing data (model=binomial, residual type=deviance)
library(dplyr)

assertthat::assert_that(length(args) > 0, msg="must give arg for specifying device eg 'Rscript <filename>.R ubergenno'")

DEVICE = args[1]
source('../PATHS.R') # load in data_dir and save_dir, depending on DEVICE value
# # location of papalexi-2021 folder
# data_dir = switch(DEVICE,
#                'laptop'='C:/Users/Cathe/Documents/School/Genetic Hypothesis Testing using Negative Controls/genData/papalexi', 
#                'desktop'='C:/Users/Catherine W/Documents/Research/genData/papalexi<UPDATE!>', 
#                'ubergenno'='/raid6/Catherine/papalexi<UPDATE!>')
# # location of intermediate save files/plots are written
# save_dir = switch(DEVICE,
#                'laptop'='C:/Users/Cathe/Documents/School/Genetic Hypothesis Testing using Negative Controls/Proximal/papalexi-2021/saves', 
#                'desktop'='C:/Users/Catherine W/Documents/Research/DoubleBridge/saves/papalexi_saves<UPDATE!>', 
#                'ubergenno'='/raid6/Catherine/papalexi/papalexi_saves<UPDATE!>')

# assertthat::assert_that(!is.null(data_dir), msg='first arg must be: laptop, desktop, or ubergenno')

assertthat::assert_that(length(args) >= 2, msg="provide device and TEST_GENENUM (currently int <=4000)")
NUM_IMPORTANT_GENES = as.integer(args[2])

# =================== Start ========================================================================
print(sprintf("[%s] START: Find Important Genes papalexi", Sys.time()))

gene_odm <- ondisc::read_odm(odm_fp      = paste0(data_dir, "/papalexi-2021/processed/gene/expression_matrix.odm"),
                             metadata_fp = paste0(data_dir, "/papalexi-2021/processed/gene/metadata.rds"))
grna_odm <- ondisc::read_odm(odm_fp      = paste0(data_dir, "/papalexi-2021/processed/grna_assignment/assignment_matrix.odm"),
                             metadata_fp = paste0(data_dir, "/papalexi-2021/processed/grna_assignment/metadata.rds"))

# load all into memory
gene <- gene_odm[[,1:ncol(gene_odm)]] # entire gene ds (18649 genes x 20729 cells)
# grna <- grna_odm[[,1:ncol(grna_odm)]] # entire grna ds (110   grnas x 20729 cells)
invisible(gc(verbose=FALSE))


# =================== Choose important features ====================================================
# Input should be: row = feature, col = cell
gene_dev = scry::devianceFeatureSelection(object=gene, fam='binomial') # < 1 min
rm(gene)
invisible(gc(verbose=FALSE))


# =================== save normalized gene expression ==============================================
print(sprintf("[%s]    - saving top genes", Sys.time()))

# save/format resulting genes
gene_dev_df = data.frame(     idx = 1:length(gene_dev),
                         deviance = gene_dev,
                        gene_name = (gene_odm |> ondisc::get_feature_covariates() |> rownames()))
# important_genes_idx = sort(gene_dev, index.return=TRUE, decreasing=TRUE)$ix[1:important_genes_num] # <- DONT USE THIS. SORT REMOVES NA BY DEFAULT
important_genes_idx = gene_dev_df |> 
                      arrange(desc(deviance)) |> 
                      head(NUM_IMPORTANT_GENES) |> 
                      pull(idx)
write.csv(gene_dev_df,    file = sprintf('%s/gene_deviance.csv',       save_dir), row.names = F)

saveRDS(important_genes_idx, file = sprintf('%s/important_genes_idx.rds', save_dir))
# save the names of the important genes
saveRDS((gene_odm |> ondisc::get_feature_covariates() |> rownames())[important_genes_idx], 
        sprintf('%s/important_genes_name.rds', save_dir))

# =================== END ==========================================================================
print(sprintf("[%s] END", Sys.time()))
