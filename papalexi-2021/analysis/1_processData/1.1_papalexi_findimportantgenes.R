# ------------------------------------------------------------------------------------------------ #
#                        Find Important Genes from papalexi-2021 dataset
#                        (and mark transcription factor genes)
#     Find top XXX genes in the papalexi dataset (processed by Tim Barry) with Will Townes' 
#     (scry::devianceFeatureSelection) method by fitting a multinomial model and using the glm 
#     deviance as the normalized values (..., fam  = "binomial", type = "deviance")
# Requires: papalexi dataset saved in ondisc format at <data_dir>/papalexi-2021/
# Outputs: (nothing) but saves genes' deviance and important genes' 
#       <save_dir>/gene_deviance.csv   # cols: gene_idx, deviance, gene_name, is_transcription_factor, importance_rank
#       <save_dir>/gene_deviance_topnoTFonly.csv # redundant, but saves space bc smaller NUM_IMPORTANT_GENES vs all genes
#       <save_dir>/important_genes_idx.rds  (CHANGE: remove, only need gene_deviance.csv, redundant)
#       <save_dir>/important_genes_name.rds (CHANGE: remove, only need gene_deviance.csv, redundant) 
# ------------------------------------------------------------------------------------------------ #
args = commandArgs(trailingOnly = TRUE)
# args = c('laptop', '4000')
    
suppressPackageStartupMessages(require(assertthat)) # for some assert statements
suppressPackageStartupMessages(require(ondisc))     # loading in data !!! USE 1.1.0 RELEASE: https://github.com/timothy-barry/ondisc/releases, install from .tar.gz file install.packages('C:/Users/Cathe/Documents/School/Genetic Hypothesis Testing using Negative Controls/genData/package_ondisc_old/ondisc-1.1.tar.gz', repos = NULL, type ='source')
suppressPackageStartupMessages(require(rhdf5))      # read/save in HDF5 format # BiocManager::install("rhdf5")
suppressPackageStartupMessages(require(HDF5Array))  # (^same)                  # BiocManager::install("HDF5Array")
suppressPackageStartupMessages(require(scry))       # normalizing data (model=binomial, residual type=deviance)
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(require(readxl)) # load in excel xlsx format for TF list


assertthat::assert_that(length(args) > 0, msg="must give arg for specifying device eg 'Rscript <filename>.R ubergenno'")

DEVICE = args[1]
source('../PATHS.R') # load in data_dir and save_dir and util_dir, depending on DEVICE value

assertthat::assert_that(length(args) >= 2, msg="provide device and TEST_GENENUM (currently int <=4000)")
NUM_IMPORTANT_GENES = as.integer(args[2])

# =================== Start ========================================================================
print(sprintf("[%s] START: Find Important Genes papalexi", Sys.time()))

# Functions changed?? read_odm not a function in ondisc library??, download older version 1.1.0 to fit
# https://github.com/Katsevich-Lab/import-papalexi-2021
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

# =================== Load Transcription Factor (TF) genes list ======================================
print(sprintf("[%s]    - loading list of TF genes", Sys.time()))
# load list of genes that are Transcription Factors (read in xlsx sheet)
tf_raw = readxl::read_xlsx(path = paste0(data_dir, "/../extra/transcriptionfactorlist.xlsx"),  # given from Kathryn over slack (suppl of a paper?)
                           sheet = 2) |> suppressMessages() # suppress messages on how they renamed columns
# clean up a bit
#      2nd col ('...2') is gene name, 4th col is TF indicator 'Is TF?'
#      First row is not data (2 rows of colnames)
tf = tf_raw[-1 , c('...2', 'Is TF?')]
colnames(tf) = c('gene_name', 'TF')
# table(tf$TF) # only levels are No, Yes
#        No  Yes 
#       1126 1639 
TF_names = tf |> dplyr::filter(TF == 'Yes') |> dplyr::pull(gene_name)


# =================== save gene deviances (='importance') ==============================================
print(sprintf("[%s]    - saving top genes", Sys.time()))

# save/format resulting genes
gene_dev_df = data.frame(gene_idx = 1:length(gene_dev),
                         deviance = gene_dev,
                        gene_name = (gene_odm |> ondisc::get_feature_covariates() |> rownames())) |>
              dplyr::mutate(is_transcription_factor = gene_name %in% TF_names) 
important_genes_ranks = gene_dev_df |> 
                      filter(!is_transcription_factor) |> # not a TF gene
                      arrange(desc(deviance)) |>          # top dev genes
                      mutate(importance_rank = 1:n()) |>
                      select(gene_idx, importance_rank)
gene_dev_df = merge(gene_dev_df, important_genes_ranks, by = 'gene_idx', all.x=TRUE) |> arrange(gene_idx)

write.csv(gene_dev_df,    file = sprintf('%s/gene_deviance.csv',       save_dir), row.names = F)

write.csv(gene_dev_df |> filter(importance_rank <= NUM_IMPORTANT_GENES) |> arrange(importance_rank),    
          file = sprintf('%s/gene_deviance_topnoTFonly.csv',       save_dir), row.names = F)



# saveRDS(important_genes_idx, file = sprintf('%s/important_genes_idx.rds', save_dir))
# # save the names of the important genes
# saveRDS((gene_odm |> ondisc::get_feature_covariates() |> rownames())[important_genes_idx], 
#         sprintf('%s/important_genes_name.rds', save_dir))

# =================== END ==========================================================================
print(sprintf("[%s] END", Sys.time()))
