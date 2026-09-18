# ------------------------------------------------------------------------------------------------ #
# Save a copy of the grna assignment in a more compatible format for other analysis
# (e.g. as a separate csv, or as a h5 file)
# (just do this locally)
#
# ------------------------------------------------------------------------------------------------ #
# args = commandArgs(trailingOnly = TRUE)
args = c('macbook')

suppressPackageStartupMessages(require(assertthat)) # for some assert statements
suppressPackageStartupMessages(require(ondisc))     # loading in data
suppressPackageStartupMessages(require(rhdf5))      # read/save in HDF5 format
suppressPackageStartupMessages(require(HDF5Array))  # (^same)

assertthat::assert_that(length(args) > 0, msg="must give arg for specifying device eg 'Rscript <filename>.R ubergenno'")
DEVICE = args[1]
source('../PATHS.R') # load in data_dir and save_dir, depending on DEVICE value


# =================== Start ========================================================================
print(sprintf("[%s] START: converge grna odm to h5 and csv format", Sys.time()))

gene_odm <- ondisc::read_odm(odm_fp      = paste0(data_dir, "/papalexi-2021/processed/gene/expression_matrix.odm"),
                             metadata_fp = paste0(data_dir, "/papalexi-2021/processed/gene/metadata.rds"))
grna_odm <- ondisc::read_odm(odm_fp      = paste0(data_dir, "/papalexi-2021/processed/grna_assignment/assignment_matrix.odm"),
                             metadata_fp = paste0(data_dir, "/papalexi-2021/processed/grna_assignment/metadata.rds"))





# load all into memory
# gene <- gene_odm[[,1:ncol(gene_odm)]] # entire gene ds (18649 genes x 20729 cells)
grna <- grna_odm[[,1:ncol(grna_odm)]] # entire grna ds (110   grnas x 20729 cells)
row.names(grna) = grna_odm@feature_covariates |> row.names()

# =================== Convert grna to HDF5 file ========================
print(sprintf("[%s]    - converting to HDF5 file", Sys.time()))
# will be saved in...
h5file   = paste0(save_dir, "/gene.h5")
HDF5Array::setHDF5DumpFile(h5file)

# delete existing 'grna' if already exists
h5f         = rhdf5::H5Fopen(h5file); 
name_exists = rhdf5::H5Lexists(h5f, 'grna')
rhdf5::h5closeAll()
if (name_exists) {  rhdf5::h5delete(file = h5file, name = 'grna') }
rm(h5f, name_exists)

# save grna expression of top genes in 'gene.h5' under 'grna'
HDF5Array::setHDF5DumpName('grna')
grna_hd5 = as(grna, "HDF5Matrix")
rhdf5::h5closeAll()
invisible(gc(verbose=FALSE))


# =================== Format papalexi-2021 data for causarray ========================
# take the tests chosen in each setting
#
# filter out for only grnas w the preselected threshold (typically 100)
# and then only select the cells for those

# AYZW_setting_name = 'A'
# 
# # load chosen AYZW names
# AY   = read.csv(sprintf('%s/AY/%s/AY.csv', save_dir, AYZW_setting_name))
# 
# # only take the grnas with a test associated (which means it passed qc) or a NT
# grna_names = AY$A |> unique()
# grna_names = grna_names[grna_names %in% row.names(grna)] |> sort()
# NT_names = grep('NTg', row.names(grna), value = TRUE)
# 
# 
# # only the cells with one of these or a nt perturbation
# cell_subset = which(colSums(grna[c(NT_names, grna_names), ]) == 1)
# gene_subset = gene[, cell_subset]
# grna_subset = grna[c(NT_names, grna_names), cell_subset]

# =================== END ==========================================================================
print(sprintf("[%s] END", Sys.time()))



'PDL1g1' %in% AY$A



