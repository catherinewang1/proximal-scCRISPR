# ------------------------------------------------------------------------------------------------ #
#                        Normalize genes in papalexi dataset
#     Normalizes the papalexi dataset (processed by Tim Barry) with Will Townes'
#     method by fitting a multinomial model and using the glm deviance as the
#     normalized values (..., fam  = "binomial", type = "deviance")
# Requires: papalexi dataset saved in ondisc format at <data_dir>/papalexi-2021/
# Outputs: (nothing) but saves loaded in and normalized gene expression at
#          <data_dir>/papalexi_saves/gene.h5    (in HDF5 format)
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

assertthat::assert_that(length(args) >= 2, msg="provide device and NUM_IMPORTANT_GENES (currently int <=4000)")
NUM_IMPORTANT_GENES = as.integer(args[2])


# =================== Start ========================================================================
print(sprintf("[%s] START: Normalize papalexi (top %d genes)", Sys.time(), NUM_IMPORTANT_GENES))

gene_odm <- ondisc::read_odm(odm_fp      = paste0(data_dir, "/papalexi-2021/processed/gene/expression_matrix.odm"),
                             metadata_fp = paste0(data_dir, "/papalexi-2021/processed/gene/metadata.rds"))
grna_odm <- ondisc::read_odm(odm_fp      = paste0(data_dir, "/papalexi-2021/processed/grna_assignment/assignment_matrix.odm"),
                             metadata_fp = paste0(data_dir, "/papalexi-2021/processed/grna_assignment/metadata.rds"))

important_genes_idx   = readRDS(sprintf('%s/important_genes_idx.rds', save_dir))

invisible(gc(verbose=FALSE))


# =================== Convert gene exp to HDF5 file by converting in memory ========================
print(sprintf("[%s]    - converting to HDF5 file", Sys.time()))
# will be saved in...
h5file   = paste0(save_dir, "/gene.h5")
HDF5Array::setHDF5DumpFile(h5file)

# delete existing 'gene' if already exists
h5f         = rhdf5::H5Fopen(h5file); 
name_exists = rhdf5::H5Lexists(h5f, 'gene')
rhdf5::h5closeAll()
if (name_exists) {  rhdf5::h5delete(file = h5file, name = 'gene') }
rm(h5f, name_exists)

# save gene expression of top genes in 'gene.h5' under 'gene'
HDF5Array::setHDF5DumpName('gene')

gene_important = gene_odm[[important_genes_idx[1:NUM_IMPORTANT_GENES],
                           1:ncol(gene_odm)]] # load top gene expr into memory (NUM_IMPORTANT_GENES x #cells)
gene_hd5 = as(gene_important, "HDF5Matrix")
rhdf5::h5closeAll()
invisible(gc(verbose=FALSE))


# =================== Perform normalization ========================================================
print(sprintf("[%s]    - performing normalization", Sys.time()))
gene_norm = scry::nullResiduals(gene_hd5, 
                                fam  = "binomial", 
                                type = "deviance") # #imp genes x 20729


# =================== save normalized gene expression ==============================================
print(sprintf("[%s]    - saving normalized gene expression", Sys.time()))

# delete existing gene_norm if already exists
h5f         = rhdf5::H5Fopen(h5file); 
name_exists = rhdf5::H5Lexists(h5f, 'gene_norm')
rhdf5::h5closeAll()
if (name_exists) {  rhdf5::h5delete(file = h5file, name = 'gene_norm') }
rm(h5f, name_exists)

# save normalized gene expression of top genes in 'gene.h5' under 'gene_norm'
if(DEVICE == 'ubergenno') { #if on ubergenno, need to realize to memory first?????
    HDF5Array::writeHDF5Array(DelayedArray::realize(gene_norm), filepath=h5file, name='gene_norm')
} else {
    HDF5Array::writeHDF5Array(gene_norm, filepath=h5file, name='gene_norm')
}

invisible(gc(verbose=FALSE))


# =================== END ==========================================================================
print(sprintf("[%s] END", Sys.time()))
