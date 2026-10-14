# ------------------------------------------------------------------------------------------------ #
#     Initial Estimates on Jin-2025 dataset
#     Mouse brain data from Xin Jin's lab ~ 2025 and downloaded from dropbox
# Email: Ergan to Jin-Hong to Catherine:
#     Please check the following two links of the data available.
#     Those data is shared in 2 batches by Seoyeon, where one is containing 10 perturbations and another one contains 48 perturbations.
#
#
# First in python script, subset for chosen celltype L45 Glut
# - Gene selection (do this on perturbed cells, or else large control sample will overwhelm)
# - Gene normalization
# - Pre-processing (gene importance and normalization)
# - Perform PCA? (or use pre calculated ones)
#
# ------------------------------------------------------------------------------------------------ #
args <- commandArgs(trailingOnly = TRUE)
# args = c('macbook')

assertthat::assert_that(length(args) > 0, msg = "must give arg for specifying device eg 'Rscript <filename>.R ubergenno'")
DEVICE <- args[1]
source("../PATHS.R") # load in data_dir and save_dir and CODE_DIR, depending on DEVICE value
assertthat::assert_that(!is.null(data_dir), msg = "first arg must be: laptop, desktop, or ubergenno")


dir.create(save_dir, showWarnings = FALSE)

suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(DelayedMatrixStats))
suppressPackageStartupMessages(library(sparseMatrixStats))
# library(Seurat)
# library(SeuratObject)
suppressPackageStartupMessages(library(anndataR))
suppressPackageStartupMessages(library(h5lite))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))

ggplot2::theme_set(theme_classic() +
                     theme(plot.title = element_text(hjust = .5)))

# params from python
PARAMS = list(
  'TAG'= 'L45Glut',
  'subfolder'='3_L45glutcells',
  'seed'=12345,
  'GRNA_MIN_CELLS'      = 20, # filter for grnas with at least this number of cells (keep low here, set higher later)
  'GENE_MIN_NONZERO_TRTCELLS' =3000, # filter for genes with at least this number of TREATED cells (not Non_target)
  'neuron_type' = "Glut",                     # subset for this neuron type
  'neuron_subgroup' = "005 L4-5 IT CTX Glut", # subset for this neuron subtype predicted_group=='005 L4-5 IT CTX Glut’
  'NUM_SUBSAMPLE_CONTROL_CELLS'= 10000, 
  'NUM_FAKE_GRNA'= 20,
  'FAKE_GRNA_SIZE'= 10, # number of nt grnas to constrruct fake grnas
  'FAKE_GRNA_MIN_NUMCELLS'= 100   # minimum size of fake grna
)

# parameters just for this 
params = list(
  NUM_IMPORTANT_GENES  = 4000,
  CELL_SUBSET_SIZE     = 25000,
  CELL_SUBSET_SIZE_PCA = 25000,
  RUN_EVEN_IF_SAVED    = TRUE,
  GENE_MIN_NONZERO_TRTCELLS=3000
)

save_dir <- paste0(save_dir, "/",  PARAMS[['subfolder']])
h5file <- sprintf("%s/jin2025_%s.h5", save_dir, PARAMS[['TAG']])

# =================== START ========================================================================
print(sprintf("[%s] START: processData in Jin2025", Sys.time()))


# ============================================================================== #
# =================== Load data from h5ad   ====================================
# ============================================================================== #
print(sprintf("[%s]    - load h5ad", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))

jin2025 <- anndataR::read_h5ad(sprintf("%s/jin2025_%s.h5ad", save_dir, PARAMS[['TAG']]))
metadata <- jin2025$obs




# ============================================================================== #
# =================== Construct grna and grnaNT matrix =========================
# ============================================================================== #
print(sprintf("[%s]    - construct grna and grnaNT", Sys.time()))


# construct grna: perturbation assignment for combined perturbations
grna <- model.matrix(~ gene_target_w_fakeNT - 1,
                     data = metadata |> mutate(gene_target_w_fakeNT = relevel(gene_target_w_fakeNT, ref = "Non_target")),
                     check.names = FALSE
)
colnames(grna) <- sub("^gene_target_w_fakeNT", "", colnames(grna))
grna <- t(grna) # feature x cell
grna <- as(grna, "dgRMatrix") # convert to sparse matrix


# # construct grnaNT: perturbation assignment for individual Non_target perturbations
# grnaNT <- model.matrix(~ guide_call - 1,
#                        data = metadata,
#                        check.names = FALSE
# )
# colnames(grnaNT) <- sub("^guide_call", "", colnames(grnaNT))
# grnaNT <- grnaNT[, grepl("Non_target", colnames(grnaNT))] # subset only NT guides
# grnaNT <- t(grnaNT) # feature x cell
# grnaNT <- as(grnaNT, "RsparseMatrix") # convert to sparse matrix dgRMatrix
# 
# 
# # some checks
# assertthat::assert_that(all(row.names(grna) %in% metadata$gene_target), msg = "grna matrix cols should match gene targets")
# assertthat::assert_that(all(row.names(grnaNT) %in% metadata$guide_call), msg = "grnaNT matrix cols should match guide calls")


# ============================================================================== #
# =================== Load TF and Housekeeping genes ===========================
# ============================================================================== #
print(sprintf("[%s]    - Load list of TF and houskeeping genes", Sys.time()))

# human TF list
# load list of genes that are Transcription Factors (read in xlsx sheet)
tf_raw <- readxl::read_xlsx(
  path = paste0(data_dir, "/../extra/transcriptionfactorlist.xlsx"),
  sheet = 2
) |> suppressMessages() # suppress messages on how they renamed columns
tf <- tf_raw[-1, c("...2", "Is TF?")] #      2nd col ('...2') is gene name, 4th col is TF indicator 'Is TF?'
colnames(tf) <- c("gene_name", "TF") #      First row is not data (2 rows of colnames)
humanTF <- tf |>
  dplyr::filter(TF == "Yes") |>
  dplyr::pull(gene_name) |>
  stringr::str_to_title()

# mouse TF list
mouseTF <- read.csv(paste0(data_dir, "/../extra/mousetranscriptionfactor/mouseTranscriptionFactors.txt"), header = FALSE, sep = "\t")
colnames(mouseTF) <- c("ensembl", "gene_name")
mouseTF <- mouseTF[, "gene_name"]

# housekeeping genes
housekeepgenes0 <- jsonlite::read_json(paste0(data_dir, "/../extra/housekeepinggenes/HOUNKPE_HOUSEKEEPING_GENES.v2026.1.Hs.json"))
housekeepgenes <- housekeepgenes0$HOUNKPE_HOUSEKEEPING_GENES$geneSymbols |> unlist()
housekeepgenes <- stringr::str_to_title(housekeepgenes)

rm(tf_raw, tf, housekeepgenes0)


# ============================================================================== #
# =================== Calculate Gene Dev/Importance ============================
# ============================================================================== #
print(sprintf("[%s]    - calculate gene dev/importance and choose top genes", Sys.time()))

# load if already saved
if (!params$RUN_EVEN_IF_SAVED & "gene_dev.csv" %in% list.files(save_dir)) {
  gene_dev_df <- read.csv(file = sprintf("%s/gene_dev.csv", save_dir))
} else {
  # -------- get counts of treated cells only
  trt_cellnames <- metadata |>
    filter(gene_target != "Non_target") |>
    row.names()
  if (length(trt_cellnames) > params$CELL_SUBSET_SIZE) {
    trt_cellnames_subsample <- sample(trt_cellnames, size = params$CELL_SUBSET_SIZE, replace = FALSE)
  } else {
    trt_cellnames_subsample <- trt_cellnames
  }
  counts_subset = jin2025$X[trt_cellnames_subsample, ] |> t() # needs to be genes x cells
  
  # -------- calculate gene dev on treated cells only
  gene_dev <- scry::devianceFeatureSelection(object = counts_subset, fam = "binomial") # Input should be: row = feature, col = cell
  rm(counts_subset)
  
  # -------- assemble gene_dev_df: deviance, ranking
  gene_dev_df <- data.frame(gene_idx = 1:length(gene_dev), deviance = gene_dev) |>
    dplyr::arrange(desc(gene_dev)) |>
    dplyr::mutate(importance_rank = 1:n()) # add rank importance
  gene_dev_df <- gene_dev_df[names(gene_dev), ] |> tibble::rownames_to_column(var = "gene_name") |> # reorder to original
    dplyr::mutate(is_grna_target = (gene_name %in% unique(metadata$gene_target))) # add T/F is a perturbation target
  
  
  # add columns: is transcription factor and is housekeeping gene
  gene_dev_df <- gene_dev_df |> mutate(
    is_transcription_factor = (gene_name %in% union(humanTF, mouseTF)),
    is_housekeeping = (gene_name %in% housekeepgenes)
  )
  
  
  # -------- assemble gene_dev_df: gene nonzero counts
  # add other gene info for gene qc/selection
  # gene_nonzero_cells: number of cells w non-zero counts
  # gene_nonzero_trtcells: number of treatment cells with non-zero counts
  
  #' Find the num nonzero cells for each gene, in this subset of selected cellnames
  #' requires jin2025_10pert and jin2025_48pert anndataR h5ad loaded objects
  #' @param cellnames (vector) of character cellnames according to the metadata (e.g. has the :10pert and :48pert suffix)
  calc_gene_nonzero_counts <- function(cellnames) {
    # counts_celltypegrna <- jin2025$X[cellnames, ]
    # rbind(
    #   jin2025_10pert$X[jin2025_10pert$obs_names[paste0(jin2025_10pert$obs_names, ":10pert") %in% cellnames], ],
    #   jin2025_48pert$X[jin2025_48pert$obs_names[paste0(jin2025_48pert$obs_names, ":48pert") %in% cellnames], ]
    # )
    gene_nonzero_counts <- colSums(jin2025$X[cellnames, ] != 0)
    return(gene_nonzero_counts)
  }
  
  gc()
  gene_nonzero_counts_trt <- calc_gene_nonzero_counts(metadata |> filter(gene_target != "Non_target") |> row.names())
  gene_nonzero_counts_celltypegrna <- calc_gene_nonzero_counts(metadata |> row.names())
  gc()
  
  
  # make into df and merge just to make sure the gene names match in case the order is not the same
  df1 <- data.frame(gene_nonzero_counts_trt) |>
    tibble::rownames_to_column() |>
    `colnames<-`(c("gene_name", "gene_nonzero_trtcells"))
  df2 <- data.frame(gene_nonzero_counts_celltypegrna) |>
    tibble::rownames_to_column() |>
    `colnames<-`(c("gene_name", "gene_nonzero_cells"))
  df3 <- merge(df1, df2, by = "gene_name")
  
  
  gene_dev_df <- merge(gene_dev_df, df3, by = "gene_name")
  
  # -------- assemble gene_dev_df: select normalized genes
  # genes to normalize (top X genes AND need min nonzero treatment cells) OR grna target genes OR (housekeeping genes AND min trt cells AND less strict rank)) x (all celltypegrna subset cells)
  gene_dev_df <- gene_dev_df |> mutate(is_normalized_gene = ((importance_rank <= params$NUM_IMPORTANT_GENES) & (gene_nonzero_trtcells >= params$GENE_MIN_NONZERO_TRTCELLS)) |
                                                                 is_grna_target |
                                         (is_housekeeping & (gene_nonzero_trtcells >= params$GENE_MIN_NONZERO_TRTCELLS) & (importance_rank <= 2 * params$NUM_IMPORTANT_GENES)))
  
  # reorder rows (genes) to the original again
  row.names(gene_dev_df) <- gene_dev_df$gene_name
  gene_dev_df <- gene_dev_df[jin2025$var_names, ]
  
  write.csv(x = gene_dev_df, file = sprintf("%s/gene_dev.csv", save_dir), row.names = FALSE) # csv
  rm(gene_dev, df1, df2, df3, gene_nonzero_counts_trt, gene_nonzero_counts_celltypegrna, trt_cellnames, trt_cellnames_subsample, calc_gene_nonzero_counts)
  
  # assertthat::assert_that(all(gene_dev_df$gene_nonzero_trtcells <= gene_dev_df$gene_nonzero_cells))
  # ggplot(gene_dev_df, aes(x = gene_nonzero_trtcells, gene_nonzero_cells)) + geom_point() + geom_abline(aes(intercept = 0, slope=1))
  # ggplot(gene_dev_df, aes(x = importance_rank, gene_nonzero_cells)) + geom_point()
  # ggplot(gene_dev_df, aes(x = importance_rank, gene_nonzero_trtcells, color = is_grna_target)) + geom_point()
  # ggplot(gene_dev_df, aes(x = importance_rank, gene_nonzero_trtcells, color = is_normalized_gene)) + geom_point()
  # ggplot(gene_dev_df, aes(x = importance_rank, gene_nonzero_trtcells, color = is_transcription_factor)) + geom_point()
  # ggplot(gene_dev_df, aes(x = deviance, y=gene_nonzero_trtcells, color = is_normalized_gene)) + geom_point() + coord_cartesian(xlim=c(0, 50000))
}


# ============================================================================== #
# =================== Save as h5 files =========================================
# ============================================================================== #
print(sprintf("[%s]    - save as h5", Sys.time()))

# save the norm and counts and metadata in a sep h5 file becuse that is easier to load



# metadata about this subtype of cell + subset grna
print(sprintf("[%s]        + cell_metadata", Sys.time()))

cell_metadata <- metadata
for (i in 1:ncol(cell_metadata)) { # some empty '' cols creating issues when saving. replace NA here
  if (all(cell_metadata[, i] == "")) {
    cell_metadata[, i] <- NA
  }
}
# colnames(cell_metadata) = sapply(colnames(cell_metadata), FUN = function(x){gsub('\\.', '_', x)}) |> unname() # replace . in colnames to _ bc character issues when saving? this is ok
h5lite::h5_write(cell_metadata, h5file, "cell_metadata/cell_metadata")


# grna assignments
print(sprintf("[%s]        + grna", Sys.time()))
h5lite::h5_write(as.matrix(grna), h5file, "grna/grna")
h5lite::h5_write(row.names(grna), h5file, "grna/row_grnanames")
h5lite::h5_write(colnames(grna), h5file, "grna/col_cellnames")

# # grna NT assignments
# print(sprintf("[%s]        + grnaNT", Sys.time()))
# h5lite::h5_write(as.matrix(grnaNT), h5file, "grnaNT/grnaNT")
# h5lite::h5_write(row.names(grnaNT), h5file, "grnaNT/row_grnanames")
# h5lite::h5_write(colnames(grnaNT), h5file, "grnaNT/col_cellnames")


# test loading in saved files
if (F) {
  mat_in_meta   <- h5_read(h5file, "cell_metadata/cell_metadata") # R reads in row and colnames
  mat_in_grna   <- h5_read(h5file, "grna/grna") # R reads in row and colnames
  mat_in_grnaNT <- h5_read(h5file, "grnaNT/grnaNT") # R reads in row and colnames
  
  assertthat::assert_that(all(colnames(mat_in_grnaNT) == colnames(mat_in_grna)))
  assertthat::assert_that(all(colnames(mat_in_grnaNT) == row.names(mat_in_meta)))
  assertthat::assert_that(all(row.names(mat_in_meta) == colnames(mat_in_grna)))
  
  rm(mat_in_meta, mat_in_grna, mat_in_grnaNT)
}







# h5_write(as.matrix(gene_counts), h5file, "counts/counts")        # it seems that python has trouble reading in with extra rownames and columns...
# h5_write(gene_counts, h5file, "counts/counts")
# h5_write(row.names(gene_counts), h5file, "counts/row_genenames") # so create hierarchy and save rownames and colnames separately
# h5_write(colnames(gene_counts), h5file, "counts/col_cellnames")  



# ============================================================================== #
# =================== gather gene_counts =======================================
# ============================================================================== #



normalized_genes_names = gene_dev_df |> filter(is_normalized_gene) |> pull(gene_name) |> sort()
# cell_names = colnames(grna) # should be the same set and ordering as colnames(grnaNT) and row.names(cell_metadata)
# 
# # subset each of the datasets and rename cells to indicate source (:10pert or :48pert)
# gene_counts_10 = jin2025_10pert$X[, normalized_genes_names]
# row.names(gene_counts_10) = paste0(row.names(gene_counts_10), ':10pert')
# gene_counts_48 = jin2025_48pert$X[, normalized_genes_names]
# row.names(gene_counts_48) = paste0(row.names(gene_counts_48), ':48pert')
# 
# gene_counts = rbind(gene_counts_10, gene_counts_48)
# gene_counts = t(gene_counts)
# gene_counts = gene_counts[, cell_names]
# 
# rm(gene_counts_10, gene_counts_48)
# gc(verbose = FALSE)


gene_counts = jin2025$X[, normalized_genes_names] |> t()

# ============================================================================== #
# =================== normalize to gene_norm ===================================
# ============================================================================== #


# remove all non-relevant objects in environment
rm(list = setdiff(objects(),
                  c('gene_counts', 'PARAMS', 'params', 'gene_norm',
                    'DEVICE', 'save_dir', 'util_dir', 'data_dir', 'h5file')))


# gc()
gene_counts <- as(gene_counts, "HDF5Matrix") # faster to first save gene_counts, and then load as delayed array. no, need to convert for gene_norm later if not here, to save gene
# gc()
# gene_counts = HDF5Array::HDF5Array(hdf5_filename, 'gene_counts/gene_counts')

# =================== normalize
# genes in the rows and samples in the columns
gene_norm = scry::nullResiduals(gene_counts,
                                fam  = "binomial",
                                type = "deviance") #  genes x  cells


# =================== save to same h5 file as metadata and grnas
h5file
rhdf5::h5ls(h5file)

# gene_counts
h5lite::h5_delete(file=h5file, name='gene_counts/', attr = NULL, warn = TRUE) # delete current h5 groupdata if already present
rhdf5::h5createGroup(h5file, group = 'gene_counts') # create group
HDF5Array::writeHDF5Array(gene_counts, filepath = h5file, 
                          name = 'gene_counts/gene_counts', 
                          verbose=TRUE, with.dimnames = TRUE)

# gene_norm
h5lite::h5_delete(file=h5file, name='gene_norm/', attr = NULL, warn = TRUE) # delete current h5 groupdata if already present
rhdf5::h5createGroup(h5file, group = 'gene_norm') # create group
HDF5Array::writeHDF5Array(gene_norm, filepath = h5file, 
                          name = 'gene_norm/gene_norm', 
                          verbose=TRUE, with.dimnames = TRUE)




assertthat::assert_that(all(colnames(gene_counts) == colnames(gene_norm)))

# ============================================================================== #
# =================== PCA ======================================================
# ============================================================================== #
# remove all non-relevant objects in environment
rm(list = setdiff(objects(),
                  c('params', 'gene_norm',
                    'DEVICE', 'save_dir', 'util_dir', 'data_dir', 'h5file')))


set.seed(12345)

# sample treated cells more
grna_h5   <- h5lite::h5_read(h5file, "grna/grna") # R reads in row and colnames
gene_norm = HDF5Array::HDF5Array(h5file, name='gene_norm/gene_norm')
row.names(gene_norm) = rhdf5::h5read(h5file, name = 'gene_norm/.gene_norm_dimnames/1') |> as.vector()
colnames( gene_norm) = rhdf5::h5read(h5file, name = 'gene_norm/.gene_norm_dimnames/2') |> as.vector()
# 
# trt_cells = which(grna_h5['Non_target', ] == 0) |> names()
# NT_cells = which(grna_h5['Non_target', ] == 1) |> names()
# subsample_cells = c(sample(trt_cells, size = min(length(trt_cells), round(.8*params$CELL_SUBSET_SIZE_PCA))) ,
#                     sample(NT_cells,  size = min(length(NT_cells), round(.2*params$CELL_SUBSET_SIZE_PCA))))

gene_dev_df = read.csv(sprintf('%s/gene_dev.csv', save_dir))
pca_genes = setdiff(row.names(gene_norm), gene_dev_df |> filter(is_grna_target | is_transcription_factor) |> pull(gene_name))



# === perform only on Non_target cells
gene_norm_PCA = gene_norm[pca_genes,  names(which(grna_h5['Non_target', ] == 1))] |> t()

# perform pca
pca_res = prcomp(x = gene_norm_PCA, rank.=50, retx = TRUE, center = TRUE, scale. = TRUE)

# construct NCs as the actual loadings: [all cells x valid genes] %*% [rotation of valid genes x PCs]
# NC_loadings = t(gene_norm) %*% pca_res$rotation
NC_loadings = scale(t(gene_norm[pca_genes, ]), center=TRUE,scale=TRUE) %*% pca_res$rotation # should center/scale first oops, but they are linearlly corr so i don't think it matters much hopefully

# or center and scale based on subsampled data? it does not matter, these lie on a line basically
# NC_loadings2 = scale(t(gene_norm), center=pca_res$center,scale=pca_res$scale) %*% pca_res$rotation # should center/scale first oops, but they are linearlly corr so i don't think it matters much hopefully
# ggplot(NULL, aes(x = NC_loadings[, 1], y = NC_loadings2[, 1])) + geom_point() + coord_cartesian(ylim = c(-10, 10), xlim = c(-10, 10))
# ggplot(NULL, aes(x = NC_loadings[, 2], y = NC_loadings2[, 2])) + geom_point() + coord_cartesian(ylim = c(-10, 10), xlim = c(-10, 10))


# save
dir.create(sprintf('%s/pca/', save_dir), showWarnings = FALSE)
saveRDS(NC_loadings, sprintf('%s/pca/NCloadings.rds', save_dir)) 
saveRDS(    pca_res, sprintf('%s/pca/pca_res.rds', save_dir)) 


# ===  perform on all cells
gene_norm_PCA = gene_norm[pca_genes, ] |> t() # just use all cells i guess (no subsampling)
pca_res = prcomp(x = gene_norm_PCA, rank.=50, retx = TRUE, center = TRUE, scale. = TRUE) # perform pca

# construct NCs as the actual loadings: [all cells x valid genes] %*% [rotation of valid genes x PCs]
NC_loadings = scale(t(gene_norm[pca_genes, ]), center=TRUE,scale=TRUE) %*% pca_res$rotation # should center/scale first oops, but they are linearlly corr so i don't think it matters much hopefully

saveRDS(NC_loadings, sprintf('%s/pca/NCloadings_mixTrtUntrtCells.rds', save_dir)) 
saveRDS(    pca_res, sprintf('%s/pca/pca_res_mixTrtUntrtCells.rds', save_dir)) 





# checking pca results
all(row.names(NC_loadings) == colnames(grna_h5))
all(row.names(NC_loadings) == colnames(gene_norm))
all(colnames(grna_h5) == colnames(gene_norm))

cell_metadata = h5lite::h5_read(h5file, "cell_metadata/cell_metadata")
all(row.names(cell_metadata) == colnames(grna_h5))

cell_metadata$PC1 = NC_loadings[, 1]
cell_metadata$PC2 = NC_loadings[, 2]
cell_metadata$PC3 = NC_loadings[, 3]

colnames(cell_metadata)
ggplot(cell_metadata |> filter(gene_target != 'Non_target'), aes(x = PC1, y = PC2, color =  predicted_subclass)) + geom_point(alpha = .3) + theme(legend.position = 'bottom')
ggplot(cell_metadata |> filter(gene_target != 'Non_target'), aes(x = PC1, y = PC2, color =  gene_target)) + geom_point(alpha = .3) + theme(legend.position = 'bottom')



umap_res = umap::umap(NC_loadings)
cell_metadata$UMAP1 = umap_res$layout[, 1]
cell_metadata$UMAP2 = umap_res$layout[, 2]
ggplot(cell_metadata, # |> filter(gene_target != 'Non_target'), 
       aes(x = UMAP1, y = UMAP2, color =  predicted_subclass)) + geom_point(alpha = .3) + theme(legend.position = 'bottom')
ggplot(cell_metadata |> filter(gene_target != 'Non_target'), 
       aes(x = UMAP1, y = UMAP2, color =  gene_target)) + geom_point(alpha = .3) + theme(legend.position = 'bottom')


ggplot(cell_metadata, aes(x = predicted_group_))










ggplot(NULL, aes(x = 1:200, y = pca_res$sdev[1:200])) + geom_point() + geom_line() + coord_cartesian(ylim = c(0, 12))

plot(NC_loadings[sample(nrow(NC_loadings), 1000), ], aes(x=PC1, y=PC3))+geom_point()






cell_metadata$PC1 = NC_loadings[, 1]
cell_metadata$PC2 = NC_loadings[, 2]
cell_metadata$PC3 = NC_loadings[, 3]

colnames(cell_metadata)
ggplot(cell_metadata |> filter(gene_target != 'Non_target'), aes(x = PC1, y = PC2, color =  predicted_subclass)) + geom_point(alpha = .3) + theme(legend.position = 'bottom')
ggplot(cell_metadata |> filter(gene_target != 'Non_target'), aes(x = PC1, y = PC2, color =  gene_target)) + geom_point(alpha = .3) + theme(legend.position = 'bottom')



# compare NC loadings
cell_metadata = h5lite::h5_read(h5file, "cell_metadata/cell_metadata")
NC_loadings_cntrlcells = readRDS(sprintf('%s/pca/NCloadings.rds', save_dir))
NC_loadings_mixcells   = readRDS(sprintf('%s/pca/NCloadings_mixTrtUntrtCells.rds', save_dir))



all(row.names(cell_metadata) == row.names(NC_loadings_cntrlcells))
all(row.names(cell_metadata) == row.names(NC_loadings_mixcells))


cor(NC_loadings_cntrlcells[, 1:5]) |> round(3) # not exactly 0s off diagonal bc subsampling
cor(NC_loadings_mixcells[, 1:5]) |> round(3)

plot_comparePCAres <- function(i) {
  ggplot(NULL, aes(x = NC_loadings_cntrlcells[, i], y = NC_loadings_mixcells[, i])) + 
    geom_abline(aes(intercept=0, slope=1)) + geom_point() + labs(title = paste0('PC ', i), x='PCA from Control Cells', y = 'PCA from Mix of Trt and Cntrl Cells')
}

g = list()
for(i in 1:12) {
  g[[i]] = plot_comparePCAres(i)
}
gridExtra::grid.arrange(grobs = g, ncol=3)




# ============================================================================== #
# =================== END ======================================================
# ============================================================================== #
print(sprintf("[%s] END", Sys.time()))






