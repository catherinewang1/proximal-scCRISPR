# ------------------------------------------------------------------------------------------------ #
#     Initial Estimates on Jin-2025 dataset
#     Mouse brain data from Xin Jin's lab ~ 2025 and downloaded from dropbox
# Email: Ergan to Jin-Hong to Catherine:
#     Please check the following two links of the data available.
#     Those data is shared in 2 batches by Seoyeon, where one is containing 10 perturbations and another one contains 48 perturbations.
#
#
# in 1.1
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
save_dir <- paste0(save_dir, "/1_initial")
assertthat::assert_that(!is.null(data_dir), msg = "first arg must be: laptop, desktop, or ubergenno")


dir.create(save_dir)

library(Matrix)
# library(Seurat)
# library(SeuratObject)
library(anndataR)
library(h5lite)
library(dplyr)
library(ggplot2)
ggplot2::theme_set(theme_classic() +
  theme(plot.title = element_text(hjust = .5)))


params <- list(
  NUM_IMPORTANT_GENES = 4000,
  CELL_SUBSET_SIZE = 25000,
  RUN_EVEN_IF_SAVED = TRUE,
  GRNA_MIN_CELLS = 250, # filter for grnas with at least this number of cells
  GENE_MIN_NONZERO_TRTCELLS = 10000, # filter for genes with at least this number of TREATED cells (not Non_target)
  neuron_type = "Glut" # subset for this neuron type
)


# =================== START ========================================================================
print(sprintf("[%s] START: processData in Jin2025", Sys.time()))


# ============================================================================== #
# =================== Load data from h5ad   ====================================
# ============================================================================== #
print(sprintf("[%s]    - load h5ad", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
# load h5ad (simpler?)
t0 <- Sys.time()
jin2025_10pert <- anndataR::read_h5ad(sprintf("%s/subset_adata_SCARF_10celltypes_10perturbations.h5ad", data_dir))
t1 <- Sys.time()
jin2025_48pert <- anndataR::read_h5ad(sprintf("%s/subset_adata_SCARF_10celltypes_48perturbations.h5ad", data_dir))
t2 <- Sys.time()
print(sprintf(
  "h5ad plain load time for is %2.2f and %2.2f mins for 10 and 48 pert, respectively",
  difftime(t1, t0, units = "mins"), difftime(t2, t1, units = "mins")
)) #  0.19 and 0.23 min
rm(t0, t1, t2)
# format(object.size(jin2025_10pert$X), units = 'Gb') # 6.9 Gb
# format(object.size(jin2025_48pert$X), units = 'Gb') # 8.9 Gb

metadata_10 <- jin2025_10pert$obs |> mutate(dataset = "10pert")
metadata_48 <- jin2025_48pert$obs |> mutate(dataset = "48pert")
row.names(metadata_10) <- paste0(row.names(metadata_10), ":10pert")
row.names(metadata_48) <- paste0(row.names(metadata_48), ":48pert")
metadata_all <- rbind(metadata_10, metadata_48)
rm(metadata_10, metadata_48)
rm(metadata_10, metadata_48)


# ============================================================================== #
# =================== Choose GRNAs/Perturbations ===============================
# ============================================================================== #
# metadata_all$gene_target |> table() |> sort()
# metadata_all |> filter(gene_target == 'Non_target') |> pull(guide_call) |> table() |> sort()
#
# metadata_subcelltype$gene_target |> table() |> sort()
# ggplot() + geom_bar(aes(x=metadata_subcelltype$gene_target)) + coord_cartesian(ylim=c(0, 800))

metadata_subcelltype <- metadata_all |> filter(neuron_type == params$neuron_type) # subset cells for this neuron type only

grna_counts <- table(metadata_subcelltype$gene_target)
selected_grnas <- names(grna_counts)[grna_counts > params$GRNA_MIN_CELLS] |>
  setdiff("Non_target") |>
  sort()


metadata_subcelltypegrna <- metadata_subcelltype |> filter(gene_target %in% selected_grnas | gene_target == "Non_target")

rm(grna_counts, metadata_subcelltype)


# ============================================================================== #
# =================== Construct grna and grnaNT matrix =========================
# ============================================================================== #
print(sprintf("[%s]    - construct grna and grnaNT", Sys.time()))


# construct grna: perturbation assignment for combined perturbations
grna <- model.matrix(~ gene_target - 1,
  data = metadata_subcelltypegrna |> mutate(gene_target = relevel(gene_target, ref = "Non_target")),
  check.names = FALSE
)
colnames(grna) <- sub("^gene_target", "", colnames(grna))
grna <- t(grna) # feature x cell
grna <- as(grna, "dgRMatrix") # convert to sparse matrix


# construct grnaNT: perturbation assignment for individual Non_target perturbations
grnaNT <- model.matrix(~ guide_call - 1,
  data = metadata_subcelltypegrna,
  check.names = FALSE
)
colnames(grnaNT) <- sub("^guide_call", "", colnames(grnaNT))
grnaNT <- grnaNT[, grepl("Non_target", colnames(grnaNT))] # subset only NT guides
grnaNT <- t(grnaNT) # feature x cell
grnaNT <- as(grnaNT, "dgRMatrix") # convert to sparse matrix


# some checks
assertthat::assert_that(all(row.names(grna) %in% metadata_subcelltypegrna$gene_target), msg = "grna matrix cols should match gene targets")
assertthat::assert_that(all(row.names(grnaNT) %in% metadata_subcelltypegrna$guide_call), msg = "grnaNT matrix cols should match guide calls")


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
# =================== Calculate Gene Dev/Importance ===========================
# ============================================================================== #
print(sprintf("[%s]    - calculate gene dev/importance and choose top genes", Sys.time()))

# load if already saved
if (!params$RUN_EVEN_IF_SAVED & "gene_dev.csv" %in% list.files(save_dir)) {
  gene_dev_df <- read.csv(file = sprintf("%s/gene_dev.csv", save_dir))
} else {
  # -------- get counts of treated cells only
  trt_cellnames <- metadata_subcelltypegrna |>
    filter(gene_target != "Non_target") |>
    row.names()
  if (length(trt_cellnames) > params$CELL_SUBSET_SIZE) {
    trt_cellnames_subsample <- sample(trt_cellnames, size = params$CELL_SUBSET_SIZE, replace = FALSE)
  } else {
    trt_cellnames_subsample <- trt_cellnames
  }
  counts_subset <- rbind(
    jin2025_10pert$X[jin2025_10pert$obs_names[paste0(jin2025_10pert$obs_names, ":10pert") %in% trt_cellnames_subsample], ],
    jin2025_48pert$X[jin2025_48pert$obs_names[paste0(jin2025_48pert$obs_names, ":48pert") %in% trt_cellnames_subsample], ]
  )
  counts_subset <- t(counts_subset) # cell x gene

  # dim(counts_subset) # 19070  genes x 25000 subsampled cells

  # -------- calculate gene dev on treated cells only
  gene_dev <- scry::devianceFeatureSelection(object = counts_subset, fam = "binomial") # Input should be: row = feature, col = cell
  rm(counts_subset)

  # -------- assemble gene_dev_df: deviance, ranking
  gene_dev_df <- data.frame(gene_idx = 1:length(gene_dev), deviance = gene_dev) |>
    dplyr::arrange(desc(gene_dev)) |>
    dplyr::mutate(importance_rank = 1:n()) # add rank importance
  gene_dev_df <- gene_dev_df[names(gene_dev), ] |> tibble::rownames_to_column(var = "gene_name") |> # reorder to original
    dplyr::mutate(is_grna_target = (gene_name %in% unique(metadata_all$gene_target))) # add T/F is a perturbation target


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
    counts_celltypegrna <- rbind(
      jin2025_10pert$X[jin2025_10pert$obs_names[paste0(jin2025_10pert$obs_names, ":10pert") %in% cellnames], ],
      jin2025_48pert$X[jin2025_48pert$obs_names[paste0(jin2025_48pert$obs_names, ":48pert") %in% cellnames], ]
    )
    gene_nonzero_counts <- colSums(counts_celltypegrna != 0)
    return(gene_nonzero_counts)
  }

  gc()
  gene_nonzero_counts_trt <- calc_gene_nonzero_counts(metadata_subcelltypegrna |> filter(gene_target != "Non_target") |> row.names())
  gene_nonzero_counts_celltypegrna <- calc_gene_nonzero_counts(metadata_subcelltypegrna |> row.names())
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
  gene_dev_df <- gene_dev_df[jin2025_48pert$var_names, ]

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
# =================== Save as h5 files ===========================
# ============================================================================== #
print(sprintf("[%s]    - save as h5", Sys.time()))


h5file <- sprintf("%s/jin2025_subset.h5", save_dir)


# metadata about this subtype of cell + subset grna
print(sprintf("[%s]        + cell_metadata", Sys.time()))

cell_metadata <- metadata_subcelltypegrna
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

# grna NT assignments
print(sprintf("[%s]        + grnaNT", Sys.time()))
h5lite::h5_write(as.matrix(grnaNT), h5file, "grnaNT/grnaNT")
h5lite::h5_write(row.names(grnaNT), h5file, "grnaNT/row_grnanames")
h5lite::h5_write(colnames(grnaNT), h5file, "grnaNT/col_cellnames")


# test loading in saved files
if (F) {
  mat_in <- h5_read(h5file, "cell_metadata/cell_metadata") # R reads in row and colnames
  mat_in <- h5_read(h5file, "grna/grna") # R reads in row and colnames
  mat_in <- h5_read(h5file, "grnaNT/grnaNT") # R reads in row and colnames
}
