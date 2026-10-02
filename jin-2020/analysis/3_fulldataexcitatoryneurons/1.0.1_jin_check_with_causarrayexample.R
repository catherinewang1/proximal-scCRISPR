# ------------------------------------------------------------------------------------------------ #
#    For the analysis on jin-2020 dataset (excitatory cells chosen from Broad website and counts from GEO)
#    check that the counts match with the subset from causarray example
# 
# ------------------------------------------------------------------------------------------------ #
args = commandArgs(trailingOnly = TRUE)
# args = c('macbook')


assertthat::assert_that(length(args) > 0, msg="must give arg for specifying device eg 'Rscript <filename>.R ubergenno'")
DEVICE = args[1]
source('../PATHS.R') # load in data_dir and save_dir and CODE_DIR, depending on DEVICE value
assertthat::assert_that(!is.null(data_dir), msg='first arg must be: laptop, desktop, or ubergenno')

library(dplyr)
library(Seurat)
library(SeuratObject)


# =================== LOAD: GEO ========================================================================
print(sprintf("[%s] LOAD: Jin Full Dataset from GEO, subselect excitatory neurons", Sys.time()))

# =================== \__grna: carried over from Broad ===========================
grna = read.csv(sprintf('%s/grna.csv', save_dir), row.name = 'X') # for cells w a pert assignment 22085 x 39
# grna = t(grna) #  30 perturbations x 4744 cells
# row.names(grna) = stringr::str_to_title(row.names(grna)) # hopefully these match broad Perturbation names now



# =================== \__counts: processed from GEO ===========================
# load in top 4000 genes (count and continuous)
gene = readRDS(sprintf('%s/gene_counts.rds', save_dir))
gene = t(gene) # so cell x gene
gene = as.matrix(gene)


# =================== \__rename: ===========================
Y_geo = gene 
A_geo = grna |> select(-GFP)

rm(gene, grna)

# =================== \__gene_dev ===========================
gene_dev = read.csv(sprintf('%s/gene_dev.csv', save_dir))

# =================== LOAD: causarray example========================================================================
print(sprintf("[%s] LOAD: causarray example data", Sys.time()))

# =================== \__from causarray tutorial by Jinhong ========================================================================
# from tutorial
# A: cell-by-perturbation 0/1 trt assignment
#    - cells that received negative control perturbation GFP have 0 in all cols
# Y: cell-by-gene matrix of gene expression counts

sc.seurat <- readRDS(sprintf('%s/causarrayexampledata/perturbseq-exneu.rds', data_dir))

# Access the counts through Seurat when its installed version recognizes the
# serialized Assay5 object; otherwise recover the same layer and dimnames
# directly. The fallback keeps this tutorial runnable with older Seurat builds.
if ("RNA" %in% Assays(sc.seurat)) {
  counts <- GetAssayData(sc.seurat, assay = "RNA", layer = "counts")
} else {
  rna.assay <- sc.seurat@assays[["RNA"]]
  counts <- rna.assay@layers[["counts"]]
  rownames(counts) <- rownames(rna.assay@features)
  colnames(counts) <- rownames(rna.assay@cells)
  rm(rna.assay)
}
Y <- data.frame(t(as.matrix(counts)), check.names = FALSE) # cell-by-gene matrix
metadata <- sc.seurat@meta.data

perturb <- metadata
colnames(perturb) <- gsub("Perturbation", "trt_", colnames(perturb))
perturb$trt_ <- relevel(as.factor(perturb$trt_), ref = "GFP")
A <- data.frame(
  model.matrix(~ trt_ - 1, data = perturb)[, -1, drop = FALSE],
  check.names = FALSE
) # cell-by-trt matrix; remove the first (GFP control) column
colnames(A) <- sub("^trt_", "", colnames(A))


# =================== \__rename: ===========================
Y_causarray = Y
A_causarray = A

rm(counts, sc.seurat, metadata, A, Y, perturb)


# =================== COMPARE: ========================================================================

# causarray: 2926 cells receiving 29 perturbations
# geo: 4744 cells receiving 29 perturbations

# =================== \__grna/perturbation assignment ===========================
dim(A_geo); dim(A_causarray)  # different number of cells somehow (more strict qc)
dim(A_geo) == dim(A_causarray) 
all(colnames(A_geo) == colnames(A_causarray)) # same colnames(perturbation names) and order



causarray_cellnames = row.names(A_causarray)

# are all the causarray cells are in geo?
all(causarray_cellnames %in% row.names(A_geo)) 


A_geo_subset  = A_geo[causarray_cellnames, ] # subset just the cells also in causarray

# make sure cells (rows) are in the same order
all(row.names(A_geo_subset) == row.names(A_causarray))

# make sure all the perturbation assignments are the same
all(A_geo_subset == A_causarray) 


# =================== \__gene counts ===========================


dim(Y_geo); dim(Y_causarray)

# causarray grna matches gene counts
all(row.names(A_causarray) == row.names(Y_causarray))

# are all the causarray cells are in geo?
all(row.names(Y_causarray) %in% row.names(Y_geo)) 

# are all the causarray genes are in geo?
all(colnames(Y_causarray) %in% colnames(Y_geo)) # no?



(colnames(Y_causarray) %in% colnames(Y_geo)) |> table() # 165 no and 3056 yes
colnames(Y_causarray)[!(colnames(Y_causarray) %in% colnames(Y_geo))] # which ones
gene_dev |> filter(gene_name %in% colnames(Y_causarray)[!(colnames(Y_causarray) %in% colnames(Y_geo))]) # ranks and other gene information of excluded
gene_dev |> filter(gene_name %in% colnames(Y_causarray)[!(colnames(Y_causarray) %in% colnames(Y_geo))]) |> pull(importance_rank) |> hist() # these genes just had more than >4000 importance rank

gene_intersect = intersect(colnames(Y_geo), colnames(Y_causarray)) # genes in both datasets
Y_geo_subset       =       Y_geo[row.names(Y_causarray), gene_intersect] # subset just the cells and genes also in causarray
Y_causarray_subset = Y_causarray[                      , gene_intersect]


dim(Y_geo_subset); dim(Y_causarray_subset)  
all(dim(Y_geo_subset) == dim(Y_causarray_subset)) # same dimensions?
all(row.names(Y_geo_subset) == row.names(Y_causarray_subset)) # same cell (row) names?
all(colnames(Y_geo_subset) == colnames(Y_causarray_subset)) # same gene (col) names?



all(Y_geo_subset == Y_causarray_subset) # same counts


# # test if i mess up a count
# Y_geo_subset[3, 8] = 999
# all(Y_geo_subset == Y_causarray_subset) 



# =================== \__which cells excluded ===========================










