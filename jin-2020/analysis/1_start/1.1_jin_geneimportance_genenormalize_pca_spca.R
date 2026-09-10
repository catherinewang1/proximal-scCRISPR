# ------------------------------------------------------------------------------------------------ #
#     Analysis on jin-2020 dataset subset from causarrayexampledata 
# 
# - Pre-processing (gene importance and normalization)
# - Perform PCA + SPCA 
# ------------------------------------------------------------------------------------------------ #
args = commandArgs(trailingOnly = TRUE)
# args = c('macbook')


assertthat::assert_that(length(args) > 0, msg="must give arg for specifying device eg 'Rscript <filename>.R ubergenno'")
DEVICE = args[1]
source('../PATHS.R') # load in data_dir and save_dir and CODE_DIR, depending on DEVICE value
assertthat::assert_that(!is.null(data_dir), msg='first arg must be: laptop, desktop, or ubergenno')

library(Seurat)
library(SeuratObject)


# =================== START ========================================================================
print(sprintf("[%s] START: Find Important Genes papalexi", Sys.time()))


# ============================================================================== #
# =================== Load Data from Jinhong's causarray tutorial ==============
# ============================================================================== #
print(sprintf("[%s]    - load subset of jin dataset that is causarray example data", Sys.time()))

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
}
# Y <- data.frame(t(as.matrix(counts)), check.names = FALSE) # cell-by-gene matrix
metadata <- sc.seurat@meta.data

# perturb <- metadata
# colnames(perturb) <- gsub("Perturbation", "trt_", colnames(perturb))
# perturb$trt_ <- relevel(as.factor(perturb$trt_), ref = "GFP")
# A <- data.frame(
#   model.matrix(~ trt_ - 1, data = perturb)[, -1, drop = FALSE],
#   check.names = FALSE
# ) # cell-by-trt matrix; remove the first (GFP control) column
# colnames(A) <- sub("^trt_", "", colnames(A))



# ============================================================================== #
# =================== Pre-processing ===========================================
# ============================================================================== #
print(sprintf("[%s]    - Pre-processing", Sys.time()))

# sometimes the names are bad for the glm formulas (e.g. AB-C). Remove non letter and number characters 
clean_name <- function(x) {  gsub(pattern='[^a-zA-Z0-9]', replacement='', x=x)  }


# =================== \__load Transcription Factor (TF) genes list ======================================
print(sprintf("[%s]        - loading list of TF genes", Sys.time()))
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



TF_names_cleaned = tolower(sapply(TF_names, clean_name))
genenames_df = data.frame(genename_jin     = row.names(counts), 
                          genename_cleaned = tolower(sapply(row.names(counts), clean_name))) |> 
               dplyr::mutate(is_transcription_factor = (genename_cleaned %in% TF_names_cleaned))

rm(tf_raw, tf, TF_names_cleaned)


# =================== \__calculate gene importance (deviance) ====================================================
print(sprintf("[%s]        - calculate gene importance (deviance)", Sys.time()))
# Input should be: row = feature, col = cell
gene_dev = scry::devianceFeatureSelection(object=counts, fam='binomial') # < 1 min
gene_dev_df = data.frame(gene_idx = 1:length(gene_dev), deviance=gene_dev)  |> 
              dplyr::arrange(desc(gene_dev)) |> dplyr::mutate(importance_rank = 1:n())  # add rank importance 
gene_dev_df = gene_dev_df[names(gene_dev), ] |> tibble::rownames_to_column(var = 'gene_name') |> # reorder to original
              dplyr::mutate(is_grna_target = (gene_name %in% unique(metadata$Perturbation))) # add T/F is a perturbation target


gene_dev_df = merge(gene_dev_df, genenames_df, 
                    by.x = 'gene_name', by.y = 'genename_jin', 
                    all.x = TRUE) |> dplyr::arrange(gene_idx)

write.csv(x = gene_dev_df, file = sprintf('%s/gene_dev.csv', save_dir), row.names = FALSE) # csv 

rm(gene_dev)




# =================== \__perform normalization ========================================================
print(sprintf("[%s]        - performing normalization", Sys.time()))
# genes in the rows and samples in the columns
gene_norm = scry::nullResiduals(counts, 
                                fam  = "binomial", 
                                type = "deviance") # 3221 genes x 2926 cells
# write.csv(   x = gene_norm, file = sprintf('%s/gene_norm.csv', save_dir)) # csv (dont, this is large)
saveRDS(object = gene_norm, file = sprintf('%s/gene_norm.rds', save_dir)) # this is smaller


# ============================================================================== #
# =================== Perform PCA + SPCA ======================================
# ============================================================================== #
print(sprintf("[%s]    - Perform PCA + SPCA", Sys.time()))

# =================== \__setup ========================================================
print(sprintf("[%s]        - prep (remove TF and grna target genes)", Sys.time()))
# select genes for PCA decomp, remove TF and grna targets which add signal
valid_gene_names = gene_dev_df |> 
  dplyr::filter((!is_transcription_factor) & !(is_grna_target)) |>    # no TF or grna target
  dplyr::arrange(gene_idx)


# PCA on all cells and genes? Or PCA on only the 'control' perturbations?

# PCA on all cells for now

# subsample to speed up computation time (e.g. not all 2926 cells)

dim(gene_norm)
# N_subsample = 5000
N_subsample = ncol(gene_norm) # don't subsample for jin dataset, there are 2926 cells here, maybe need to subsample later
set.seed(12345)
gene_norm_validgenes = gene_norm[valid_gene_names$gene_idx, ] # 3221 to 3103 genes

if(N_subsample == ncol(gene_norm)) {
  gene_norm_SPC = t(gene_norm_validgenes) # use all cells
} else {
  gene_norm_SPC = t(gene_norm_validgenes[, sample(1:ncol(gene_norm), N_subsample)]) # sample X out of XX cells
}
gene_norm_SPC = as.matrix(gene_norm_SPC)
# gene_norm_SPC = t(gene_norm[myGenenames_df |> dplyr::filter((!TF) & (!grna_target)) |> dplyr::pull(importance_rank),
#                             sample(1:ncol(gene_norm), N_subsample)]) # sample 5k out of 21k cells
dim(gene_norm_SPC)


# =================== \__PCA ========================================================
print(sprintf("[%s]        - PCA", Sys.time()))
# perform pca
pca_res = prcomp(x = gene_norm_SPC, rank.=100, retx = TRUE, center = TRUE, scale. = TRUE)

# construct NCs as the actual loadings
NC_loadings = t(gene_norm_validgenes) %*% pca_res$rotation

# save
dir.create(sprintf('%s/pca/', save_dir), showWarnings = FALSE)
saveRDS(NC_loadings, sprintf('%s/pca/NCloadings.rds', save_dir)) 


# =================== \__Sparse PCA ===========================================================
print(sprintf("[%s]        - Sparse PCA", Sys.time()))
dir.create(sprintf('%s/spca/', save_dir), showWarnings = FALSE)


# =================== \__\__Use SPC.cv to choose tuning parameters ======================================
print(sprintf("[%s]        - Sparse PCA: cv tuning params", Sys.time()))
# Use SPC.cv to choose tuning parameters:
if(F) {
  cv.out <- PMA::SPC.cv(gene_norm_SPC, 
                   sumabsvs = c(seq(1.2, 5, len = 5), seq(6, floor(sqrt(ncol(gene_norm_SPC))), len = 5)),
                   orth=TRUE)
  print(cv.out)
  plot(cv.out)
  print(cv.out$bestsumabsv)    # 42.75
  print(cv.out$bestsumabsv1se) # 30.5  =smallest sumabsv value that has CV error within 1 SE of best CV error:  30.5 
 
  saveRDS(cv.out, 
          sprintf('%s/spca/cvout.rds', save_dir)) # save cv res
}

# =================== \__\__perform SPCA ======================================
print(sprintf("[%s]        - Sparse PCA: fit", Sys.time()))
my_sumabsv = 30.5
my_K       = 100   # number of factors in the PMD to be returned
out.orth <- PMA::SPC(gene_norm_SPC,
                sumabsv=my_sumabsv, # tuning parameter
                # sumabsv=cv.out$bestsumabsv1se, # 33.5 not strong enough. v's not sparse and number of nonzero coefs too large
                K=my_K, 
                orth=TRUE)
saveRDS(out.orth, 
        sprintf('%s/spca/outorth_sumabs=%.1f_K=%d_N=%d.rds', save_dir, my_sumabsv, my_K, N_subsample)) # save res

# gene_norm_noTFTargets = gene_norm[myGenenames_df |> dplyr::filter((!TF) & (!grna_target)) |> dplyr::pull(importance_rank), ]
# construct NCs as the actual loadings
NC_loadings = t(gene_norm_validgenes) %*% out.orth$v


saveRDS(NC_loadings, 
        sprintf('%s/spca/NCloadings_sumabs=%.1f_K=%d_N=%d.rds', save_dir, my_sumabsv, my_K, N_subsample)) # save 


# gene_norm_noTFTargets |> dim()
# out.orth$v|> dim()
# NC_loadings |> dim()
# 
# # check orthogonality (not exactly Identity matrix bc SPC used sample)
# t(NC_loadings) %*% NC_loadings / nrow(NC_loadings)
# 
# 
# 
# heatmap( t(NC_loadings) %*% NC_loadings / nrow(NC_loadings), 
#          Rowv=NA, Colv=NA, col = heat.colors(256),  margins=c(5,10))

# =================== END ==========================================================================
print(sprintf("[%s] END", Sys.time()))






# ============================================================================== #
# =================== chromosome information ====================================
# ============================================================================== #
if(F) {
  # chr information doesn't actually matter scientifically, but we can just do it to be safe
  
  library(biomaRt)
  
  # make dir for saving
  dir.create(sprintf('%s/chromosome', save_dir), showWarnings = FALSE)
  
  # attributes to pull
  myAttributes = c("wikigene_name", "wikigene_id", "chromosome_name", 
                   "ensembl_gene_id", "external_gene_name")
  myGrnaGenenames = c(row.names(counts), unique(metadata$Perturbation)) |> unique() |> sort()
  myGrnaGenenames[100:105]
  myGrnaGenenames[500:510]
  
  # make a mart object
  ensembl = useMart("ensembl") #listDatasets(ensembl)
  has_mouse <- function(x) {grepl(pattern='mouse', x=x)}
  
  listDatasets(ensembl)[sapply(X = listDatasets(ensembl)$description, FUN = has_mouse), ]
  
  
  ensembl = useDataset("hsapiens_gene_ensembl", mart=ensembl)
  
  
  # get chromosome info ======================================
  print(sprintf("[%s]    - get GENE chromosome info from biomaRt", Sys.time()))
  
  g = getBM(attributes = myAttributes,
            filters = "wikigene_name",
            values = myGrnaGenenames[1:10],
            mart = ensembl)
  
  
  # format GENE chr df and save ===================================
  print(sprintf("[%s]        - format GENE chr df and save", Sys.time()))
  
  # some chr names are like CHR_HSCHR6_MHC_SSTO_CTG1??
  # g |> filter(!chromosome_name %in% c(as.character(1:100), 'A', 'Y'))
  
  chr_df = merge(data.frame('wikigene_name' = myGenenames,
                            importance_rank = 1:length(myGenenames),
                            gene_idx        = myGenenames_idx),
                 g |> filter(chromosome_name %in% c(as.character(1:100), 'A', 'Y')),
                 all.x = TRUE) |> arrange(importance_rank)
  
  # some chr names have multiple ensembl names, just choose smallest ensembl_gene_id
  chr_df = chr_df |> group_by(wikigene_name) |> slice_min(ensembl_gene_id)
  
  write.csv(chr_df, file = sprintf('%s/chromosome/gene_chromosome.csv', save_dir), row.names = F)
}








