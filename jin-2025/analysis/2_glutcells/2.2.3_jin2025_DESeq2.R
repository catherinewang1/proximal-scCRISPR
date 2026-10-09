# ------------------------------------------------------------------------------------------------ #
#     Initial Estimates on Jin-2025 dataset on glut cell type
#     Mouse brain data from Xin Jin's lab ~ 2025 and downloaded from dropbox 
# Email: Ergan to Jin-Hong to Catherine:
#     Please check the following two links of the data available. 
#     Those data is shared in 2 batches by Seoyeon, where one is containing 10 perturbations and another one contains 48 perturbations.
#     
# 
# 2.1 preprocessing
# - Pre-processing (gene importance and normalization)
# - Perform PCA 
# 2.2 estimates
# - 2.2.1 glms + proximal

# - 2.2.3 causarray
# - 2.2.4 deseq2 (run this after causarray to used the saved data)
# - 2.2.5 sceptre
# ------------------------------------------------------------------------------------------------ #
args = commandArgs(trailingOnly = TRUE)
# args = c('macbook')


assertthat::assert_that(length(args) > 0, msg="must give arg for specifying device eg 'Rscript <filename>.R ubergenno'")
DEVICE = args[1]
source('../PATHS.R') # load in data_dir and save_dir and CODE_DIR, depending on DEVICE value
save_dir = paste0(save_dir, "/2_glutcells")
assertthat::assert_that(!is.null(data_dir), msg='first arg must be: laptop, desktop, or ubergenno')

library(dplyr)
library(DESeq2)
library(BiocParallel)
register(MulticoreParam(4))



# =================== START ========================================================================
print(sprintf("[%s] START: deseq2 for Jin2025", Sys.time()))
dir.create(sprintf('%s/deseq2/', save_dir))


# h5ad save file from causarray analysis to be reused so that we have the fake NT treatments
assertthat::assert_that('jin2025_glut.h5ad' %in% list.files(save_dir), msg = 'prev saved cleaned and subsetted data should be saved beforehand (subset glut cells + create fake NT treatments + subsample control NT cells')

# ============================================================================== #
# ================================ Load ========================================
# ============================================================================== #
print(sprintf("[%s]    - load h5ad (prev saved from causarray analysis)", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))


h5adfile   = sprintf('%s/jin2025_glut.h5ad', save_dir)
jin2025_glut = anndataR::read_h5ad(h5adfile)

# --------------------- Load gene_counts
gene_counts = jin2025_glut$X |> t()

# --------------------- Load grna
grna_A = jin2025_glut$obs$gene_target_w_fakeNT |> factor() |> relevel('Non_target')
grna_names = setdiff(levels(grna_A), 'Non_target') # non Non_target (control) grna names


all(jin2025_glut$obs |> row.names() == jin2025_glut$X |> row.names()) # check cell names match
all(jin2025_glut$obs |> row.names() == gene_counts |> colnames()) # check cell names match
rm(jin2025_glut)


# ============================================================================== #
# ================================ DESeq2 ======================================
# ============================================================================== #
dds <- DESeqDataSetFromMatrix(countData = gene_counts,
                              colData = data.frame(grna = grna_A),
                              design = ~grna)


# --------------------- DESeq2 run  
# modify from JinHong's causarray https://github.com/JinHongDu-Lab/causarray/blob/2f6535c406079db72457edb8afbf03eb6ae363b7/paper/methods/R_functions.R#L126
# dds <- DESeq(dds)
# to avoid NA p-values
geoMeans <- exp(MatrixGenerics::rowMeans(log(counts(dds))))
if(all(geoMeans==0)){
  geoMeans <- apply(counts(dds), 1, function(row) if (all(row == 0)) 0 else exp(mean(log(row[row != 0]))))
}
sizeFactors(dds) <- estimateSizeFactorsForMatrix(counts(dds), geoMeans=geoMeans)
# dds <- DESeq(dds, minReplicatesForReplace=Inf)
dds <- DESeq(dds, test = 'Wald', minReplicatesForReplace=Inf, parallel=TRUE, BPPARAM=MulticoreParam(4))
# saveRDS(dds, sprintf('%s/deseq2/deseq2_obj.rds', save_dir)) # no this is large



# --------------------- DESeq gather results  

# The & symbols were replaced with . 


ddsResNames = setdiff(resultsNames(dds),"Intercept")
deseq_res = NULL
for(resName in ddsResNames) {
  res = results(dds, name = resName, cooksCutoff=FALSE, independentFiltering=FALSE) |> data.frame()
  res = res |> tibble::rownames_to_column('gene') |> mutate(ddsResName = resName, .before='gene')
  
  deseq_res = rbind(deseq_res, res)
}

#' clean ddsResName so that they match the original grna names
#' remove the extra grna_ and remove the vs_Non_taarget
#' @example sapply(X=deseq_res$ddsResName[c(1, 10^4, 10^5)], FUN=clean_ddsResName)
clean_ddsResName <- function(x) {
  x = gsub(pattern='grna_', replacement='', x=x)
  x = gsub(pattern='_vs_Non_target', replacement='', x=x)
  x = gsub(pattern='\\.', replacement='&', x=x)
  return(x)
}

deseq_res$grna = sapply(X=deseq_res$ddsResName, FUN=clean_ddsResName)
deseq_res = deseq_res |> relocate(grna, .before ='gene')


# check all grna names are in original grna names
all(unique(deseq_res$grna) %in% grna_names)


# less save if names don't match
# # change the grna_names the loop loops over
# grna_names2 = sapply(resultsNames(dds), FUN = function(x){strsplit(x=x, split = '_vs_')[[1]][1]}) |> unname()
# 
# deseq_res = NULL
# for(grna_name_ in grna_names2) {
#   cur_name = paste0('grna_', grna_name_, '_vs_Non_target')
#   res = results(dds, name = cur_name, cooksCutoff=FALSE, independentFiltering=FALSE) |> data.frame()
#   res = res |> tibble::rownames_to_column('gene') |> mutate(grna = grna_name_, .before='gene')
#   
#   deseq_res = rbind(deseq_res, res)
# }


# --------------------- DESeq write results  
write.csv(deseq_res, sprintf('%s/deseq2/deseq_effects.csv', save_dir), row.names = FALSE)








# 
# # test
# df = read.csv(sprintf('%s/deseq2/deseq_effects.csv', save_dir))












# ============================================================================== #
# ========= Older Version: choosing genes from all normalized genes ============
# =========                and does not perform Null tests          ============
# ============================================================================== #
# 
# # =================== START ========================================================================
# print(sprintf("[%s] START: deseq2 for Jin2025", Sys.time()))
# dir.create(sprintf('%s/deseq2/', save_dir))
# 
# 
# # ============================================================================== #
# # ================================ Load ========================================
# # ============================================================================== #
# print(sprintf("[%s]    - load h5 (prev saved from processData)", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
# h5file   = sprintf('%s/jin2025_glut.h5', save_dir)
# 
# # --------------------- Load gene_counts
# gene_counts = HDF5Array::HDF5Array(h5file, name='gene_counts/gene_counts')
# row.names(gene_counts) = rhdf5::h5read(h5file, name = 'gene_counts/.gene_counts_dimnames/1') |> as.vector()
# colnames( gene_counts) = rhdf5::h5read(h5file, name = 'gene_counts/.gene_counts_dimnames/2') |> as.vector()
# 
# # --------------------- Load grna
# grna     = h5lite::h5_read(h5file, name = 'grna/grna')
# 
# # --------------------- Load gene_dev_df
# gene_dev_df = read.csv(file = sprintf('%s/gene_dev.csv', save_dir)) 
# 
# 
# 
# 
# 
# # ============================================================================== #
# # ================================ Choose Y (responses) ========================
# # ============================================================================== #
# 
# # --------------------- choose Y (responses)
# genes_nonzero = row.names(gene_counts)[(rowSums(gene_counts == 0) >= GENE_NUMNONZEROCELLS_MIN)] # genes that meet the nonzero threshold
# all_Ys_df = gene_dev_df |>  
#   filter(gene_name %in% genes_nonzero) |>  
#   filter(is_normalized_gene) |> 
#   # filter(!is_grna_target & !is_transcription_factor) |>
#   arrange(importance_rank)
# 
# all_Ys = all_Ys_df |> filter(importance_rank <= GENE_IMPORTANCERANK_MAX | is_grna_target) |> pull(gene_name)
# 
# 
# 
# # ============================================================================== #
# # ================================ Subsample control cells =====================
# # ============================================================================== #
# # --------------------- subsample control cells  e.g. subsample 10,000 control cells out of ~135k
# set.seed(123456789)
# control_cell_names = which(grna['Non_target', ] == 1) |> names()
# if(!is.na(NUM_SUBSAMPLE_CONTROL_CELLS) & (NUM_SUBSAMPLE_CONTROL_CELLS <= length(control_cell_names))) {
#   control_cell_names = sample(x=control_cell_names, size=NUM_SUBSAMPLE_CONTROL_CELLS, replace=FALSE)
# }
# trt_cell_names = which(grna['Non_target', ] == 0) |> names()
# selected_cell_names = c(trt_cell_names, control_cell_names)
# 
# 
# gene_counts = gene_counts[all_Ys, selected_cell_names]
# grna        =        grna[      , selected_cell_names]
# 
# 
# rm(control_cell_names, trt_cell_names)
# 
# 
# 
# # ============================================================================== #
# # ================================ DESeq2 ======================================
# # ============================================================================== #
# # --------------------- DESeq2 setup
# grna_A = row.names(grna)[apply(X = (grna==1), FUN = which, MARGIN = 2)] |> factor() |> relevel('Non_target')
# 
# grna_names = setdiff(row.names(grna), 'Non_target') # non GFP grna names
# 
# dds <- DESeqDataSetFromMatrix(countData = gene_counts,
#                               colData = data.frame(grna = grna_A),
#                               design = ~grna)
# 
# 
# # --------------------- DESeq2 run  
# # modify from JinHong's causarray https://github.com/JinHongDu-Lab/causarray/blob/2f6535c406079db72457edb8afbf03eb6ae363b7/paper/methods/R_functions.R#L126
# # dds <- DESeq(dds)
# # to avoid NA p-values
# geoMeans <- exp(MatrixGenerics::rowMeans(log(counts(dds))))
# if(all(geoMeans==0)){
#   geoMeans <- apply(counts(dds), 1, function(row) if (all(row == 0)) 0 else exp(mean(log(row[row != 0]))))
# }
# sizeFactors(dds) <- estimateSizeFactorsForMatrix(counts(dds), geoMeans=geoMeans)
# # dds <- DESeq(dds, minReplicatesForReplace=Inf)
# dds <- DESeq(dds, test = 'Wald', minReplicatesForReplace=Inf, parallel=TRUE, BPPARAM=MulticoreParam(4))
# # saveRDS(dds, sprintf('%s/deseq2/deseq2_obj.rds', save_dir)) # no this is large
# 
# # --------------------- DESeq gather results  
# deseq_res = NULL
# for(grna_name_ in grna_names) {
#   cur_name = paste0('grna_', grna_name_, '_vs_Non_target')
#   res = results(dds, name = cur_name, cooksCutoff=FALSE, independentFiltering=FALSE) |> data.frame()
#   res = res |> tibble::rownames_to_column('gene') |> mutate(grna = grna_name_, .before='gene')
#   
#   deseq_res = rbind(deseq_res, res)
# }
# 
# 
# # --------------------- DESeq write results  
# write.csv(deseq_res, sprintf('%s/deseq2/deseq_effects.csv', save_dir), row.names = FALSE)
# 
# 
# # # test
# # df = read.csv(sprintf('%s/deseq2/deseq_effects.csv', save_dir))
# # 
# # df
# # dim(df)
# # 
# # nrow(df) / 60
# # df$gene |> unique() |> length()
# 
# 
