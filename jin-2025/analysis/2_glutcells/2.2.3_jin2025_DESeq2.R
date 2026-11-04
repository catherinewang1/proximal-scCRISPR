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
assertthat::assert_that(!is.null(data_dir), msg='first arg must be: laptop, desktop, or ubergenno')

library(dplyr)
library(DESeq2)
library(BiocParallel)
register(MulticoreParam(4))




PARAMS = list( 
  'TAG'       = 'L45Glut',
  'subfolder' = '3_L45glutcells',
  'seed'      = 12345,
  'A_MIN_TRTCELLS'  = 40, # filter for grnas with at least this number of cells (low prev, can set higher here))
  'Y_MAX_GENE_IMPORTANCE' = 1500
)

save_dir = paste0(save_dir, '/', PARAMS[['subfolder']])
set.seed(PARAMS[['seed']])

h5adfilename = paste0('jin2025_', PARAMS['TAG'], '.h5ad')

# =================== START ========================================================================
print(sprintf("[%s] START: deseq2 for Jin2025 for L45 Glut Cells", Sys.time()))
dir.create(sprintf('%s/deseq2/', save_dir))


# h5ad save file from causarray analysis to be reused so that we have the fake NT treatments
assertthat::assert_that(h5adfilename %in% list.files(save_dir), msg = 'prev saved cleaned and subsetted data should be saved beforehand (subset glut cells + create fake NT treatments + subsample control NT cells')

# ============================================================================== #
# ================================ Load ========================================
# ============================================================================== #
print(sprintf("[%s]    - load h5ad (prev saved from causarray analysis)", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))


jinn2025 = anndataR::read_h5ad(sprintf('%s/%s', save_dir, h5adfilename))


# --------------------- Load grna
grna_A_orig = jinn2025$obs$gene_target_w_fakeNT |> factor() |> relevel('Non_target')

# --------------------- Load gene_dev_df
gene_dev_df = read.csv(sprintf('%s/gene_dev.csv', save_dir))

# --------------------- Subset grna and gene 
# grnas with enough cells 
grna_counts = table(grna_A_orig)
A_treatments = names(grna_counts)[grna_counts>=PARAMS[['A_MIN_TRTCELLS']]]
subset_cells = jinn2025$obs |> filter(gene_target_w_fakeNT %in% A_treatments) |> row.names()

# genes that are important/high deviance
Y_responses = gene_dev_df |> filter((gene_name %in% A_treatments) | importance_rank <= PARAMS[['Y_MAX_GENE_IMPORTANCE']]) |> arrange(gene_idx) |> pull(gene_name)


# --------------------- construct grna
grna_A = jinn2025$obs[subset_cells, 'gene_target_w_fakeNT'] |> factor() |> relevel('Non_target')

# --------------------- Load gene_counts
gene_counts = jinn2025$X[subset_cells, Y_responses] |> t()
grna_names = setdiff(levels(grna_A), 'Non_target') # non Non_target (control) grna names



all(jinn2025$obs |> row.names() == jinn2025$X |> row.names()) # check cell names match
all(jinn2025$obs |> row.names() == gene_counts |> colnames()) # check cell names match
rm(jinn2025, grna_counts, A_treatments, Y_responses, subset_cells)


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





# --------------------- DESeq write results  
write.csv(deseq_res, sprintf('%s/deseq2/deseq_effects.csv', save_dir), row.names = FALSE)















