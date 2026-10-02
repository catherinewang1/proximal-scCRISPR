# ------------------------------------------------------------------------------------------------ #
#     Estimate ATE/Log FC
# For only 1/few choices of Proximal....
# - Positive Tests: perform all
# - Discovery Tests: 
#    + Perform initial Linear Reg for ALL/A LOT OF tests
#    + For those with significance (small p-values), perform proximal
# ------------------------------------------------------------------------------------------------ #
args = commandArgs(trailingOnly = TRUE)
# args = c('macbook')


assertthat::assert_that(length(args) > 0, msg="must give arg for specifying device eg 'Rscript <filename>.R ubergenno'")
DEVICE = args[1]
source('../PATHS.R') # load in data_dir and save_dir and CODE_DIR, depending on DEVICE value
save_dir = paste0(save_dir, '_exneu_unpert') # FOR EXCITATORY NEURONS + USE UNPERTURBED CELLS AS CONTROLS
assertthat::assert_that(!is.null(data_dir), msg='first arg must be: laptop, desktop, or ubergenno')

library(dplyr)
# library(future.apply)
# options(future.globals.maxSize= 850*1024^2) #1st num is MB
# options(future.globals.maxSize= 2500*1024^2) #1st num is MB
# plan(multisession, workers = 24)
# plan(multisession, workers = 15)
# plan(sequential)

library(DESeq2)



dir.create(sprintf('%s/deseq2/', save_dir))




# --------------------- load data (old)
# gene = readRDS(sprintf('%s/gene_counts.rds', save_dir))
# 
# grna = read.csv(sprintf('%s/grna.csv', save_dir), row.name = 'X') # for cells w a pert assignment 22085 x 39
# grna = t(grna) #  30 perturbations x 4744 cells
# row.names(grna) = stringr::str_to_title(row.names(grna)) # hopefully these match broad Perturbation names now



# --------------------- load data 
h5file <- sprintf('%s/jin_exneu.h5', save_dir)

grna      = h5lite::h5_read(h5file, "grna/grna")      # 31 8215
gene      = h5lite::h5_read(h5file, "counts/counts")  # 4017 8215


# --------------------- DESeq2 setup
grna_A = row.names(grna)[apply(X = (grna==1), FUN = which, MARGIN = 2)] |> factor() |> relevel('GFP') |> relevel('UnPert')
grna_names = setdiff(row.names(grna), 'UnPert') # non GFP grna names

dds <- DESeqDataSetFromMatrix(countData = gene,
                              colData = data.frame(grna = grna_A),
                              design = ~grna)


# --------------------- DESeq run  
# modify from JinHong's causarray https://github.com/JinHongDu-Lab/causarray/blob/2f6535c406079db72457edb8afbf03eb6ae363b7/paper/methods/R_functions.R#L126
# dds <- DESeq(dds)
# to avoid NA p-values
geoMeans <- exp(MatrixGenerics::rowMeans(log(counts(dds))))
if(all(geoMeans==0)){
  geoMeans <- apply(counts(dds), 1, function(row) if (all(row == 0)) 0 else exp(mean(log(row[row != 0]))))
}
sizeFactors(dds) <- estimateSizeFactorsForMatrix(counts(dds), geoMeans=geoMeans)
dds <- DESeq(dds, minReplicatesForReplace=Inf)


deseq_res = NULL
for(grna_name_ in grna_names) {
  cur_name = paste0('grna_', grna_name_, '_vs_UnPert')
  res = results(dds, name = cur_name, cooksCutoff=FALSE, independentFiltering=FALSE) |> data.frame()
  res = res |> tibble::rownames_to_column('gene') |> mutate(grna = grna_name_, .after='gene')
  
  deseq_res = rbind(deseq_res, res)
}
# --------------------- 


write.csv(deseq_res, sprintf('%s/deseq2/deseq_effects.csv', save_dir), row.names = FALSE)








# ============================================================================== #
# =================== TRASH ==========================================
# ============================================================================== #


# process to #features x #cells to be consistent
# 
# grna: perturbation-by-cell 0/1 trt assignment
# gene:      gene-by-cell count gene expression 
# gene_norm: gene-by-cell normalized gene expression
# 
# but others may be #cells x #features
# metadata: #cells x 13 variables




# # =================== \__Load metadata: carried over from Broad ===========================
# # metadata = read.csv(sprintf('%s/processGEO/grna.csv', data_dir)) # for all
# metadata0 = read.csv(sprintf('%s/cell_metadata.csv', save_dir), row.name = 'X') # excitatory neuron cells w/ a pert assignment w/ > 80(?) sample size
# gene_dev = read.csv(sprintf('%s/gene_dev.csv', save_dir))
# # Add TF info, but does not have for GFP (will be NA cols for is_transcription_factor and importance_rank)
# # merge(metadata, gene_dev |> select(gene_name, is_transcription_factor, importance_rank), 
# #       by.x = 'Perturbation_alias', by.y = 'gene_name') |> dim()
# # all(metadata$Perturbation %in% gene_dev$gene_name)
# # all(metadata$Perturbation_alias %in% gene_dev$gene_name)
# # unique(metadata$Perturbation_alias)[which(!unique(metadata$Perturbation_alias) %in% gene_dev$gene_name)]
# # metadata |> filter(Perturbation_alias %in% c('Mll1', 'Myst4'))
# # metadata |> filter(Perturbation_alias %in% c('Kmt2a'))
# # 'Kmt2a' %in% gene_dev$gene_name
# metadata = merge(metadata0, gene_dev |> select(gene_name, is_transcription_factor, importance_rank), 
#                  by.x = 'Perturbation_alias', by.y = 'gene_name', all.x = TRUE)
# row.names(metadata) = metadata$NAME
# metadata = metadata[row.names(metadata0), ] # reorder to same as original
# 
# 
# # dim(metadata); dim(metadata0)
# # row.names(metadata)[1:5];row.names(metadata0)[1:5]
# # all(row.names(metadata) %in% row.names(metadata0))
# # all(row.names(metadata0) %in% row.names(metadata))
# # all(row.names(metadata0) == row.names(metadata))
# rm(metadata0, gene_dev)
# 
# # =================== \__Load grna: carried over from Broad ===========================
# grna = read.csv(sprintf('%s/grna.csv', save_dir), row.name = 'X') # for cells w a pert assignment 22085 x 39
# grna = t(grna) #  30 perturbations x 4744 cells
# row.names(grna) = stringr::str_to_title(row.names(grna)) # hopefully these match broad Perturbation names now
# 
# # =================== \__Load counts: processed from GEO ===========================
# # load in top 4000 genes (count and continuous)
# gene = readRDS(sprintf('%s/gene_counts.rds', save_dir))
# gene = as.matrix(gene)
# 
# # =================== \__Load continuous: normalized counts ==============================================
# gene_norm = readRDS(sprintf('%s/gene_norm.rds', save_dir))
# gene_norm = as.matrix(gene_norm)
# 
# (dim(gene_norm) == dim(gene)) |> all() # same dimensions (all genes)
# (row.names(gene_norm) == row.names(gene)) |> all() # same genes and order
# (row.names(gene_norm) == row.names(gene)) |> table() # same genes and order
# (colnames(gene_norm) == colnames(gene)) |> all()  # same cells and order
# 
# # =================== \__Reorder cells in data to match gene/gene_norm order ===========================
# grna      = grna[, colnames(gene)]
# metadata  = metadata[colnames(gene), ]
# 
# 
# (colnames(gene_norm) == colnames(grna)) |> table() # same cells
# (colnames(gene_norm) == row.names(metadata)) |> table() # same cells
# 
# format(object.size(gene), units = 'Gb') # less when loaded as original sparse matrix form. but then transformed to matrix
# format(object.size(gene_norm), units = 'Gb')
# format(object.size(grna), units = 'Gb')
# format(object.size(metadata), units = 'Gb')
# 
# 
# 
# # =================== \__Load gene info (dev, TF, ...) ==============================================
# 
# gene_dev_df = read.csv(file = sprintf('%s/gene_dev.csv', save_dir)) # load gene dev info
# gene_dev_df = gene_dev_df |> filter(gene_name %in% row.names(gene)) # only choose top 4000 and grna target
# gene_dev_df$gene_name_upper = toupper(gene_dev_df$gene_name)
# # gene_dev_df |> filter(gene_name_upper == 'PISD') # there is Pisd and PISD apparently??
# row.names(gene_dev_df) = gene_dev_df$gene_name
# gene_dev_df = gene_dev_df[row.names(gene), ] # select the rows=genes that are top or grna target (the selected gene and gene_norm)
# 
# 
# (row.names(gene_norm) %in% gene_dev_df$gene_name) |> all() 
# (row.names(gene) %in% gene_dev_df$gene_name) |> all() 
# 
# (row.names(gene_norm) == gene_dev_df$gene_name) |> table()
# (row.names(gene) == gene_dev_df$gene_name) |> table()
# 
# gene_dev_df |> head()
# # gene_name gene_idx deviance importance_rank is_grna_target genename_cleaned is_transcription_factor gene_name_upper
# # Actb        Actb     9932 579441.3               1          FALSE             actb                   FALSE            ACTB
# # Nrgn        Nrgn    17935 242049.9               2          FALSE             nrgn                   FALSE            NRGN
# # Malat1    Malat1    27229 232525.0               3          FALSE           malat1                   FALSE          MALAT1
# # Tuba1a    Tuba1a    23642 187071.7               4          FALSE           tuba1a                   FALSE          TUBA1A
# # Pcp4        Pcp4    24555 179699.0               5          FALSE             pcp4                   FALSE            PCP4
# # Calm2      Calm2    25756 136040.4               6          FALSE            calm2                   FALSE           CALM2
# 
# 
# 
# 
# 
# # ======
# 
# 
# library()
# 
# 
# library(Seurat)
# library(SeuratObject)
# 
# a = CreateSeuratObject(gene, assay = 'counts', meta.data = data.frame( t(grna)))
# 
# gfp_idx = which(grna['Gfp', ]==1) |> names()
# trt_idx = which(grna['Satb2', ]==1) |> names()
# 
# Seurat::FindMarkers(a, slot = 'counts', cells.1 = gfp_idx, cells.2 = trt_idx, 
#                     test.use = 'DESeq2')
# 
# BiocManager::install("DESeq2")
# library("DESeq2")
# 
# 
# grna_A = row.names(grna)[apply(X = (grna==1), FUN = which, MARGIN = 2)] |> factor()
# grna_A = relevel(grna_A, 'Gfp')
# 
# grna_df = data.frame(t(grna)) |> select(-Gfp)
# 
# deseq_formula = as.formula(paste0('~ ', paste0(colnames(grna_df), collapse = '+')))
# 
# dds <- DESeqDataSetFromMatrix(countData = gene[1:10, ],
#                               # colData = data.frame(grna = grna_A),
#                               colData = grna_df,
#                               # design = ~ grna
#                               design = deseq_formula
#                               )
# dds <- DESeq(dds)
# res <- results(dds)
# res
# 
# grna_names = setdiff(row.names(grna), 'Gfp') # non GFP grna names
# 
# t0 <- Sys.time()
# dds <- DESeqDataSetFromMatrix(countData = gene[1:10, ],
#                               colData = data.frame(grna = grna_A),
#                               # colData = grna_df,
#                               design = ~grna)
# # dds <- DESeq(dds)
# # to avoid NA p-values
# geoMeans <- exp(MatrixGenerics::rowMeans(log(counts(dds))))
# if(all(geoMeans==0)){
#   geoMeans <- apply(counts(dds), 1, function(row) if (all(row == 0)) 0 else exp(mean(log(row[row != 0]))))
# }
# sizeFactors(dds) <- estimateSizeFactorsForMatrix(counts(dds), geoMeans=geoMeans)
# dds <- DESeq(dds, minReplicatesForReplace=Inf)
# # equivalently
# 
# dds
# res <- results(dds, name="trt_1_vs_0", cooksCutoff=FALSE, independentFiltering=FALSE)
# 
# 
# i = 4
# 
# res <- results(dds, name = resultsNames(dds)[i], cooksCutoff=FALSE, independentFiltering=FALSE)
# res
# 
# 
# resultsNames(dds)




