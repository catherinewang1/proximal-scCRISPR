# ------------------------------------------------------------------------------------------------ #
#     Analysis on jin-2020 dataset- excitatory cells chosen from Broad website
#                                   and counts from GEO
# 
# - Pre-processing (gene importance and normalization)
# - Perform PCA + SPCA 
# ------------------------------------------------------------------------------------------------ #
args = commandArgs(trailingOnly = TRUE)
# args = c('macbook')


assertthat::assert_that(length(args) > 0, msg="must give arg for specifying device eg 'Rscript <filename>.R ubergenno'")
DEVICE = args[1]
source('../PATHS.R') # load in data_dir and save_dir and CODE_DIR, depending on DEVICE value
save_dir = paste0(save_dir, '_exneu_unpert')
assertthat::assert_that(!is.null(data_dir), msg='first arg must be: laptop, desktop, or ubergenno')



library(dplyr)
library(Seurat)
library(SeuratObject)

NUM_IMPORTANT_GENES = 4000


dir.create(save_dir)

# =================== START ========================================================================
print(sprintf("[%s] START: Find Important Genes on Jin Full Dataset from GEO, subselect excitatory neurons", Sys.time()))

# ============================================================================== #
# =================== Load Data from combineBroadGEO download    ==============
# ============================================================================== #
print(sprintf("[%s]    - load combineBroadGEO", Sys.time()))

# From Broad
metadata = read.delim(sprintf('%s/broad/metadata/meta_PertCortex.txt', data_dir), header = T, stringsAsFactors = F)
metadata = metadata[-1, ] # 1st row (after header) is col type
rownames(metadata) <- metadata[, 1]

metadata = metadata |> mutate(isKey = (isKey == 'TRUE'), # string to boolean
                              isAnalysed = (isAnalysed == 'TRUE'),# string to boolean
                              SCRUBLET = as.numeric(SCRUBLET)) # string to numeric

# From combineBroadGEO
counts = readRDS(sprintf('%s/combineBroadGEO/gene.rds', data_dir)) #  27998 genes x  49067 cells


# =================== \__subset cells ==================================================
# subset cells for:
#     Excitatory Neuron cell type 
#   & has a perturbatin assignment (include unperturbed cells)
#   & is a Key??
#   & assigned a perturbation with at least sample size 50

# filter for Excitatory neurons, and has perturbation, and should remove if sample 17/or is it batch 17?
metadata_excitatory = metadata |> filter(CellType == 'Excitatory') #  & !is.na(Perturbation)) 

metadata_excitatory |> select(isKey, isAnalysed) |> table()
metadata_excitatory$Perturbation |> table() |> sort() # what is 'NotKey'
metadata_excitatory |> select(isKey, Perturbation) |> table() # removing bc these don't seem to be a 'valid' Pert

# filter for isKey
metadata_excitatory = metadata_excitatory |> filter(isKey)  # 8450   10
dim(metadata_excitatory)


# perturbations with at least 50 sample size.
perturbations_with_enoughcells = metadata_excitatory |> group_by(Perturbation) |> summarize(count = n()) |> filter(count >= 50) |> pull(Perturbation)


# subset cells with only these perturbations or unperturbed
metadata_excitatory = metadata_excitatory |> filter(is.na(Perturbation) | Perturbation %in% perturbations_with_enoughcells) # 8215 cells x  10
metadata_excitatory$Perturbation |> table(useNA = "ifany")

# replace unperturbed cell's label from NA to 'UnPert'
metadata_excitatory[is.na(metadata_excitatory$Perturbation), 'Perturbation'] = 'UnPert' 
# sort cell names a bit- Batch, then cbc
metadata_excitatory = metadata_excitatory |> arrange(Batch, NAME)


# plot(x = metadata_excitatory$isAnalysed, y = metadata_excitatory$SCRUBLET, ylim = c(0, .07))
# metadata_excitatory |> group_by(isAnalysed) |> summarize(scrublet_mean = mean(SCRUBLET))


# =================== \__subset cells in counts ==================================================
counts_excitatory = counts[, metadata_excitatory$NAME] # 27998 genes x  8215 cells


# create Perturbation_alias which is the perturbation name if the gene is measured or not present in counts (GEO h5 files)
# but is changed to another alias if the perturbation's alias is present in counts
# e.g. Perturbation SCN2A is not measured in counts, but the alias SCN2A1 is
# e.g. Perturbation MLL1 is not measured in counts, but alias KMT2A is
# e.g. Perturbation MYST4 is not measured in counts, but alias KAT6B is
# e.g. Perturbation GFP is not measured in counts, but is also not in counts and cannot find aliases 
# now, only missing 1 grna target measurement
# Perturbation_alias_dict = list('SCN2A'='SCN2A1',  # capital if using GEO + supplemental file 5
#                                'MLL1'='KMT2A',    # for processed straight from GEO: it uses all 3, 
#                                'MYST4'='KAT6B')   # for using Broad Metadata, use MLL1 and MYST4... they mixed these...
Perturbation_alias_dict = list('Scn2a'='Scn2a1',  # title case if using GEO + Broad
                               'Mll1'='Kmt2a',    # for processed straight from GEO: it uses all 3, 
                               'Myst4'='Kat6b')   # for using Broad Metadata, use MLL1 and MYST4... they mixed these...


perturbation_target_names = unique(metadata$Perturbation) |> sort() # perturbation names and gene names don't always match??
perturbation_target_names = perturbation_target_names[!is.na(perturbation_target_names)] 
perturbation_target_alias = sapply(perturbation_target_names, FUN=function(x){if(x %in% names(Perturbation_alias_dict)){return(Perturbation_alias_dict[[x]])} else {x}})

perturbation_alias_df = data.frame('Perturbation_alias'=perturbation_target_alias) |> tibble::rownames_to_column(var = 'Perturbation') |> mutate(Perturbation_alias_lower = tolower(Perturbation_alias))
write.csv(perturbation_alias_df, sprintf('%s/perturbation_alias.csv', save_dir), row.names = F) 



# =================== \__save gene and grna subset ==================================================

cell_name_order = metadata_excitatory$NAME # preserve original order
metadata_excitatory = merge(metadata_excitatory, perturbation_alias_df, all.x = TRUE) # add pert target's alias
row.names(metadata_excitatory) = metadata_excitatory$NAME
metadata_excitatory = metadata_excitatory[cell_name_order, ]
rm(cell_name_order)
# write.csv(metadata_excitatory, sprintf('%s/cell_metadata.csv', save_dir), row.names = TRUE)
# saveRDS(    counts_excitatory, sprintf('%s/gene.rds', save_dir))
# write.csv(  counts_excitatory, sprintf('%s/gene.csv', save_dir), row.names = TRUE)





dim(counts_excitatory) # 27998  8215
format(object.size(counts_excitatory), units = 'Gb') # .4 Gb


rm(counts, metadata) # original counts and metadata with all cells
rm(Perturbation_alias_dict) # should be in metadata/grna now

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
genenames_df = data.frame(genename_jin     = row.names(counts_excitatory), 
                          genename_cleaned = tolower(sapply(row.names(counts_excitatory), clean_name))) |> 
               dplyr::mutate(is_transcription_factor = (genename_cleaned %in% TF_names_cleaned))

rm(tf_raw, tf, TF_names_cleaned)

# =================== \__subset cells w pert assignment ====================================================

# create grna assignment matrix: Keep GFP col! each row should sum=1
perturb_ = relevel(as.factor(metadata_excitatory$Perturbation), ref = "GFP") # levels with Unpert, GFP, other pert
perturb_ = relevel(perturb_, ref = "UnPert") 
grna = data.frame(
    model.matrix(~ perturb_-1, data = NULL),
    check.names = FALSE
  )
colnames(grna) <- sub("^perturb_", "", colnames(grna))
row.names(grna) = row.names(metadata_excitatory)

# save 
# write.csv(x = grna, file = sprintf('%s/grna.csv', save_dir), row.names = TRUE) # csv 
rm(perturb_)


# =================== \__calculate gene importance (deviance) ====================================================
print(sprintf("[%s]        - calculate gene importance (deviance)", Sys.time()))
# Input should be: row = feature, col = cell
gene_dev = scry::devianceFeatureSelection(object=counts_excitatory, fam='binomial') # < 1 min
gene_dev_df = data.frame(gene_idx = 1:length(gene_dev), deviance=gene_dev)  |> 
              dplyr::arrange(desc(gene_dev)) |> dplyr::mutate(importance_rank = 1:n())  # add rank importance 
gene_dev_df = gene_dev_df[names(gene_dev), ] |> tibble::rownames_to_column(var = 'gene_name') |> # reorder to original
              dplyr::mutate(is_grna_target = (tolower(gene_name) %in%  perturbation_alias_df$Perturbation_alias_lower  )) # add T/F is a perturbation target


gene_dev_df = merge(gene_dev_df, genenames_df, 
                    by.x = 'gene_name', by.y = 'genename_jin', 
                    all.x = TRUE) |> dplyr::arrange(gene_idx)

write.csv(x = gene_dev_df, file = sprintf('%s/gene_dev.csv', save_dir), row.names = FALSE) # csv 

rm(gene_dev)


top_genenames_orgrnatarget = gene_dev_df |> filter(importance_rank <= NUM_IMPORTANT_GENES | is_grna_target) |> arrange(importance_rank) |> pull(gene_name)

gene_counts = counts_excitatory[top_genenames_orgrnatarget, ]
# row.names(gene_counts) = toupper(row.names(gene_counts))
# saveRDS(object = gene_counts, file = sprintf('%s/gene_counts.rds', save_dir)) # quickly get gene counts of top genes




# =================== \__perform normalization ========================================================
print(sprintf("[%s]        - performing normalization", Sys.time()))
# genes in the rows and samples in the columns
gene_norm = scry::nullResiduals(gene_counts, 
                                fam  = "binomial", 
                                type = "deviance") #  genes x  cells
# write.csv(   x = gene_norm, file = sprintf('%s/gene_norm.csv', save_dir)) # csv (dont, this is large)
# saveRDS(object = gene_norm, file = sprintf('%s/gene_norm.rds', save_dir)) # this is smaller (is rds smaller or h5?)
# format(object.size(gene_norm), units = 'Gb') # .2 Gb
# format(object.size(gene_counts), units = 'Gb') # .2 Gb





# =================== \__save all as h5 file =========================================================
print(sprintf("[%s]    - Saving as h5", Sys.time()))
library(h5lite)

h5file <- sprintf('%s/jin_exneu.h5', save_dir)

# metadata about the excitatory neuron cells
print(sprintf("[%s]        + cell_metadata", Sys.time()))
h5_write(metadata_excitatory , h5file, "cell_metadata/cell_metadata")


# grna assignments
print(sprintf("[%s]        + grna", Sys.time()))
grna_wide = t(as.matrix(grna)) # format so consistent with genes: (feature) x cell
h5_write(          grna_wide , h5file, "grna/grna")
h5_write(row.names(grna_wide), h5file, "grna/row_grnanames")
h5_write( colnames(grna_wide), h5file, "grna/col_cellnames")

# gene counts
print(sprintf("[%s]        + counts", Sys.time()))
h5_write(as.matrix(gene_counts), h5file, "counts/counts")        # it seems that python has trouble reading in with extra rownames and columns... 
h5_write(row.names(gene_counts), h5file, "counts/row_genenames") # so create hierarchy and save rownames and colnames separately
h5_write(colnames(gene_counts), h5file, "counts/col_cellnames")  

# gene normalized
print(sprintf("[%s]        + gene_norm", Sys.time()))
h5_write(as.matrix(gene_norm), h5file, "gene_norm/gene_norm")
h5_write(row.names(gene_norm), h5file, "gene_norm/row_genenames")
h5_write( colnames(gene_norm), h5file, "gene_norm/col_cellnames")

# test
# mat_in <- h5_read(h5file, "counts/counts") # R reads in row and colnames
# colnames(gene_counts) == rownames(grna)

# check cells are in the same order
all(row.names(metadata_excitatory) == metadata_excitatory$NAME)
all(row.names(metadata_excitatory) == colnames(gene_counts))
all(row.names(metadata_excitatory) == colnames(gene_norm))
all(row.names(metadata_excitatory) == colnames(grna_wide))


# ============================================================================== #
# =================== Perform PCA + SPCA ======================================
# ============================================================================== #
gc()
# table(perturbation_alias_df$Perturbation_alias_lower %in% tolower(row.names(gene_norm)))  # only missing 1 now
# table(perturbation_alias_df$Perturbation_alias_lower %in% tolower(row.names(gene_counts))) # GFP gene

# # can start here and load in prev saves
# gene_dev_df = read.csv(sprintf('%s/gene_dev.csv', save_dir))
# gene_norm = readRDS(sprintf('%s/gene_norm.rds', save_dir))
# gene_norm <- h5lite::h5_read(h5file, "gene_norm/gene_norm") # R reads in row and colnames


print(sprintf("[%s]    - Perform PCA + SPCA", Sys.time()))

# =================== \__setup ========================================================
print(sprintf("[%s]        - prep (remove TF and grna target genes)", Sys.time()))
# select genes for PCA decomp, remove TF and grna targets which add signal
valid_gene_names = gene_dev_df |> 
  dplyr::filter((!is_transcription_factor) & (!is_grna_target) & importance_rank <= NUM_IMPORTANT_GENES) |>    # no TF or grna target
  dplyr::arrange(gene_idx)


# PCA on all cells and genes? Or PCA on only the 'control' perturbations?
# PCA on all cells for now

dim(gene_norm)
# N_subsample = 5000
N_subsample = ncol(gene_norm) 
set.seed(12345)
gene_norm_validgenes = gene_norm[valid_gene_names$gene_name, ] # 4014 to 3833 genes

if(N_subsample == ncol(gene_norm)) {
  gene_norm_SPC = t(gene_norm_validgenes) # use all cells
} else {
  gene_norm_SPC = t(gene_norm_validgenes[, sample(1:ncol(gene_norm), N_subsample)]) # sample X out of XX cells
}
gene_norm_SPC = as.matrix(gene_norm_SPC)
# gene_norm_SPC = t(gene_norm[myGenenames_df |> dplyr::filter((!TF) & (!grna_target)) |> dplyr::pull(importance_rank),
#                             sample(1:ncol(gene_norm), N_subsample)]) # sample 5k out of 21k cells
dim(gene_norm_SPC) # 8215 cells x  3830 genes


# =================== \__PCA ========================================================
print(sprintf("[%s]        - PCA", Sys.time()))
# perform pca
pca_res = prcomp(x = gene_norm_SPC, rank.=100, retx = TRUE, center = TRUE, scale. = TRUE)

# construct NCs as the actual loadings: [all cells x valid genes] %*% [rotation of valid genes x PCs]
NC_loadings = t(gene_norm_validgenes) %*% pca_res$rotation

# save
dir.create(sprintf('%s/pca/', save_dir), showWarnings = FALSE)
saveRDS(NC_loadings, sprintf('%s/pca/NCloadings.rds', save_dir)) 


# =================== \__Sparse PCA ===========================================================
if(FALSE) { # don't do spca, it takes a long time
print(sprintf("[%s]        - Sparse PCA", Sys.time()))
dir.create(sprintf('%s/spca/', save_dir), showWarnings = FALSE)


# =================== \__\__Use SPC.cv to choose tuning parameters ======================================
print(sprintf("[%s]        - Sparse PCA: cv tuning params", Sys.time()))
# Use SPC.cv to choose tuning parameters:
if(F) {
  # maybe subset here
  # gene_norm_SPC_cv = gene_norm_SPC[sample(nrow(gene_norm_SPC), size=10000, replace=FALSE), ]
  gene_norm_SPC_cv = gene_norm_SPC
  
  cv.out <- PMA::SPC.cv(gene_norm_SPC_cv, 
                   sumabsvs = c(seq(1.2, 5, len = 5), seq(6, floor(sqrt(ncol(gene_norm_SPC_cv))), len = 5)),
                   orth=TRUE)
  print(cv.out)
  plot(cv.out)
  print(cv.out$bestsumabsv)    # 42.75
  print(cv.out$bestsumabsv1se) # 47.25  =smallest sumabsv value that has CV error within 1 SE of best CV error:  30.5 
 
  saveRDS(cv.out, 
          sprintf('%s/spca/cvout.rds', save_dir)) # save cv res
}

# =================== \__\__perform SPCA ======================================

print(sprintf("[%s]        - Sparse PCA: fit", Sys.time()))
my_sumabsv = 47.25
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
}

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








