# ------------------------------------------------------------------------------------------------ #
#     Test basic analysis on jin-2020 dataset subset from causarrayexampledata 
# this file tests code
# ------------------------------------------------------------------------------------------------ #
args = commandArgs(trailingOnly = TRUE)
args = c('macbook')




assertthat::assert_that(length(args) > 0, msg="must give arg for specifying device eg 'Rscript <filename>.R ubergenno'")
DEVICE = args[1]
source('../PATHS.R') # load in data_dir and save_dir and CODE_DIR, depending on DEVICE value
assertthat::assert_that(!is.null(data_dir), msg='first arg must be: laptop, desktop, or ubergenno')

library(dplyr)
library(Seurat)
library(SeuratObject)

save_dir
data_dir

list.files(data_dir)
list.files()




# =================== causarray tutorial by Jinhong ========================================================================

# perturbseq_exneu_rds = readRDS(sprintf('%s/causarrayexampledata/perturbseq-exneu.rds', data_dir))


# data_dir = "/Users/catherinewang/Documents/School/genData/papalexi"

# jincounts = read.delim(sprintf('%s/../jin/other/Counts.PBC.txt', data_dir))
# jincounts$gene |> unique() |> length()
# summary(jincounts)

# sc.seurat <- readRDS(sprintf('%s/../jin/causarrayexampledata/perturbseq-exneu.rds', data_dir))
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



# A is pert assignment of 2926 x 29  #cells x #perturbations?
# Y is counts of 2926 x 3221 of cells x genes?

# A_: include the GFP column
A_ = data.frame(
  model.matrix(~ trt_ - 1, data = perturb),
  check.names = FALSE
)
A_


# ============================================================================== #
# =================== Negative Binomial GLM ====================================
# ============================================================================== #


# A: cell-by-perturbation 0/1 trt assignment, cells that received negative control perturbation GFP have 0
# Y: cell-by-gene matrix of gene expression counts


# sometimes the names are bad for the glm formulas (e.g. AB-C). Remove non letter and number characters 
clean_name <- function(x) {  gsub(pattern='[^a-zA-Z0-9]', replacement='', x=x)  }



# Neg Bin fits without any covariates
libsize = data.frame(libsize_log=log(metadata$nUMI |> as.numeric() + 1e-20))
# libsize_log = log((metadata$nUMI |> as.numeric())+ 1e-20)




# which_genes = colnames(Y) # all genes
N_genes = 1000   # or randomly choose some
set.seed(12345)

perturbations_names = unique(metadata$Perturbation) |> sort()
pert_targetgene_in_subsetdata = intersect(perturbations_names, colnames(Y)) # always test these bc they are pert targets
which_genes = union(pert_targetgene_in_subsetdata, 
                    sample(setdiff(colnames(Y), perturbations_names),  replace = F, size = N_genes - length(pert_targetgene_in_subsetdata)))

# fit and save res: 1 glm per gene (vs individual pert-gene fits)
res = NULL
for(Y_name in which_genes) {
  
  df = cbind(Y[, Y_name, drop=FALSE] |> `colnames<-`(paste0('gene_', clean_name(Y_name))), 
             libsize,
             A
             )
  
  
  # formula = paste0('gene_', Y_name, ' ~ offset(libsize_log) + ', paste0(colnames(A), collapse = ' + '))
  formula = paste0('gene_', clean_name(Y_name), 
                   ' ~ offset(libsize_log) + ', 
                   paste0(sapply(colnames(A), clean_name), collapse = ' + '))
  nbfit = MASS::glm.nb(formula=formula, data = df)
  
  res = rbind(res, 
              summary(nbfit)$coefficients[-1, ] |> 
                `colnames<-`(c('estimate', 'se', 'zvalue', 'pvalue')) |> 
                data.frame() |> 
                dplyr::mutate(grna=colnames(A), gene = Y_name, .after="estimate") )
  
  
  # summary(nbfit)$coefficients[colnames(A), ] |> 
  #   `colnames<-`(c('estimate', 'se', 'zvalue', 'pvalue')) |> 
  #   data.frame() |> 
  #   tibble::rownames_to_column("grna") |> 
  #   dplyr::mutate(gene = Y_name, .after="grna")
}




# ============================================================================== #
# =================== Pre-processing ===========================================
# ============================================================================== #


# =================== Load Transcription Factor (TF) genes list ======================================
print(sprintf("[%s]    - loading list of TF genes", Sys.time()))
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


# =================== calculate gene importance (deviance) ====================================================
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




# =================== Perform normalization ========================================================
print(sprintf("[%s]    - performing normalization", Sys.time()))
# genes in the rows and samples in the columns
gene_norm = scry::nullResiduals(counts, 
                                fam  = "binomial", 
                                type = "deviance") # 3221 genes x 2926 cells
# write.csv(   x = gene_norm, file = sprintf('%s/gene_norm.csv', save_dir)) # csv (dont, this is large)
saveRDS(object = gene_norm, file = sprintf('%s/gene_norm.rds', save_dir)) # this is smaller



# =================== Perform PCA + SPCA ========================================================


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
# rm(gene_norm_validgenes); gc()
# gene_norm_SPC = t(gene_norm[myGenenames_df |> dplyr::filter((!TF) & (!grna_target)) |> dplyr::pull(importance_rank),
#                             sample(1:ncol(gene_norm), N_subsample)]) # sample 5k out of 21k cells
dim(gene_norm_SPC)


# =================== Perform PCA ========================================================
# perform pca
pca_res = prcomp(x = gene_norm_SPC, rank.=100, retx = TRUE, center = TRUE, scale. = TRUE)

# construct NCs as the actual loadings
NC_loadings = t(gene_norm_validgenes) %*% pca_res$rotation

# save
dir.create(sprintf('%s/pca/', save_dir), showWarnings = FALSE)
saveRDS(NC_loadings, sprintf('%s/pca/NCloadings.rds', save_dir)) 


# =================== Perform Sparse PCA ===========================================================
dir.create(sprintf('%s/spca/', save_dir), showWarnings = FALSE)


# =================== +-- Use SPC.cv to choose tuning parameters: ======================================
# Use SPC.cv to choose tuning parameters:
if(F) {
  cv.out <- PMA::SPC.cv(gene_norm_SPC, 
                   sumabsvs = c(seq(1.2, 5, len = 5), seq(6, floor(sqrt(ncol(gene_norm_SPC))), len = 5)),
                   orth=TRUE)
  print(cv.out)
  plot(cv.out)
  print(cv.out$bestsumabsv) # testing out on top 1000 genes sampled 5000 cells: 5 and 4.58 but default values are 1.2-5. Should be between 1 and sqrt(p). which is 1000 here... 31.6
  print(cv.out$bestsumabsv1se)
  
  saveRDS(cv.out, 
          sprintf('%s/spca/cvout.rds', save_dir)) # save cv res
}

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
NC_loadings |> dim()

# check orthogonality (not exactly Identity matrix bc SPC used sample)
t(NC_loadings) %*% NC_loadings / nrow(NC_loadings)



heatmap( t(NC_loadings) %*% NC_loadings / nrow(NC_loadings), 
         Rowv=NA, Colv=NA, col = heat.colors(256),  margins=c(5,10))


# ========================== Raw Counts to Seurat Object ======================
# download raw counts from 
# ========================== 
# Datasets from Jin-2020 experiment downloaded from different sites:
#   
#   broad/ 
#   - from the broad institute portal
#   - only has the log TPM gene expression for all the cells (not raw counts)
# 
#   causarrayexample/
#   - from Jin-Hong's tutorial for causarray package
#   - it is in count form, but only has a subset of the data
# 
#   GEO/
#   - data from GEO website
#   - needs cleaning. And I don't quite understand the formats
# 
# 
# I just want gene expression (in raw counts) and perturbation assignment for all the qc cells...

# ========================== \__causarrayexample  ========================== 
sc.seurat <- readRDS(sprintf('%s/causarrayexampledata/perturbseq-exneu.rds', data_dir))
jin_causarray_gene = GetAssayData(sc.seurat, assay = "RNA", layer = "counts")
jin_causarray_metadata = sc.seurat@meta.data
rm(sc.seurat)
dim(jin_causarray_gene) # 3221 genes x 2926 cells
dim(jin_causarray_grna) # 2926 cells x 13 metadata
head(jin_causarray_gene);head(jin_causarray_grna)

jin_causarray_grna$Batch |> table() # these are the samples 1-18 \ 17


# does NAME always have PertCortex_1? No, the number is the batch#/Sample#
test = sapply(X = jin_causarray_grna$NAME, FUN = function(x) {substring(x, first = 1, last = 13)})
test |> table()

# take an example gene to see how the counts might have been added up
jin_causarray_Cd24a = jin_causarray_gene['Cd24a', ]

  
# ========================== \__broad  ========================== 
library(Matrix)

list.dirs(sprintf('%s/broad/', data_dir))

# no... this is log (or some continuous expression...)
# load in all counts in txt.gz format (this takes a while)
jin <- read.delim(sprintf("%s/broad/expression/expression_PertCortex.new.txt.gz", data_dir), header = T, stringsAsFactors = F)

dim(jin) # 27998 genes x  49068 cells

# first col is gene name
gene_names = jin[, 1]  # should be unique: max(table(gene_names))
jin_sparse = as(jin[, -1], "sparseMatrix") 
rownames(jin_sparse) = gene_names

jin_logTPM_cellnames = colnames(jin_sparse)

sapply(jin_logTPM_cellnames, FUN = function(x) {strsplit(x, split = '_')[[1]][1]}) |> table()
sapply(jin_logTPM_cellnames, FUN = function(x) {strsplit(x, split = '_')[[1]][3]})



# load metadata
jin_meta = read.delim(sprintf('%s/broad/metadata/meta_PertCortex.txt', data_dir), header = T, stringsAsFactors = F) 
jin_meta = jin_meta[-1, ] # 1st row (after header) is col type
rownames(jin_meta) <- jin_meta[, 1]


format(object.size(jin), units = 'Gb') # 10 Gb
format(object.size(jin_sparse), units = 'Gb') # 1 Gb


# try again... other file? this is for perturbations...
jin2 <- read.delim(sprintf("%s/broad/other/Counts.PBC.txt", data_dir), header = T, stringsAsFactors = F)
dim(jin2) #  301349      7
head(jin2)
#                          Name              cbc                PBC Count         batch Num gene
# PertCortex_4_AAACCTGCATGCGCAC AAACCTGCATGCGCAC ACTAAAGCTGCATCGCGG     4 ctx170123.txt   4 Mll1
# PertCortex_4_AAACCTGTCCGCGGTA AAACCTGTCCGCGGTA ACTAAAGCTGCATCGCGG     1 ctx170123.txt   4 Mll1


jin2$gene |> table() |> length()

# ========================== \__GEO  ==========================

# try again... download from GEO this time... but this is probably a small sample... I do not know what this is...
jin3 <- read.delim(sprintf("%s/GEO/GSE157977_RAW/GSM4782532_counts.dialout.sample1.UMI.Counts.csv.gz", data_dir), 
                   header = T, stringsAsFactors = F, sep = ",")



# From https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE157977
#
GEO_data_folder = sprintf('%s/GEO/GSE157977_RAW/', data_dir)
GEO_files = list.files(GEO_data_folder)

for(fn in list.files(GEO_data_folder)) {
  
}

fn = GEO_files[1]

jin4 <- read.delim(sprintf('%s/%s', GEO_data_folder, fn), 
                   header = T, stringsAsFactors = F, sep = ",")


# BiocManager::install("rhdf5")
library(rhdf5)
fn = "GSM4782559_UMI.counts.GFP.h5" 
h5ls(sprintf('%s/%s', GEO_data_folder, fn))

h5_group_names = h5ls(sprintf('%s/%s', GEO_data_folder, fn))
paste0(h5_group_names$group, '/', h5_group_names$name)

mydata <- h5read(sprintf('%s/%s', GEO_data_folder, fn), "//mm10")
names(mydata)

mydata$barcodes
mydata$shape


mydata$data
mydata$gene_names
mydata$indices
mydata$indptr


fn = "GSM4782566_UMI.counts.sample3.h5"
h5ls(sprintf('%s/%s', GEO_data_folder, fn))
# group       name       otype  dclass     dim
# 0     /       mm10   H5I_GROUP                
# 1 /mm10   barcodes H5I_DATASET  STRING  737280
# 2 /mm10       data H5I_DATASET INTEGER 3061611
# 3 /mm10 gene_names H5I_DATASET  STRING   27998
# 4 /mm10      genes H5I_DATASET  STRING   27998
# 5 /mm10    indices H5I_DATASET INTEGER 3061611
# 6 /mm10     indptr H5I_DATASET INTEGER  737281
# 7 /mm10      shape H5I_DATASET INTEGER       2

library(Seurat)
library(SeuratObject)

expression_matrix = Seurat::Read10X_h5(sprintf('%s/%s', GEO_data_folder, fn))

# seu = CreateSeuratObject(counts = expression_matrix, assay = "RNA")
dim(expression_matrix) # 27998 737280 seems to be genes by cells?
expression_matrix |> colnames() |> head()  # idk what these are, cell barcodes?
expression_matrix |> row.names() |> head() # seem to be gene names




# compare 2 samples
fn1 = "GSM4782566_UMI.counts.sample3.h5"
fn2 = "GSM4782567_UMI.counts.sample4.h5"
expression_matrix1 = Seurat::Read10X_h5(sprintf('%s/%s', GEO_data_folder, fn1))
expression_matrix2 = Seurat::Read10X_h5(sprintf('%s/%s', GEO_data_folder, fn2))


dim(expression_matrix1); dim(expression_matrix2) # the same  27998 737280
colnames(expression_matrix1)[1:5] == colnames(expression_matrix2)[1:5]

all(colnames(expression_matrix1) == colnames(expression_matrix2)) # same colnames
all(row.names(expression_matrix1) == row.names(expression_matrix2))  # same rownames

all(expression_matrix1[1:50, 1:50] == expression_matrix2[1:50, 1:50])
(colSums(expression_matrix1) == colSums(expression_matrix2)) |> table()


cur_fns = GEO_files[19:length(GEO_files)] 
for(fn in cur_fns) {
  filepath = sprintf('%s/%s', GEO_data_folder, fn)
  print(filepath)
  expression_matrix = Seurat::Read10X_h5(filepath)
  print(dim(expression_matrix))
  
  sum(expression_matrix != 0)
}




jin_logTPM_cellnames_ATCG = sapply(jin_logTPM_cellnames, FUN = function(x) {strsplit(x, split = '_')[[1]][3]})
expression_matrix_cols_ATCG = sapply(colnames(expression_matrix), FUN = function(x) {strsplit(x, split = '-')[[1]][1]})
all(jin_logTPM_cellnames_ATCG %in% expression_matrix_cols_ATCG) # all the ATCG names of jin_logTPM are in expression_matrix's


# I don't understand... there are 18 samples, + a file for GFP and a file for Ank2
# across these h5 files, they all seem to have the same colnames (cell barcodes)
# Do they re-use the cell barcodes? 


# [1]  "GSM4782532_counts.dialout.sample1.UMI.Counts.csv.gz"  "GSM4782534_counts.dialout.sample2.UMI.Counts.csv.gz" 
# .... 
# [19] "GSM4782559_UMI.counts.GFP.h5"                         "GSM4782561_UMI.counts.Ank2.h5"                        "GSM4782562_UMI.counts.sample1.h5"                    
# [22] "GSM4782564_UMI.counts.sample2.h5"                     "GSM4782566_UMI.counts.sample3.h5"                     "GSM4782567_UMI.counts.sample4.h5"                    
# [25] "GSM4782569_UMI.counts.sample5.h5"                     "GSM4782570_UMI.counts.sample6.h5"                     "GSM4782572_UMI.counts.sample7.h5"                    
# [28] "GSM4782573_UMI.counts.sample8.h5"                     "GSM4782575_UMI.counts.sample9.h5"                     "GSM4782577_UMI.counts.sample10.h5"                   
# [31] "GSM4782578_UMI.counts.sample11.h5"                    "GSM4782580_UMI.counts.sample12.h5"                    "GSM4782582_UMI.counts.sample13.h5"                   
# [34] "GSM4782583_UMI.counts.sample14.h5"                    "GSM4782584_UMI.counts.sample15.h5"                    "GSM4782586_UMI.counts.sample16.h5"                   
# [37] "GSM4782588_UMI.counts.sample17.h5"                    "GSM4782589_UMI.counts.sample18.h5"   


# some criteria for filtering out cells
QCCELL = list(MIN_NONZEROGENES = 500) 



sample_numbers = setdiff(1:18, 17) # somehow sample 17 did not pass qc by them so I will just exclude?


dialouts = list()
dialout_colnames = list()
for(sample_number in sample_numbers) {
  sample_number_filenames = grep(pattern = sprintf('sample%s\\.', sample_number), x = GEO_files, value = TRUE)
  # sample_number_filenames = GEO_files[grepl(pattern = sprintf('sample%s', sample_number), x = GEO_files)]
  dialout_filename = grep('dialout', sample_number_filenames, value = TRUE) # something like GSM4782536_counts.dialout.sample4.UMI.Counts.csv.gz
  counts_filename = grep('h5', sample_number_filenames, value = TRUE)       # something like GSM4782567_UMI.counts.sample4.h5
  # dialout_filename = sprintf('%s/GSM4782532_counts.dialout.sample%s.UMI.Counts.csv.gz', GEO_data_folder, sample_number)
  # counts_filename  = sprintf('%s/GSM4782588_UMI.counts.sample%s.h5', GEO_data_folder, sample_number)
  
  dialout_ = read.delim(sprintf('%s/%s', GEO_data_folder, dialout_filename), header = T, stringsAsFactors = F, sep = ",")
  # counts_  = Seurat::Read10X_h5(sprintf('%s/%s', GEO_data_folder, counts_filename))
  
  dialouts[[sprintf('sample%02.f', sample_number)]] = dialout_
  dialout_colnames[[sprintf('sample%02.f', sample_number)]] = colnames(dialout_)[-1]
  print(dim(dialout_))
  print(colnames(dialout_))
  # rowSums(dialout_[, -1])
  # colnames(counts_) |> table() |> max() # unique colnames (the cell barcodes supposedly)
  
}


all_dialoutcolnames = c()
for(dn in names(dialout_colnames)) {
  all_dialoutcolnames = c(all_dialoutcolnames, dialout_colnames[[dn]])
}



## or... I could just try to convert log(count) --> counts...
# The data from both datasets were loaded into Seurat separately and transformed to log counts per million.
# gene exp is log of TPM (but there are 0s, so cannot be just log(count) and also dont know what is the base)
# ?? what is this. counts --> log(counts)
# log counts per million
# https://www.reneshbedre.com/blog/expression_units.html


GEO_Cd24a = c()

cur_fns = GEO_files[21:length(GEO_files)] 
for(fn in cur_fns) {
  filepath = sprintf('%s/%s', GEO_data_folder, fn)
  print(filepath)
  expression_matrix = Seurat::Read10X_h5(filepath)
  
  print(dim(expression_matrix))
  
  # expression_matrix['Cd24a', 'AAAGCAAAGATGTGTA-1']
  # colnames(expression_matrix)[1:5]
  sapply(colnames(expression_matrix), FUN = function(x){strsplit(x, '-')[[1]][2]}) |> table()  |> print() # seems to be just 1? '<cbc>-1'
  
  # sum(expression_matrix != 0)
}

cbc_causarray = sapply(names(jin_causarray_Cd24a[1:10]), FUN = function(x){strsplit(x=x, split='_')[[1]][3]})


expression_matrix['Cd24a', paste0(cbc_causarray, '-1')]
jin_causarray_Cd24a[1:10]


# ================= \__ combine ==================

# from the broad institute data, take the metadata with the cell names (and perturbation assignment?)
# take NAMES (should have batch/sample info),
# then, extract the counts from the GEO website...
# surely this is not the best way, but I think it works...


# From GEO
GEO_data_folder = sprintf('%s/GEO/GSE157977_RAW/', data_dir)
GEO_files = list.files(GEO_data_folder)
# From Broad
jin_meta = read.delim(sprintf('%s/broad/metadata/meta_PertCortex.txt', data_dir), header = T, stringsAsFactors = F) 
jin_meta = jin_meta[-1, ] # 1st row (after header) is col type
rownames(jin_meta) <- jin_meta[, 1]

# prepare cbc names- separate the names into batch/sample number and cbc
cell_names = jin_meta$NAME
cell_batch = jin_meta$Batch # should be same
cell_samplenumber = sapply(X=cell_names, FUN=function(x){strsplit(x=x,split='_')[[1]][2]}) 
cell_cbc          = sapply(X=cell_names, FUN=function(x){strsplit(x=x,split='_')[[1]][3]}) 
cell_names_df = data.frame(NAME=cell_names, Batch=cell_batch, 
                           samplenumber=cell_samplenumber, barcode=cell_cbc)

(cell_names_df$Batch == cell_names_df$samplenumber) |> table() # check all same


jin_counts_assembled = NULL
for(cur_samplenumber in c(1:16, 17, 18)) { 
  # filter for current sample number
  cur_cell_names_df = cell_names_df |> filter(samplenumber == cur_samplenumber)
  cbc_formatted = paste0(cur_cell_names_df$barcode, '-1') # it is formatted differently
  
  
  # find the filename: has 'sampleXX.h5' in the name
  fn = grep(paste0('sample', cur_samplenumber, '.h5'), GEO_files, value = TRUE)
  print(fn)
  assertthat::assert_that(length(fn) == 1, msg = 'GEO h5 file found too many or few files (only should have 1 match)')
  
  # load in raw counts from GEO 
  expression_matrix = Seurat::Read10X_h5( sprintf('%s/%s', GEO_data_folder, fn))
  
  assertthat::assert_that(all(cbc_formatted %in% colnames(expression_matrix)), 
                          msg = 'selected cbc must be in the loaded GEO h5 file')
  
  counts_cur = expression_matrix[, cbc_formatted]
  colnames(counts_cur) = cur_cell_names_df$NAME
  
  jin_counts_assembled = cbind(jin_counts_assembled, counts_cur)
  # rm(cur_cell_names_df, cbc_formatted, fn, expression_matrix, counts_cur)
}

saveRDS(jin_counts_assembled, sprintf('%s/combineBroadGEO/gene_useCellNames.rds', data_dir))









jin_counts_assembled
dim(jin_counts_assembled)
jin_meta

jin_counts_assembled['Srgap2', jin_meta$NAME[1:15]]



jin_causarray_gene['Cd24a', colnames(jin_causarray_gene)[1:15]]
jin_counts_assembled['Cd24a', colnames(jin_causarray_gene)[1:15]]

test_cells_compare = intersect(colnames(jin_causarray_gene), colnames(jin_counts_assembled))
(jin_causarray_gene['Cd24a', test_cells_compare] == jin_counts_assembled['Cd24a', test_cells_compare]) |> table()
(jin_causarray_gene['Cd24a', test_cells_compare] == jin_counts_assembled['Cd24a', test_cells_compare]) |> which()

test_cells_compare[(jin_causarray_gene['Cd24a', test_cells_compare] != jin_counts_assembled['Cd24a', test_cells_compare])]
test_cells_bad = test_cells_compare[(jin_causarray_gene['Cd24a', test_cells_compare] != jin_counts_assembled['Cd24a', test_cells_compare])]

# don't seem to match on sample 3. load in expression matrix for sample 3
expression_matrix = Seurat::Read10X_h5( sprintf('%s/%s', GEO_data_folder, "GSM4782566_UMI.counts.sample3.h5"))
expression_matrix['Cd24a', "TACCTTACACCGAAAG-1"]
jin_causarray_gene['Cd24a', "PertCortex_3_TACCTTACACCGAAAG"]

test_cells_bad_formatted = paste0(sapply(test_cells_bad, FUN=function(x){strsplit(x, split='_')[[1]][3]}), '-1')
expression_matrix['Cd24a', test_cells_bad_formatted]
jin_causarray_gene['Cd24a', test_cells_bad]
jin_counts_assembled['Cd24a', test_cells_bad]



(expression_matrix['Cd24a', test_cells_bad_formatted] == jin_causarray_gene['Cd24a', test_cells_bad]) |> table()


# later
all(colnames(jin_causarray_gene) %in% colnames(jin_counts_assembled))



# GFP needs to be added differently...? oh it seems to be the same? why is there a separate GFP (and Ank2 file?)
jin_meta |> filter(Perturbation == 'GFP')




# ================= \__different? ==================

# There is some discrepancy between the assembled data and jinhong causarray example??


# clean environment and then test


# load causarray example
sc.seurat <- readRDS(sprintf('%s/causarrayexampledata/perturbseq-exneu.rds', data_dir))
jin_causarray_gene = GetAssayData(sc.seurat, assay = "RNA", layer = "counts")
jin_causarray_metadata = sc.seurat@meta.data
rm(sc.seurat)

# load causarray example
jin_assembled = readRDS(sprintf('%s/combineBroadGEO/counts.rds', data_dir))


# causarray example should be a subset

all(row.names(jin_causarray_gene) %in% row.names(jin_assembled)) # all genes contained
all(colnames(jin_causarray_gene) %in% colnames(jin_assembled)) # all cells contained


jin_assembled_subset = jin_assembled[row.names(jin_causarray_gene), colnames(jin_causarray_gene) ]

dim(jin_causarray_gene); dim(jin_assembled_subset)
(jin_causarray_gene == jin_assembled_subset) |> dim()

(jin_causarray_gene == jin_assembled_subset) |> as.vector() |> table()
#     FALSE      TRUE 
# 5,949,472 3,475,174 
# a lot of discrepancies...



# check matches
match_idx = (jin_causarray_gene == jin_assembled_subset) |> as.vector() |> which()
jin_causarray_gene[match_idx]
jin_assembled_subset[match_idx]
all(jin_causarray_gene[match_idx] == jin_assembled_subset[match_idx])

# check mismatches
mismatch_idx = (jin_causarray_gene != jin_assembled_subset) |> as.vector() |> which()
jin_causarray_gene[mismatch_idx]
jin_assembled_subset[mismatch_idx]



jin_causarray_gene[mismatch_idx[1:15]]
jin_assembled_subset[mismatch_idx[1:15]]


table(jin_assembled_subset[mismatch_idx] == 0)
table(jin_assembled_subset[match_idx] == 0)

mismatch_idx = which(as.matrix(jin_causarray_gene != jin_assembled_subset), arr.ind = TRUE) |> data.frame()
mismatch_idx$col |> hist()
mismatch_idx$row |> hist()

# take a cell with a large number of mismatches
mismatch_idx |> group_by(col) |> summarize(count = n()) |> arrange(desc(count)) |> head()

colnames(jin_causarray_gene)[492]

test_badcell = "PertCortex_4_CCCAATCCAAGTTGTC"

# should be from sample 4, so load in:
expression_matrix = Seurat::Read10X_h5( sprintf('%s/%s', data_dir, "GEO/GSE157977_RAW/GSM4782567_UMI.counts.sample4.h5"))

# (expression_matrix[row.names(jin_causarray_gene), paste0(strsplit(test_badcell, '_')[[1]][3], '-1')] == 
# jin_causarray_gene[row.names(jin_causarray_gene), test_badcell]) |> table()


badcell_expr = data.frame(geo       = expression_matrix[row.names(jin_causarray_gene), paste0(strsplit(test_badcell, '_')[[1]][3], '-1')], 
                          causarray = jin_causarray_gene[row.names(jin_causarray_gene), test_badcell])

hist(badcell_expr$geo)
range(badcell_expr$geo)
expression_matrix[, ]
# in GEO, this cell "PertCortex_4_CCCAATCCAAGTTGTC" has all 0 expression counts, so it should not have passed quality control
# (cells with at least XX > Non-zero gene expression counts)
# (this is true both in the subselection of genes and all the genes)
expression_matrix[row.names(jin_causarray_gene), paste0(strsplit(test_badcell, '_')[[1]][3], '-1')] |> range()
expression_matrix[, paste0(strsplit(test_badcell, '_')[[1]][3], '-1')] |> range()

# see if this cell is selected in the broad data, it is and it has a lot of expression...
jin_meta = read.delim(sprintf('%s/broad/metadata/meta_PertCortex.txt', data_dir), header = T, stringsAsFactors = F) 
jin_meta = jin_meta[-1, ] # 1st row (after header) is col type
rownames(jin_meta) <- jin_meta[, 1]

jin_meta[test_badcell, ]
#                                                        NAME nGene  nUMI Cluster Batch   CellType Perturbation isKey isAnalysed           SCRUBLET
# PertCortex_4_CCCAATCCAAGTTGTC PertCortex_4_CCCAATCCAAGTTGTC  7913 55494      16     4 Excitatory         Mll1  TRUE       TRUE 0.0328883099967437

test =  read.delim( sprintf('%s/%s', data_dir, "GEO/GSE157977_RAW/GSM4782536_counts.dialout.sample4.UMI.Counts.csv.gz"), 
                    header = T, stringsAsFactors = F, sep = ",")

test |> dim()




jin_causarray_metadata[test_badcell, ]




rowSums(expression_matrix) |> hist(xlim = c(0, 20), breaks = seq(0, 1000000, by = 1)) 
colSums(expression_matrix) |> hist(xlim = c(0, 20), breaks = seq(0, 1000000, by = 1)) # seems really lowly expressed


GEO_h5_filenames = c("GSM4782562_UMI.counts.sample1.h5",                     "GSM4782564_UMI.counts.sample2.h5",
                     "GSM4782566_UMI.counts.sample3.h5",                     "GSM4782567_UMI.counts.sample4.h5",
                     "GSM4782569_UMI.counts.sample5.h5",                     "GSM4782570_UMI.counts.sample6.h5",
                     "GSM4782572_UMI.counts.sample7.h5",                     "GSM4782573_UMI.counts.sample8.h5",
                     "GSM4782575_UMI.counts.sample9.h5",                     "GSM4782577_UMI.counts.sample10.h5",
                     "GSM4782578_UMI.counts.sample11.h5",                    "GSM4782580_UMI.counts.sample12.h5",
                     "GSM4782582_UMI.counts.sample13.h5",                    "GSM4782583_UMI.counts.sample14.h5",
                     "GSM4782584_UMI.counts.sample15.h5",                    "GSM4782586_UMI.counts.sample16.h5",
                     # "GSM4782588_UMI.counts.sample17.h5",                   
                     "GSM4782589_UMI.counts.sample18.h5")

badcell_total = matrix(nrow = length(GEO_h5_filenames), ncol = 27998) # I don't think these cells are comparable, they have the same cell barcode but they are different batches??
i = 1
for(fn in GEO_h5_filenames) {
  filepath = sprintf('%s/GEO/GSE157977_RAW/%s', data_dir, fn)
  print(filepath)
  expression_matrix = Seurat::Read10X_h5(filepath)
  
  # print(dim(expression_matrix))
  # colSums(expression_matrix) |> hist(xlim = c(0, 20), breaks = seq(0, 1000000, by = 1), main = fn)
  badcell_total[i, ] = expression_matrix[, paste0(strsplit(test_badcell, '_')[[1]][3], '-1')]
  i = i + 1
}





badcell_total[, 1:5]


rowSums(badcell_total)
jin_causarray_metadata$Batch |> unique()
# the nUMI 55494 matches with batch/sample 6?? which does not match with the numbering? 
# e.g. jin_meta[test_badcell, ] has
#                          NAME nGene  nUMI Cluster Batch   CellType Perturbation isKey isAnalysed           SCRUBLET
# PertCortex_4_CCCAATCCAAGTTGTC  7913 55494      16     4 Excitatory         Mll1  TRUE       TRUE 0.0328883099967437
# 

expression_matrix6 = Seurat::Read10X_h5( sprintf('%s/%s', data_dir, "GEO/GSE157977_RAW/GSM4782570_UMI.counts.sample6.h5"))


badcell_expr = data.frame(geo6     = expression_matrix6[row.names(jin_causarray_gene), paste0(strsplit(test_badcell, '_')[[1]][3], '-1')], 
                          causarray = jin_causarray_gene[row.names(jin_causarray_gene), test_badcell])



badcell_expr



# so it seems that the batch numbers/sample numbers are not exactly aligned...?\
# but Jin-Hong got the correct values?


jin_meta$Cluster |> unique() |> as.numeric() |> sort()
jin_meta$Batch |> unique() |> as.numeric() |> sort()
jin_meta  |> pull(Batch) |> as.numeric() |>  hist(breaks = seq(0, 20, by = 1))

jin_meta |> filter(isAnalysed=='TRUE')  |> pull(Batch) |> as.numeric() |>  hist(breaks = seq(0, 20, by = 1)-.5)

# test if calculated nUMIs matches the metadata
jin_assembled_nUMIs = jin_assembled |> colSums()
jin_meta
jin_assembled_nUMIs

table(jin_meta[names(jin_assembled_nUMIs), 'nUMI'] == jin_assembled_nUMIs)


jin_meta[which(!( row.names(jin_meta)  %in% colnames(jin_assembled))), ]



# Or just go through each sample and filter for cells ourselves...



GEO_dialout_filenames = c( "GSM4782532_counts.dialout.sample1.UMI.Counts.csv.gz"  ,"GSM4782534_counts.dialout.sample2.UMI.Counts.csv.gz" ,
                           "GSM4782535_counts.dialout.sample3.UMI.Counts.csv.gz"  ,"GSM4782536_counts.dialout.sample4.UMI.Counts.csv.gz" ,
                           "GSM4782538_counts.dialout.sample5.UMI.Counts.csv.gz"  ,"GSM4782539_counts.dialout.sample6.UMI.Counts.csv.gz" ,
                           "GSM4782541_counts.dialout.sample7.UMI.Counts.csv.gz"  ,"GSM4782542_counts.dialout.sample8.UMI.Counts.csv.gz" ,
                           "GSM4782544_counts.dialout.sample9.UMI.Counts.csv.gz"  ,"GSM4782545_counts.dialout.sample10.UMI.Counts.csv.gz",
                            "GSM4782546_counts.dialout.sample11.UMI.Counts.csv.gz", "GSM4782548_counts.dialout.sample12.UMI.Counts.csv.gz",
                            "GSM4782549_counts.dialout.sample13.UMI.Counts.csv.gz", "GSM4782551_counts.dialout.sample14.UMI.Counts.csv.gz",
                            "GSM4782553_counts.dialout.sample15.UMI.Counts.csv.gz", "GSM4782554_counts.dialout.sample16.UMI.Counts.csv.gz",
                            "GSM4782556_counts.dialout.sample17.UMI.Counts.csv.gz", "GSM4782558_counts.dialout.sample18.UMI.Counts.csv.gz")

GEO_h5_filenames = c("GSM4782562_UMI.counts.sample1.h5",                     "GSM4782564_UMI.counts.sample2.h5",
                     "GSM4782566_UMI.counts.sample3.h5",                     "GSM4782567_UMI.counts.sample4.h5",
                     "GSM4782569_UMI.counts.sample5.h5",                     "GSM4782570_UMI.counts.sample6.h5",
                     "GSM4782572_UMI.counts.sample7.h5",                     "GSM4782573_UMI.counts.sample8.h5",
                     "GSM4782575_UMI.counts.sample9.h5",                     "GSM4782577_UMI.counts.sample10.h5",
                     "GSM4782578_UMI.counts.sample11.h5",                    "GSM4782580_UMI.counts.sample12.h5",
                     "GSM4782582_UMI.counts.sample13.h5",                    "GSM4782583_UMI.counts.sample14.h5",
                     "GSM4782584_UMI.counts.sample15.h5",                    "GSM4782586_UMI.counts.sample16.h5",
                     # "GSM4782588_UMI.counts.sample17.h5",                   
                     "GSM4782589_UMI.counts.sample18.h5")




# =============== \__QC cells from GEO ========


GEO_data_folder = sprintf('%s/GEO/GSE157977_RAW/', data_dir)
GEO_files = list.files(GEO_data_folder)


GENE_MIN_NONZERO = 500

gene_matrix_assembled = NULL
grna_meta = NULL

for(cur_samplenumber in c(1:16, 18)) {
# for(cur_samplenumber in c(1:3)) { 
  # find the filename: 
  fn_gene = grep(paste0('sample', cur_samplenumber, '\\.h5'), GEO_files, value = TRUE) # has 'sampleXX.h5' in the name
  fn_grna = grep(paste0('dialout\\.sample', cur_samplenumber, '\\.UMI'), GEO_files, value = TRUE)
  
  assertthat::assert_that(length(fn_gene) == 1, msg = 'GEO h5 file found too many or few files (only should have 1 match)')
  assertthat::assert_that(length(fn_grna) == 1, msg = 'GEO dialout file found too many or few files (only should have 1 match)')
  
  
  print(sprintf('%s %s', fn_gene, fn_grna))
  # load in raw counts from GEO 
  gene_matrix = Seurat::Read10X_h5( sprintf('%s/%s', GEO_data_folder, fn_gene))
  grna_matrix = read.delim(sprintf('%s/%s', GEO_data_folder, fn_grna), header = T, stringsAsFactors = F, sep = ",")
  row.names(grna_matrix) = grna_matrix$cbc
  grna_matrix = grna_matrix |> dplyr::select(-cbc)
  
  # from genes, filter cells with <500 genes expressed 
  nUMIs = colSums(gene_matrix)
  n_nonzero = colSums(gene_matrix != 0)
  
  chosen_cells = names(n_nonzero)[n_nonzero >= GENE_MIN_NONZERO]
  chosen_cells_cbc = sapply(chosen_cells, FUN=function(x){strsplit(x, '-')[[1]][1]})
  chosen_cells_formatted = paste0('sample_', cur_samplenumber, '_', chosen_cells_cbc)
  
  
  gene_matrix_subset = gene_matrix[, chosen_cells]
  colnames(gene_matrix_subset) = chosen_cells_formatted
  
  gene_matrix_assembled = cbind(gene_matrix_assembled, gene_matrix_subset)
  
  # cells were assigned to a pertubation if they had UMI supporting that perturbation 
  # and the number of UMI supporting it were >1.3 times the next highest.
  
  grna_matrix_subset = grna_matrix[chosen_cells_cbc, ] |> as.matrix()
  
  assign_perturbation <- function(pert_umi) {
    top_two = sort(pert_umi, decreasing = TRUE, na.last = TRUE)[c(1, 2)]
    if(is.na(top_two[1]) | is.na(top_two[2])) {
      return(NA)
    } else if(top_two[1] > 1.3*top_two[2]) {
      return(names(top_two[1]))
    } else {
      return(NA)
    }
  }
  # test = apply(X = grna_matrix_subset[1:3, ], MARGIN = 1, FUN = assign_perturbation)
  # names(test) == row.names(grna_matrix_subset[1:3, ])
  grna_meta_cur = data.frame(cbc = chosen_cells_cbc, 
                             NAME = chosen_cells_formatted, 
                             samplenumber = cur_samplenumber,
                             Perturbation_PBC = apply(X = grna_matrix_subset, MARGIN = 1, FUN = assign_perturbation),
                             nUMI = nUMIs[chosen_cells], 
                             n_nonzero = n_nonzero[chosen_cells], row.names = NULL)

  
  grna_meta = rbind(grna_meta, grna_meta_cur)
}


# get PBC to Gene map
jin_counts_PBC = read.delim(sprintf('%s/broad/other/Counts.PBC.txt', data_dir), 
                            header = T, stringsAsFactors = F, sep = "\t")

# add this Perturbation gene name (nice name) to grna_meta
PBC_map = jin_counts_PBC |> select(PBC, gene) |> distinct()


grna_meta = merge(grna_meta, PBC_map, by.x = 'Perturbation_PBC', by.y = 'PBC', all.x =TRUE) |> 
       rename(Perturbation=gene) |> arrange(samplenumber, cbc) |> relocate(Perturbation_PBC, .before=Perturbation)

# this is more QC GEO, and then add in info about PBC from Broad other data
write.csv(grna_meta, sprintf('%s/combineBroadGEO/grna.csv', data_dir), row.names = FALSE)
write.csv(gene_matrix_assembled, sprintf('%s/combineBroadGEO/gene.csv', data_dir), row.names = FALSE)
saveRDS(gene_matrix_assembled, sprintf('%s/combineBroadGEO/gene.rds', data_dir))






