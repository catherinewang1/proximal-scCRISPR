# ------------------------------------------------------------------------------------------------ #
#     Process the downloaded gene counts and perturbation counts from GEO website
# (also uses the file /other/Counts.PBC.txt from BROAD site to match perturbation barcode to nice gene name)
# 
# Requires:
#   - GEO: downloaded and unzipped files from Jin et al., 2020 GEO files
#     https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE157977
# 
#     at location <data_dir>/GEO/GSE157977_RAW/...
#   - Paper Supplemental File:
#     https://www.science.org/doi/suppl/10.1126/science.aaz6063/suppl_file/aaz6063_jin_table-s5.csv
#     at location <data_dir>/paper_supplement/aaz6063_jin_table-s5.csv'
# 
#   - NO CHANGED: BROAD: /other/Counts.PBC.txt
#     https://singlecell.broadinstitute.org/single_cell/study/SCP1184/in-vivo-perturb-seq-reveals-neuronal-and-glial-abnormalities-associated-with-asd-risk-genes#study-download
# 
#     at location <data_dir>/broad/other/Counts.PBC.txt
# 
# Constructs and saves:
#   - <data_dir>/processGEO/gene.csv
#   - <data_dir>/processGEO/gene.rds
#   - <data_dir>/processGEO/grna.csv
# 
# Structure:
#   - Load in cell-level gene expression in the samples {1-16, 18} in the h5 files like:
#          e.g. sample 8: GSM4782573_UMI.counts.sample8.h5
#   - Only keep cells with at least 500 genes expressed
#   - For those cells, get the perturbation UMI information in the csv.gz files like:
#          e.g. sample 8: GSM4782542_counts.dialout.sample8.UMI.Counts.csv.gz
#   - Assign a Perturbation if the largest perturbation UMI is >=1.3x the next highest count
#   - Assemble a dataframe of
#       - gene: #genes x #cells
#       - grna: #cells x #meta data (including Perturbation assignment)
# 
# From the paper: After quality control, we retained for further analysis a total of
#         46,770 neocortical cells across 17 high quality experimental batches.
# 
# Problems that occurred when performing this:
# 
# The Broad Institute data had the gene x cell and grna assignment information, but the gene information was available
# as log TPM (continuous) values and not raw counts. The GEO website seemed to have the original data of gene counts and perturbation
# counts across 18 samples (they remove sample 17 for quality control). So I wanted to get the cells selected from Broad website,
# and then extract the count gene expression from the GEO data. The Broad data seems to label the cells (and also has meta data about batches)
# based on the sample number (cell names and batches also range from 1-18). However, when trying to select these cells based on these
# batch/sample numbers, they do not always match. They often do match (checking nUMIs and causarray example data), but also do not match for chunks at a time.
# 
# e.g. cell name PertCortex_4_CCCAATCCAAGTTGTC has nUMI 55494. But from GEO, sample 4's cell wtih
# cell barcode CCCAATCCAAGTTGTC has all 0 values. However, sample 6's cell with this barcode has total UMI count 55494.
# This makes it seem like there are 1-18 batches and 1-18 samples, but the numbers do not match up.
# 
# This means I can't use the Broad data's perturbation assignment because the cell names will not match.
# 
# So, this is why I just take the GEO data's cell expression by sample, filter for cells, and assign perturbation again.
# 
# Somehow, Jin-Hong's causarray data with counts does match with Broad's cell name. So there is something I am missing.
# 
# ------------------------------------------------------------------------------------------------ #

# =============== Setup directories ========
args = commandArgs(trailingOnly = TRUE)
args = c('macbook')




assertthat::assert_that(length(args) > 0, msg="must give arg for specifying device eg 'Rscript <filename>.R ubergenno'")
DEVICE = args[1]
source('../PATHS.R') # load in data_dir and save_dir and CODE_DIR, depending on DEVICE value
assertthat::assert_that(!is.null(data_dir), msg='first arg must be: laptop, desktop, or ubergenno')

library(dplyr)
library(Seurat)
library(SeuratObject)

dir.create(sprintf('%s/processGEO/', data_dir)) # create dir to save grna.csv, gene.csv, gene.rds


# =============== Load Extra/Helper files ========

# no?? for some reason, in the data from broad, a perturbation barcode maps to 2 genes??
#   Perturbation_barcode gene
# 1   ACTAAAGCTGCATCGCGG Mll1
# 2   ACTAAAGCTGCATCGCGG Ank2
if(F) {
  # From Broad data, use this Counts.PBC.txt file to map the 
  # Perturbation BarCode (PBC) to the grna target gene name (nice format)
  
  # get PBC to Gene map
  jin_counts_PBC = read.delim(sprintf('%s/broad/other/Counts.PBC.txt', data_dir), 
                              header = T, stringsAsFactors = F, sep = "\t")
  
  # add this Perturbation gene name (nice name) to grna_meta
  PBC_map = jin_counts_PBC |> rename(Perturbation_barcode = PBC) |> select(Perturbation_barcode, gene) |> distinct()
  rm(jin_counts_PBC)
  
  
  PBC_map |> filter(Perturbation_barcode == 'ACTAAAGCTGCATCGCGG') # there are 2 genes mapped to this PBC
  
  
  # 
  # jin_counts_PBC |> filter(PBC == 'ACTAAAGCTGCATCGCGG') |> select(PBC, gene) |> table()
  # 
  # jin_counts_PBC |> filter(PBC == 'ACTAAAGCTGCATCGCGG' & gene == 'Ank2') |> head(10)
  # 
  # jin_counts_PBC |> filter(PBC == 'ACTAAAGCTGCATCGCGG' & gene == 'Mll1') |> head(10)
}


# Instead, use this supplemental file 5 (from published paper's supplemental files in Science)
# https://www.science.org/doi/suppl/10.1126/science.aaz6063/suppl_file/aaz6063_jin_table-s5.csv
# Supplemental File 5 seems to have the perturbation barcodes
supp_file_5 = read.csv(sprintf('%s/paper_supplement/aaz6063_jin_table-s5.csv', data_dir), )
PBC_map = supp_file_5

colnames(PBC_map) = PBC_map[1, ]
PBC_map = PBC_map[-1, ]
PBC_map = PBC_map |> rename(Perturbation_barcode = `Perturbation barcode`, Perturbation = gene) |> select(Perturbation, Perturbation_barcode)
PBC_map$Perturbation = sapply(PBC_map$Perturbation, FUN=function(x){ if(x=='GFP (control)') {'GFP'} else {x}}) # clean the GFP name
row.names(PBC_map) = NULL


assertthat::assert_that(PBC_map$Perturbation_barcode |> table() |> max() == 1) # make sure barcode to gene name unique




# =============== QC cells from GEO ========
GEO_data_folder = sprintf('%s/GEO/GSE157977_RAW/', data_dir)
GEO_files = list.files(GEO_data_folder)
GENE_MIN_NONZERO = 500 # cells must have at least this many genes expressed 


gene_matrix_assembled = NULL
grna_meta = NULL
for(cur_samplenumber in c(1:16, 18)) {
  # for(cur_samplenumber in c(1:3)) { 
  
  # find the filenames of gene UMI and pert UMIs: 
  fn_gene = grep(paste0('sample', cur_samplenumber, '\\.h5'), GEO_files, value = TRUE) # has 'sampleXX.h5' in the name
  fn_grna = grep(paste0('dialout\\.sample', cur_samplenumber, '\\.UMI'), GEO_files, value = TRUE)
  
  assertthat::assert_that(length(fn_gene) == 1, msg = 'GEO h5 file found too many or few files (only should have 1 match)')
  assertthat::assert_that(length(fn_grna) == 1, msg = 'GEO dialout file found too many or few files (only should have 1 match)')
  
  
  print(sprintf('%s %s', fn_gene, fn_grna))
  
  
  # load in raw gene counts from GEO 
  gene_matrix = Seurat::Read10X_h5( sprintf('%s/%s', GEO_data_folder, fn_gene))
  grna_matrix = read.delim(sprintf('%s/%s', GEO_data_folder, fn_grna), header = T, stringsAsFactors = F, sep = ",")
  row.names(grna_matrix) = grna_matrix$cbc
  grna_matrix = grna_matrix |> dplyr::select(-cbc)
  
  # from genes, select cells with >=500 genes expressed 
  nUMIs     = colSums(gene_matrix)
  n_nonzero = colSums(gene_matrix != 0)
  
  # chosen_cells = names(    nUMIs)[    nUMIs >= GENE_MIN_NONZERO] # filter on nUMIs
  chosen_cells = names(n_nonzero)[n_nonzero >= GENE_MIN_NONZERO]  # filter on n non-zero genes
  
  chosen_cells_cbc = sapply(chosen_cells, FUN=function(x){strsplit(x, '-')[[1]][1]})  # extract just cbc 
  chosen_cells_formatted = paste0('sample_', cur_samplenumber, '_', chosen_cells_cbc) # format cell name
  
  grna_matrix_subset = grna_matrix[chosen_cells_cbc, ] |> as.matrix()
  gene_matrix_subset = gene_matrix[, chosen_cells]
  colnames(gene_matrix_subset) = chosen_cells_formatted
  
  gene_matrix_assembled = cbind(gene_matrix_assembled, gene_matrix_subset)
  
  
  # cells were assigned to a perturbation if they had UMI supporting that perturbation 
  # and the number of UMI supporting it were >1.3 times the next highest.
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
                             Perturbation_barcode = apply(X = grna_matrix_subset, MARGIN = 1, FUN = assign_perturbation),
                             nUMI = nUMIs[chosen_cells], 
                             n_nonzero = n_nonzero[chosen_cells], row.names = NULL)
  
  
  grna_meta = rbind(grna_meta, grna_meta_cur)
}

# different number of cells than paper... maybe it is nUMIs and not num nonzero genes (e.g. sum(vec != 0) vs sum(vec))

dim(gene_matrix_assembled); dim(grna_meta)
# if QC on num non-zero:
# [1] 27998 43954
# [1] 43954     6
# if QC on nUMIs
# [1] 27998 51116
# [1] 51116     6


# grna_meta_ = 
grna_meta = merge(grna_meta, PBC_map, by.x = 'Perturbation_barcode', by.y = 'Perturbation_barcode', all.x =TRUE) |> 
  arrange(samplenumber, cbc) |> relocate(Perturbation_barcode, .before=Perturbation)


# save grna as csv and gene as both csv (more compatible) and rds (smaller and faster to load), and then add in info about PBC from Broad other data
write.csv(grna_meta, sprintf('%s/processGEO/grna.csv', data_dir), row.names = FALSE)
saveRDS(gene_matrix_assembled, sprintf('%s/processGEO/gene.rds', data_dir))
write.csv(gene_matrix_assembled, sprintf('%s/processGEO/gene.csv', data_dir), row.names = TRUE)



# number of cells with a perturbation assignment
grna_meta |> filter(!is.na(Perturbation)) |> nrow() 
# if QC by n_nonzero genes: 22085 of 43954 cells have perturbation
# if QC by nUMIs: 24339 of 51116 cells have perturbation
# none of these match...the papers' 46,770


# dim(grna_meta)
# dim(grna_meta_)
# 
# 
# head(grna_meta)
# head(grna_meta_)
# 
# # check for uniqueness
# grna_meta$NAME |> table() |> max()
# grna_meta_$NAME |> table() |> max()
# grna_meta_$NAME |> table()  |> hist()
# grna_meta_$NAME |> table() |> sort(decreasing = TRUE)
# 
# grna_meta_ |> group_by(NAME, Perturbation_barcode) |> summarize(count = n()) |> filter(count == 2)|> arrange(NAME)
# 
# 
# grna_meta |> filter(NAME == 'sample_10_ACGGCCAAGGTGATTA')
# grna_meta_ |> filter(NAME == 'sample_10_ACGGCCAAGGTGATTA')
# PBC_map |> filter(PBC == 'ACTAAAGCTGCATCGCGG')


# =============== Combine Broad and GEO =================

dir.create(sprintf('%s/combineBroadGEO/', data_dir))

# From GEO
GEO_data_folder = sprintf('%s/GEO/GSE157977_RAW/', data_dir)
GEO_files = list.files(GEO_data_folder)


# From Broad
jin_meta = read.delim(sprintf('%s/broad/metadata/meta_PertCortex.txt', data_dir), header = T, stringsAsFactors = F) 
jin_meta = jin_meta[-1, ] # 1st row (after header) is col type
rownames(jin_meta) <- jin_meta[, 1]

# add cbc to jin_meta by parsing the NAME
jin_meta = jin_meta |> mutate(cbc = sapply(X=NAME, FUN=function(x){strsplit(x=x,split='_')[[1]][3]})) 



# =============== \__Match Batch Number and Sample Number =================

# Try to match the batch (from Broad) to the sample number (from GEO) but checking the nUMIs
# Easier to go through sample number first (because it takes longer to load gene expr)
sample_batch_map = list() # x[[sample number]] = batch number
for(cur_samplenumber in 1:18) {
  # cur_samplenumber = 1 # debug
  # find the filename: has 'sampleXX.h5' in the name
  fn = grep(paste0('sample', cur_samplenumber, '.h5'), GEO_files, value = TRUE)
  print(fn)
  assertthat::assert_that(length(fn) == 1, msg = 'GEO h5 file found too many or few files (only should have 1 match)')
  
  # load in raw counts from GEO 
  expression_matrix = Seurat::Read10X_h5( sprintf('%s/%s', GEO_data_folder, fn))
  
  
  for(cur_batchnumber in 1:18) {
    batch_cells = jin_meta |> filter(Batch == cur_batchnumber) |> select(cbc, nUMI)
    batch_cells_cbc_formatted = paste0(batch_cells$cbc, '-1') # it is formatted differently
    
    
    if(all(batch_cells_cbc_formatted %in% colnames(expression_matrix))) {
      # count the total nUMIs in the loaded gene expr for this subset of cells
      cur_nUMIs = colSums(expression_matrix[, batch_cells_cbc_formatted])
      
      # if nUMIs match, then match this sample number with this batch number
      if(all(cur_nUMIs == batch_cells$nUMI)) {
        sample_batch_map[[as.character(cur_samplenumber)]] = cur_batchnumber
      } # keep looping through in case something else goes wrong (this means that there may be one to many)
    } # if not all the cells' barcode are present, then skip this batch
    
    
    rm(batch_cells, batch_cells_cbc_formatted)
  }
  rm(fn, expression_matrix, cur_batchnumber)
}
rm(cur_samplenumber)


batch_sample_df = data.frame(batch = unlist(sample_batch_map),
                             sample = as.integer(names(sample_batch_map)))


write.csv(x = batch_sample_df, file = sprintf('%s/combineBroadGEO/broadbatch_geosample_match.csv', data_dir), row.names = FALSE)



# =============== \__Gather raw cell counts =================

if(!'batch_sample_df' %in% ls()) { # load sample-batch if not in environment
  batch_sample_df = read.csv(sprintf('%s/combineBroadGEO/broadbatch_geosample_match.csv', data_dir))
}


# prepare cbc names- separate the names into batch/sample number and cbc
# cell_names = jin_meta$NAME
# cell_batch = jin_meta$Batch # should be same
# cell_samplenumber = sapply(X=cell_names, FUN=function(x){strsplit(x=x,split='_')[[1]][2]})
# cell_cbc          = sapply(X=cell_names, FUN=function(x){strsplit(x=x,split='_')[[1]][3]})
# cell_names_df = data.frame(NAME=cell_names, Batch=cell_batch,
#                            samplenumber=cell_samplenumber, barcode=cell_cbc)
# 
# (cell_names_df$Batch == cell_names_df$samplenumber) |> table() # check all same


jin_counts_assembled = NULL
for(cur_batchnumber in 1:18) { 
  # filter for current sample number
  batch_cells = jin_meta |> filter(Batch == cur_batchnumber)
  batch_cbc_formatted = paste0(batch_cells$cbc, '-1') # format differently for h5
  
  # map the Broad Batch number to the GEO sample number
  cur_samplenumber = batch_sample_df |> filter(batch == cur_batchnumber) |> pull(sample)
  assertthat::assert_that(length(cur_samplenumber) == 1, msg = 'For each Broad Batch Number, only should have 1 GEO Sample Number')
  
  
  # find the filename: has 'sampleXX.h5' in the name
  fn = grep(paste0('sample', cur_samplenumber, '.h5'), GEO_files, value = TRUE)
  print(fn)
  assertthat::assert_that(length(fn) == 1, msg = 'GEO h5 file found too many or few files (only should have 1 match)')
  
  # load in raw counts from GEO 
  expression_matrix = Seurat::Read10X_h5( sprintf('%s/%s', GEO_data_folder, fn))
  
  assertthat::assert_that(all(batch_cbc_formatted %in% colnames(expression_matrix)), 
                          msg = 'selected cbc must be in the loaded GEO h5 file')
  
  cur_counts = expression_matrix[, batch_cbc_formatted]
  colnames(cur_counts) = batch_cells$NAME
  
  jin_counts_assembled = cbind(jin_counts_assembled, cur_counts)
  # rm(cur_cell_names_df, cbc_formatted, fn, expression_matrix, counts_cur)
}


saveRDS(jin_counts_assembled, sprintf('%s/combineBroadGEO/gene.rds', data_dir))





# =============== \__Check created gene counts =================
# (1) check that the names match with Broad gene TPM  and (2) that the nUMIs match

# load in all counts as logTPM in txt.gz format (this takes a while)
Broad_gene = read.delim(sprintf("%s/broad/expression/expression_PertCortex.new.txt.gz", data_dir), header = T, stringsAsFactors = F)
row.names(Broad_gene) = Broad_gene$GENE # 1st col is gene name
Broad_gene = Broad_gene |> select(-GENE)# assign as rowname and remove

# load in created counts
jin_counts_assembled = readRDS(sprintf('%s/combineBroadGEO/gene.rds', data_dir))

# (1) check the names match with Broad gene

# dimensions
print(sprintf('%s - CHECK: dims match', all(dim(Broad_gene)==dim(jin_counts_assembled))))

# indices of 0 should match (this can take a little)
zero_idx_match = (which(Broad_gene == 0) == which(as.matrix(jin_counts_assembled) == 0)) 
print(sprintf('%s - CHECK: zero idx match', all(zero_idx_match)))


print(sprintf('%s - CHECK: cell names in Broad are in Combined', all(colnames(Broad_gene) %in% colnames(jin_counts_assembled))))
print(sprintf('%s - CHECK: cell names in Combined are in Broad', all(colnames(jin_counts_assembled) %in% colnames(Broad_gene))))
print(sprintf('%s - CHECK: cell names match (same order)', all(colnames(Broad_gene) == colnames(jin_counts_assembled))))

# (2) that the nUMIs match

# From Broad, load in again 
jin_meta = read.delim(sprintf('%s/broad/metadata/meta_PertCortex.txt', data_dir), header = T, stringsAsFactors = F) 
jin_meta = jin_meta[-1, ] # 1st row (after header) is col type
rownames(jin_meta) <- jin_meta[, 1]

jin_counts_assembled_nUMIs = colSums(jin_counts_assembled)
print(sprintf('%s - CHECK: nUMIs match', all(jin_meta[names(jin_counts_assembled_nUMIs), 'nUMI'] == jin_counts_assembled_nUMIs)))




# =============== END =================



# =============== TRASH =================
# =============== TRASH =================
# =============== TRASH =================
# =============== TRASH =================
# =============== TRASH =================
# =============== TRASH =================



if(F) {
  
  
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
  
  
  
}