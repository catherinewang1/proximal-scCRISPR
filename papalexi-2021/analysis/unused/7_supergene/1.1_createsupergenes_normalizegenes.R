# ---------------------------------------------------------------------------- #
#                   desc                                                       #
# desc more                                                                    #
# ---------------------------------------------------------------------------- #
# setup script parameters (e.g. paths)
args = commandArgs(trailingOnly = TRUE)
# args = c('laptop', 100)



suppressPackageStartupMessages(library(crayon))
cat(crayon::blue('\n')) # for colors, sets some setting in terminal to display ANSI colors correctly





assertthat::assert_that(length(args) > 0, msg="must give arg for specifying device eg 'Rscript <filename>.R ubergenno'")
DEVICE = args[1]
source('../PATHS.R') # load in data_dir and save_dir and CODE_DIR, depending on DEVICE value

assertthat::assert_that(length(args) > 1, msg="must give arg for specifying number of important non-Transcription Factor genes'")
NUM_IMPORTANT_GENES = as.integer(args[2])



myPrintColor <- function(txt, col=36) {
  # change number from 29:47 # for(col in 29:47){ cat(paste0("\033[0;", col, "m",txt,"\033[0m","\n"))}
  cat(paste0("\033[0;", col, "m",txt,"\033[0m","\n"))
}

# folder to save results in
dir.create(sprintf('%s/supergene/', save_dir), showWarnings = FALSE)

cat(' \n ')

# ==================================================================================================
# =================== Start ========================================================================
# ==================================================================================================
# print(sprintf("[%s] START: 7_supergenes/1.0_createsupergenes.R", Sys.time()))
myPrintColor(sprintf("[%s] START: 7_supergenes/1.0_createsupergenes.R", Sys.time()) |> as.character())


# load in list of all genes measured from papalexi-2021
gene_odm <- ondisc::read_odm(odm_fp      = paste0(data_dir, "/papalexi-2021/processed/gene/expression_matrix.odm"),
                             metadata_fp = paste0(data_dir, "/papalexi-2021/processed/gene/metadata.rds"))
grna_odm <- ondisc::read_odm(odm_fp      = paste0(data_dir, "/papalexi-2021/processed/grna_assignment/assignment_matrix.odm"),
                             metadata_fp = paste0(data_dir, "/papalexi-2021/processed/grna_assignment/metadata.rds"))

# load all into memory
gene <- gene_odm[[,1:ncol(gene_odm)]] # entire gene ds (18649 genes x 20729 cells)
# grna <- grna_odm[[,1:ncol(grna_odm)]] # entire grna ds (110   grnas x 20729 cells)
invisible(gc(verbose=FALSE))

gene_names_original = gene_odm |> ondisc::get_feature_ids()


# ==================================================================================================
# =================== Choose Top XX nonTF genes ====================================================
# ==================================================================================================
# print(sprintf("[%s]    - Choose Top %d nonTF genes", Sys.time(), NUM_IMPORTANT_GENES))
myPrintColor(sprintf("[%s]    - Choose Top %d nonTF genes", Sys.time(), NUM_IMPORTANT_GENES))


# =================== load list of genes that are Transcription Factors ===========================
# print(sprintf("[%s]          loading list of TF genes", Sys.time()))
myPrintColor(sprintf("[%s]          loading list of TF genes", Sys.time()))
# read in xlsx sheet
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



# =================== get gene importance (deviance) ==============================================
# Load in top genes in order of importance (output from 1_processData 1.1_findimportantgenes_papalex.R)
# print(sprintf("[%s]          loading gene deviances/importances", Sys.time()))
myPrintColor(sprintf("[%s]          loading gene deviances/importances", Sys.time()))

gene_deviance = read.csv(sprintf('%s/gene_deviance.csv', save_dir))


nonTF_names = gene_deviance |> 
              dplyr::filter(!gene_name %in% TF_names) |> 
              dplyr::arrange(desc(deviance))
# Check removed the right amount of genes
#                        number of genes
# dim(gene_deviance)[1] # before: 18649
# length(TF_names)      # TF: 1639
# dim(gene_deviance)[1] - length(TF_names) # difference: 17010
# dim(nonTF_names)[1]   # remaining: 17370
#
# hmm there are some TF genes listed that aren't measured, is that okay?
# which(!TF_names %in% gene_deviance$gene_name) # |> length()# names from importance values
# which(!TF_names %in% gene_names_original) # |> length()# names from original 
# hmm there are some names in this sheet that were not measured in papalexi-2021



# =================== take top XX non-TF gene names ==============================================
# take top XX non-TF gene names (#top specified as script parameter)
# print(sprintf("[%s]          loading list of TF genes", Sys.time()))
myPrintColor(sprintf("[%s]          loading list of TF genes", Sys.time()))
topnonTF_gene_names = nonTF_names[1:NUM_IMPORTANT_GENES, ] |> dplyr::pull(gene_name)

# Save top NUM_IMPORTANT_GENES with their deviances (in order of most to least imp)
write.csv(nonTF_names[1:NUM_IMPORTANT_GENES, ], 
         sprintf('%s/supergene/important_genes_nonTF.csv', save_dir),
         row.names = FALSE)


# perform sparse PCA, take top YY loadings, take genes with non-zero coefs as members of each PC loading
# perform PCA on top XX non-TF genes + expressions (expensive)
# ==================================================================================================
# =================== Perform Normalization! Expensive! ============================================
# ==================================================================================================
# code changed from 1_processData/1.2_normalize_papalexi.R
# print(sprintf("[%s]    - Normalize papalexi (top %d non TF genes)", Sys.time(), NUM_IMPORTANT_GENES))
myPrintColor(sprintf("[%s]    - Normalize papalexi (top %d non TF genes)", Sys.time(), NUM_IMPORTANT_GENES))

# =================== Convert gene exp to HDF5 file by converting in memory ========================
# print(sprintf("[%s]          loading list of TF genes", Sys.time()))
myPrintColor(sprintf("[%s]          loading list of TF genes", Sys.time()))

# topnonTF_gene_names # instead of  important_genes_idx   = readRDS(sprintf('%s/important_genes_idx.rds', save_dir))

invisible(gc(verbose=FALSE))

# print(sprintf("[%s]        - converting to HDF5 file", Sys.time()))
myPrintColor(sprintf("[%s]        - converting to HDF5 file", Sys.time()))
# will be saved in...
h5file   = paste0(save_dir, "/gene_nonTF.h5")
HDF5Array::setHDF5DumpFile(h5file)

# delete existing 'gene_nonTF' if already exists
h5f         = rhdf5::H5Fopen(h5file);
name_exists = rhdf5::H5Lexists(h5f, 'gene_nonTF')
rhdf5::h5closeAll()
if (name_exists) {  rhdf5::h5delete(file = h5file, name = 'gene_nonTF') }
rm(h5f, name_exists)

# save gene expression of top genes in 'gene.h5' under 'gene'
HDF5Array::setHDF5DumpName('gene_nonTF')

gene_important = gene_odm[[topnonTF_gene_names,
                           1:ncol(gene_odm)]] # load top gene expr into memory (NUM_IMPORTANT_GENES x #cells)
gene_hd5 = as(gene_important, "HDF5Matrix")
rhdf5::h5closeAll()
invisible(gc(verbose=FALSE))


# =================== Perform normalization ========================================================
# print(sprintf("[%s]        - performing normalization", Sys.time()))
myPrintColor(sprintf("[%s]        - performing normalization", Sys.time()))
gene_norm = scry::nullResiduals(gene_hd5, 
                                fam  = "binomial", 
                                type = "deviance") # #imp genes x 20729


# =================== save normalized gene expression ==============================================
# print(sprintf("[%s]        - saving normalized gene expression", Sys.time()))
myPrintColor(sprintf("[%s]        - saving normalized gene expression", Sys.time()))
# delete existing gene_norm_nonTF if already exists
h5f         = rhdf5::H5Fopen(h5file); 
name_exists = rhdf5::H5Lexists(h5f, 'gene_norm_nonTF')
rhdf5::h5closeAll()
if (name_exists) {  rhdf5::h5delete(file = h5file, name = 'gene_norm_nonTF') }
rm(h5f, name_exists)

# save normalized gene expression of top genes in 'gene.h5' under 'gene_norm_nonTF'
if(DEVICE == 'ubergenno') { #if on ubergenno, need to realize to memory first?????
  HDF5Array::writeHDF5Array(DelayedArray::realize(gene_norm), filepath=h5file, name='gene_norm_nonTF')
} else {
  HDF5Array::writeHDF5Array(gene_norm, filepath=h5file, name='gene_norm_nonTF')
}

invisible(gc(verbose=FALSE))

# 



# ==================================================================================================
# =================== PCA (normalized genes) =======================================================
# ==================================================================================================











# =================== END ==========================================================================
print(sprintf("[%s] END", Sys.time()))






