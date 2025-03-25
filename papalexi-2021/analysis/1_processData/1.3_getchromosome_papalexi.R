# ---------------------------------------------------------------------------- #
#                   Do CB using genes as NCE/NCO                               #
# 1.3  - add chromosome number information                                     #
# Requires: prev saved normalized gene expression (HDF5)                       #
# Ouputs: (nothing) but saves                                                  #
#         gene_chr in the csv file                                             #
#                         "<save_dir>/chromosome/gene_chromosome.csv"          #
#         grna_chr in the csv file                                             #
#                         "<save_dir>/chromosome/grna_chromosome.csv"          #
# ---------------------------------------------------------------------------- #
args = commandArgs(trailingOnly = TRUE)

require(assertthat) # for some assert statements
require(ondisc)
require(tibble)
library(dplyr)
library(biomaRt)    # for pulling gene-chr info

assertthat::assert_that(length(args) > 0, msg="must give arg for specifying device eg 'Rscript <filename>.R ubergenno'")

DEVICE = args[1]
source('../PATHS.R') # load in data_dir and save_dir and CODE_DIR, depending on DEVICE value


# =================== Start ====================================================
print(sprintf("[%s] START: CB Chromosome", Sys.time()))

# make dir for saving
dir.create(sprintf('%s/chromosome', save_dir), showWarnings = FALSE)

# attributes to pull
myAttributes = c("wikigene_name", "wikigene_id", "chromosome_name", 
                 "ensembl_gene_id", "external_gene_name")

# gene names to pull (using wikigene_id)
# #      to get all genes
# gene_odm <- ondisc::read_odm(odm_fp      = paste0(data_dir, "/papalexi-2021/processed/gene/expression_matrix.odm"),
#                              metadata_fp = paste0(data_dir, "/papalexi-2021/processed/gene/metadata.rds"))
# 
# myGenenames = gene_odm |> ondisc::get_feature_covariates() |> row.names() 

# OR
#      to get top genes
myGenenames     = readRDS(sprintf('%s/important_genes_name.rds', save_dir))
myGenenames_idx = readRDS(sprintf('%s/important_genes_idx.rds',  save_dir))


# grna target gene names to pull
grna_odm <- ondisc::read_odm(odm_fp      = paste0(data_dir, "/papalexi-2021/processed/grna_assignment/assignment_matrix.odm"),
                             metadata_fp = paste0(data_dir, "/papalexi-2021/processed/grna_assignment/metadata.rds"))
myGrnaGenenames = grna_odm |> ondisc::get_feature_covariates() |> 
                  filter(target != 'non-targeting') |> 
                  pull(target) |> unique()

# make a mart object
ensembl = useMart("ensembl")
ensembl = useDataset("hsapiens_gene_ensembl", mart=ensembl)


# =================== get GENE chromosome info ======================================
print(sprintf("[%s]    - get GENE chromosome info from biomaRt", Sys.time()))

g = getBM(attributes = myAttributes,
          filters = "wikigene_name",
          values = myGenenames,
          mart = ensembl)


# =================== format GENE chr df and save ===================================
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


# =================== get GRNA chromosome info ======================================
print(sprintf("[%s]    - get GRNA chromosome info from biomaRt", Sys.time()))

g_grna = getBM(attributes = myAttributes, 
               filters = "wikigene_name", 
               values = myGrnaGenenames, 
               mart = ensembl)


# =================== format GRNA chr df and save ===================================
print(sprintf("[%s]        - format GRNA chr df and save", Sys.time()))

# # ???? shows up on website but not using biomaRt (ignore this gene for now)
# getBM(attributes = myAttributes, 
#       filters = 'hgnc_symbol', #"wikigene_name", 
#       values = 'MARCH8', 
#       mart = ensembl)
grna_df = merge(grna_odm |> 
                  ondisc::get_feature_covariates() |> 
                  tibble::rownames_to_column("grna"),
                g_grna |> 
                  filter(chromosome_name %in% c(as.character(1:100), 'A', 'Y')) |>
                  mutate(target_wikigene_name = wikigene_name,
                         target_chromosome_name = chromosome_name) |> 
                  dplyr::select(target_wikigene_name, target_chromosome_name),
                by.x = 'target', by.y = 'target_wikigene_name',
                all.x = TRUE) |> dplyr::select(grna, target, target_chromosome_name, 
                                        n_nonzero, known_protein_effect)
write.csv(grna_df, file = sprintf('%s/chromosome/grna_chromosome.csv', save_dir), row.names = F)


# =================== END ======================================================
print(sprintf("[%s] END", Sys.time()))
