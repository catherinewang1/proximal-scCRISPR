# ------------------------------------------------------------------------------------------------ #
#    VCC data subsampled
# ------------------------------------------------------------------------------------------------ #

# =============== Setup directories ========
args = commandArgs(trailingOnly = TRUE)
args = c('macbook')




assertthat::assert_that(length(args) > 0, msg="must give arg for specifying device eg 'Rscript <filename>.R ubergenno'")
DEVICE = args[1]
source('../PATHS.R') # load in data_dir and save_dir and CODE_DIR, depending on DEVICE value
save_dir = paste0(save_dir, '_subsampled') # FOR EXCITATORY NEURONS + USE UNPERTURBED CELLS AS CONTROLS
assertthat::assert_that(!is.null(data_dir), msg='first arg must be: laptop, desktop, or ubergenno')

# library(dplyr)
library(Seurat)
# library(SeuratObject)
library(anndataR)



# load h5ad
h5ad_file  = sprintf('%s/VCC-subsampled.h5ad', data_dir)
vcc_seurat = anndataR::read_h5ad(h5ad_file, as = "Seurat")

gene_continuous = vcc_seurat[['RNA']]['X']
cell_metadata   = vcc_seurat[[]]


gene_continuous[1:2, 1:4]
exp(0.5705448 )
2^0.5705448 
exp(0.5705448/1e6)

test_y = gene_continuous[1, ] |> unique() |> sort()

exp(test_y)
exp(test_y)
exp(test_y - 1)
exp(test_y) - 1


Sys.setenv(UV_OFFLINE=1)
vcc_full = anndata::read_h5ad(sprintf('%s/VCC.h5ad', data_dir))

names(vcc_full)
vcc_full$layers
vcc_full$
vcc_full$shape() # 17k genes, subset has 15k
a = vcc_full$as_Seurat()
vcc_full_seurat = anndataR::read_h5ad(sprintf('%s/VCC.h5ad', data_dir), as = "Seurat")


gene_continuous



cell_metadata |> head()

sum(gene_continuous[, 'train__AAACAAGCATTGCCGCACTTTAGG-Flex_1_01'])

cell_metadata['train__AAACAAGCATTGCCGCACTTTAGG-Flex_1_01', 'nCount_RNA']

# nCount_RNA are non-integers
# is this the sum here or overall...just here for the 15k genes out of 17k total genes measured
colSums(gene_continuous[, 1:1000]) == cell_metadata[1:1000, 'nCount_RNA']
all(colSums(gene_continuous) == cell_metadata[, 'nCount_RNA'])



gene_continuous[, 'train__AAACAAGCATTGCCGCACTTTAGG-Flex_1_01'] |> hist()
exp(gene_continuous[, 'train__AAACAAGCATTGCCGCACTTTAGG-Flex_1_01']) |> hist()


log_expr = gene_continuous[, 'train__AAACAAGCATTGCCGCACTTTAGG-Flex_1_01']

temp <- function(x) {
  if(x==0) {return(0)} 
  else{ return(exp(x))}
}
sapply(X=log_expr, FUN = temp)

sapply(X=log_expr, FUN = temp) |> hist()
sapply(X=log_expr, FUN = temp)



# ==== testing out loading data

# Sys.setenv(UV_OFFLINE=1)
# vcc = anndata::read_h5ad(sprintf('%s/VCC-subsampled.h5ad', data_dir))
# vcc_seurat = vcc$as_Seurat()


h5ad_file = sprintf('%s/VCC-subsampled.h5ad', data_dir)
vcc_seurat = anndataR::read_h5ad(h5ad_file, as = "Seurat")


Seurat::Assays(vcc_seurat)

names(vcc_seurat)
vcc_seurat$cell_type
vcc_seurat$nCount_RNA
vcc_seurat$perturbation
vcc_seurat[['RNA']]

gene_continuous = vcc_seurat[['RNA']]['X']
cell_metadata = vcc_seurat[[]]


gene_continuous[, 1] |> hist()




vcc_seurat$dataset
vcc_seurat$cell_type
vcc_seurat$is_control
vcc_seurat$perturbation
vcc_seurat[['perturbation']]


cell_metadata = data.frame(temp = rep(NA, times = ncol(gene_continuous)))
row.names(cell_metadata) = colnames(gene_continuous)

for(vcc_name in names(vcc_seurat[[]])) {
  cur_vcc_data = vcc_seurat[[vcc_name]]
  assertthat::assert_that(all(row.names(cell_metadata) == row.names(cur_vcc_data)), msg = sprintf('%s must have same cellnames and cellname order', vcc_name))
  cell_metadata = cbind(cell_metadata, cur_vcc_data)
}





Seurat::ReadH5AD

dim(vcc)

vcc$target_gene


names(vcc)
vcc[['layers']]


gene_continuous = vcc$X # 67179 cells x 15565 genes
dim(gene_continuous) # 67179 15565
vcc$n_obs # 67179
vcc$n_vars # 15565
temp = vcc$obs_names
temp

row.names(gene_continuous) = as.vector(vcc$obs_names)




expm1 <- function(x) {exp(x) - 1}

test_transform = expm1(gene_continuous[1, ]) / 1e6 * cell_metadata[1, 'nCount_RNA']
hist(test_transform)



plot(gene_continuous[1, ], test_transform)


test_transform = expm1(gene_continuous[1, ]) / 1e6 * 1000

plot(gene_continuous[1, ], expm1(gene_continuous[1, ]) / 1e6 * 1e6 * 1.3)




exp(cell_metadata[1, 'nCount_RNA'])


continuous_unique = gene_continuous[1, ] |> round(3) |> unique() |> sort()
plot(continuous_unique)



plot(gene_continuous[1, ], exp(gene_continuous[1, ]))


head(expm1(gene_continuous[1, ])/1e6)


head(2^gene_continuous[1, ])

2^(gene_continuous[1, ])-1
log(tpm/10+1, 2)








