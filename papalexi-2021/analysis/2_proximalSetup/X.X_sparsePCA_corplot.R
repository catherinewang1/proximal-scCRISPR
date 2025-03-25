# Make a correlation plot of the sparse PCA result




library(dplyr)
library(ggplot2)
library(RColorBrewer)




# === Parameter Settings from normalizing genes 
NUM_IMPORTANTGENES = 4000
# === Parameter Settings from SPCA
my_sumabsv = 5
my_K = 60
N_subsample = 'all' # subsample size, or 'all' if using all cells
# === Parameter Settings for device (which device is being used)
DEVICE = 'laptop'


# first set working dir to current file save location (should be in some sort of papalexi-2021/analysis/)
source('../PATHS.R') # load in data_dir and save_dir and CODE_DIR, depending on DEVICE value

dir.create(sprintf('%s/spca/spcaplots/', save_dir))



# create gene and grna ondisc managers
gene_odm <- ondisc::read_odm(odm_fp      = paste0(data_dir, "/papalexi-2021/processed/gene/expression_matrix.odm"),
                             metadata_fp = paste0(data_dir, "/papalexi-2021/processed/gene/metadata.rds"))
grna_odm <- ondisc::read_odm(odm_fp      = paste0(data_dir, "/papalexi-2021/processed/grna_assignment/assignment_matrix.odm"),
                             metadata_fp = paste0(data_dir, "/papalexi-2021/processed/grna_assignment/metadata.rds"))

# load grna assignments (load all into memory)
grna = grna_odm[[,1:ncol(grna_odm)]] # |> as.matrix() # ~110 x 20729 = #grnas x #cells
grna_rownames = grna_odm |> ondisc::get_feature_covariates() |> rownames()

# load normalized gene exp
h5file      = paste0(save_dir, "/gene.h5"); print(h5file)
reading_hd5file  = rhdf5::H5Fopen(name = h5file)
readin_gene_norm = reading_hd5file&'gene_norm'
gene_norm = readin_gene_norm[1:NUM_IMPORTANTGENES, ] # eg dim = 4000 x 20729 = #important x #cells
rhdf5::h5closeAll()
invisible(gc(verbose=FALSE))


# use Sparse PCA Loadings
if(N_subsample == 'all') {  N_subsample = ncol(gene_odm) }
NCs     = readRDS(sprintf('%s/spca/NCloadings_sumabs=%.1f_K=%d_N=%d.rds', save_dir, my_sumabsv, my_K, N_subsample))  
outorth = readRDS(sprintf('%s/spca/outorth_sumabs=%.1f_K=%d_N=%d.rds', save_dir, my_sumabsv, my_K, N_subsample))  



# Overall: PCi vs PCj
NCcor = cor(NCs)
NCcor_melted = reshape2::melt(NCcor)
p1 = ggplot(NCcor_melted) +
  geom_tile(aes(x = Var1, y = Var2, fill = value)) +
  scale_fill_gradient2(limits = c(-1, 1), 
                       midpoint = 0, 
                       high = RColorBrewer::brewer.pal(11, 'RdBu')[ 2], 
                       low  = RColorBrewer::brewer.pal(11, 'RdBu')[10]) +
  # scale_fill_brewer(type = 'div',  limits = c(-1, 1)) +
  # scale_fill_distiller(palette = 'RdBu', type = 'div', limits = c(-1, 1)) +
  scale_x_continuous(expand = c(0, 0), breaks = c(1, 15, 30, 45, 60)) +
  scale_y_reverse(expand = c(0, 0), breaks = c(1, 15, 30, 45, 60)) +
  labs(x = 'Sparse PC', y = 'Sparse PC', fill = 'correlation',
       title = 'Correlation between Top 60 PCs') +
  theme_classic() +
  theme(axis.line = element_blank(),
        axis.ticks = element_blank(),
        legend.position = 'right',
        plot.title = element_text(hjust = .5))
p1

ggsave(sprintf('%s/spca/spcaplots/SPC_cor.pdf', save_dir), height = 5, width = 6)

# limit to top 20 PCs
topXPCs = 20
NCcor = cor(NCs[,1:20])
NCcor_melted = reshape2::melt(NCcor)
p1 = ggplot(NCcor_melted) +
  geom_tile(aes(x = Var1, y = Var2, fill = value)) +
  scale_fill_gradient2(limits = c(-1, 1), 
                       midpoint = 0, 
                       high = RColorBrewer::brewer.pal(11, 'RdBu')[ 2], 
                       low  = RColorBrewer::brewer.pal(11, 'RdBu')[10]) +
  # scale_fill_brewer(type = 'div',  limits = c(-1, 1)) +
  # scale_fill_distiller(palette = 'RdBu', type = 'div', limits = c(-1, 1)) +
  scale_x_continuous(expand = c(0, 0), 
                     # breaks = c(1, 15, 30, 45, 60),
                     breaks = c(1, 5, 10, 15, 20)) +
  scale_y_reverse(expand = c(0, 0), 
                  # breaks = c(1, 15, 30, 45, 60),
                  breaks = c(1, 5, 10, 15, 20)) +
  labs(x = 'Sparse PC', y = 'Sparse PC', fill = 'correlation',
       title = sprintf('Correlation between Top %d PCs', topXPCs)) +
  theme_classic() +
  theme(axis.line = element_blank(),
        axis.ticks = element_blank(),
        legend.position = 'right',
        plot.title = element_text(hjust = .5))
p1

ggsave(sprintf('%s/spca/spcaplots/SPC_cor_%d.pdf', save_dir, topXPCs), height = 5, width = 6)


# Individual Genes of top 10 PCs:


# setup the data
myGenenames_df = read.csv(sprintf('%s/important_genes_name_TF_Target.csv', save_dir))
myGenenames_df_noTFTargets = myGenenames_df |> dplyr::filter((!TF) & (!grna_target))
gene_norm_noTFTargets = gene_norm[myGenenames_df_noTFTargets |> dplyr::pull(importance_rank), ]

# Get the top 10 genes of top 10 PCs
PC10 = outorth$v[,1:10] |> data.frame()
colnames(PC10) = paste0('PC', 1:10)
PC10$PCgeneidx = 1:nrow(PC10)

# save which genes are top in which PCs
PCgenemap = NULL
# save the genes' expressions 
gene_norm_PC10 = NULL
gene_norm_PC10_coef = NULL # gene exp mult by SPC coef

# do for top 10 PCs
for(PCnumber in 1:10) {
  # PCnumber = 1
  # find the genes w the largest abs coef
  PCgeneidxs = PC10 |> 
    arrange(desc(abs(eval(parse(text = paste0('PC', PCnumber)))))) |>
    slice(1:50) |>
    pull(PCgeneidx)
  
  PCgenecoefs = PC10[PCgeneidxs, paste0('PC', PCnumber)]
  # get those genes' expression levels
  PCgeneidxs_gene = gene_norm_noTFTargets[PCgeneidxs, ] |> t()
  # colnames(PCgeneidxs_gene) = paste0('PC', PCnumber, '_', 'gene', 1:ncol(PCgeneidxs_gene))
  gene_norm_PC10      = cbind(gene_norm_PC10,      PCgeneidxs_gene) # append to overall df
  gene_norm_PC10_coef = cbind(gene_norm_PC10_coef, PCgeneidxs_gene %*% diag(PCgenecoefs))

  PCgenemap = rbind(PCgenemap,
                    data.frame(PC              = PCnumber, 
                               PCgeneidx       = PCgeneidxs, 
                               importance_rank = myGenenames_df_noTFTargets[PCgeneidxs, 'importance_rank'])   )
}

# Plot: Cor between genes
PC10_cor = cor(gene_norm_PC10)
PC10_cor_melted = reshape2::melt(PC10_cor)

p2 = ggplot(PC10_cor_melted) +
  geom_tile(aes(x = Var1, y = Var2, fill = value)) +
  scale_fill_gradient2(limits = c(-1, 1), 
                       midpoint = 0, 
                       high = RColorBrewer::brewer.pal(11, 'RdBu')[ 2], 
                       low  = RColorBrewer::brewer.pal(11, 'RdBu')[10]) +
  # scale_fill_brewer(type = 'div',  limits = c(-1, 1)) +
  # scale_fill_distiller(palette = 'RdBu', type = 'div', limits = c(-1, 1)) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_reverse(expand = c(0, 0)) +
  labs(x = 'gene', y = 'gene', fill = 'correlation',
       title = 'Correlation between Top PCs Top Genes') +
  theme_classic() +
  theme(axis.line = element_blank(),
        axis.ticks = element_blank(),
        legend.position = 'right',
        plot.title = element_text(hjust = .5))
p2
ggsave(sprintf('%s/spca/spcaplots/SPC_gene_cor.pdf', save_dir), height = 5, width = 6)

# Plot: Cor between genes mult by SPC coefs
PC10_cor_coef = cor(gene_norm_PC10_coef)
PC10_cor_coef_melted = reshape2::melt(PC10_cor_coef)

p3 = ggplot(PC10_cor_coef_melted) +
  geom_tile(aes(x = Var1, y = Var2, fill = value)) +
  scale_fill_gradient2(limits = c(-1, 1), 
                       midpoint = 0, 
                       high = RColorBrewer::brewer.pal(11, 'RdBu')[ 2], 
                       low  = RColorBrewer::brewer.pal(11, 'RdBu')[10]) +
  # scale_fill_brewer(type = 'div',  limits = c(-1, 1)) +
  # scale_fill_distiller(palette = 'RdBu', type = 'div', limits = c(-1, 1)) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_reverse(expand = c(0, 0)) +
  labs(x = 'gene', y = 'gene', fill = 'correlation',
       title = 'Correlation between Top PCs Top Genes*PC Coef') +
  theme_classic() +
  theme(axis.line = element_blank(),
        axis.ticks = element_blank(),
        legend.position = 'right',
        plot.title = element_text(hjust = .5))
p3
ggsave(sprintf('%s/spca/spcaplots/SPC_genecoef_cor.pdf', save_dir), height = 5, width = 6)




pdf(sprintf('%s/spca/spcaplots/SPC_allcors.pdf', save_dir), height = 5, width = 18)
gridExtra::grid.arrange(p1, p2, p3, nrow = 1)
dev.off()



