# ---------------------------------------------------------------------------- #
#                   Perform PCA then PCA/UMAP on papalexi dataset              #                            #
# Requires: prev saved normalized gene expression (HDF5)                       #
# Ouputs: (nothing) but saves dim red and plots at                             #
#         '<save_dir>/pca/'                                                    #
#         '<save_dir>/pca_umap/'                                               #
# ---------------------------------------------------------------------------- #
# args = commandArgs(trailingOnly = TRUE)
args = c('laptop', '2000')

require(assertthat) # for some assert statements
require(rhdf5)      # read/write HDF5 format
library(umap)       # perform umap
library(dplyr)
library(ggplot2)
library(cowplot)
library(viridis)
theme_set(theme_cowplot() +
            theme(plot.title = element_text(hjust = .5),
                  plot.subtitle = element_text(hjust = .5),
                  strip.background = element_rect(fill = "white", color = 'black'),
                  panel.grid.minor = element_blank(),
                  panel.grid.major = element_blank()))

assertthat::assert_that(length(args) > 0, msg="must give arg for specifying device eg 'Rscript <filename>.R ubergenno'")
DEVICE = args[1]
source('../PATHS.R') # load in data_dir and save_dir, depending on DEVICE value

assertthat::assert_that(length(args) >= 2, msg="provide device and NUM_IMPORTANT_GENES (currently int <=4000)")
NUM_IMPORTANT_GENES = as.integer(args[2])


PCA_RANK = 50

# =================== START ====================================================
print(sprintf("[%s] START: PCA and PCA/UMAP papalexi (top %d genes)", Sys.time(), NUM_IMPORTANT_GENES))
h5file      = paste0(save_dir, "/gene.h5")

# ==================================================================================================
# =================== loading gene_norm ========================================
print(sprintf("[%s]    - loading normalized gene exp", Sys.time()))
reading_hd5file  = rhdf5::H5Fopen(name = h5file)
readin_gene_norm = reading_hd5file&'gene_norm'
gene_norm = readin_gene_norm[1:NUM_IMPORTANT_GENES, ] # dim = 4000 x 20729 = #important x #cells
rhdf5::h5closeAll()
invisible(gc(verbose=FALSE))


# =================== scaling gene_norm ========================================
print(sprintf("[%s]    - scaling normalized gene exp", Sys.time()))
gene_norm_scaled = apply(gene_norm, MARGIN = 1, FUN = scale)
invisible(gc(verbose=FALSE))

# filling NAs in gene_norm_scaled
# anyNA(gene_norm_scaled) # but after feat sel, shouldn't have any?? some still 0sd
# gene_norm_scaled[is.na(gene_norm_scaled)] = 0 
# invisible(gc(verbose=FALSE))



# =================== performing PCA ===========================================
print(sprintf("[%s]    - performing PCA (top %s PCs)", Sys.time(), PCA_RANK))
t0 = Sys.time()
pca_fit = gene_norm_scaled[,] |> prcomp(rank. = PCA_RANK)
t1 = Sys.time(); print(sprintf('[%s]         %.2f mins', Sys.time(), difftime(t1, t0, units ="mins")))
invisible(gc(verbose=FALSE))


# =================== saving PCA result ========================================
print(sprintf("[%s]    - saving PCA", Sys.time()))
dir.create(sprintf("%s/eda/objects", save_dir), recursive=TRUE, showWarnings=FALSE)
saveRDS(pca_fit,
        sprintf("%s/eda/objects/pca_NUMGENES=%s_PCARANK=%s.rds", save_dir, NUM_IMPORTANT_GENES, PCA_RANK))


# =================== saving PCA plot ==========================================
print(sprintf("[%s]    - saving PCA plot", Sys.time()))

p = ggplot(NULL,
           aes(x = pca_fit$x[,'PC1'], y = pca_fit$x[,'PC2'])) +
  geom_point(alpha=.05, size = .5) +
  labs(x = 'PC1', y = 'PC2', 
       title = sprintf('PCA #Genes=%s #PCs=%s', NUM_IMPORTANT_GENES, PCA_RANK)) +
  theme(plot.title = element_text(hjust = .5))

ggsave(plot=p, 
       filename=sprintf('%s/eda/pca_NUMGENES=%s_PCARANK=%s.png', 
                        save_dir, NUM_IMPORTANT_GENES, PCA_RANK),
       height = 5, width = 6, dpi = 320)

ggsave(plot=p, 
       filename=sprintf('%s/eda/pca_NUMGENES=%s_PCARANK=%s.pdf', 
                        save_dir, NUM_IMPORTANT_GENES, PCA_RANK),
       height = 5, width = 6)


# =================== performing UMAP ==========================================
print(sprintf("[%s]    - performing UMAP", Sys.time()))
t0 = Sys.time()
umap_fit = pca_fit$x |> umap(n_components = 2)
t1 = Sys.time(); print(sprintf('[%s]         %.2f mins', Sys.time(), difftime(t1, t0, units ="mins")))
invisible(gc(verbose=FALSE))


# =================== saving UMAP result =======================================
print(sprintf("[%s]    - saving UMAP", Sys.time()))
umap_fit$knn  = 'removed for space'

# dir.create(sprintf("%s/eda/objects/", save_dir), recursive=TRUE, showWarnings=FALSE)
saveRDS(umap_fit,
        sprintf("%s/eda/objects/pca_umap_NUMGENES=%s_PCARANK=%s.rds", save_dir, NUM_IMPORTANT_GENES, PCA_RANK))


# =================== saving UMAP plot =========================================
print(sprintf("[%s]    - saving UMAP plot", Sys.time()))

p = ggplot(NULL,
           aes(x = umap_fit$layout[,1], y = umap_fit$layout[,2])) +
  geom_point(alpha=.05, size = .5) +
  labs(x = 'pca-umap1', y = 'pca-umap2', 
       title = sprintf('PCA-UMAP #Genes=%s #PCs=%s', NUM_IMPORTANT_GENES, PCA_RANK)) +
  theme(plot.title = element_text(hjust = .5))

ggsave(plot=p, 
       filename=sprintf('%s/eda/pca_umap_NUMGENES=%s_PCARANK=%s.png', 
                        save_dir, NUM_IMPORTANT_GENES, PCA_RANK),
       height = 5, width = 6, dpi = 320)

ggsave(plot=p, 
       filename=sprintf('%s/eda/pca_umap_NUMGENES=%s_PCARANK=%s.pdf', 
                        save_dir, NUM_IMPORTANT_GENES, PCA_RANK),
       height = 5, width = 6)





# ==================================================================================================
# =================== Plots by Variable (e.g. confounders or treatment) ============================
print(sprintf("[%s]    - Plots by Variable", Sys.time()))


# loads in prev saved pca and umap results (so you can start the script here)
source(sprintf('%s/dimred_plots_utils.R', util_dir))


papalexi_dir = paste0(data_dir, '/papalexi-2021')
gene_odm <- ondisc::read_odm(odm_fp      = paste0(papalexi_dir, "/processed/gene/expression_matrix.odm"),
                     metadata_fp = paste0(papalexi_dir, "/processed/gene/metadata.rds"))
grna_odm <- ondisc::read_odm(odm_fp      = paste0(papalexi_dir, "/processed/grna_assignment/assignment_matrix.odm"),
                     metadata_fp = paste0(papalexi_dir, "/processed/grna_assignment/metadata.rds"))

chosen_genes = c('FTH1', 'EEF1A1', 'FAU') 
chosen_grna_T  = c('IFNGR2g1', 'ATF2g1', 'CD86g1')
all_grna_NT = grna_odm |> ondisc::get_feature_covariates() |> filter(target_type == 'non-targeting') |> rownames()





# =================== Constructing Dataframes of Dim Reds ========================
print(sprintf("[%s]    - Constructing dataframes of dim reds", Sys.time()))



papalexi_pca      = readRDS(sprintf(     '%s/eda/objects/pca_NUMGENES=%s_PCARANK=%s.rds', save_dir, NUM_IMPORTANT_GENES, PCA_RANK))
papalexi_pca_umap = readRDS(sprintf('%s/eda/objects/pca_umap_NUMGENES=%s_PCARANK=%s.rds', save_dir, NUM_IMPORTANT_GENES, PCA_RANK))

# dataframe indicating gRNA each cell received
grna_assignment = get_grna_assignment(grna_odm = grna_odm)

# function for getting dim red for same (non-dimred method) params 
my_get_dimred = get_dimred_bydimredmethod(grna_assignment = grna_assignment, 
                                          cell_covariates = gene_odm |> get_cell_covariates(),
                                          gRNA_names = chosen_grna_T,
                                          NT_names = all_grna_NT,
                                          NT_exposed_size=NULL,
                                          seed = 12345) 

pca_cells      = my_get_dimred(df_dimred = as.data.frame(papalexi_pca$x[,1:2]))
# umap_cells     = my_get_dimred(df_dimred = as.data.frame(papalexi_umap$layout)) # no? only umap after pca
pca_umap_cells = my_get_dimred(df_dimred = as.data.frame(papalexi_pca_umap$layout))


# =================== Plotting Dim Reductions ========================
print(sprintf("[%s]    - Plotting dim red (all at once)", Sys.time()))
# function for plotting dim red for same (non-dimred method) params

my_plot_dimred = plot_dimred_bydimredmethod(save_dir=sprintf('%s/eda/', save_dir), 
                                            num_imp_genes = NUM_IMPORTANT_GENES, 
                                            gRNA_names = chosen_grna_T,
                                            NT_names = all_grna_NT)

my_plot_dimred(dimred_type='pca',      gRNAcells=pca_cells$gRNAcells,      allcells=pca_cells$allcells)
# my_plot_dimred(dimred_type='umap',     gRNAcells=umap_cells$gRNAcells,     allcells=umap_cells$allcells)
my_plot_dimred(dimred_type='pca_umap', gRNAcells=pca_umap_cells$gRNAcells, allcells=pca_umap_cells$allcells)


# =================== Plotting Dim Reductions ========================
print(sprintf("[%s]    - Plotting dim red (choose some to customize by hand)", Sys.time()))


dir.create(sprintf('%s/eda/custom/', save_dir))
dir.create(sprintf('%s/eda/custom/%s_NUMGENES=%s_PCARANK=%s/', 
                   save_dir, 'pca_umap', NUM_IMPORTANT_GENES, PCA_RANK))




# gRNA i = IFNGR2g1 vs NT
dimred_type = 'pca_umap'
gRNA_i = 'IFNGR2g1'
gRNAcells = pca_umap_cells$gRNAcells
x_label = 'UMAP 1'
y_label = 'UMAP 2'

ggplot(gRNAcells |> filter(gRNA %in% c(gRNA_i, all_grna_NT)), 
       aes(x = dim1, y = dim2, color = gRNA_label)) +
  geom_point(alpha = .6, size = 1.5, shape = 19, stroke=NA) +
  labs(title = paste0(gRNA_i, ' vs Non-Targeting'), color = 'gRNA',
       x = x_label, y = y_label) +
  scale_color_brewer(palette = 'Dark2') +
  guides(color = guide_legend(override.aes=list(shape = 15, size = 3.5, alpha = 1))) +
  theme(plot.title = element_text(hjust = .5),
        legend.position = c(.73, .86),
        legend.box.background = element_rect(colour = "black", linejoin = "bevel", linewidth = .5), 
        legend.margin = margin(t = 2, r = 4, b = 2, l = 4, unit = "pt"), 
        legend.key.size = unit(.45, 'cm'),
        panel.background = element_rect(color = 'black', linewidth = 1), 
        axis.line = element_blank())

ggsave(filename=sprintf('%s/eda/custom/%s_NUMGENES=%s_PCARANK=%s/gRNA_%s.png', 
                        save_dir, dimred_type, NUM_IMPORTANT_GENES, PCA_RANK, gRNA_i),
       height = 5, width = 6, dpi = 320)

ggsave(filename=sprintf('%s/eda/custom/%s_NUMGENES=%s_PCARANK=%s/gRNA_%s.pdf', 
                        save_dir, dimred_type, NUM_IMPORTANT_GENES, PCA_RANK, gRNA_i),
       height = 5, width = 6)





# Plot all chosen gRNA vs Non-Targeting gRNA
ggplot() +
  geom_point(data = gRNAcells |> filter(gRNA_label != 'NT'),
             aes(x = dim1, y = dim2, color = gRNA_label),
             alpha = .7, size = 1.5, shape = 19, stroke=NA) +
  geom_point(data = gRNAcells |> filter(gRNA_label == 'NT'),  # NT cells more transparent and smaller
             aes(x = dim1, y = dim2, color = gRNA_label),
             alpha = .3, size = 1.25, shape = 19, stroke=NA) +
  labs(title = paste0('Some gRNA vs Non-Targeting gRNA'), color = 'gRNA',
       x = x_label, y = y_label) +
  scale_color_brewer(palette = 'Dark2') +
  guides(color = guide_legend(override.aes=list(shape = 15, size = 3.5, alpha = c(1, 1, 1, .5)))) +
  theme(plot.title = element_text(hjust = .5),
        legend.position = c(.72, .81),
        legend.box.background = element_rect(colour = "black", linejoin = "bevel", linewidth = .5), 
        legend.margin = margin(t = 2, r = 4, b = 2, l = 4, unit = "pt"), 
        legend.key.size = unit(.45, 'cm'),
        panel.background = element_rect(color = 'black', linewidth = 1), 
        axis.line = element_blank())


ggsave(filename=sprintf('%s/eda/custom/%s_NUMGENES=%s_PCARANK=%s/gRNA_all.png', 
                        save_dir, dimred_type, NUM_IMPORTANT_GENES, PCA_RANK),
       height = 5, width = 6, dpi = 320)

ggsave(filename=sprintf('%s/eda/custom/%s_NUMGENES=%s_PCARANK=%s/gRNA_all.pdf', 
                        save_dir, dimred_type, NUM_IMPORTANT_GENES, PCA_RANK),
       height = 5, width = 6)



# UMAP colored by potential confounders

# overall_subsamplesize = 20729
overall_subsamplesize = 7000

# bio_rep
set.seed(12345)
subsamplesize = overall_subsamplesize
dimred_type = 'pca_umap'
covariate_name = 'bio_rep'
allcells = pca_umap_cells$allcells # 20729 rows
allcells = allcells[sample(1:nrow(allcells), size = subsamplesize, replace = FALSE), ]
x_label = 'UMAP 1'
y_label = 'UMAP 2'



ggplot(allcells,
       aes(x = dim1, y = dim2, 
           color=eval(parse(text = covariate_name)))) +
  geom_point(alpha = .6, size = 1.5, shape = 19, stroke=NA) +
  labs(title = paste0('Bio Rep'),
       x = x_label, y = y_label,
       color = 'Bio Rep') +
  # scale_color_discrete(labels = c('rep 1', 'rep 2', 'rep 3')) +
  scale_color_brewer(palette = 'Dark2', labels = c('rep 1', 'rep 2', 'rep 3')) +
  guides(color = guide_legend(override.aes=list(shape = 15, size = 3.5, alpha = 1))) +
  theme(plot.title = element_text(hjust = .5),
        legend.position = c(.8, .85),
        legend.box.background = element_rect(colour = "black", linejoin = "bevel", linewidth = .5), 
        legend.margin = margin(t = 2, r = 4, b = 2, l = 4, unit = "pt"), 
        panel.background = element_rect(color = 'black', linewidth = 1), 
        axis.line = element_blank(),
        legend.key.size = unit(.45, 'cm'))

ggsave(filename=sprintf('%s/eda/custom/%s_NUMGENES=%s_PCARANK=%s/covariate_%s.png', 
                        save_dir, dimred_type, NUM_IMPORTANT_GENES, PCA_RANK, covariate_name),
       height = 5, width = 6, dpi = 320)

ggsave(filename=sprintf('%s/eda/custom/%s_NUMGENES=%s_PCARANK=%s/covariate_%s.pdf', 
                        save_dir, dimred_type, NUM_IMPORTANT_GENES, PCA_RANK, covariate_name),
       height = 5, width = 6)




# lane
set.seed(12345)
subsamplesize = overall_subsamplesize
dimred_type = 'pca_umap'
covariate_name = 'lane'
allcells = pca_umap_cells$allcells # 20729 rows
allcells = allcells[sample(1:nrow(allcells), size = subsamplesize, replace = FALSE), ]
x_label = 'UMAP 1'
y_label = 'UMAP 2'



ggplot(allcells,
       aes(x = dim1, y = dim2, 
           color=eval(parse(text = covariate_name)))) +
  geom_point(alpha = .4, size = 1.5, shape = 19, stroke=NA) +
  labs(title = paste0('Lane'),
       x = x_label, y = y_label,
       color = 'Lane') +
  scale_color_discrete(labels = paste0('Lane ', 1:8)) +
  guides(color = guide_legend(override.aes=list(shape = 15, size = 3.5, alpha = 1))) +
  theme(plot.title = element_text(hjust = .5),
        legend.position = c(.8, .75),
        legend.box.background = element_rect(colour = "black", linejoin = "bevel", linewidth = .5), 
        legend.margin = margin(t = 2, r = 4, b = 2, l = 4, unit = "pt"), 
        panel.background = element_rect(color = 'black', linewidth = 1), 
        axis.line = element_blank(),
        legend.key.size = unit(.45, 'cm'))

ggsave(filename=sprintf('%s/eda/custom/%s_NUMGENES=%s_PCARANK=%s/covariate_%s.png', 
                        save_dir, dimred_type, NUM_IMPORTANT_GENES, PCA_RANK, covariate_name),
       height = 5, width = 6, dpi = 320)

ggsave(filename=sprintf('%s/eda/custom/%s_NUMGENES=%s_PCARANK=%s/covariate_%s.pdf', 
                        save_dir, dimred_type, NUM_IMPORTANT_GENES, PCA_RANK, covariate_name),
       height = 5, width = 6)





# cell phase
set.seed(12345)
subsamplesize = overall_subsamplesize
dimred_type = 'pca_umap'
covariate_name = 'phase'
allcells = pca_umap_cells$allcells # 20729 rows
allcells = allcells[sample(1:nrow(allcells), size = subsamplesize, replace = FALSE), ]
x_label = 'UMAP 1'
y_label = 'UMAP 2'



ggplot(allcells,
       aes(x = dim1, y = dim2, 
           color=eval(parse(text = covariate_name)))) +
  geom_point(alpha = .5, size = 1.5, shape = 19, stroke=NA) +
  labs(title = paste0('Cell Phase'),
       x = x_label, y = y_label,
       color = 'Cell Phase') +
  scale_color_brewer(palette = 'Dark2') +
  guides(color = guide_legend(override.aes=list(shape = 15, size = 3.5, alpha = 1))) +
  theme(plot.title = element_text(hjust = .5),
        legend.position = c(.75, .83),
        legend.box.background = element_rect(colour = "black", linejoin = "bevel", linewidth = .5), 
        legend.margin = margin(t = 2, r = 4, b = 2, l = 4, unit = "pt"), 
        panel.background = element_rect(color = 'black', linewidth = 1), 
        axis.line = element_blank(),
        legend.key.size = unit(.45, 'cm'))

ggsave(filename=sprintf('%s/eda/custom/%s_NUMGENES=%s_PCARANK=%s/covariate_%s.png', 
                        save_dir, dimred_type, NUM_IMPORTANT_GENES, PCA_RANK, covariate_name),
       height = 5, width = 6, dpi = 320)

ggsave(filename=sprintf('%s/eda/custom/%s_NUMGENES=%s_PCARANK=%s/covariate_%s.pdf', 
                        save_dir, dimred_type, NUM_IMPORTANT_GENES, PCA_RANK, covariate_name),
       height = 5, width = 6)


# n_nonzero
set.seed(12345)
subsamplesize = overall_subsamplesize
dimred_type = 'pca_umap'
covariate_name = 'n_nonzero'
allcells = pca_umap_cells$allcells # 20729 rows
allcells = allcells[sample(1:nrow(allcells), size = subsamplesize, replace = FALSE), ]
x_label = 'UMAP 1'
y_label = 'UMAP 2'




ggplot(allcells,
       aes(x = dim1, y = dim2, 
           color=eval(parse(text = sprintf('(%s)', covariate_name))))) +
  geom_point(alpha = .7, size = 1.5, shape = 19, stroke=NA) +
  labs(title = paste0('# Nonzero Genes'),
       x = x_label, y = y_label,
       color = '# nonzero') +
  scale_color_viridis(discrete = FALSE, option = "magma", direction = -1,
                      begin = 0, breaks = c(0, 2000, 4000, 6000)) +
  theme(plot.title = element_text(hjust = .5),
        legend.position = c(.76, .78),
        legend.box.background = element_rect(colour = "black", linejoin = "bevel", linewidth = .5), 
        legend.margin = margin(t = 2, r = 4, b = 4, l = 4, unit = "pt"), 
        panel.background = element_rect(color = 'black', linewidth = 1), 
        axis.line = element_blank(),
        legend.key.size = unit(.45, 'cm'))

ggsave(filename=sprintf('%s/eda/custom/%s_NUMGENES=%s_PCARANK=%s/covariate_%s.png', 
                        save_dir, dimred_type, NUM_IMPORTANT_GENES, PCA_RANK, covariate_name),
       height = 5, width = 6, dpi = 320)

ggsave(filename=sprintf('%s/eda/custom/%s_NUMGENES=%s_PCARANK=%s/covariate_%s.pdf', 
                        save_dir, dimred_type, NUM_IMPORTANT_GENES, PCA_RANK, covariate_name),
       height = 5, width = 6)



# p_mito
set.seed(12345)
subsamplesize = overall_subsamplesize
dimred_type = 'pca_umap'
covariate_name = 'p_mito'
allcells = pca_umap_cells$allcells # 20729 rows
allcells = allcells[sample(1:nrow(allcells), size = subsamplesize, replace = FALSE), ]
x_label = 'UMAP 1'
y_label = 'UMAP 2'




ggplot(allcells,
       aes(x = dim1, y = dim2, 
           color=eval(parse(text = sprintf('(%s)', covariate_name))))) +
  geom_point(alpha = .7, size = 1.5, shape = 19, stroke=NA) +
  labs(title = paste0('Proportion Mitochondria'),
       x = x_label, y = y_label,
       color = 'prop mito') +
  scale_color_viridis(discrete = FALSE, option = "magma", direction = -1,
                      breaks = c(.025, .050, .075)) +
  theme(plot.title = element_text(hjust = .5),
        legend.position = c(.77, .79),
        legend.box.background = element_rect(colour = "black", linejoin = "bevel", linewidth = .5), 
        legend.margin = margin(t = 2, r = 4, b = 4, l = 4, unit = "pt"), 
        panel.background = element_rect(color = 'black', linewidth = 1), 
        axis.line = element_blank(),
        legend.key.size = unit(.45, 'cm'))

ggsave(filename=sprintf('%s/eda/custom/%s_NUMGENES=%s_PCARANK=%s/covariate_%s.png', 
                        save_dir, dimred_type, NUM_IMPORTANT_GENES, PCA_RANK, covariate_name),
       height = 5, width = 6, dpi = 320)

ggsave(filename=sprintf('%s/eda/custom/%s_NUMGENES=%s_PCARANK=%s/covariate_%s.pdf', 
                        save_dir, dimred_type, NUM_IMPORTANT_GENES, PCA_RANK, covariate_name),
       height = 5, width = 6)











# =================== END ======================================================
print(sprintf("[%s] END", Sys.time()))






