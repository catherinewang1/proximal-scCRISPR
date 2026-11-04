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
# - 2.2.2 deseq2
# - 2.2.3 causarray
# - 2.2.4 sceptre
# ------------------------------------------------------------------------------------------------ #
args = commandArgs(trailingOnly = TRUE)
# args = c('macbook')


assertthat::assert_that(length(args) > 0, msg="must give arg for specifying device eg 'Rscript <filename>.R ubergenno'")
DEVICE = args[1]
source('../PATHS.R') # load in data_dir and save_dir and CODE_DIR, depending on DEVICE value
assertthat::assert_that(!is.null(data_dir), msg='first arg must be: laptop, desktop, or ubergenno')


suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggh4x))
suppressPackageStartupMessages(library(cowplot))
theme_set(theme_cowplot() +
            theme(plot.title = element_text(hjust = .5),
                  plot.subtitle = element_text(hjust = .5),
                  strip.background = element_rect(fill = "white", color = 'black'),
                  panel.grid.minor = element_blank(),
                  panel.grid.major = element_blank()))
                                                                              

PARAMS = list( 
  'TAG'       = 'L45Glut',
  'subfolder' = '3_L45glutcells',
  'seed'      = 12345,
  'A_MIN_TRTCELLS'  = 40, # filter for grnas with at least this number of cells (keep low here, set higher later))
  'Y_MAX_GENE_IMPORTANCE' = 1500
)
save_dir = paste0(save_dir, '/', PARAMS[['subfolder']])
set.seed(PARAMS[['seed']])


h5ad_folder = sprintf('%s/causarray/jin2025_%s/', save_dir, PARAMS[['TAG']])
dir.create(sprintf('%s/plots/', h5ad_folder), showWarnings=FALSE) 



get_AY_test_type2 <- function(grna_name, gene_name) {
  # negative if 'Non_target' is in grna_name
  if(grepl('Non_target', grna_name)) {return('negative')} 
  else if(grna_name == gene_name) {return('positive')} 
  else {return('discovery')}
}

# ============================================================================== #
# =================== Wilcoxan ========================================
# ============================================================================== #
# load in from causarray's comparison analysis
wilc = anndataR::read_h5ad(sprintf('%s/wilcoxon.h5ad', h5ad_folder))

wilc_pvals = reshape2::melt(wilc$layers[['pvalue']] |> as.data.frame() |> tibble::rownames_to_column('grna'), id='grna', value.name = 'pval') |> dplyr::rename(gene=variable)
wilc_pvals$type_raw = mapply(FUN = get_AY_test_type2, grna_name = wilc_pvals$grna, gene_name = wilc_pvals$gene)
wilc_pvals$type = factor(wilc_pvals$type_raw,
                      levels = c('negative', 'discovery', 'positive'),
                      labels = c(    'Null', 'Discovery', 'Alternative'))
wilc_pvals = wilc_pvals |> dplyr::relocate(type, type_raw, .after = gene)
head(wilc_pvals)


ggplot(wilc_pvals, aes(x = pval, fill = type)) + 
  geom_histogram(aes(y=after_stat(..density..)), binwidth = .020) +
  # facet_grid(cols = vars(type), scales = 'free_y') +
  facet_grid2(cols = vars(type), independent = 'y', shrink = TRUE, scales='free_y') +
  scale_x_continuous(breaks = c(0, .25, .5, .75, 1), 
                     labels = c('0', '.25', '.5', '.75', '1')) +
  labs(title = 'Hist of p-values from Wilcoxan Rank Sum Test')



ggsave(sprintf('%s/plots/hist_pvals_wilcoxan.pdf', h5ad_folder), height = 3, width = 6)

# ============================================================================== #
# =================== nbglm ========================================
# ============================================================================== #

# load in from causarray's comparison analysis
# nbglm = anndataR::read_h5ad(sprintf('%s/causarray/jin2025_L45Glut/nbglm.h5ad', save_dir))
nbglm = anndataR::read_h5ad(sprintf('%s/nbglm.h5ad', h5ad_folder))

nbglm_pvals = reshape2::melt(nbglm$layers[['pvalue']] |> as.data.frame() |> tibble::rownames_to_column('grna'), id='grna', value.name = 'pval') |> dplyr::rename(gene=variable)
nbglm_pvals$type_raw = mapply(FUN = get_AY_test_type2, grna_name = nbglm_pvals$grna, gene_name = nbglm_pvals$gene)
nbglm_pvals$type = factor(nbglm_pvals$type_raw,
                          levels = c('negative', 'discovery', 'positive'),
                          labels = c(    'Null', 'Discovery', 'Alternative'))
nbglm_pvals = nbglm_pvals |> dplyr::relocate(type, type_raw, .after = gene)
head(nbglm_pvals)


ggplot(nbglm_pvals, aes(x = pval, fill = type)) + 
  geom_histogram(aes(y=after_stat(..density..)), binwidth = .020) +
  # facet_grid(cols = vars(type), scales = 'free_y') +
  facet_grid2(cols = vars(type), independent = 'y', shrink = TRUE, scales='free_y') +
  scale_x_continuous(breaks = c(0, .25, .5, .75, 1), 
                     labels = c('0', '.25', '.5', '.75', '1')) +
  labs(title = 'Hist of p-values from NB GLM Fit')

ggsave(sprintf('%s/plots/hist_pvals_nbglm.pdf', h5ad_folder), height = 3, width = 6)



# ============================================================================== #
# =================== causarray/min_counts40 ========================================
# ============================================================================== #

# load in from causarray
effects_causarray = read.csv(sprintf('%s/min_counts40/lfc.csv', h5ad_folder))


effects_causarray = effects_causarray |>
  dplyr::rename(pval = pvalue, ATE = log2fc, se = log2fc_se, gene = gene_names, grna=trt, n_trtmnt=n_treated, n_cntrl=n_control) |>
  dplyr::select(grna, gene, ATE, se, pval, padj, rej, n_trtmnt, n_cntrl, count_treated, count_control) |> 
  dplyr::mutate(time_sec = NA, NC_type = NA,
                numNC = NA, method='causarray'
  )

effects_causarray$type = mapply(FUN = get_AY_test_type2, grna_name = effects_causarray$grna, gene_name = effects_causarray$gene)
effects_causarray = effects_causarray |> relocate(type, .after=gene)
effects_causarray$type = factor(effects_causarray$type,
                          levels = c('negative', 'discovery', 'positive'),
                          labels = c(    'Null', 'Discovery', 'Alternative'))

ggplot(effects_causarray, aes(x = pval, fill = type)) + 
  geom_histogram(aes(y=after_stat(..density..)), binwidth = .020) +
  # facet_grid(cols = vars(type), scales = 'free_y') +
  facet_grid2(cols = vars(type), independent = 'y', shrink = TRUE, scales='free_y') +
  scale_x_continuous(breaks = c(0, .25, .5, .75, 1), 
                     labels = c('0', '.25', '.5', '.75', '1')) +
  labs(title = 'Hist of p-values from Causarray')

ggsave(sprintf('%s/plots/hist_pvals_causarray.pdf', h5ad_folder), height = 3, width = 6)


# ============================================================================== #
# =================== deseq2 ========================================
# ============================================================================== #

effects_deseq = read.csv(sprintf('%s/deseq2/deseq_effects.csv', save_dir)) 

effects_deseq = effects_deseq |>
  dplyr::rename(pval = pvalue, ATE = log2FoldChange, se = lfcSE) |>
  dplyr::select(grna, gene, ATE, se, pval) |> 
  dplyr::mutate(time_sec = NA, NC_type = NA,
                numNC = NA, method='deseq2'
  )
effects_deseq$type = mapply(FUN = get_AY_test_type2, grna_name = effects_deseq$grna, gene_name = effects_deseq$gene)
effects_deseq = effects_deseq |> relocate(type, .after=gene)
effects_deseq$type = factor(effects_deseq$type,
                                levels = c('negative', 'discovery', 'positive'),
                                labels = c(    'Null', 'Discovery', 'Alternative'))

ggplot(effects_deseq, aes(x = pval, fill = type)) + 
  geom_histogram(aes(y=after_stat(..density..)), binwidth = .020) +
  # facet_grid(cols = vars(type), scales = 'free_y') +
  facet_grid2(cols = vars(type), independent = 'y', shrink = TRUE, scales='free_y') +
  scale_x_continuous(breaks = c(0, .25, .5, .75, 1), 
                     labels = c('0', '.25', '.5', '.75', '1')) +
  labs(title = 'Hist of p-values from DESeq2')

ggsave(sprintf('%s/plots/hist_pvals_deseq2.pdf', h5ad_folder), height = 3, width = 6)





