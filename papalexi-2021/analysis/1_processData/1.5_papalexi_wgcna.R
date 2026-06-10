# ------------------------------------------------------------------------------------------------ #
#                        Construct Negative Controls using WGCNA Clusters 
#     Performs WGCNA Clustering (package: WGCNA, Langfelder and Horvath 2008 doi.org/10.1186/1471-2105-9-559). 
#     Saves the loadings to be used as Negative Controls (Exposures and Outcomes) in Prox Inf
# Requires: normalized values
# Outputs: saves clusters and averages at 
#     sprintf('%s/spca/NCloadings_sumabs=%.1f_K=%d_N=%d.rds', save_dir, my_sumabsv, my_K, N_subsample)
#     sprintf('%s/spca/NCavg_sumabs=%.1f_K=%d_N=%d.rds', save_dir, my_sumabsv, my_K, N_subsample)
# ------------------------------------------------------------------------------------------------ #
args = commandArgs(trailingOnly = TRUE)
args = c('laptop', '4000')

POWERS = c(1, 2, 3, 4, 6)



suppressWarnings(suppressMessages(require(assertthat))) # for some assert statements
suppressWarnings(suppressMessages(library(WGCNA)))      # WGCNA analysis
suppressWarnings(suppressMessages(require(rhdf5)))      # read/write HDF5 format
suppressWarnings(suppressMessages(require(cowplot)))
suppressWarnings(suppressMessages(library(ggplot2)))
suppressWarnings(suppressMessages(library(RColorBrewer)))

options(stringsAsFactors = FALSE) # The following setting is important, do not omit.
WGCNA::enableWGCNAThreads()



assertthat::assert_that(length(args) > 0, msg="must give arg for specifying device eg 'Rscript <filename>.R ubergenno'")
DEVICE = args[1]
source('../PATHS.R') # load in data_dir and save_dir, depending on DEVICE value

assertthat::assert_that(length(args) >= 2, msg="provide device and NUM_IMPORTANT_GENES (likely int <=4000)") 
NUM_IMPORTANT_GENES = as.integer(args[2]) # should be <= 1.1 previous script param

# assertthat::assert_that(length(args) > 0, msg="must give arg for specifying power(s) 'Rscript <filename>.R ubergenno 4000 1 3'")
# POWERS = as.numeric(args[3:length(args)])
# assertthat::assert_that(all(POWERS >= 1), msg="specified POWERS must be >= 1 for default WGCNA::blockwiseModules function")

# =================== Start ========================================================================
print(sprintf("[%s] START: WGCNA gene modules", Sys.time()))



# ==================================================================================================
# =================== Load Normalized Gene Expr ====================================================
# ==================================================================================================
print(sprintf("[%s]    - Loading Normalized Gene Expression w/ %d Genes", Sys.time(), NUM_IMPORTANT_GENES))

# =================== loading top genes ===========================================================
h5file      = paste0(save_dir, "/gene.h5"); print(h5file)
reading_hd5file  = rhdf5::H5Fopen(name = h5file)
# readin_gene_norm = reading_hd5file&'gene_norm_removeconfounders' # prev version removed known confounders, but we are pretending we don't know confounders
readin_gene_norm = reading_hd5file&'gene_norm'
gene_norm = readin_gene_norm[1:NUM_IMPORTANT_GENES, ] # dim = 4000 x 20729 = #important x #cells
rhdf5::h5closeAll()
invisible(gc(verbose=FALSE))

datExpr = t(gene_norm) # follow tutorial's naming scheme
rm(gene_norm)


# ==================================================================================================
# =================== Automatically Construct Gene Clusters ========================================
# ==================================================================================================

# =================== Automatically Construct Gene Clusters ========================================================================
print(sprintf("[%s]    - WGCNA cluster genes automatically", Sys.time()))
# print(sprintf("[%s]    - automatically construct gene clusters", Sys.time()))

assertthat::assert_that(dir.exists(save_dir), msg = 'save_dir must exist, make sure you are in the right folder')
WGCNA_results_savepath = sprintf('%s/WGCNA/autoWGCNA', save_dir)
# create folders for saving results
dir.create(sprintf('%s/modules', WGCNA_results_savepath), recursive = T, showWarnings = F)

dir.create(sprintf('%s/modulesDendrograms/',    WGCNA_results_savepath), showWarnings = F)
dir.create(sprintf('%s/modulesDendrograms/pdf', WGCNA_results_savepath), showWarnings = F)
dir.create(sprintf('%s/modulesDendrograms/png', WGCNA_results_savepath), showWarnings = F)
dir.create(sprintf('%s/modulesDendrograms/svg', WGCNA_results_savepath), showWarnings = F)

dir.create(sprintf('%s/modulesHeatmaps/',    WGCNA_results_savepath), showWarnings = F)
dir.create(sprintf('%s/modulesHeatmaps/pdf', WGCNA_results_savepath), showWarnings = F)
dir.create(sprintf('%s/modulesHeatmaps/png', WGCNA_results_savepath), showWarnings = F)
dir.create(sprintf('%s/modulesHeatmaps/svg', WGCNA_results_savepath), showWarnings = F)

dir.create(sprintf('%s/modulesEigengeneCor/',    WGCNA_results_savepath), showWarnings = F)
dir.create(sprintf('%s/modulesEigengeneCor/pdf', WGCNA_results_savepath), showWarnings = F)
dir.create(sprintf('%s/modulesEigengeneCor/png', WGCNA_results_savepath), showWarnings = F)
dir.create(sprintf('%s/modulesEigengeneCor/svg', WGCNA_results_savepath), showWarnings = F)

dir.create(sprintf('%s/savedTOMs/', WGCNA_results_savepath), showWarnings = F)

# get gene modules for each power/beta
num_modules = c()
for(pow in POWERS) {
  dir.create(sprintf('%s/savedTOMs/TOM_NUMGENES=%04.0f_POWER=%01.1f', WGCNA_results_savepath, NUM_IMPORTANT_GENES, pow), 
             recursive = TRUE, showWarnings = FALSE)
  
  
  net = WGCNA::blockwiseModules(datExpr, power = pow,
                                TOMType = "unsigned", 
                                # Basic tree cut options
                                # try to create many but potentially smaller modules
                                deepSplit = 4,  # (0-4 where 4 is more sensitive and creates smaller )
                                detectCutHeight = .9999,
                                minModuleSize = 20,
                                # Gene reassignment, module trimming, and module "significance" criteria
                                reassignThreshold = 1e-6,
                                minCoreKME = 0.2, 
                                minCoreKMESize = 5,
                                minKMEtoStay = 0.15,
                                # Module merging options
                                mergeCutHeight = .10, # mergeCutHeight = 0.25, # keep low to make fewer merges
                                # rest
                                numericLabels = TRUE, pamRespectsDendro = FALSE,
                                saveTOMs = TRUE,
                                saveTOMFileBase = sprintf('%s/savedTOMs/TOM_NUMGENES=%04.0f_POWER=%01.1f/', WGCNA_results_savepath, NUM_IMPORTANT_GENES, pow), 
                                verbose = 3, nThreads = 8)
  
  saveRDS(net, sprintf('%s/modules/WGCNAmodules_NUMGENES=%04.0f_POWER=%01.1f.rds', WGCNA_results_savepath, NUM_IMPORTANT_GENES, pow))
  
  
  num_modules = c(num_modules, (net$colors |> table() |> length()) - 1) # keep track of how many modules found
  
  ########################################################
  # Plot the dendrogram and the module colors underneath
  mergedColors = WGCNA::labels2colors(net$colors) # Convert labels to colors for plotting
  plot_dendro <- function() {
    plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                        "Module colors",
                        dendroLabels = FALSE, hang = 0.03,
                        addGuide = TRUE, guideHang = 0.05,
                        main = sprintf('Dendrogram (default WGCNA, power=%01.1f)', pow))
  }
  
  pdf(sprintf('%s/modulesDendrograms/pdf/wgcnaDendro_NUMGENES=%04.0f_POWER=%01.1f.pdf', WGCNA_results_savepath, NUM_IMPORTANT_GENES, pow),
      width = 10, height = 6)
  plot_dendro()
  dev.off()
  png(sprintf('%s/modulesDendrograms/png/wgcnaDendro_NUMGENES=%04.0f_POWER=%01.1f.png', WGCNA_results_savepath, NUM_IMPORTANT_GENES, pow),
      width = 1280, height = 720, units = 'px', pointsize = 30)
  plot_dendro()
  dev.off()
  svg(sprintf('%s/modulesDendrograms/svg/wgcnaDendro_NUMGENES=%04.0f_POWER=%01.1f.svg', WGCNA_results_savepath, NUM_IMPORTANT_GENES, pow),
      width = 10, height = 6)
  plot_dendro()
  dev.off()
  
  
  ########################################################
  # Plot correlation heatmap
  
  # reorder based on gene module number (sort by module num, then gene import num)
  sorted_cols = sort(paste0(net$colors, '_', sprintf('%05.f', 1:NUM_IMPORTANT_GENES)), index.return=TRUE)
  # correlation matrix
  datExprCor = cor(datExpr[,sorted_cols$ix])
  # function makes plot of Correlations btw top genes ordered by module
  plot_corheatmap <- function() {
    heatmap(abs(datExprCor), Rowv = NA, Colv = NA, scale = 'none',
            labRow = NA, labCol = NA, margins = c(0, 0),
            ColSideColors = mergedColors[sorted_cols$ix],
            RowSideColors = mergedColors[sorted_cols$ix],
            col = colorRampPalette(c('white', 'darkred'))(n=256),
            zlim = c(0, 1), # color range between 0 and 1 (correlations)
            main=sprintf('abs(Correlation)\n(%04.f genes, sorted by modules using power=%01.1f)', NUM_IMPORTANT_GENES, pow))     
  }
  
  pdf(sprintf('%s/modulesHeatmaps/pdf/wgcnaHeatmap_NUMGENES=%04.0f_POWER=%01.1f.pdf', WGCNA_results_savepath, NUM_IMPORTANT_GENES, pow),
      width = 7, height = 8)
  plot_corheatmap()
  dev.off()
  png(sprintf('%s/modulesHeatmaps/png/wgcnaHeatmap_NUMGENES=%04.0f_POWER=%01.1f.png', WGCNA_results_savepath, NUM_IMPORTANT_GENES, pow), 
      width = 2560, height = 2800, pointsize = 75, 
      units = "px")
  plot_corheatmap()
  dev.off()
  svg(sprintf('%s/modulesHeatmaps/svg/wgcnaHeatmap_NUMGENES=%04.0f_POWER=%01.1f.svg', WGCNA_results_savepath, NUM_IMPORTANT_GENES, pow),
      width = 7, height = 8)
  plot_corheatmap()
  dev.off()
  
  
  ########################################################
  # Plot correlation between Module EigenGenes
  p = ggplot(reshape2::melt(cor(net$MEs))) +
    geom_tile(aes(x = Var1, y = Var2, fill= value)) +
    geom_text(aes(x = Var1, y = Var2, label = sprintf('%.2f', value)),
              size = 10, color = 'black') +
    labs(title = sprintf('Correlation between Module Eigengenes'),
         subtitle = sprintf('%04.f genes, using power=%01.1f', NUM_IMPORTANT_GENES, pow)) +
    scale_x_discrete(expand = c(0, 0, 0, 0), limits=rev) +
    scale_y_discrete(expand = c(0, 0, 0, 0), limits=rev) +
    scale_fill_distiller(palette = "Spectral", limits = c(-1, 1)) +
    cowplot::theme_cowplot() +
    theme(axis.line = element_blank(), axis.title = element_blank(),
          legend.position = 'right',
          plot.title = element_text(hjust = .5),
          plot.subtitle = element_text(hjust = .5))
  
  ggsave(sprintf('%s/modulesEigengeneCor/png/wgcnaEigengeneCor_NUMGENES=%04.0f_POWER=%01.1f.png', 
                 WGCNA_results_savepath, NUM_IMPORTANT_GENES, pow), 
         plot = p, device = 'png',
         width = 1080, height = 1080, units = 'px',
         scale = 1.5,
         bg = 'white')
  ggsave(sprintf('%s/modulesEigengeneCor/pdf/wgcnaEigengeneCor_NUMGENES=%04.0f_POWER=%01.1f.pdf', 
                 WGCNA_results_savepath, NUM_IMPORTANT_GENES, pow), 
         plot = p, device = 'pdf',
         width = 6, height = 6,
         scale = 1,
         bg = 'white')
  ggsave(sprintf('%s/modulesEigengeneCor/svg/wgcnaEigengeneCor_NUMGENES=%04.0f_POWER=%01.1f.svg', 
                 WGCNA_results_savepath, NUM_IMPORTANT_GENES, pow), 
         plot = p, device = 'svg',
         width = 6, height = 6,
         scale = 1,
         bg = 'white')
}


# =================== END ======================================================
print(sprintf("[%s] END", Sys.time()))


