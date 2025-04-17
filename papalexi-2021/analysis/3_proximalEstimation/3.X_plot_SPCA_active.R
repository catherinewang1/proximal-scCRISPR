# ---------------------------------------------------------------------------- #
#              Plot Active p-values (using sPCA of genes as NCE/NCO)           #
# Requires: prev saved active p-values
# Ouputs: (nothing) but saves                                                  #
# ---------------------------------------------------------------------------- #
args = commandArgs(trailingOnly = TRUE)
# args = c('laptop', 'spca/cbgenes', 'allPerturbations', 'simple', '0')
args = c('laptop', 'spca/cbgenes', 'A', 'simple', 'NA')

GAMMA = 0


suppressPackageStartupMessages(require(assertthat)) # for some assert statements
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(ggplot2))    # plotting
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(RColorBrewer))




theme_set(theme_cowplot() +
            theme(plot.title = element_text(hjust = .5),
                  plot.subtitle = element_text(hjust = .5),
                  strip.background = element_rect(color = 'black', fill = 'white')))



assertthat::assert_that(length(args) > 0, msg="must give arg for specifying device eg 'Rscript <filename>.R ubergenno'")

DEVICE = args[1]
source('../PATHS.R') # load in data_dir and save_dir and CODE_DIR, depending on DEVICE value
# # location of papalexi-2021 folder
# data_dir = switch(DEVICE,
#                   'laptop'='C:/Users/Cathe/Documents/School/Genetic Hypothesis Testing using Negative Controls/genData/papalexi', 
#                   'desktop'='C:/Users/Catherine W/Documents/Research/genData/papalexi', 
#                   'ubergenno'='/raid6/Catherine/papalexi')
# # location of intermediate save files/plots are written
# save_dir = switch(DEVICE,
#                   'laptop'='C:/Users/Cathe/Documents/School/Genetic Hypothesis Testing using Negative Controls/DoubleBridge/saves/papalexi_saves', 
#                   'desktop'='C:/Users/Catherine W/Documents/Research/DoubleBridge/saves/papalexi_saves', 
#                   'ubergenno'='/raid6/Catherine/papalexi/papalexi_saves')
# # location of utils code
# util_dir = switch(DEVICE,
#                   'laptop'='C:/Users/Cathe/Documents/School/Genetic Hypothesis Testing using Negative Controls/DoubleBridge/code/utils', 
#                   'desktop'='C:/Users/Catherine W/Documents/Research/DoubleBridge/code/utils', 
#                   'ubergenno'='/home/catheri2/DoubleBridge/code/utils')
# 
# assertthat::assert_that(!is.null(data_dir), msg='first arg must be: laptop, desktop, or ubergenno')
assertthat::assert_that(length(args) > 1, msg="must give arg for specifying cbgenes or cbgenes_pooled name 'Rscript <filename>.R ubergenno cbgenes C'")
cbgenes_setting_name = args[2]
assertthat::assert_that(cbgenes_setting_name == 'spca/cbgenes', msg="this plotting script is ONLY for SPCA results, with 1 ATE est per AY pair")

assertthat::assert_that(length(args) > 1, msg="must give arg for specifying chosen AYZW name 'Rscript <filename>.R ubergenno C'")
AYZW_setting_name = args[3]

assertthat::assert_that(length(args) > 2, msg="must give arg for specifying chosen CB Setting name 'Rscript <filename>.R ubergenno cbgenes C default'")
CB_setting_name = args[4]

assertthat::assert_that(length(args) > 3, msg="must give arg for specifying lmYA threshold 'Rscript <filename>.R ubergenno cbgenes C default .1'")
lmYA_threshold = as.numeric(args[5])

# =================== Start ====================================================
print(sprintf("[%s] START: CB Plots", Sys.time()))

# load chosen AYZW names
AY   = read.csv(sprintf('%s/%s/%s/AY.csv', save_dir, cbgenes_setting_name, AYZW_setting_name))
# AYZW = readRDS(sprintf('%s/%s/%s/AYZW.rds', save_dir, cbgenes_setting_name, AYZW_setting_name)) # none for spca

# load gene info (top XX genes, removing TF and Targets)
gene_TF_Target = read.csv(sprintf('%s/important_genes_name_TF_Target.csv', save_dir)) # top important genes noting the TF and gRNA targets






pvals = read.csv(sprintf('%s/spca/cbgenes/%s/%s/ATE_activearbdep_gamma=%.2f.csv', save_dir, AYZW_setting_name, CB_setting_name, GAMMA))

pvals = readRDS(sprintf('%s/spca/cbgenes/%s/%s/ATE_activearbdep_gamma=%.2f.rds', save_dir, AYZW_setting_name, CB_setting_name, GAMMA))

nrow(AY); AY$type|>table()
nrow(pvals)
AY$AY_idx = 1:nrow(AY)
pvals = merge(AY, pvals, by = 'AY_idx')


# Nicer Labels for type
type_map = c("negative" = "Non-Causal A-Y",
             "maybe"    = "Candidate A-Y",
             "positive" = "Causal A-Y")
type_labeller = labeller(type = type_map)



# =================== Make Plots ======================================
print(sprintf("[%s]    - Make Plots", Sys.time()))

plot_savepath = sprintf('%s/%s/%s/%s/plotsActive/lmYAthresh=%05.02f', save_dir, cbgenes_setting_name, AYZW_setting_name, CB_setting_name, lmYA_threshold)
dir.create(plot_savepath, showWarnings = FALSE, recursive = TRUE)



# active pvals
p_active = ggplot(pvals,
       aes(x = pval_active)) +
  geom_histogram() +
  labs(title = "Active") +
  facet_grid(vars(type))


# proxy pvals
p_proxy = ggplot(pvals,
       aes(x = pval_proxy)) +
  geom_histogram() +
  labs(title = "Proxy") +
  facet_grid(vars(type))

# true pvals
p_true = ggplot(pvals,
       aes(x = pval_true)) +
  geom_histogram() +
  labs(title = "True") +
  facet_grid(vars(type))


pdf(file = sprintf('%s/histogram.pdf', plot_savepath), width = 7, height = 4); grid.arrange(p_active, p_proxy, p_true, nrow = 1); dev.off()



pvals_tall = pvals |> 
             select(type, pval_active, pval_proxy, pval_true) |> 
             reshape::melt(id = 'type', variable_name = 'pval_type') |>
             mutate(pval_type = mapply(function(x) {substring(x, first = 6) |> stringr::str_to_title()},
                                      pval_type)) |>
             rename(pval = value)


# histogram
ggplot(pvals_tall, aes(x = pval)) +
  geom_histogram(aes(y = after_stat(density)), binwidth = .05) +
  labs(title = 'Histogram') +
  facet_grid(rows = vars(type), cols = vars(pval_type), scales = 'free_y')

ggsave(file = sprintf('%s/histogram.pdf', plot_savepath), width = 7, height = 4)



# QQ Uniform
ggplot(pvals_tall, aes(sample = pval, color = pval_type)) +
  geom_qq(aes(), distribution = stats::qunif, alpha = .4) +
  geom_abline(aes(intercept = 0, slope = 1)) +
  labs(title = 'QQ Uniform') +
  facet_grid(rows = vars(type), 
             cols = vars(pval_type),
             scales = 'free_y')
ggsave(file = sprintf('%s/qqunif.pdf', plot_savepath), width = 7, height = 4)











