# ---------------------------------------------------------------------------- #
#                   Active p-values (using sPCA of genes as NCE/NCO)           #
# using already saved proxy and true pvalues
# Requires: prev saved proxy and true pvals                                    #
# Ouputs: (nothing) but saves                                                  #
# ---------------------------------------------------------------------------- #
args = commandArgs(trailingOnly = TRUE)
# args = c('laptop', 'spca/cbgenes', 'A', 'simple', 'NA')

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
assertthat::assert_that(!is.null(data_dir), msg='first arg must be: laptop, desktop, or ubergenno')

assertthat::assert_that(length(args) > 1, msg="must give arg for specifying cbgenes or cbgenes_pooled name 'Rscript <filename>.R ubergenno cbgenes C'")
cbgenes_setting_name = args[2]
assertthat::assert_that(cbgenes_setting_name == 'spca/cbgenes', msg="this plotting script is ONLY for SPCA results, with 1 ATE est per AY pair")

assertthat::assert_that(length(args) > 1, msg="must give arg for specifying chosen AYZW name 'Rscript <filename>.R ubergenno C'")
AYZW_setting_name = args[3]

assertthat::assert_that(length(args) > 2, msg="must give arg for specifying chosen CB Setting name 'Rscript <filename>.R ubergenno cbgenes C default'")
CB_setting_name = args[4]

assertthat::assert_that(length(args) > 3, msg="must give arg for specifying lmYA threshold 'Rscript <filename>.R ubergenno cbgenes C default .1'")
lmYA_threshold = as.numeric(args[5])


# dir.create(sprintf('%s/spca/cbgenes/%s/%s', save_dir, AYZW_setting_name, CB_setting_name), recursive = TRUE, showWarnings = FALSE)


# =================== Start ====================================================
print(sprintf("[%s] START: CB Estimate", Sys.time()))

# source CB utility functions (mainly CB fn for estimating)
source(sprintf('%s/CBEstSPCAActive.R', util_dir)) # for active arbdep with SPCA

# load chosen AYZW names and saved pvals
AY   = read.csv(sprintf('%s/%s/%s/AY.csv', save_dir, cbgenes_setting_name, AYZW_setting_name))
pvals = read.csv(sprintf('%s/spca/cbgenes/%s/%s/ATE_activearbdep_gamma=%.2f.csv', save_dir, AYZW_setting_name, CB_setting_name, 0)) # GAMMA = 0



nrow(AY); AY$type|>table()
nrow(pvals)
AY$AY_idx = 1:nrow(AY)
pvals = merge(AY, pvals, by = 'AY_idx')





# get active p-values for different gamma values
GAMMAs = seq(from = 0, to = 1, length.out = 5)[1:4]
pvals_GAMMAs = NULL
for(GAMMA in GAMMAs) {
  get_active_pval_time = get_active_arbdep_pval_using_saved_make(pvals_df = pvals, gamma = GAMMA)
  
  # mapply output format is weird
  # pvals_GAMMA = mapply(get_active_pval_time,
  #                        AY_idx = pvals$AY_idx, SIMPLIFY = TRUE)
  # pvals_GAMMAs = rbind(pvals_GAMMAs,
  #                      t(pvals_GAMMA) |> as.data.frame())
  
  
  # format into data.frame properly
  for(AY_idx in pvals$AY_idx) {
    pvals_GAMMAs = rbind(pvals_GAMMAs, 
                         get_active_pval_time(AY_idx = AY_idx) |> data.frame())
  }
}
pvals_GAMMAs = merge(AY, pvals_GAMMAs, by = 'AY_idx') # add AY info (e.g. test type)
pvals_GAMMAs = pvals_GAMMAs |> group_by(type, gamma) |> mutate(observed = sort(pval_active), expected = ppoints(n())) |> ungroup()




# =================== Make Plots ======================================
print(sprintf("[%s]    - Make Plots", Sys.time()))

plot_savepath = sprintf('%s/%s/%s/%s/plotsActive/lmYAthresh=%05.02f', save_dir, cbgenes_setting_name, AYZW_setting_name, CB_setting_name, lmYA_threshold)
dir.create(plot_savepath, showWarnings = FALSE, recursive = TRUE)



# prob of querying true
ggplot(pvals_GAMMAs, 
       aes(x = queried_true_prob)) +
  geom_histogram(aes(y = after_stat(density)), bins = 25)  +
  facet_grid(type ~ gamma, scales = 'free_y')
ggsave(file = sprintf('%s/gamma_queryprob.pdf', plot_savepath), width = 7, height = 4)


# histogram of active
ggplot(pvals_GAMMAs, 
       aes(x = pval_active)) +
  geom_histogram(aes(y = after_stat(density)), bins = 25) +
  facet_grid(type ~ gamma, scales = 'free_y')
ggsave(file = sprintf('%s/gamma_histogram.pdf', plot_savepath), width = 7, height = 4)


# plot QQ-plot (against Unif(0,1)) of ACTIVE pvals across differnt gamma values
ggplot(pvals_GAMMAs, 
       aes(x = expected,
           y = observed, group = gamma, color = GAMMA)) +
  geom_point() +
  facet_grid(type ~ gamma)
ggsave(file = sprintf('%s/gamma_qqunif_faceted.pdf', plot_savepath), width = 7, height = 4)

ggplot(pvals_GAMMAs, 
       aes(x = expected,
           y = observed, group = gamma, color = gamma)) +
  geom_abline(aes(intercept = 0, slope = 1)) +
  geom_point(alpha = .8) +
  scale_color_viridis_b(begin = 0, end = 1)+
  facet_grid(~type)
ggsave(file = sprintf('%s/gamma_qqunif.pdf', plot_savepath), width = 7, height = 4)

