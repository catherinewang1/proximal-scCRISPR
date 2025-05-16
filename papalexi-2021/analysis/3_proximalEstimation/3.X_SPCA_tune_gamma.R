# ---------------------------------------------------------------------------- #
#                   Active p-values (using sPCA of genes as NCE/NCO)           #
# using already saved proxy and true pvalues                                   #
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
suppressPackageStartupMessages(library(latex2exp))

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


# get indep active p-values
set.seed(1234567)
n_train = 1000
train_idx = sample(which(pvals$type == 'negative'))[1:n_train]  # split to train/test


pvals_indep_active = get_active_indep_pval_using_saved(pvals_df = pvals[-train_idx,], pval_proxy_null = pvals[train_idx, 'pval_proxy'])

pvals_indep_active = cbind(pvals[-train_idx, c('AY_idx', 'type', 'pval_proxy', 'pval_true')],
                           pvals_indep_active)


pvals_indep_active = pvals_indep_active |> 
                      rename(pval_active = pval_active_indep, 
                             time_active = time_active_indep,
                             time_active_expect = time_active_indep_expect)
# add expected under unif(0,1) for qqplots
pvals_indep_active = pvals_indep_active |> group_by(type) |> mutate(observed = sort(pval_active), expected = ppoints(n())) |> ungroup()
# Nicer Labels for type
pvals_indep_active = merge(pvals_indep_active, 
                     data.frame('type'      = c('negative', 'maybe', 'positive'),
                                # 'type_long' = c(TeX('Non-Causal AY'), TeX('Candidate AY'), TeX('Causal AY'))), 
                                # 'type_long' = c('Non-Causal', 'Candidate', 'Causal')), 
                                'type_long' = c('Null', 'Candidate', 'Alternative')), 
                     by = 'type')
pvals_indep_active$type_long = factor(pvals_indep_active$type_long, ordered = TRUE, levels = c('Null', 'Alternative'))

head(pvals_indep_active)



# get arbdep active p-values for different gamma values
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


# Nicer Labels for type
pvals_GAMMAs = merge(pvals_GAMMAs, 
                     data.frame('type'      = c('negative', 'maybe', 'positive'),
                                # 'type_long' = c(TeX('Non-Causal AY'), TeX('Candidate AY'), TeX('Causal AY'))), 
                                # 'type_long' = c('Non-Causal', 'Candidate', 'Causal')), 
                                'type_long' = c('Null', 'Candidate', 'Alternative')), 
                     by = 'type')
pvals_GAMMAs$type_long = factor(pvals_GAMMAs$type_long, ordered = TRUE, levels = c('Null', 'Alternative'))

# type_map = c("negative" = "Non-Causal AY",
#              "maybe"    = "Candidate AY",
#              "positive" = "Causal AY")
# type_labeller = labeller(type = type_map)
# pvals_GAMMAs$type = factor(pvals_GAMMAs$type)
# levels(pvals_GAMMAs$type) <- c("negative" = TeX("$$Non-Causal AY"),
#                                "maybe"    = TeX("Candidate AY"),
#                                "positive" = TeX("Causal AY"))

# nicer labels for gamma
pvals_GAMMAs$gamma_numeric = pvals_GAMMAs$gamma 
pvals_GAMMAs$gamma = factor(pvals_GAMMAs$gamma)
levels(pvals_GAMMAs$gamma) <- c('0' = TeX("$\\gamma=0$"), '.25' = TeX("$\\gamma=.25$"), '.5' = TeX("$\\gamma=.5$"), '.75' = TeX("$\\gamma=.75$"))


# =================== Make Plots ===============================================
print(sprintf("[%s]    - Make Plots", Sys.time()))

plot_savepath = sprintf('%s/%s/%s/%s/plotsActive/lmYAthresh=%05.02f', save_dir, cbgenes_setting_name, AYZW_setting_name, CB_setting_name, lmYA_threshold)
dir.create(plot_savepath, showWarnings = FALSE, recursive = TRUE)

# custom ggproto- to limit axis to 0,1
cutoff01 <- list(scale_y_continuous(limits = c(-.01, 1.01)),
                 scale_x_continuous(limits = c(-.01, 1.01)))

# prob of querying true ========================================================
p_histogram = ggplot(pvals_GAMMAs, 
       aes(x = queried_true_prob)) +
  geom_histogram(bins = 25)  + # aes(y = after_stat(density))
  labs(title = 'Probability of Querying True',
       x = 'Probability of Querying True p-value') +
  facet_grid(type_long ~ gamma, scales = 'free_y', labeller = label_parsed)

p_histogram
ggsave(file = sprintf('%s/gamma_queryprob.pdf', plot_savepath), width = 7, height = 4)

p_histogram + scale_x_continuous(limits = c(-.01, 1.01))
ggsave(file = sprintf('%s/gamma_queryprob01.pdf', plot_savepath), width = 7, height = 4)

# histogram of active ==========================================================
p_histogram = ggplot(pvals_GAMMAs, 
       aes(x = pval_active)) +
  geom_histogram(aes(y = after_stat(density)), binwidth = .1) +
  labs(title = 'Histogram of Active p-values',
       x = 'Active p-value') +
  facet_grid(type_long ~ gamma, scales = 'free_y', labeller = label_parsed)

p_histogram
ggsave(file = sprintf('%s/gamma_histogram.pdf', plot_savepath), width = 7, height = 4)

p_histogram + scale_x_continuous(limits = c(-.01, 1.01))
ggsave(file = sprintf('%s/gamma_histogram01.pdf', plot_savepath), width = 7, height = 4)



# plot QQ-plot (against Unif(0,1))  ============================================
# of ACTIVE pvals across differnt gamma values

# ggplot(pvals_GAMMAs, 
#        aes(x = expected,
#            y = observed, group = gamma, color = GAMMA)) +
#   geom_point() +
#   labs(title = 'QQ-plot of Active p-values vs Unif(0,1)',
#        x = 'Expected', y = 'Observed') +
#   facet_grid(type_long ~ gamma, labeller = label_parsed)
# ggsave(file = sprintf('%s/gamma_qqunif_faceted.pdf', plot_savepath), width = 7, height = 4)

p_qqplot <- ggplot(pvals_GAMMAs, 
       aes(x = expected,
           y = observed, color = gamma_numeric)) +
  geom_abline(aes(intercept = 0, slope = 1)) +
  geom_point(alpha = .8) +
  scale_color_viridis_b(begin = 0, end = 1, breaks = c(0, .2, .4, .6), labels = c(0, .25, .5, .75)) +
  labs(title = 'QQ-plot of Active p-values vs Unif(0,1)',
       x = 'Expected', y = 'Observed',
       color = TeX("$\\gamma$")) +
  facet_grid(~type_long, labeller = label_parsed) +
  theme(legend.text = element_text(vjust = -.7))

p_qqplot
ggsave(file = sprintf('%s/gamma_qqunif.pdf', plot_savepath), width = 7, height = 4)

p_qqplot + cutoff01
ggsave(file = sprintf('%s/gamma_qqunif01.pdf', plot_savepath), width = 7, height = 4)


# plot runtimes  ============================================

pvals_time = merge(pvals_GAMMAs,
                   pvals |> select(AY_idx, time_proxy, time_true),
                   by = 'AY_idx')

pvals_time_tall = reshape::melt(pvals_time, 
              id.vars = c('AY_idx', 'type_long', 'gamma'), 
              measure.vars = c('time_proxy', 'time_true', 'time_active_expect'),
              variable_name = 'time_type') |> rename(time = value)


p_runtime = ggplot(pvals_time_tall, 
       aes(x = time, fill = time_type)) +
  geom_histogram(aes(y = after_stat(density)), 
                 color = 'gray20',  bins = 18, alpha = .7, 
                 position = 'nudge')  + 
  scale_fill_brewer(palette = "Dark2", label = c("time_proxy" = "proxy",
                                                 "time_true"    = "true",
                                                 "time_active_expect" = "active")) +
  labs(title = 'Runtimes',
       x = 'time (s)',
       fill = 'test') +
  facet_grid(type_long ~ gamma, scales = 'free_y', labeller = label_parsed)

p_runtime
ggsave(file = sprintf('%s/gamma_runtime.pdf', plot_savepath), width = 7, height = 4)


labeller(type = c("time_proxy" = "proxy",
                                  "time_true"    = "true",
                                  "time_active_expect" = "active"))

# ==============================================================================
# Now include indep active p-values
# ==============================================================================
pvals_indep_active$gamma_numeric = -1
pvals_all = rbind(pvals_GAMMAs |> select(AY_idx, type_long, gamma_numeric, pval_active, time_active_expect, expected, observed),
                  pvals_indep_active |> select(AY_idx, type_long, gamma_numeric, pval_active, time_active_expect, expected, observed))
pvals_all$gamma = factor(pvals_all$gamma_numeric)
levels(pvals_all$gamma) <- c('-1' = TeX("$\\gamma=-1$"), '0' = TeX("$\\gamma=0$"), '.25' = TeX("$\\gamma=.25$"), '.5' = TeX("$\\gamma=.5$"), '.75' = TeX("$\\gamma=.75$"))
levels(pvals_all$gamma) <- c('-1' = TeX("Independent"), 
                             '0' = TeX("Arb Dep, $\\gamma=0$"), 
                             '.25' = TeX("Arb Dep, $\\gamma=.25$"), 
                             '.5' = TeX("Arb Dep, $\\gamma=.5$"), 
                             '.75' = TeX("Arb Dep, $\\gamma=.75$"))
gamma_colors = viridis::rocket(6)[2:5]

other_colors = viridis::viridis(4, begin = 0, end = 1)[2:4]

# form pvals_all2 that has proxy and true pvalues too
pvals1 = pvals
pvals1 = merge(pvals1, 
               data.frame('type'      = c('negative', 'maybe', 'positive'),
                          # 'type_long' = c(TeX('Non-Causal AY'), TeX('Candidate AY'), TeX('Causal AY'))), 
                          # 'type_long' = c('Non-Causal', 'Candidate', 'Causal')), 
                          'type_long' = c('Null', 'Candidate', 'Alternative')), 
               by = 'type')
pvals1$type_long = factor(pvals1$type_long, ordered = TRUE, levels = c('Null', 'Alternative'))

pvals1$gamma_numeric = -3
pvals_proxy = pvals1 |> rename(pval = pval_proxy, time = time_proxy) |> select(AY_idx, type_long, gamma_numeric, pval, time)
pvals1$gamma_numeric = -2
pvals_true  = pvals1 |> rename(pval = pval_true,  time = time_true)  |> select(AY_idx, type_long, gamma_numeric, pval, time)



pvals_all2 = rbind(pvals_proxy,
                   pvals_true,
                   pvals_all |> rename(pval = pval_active, time = time_active_expect) |> select(AY_idx, type_long, gamma_numeric, pval, time))
pvals_all2$gamma = factor(pvals_all2$gamma_numeric)
levels(pvals_all2$gamma) <- c('-3' = TeX("Proxy"),
                              '-2' = TeX("True"),
                             '-1' = TeX("Independent"), 
                             '0' = TeX("Arb Dep, $\\gamma=0$"), 
                             '.25' = TeX("Arb Dep, $\\gamma=.25$"), 
                             '.5' = TeX("Arb Dep, $\\gamma=.5$"), 
                             '.75' = TeX("Arb Dep, $\\gamma=.75$"))

pvals_all2 = pvals_all2 |> group_by(type_long, gamma) |> mutate(observed = sort(pval), expected = ppoints(n())) |> ungroup()

# qqplot of active pvals vs Unif(0,1) ==========================================
qqplot_both = ggplot(pvals_all, 
                   aes(x = expected,
                       y = observed,
                       color = gamma)) +
  geom_abline(aes(intercept = 0, slope = 1)) +
  geom_point(alpha = .8, size = .8) +
  # scale_color_viridis_b(begin = 0, end = 1, breaks = c(-1, 0, .2, .4, .6), labels = c('Indep', 0, .25, .5, .75)) +
  scale_color_manual(values = c('gray32', gamma_colors), 
                     labels = c(TeX('Independent'), TeX('Arb Dep, \\gamma=0'), TeX('Arb Dep, \\gamma=.25'), TeX('Arb Dep, \\gamma=.50'), TeX('Arb Dep, \\gamma=.75'))) +
  labs(title = 'QQ-plot of Active p-values vs Unif(0,1)',
       x = 'Expected', y = 'Observed',
       color = 'active pval type') +
  facet_grid(~type_long, labeller = label_parsed) +
  theme(legend.text = element_text(vjust = .7, size = 10),
        legend.title = element_text(size = 12))


qqplot_both
ggsave(file = sprintf('%s/indep_arbdep_qqunif.pdf', plot_savepath), width = 7, height = 4)

qqplot_both + cutoff01
ggsave(file = sprintf('%s/indep_arbdep_qqunif01.pdf', plot_savepath), width = 7, height = 4)


qqplot_all = ggplot(pvals_all2, 
                     aes(x = expected,
                         y = observed,
                         color = gamma)) +
  geom_abline(aes(intercept = 0, slope = 1)) +
  geom_point(alpha = .9, size = .8) +
  # scale_color_viridis_b(begin = 0, end = 1, breaks = c(-1, 0, .2, .4, .6), labels = c('Indep', 0, .25, .5, .75)) +
  scale_color_manual(values = c(other_colors, gamma_colors), 
                     labels = c(TeX('Proxy'), TeX('True'), TeX('Independent'), TeX('Arb Dep, \\gamma=0'), TeX('Arb Dep, \\gamma=.25'), TeX('Arb Dep, \\gamma=.50'), TeX('Arb Dep, \\gamma=.75'))) +
  labs(title = 'QQ-plot of p-values vs Unif(0,1)',
       x = 'Expected', y = 'Observed',
       color = 'pval type') +
  facet_grid(type_long~., labeller = label_parsed) +
  theme(legend.text = element_text(vjust = .7, size = 10),
        legend.title = element_text(size = 12))  +
  guides(color = guide_legend(override.aes = list(size = 4)))

qqplot_all
ggsave(file = sprintf('%s/all_qqunif.pdf', plot_savepath), width = 5, height = 5)

qqplot_all + cutoff01
ggsave(file = sprintf('%s/all_qqunif01.pdf', plot_savepath), width = 5, height = 5)

# Histogram of active pvals ====================================================
hist_both = ggplot(pvals_all, 
                     aes(x = pval_active,
                         fill = gamma)) +
  geom_histogram(aes(y = after_stat(density)), alpha = .8, size = .8, breaks = seq(0, 4, length.out = 4*10), color = 'black') +
  # scale_color_viridis_b(begin = 0, end = 1, breaks = c(-1, 0, .2, .4, .6), labels = c('Indep', 0, .25, .5, .75)) +
  scale_fill_manual(values = c('gray32', gamma_colors), 
                     labels = c(TeX('Independent'), TeX('Arb Dep, \\gamma=0'), TeX('Arb Dep, \\gamma=.25'), TeX('Arb Dep, \\gamma=.50'), TeX('Arb Dep, \\gamma=.75'))) +
  labs(title = 'Histogram of Active p-values',
       x = 'p-values', y = 'density',
       fill = 'active pval type') +
  facet_grid(type_long~gamma, labeller = label_parsed, scales = 'free') +
  theme(legend.text = element_text(vjust = .7, size = 10),
        legend.title = element_text(size = 12),
        legend.position = 'none',
        strip.text = element_text(size = 10),
        axis.text.x = element_text(size = 7))

hist_both + scale_x_continuous(limits = c(0, 1))
ggsave(file = sprintf('%s/indep_arbdep_hist01.pdf', plot_savepath), width = 7, height = 4)


hist_all = ggplot(pvals_all2, 
       aes(x = pval,
           fill = gamma)) +
  geom_histogram(aes(y = after_stat(density)), alpha = .8, size = .8, breaks = seq(0, 4, length.out = 4*10), color = 'black') +
  # scale_color_viridis_b(begin = 0, end = 1, breaks = c(-1, 0, .2, .4, .6), labels = c('Indep', 0, .25, .5, .75)) +
  scale_fill_manual(values = c(other_colors, gamma_colors), 
                    labels = c(TeX('Proxy'), TeX('True'), TeX('Independent'), TeX('Arb Dep, \\gamma=0'), TeX('Arb Dep, \\gamma=.25'), TeX('Arb Dep, \\gamma=.50'), TeX('Arb Dep, \\gamma=.75'))) +
  labs(title = 'Histogram of p-values',
       x = 'p-values', y = 'density',
       fill = 'active pval type') +
  facet_grid(type_long~gamma, labeller = label_parsed, scales = 'free') +
  theme(legend.text = element_text(vjust = .7, size = 10),
        legend.title = element_text(size = 12),
        legend.position = 'none',
        strip.text = element_text(size = 10),
        axis.text.x = element_text(size = 8))

hist_all + scale_x_continuous(limits = c(0, 1))
ggsave(file = sprintf('%s/all_hist01.pdf', plot_savepath), width = 10, height = 5)



# plot runtimes  ============================================
runtime_both = ggplot(pvals_all, 
       aes(x = time_active_expect,
           fill = gamma)) +
  geom_histogram(aes(y = after_stat(density)), 
                 color = 'gray20',  bins = 18, alpha = .7, 
                 position = 'nudge')  + 
  # scale_color_viridis_b(begin = 0, end = 1, breaks = c(-1, 0, .2, .4, .6), labels = c('Indep', 0, .25, .5, .75)) +
  scale_fill_manual(values = c('gray32', gamma_colors), 
                     labels = c(TeX('Independent'), TeX('Arb Dep, \\gamma=0'), TeX('Arb Dep, \\gamma=.25'), TeX('Arb Dep, \\gamma=.50'), TeX('Arb Dep, \\gamma=.75'))) +
  labs(title = 'Runtimes',
       x = 'time (s)',
       fill = 'active pval type') +
  facet_grid(type_long~gamma, labeller = label_parsed) +
  theme(legend.text = element_text(vjust = .7, size = 10),
        legend.title = element_text(size = 12),
        legend.position = 'none',
        strip.text = element_text(size = 9))

runtime_both
ggsave(file = sprintf('%s/indep_arbdep_runtime_hist.pdf', plot_savepath), width = 7, height = 4)


# pvals |> group_by(type) |> summarize(time_proxy_mean = mean(time_proxy),
#                                      time_proxy_sd   = sd(time_proxy, na.rm=TRUE),
#                                      time_true_mean  = mean(time_true),
#                                      time_true_sd    = sd(time_true, na.rm=TRUE))

# library(ggpattern)


ggplot(pvals_all, 
       aes(y = time_active_expect,
           # x = paste0(gamma_numeric, type_long),
           fill = gamma)) +
  # geom_boxplot_pattern(aes(pattern = type_long), pattern_fill='black', pattern_fill2='black', orientation = 'x') +
  # scale_pattern_manual(values = c('magick', 'stripe')) +
  geom_boxplot(aes(alpha = type_long), outliers = FALSE) +
  # geom_violin(aes(x = type_long, y = time_active_expect, alpha = type_long), orientation = 'x') +
  scale_alpha_manual(values = c(.5, 1)) +
  scale_fill_manual(values = c('gray32', gamma_colors), guide = 'none', 
                    labels = c(TeX('Independent'), TeX('Arb Dep, \\gamma=0'), TeX('Arb Dep, \\gamma=.25'), TeX('Arb Dep, \\gamma=.50'), TeX('Arb Dep, \\gamma=.75'))) +
  labs(title = 'Expected Runtimes',
       y = 'time (s)',
       fill = 'active pval type',
       alpha = 'test type:',
       x = '') +
  facet_grid(~gamma, labeller = label_parsed, switch = 'both') +
  theme(legend.position = c(.5,.1),
        legend.text = element_text(vjust = .7, size = 10),
        legend.direction = 'horizontal',
        legend.justification = 'center',
        legend.title = element_text(size = 12),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.text = element_text(size = 8)) +
  guides(alpha = guide_legend(override.aes = list(fill = 'gray10')))

ggsave(file = sprintf('%s/indep_arbdep_runtime_boxplot.pdf', plot_savepath), width = 7, height = 4)



ggplot(pvals_all, 
       aes(y = time_active_expect,
           # x = paste0(gamma_numeric, type_long),
           fill = gamma)) +
  # geom_boxplot_pattern(aes(pattern = type_long), pattern_fill='black', pattern_fill2='black', orientation = 'x') +
  # scale_pattern_manual(values = c('magick', 'stripe')) +
  # geom_boxplot(aes(alpha = type_long), outliers = FALSE) +
  geom_violin(aes(x = type_long, y = time_active_expect, alpha = type_long), orientation = 'x') +
  scale_alpha_manual(values = c(.5, 1)) +
  scale_fill_manual(values = c('gray32', gamma_colors), guide = 'none',
                    labels = c(TeX('Independent'), TeX('Arb Dep, \\gamma=0'), TeX('Arb Dep, \\gamma=.25'), TeX('Arb Dep, \\gamma=.50'), TeX('Arb Dep, \\gamma=.75'))) +
  labs(title = 'Expected Runtimes',
       y = 'time (s)',
       fill = 'active pval type',
       alpha = 'test type',
       x = '') +
  facet_grid(~gamma, labeller = label_parsed, switch = 'both') +
  theme(legend.position = c(.5,.1),
        legend.text = element_text(vjust = .7, size = 10),
        legend.direction = 'horizontal',
        legend.justification = 'center',
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.text = element_text(size = 8)) +
  guides(alpha = guide_legend(override.aes = list(fill = 'gray10')))

ggsave(file = sprintf('%s/indep_arbdep_runtime_violin.pdf', plot_savepath), width = 7, height = 4)


ggplot(pvals_all2, 
       aes(y = time,
           fill = gamma)) +
  geom_boxplot(aes(alpha = type_long), outliers = FALSE) +
  scale_alpha_manual(values = c(.5, 1)) +
  scale_fill_manual(values = c(other_colors, gamma_colors), guide = 'none', 
                    labels = c(TeX('Independent'), TeX('Arb Dep, \\gamma=0'), TeX('Arb Dep, \\gamma=.25'), TeX('Arb Dep, \\gamma=.50'), TeX('Arb Dep, \\gamma=.75'))) +
  labs(title = 'Expected Runtimes',
       y = 'time (s)',
       fill = 'active pval type',
       alpha = 'test type:',
       x = '') +
  facet_grid(~gamma, labeller = label_parsed, switch = 'both') +
  theme(legend.position = c(.5,.1),
        legend.text = element_text(vjust = .7, size = 10),
        legend.direction = 'horizontal',
        legend.justification = 'center',
        legend.title = element_text(size = 12),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.text = element_text(size = 8)) +
  guides(alpha = guide_legend(override.aes = list(fill = 'gray10')))

ggsave(file = sprintf('%s/all_runtime_boxplot.pdf', plot_savepath), width = 10, height = 5)



ggplot(pvals_all2, 
       aes(fill = gamma)) +
  geom_violin(aes(x = type_long, y = time, alpha = type_long), orientation = 'x') +
  scale_alpha_manual(values = c(.5, 1)) +
  scale_fill_manual(values = c(other_colors, gamma_colors), guide = 'none', 
                    labels = c(TeX('Independent'), TeX('Arb Dep, \\gamma=0'), TeX('Arb Dep, \\gamma=.25'), TeX('Arb Dep, \\gamma=.50'), TeX('Arb Dep, \\gamma=.75'))) +
  labs(title = 'Expected Runtimes',
       y = 'time (s)',
       fill = 'active pval type',
       alpha = 'test type:',
       x = '') +
  facet_grid(~gamma, labeller = label_parsed, switch = 'both') +
  theme(legend.position = c(.5,.1),
        legend.text = element_text(vjust = .7, size = 10),
        legend.direction = 'horizontal',
        legend.justification = 'center',
        legend.title = element_text(size = 12),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.text = element_text(size = 8)) +
  guides(alpha = guide_legend(override.aes = list(fill = 'gray10')))

ggsave(file = sprintf('%s/all_runtime_violin.pdf', plot_savepath), width = 10, height = 5)


