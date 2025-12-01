# ---------------------------------------------------------------------------- #
#  plot p-val from lmYA estimates
# ---------------------------------------------------------------------------- #
args = commandArgs(trailingOnly = TRUE)
# args = c('laptop', 'C_ubergenno')
# args = c('laptop', 'A')
# args = c('laptop', 'D') # lots of negative AY tests
# args = c('laptop', 'E_PCA')

require(assertthat) # for some assert statements
library(dplyr)
library(tibble)
library(ggplot2)    # plotting
library(cowplot)


theme_set(theme_cowplot() +
            theme(plot.title = element_text(hjust = .5),
                  plot.subtitle = element_text(hjust = .5),
                  strip.background = element_rect(color = 'black', fill = 'white')))

assertthat::assert_that(length(args) > 0, msg="must give arg for specifying device eg 'Rscript <filename>.R ubergenno'")



DEVICE = args[1]
source('../PATHS.R') # load in data_dir and save_dir and CODE_DIR, depending on DEVICE value
assertthat::assert_that(!is.null(data_dir), msg='first arg must be: laptop, desktop, or ubergenno')

assertthat::assert_that(length(args) > 1, msg="must give arg for specifying chosen AYZW name 'Rscript <filename>.R ubergenno C'")
AYZW_setting_name = args[2]





# =================== Start ====================================================
print(sprintf("[%s] START: lmYA plots", Sys.time()))


# load in prev saved ATE estimates and pvals
ATE_df = read.csv(sprintf('%s/cbgenes/%s/lmYA.csv', save_dir, AYZW_setting_name))
# location for saving the plots
plot_savepath = sprintf('%s/cbgenes/%s/plots_lmYA/', save_dir, AYZW_setting_name)
dir.create(plot_savepath, showWarnings = FALSE, recursive = TRUE)

# ATE_df$A |> unique() |> length()               # number of grnas
# ATE_df |> filter(type == 'positive') |> nrow() # number of causal AY tests



# Nicer Labels for type
type_labeller = labeller(type = 
                           c("negative" = "Non-Causal A-Y Pairs",
                             "maybe"    = "Candidate A-Y Pairs",
                             "positive" = "Causal A-Y Pairs"))

# =================== some plots ====================================================
print(sprintf("[%s]    - Make Plots", Sys.time()))

# =============== QQ Unif of pvals ===========================================
# original pvals
ggplot(ATE_df,
       aes(sample = pval)) +
  geom_qq(distribution = stats::qunif) +
  geom_abline(aes(slope = 1, intercept = 0)) +
  labs(title = 'QQ Unif of pvals') +
  facet_wrap(vars(type), ncol = 1,
             labeller = type_labeller)
ggsave(sprintf('%s/pval_qq.pdf', plot_savepath), height = 5, width = 5)
ggsave(sprintf('%s/pval_qq.svg', plot_savepath), height = 5, width = 5)
ggsave(sprintf('%s/pval_qq.jpg', plot_savepath), height = 5, width = 5)



# transformed pvals
transform_pval <- function(x){-log(x, base = 10)}

unifpointsref = runif(n = 10000, min = 0, max = 1)

png(sprintf('%s/pval_qq_transformed.png', plot_savepath))
qqplot(unifpointsref |> transform_pval(),
       ATE_df |> filter(type == 'negative') |> pull(pval) |> transform_pval(),
       main= 'QQplot of pvals transformed vs Unif transformed');
qqline(unifpointsref |> transform_pval(), col = 'red');
dev.off()

# # same distribution of unif transformed
# qqplot(runif(n =  1000, min = 0, max = 1) |> transform_pval(),
#        runif(n =  nrow(ATE_df), min = 0, max = 1) |> transform_pval());
# qqline(runif(n =  1000, min = 0, max = 1) |> transform_pval())


# =============== Histogram of pvals ===========================================
ggplot(ATE_df) +
  geom_histogram(aes(x = pval), binwidth = .01, fill = 'orange', alpha = .7) +
  # geom_histogram(aes(x = pval |> transform_pval())) +
  geom_vline(aes(xintercept = 0.05)) +
  labs(title = 'Histogram of pvals') +
  facet_wrap(vars(type), ncol = 1, scales = 'free_y',
             labeller = type_labeller) +
  scale_x_continuous(expand = c(0.005, 0.005),
                     breaks = seq(from=0, to=1, by=.1)) +
  scale_y_continuous(expand = c(0.005, 0.005)) 
ggsave(sprintf('%s/pval_hist.pdf', plot_savepath), height = 5, width = 5)
ggsave(sprintf('%s/pval_hist.svg', plot_savepath), height = 5, width = 5)
ggsave(sprintf('%s/pval_hist.jpg', plot_savepath), height = 5, width = 5)

# ggplot(ATE_df |> filter(type == 'negative')) +
#   geom_histogram(aes(x = pval |> transform_pval()))
# 
# 
# ggplot(ATE_df |> filter(type == 'positive')) +
#   geom_histogram(aes(x = pval))


# =============== Histogram of ATEs ===========================================
ggplot(ATE_df) +
  geom_histogram(aes(x = ATE), binwidth = .1, fill = 'orange', alpha = .7) +
  # geom_histogram(aes(x = pval |> transform_pval())) +
  geom_vline(aes(xintercept = 0)) +
  labs(title = 'Histogram of ATEs') +
  facet_wrap(vars(type), ncol = 1, scales = 'free_y',
             labeller = type_labeller) +
  scale_x_continuous(expand = c(0.005, 0.005),
                     limits = c(-5, 2.5)) +
  scale_y_continuous(expand = c(0.005, 0.005)) 
ggsave(sprintf('%s/ATE_hist.pdf', plot_savepath), height = 5, width = 5)
ggsave(sprintf('%s/ATE_hist.svg', plot_savepath), height = 5, width = 5)
ggsave(sprintf('%s/ATE_hist.jpg', plot_savepath), height = 5, width = 5)




