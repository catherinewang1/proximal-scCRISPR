# ---------------------------------------------------------------------------- #
#                   Plot Estimates of ATE (lmYA + CB)                          #
# e.g. Rscript <filename>.R ubergenno cbgenes C default .1
# 3.1 - add chromosome number information                                      #
# 3.2 - choose A,Y,Z,W combinations                                            #
# 3.3 - CB Effect Estimate                                                     #
# 3.4 - Plot ATE Estimates 
#       lmYA_threshold =  
#       lmYA_threshold_plot- for plotting only AY pairs lmYA above specified   #
#                            threshold                                         #
#       e.g. Rscript 3.4_CB_plots.R laptop .2                                  #
# Requires: prev saved ATEs in                                                 #
#    "<save_dir>/cbgenes/<AYZW_setting_name>/ATE.csv"                          #
# Ouputs: (nothing) but saves plots in                                         #
#    "<save_dir>/cbgenes/AYZW_setting_name/plots/<lmYA_threshold_plot>"        #
# ---------------------------------------------------------------------------- #
args = commandArgs(trailingOnly = TRUE)
# args = c('laptop', 'E_PCA', 'default', '.01')

# args = c('laptop', 'cbgenes_pooled', 'B_ubergenno', 'default', '.01')


# device, cbgenes or cbgenes_pooled or cbgenes_pca, AYZW setting name, CB setting name, threshold
# args = c('laptop', 'cbgenes_pooled', 'C_ubergenno', 'default', '.10')
# args = c('laptop', 'cbgenes_pca', 'E_PCA', 'simple', '.01')
# args = c('laptop', 'cbgenes', 'E_PCA', 'simple', '.1')
# args = c('laptop', 'cbgenes_cca', 'E_PCA', 'simple', '.01')
args = c('laptop', 'cbgenes_count', 'E_PCA', '.00')
args = c('laptop', 'cbgenes_count', 'A', '.00')


suppressPackageStartupMessages(require(assertthat)) # for some assert statements
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(ggplot2))    # plotting
suppressPackageStartupMessages(library(cowplot))
# suppressPackageStartupMessages(library(ggpattern))  # patterns



suppressPackageStartupMessages(require(assertthat)) # for some assert statements
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(ggplot2))    # plotting
suppressPackageStartupMessages(library(cowplot))
# suppressPackageStartupMessages(library(ggpattern))  # patterns

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

assertthat::assert_that(length(args) > 2, msg="must give arg for specifying chosen AYZW name 'Rscript <filename>.R ubergenno C'")
AYZW_setting_name = args[3]


assertthat::assert_that(length(args) > 3, msg="must give arg for specifying lmYA threshold 'Rscript <filename>.R ubergenno cbgenes C default .1'")
lmYA_threshold = as.numeric(args[4])

# =================== Start ====================================================
print(sprintf("[%s] START: CB Plots", Sys.time()))

# load chosen AYZW names
AY   = read.csv(sprintf('%s/%s/%s/AY.csv', save_dir, cbgenes_setting_name, AYZW_setting_name))
AYZW = readRDS(sprintf('%s/%s/%s/AYZW.rds', save_dir, cbgenes_setting_name, AYZW_setting_name))

# load in ATE if there, or else construct df from intermediate saved ATEs .csv
if('ATE.csv' %in% list.files(sprintf('%s/%s/%s/', save_dir, cbgenes_setting_name, AYZW_setting_name))) {
  ATE = read.csv(sprintf('%s/%s/%s/ATE.csv', save_dir, cbgenes_setting_name, AYZW_setting_name))
} else {
  ate_filenames = list.files(sprintf('%s/%s/%s/%s/intermediateATEs/', save_dir, cbgenes_setting_name, AYZW_setting_name, CB_setting_name))
  get_AY_idx <- function(x) {strsplit(x, '_|\\.')[[1]][2]}
  get_ZW_idx <- function(x) {strsplit(x, '_|\\.')[[1]][3]}
  
  ATE = NULL
  for(fn in ate_filenames) {
    ate_fn = read.csv(sprintf('%s/%s/%s/%s/intermediateATEs/%s', 
                              save_dir, cbgenes_setting_name, AYZW_setting_name, CB_setting_name, fn))
    ate_fn$AY_idx = get_AY_idx(fn)
    ate_fn$ZW_idx = get_ZW_idx(fn)
    ATE = rbind(ATE, ate_fn)
  }
}

ATE$EstMethod = ATE$method

# long EstMethod name (method, basis, numNC)
ATE$EstMethod_l = apply(cbind(ATE$EstMethod, ATE$numNC), 1,                  # EstMethod col is long name
                        function(x) paste(x[!is.na(x)], collapse = ""))

head(ATE); dim(ATE)
nrow(ATE |> distinct(AY_idx))
nrow(ATE |> distinct(AY_idx, ZW_idx))

AY$AY_idx = 1:nrow(AY)
df = merge(AY,
           ATE,
           by = 'AY_idx')

# df = reshape2::melt(df, id = c('AY_idx', 'ZW_idx', 'type', 'A', 'A_chr', 'Y', 'Y_chr'),
#                 value.name = 'ATE') 
# df = df |> rename(EstMethod = variable)


# convert AY type to factor with order
df$type = factor(df$type, levels = c('maybe', 'negative', 'positive'))

# AY_idx with |lmYA estimate| > lmYA_threshold
AY_idx_threshold = df |> filter(EstMethod == 'glmYA' & 
                                  abs(ATE) >= lmYA_threshold
                                  # ATE <= -1 * abs(lmYA_threshold) # estimating too negative
                                ) |>
  pull(AY_idx) |> unique()

df = df |> filter(AY_idx %in% AY_idx_threshold)
# dim(df); head(df)
# df |> group_by(type, EstMethod) |> summarize(ATE = mean(ATE, na.rm = T))

# Nicer Labels for type
type_labeller = labeller(type = 
                           c("negative" = "Non-Causal A-Y Pairs",
                             "maybe"    = "Candidate A-Y Pairs",
                             "positive" = "Causal A-Y Pairs"))

# order methods for display
xorder = c("glmYAX", "glmYA", "Pois2SLS")

xorder = xorder[xorder %in% unique(df$EstMethod)]

# order methods for display (more detail, includes numNC)
xorder_l = c("glmYAX", "glmYA")
for(twoSLSname in c('Pois2SLS')) {
  for(numNC in sort(unique(df$numNC))) {
    xorder_l = c(xorder_l,
                 paste0(twoSLSname, ""   , numNC))
  }
}
# for(m in c('OCBGMM', 'OCBGMMRW', 'OCBGMMRWReg', 
#            'OCBLinPI', 'OCBLinOS', 'OCBLinOStrim')) {
#   for(b in c(1, 2)) {
#     for(numNC in sort(unique(df$numNC))) {
#       xorder_l = c(xorder_l, paste0(m, 'basis', b, numNC))
#     }
#   }
# }

xorder_l = xorder_l[xorder_l %in% unique(df$EstMethod_l)]


xorder_l_breaks = c('glmYAX', 'glmYA', 'Pois2SLS')
xorder_l_breaks = xorder_l_breaks[xorder_l_breaks %in% unique(df$EstMethod_l)]

ylim_ATE = c(-3, 3)
ylim_FC  = c(0, 2)

# # lmYA estimates and p-values only (saved in cbgenes folder)
# df_lmYA = tryCatch({read.csv(sprintf('%s/cbgenes/%s/lmYA.csv', 
#                                      save_dir, AYZW_setting_name))},
#                    error = function(cond) { return(NULL) })


# =================== Make Plots ======================================
print(sprintf("[%s]    - Make Plots", Sys.time()))

plot_savepath = sprintf('%s/%s/%s/plots/lmYAthresh=%05.02f', save_dir, cbgenes_setting_name, AYZW_setting_name, lmYA_threshold)
dir.create(plot_savepath, showWarnings = FALSE, recursive = TRUE)


# =================== Make Plots: all ATEs ======================================
print(sprintf("[%s]        - all ATEs", Sys.time()))

# all ATE Estimates: Histogram by type
ggplot(df,
       aes(x = ATE, y = after_stat(density))) +
  geom_histogram(fill = 'orange', alpha = .6, binwidth = .1) +
  geom_density(color = 'orange') +
  geom_vline(aes(xintercept = 0), color = 'darkorange') +
  scale_x_continuous(limits = c(-4, 4)) +
  facet_wrap(vars(type), nrow = 3, labeller = type_labeller) +
  # facet_grid(type ~ . , labeller = type_labeller) +
  labs(title = 'Histogram of All ATE Estimates',  
       x = 'ATE Estimate') 
ggsave(sprintf('%s/ATE_hist.pdf', plot_savepath), height = 5, width = 5)
ggsave(sprintf('%s/ATE_hist.svg', plot_savepath), height = 5, width = 5)
ggsave(sprintf('%s/ATE_hist.jpg', plot_savepath), height = 5, width = 5)


# all ATE Estimates: Boxplot by type
# ggplot(df[sample(1:nrow(df), size = 1000, ), ]) +
ggplot(df) +
  geom_hline(aes(yintercept = 0), alpha = 1, linetype = 'solid') +
  geom_boxplot(aes(x = EstMethod_l, y = ATE, 
                   # pattern_fill  = basis,
                   # alpha = basis, 
                   fill = method),
               # fill = 'orange', 
               # alpha = .6,
               linewidth = .1,
               outlier.shape = NA, # don't show outliers
               outlier.size = .1, outlier.alpha = .3) +
  coord_cartesian(ylim = c(-2.25, 1.75)) +
  scale_x_discrete(limits = xorder_l) +
  # scale_alpha_discrete(guide = 'none') +
  facet_wrap(vars(type), nrow = 3) +
  labs(title = 'Boxplots of All ATE Estimates',
       x = 'Estimation Method', y = 'ATE Estimate') +
  theme(axis.text.x = element_text(angle = 45, vjust = 1.2, hjust = 1.3, size = 10))
ggsave(sprintf('%s/ATE_boxplot.pdf', plot_savepath), height = 5, width = 8)
ggsave(sprintf('%s/ATE_boxplot.svg', plot_savepath), height = 5, width = 8)
ggsave(sprintf('%s/ATE_boxplot.jpg', plot_savepath), height = 5, width = 8)



# ATE Violin Plot by Estimation Method
ggplot(df,
       aes(x = EstMethod_l, y = ATE, 
           fill = method, color = method)) +
  geom_hline(aes(yintercept = 0), alpha = .8, linetype = 'solid', color = 'red') +
  geom_jitter(shape=16, position=position_jitter(0.2),
              alpha = .6, 
              # color = 'darkorange2', 
              size = 1) +
  geom_violin(alpha = .3, 
              color = 'black',
              # fill = 'orange',
              draw_quantiles = c(.05, .95),
              trim = TRUE,
              # bounds = c(0, 1),
              scale = TRUE) +
  # geom_boxplot(# fill = 'orange', 
  #              # alpha = .6,
  #              linewidth = .1,
  #              outlier.shape = NA, # don't show outliers
  #              outlier.size = .1, outlier.alpha = .3) +
  # coord_cartesian(ylim = ylim*2) +
  # scale_x_discrete(limits = c('lmYA0', xorder_l), #) + #,
  #                  breaks = c('lmYA0', xorder_l_breaks)) +
  # scale_alpha_discrete(guide = 'none') +
  scale_y_continuous(limits = ylim_ATE) +
  facet_wrap(vars(type), nrow = 3, labeller = type_labeller) +
  labs(title = 'Violin plot ATE Estimates by Estimation Methods',
       x = 'Estimation Method', y = 'ATE Estimate') +
  theme(axis.text.x = element_text(angle = 45, vjust = 1.2, hjust = 1.3, size = 10))
ggsave(sprintf('%s/ATE_violinplot.pdf', plot_savepath), height = 10, width = 20)
ggsave(sprintf('%s/ATE_violinplot.svg', plot_savepath), height = 10, width = 20)
ggsave(sprintf('%s/ATE_violinplot.jpg', plot_savepath), height = 10, width = 20)





# FC Violin Plot by Estimation Method
ggplot(df,
       aes(x = factor(EstMethod_l, levels = xorder_l), y = FC, 
           fill = method, color = method)) +
  geom_hline(aes(yintercept = 1), alpha = .8, linetype = 'solid', color = '#EE2C2C') +
  geom_jitter(shape=16, position=position_jitter(0.2),
              alpha = .6, 
              size = 1) +
  geom_violin(alpha = .3, 
              color = 'black',
              draw_quantiles = c(.05, .95),
              trim = FALSE,
              # bounds = c(0, 1),
              scale = TRUE) +
  scale_x_discrete(breaks = xorder_l) +
  scale_y_continuous(limits = ylim_FC) +
  facet_wrap(vars(type), nrow = 3, labeller = type_labeller) +
  labs(title = 'Violin plot FC Estimates by Estimation Methods',
       x = 'Estimation Method', y = 'FC Estimate') +
  theme(axis.text.x = element_text(angle = 45, vjust = 1.2, hjust = 1.3, size = 10))
ggsave(sprintf('%s/ATE_violinplot.pdf', plot_savepath), height = 10, width = 20)
ggsave(sprintf('%s/ATE_violinplot.svg', plot_savepath), height = 10, width = 20)
ggsave(sprintf('%s/ATE_violinplot.jpg', plot_savepath), height = 10, width = 20)

# log FC Violin Plot by Estimation Method
ggplot(df,
       aes(x = factor(EstMethod_l, levels = xorder_l), y = log(FC), 
           fill = method, color = method)) +
  geom_hline(aes(yintercept = 0), alpha = .8, linetype = 'solid', color = '#EE2C2C') +
  geom_jitter(shape=16, position=position_jitter(0.2),
              alpha = .6, 
              size = 1) +
  geom_violin(alpha = .3, 
              color = 'black',
              draw_quantiles = c(.05, .95),
              trim = FALSE,
              # bounds = c(0, 1),
              scale = TRUE) +
  scale_x_discrete(breaks = xorder_l) +
  scale_y_continuous(limits = log(ylim_FC + .1)) +
  facet_wrap(vars(type), nrow = 3, labeller = type_labeller) +
  labs(title = 'Violin plot FC Estimates by Estimation Methods',
       x = 'Estimation Method', y = 'FC Estimate') +
  theme(axis.text.x = element_text(angle = 45, vjust = 1.2, hjust = 1.3, size = 10))
ggsave(sprintf('%s/ATE_violinplot.pdf', plot_savepath), height = 10, width = 20)
ggsave(sprintf('%s/ATE_violinplot.svg', plot_savepath), height = 10, width = 20)
ggsave(sprintf('%s/ATE_violinplot.jpg', plot_savepath), height = 10, width = 20)






# pval Violin Plot by Estimation Method
ggplot(df,
       aes(x = factor(EstMethod_l, levels = xorder_l), y = pval, 
           fill = method, color = method)) +
  geom_hline(aes(yintercept = 0.05), alpha = .8, linetype = 'solid', color = '#EE2C2C') +
  geom_jitter(shape=16, position=position_jitter(0.2),
              alpha = .6, 
              size = 1) +
  geom_violin(alpha = .3, 
              color = 'black',
              draw_quantiles = c(.05, .95),
              trim = FALSE,
              # bounds = c(0, 1),
              scale = TRUE) +
  scale_x_discrete(breaks = xorder_l) +
  scale_y_continuous(limits = c(0,1)) +
  facet_wrap(vars(type), nrow = 3, labeller = type_labeller) +
  labs(title = 'Violin plot FC Estimates by Estimation Methods',
       x = 'Estimation Method', y = 'FC Estimate') +
  theme(axis.text.x = element_text(angle = 45, vjust = 1.2, hjust = 1.3, size = 10))
ggsave(sprintf('%s/ATE_violinplot.pdf', plot_savepath), height = 10, width = 20)
ggsave(sprintf('%s/ATE_violinplot.svg', plot_savepath), height = 10, width = 20)
ggsave(sprintf('%s/ATE_violinplot.jpg', plot_savepath), height = 10, width = 20)




# FAKE P-VALUE CONSTRUCTION (COMPARE TO EMPIRICAL DISTN)
# =================== Make Plots: looking at pvalues ======================================
# folder for p-value plots to be saved in 
dir.create(sprintf('%s/pvals/', plot_savepath), showWarnings = FALSE)

# ---- p-values from z-tests on each AY pair means of ZW choices, estimate variance from
# all means for particular estimation method setting (e.g. OCB2SLSReg with 5 NCs)

# using estimated standard deviation of the estimates for non-causal AY for each EstMethod_l
# (Approx with a normal... but could just use the quantiles... but then REALLY need to use separate sample)
EstMethod_sd = df |> 
  filter(type == 'negative') |>
  # filter(EstMethod_l != 'lmYA') |>
  group_by(AY_idx, EstMethod_l) |> # first find ATEmean for each AY, EstMethod_l 
  summarize(ATEmean = mean(ATE),
            .groups = 'drop') |>
  group_by(EstMethod_l) |> 
  summarize(EstMethod_l_sd = sd(ATEmean))

# dim(EstMethod_sd); head(EstMethod_sd)

# should be faster access to sd than filtering
EstMethod_sd_list = setNames(EstMethod_sd$EstMethod_l_sd, EstMethod_sd$EstMethod_l)
mypvalfunc <- function(x, m)  {
  # m =  'OCB2SLS10'; x = -.02
  # pnorm(0, mean = 0, sd = .03)
  pnorm(q=x, 
        mean = 0, 
        sd = EstMethod_sd_list[m])
  # sd = EstMethod_sd |> filter(EstMethod_l == m) |> pull(EstMethod_l_sd))
}

df_ATEmean = df |> 
  group_by(AY_idx, type, numNCO, method, EstMethod_l) |>
  summarize(ATEmean = mean(ATE),
            .groups = 'drop')

# dim(df_ATEmean); head(df_ATEmean)
# mapply(FUN = mypvalfunc, 
#        x = c(-1, 0, .2),
#        m = c('OCB2SLS10', 'OCB2SLS10', 'OCB2SLS50'))
df_ATEmean$zpvals = mapply(FUN = mypvalfunc, 
                           x = df_ATEmean$ATEmean,
                           m = df_ATEmean$EstMethod_l)

# head(df_ATEmean)
# df_ATEmean |> filter(is.na(zpvals))
# ggplot(df_ATEmean, 
#        aes(x = zpvals)) +
#   geom_histogram() +
#   facet_wrap(vars(type), ncol = 1)
# df_ATEmean$zpvals |> hist()


# # z pvalues violin plot per EstMethod
# ggplot(df_ATEmean |> filter(numNC > 0 | is.na(numNC)),
#        aes(x = EstMethod_l, y = zpvals, 
#            fill = method, color = method)) +
#   geom_hline(aes(yintercept = 0.05), alpha = .8, linetype = 'solid', color = 'red') +
#   geom_jitter(shape=16, position=position_jitter(0.2),
#               alpha = .6, 
#               # color = 'darkorange2', 
#               size = 1) +
#   geom_violin(alpha = .3, 
#               color = 'black',
#               # fill = 'orange',
#               draw_quantiles = c(.05),
#               trim = TRUE,
#               bounds = c(0, 1),
#               scale = TRUE) +
#   # geom_boxplot(# fill = 'orange', 
#   #              # alpha = .6,
#   #              linewidth = .1,
#   #              outlier.shape = NA, # don't show outliers
#   #              outlier.size = .1, outlier.alpha = .3) +
#   # coord_cartesian(ylim = ylim*2) +
#   scale_x_discrete(limits = xorder_l,
#                    breaks = xorder_l_breaks) +
#   # scale_alpha_discrete(guide = 'none') +
#   facet_wrap(vars(type), nrow = 3) +
#   labs(title = 'Violin plot of z-test pvals (sd est neg AY pairs)',
#        x = 'Estimation Method', y = 'pval') +
#   theme(axis.text.x = element_text(angle = 45, vjust = 1.2, hjust = 1.3, size = 10))
# ggsave(sprintf('%s/pvals/zpvals_violinplot.pdf', plot_savepath), height = 10, width = 16)
# ggsave(sprintf('%s/pvals/zpvals_violinplot.svg', plot_savepath), height = 10, width = 16)
# ggsave(sprintf('%s/pvals/zpvals_violinplot.jpg', plot_savepath), height = 10, width = 16)



# Include the lmYA pvalues (pvalues from linear model of Y ~ A)
df_lmYA$numNCO = NA
df_lmYA$method = 'lmYA0'
df_lmYA$EstMethod_l = 'lmYA0'


df_pvals = dplyr::bind_rows(df_ATEmean |> dplyr::rename(pval = zpvals), 
                            df_lmYA |> dplyr::select(-A, -A_chr, -Y, -Y_chr, -ZW_idx, -ATE))
head(df_pvals)

# z pvalues violin plot per EstMethod incl lmYA pvals 
ggplot(df_pvals |> filter(numNC > 0 | is.na(numNC)),
       aes(x = EstMethod_l, y = pval, 
           fill = method, color = method)) +
  geom_hline(aes(yintercept = 0.05), alpha = .8, linetype = 'solid', color = 'red') +
  geom_jitter(shape=16, position=position_jitter(0.2),
              alpha = .6, 
              # color = 'darkorange2', 
              size = 1) +
  geom_violin(alpha = .3, 
              color = 'black',
              # fill = 'orange',
              draw_quantiles = c(.05),
              trim = TRUE,
              bounds = c(0, 1),
              scale = TRUE) +
  # geom_boxplot(# fill = 'orange', 
  #              # alpha = .6,
  #              linewidth = .1,
  #              outlier.shape = NA, # don't show outliers
  #              outlier.size = .1, outlier.alpha = .3) +
  # coord_cartesian(ylim = ylim*2) +
  scale_x_discrete(limits = c('lmYA0', xorder_l), #) + #,
                   breaks = c('lmYA0', xorder_l_breaks)) +
  # scale_alpha_discrete(guide = 'none') +
  facet_wrap(vars(type), nrow = 3) +
  labs(title = 'Violin plot of lmYA pval/z-test pvals (sd est neg AY pairs)',
       x = 'Estimation Method', y = 'pval') +
  theme(axis.text.x = element_text(angle = 45, vjust = 1.2, hjust = 1.3, size = 10))
ggsave(sprintf('%s/pvals/zpvals_violinplot.pdf', plot_savepath), height = 10, width = 20)
ggsave(sprintf('%s/pvals/zpvals_violinplot.svg', plot_savepath), height = 10, width = 20)
ggsave(sprintf('%s/pvals/zpvals_violinplot.jpg', plot_savepath), height = 10, width = 20)


# ---- p-values from t-tests on each AY pair over ZW choices 
# e.g. for A5 --> Y4, there were 5 sets of ZW choices. 
# Perform t-test on these 5 estimates, to test against ATE=0
# No, the estimates for each AY is too correlated, so underestimates the variance
# simple way to get a quick 'p-value'? 
# t-test
myttestfunc <- function(x) {
  if(length(x) <= 1 | all(is.na(x))) {
    return(NA)
  }
  t.test(x, mu=0)$p.value
}
# res = t.test(3:10, mu=0)

ttest_pval = df |> 
  # filter(type == 'negative') |>
  filter(EstMethod_l != 'lmYA') |>
  group_by(AY_idx, EstMethod_l) |>
  summarize(pval = myttestfunc(ATE),
            ATEmean = mean(ATE),
            ATEmedian = median(ATE),
            ATEsd = sd(ATE),
            count = n(),
            .groups = 'drop')


# head(ttest_pval)

df_tpval = merge(df |> select(AY_idx, type, A, A_chr, Y, Y_chr, 
                              method, method_type, numNC, basis, 
                              EstMethod, EstMethod_l) |> distinct(),
                 ttest_pval)
# add lmYA0 pvalues (from linear model)
df_tpval = dplyr::bind_rows(df_tpval,
                            df_lmYA |> dplyr::select(-A, -A_chr, -Y, -Y_chr, -ZW_idx, -ATE))
# head(df_tpval)
# dim(df_tpval); dim(df)
# df |> filter(AY_idx == 2 & EstMethod_l == 'OCB2SLS1')
# hist(ttest_pval$pval)

ggplot(df_tpval, 
       aes(x = EstMethod_l, y = pval,
           fill = method, color = method)) +
  geom_hline(aes(yintercept = 0), alpha = 1, linetype = 'solid') +
  geom_jitter(shape=16, position=position_jitter(0.2),
              alpha = .6, 
              # color = 'darkorange2', 
              size = 1) +
  geom_violin(alpha = .3, 
              color = 'black',
              # fill = 'orange',
              # draw_quantiles = c(.05),
              trim = TRUE,
              # bounds = c(0, 1),
              scale = TRUE) +
  # geom_boxplot(aes( 
  #                  # pattern_fill  = basis,
  #                  alpha = basis, 
  #                  fill = method),
  #              # fill = 'orange', 
  #              # alpha = .6,
  #              linewidth = .1,
  #              # outlier.shape = NA, # don't show outliers
  #              outlier.size = 1, outlier.alpha = 1) +
  # coord_cartesian(ylim = ylim*2) +
  # scale_y_log10() +
scale_x_discrete(limits = c('lmYA0', xorder_l), #) + #,
                 breaks = c('lmYA0', xorder_l_breaks)) +
  # scale_y_continuous(limits = c(0, 1)) +
  # scale_alpha_discrete(guide = 'none') +
  facet_wrap(vars(type), nrow = 3) +
  labs(title = 'Violin plot of lmYA pval/t-test pvals',
       x = 'Estimation Method', y = 'pval') +
  theme(axis.text.x = element_text(angle = 45, vjust = 1.2, hjust = 1.3, size = 10))
ggsave(sprintf('%s/pvals/tpvals_violinplot.pdf', plot_savepath), height = 10, width = 20)
ggsave(sprintf('%s/pvals/tpvals_violinplot.svg', plot_savepath), height = 10, width = 20)
ggsave(sprintf('%s/pvals/tpvals_violinplot.jpg', plot_savepath), height = 10, width = 20)







