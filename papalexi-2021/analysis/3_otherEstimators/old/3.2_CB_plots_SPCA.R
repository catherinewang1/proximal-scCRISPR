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
# args = c('laptop', 'spca/cbgenes', 'allPerturbations', 'simple', '0')
args = c('laptop', 'spca/cbgenes', 'allPerturbations', 'simple', 'NA')
# args = c('laptop', 'cbgenes_pooled', 'B_ubergenno', 'default', '.01')


# device, cbgenes or cbgenes_pooled or cbgenes_pca, AYZW setting name, CB setting name, threshold
# args = c('laptop', 'cbgenes_pooled', 'C_ubergenno', 'default', '.10')
# args = c('laptop', 'cbgenes_pca', 'E_PCA', 'simple', '.01')
# args = c('laptop', 'cbgenes', 'E_PCA', 'simple', '.1')
# args = c('laptop', 'cbgenes_cca', 'E_PCA', 'simple', '.01')
# args = c('laptop', 'cbgenes_cca', 'E_PCA', 'simple', '.00')


suppressPackageStartupMessages(require(assertthat)) # for some assert statements
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(ggplot2))    # plotting
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(RColorBrewer))
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


# load in ATE if there, or else construct df from intermediate saved ATEs .csv
if('ATE.csv' %in% list.files(sprintf('%s/%s/%s/%s/', save_dir, cbgenes_setting_name, AYZW_setting_name, CB_setting_name))) {
  ATE = read.csv(sprintf('%s/%s/%s/%s/ATE.csv', save_dir, cbgenes_setting_name, AYZW_setting_name, CB_setting_name))
} else {
  ate_filenames = list.files(sprintf('%s/%s/%s/%s/intermediateATEs/', save_dir, cbgenes_setting_name, AYZW_setting_name, CB_setting_name))
  if(length(ate_filenames) < nrow(AY)) {print('   CAUTION: fewer intermedate ATEs than requested AY pairs')}
  
  get_AY_idx <- function(x) {strsplit(x, '_|\\.')[[1]][2]}
  if(!grepl('spca', cbgenes_setting_name)) {get_ZW_idx <- function(x) {strsplit(x, '_|\\.')[[1]][3]}}
  
  
  ATE = NULL
  for(fn in ate_filenames) {
    ate_fn = read.csv(sprintf('%s/%s/%s/%s/intermediateATEs/%s', 
                              save_dir, cbgenes_setting_name, AYZW_setting_name, CB_setting_name, fn))
    ate_fn$AY_idx = get_AY_idx(fn)
    if(!grepl('spca', cbgenes_setting_name)) {ate_fn$ZW_idx = get_ZW_idx(fn)}
    
    ATE = rbind(ATE, ate_fn)
  }
  
  write.csv(ATE, sprintf('%s/%s/%s/%s/ATE.csv', save_dir, cbgenes_setting_name, AYZW_setting_name, CB_setting_name))
}




ATE$EstMethod = apply(cbind(ATE$method, ATE$basis), 1,                  # EstMethod col is long name
                      function(x) paste(x[!is.na(x)], collapse = ""))
# ATE = ATE |> mutate(EstMethod = paste0(method, '_', 'basis', basis))

# long EstMethod name (method, basis, numNC)
ATE$EstMethod_l = apply(cbind(ATE$EstMethod, ATE$numNC), 1,                  # EstMethod col is long name
      function(x) paste(x[!is.na(x)], collapse = ""))

head(ATE); dim(ATE)
nrow(ATE |> distinct(AY_idx))
# nrow(ATE |> distinct(AY_idx, ZW_idx))

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
if(!is.na(lmYA_threshold)) {
  AY_idx_threshold = df |> filter(EstMethod == 'lmYA' & 
                                    # abs(ATE) >= lmYA_threshold
                                    ATE <= -1 * abs(lmYA_threshold) # estimating too negative
  ) |>
    pull(AY_idx) |> unique()
  
  df = df |> filter(AY_idx %in% AY_idx_threshold) 
}

# dim(df); head(df)
# df |> group_by(type, EstMethod) |> summarize(ATE = mean(ATE, na.rm = T))

# Nicer Labels for type
type_map = c("negative" = "Non-Causal A-Y",
             "maybe"    = "Candidate A-Y",
             "positive" = "Causal A-Y")
type_labeller = labeller(type = type_map)

EstMethod_map = c("lmYA" = "Naive",
               "OCB2SLS"    = "2SLS",
               "OCB2SLSpci2s" = "2SLS (pci2s)",
               "OCB2SLSReg"    = "2SLS (Reg)",
               "OCBMEstbasis1"    = "M-Est",
               "OCBMEstRWbasis1"    = "M-Est (RW)",
               "OCBMEstRWRegbasis1"    = "M-Est (RW & Reg)")

EstMethod_labeller = labeller(type =  EstMethod_map)

# order methods for display
xorder = c(# "Oracle", 
           "lmYA", "OCB2SLSpci2s", "OCB2SLS", "OCB2SLSReg")
for(m in c('OCBMEst', 'OCBMEstRW', 'OCBMEstRWReg', 
           'OCBLinPI', 'OCBLinOS', 'OCBLinOStrim')) {
  for(b in c(1, 2)) {
    xorder = c(xorder, paste0(m, 'basis', b))
  }
}
xorder = xorder[xorder %in% unique(df$EstMethod)]

# order methods for display (more detail, includes numNC)
xorder_l = c(# "Oracle", 
             "lmYA")
for(twoSLSname in c('OCB2SLS', "OCB2SLSpci2s", 'OCB2SLSReg')) {
  for(numNC in sort(unique(df$numNC))) {
    xorder_l = c(xorder_l,
                 paste0(twoSLSname, ""   , numNC))
  }
}
for(m in c('OCBMEst', 'OCBMEstRW', 'OCBMEstRWReg', 
           'OCBLinPI', 'OCBLinOS', 'OCBLinOStrim')) {
  for(b in c(1, 2)) {
    for(numNC in sort(unique(df$numNC))) {
      xorder_l = c(xorder_l, paste0(m, 'basis', b, numNC))
    }
  }
}

xorder_l = xorder_l[xorder_l %in% unique(df$EstMethod_l)]


xorder_l_breaks = c('lmYA', 'OCB2SLS1', 'OCB2SLSReg1', 
                    'OCBMEstbasis11', 'OCBMEstbasis21', 
                    'OCBMEstRWbasis11', 'OCBMEstRWbasis21', 
                    'OCBMEstRWRegbasis11', 'OCBMEstRWRegbasis21', 
                    'OCBLinPIbasis11', 'OCBLinPIbasis21')
xorder_l_breaks = xorder_l_breaks[xorder_l_breaks %in% unique(df$EstMethod_l)]

ylim = c(-1, .5)


# lmYA estimates and p-values only (saved in cbgenes folder)
# df_lmYA = tryCatch({read.csv(sprintf('%s/cbgenes/%s/lmYA.csv', 
#                                      save_dir, AYZW_setting_name))},
#                   error = function(cond) { return(NULL) })

# !! Removing TF AY pairs and cleaning up for Active p-values!! 
if(F) {
  # get Y TFs
  df3 = merge(df, gene_TF_Target |> mutate(TFY = TF) |> select(gene, TFY), 
              by.x = 'Y', by.y = 'gene',
              all.x = TRUE, all.y = FALSE)
  
  # get A TFs
  # assuming that A's target is all but the last to chars
  df3$A_target = sapply(df3$A, FUN = function(x) {substring(text=x, first=1, last=nchar(x)-2)}) 
  df4 = merge(df3, gene_TF_Target |> mutate(TFA = TF) |> select(gene, TFA), 
              by.x = 'A_target', by.y = 'gene',
              all.x = TRUE, all.y = FALSE)
  
  # checking
  df4$TFY |> table() # Y is TF
  df4$TFA |> table(); sum(is.na(df4$TFA)) # A is TF
  
  
  df5 = df4 |> filter((!TFY) & (is.na(TFA) | (!TFA))) # remove TFs
  df5 = df5 |> filter(type != 'maybe') # remove maybe types
  
  
  # cleaning up to save for active p-values
  rename_to_active <- function(x) {
    if(x == 'OCB2SLSpci2s20'){return('true')}
    else if(x == 'lmYA'){return('proxy')} 
    else{return('')}
    
  }
  df6 = df5 |> select(AY_idx, type, A, Y, EstMethod, EstMethod_l, pval) |> 
    filter(EstMethod_l == 'OCB2SLSpci2s20' | EstMethod_l == 'lmYA') |>
    mutate(pval_type = mapply(rename_to_active, EstMethod_l)) |>
    arrange(AY_idx, pval_type) |>
    select(-EstMethod_l) |>
    select(AY_idx, type, pval_type, EstMethod, pval)
  write.csv(df6, 
            sprintf('%s/%s/%s/pvalues/pvals_papalexi_noTFs.csv', save_dir, cbgenes_setting_name, AYZW_setting_name),
            row.names = FALSE)
  # checking
  df6 |> group_by(type, pval_type) |> summarize(count = n())
}




# !! Removing Candidate/Maybe AY pairs !!
df = df |> filter(type != 'maybe')

# =================== Make Plots ======================================
print(sprintf("[%s]    - Make Plots", Sys.time()))

plot_savepath = sprintf('%s/%s/%s/%s/plots/lmYAthresh=%05.02f', save_dir, cbgenes_setting_name, AYZW_setting_name, CB_setting_name, lmYA_threshold)
dir.create(plot_savepath, showWarnings = FALSE, recursive = TRUE)


# =================== Make Plots: all ATEs ======================================
print(sprintf("[%s]        - all ATEs", Sys.time()))

# # all ATE Estimates: Histogram by type
# ggplot(df,
#        aes(x = ATE, y = after_stat(density))) +
#   geom_histogram(fill = 'orange', alpha = .6, binwidth = .1) +
#   geom_density(color = 'orange') +
#   geom_vline(aes(xintercept = 0), color = 'darkorange') +
#   scale_x_continuous(limits = c(-4, 4)) +
#   facet_wrap(vars(type), nrow = 3, labeller = type_labeller) +
#   # facet_grid(type ~ . , labeller = type_labeller) +
#   labs(title = 'Histogram of All ATE Estimates',  
#        x = 'ATE Estimate') 
# ggsave(sprintf('%s/ATE_hist.pdf', plot_savepath), height = 5, width = 5)
# ggsave(sprintf('%s/ATE_hist.svg', plot_savepath), height = 5, width = 5)
# ggsave(sprintf('%s/ATE_hist.jpg', plot_savepath), height = 5, width = 5)


# all ATE Estimates: Boxplot by type
# ggplot(df[sample(1:nrow(df), size = 1000, ), ]) +
# ggplot(df) +
#   geom_hline(aes(yintercept = 0), alpha = 1, linetype = 'solid') +
#   geom_boxplot(aes(x = EstMethod, y = ATE, 
#                    # pattern_fill  = basis,
#                    # alpha = basis, 
#                    fill = method),
#                  # fill = 'orange', 
#                  linewidth = .1,
#                  outlier.shape = NA, # don't show outliers
#                  outlier.size = .1, outlier.alpha = .3, alpha = .8) +
#   coord_cartesian(ylim = ylim) +
#   scale_x_discrete(limits = xorder) +
#   scale_alpha_discrete(guide = 'none') +
#   facet_wrap(vars(type), nrow = 3) +
#   labs(title = 'Boxplots of All ATE Estimates',
#        x = 'Estimation Method', y = 'ATE Estimate') +
#   theme(axis.text.x = element_text(angle = 45, vjust = 1.2, hjust = 1.3, size = 10))
# ggsave(sprintf('%s/ATE_boxplot.pdf', plot_savepath), height = 5, width = 8)
# ggsave(sprintf('%s/ATE_boxplot.svg', plot_savepath), height = 5, width = 8)
# ggsave(sprintf('%s/ATE_boxplot.jpg', plot_savepath), height = 5, width = 8)

 
# =================== Make Plots: AY pairs' ATEs ======================================
# print(sprintf("[%s]        - AY pairs' ATEs", Sys.time()))
# 
# # AY ATE Estimates (average over ZW selections): Boxplot by type
# # mean
# ggplot(df |> group_by(type, AY_idx, EstMethod, method, basis) |> 
#          summarize(ATE = mean(ATE, na.rm = TRUE))) +
#   geom_hline(aes(yintercept = 0), alpha = 1, linetype = 'solid') +
#   geom_boxplot(aes(x = EstMethod, y = (ATE), alpha = basis, fill = method),
#                # fill = 'orange', alpha = .6, 
#                outlier.size = .1, outlier.alpha = .3) +
#   coord_cartesian(ylim = ylim) +
#   scale_x_discrete(limits = xorder) +
#   scale_alpha_discrete(guide = 'none') +
#   facet_wrap(vars(type), nrow = 3,
#              labeller = type_labeller) +
#   labs(title = 'Mean ATE for AY pairs',
#        y = 'Mean ATE', x = 'Estimation Method') +
#   theme(axis.text.x = element_text(angle = 45, vjust = 1.2, hjust = 1.1, size = 10))
# ggsave(sprintf('%s/AY_meanATE_boxplot.pdf', plot_savepath), height = 5, width = 8)
# ggsave(sprintf('%s/AY_meanATE_boxplot.svg', plot_savepath), height = 5, width = 8)
# ggsave(sprintf('%s/AY_meanATE_boxplot.jpg', plot_savepath), height = 5, width = 8)
# 
# # median
# ggplot(df |> group_by(type, AY_idx, EstMethod, method, basis) |> 
#          summarize(ATE = median(ATE, na.rm = TRUE))) +
#   geom_hline(aes(yintercept = 0), alpha = 1, linetype = 'solid', color = 'gray') +
#   geom_boxplot(aes(x = EstMethod, y = (ATE), 
#                    alpha = basis, fill = method),
#                # fill = 'orange', alpha = .6, 
#                outlier.size = .1, outlier.alpha = .3) +
#   coord_cartesian(ylim = ylim) +
#   scale_x_discrete(limits = xorder) +
#   scale_alpha_discrete(guide = 'none') +
#   facet_wrap(vars(type), nrow = 3,
#              labeller = type_labeller) +
#   labs(title = 'Median ATE for AY pairs',
#        y = 'Median ATE', x = 'Estimation Method') +
#   theme(axis.text.x = element_text(angle = 45, vjust = 1.2, hjust = 1.1, size = 10))
# ggsave(sprintf('%s/AY_medianATE_boxplot.pdf', plot_savepath), height = 5, width = 8)
# ggsave(sprintf('%s/AY_medianATE_boxplot.svg', plot_savepath), height = 5, width = 8)
# ggsave(sprintf('%s/AY_medianATE_boxplot.jpg', plot_savepath), height = 5, width = 8)
# 
# 
# # AY ATE Estimates (average over ZW selections): Lines 
# # median
# ggplot(df |> filter(type == 'positive' | type == 'negative') |> 
#     group_by(type, AY_idx, EstMethod, method, basis) |> 
#     mutate(AY_idx = as.factor(AY_idx)) |>
#     summarize(ATE_median = median(ATE),
#               .groups = 'drop'),
#   aes(x = EstMethod,  y = ATE_median,
#       color = type,
#       fill = type,
#       group = AY_idx)) +
#   geom_hline(aes(yintercept = 0), alpha = 1, linetype = 'solid', color = 'black', linewidth = .7) +
#   geom_line(alpha = .5, linewidth = .3) +
#   coord_cartesian(ylim = c(-1.5, .5)) +
#   scale_x_discrete(limits = xorder) +
#   facet_wrap(vars(type), ncol=1, scales = 'free_y',
#              labeller = type_labeller) +
#   labs(title = 'Median ATE for AY pairs',
#        y = 'Median ATE', x = 'Estimation Method') +
#   theme(legend.position = 'none',
#         axis.text.x = element_text(angle = 45, vjust = 1.2, hjust = 1.1, size = 10))
# ggsave(sprintf('%s/AY_medianATE_lineplot.pdf', plot_savepath), height = 8, width = 8)
# ggsave(sprintf('%s/AY_medianATE_lineplot.svg', plot_savepath), height = 8, width = 8)
# ggsave(sprintf('%s/AY_medianATE_lineplot.jpg', plot_savepath), height = 8, width = 8)

# ============================================================================
# df |> filter(type == 'positive' | type == 'negative') |> 
#   # group_by(type, AY_idx, EstMethod) |> 
#   mutate(AY_idx = as.factor(AY_idx)) |>
#   reframe(ATE_median = median(ATE),
#             ATE_upper  = sort(ATE, na.last = TRUE)[ceiling(n()*.5 + 1.96*sqrt(n()*.5*(1-.5)))],
#             ATE_lower  = sort(ATE, na.last = TRUE)[  floor(n()*.5 - 1.96*sqrt(n()*.5*(1-.5)))],
#           count = n(),
#           # .by = c(AY_idx, numNC, method_type, basis)
#           .by = c(AY_idx, EstMethod))
# 
# df |> filter(type == 'positive' | type == 'negative') |> 
#   group_by(type, AY_idx, EstMethod) |>
#   mutate(AY_idx = as.factor(AY_idx)) |>
#   summarize(count = n(), meanATE = mean(ATE),
#           .groups = 'drop')
# 
# df |> filter(AY_idx == 4 & EstMethod == 'OCB2SLS')

# ============================================================================
# plot_df = df |> filter(type == 'positive' | type == 'negative') |> 
#   # group_by(type, AY_idx, EstMethod) |> 
#   mutate(AY_idx = as.factor(AY_idx)) |>
#   reframe(ATE_median = median(ATE),
#           ATE_upper  = sort(ATE, na.last = TRUE)[ceiling(n()*.5 + 1.96*sqrt(n()*.5*(1-.5)))],
#           ATE_lower  = sort(ATE, na.last = TRUE)[  floor(n()*.5 - 1.96*sqrt(n()*.5*(1-.5)))],
#           count = n(),
#           # .by = c(AY_idx, numNC, method_type, basis)
#           .by = c(AY_idx, EstMethod, type))
# # median w CI 
# ggplot(df |> filter(type == 'positive' | type == 'negative') |> 
#          # group_by(type, AY_idx, EstMethod) |> 
#          mutate(AY_idx = as.factor(AY_idx)) |>
#          reframe(ATE_median = median(ATE),
#                  ATE_upper  = sort(ATE, na.last = TRUE)[ceiling(n()*.5 + 1.96*sqrt(n()*.5*(1-.5)))],
#                  ATE_lower  = sort(ATE, na.last = TRUE)[  floor(n()*.5 - 1.96*sqrt(n()*.5*(1-.5)))],
#                  count = n(),
#                  # .by = c(AY_idx, numNC, method_type, basis)
#                  .by = c(AY_idx, EstMethod, type)),
#        aes(x = EstMethod, 
#            # y = ATE, #|> abs() |> log(), 
#            color = type,
#            fill = type,
#            group = AY_idx)) +
#   geom_hline(aes(yintercept = 0), alpha = 1, linetype = 'solid', color = 'black') +
#   geom_ribbon(aes(ymin = ATE_lower, ymax = ATE_upper),
#               alpha = .1,
#               linewidth = 0, color = NA) +
#   geom_line(aes(y = ATE_median), alpha = .3, linewidth = .5) +
#   coord_cartesian(ylim = c(-1.0, .5)) +
#   scale_x_discrete(limits = xorder) +
#   facet_wrap(vars(type), ncol=1, scales = 'free_y', 
#              labeller = type_labeller) +
#   labs(title = 'Median ATE for AY pairs',
#        y = 'Median ATE', x = 'Estimation Method') +
#   theme(legend.position = 'none',
#         axis.text.x = element_text(angle = 45, vjust = 1.2, hjust = 1.1, size = 10))
# ggsave(sprintf('%s/AY_medianATECI_lineplot.pdf', plot_savepath), height = 8, width = 8)
# ggsave(sprintf('%s/AY_medianATECI_lineplot.svg', plot_savepath), height = 8, width = 8)
# ggsave(sprintf('%s/AY_medianATECI_lineplot.jpg', plot_savepath), height = 8, width = 8)



# # Remove Outlier: ATE Estimates (average over ZW selections): Lines 
# # median
# # outlier AY_idx is 40: df |> filter(EstMethod == 'lmYA' & abs(ATE) >= 2 & type == 'negative') |> distinct(AY_idx)
# 
# ggplot(df |> filter(AY_idx != 40) |> # remove AY_idx == 40 (outlier)
#          filter(EstMethod %in% c('lmYA', paste0('CB', 1:20))) |># stop showing CB early
#          filter(type == 'positive' | type == 'negative') |> 
#          group_by(type, AY_idx, EstMethod) |> 
#          mutate(AY_idx = as.factor(AY_idx)) |>
#          summarize(ATE_median = median(ATE),
#                    .groups = 'drop'),
#        aes(x = EstMethod,  y = ATE_median,
#            color = type,
#            fill = type,
#            group = AY_idx)) +
#   geom_hline(aes(yintercept = 0), alpha = 1, linetype = 'solid', color = 'black', linewidth = .7) +
#   geom_line(alpha = .5, linewidth = .3) +
#   # scale_y_continuous(breaks = seq(-20, 20, by = 1)) +
#   scale_color_manual(values = c('orangered2', 'steelblue4')) +
#   scale_fill_manual(values = c('orangered2', 'steelblue4')) +
#   facet_wrap(vars(type), ncol=1, scales = 'free_y',
#              labeller = type_labeller) +
#   labs(title = 'Median ATE for A-Y pairs',
#        y = 'Median ATE', x = 'Estimation Method') +
#   theme(legend.position = 'none',
#         axis.text.x = element_text(angle = 45, vjust = 1.2, hjust = 1.1, size = 10))
# ggsave(sprintf('%s/AY_medianATE_lineplot_removeoutlier.pdf', plot_savepath), height = 8, width = 9.5)
# ggsave(sprintf('%s/AY_medianATE_lineplot_removeoutlier.svg', plot_savepath), height = 8, width = 8)
# ggsave(sprintf('%s/AY_medianATE_lineplot_removeoutlier.jpg', plot_savepath), height = 8, width = 8)
# 
# 
# # median w CI 
# ggplot(df |> filter(AY_idx != 40) |> 
#          filter(EstMethod %in% c('lmYA', paste0('CB', 1:20))) |># stop showing CB early
#          filter(type == 'positive' || type == 'negative') |> 
#          group_by(type, AY_idx, EstMethod) |> 
#          mutate(AY_idx = as.factor(AY_idx)) |>
#          summarize(ATE_median = median(ATE),
#                    ATE_upper  = sort(ATE, na.last = TRUE)[ceiling(n()*.5 + 1.96*sqrt(n()*.5*(1-.5)))],
#                    ATE_lower  = sort(ATE, na.last = TRUE)[  floor(n()*.5 - 1.96*sqrt(n()*.5*(1-.5)))],
#                    .groups = 'drop'),
#        aes(x = EstMethod, 
#            # y = ATE, #|> abs() |> log(), 
#            color = type,
#            fill = type,
#            group = AY_idx)) +
#   geom_hline(aes(yintercept = 0), alpha = 1, linetype = 'solid', color = 'black') +
#   geom_ribbon(aes(ymin = ATE_lower, ymax = ATE_upper),
#               alpha = .1,
#               linewidth = 0, color = NA) +
#   geom_line(aes(y = ATE_median), alpha = .3, linewidth = .5) +
#   scale_y_continuous(breaks = seq(-20, 20, by = 1)) +
#   scale_color_manual(values = c('orangered2', 'steelblue4')) +
#   scale_fill_manual(values = c('orangered2', 'steelblue4')) +
#   facet_wrap(vars(type), ncol=1, scales = 'free_y',
#              labeller = type_labeller) +
#   labs(title = 'Median ATE for AY pairs',
#        y = 'Median ATE', x = 'Estimation Method') +
#   theme(legend.position = 'none',
#         axis.text.x = element_text(angle = 45, vjust = 1.2, hjust = 1.1, size = 10))
# ggsave(sprintf('%s/AY_medianATECI_lineplot_removeoutlier.pdf', plot_savepath), height = 8, width = 9.5)
# ggsave(sprintf('%s/AY_medianATECI_lineplot_removeoutlier.svg', plot_savepath), height = 8, width = 8)
# ggsave(sprintf('%s/AY_medianATECI_lineplot_removeoutlier.jpg', plot_savepath), height = 8, width = 8)





# =================== Make Plots: max NCs ATEs ======================================
print(sprintf("[%s]        - max ATEs", Sys.time()))
maxNC = max(df$numNC, na.rm = TRUE)
plot_savepath_subfolder = sprintf('%s/maxNC=%2d/', plot_savepath, maxNC)
dir.create(plot_savepath_subfolder, showWarnings = FALSE)


# # all ATE Estimates: Histogram by type
# ggplot(df |> filter(numNC == maxNC),
#        aes(x = ATE, y = after_stat(density))) +
#   geom_histogram(fill = 'orange', alpha = .6, binwidth = .1) +
#   geom_density(color = 'orange') +
#   geom_vline(aes(xintercept = 0), color = 'darkorange') +
#   scale_x_continuous(limits = c(-4, 4)) +
#   facet_wrap(vars(type), nrow = 3, labeller = type_labeller) +
#   # facet_grid(type ~ . , labeller = type_labeller) +
#   labs(title = 'Histogram of All ATE Estimates w/ Max NC',  
#        x = 'ATE Estimate') 
# ggsave(sprintf('%s/ATE_hist_maxNC.pdf', plot_savepath_subfolder), height = 5, width = 5)
# ggsave(sprintf('%s/ATE_hist_maxNC.svg', plot_savepath_subfolder), height = 5, width = 5)
# ggsave(sprintf('%s/ATE_hist_maxNC.jpg', plot_savepath_subfolder), height = 5, width = 5)


# all ATE Estimates: Boxplot by type
# ggplot(df[sample(1:nrow(df), size = 1000, ), ]) +
ggplot(df |> filter(numNC == maxNC | is.na(numNC))) +
  geom_hline(aes(yintercept = 0), alpha = 1, linetype = 'solid') +
  geom_boxplot(aes(x = EstMethod, y = ATE, 
                   # pattern_fill  = basis,
                   # alpha = basis, 
                   fill = method),
               # fill = 'orange', 
               # alpha = .6,
               linewidth = .1,
               outlier.shape = NA, # don't show outliers
               outlier.size = .1, outlier.alpha = .3) +
  coord_cartesian(ylim = ylim*2) +
  scale_x_discrete(limits = xorder, labels = EstMethod_map) +
  # scale_alpha_discrete(guide = 'none') +
  facet_wrap(vars(type), nrow = 3, labeller = type_labeller) +
  labs(title = 'Boxplots of All ATE Estimates (w/ max NC)',
       x = 'Estimation Method', y = 'ATE Estimate') +
  theme(axis.text.x = element_text(angle = -65, 
                                   vjust = 1, hjust = 0,
                                   size = 10),
        axis.ticks.x = element_blank(),
        legend.position = 'none')
ggsave(sprintf('%s/ATE_boxplot_maxNC.pdf', plot_savepath_subfolder), height = 5, width = 7, scale = 1)
ggsave(sprintf('%s/ATE_boxplot_maxNC.svg', plot_savepath_subfolder), height = 5, width = 7, scale = 1)
ggsave(sprintf('%s/ATE_boxplot_maxNC.jpg', plot_savepath_subfolder), height = 5, width = 7, scale = 1)


# # all ATE Estimates: Violinplot by type
# ggplot(df |> filter(numNC == maxNC | is.na(numNC))) +
#   geom_hline(aes(yintercept = 0), alpha = 1, linetype = 'solid') +
#   geom_violin(aes(x = EstMethod, y = ATE, 
#                    # pattern_fill  = basis,
#                    alpha = basis, 
#                    fill = method),
#                # fill = 'orange', 
#                # alpha = .6,
#                color = 'black', scale = 'count') +
#   coord_cartesian(ylim = ylim*2) +
#   scale_x_discrete(limits = xorder) +
#   scale_alpha_discrete(guide = 'none') +
#   facet_wrap(vars(type), nrow = 3) +
#   labs(title = 'Boxplots of All ATE Estimates (w/ max NC)',
#        x = 'Estimation Method', y = 'ATE Estimate') +
#   theme(axis.text.x = element_text(angle = 45, vjust = 1.2, hjust = 1.3, size = 10))

# library(scales) # for prettier axis labels/scales
# vertical histograms
ggplot(df |> 
         mutate(EstMethod = factor(EstMethod,  levels = names(EstMethod_map), labels = EstMethod_map)) |>
         mutate(type = factor(type, levels = names(type_map), labels = type_map)) |> 
         # filter(type == 'maybe') |>
         filter(numNC == maxNC | is.na(numNC)) ) +
  geom_hline(aes(yintercept = 0), alpha = 1, linetype = 'solid') + 
  geom_histogram(aes(y= ATE, x = after_stat(density), group = type, fill = method), 
                 orientation = 'y', position = 'identity', binwidth = .1,
                 alpha = .8) +
  scale_y_continuous(limits = c(-1.25, 1.25)) +
  # scale_x_continuous(breaks = breaks_pretty()) +
  scale_x_continuous(breaks = function(x) unique(floor(pretty(seq(min(x), (max(x) + 1) * 1.1))))) + # https://stackoverflow.com/questions/15622001/how-to-display-only-integer-values-on-an-axis-using-ggplot2
  facet_grid(# rows = vars(type), cols = vars(EstMethod), 
             type ~ EstMethod,
             scales = 'free' # , 
             # labeller = labeller(EstMethod = EstMethod_map, type   = type_map)
             ) +
  theme(axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 10), 
        # legend.text = element_text(size = 8),
        legend.position = 'none', 
        strip.text.x = element_text(size = 7),
        strip.text.y = element_text(size = 8))

ggsave(sprintf('%s/ATE_hist_maxNC_allEstMethod.pdf', plot_savepath_subfolder), height = 5, width = 8)
ggsave(sprintf('%s/ATE_hist_maxNC_allEstMethod.svg', plot_savepath_subfolder), height = 5, width = 8)
ggsave(sprintf('%s/ATE_hist_maxNC_allEstMethod.jpg', plot_savepath_subfolder), height = 5, width = 8)



# vertical histograms: subset of estimation methods
ggplot(df |> filter(EstMethod %in% c('lmYA', 'OCB2SLSpci2s', 'OCB2SLSReg', 'OCBMEstRWbasis1', 'OCBMEstRWRegbasis1')) |>
         mutate(EstMethod = factor(EstMethod,  levels = names(EstMethod_map), labels = EstMethod_map)) |>
         mutate(type = factor(type, levels = names(type_map), labels = type_map)) |> 
         # filter(type == 'maybe') |>
         filter(numNC == maxNC | is.na(numNC)) ) +
  geom_hline(aes(yintercept = 0), alpha = 1, linetype = 'solid') + 
  geom_histogram(aes(y= ATE, x = after_stat(density), group = type, fill = method), 
                 orientation = 'y', position = 'identity', binwidth = .125,
                 alpha = .8) +
  scale_y_continuous(limits = c(-1.25, 1.25)) +
  # scale_x_continuous(breaks = breaks_pretty()) +
  scale_x_continuous(breaks = function(x) unique(floor(pretty(seq(min(x), (max(x) + 1) * 1.1))))) + # https://stackoverflow.com/questions/15622001/how-to-display-only-integer-values-on-an-axis-using-ggplot2
  labs(y = 'ATE Estimate') +
  facet_grid(# rows = vars(type), cols = vars(EstMethod), 
    type ~ EstMethod,
    scales = 'free' # , 
    # labeller = labeller(EstMethod = EstMethod_map, type   = type_map)
  ) +
  theme(axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 10), 
        # legend.text = element_text(size = 8),
        legend.position = 'none', 
        strip.text.x = element_text(size = 10),
        strip.text.y = element_text(size = 10))

ggsave(sprintf('%s/ATE_hist_maxNC_someEstMethod.pdf', plot_savepath_subfolder), height = 5, width = 8)
ggsave(sprintf('%s/ATE_hist_maxNC_someEstMethod.svg', plot_savepath_subfolder), height = 5, width = 8)
ggsave(sprintf('%s/ATE_hist_maxNC_someEstMethod.jpg', plot_savepath_subfolder), height = 5, width = 8)

# # look at naive estimates of ATE 
# plot_df = df |> 
#   filter(type == 'positive' & method_type %in% c('naive')) |>
#   mutate(AY_name = paste0(A, ":", Y))
# 
# ggplot(plot_df, aes(x = ATE)) +
#   geom_histogram(fill = 'orange', alpha = .6, binwidth = .1)
# 
# =================== Make Plots: pvals vs Unif ======================================
print(sprintf("[%s]        - pvals vs unif", Sys.time()))
plot_savepath_subfolder = sprintf('%s/pvals/', plot_savepath)
dir.create(plot_savepath_subfolder, showWarnings = FALSE)

# =================== Make Plots: pvals vs Unif: maxNC ======================================
plot_savepath_subsubfolder = sprintf('%smaxNC', plot_savepath_subfolder)
dir.create(plot_savepath_subsubfolder, showWarnings = FALSE)

lmYApvals  = df |> filter(type == 'negative' & method == 'lmYA')                          |> pull(pval) |> sort()
pci2spvals = df |> filter(type == 'negative' & method == 'OCB2SLSpci2s' & numNC == maxNC) |> pull(pval) |> sort()
mestpvals  = df |> filter(type == 'negative' & method == 'OCBMEst'      & numNC == maxNC) |> pull(pval) |> sort()

# lmYAxunif  = runif(length( lmYApvals), min = 0, max = 1) |> sort()
# pci2sxunif = runif(length(pci2spvals), min = 0, max = 1) |> sort()
# mestxunif  = runif(length( mestpvals), min = 0, max = 1) |> sort()
# lmYAppoints  = ppoints(length( lmYApvals))
# pci2sppoints = ppoints(length(pci2spvals))
# mestppoints  = ppoints(length( mestpvals))


plot_df = rbind(data.frame(EstMethod = 'lmYA',
                           pvals = lmYApvals,
                           theoretical = ppoints(length(lmYApvals))),
                data.frame(EstMethod = 'pci2s',
                           pvals = pci2spvals,
                           theoretical = ppoints(length(pci2spvals))),
                data.frame(EstMethod = 'MEst',
                           pvals = mestpvals,
                           theoretical = ppoints(length(mestpvals)))) |> 
  data.frame() |>
  mutate(EstMethod = factor(EstMethod, levels = c('lmYA', 'pci2s', 'MEst'))) |>
  mutate(pvals_trans = -log(pvals, base = 10),
         theoretical_trans = -log(theoretical, base = 10))


method_colors = c("#66C2A5", "#e93ba6", "#60ad25")

# qqplot of pvals vs unif
ggplot(plot_df, aes(sample=pvals, group = EstMethod, color = EstMethod)) +
  stat_qq_line(distribution = stats::qunif, 
               # color = 'purple', 
               color = 'black',
               line.p = c(0, 1)) +
  stat_qq(distribution = stats::qunif, size = 1, alpha = .8) +
  # scale_color_brewer(palette = "Set2") +
  scale_color_discrete(type = method_colors) +
  labs(title = 'QQ Plot of p-values against Unif(0,1)',
       x = 'theoretical Unif(0,1)', y = 'sample p-values') +
  theme(legend.position = 'inside', 
        legend.position.inside = c(.1, .75))

ggsave(sprintf('%s/pval_qqunif_maxNC.pdf', plot_savepath_subsubfolder), height = 5, width = 5)
ggsave(sprintf('%s/pval_qqunif_maxNC.svg', plot_savepath_subsubfolder), height = 5, width = 5)
ggsave(sprintf('%s/pval_qqunif_maxNC.jpg', plot_savepath_subsubfolder), height = 5, width = 5)


# qqplot of transformed pvals vs unif
p1 = ggplot(plot_df, 
       aes(# x = theoretical, y = pvals,
           x = theoretical_trans, y = pvals_trans,
           group = EstMethod, color = EstMethod)) +
  geom_point() +
  # stat_qq_line(distribution = stats::qunif, color = 'purple', line.p = c(0, 1)) +
  # stat_qq(distribution = stats::qunif, size = 1, alpha = .8) +
  geom_abline(slope = 1, intercept = 0, 
              # color = 'purple',
              color = 'black') +
  # scale_color_brewer(palette = "Set2") +
  scale_color_discrete(type = method_colors) +
  labs(title = '-log(p-values) vs -log(Unif(0,1))',
       x = 'theoretical -log(Unif(0,1))', y = 'sample -log(pval)') +
  theme(legend.position = 'inside', 
        legend.position.inside = c(.1, .75))

p1 
ggsave(sprintf('%s/pval_qqunif_maxNC_neglog.pdf', plot_savepath_subsubfolder), height = 5, width = 5)
ggsave(sprintf('%s/pval_qqunif_maxNC_neglog.svg', plot_savepath_subsubfolder), height = 5, width = 5)
ggsave(sprintf('%s/pval_qqunif_maxNC_neglog.jpg', plot_savepath_subsubfolder), height = 5, width = 5)


p1 + scale_y_continuous(limits = c(0,50))
ggsave(sprintf('%s/pval_qqunif_maxNC_neglog_cutoff50.pdf', plot_savepath_subsubfolder), height = 5, width = 5)
ggsave(sprintf('%s/pval_qqunif_maxNC_neglog_cutoff50.svg', plot_savepath_subsubfolder), height = 5, width = 5)
ggsave(sprintf('%s/pval_qqunif_maxNC_neglog_cutoff50.jpg', plot_savepath_subsubfolder), height = 5, width = 5)

p1 + scale_y_continuous(limits = c(0,25))
ggsave(sprintf('%s/pval_qqunif_maxNC_neglog_cutoff25.pdf', plot_savepath_subsubfolder), height = 5, width = 5)
ggsave(sprintf('%s/pval_qqunif_maxNC_neglog_cutoff25.svg', plot_savepath_subsubfolder), height = 5, width = 5)
ggsave(sprintf('%s/pval_qqunif_maxNC_neglog_cutoff25.jpg', plot_savepath_subsubfolder), height = 5, width = 5)



# Histogram of p-values
ggplot(df |> filter(numNC == maxNC | is.na(numNC)) |>
             filter(method %in% c('lmYA', 'OCB2SLSpci2s', 'OCBMEst')), 
       aes(x = pval, fill = EstMethod)) +
  geom_histogram(aes(), breaks = seq(from = 0, to = 1, by = .025), 
                 # fill = "#8DA0CB",
                 position = 'identity') +
  scale_fill_discrete(type = method_colors) +
  # facet_grid(rows = vars(type), cols = vars(method), scales = 'free_y', axes = 'margins') +
  ggh4x::facet_grid2(rows = vars(type), cols = vars(method), scales = 'free_y', independent = 'y') +
  labs(x = 'p-value')
ggsave(sprintf('%s/pval_hist_maxNC.pdf', plot_savepath_subsubfolder), height = 5, width = 8)
ggsave(sprintf('%s/pval_hist_maxNC.svg', plot_savepath_subsubfolder), height = 5, width = 8)
ggsave(sprintf('%s/pval_hist_maxNC.jpg', plot_savepath_subsubfolder), height = 5, width = 8)

# Histogram of p-values (only noncausal)
ggplot(df |> filter(numNC == maxNC | is.na(numNC)) |>
         filter(method %in% c('lmYA', 'OCB2SLSpci2s', 'OCBMEst')) |>
         filter(type == 'negative'), 
       aes(x = pval, fill = EstMethod)) +
  geom_histogram(aes(y = after_stat(density)), breaks = seq(from = 0, to = 1, by = .025), 
                 # fill = "#8DA0CB",
                 position = 'identity') +
  scale_fill_discrete(type = c("#66C2A5", "#d354a3", "#72ac45")) +
  # facet_grid(rows = vars(type), cols = vars(method), scales = 'free_y', axes = 'margins') +
  ggh4x::facet_grid2(cols = vars(method), scales = 'free_y', independent = 'y') +
  labs(x = 'p-value') +
  theme(strip.background = element_blank())
ggsave(sprintf('%s/pval_hist_maxNC_negative.pdf', plot_savepath_subsubfolder), height = 3, width = 10)
ggsave(sprintf('%s/pval_hist_maxNC_negative.svg', plot_savepath_subsubfolder), height = 3, width = 10)
ggsave(sprintf('%s/pval_hist_maxNC_negative.jpg', plot_savepath_subsubfolder), height = 3, width = 10)


# =================== Make Plots: pvals vs Unif: individual type and NC ======================================
plot_savepath_subsubfolder = sprintf('%stype_numNC', plot_savepath_subfolder)
dir.create(plot_savepath_subsubfolder, showWarnings = FALSE)

#' @param myType (character) 'positive' 'negative' 'maybe'
#' @param myNumNC (integer) number of NC 
#' @example make_pval_plot('negative', 3)
make_pval_plot <- function(myType, myNumNC) {
  
  lmYApvals  = df |> filter(type == myType & method == 'lmYA')                          |> pull(pval) |> sort()
  pci2spvals = df |> filter(type == myType & method == 'OCB2SLSpci2s' & numNC == myNumNC) |> pull(pval) |> sort()
  mestpvals  = df |> filter(type == myType & method == 'OCBMEst'      & numNC == myNumNC) |> pull(pval) |> sort()

  
  
  plot_df = rbind(data.frame(EstMethod = 'lmYA',
                             pvals = lmYApvals,
                             theoretical = ppoints(length(lmYApvals))),
                  data.frame(EstMethod = 'pci2s',
                             pvals = pci2spvals,
                             theoretical = ppoints(length(pci2spvals))),
                  data.frame(EstMethod = 'MEst',
                             pvals = mestpvals,
                             theoretical = ppoints(length(mestpvals)))) |> 
    data.frame() |>
    mutate(EstMethod = factor(EstMethod, levels = c('lmYA', 'pci2s', 'MEst'))) |>
    mutate(pvals_trans = -log(pvals, base = 10),
           theoretical_trans = -log(theoretical, base = 10))
  
  # qqplot of pvals vs unif
  p1 = ggplot(plot_df, aes(sample=pvals, group = EstMethod, color = EstMethod)) +
    geom_abline(slope = 1, intercept = 0, color = 'purple3', linewidth = 1, alpha = .8) +
    # stat_qq_line(distribution = stats::qunif, color = 'purple', line.p = c(0, 1)) +
    stat_qq(distribution = stats::qunif, size = 1, alpha = .8) +
    # scale_color_brewer(palette = "Set2") +
    scale_color_discrete(type = method_colors) +
    labs(title = 'QQ Plot of p-values against Unif(0,1)',
         x = 'theoretical Unif(0,1)', y = 'sample p-values') +
    theme(legend.position = 'inside', 
          legend.position.inside = c(.1, .75))
  
  ggsave(sprintf('%s/pval_qqunif_type=%s_numNC=%02d.pdf', plot_savepath_subsubfolder, myType, myNumNC), 
         plot = p1, height = 5, width = 5)
  ggsave(sprintf('%s/pval_qqunif_type=%s_numNC=%02d.svg', plot_savepath_subsubfolder, myType, myNumNC), 
         plot = p1, height = 5, width = 5)
  ggsave(sprintf('%s/pval_qqunif_type=%s_numNC=%02d.jpg', plot_savepath_subsubfolder, myType, myNumNC), 
         plot = p1, height = 5, width = 5)
  
  
  # qqplot of transformed pvals vs unif
  p2 = ggplot(plot_df, 
         aes(# x = theoretical, y = pvals,
           x = theoretical_trans, y = pvals_trans,
           group = EstMethod, color = EstMethod)) +
    geom_point() +
    # stat_qq_line(distribution = stats::qunif, color = 'purple', line.p = c(0, 1)) +
    # stat_qq(distribution = stats::qunif, size = 1, alpha = .8) +
    geom_abline(slope = 1, intercept = 0, color = 'purple3', linewidth = 1, alpha = .8) +
    # scale_color_brewer(palette = "Set2") +
    scale_color_discrete(type = method_colors) +
    labs(title = '-log(p-values) vs -log(Unif(0,1))',
         x = 'theoretical -log(Unif(0,1))', y = 'sample -log(pval)') +
    theme(legend.position = 'inside', 
          legend.position.inside = c(.1, .75))
  
  ggsave(sprintf('%s/pval_qqunif_type=%s_numNC=%02d_neglog.pdf', plot_savepath_subsubfolder, myType, myNumNC),
         plot = p2, height = 5, width = 5)
  ggsave(sprintf('%s/pval_qqunif_type=%s_numNC=%02d_neglog.svg', plot_savepath_subsubfolder, myType, myNumNC),
         plot = p2, height = 5, width = 5)
  ggsave(sprintf('%s/pval_qqunif_type=%s_numNC=%02d_neglog.jpg', plot_savepath_subsubfolder, myType, myNumNC),
         plot = p2, height = 5, width = 5)
  
  # cutoff -log(p) axis [0, c]
  ggsave(sprintf('%s/pval_qqunif_type=%s_numNC=%02d_neglog_cutoff50.pdf', plot_savepath_subsubfolder, myType, myNumNC),
         plot = p2 + scale_y_continuous(limits = c(0,50)), height = 5, width = 5)
  ggsave(sprintf('%s/pval_qqunif_type=%s_numNC=%02d_neglog_cutoff50.svg', plot_savepath_subsubfolder, myType, myNumNC),
         plot = p2 + scale_y_continuous(limits = c(0,50)), height = 5, width = 5)
  ggsave(sprintf('%s/pval_qqunif_type=%s_numNC=%02d_neglog_cutoff50.jpg', plot_savepath_subsubfolder, myType, myNumNC),
         plot = p2 + scale_y_continuous(limits = c(0,50)), height = 5, width = 5)
  
  ggsave(sprintf('%s/pval_qqunif_type=%s_numNC=%02d_neglog_cutoff25.pdf', plot_savepath_subsubfolder, myType, myNumNC),
         plot = p2 + scale_y_continuous(limits = c(0,25)), height = 5, width = 5)
  ggsave(sprintf('%s/pval_qqunif_type=%s_numNC=%02d_neglog_cutoff25.svg', plot_savepath_subsubfolder, myType, myNumNC),
         plot = p2 + scale_y_continuous(limits = c(0,25)), height = 5, width = 5)
  ggsave(sprintf('%s/pval_qqunif_type=%s_numNC=%02d_neglog_cutoff25.jpg', plot_savepath_subsubfolder, myType, myNumNC),
         plot = p2 + scale_y_continuous(limits = c(0,25)), height = 5, width = 5)
  
  
}


for(myType in c('negative', 'positive')) {
  for(myNumNC in sort(unique(df$numNC))) {
    # print(sprintf('%s : %s', myType, myNumNC))
    make_pval_plot(myType, myNumNC)
  }
}


# =================== Make Plots: pvals vs Unif: individual type and est method over numNC ======================================
# Compare qqplots by NumNC (to see the improvement as numNC increases)

plot_savepath_subsubfolder = sprintf('%sbyNumNC', plot_savepath_subfolder)
dir.create(plot_savepath_subsubfolder, showWarnings = FALSE)



#' @param myMethod (character) Estimation method 
#' @param myType (character) 'positive' 'negative' 'maybe'
#' @example make_pvalNumNC_plot(myMethod = 'OCB2SLSpci2s', myType = 'positive')
make_pvalNumNC_plot <- function(myMethod, myType) {
  # myType = 'negative'
  # myMethod = 'OCB2SLSpci2s'
  df_subset = df |> filter(type == myType & method == myMethod)
  
  plot_df = NULL
  for(curNumNC in sort(unique(df_subset$numNC))) {
    pvals = df_subset |> filter(numNC == curNumNC) |> pull(pval) |> sort()
    plot_df = rbind(plot_df, 
                    data.frame(numNC = curNumNC,
                               pvals = pvals,
                               theoretical = ppoints(length(pvals))))
  }
  plot_df = plot_df  |>
            mutate(pvals_trans = -log(pvals, base = 10),
                   theoretical_trans = -log(theoretical, base = 10)) 
  
  # qqplot of pvals vs unif
  p1 = ggplot(plot_df, aes(sample=pvals, group = numNC, color = numNC)) +
    geom_abline(slope = 1, intercept = 0, 
                # color = 'purple3', 
                color = 'black',
                linewidth = 1, alpha = .8) +
    # stat_qq_line(distribution = stats::qunif, color = 'purple', line.p = c(0, 1)) +
    stat_qq(distribution = stats::qunif, 
            # size = .8, 
            alpha = .8, 
            # geom = 'point',
            geom = 'line', linewidth = 1.3) +
    # scale_color_brewer(palette = "Greens") +
    scale_color_continuous(type = 'viridis') +
    labs(title = 'QQ Plot of p-values against Unif(0,1)',
         x = 'theoretical Unif(0,1)', y = 'sample p-values', color = '# NC') +
    theme(legend.position = 'inside', 
          legend.position.inside = c(.1, .75), 
          legend.title = element_text(hjust = 0))
  
  ggsave(sprintf('%s/pval_qqunif_method=%s_type=%s.pdf', plot_savepath_subsubfolder, myMethod, myType), 
         plot = p1, height = 5, width = 5)
  ggsave(sprintf('%s/pval_qqunif_method=%s_type=%s.svg', plot_savepath_subsubfolder, myMethod, myType), 
         plot = p1, height = 5, width = 5)
  ggsave(sprintf('%s/pval_qqunif_method=%s_type=%s.jpg', plot_savepath_subsubfolder, myMethod, myType), 
         plot = p1, height = 5, width = 5)
  
  
  # qqplot of transformed pvals vs unif
  p2 = ggplot(plot_df, 
              aes(# x = theoretical, y = pvals,
                x = theoretical_trans, y = pvals_trans,
                group = numNC, color = numNC)) +
    geom_abline(slope = 1, intercept = 0, 
                # color = 'purple3', 
                color = 'black',
                linewidth = 1, alpha = .8) +
    # stat_qq_line(distribution = stats::qunif, color = 'purple', line.p = c(0, 1)) +
    # stat_qq(distribution = stats::qunif, size = 1, alpha = .8) +
    geom_point(alpha = .8) +
    # geom_line(linewidth = 1.3) +
    # scale_color_brewer(palette = "Set2") +
    scale_color_continuous(type = 'viridis') +
    labs(title = '-log(p-values) vs -log(Unif(0,1))',
         x = 'theoretical -log(Unif(0,1))', y = 'sample -log(pval)', color = '# NC') +
    theme(legend.position = 'inside', 
          legend.position.inside = c(.1, .75),
          legend.title = element_text(hjust = 0))
  
  ggsave(sprintf('%s/pval_qqunif_method=%s_type=%s_neglog.pdf', plot_savepath_subsubfolder, myMethod, myType), 
         plot = p2, height = 5, width = 5)
  ggsave(sprintf('%s/pval_qqunif_method=%s_type=%s_neglog.svg', plot_savepath_subsubfolder, myMethod, myType), 
         plot = p2, height = 5, width = 5)
  ggsave(sprintf('%s/pval_qqunif_method=%s_type=%s_neglog.jpg', plot_savepath_subsubfolder, myMethod, myType), 
         plot = p2, height = 5, width = 5)

}



for(myMethod in c('OCB2SLSpci2s', 'OCBMEst')) {
  for(myType in c('negative', 'positive')) {
    make_pvalNumNC_plot(myMethod = myMethod, myType = myType)
  }
}


# qqplot(-log10(ppoints(length(lmYApvals))), 
#                       -log10(lmYApvals), 
#        xlab = "theoretical", ylab = "obs'd", main = "Q-Q Plot for -log10 Pval", 
#        cex = 0.5, pch = 3) 
# abline(0, 1, col = "red")
# 
# plot(-log(sort(runif(length(pci2spvals), min = 0, max = 1))), 
#      -log(sort(lmYApvals)))
# 
# 
# # qqplot of pvals vs unif
# ggplot(plot_df, aes(sample=pvals, group = EstMethod, color = EstMethod)) +
#   # stat_qq_line(distribution = stats::qunif, color = 'purple', line.p = c(0, 1)) +
#   stat_qq(distribution = stats::qunif, size = 1, alpha = .8) +
#   scale_color_brewer(palette = "Set2") +
#   labs(title = 'QQ Plot of ')
# 
# hist(df |> filter(type == 'negative' & method == 'lmYA') |> pull(pval), breaks = seq(0, 1, by=.005))
# hist(df |> filter(type == 'negative' & method == 'OCB2SLSpci2s' & numNC == maxNC) |> pull(pval), breaks = seq(0, 1, by=.005))
# hist(df |> filter(type == 'negative' & method == 'OCBMEst' & numNC == maxNC) |> pull(pval), breaks = seq(0, 1, by=.005))
# 
# 
# 
# 
# qqplot( lmYApvals, xunif); abline(a = 0, b = 1)
# qqplot( pci2spvals, xunif); abline(a = 0, b = 1)
# qqplot( mestpvals, xunif); abline(a = 0, b = 1)
# qqplot( -log( lmYApvals), -log(xunif)); abline(a = 0, b = 1)
# qqplot(qunif(lmYApvals, min = 0, max = 1), lmYApvals)
# 
# make_pval_plot <- function(pvals, title = 'pvals vs Unif(0,1)') {
#   ggplot(NULL, aes(sample=pvals)) +
#     stat_qq_line(distribution = stats::qunif, color = 'orange', line.p = c(0, 1)) +
#     stat_qq(distribution = stats::qunif, color = 'purple') + 
#     labs(title = title)
# }
# 
# 
# make_pval_plot(lmYApvals)
# 
# 
# 
# 

# =================== Make Plots: by numNC ======================================
print(sprintf("[%s]        - by numNC", Sys.time()))
# Save in separate folder
plot_savepath_subfolder = sprintf('%s/numNC/', plot_savepath)
dir.create(plot_savepath_subfolder, showWarnings = FALSE)

plot_df = df |> 
  # filter(EstMethod %in% c('lmYA', 'OCB2SLSpci2s', 'OCB2SLSReg', 'OCBMEstRWbasis1', 'OCBMEstRWRegbasis1'))
  filter(EstMethod %in% c('lmYA', 'OCB2SLSpci2s', 'OCBMEstRWbasis1'))
xorder_l_cur = xorder_l[xorder_l %in% unique(plot_df$EstMethod_l)]
ggplot(plot_df       ) +
  geom_hline(aes(yintercept = 0), alpha = 1, linetype = 'solid', color = 'gray') +
  geom_boxplot(aes(x = EstMethod_l, y = (ATE),
                   alpha = basis, fill = method),
               # fill = 'orange', alpha = .6,
               outlier.size = .1, outlier.alpha = .3) +
  coord_cartesian(ylim = c(-1.5, .5)) +
  # coord_cartesian(ylim = c(-10, 10)) +
  scale_x_discrete(limits = xorder_l_cur) +
  scale_fill_brewer(palette = 'Set2') +
  scale_alpha_discrete(guide = 'none') +
  facet_wrap(vars(type), nrow = 3,
             labeller = type_labeller) +
  # labs(title = 'Median ATE for AY pairs (w/ max NC)',
  #      y = 'Median ATE', x = 'Estimation Method') +
  labs(title = 'ATE by Estimation Method',
       y = 'ATE', x = 'Estimation Method') +
  theme(axis.text.x = element_text(angle = 30, vjust = 1.1, hjust = 1.2, size = 14),
        legend.position = 'none')
ggsave(sprintf('%s/ATE_boxplot_numNC.pdf', plot_savepath_subfolder), height = 6, width = 8)
ggsave(sprintf('%s/ATE_boxplot_numNC.svg', plot_savepath_subfolder), height = 6, width = 8)
ggsave(sprintf('%s/ATE_boxplot_numNC.jpg', plot_savepath_subfolder), height = 6, width = 8)


# vertical histograms: subset of estimation methods
ggplot(plot_df) +
  geom_hline(aes(yintercept = 0), alpha = 1, linetype = 'solid') + 
  geom_histogram(aes(y= ATE, x = after_stat(density), group = type, fill = method), 
                 orientation = 'y', position = 'identity', binwidth = .125,
                 alpha = .8) +
  scale_y_continuous(limits = c(-1.25, 1.25)) +
  # scale_x_continuous(breaks = breaks_pretty()) +
  scale_x_continuous(breaks = function(x) unique(floor(pretty(seq(min(x), (max(x) + 1) * 1.1))))) + # https://stackoverflow.com/questions/15622001/how-to-display-only-integer-values-on-an-axis-using-ggplot2
  facet_grid(# rows = vars(type), cols = vars(EstMethod), 
    type ~ EstMethod_l,
    scales = 'free' # , 
    # labeller = labeller(EstMethod = EstMethod_map, type   = type_map)
  ) +
  theme(axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 10), 
        # legend.text = element_text(size = 8),
        legend.position = 'none', 
        strip.text.x = element_text(size = 10),
        strip.text.y = element_text(size = 10))
ggsave(sprintf('%s/ATE_hist_numNC.pdf', plot_savepath_subfolder), height = 5, width = 8)
ggsave(sprintf('%s/ATE_hist_numNC.svg', plot_savepath_subfolder), height = 5, width = 8)
ggsave(sprintf('%s/ATE_hist_numNC.jpg', plot_savepath_subfolder), height = 5, width = 8)

# =================== Make Plots: individual AY pairs (Negative then Positive) ======================================
print(sprintf("[%s]        - individual AY pairs", Sys.time()))
# Save individual AY pairs in separate folder
dir.create(sprintf('%s/individualAY/', plot_savepath), showWarnings = FALSE)

# Choosing some NEGATIVE AY pairs to look at
some_pair_id = df |> 
  filter(type == 'negative' & EstMethod == 'lmYA' & abs(ATE) > .03) |> 
  distinct(AY_idx, A, Y, ATE) |>
  arrange(desc(abs(ATE)))
ggplot(df |> 
         filter(type == 'negative' & 
                  AY_idx %in% some_pair_id[1:9, 'AY_idx']) |>
         # filter(is.na(numNC) | numNC == 20) |>
         mutate(AY_name = paste0(A, ":", Y)),
       aes(x = EstMethod, y = ATE)) +
  geom_hline(aes(yintercept = 0), alpha = 1, linetype = 'solid', color = 'gray26') +
  geom_boxplot(alpha = .3, linewidth = .3,
               outlier.shape = '.', outlier.size = .3,
               outlier.stroke = 0, outlier.alpha = .5,
               fill = 'darkorchid1', color = 'darkorchid3') +
  coord_cartesian(ylim = c(-10, 10)) +
  scale_x_discrete(limits = xorder) +
  # scale_y_continuous(limits = quantile(ATE, c(0.1, 0.9), na.rm=TRUE)) +
  facet_wrap(vars(AY_name), nrow = 3, scales = 'free_y') +
  labs(x = 'Estimation Method', y = 'ATE Estimate') +
  theme(axis.text.x = element_text(angle = 45, vjust = 1.2, hjust = 1.1, size = 6))
ggsave(sprintf('%s/individualAY/AY_allATE_boxplot_negative.pdf', plot_savepath), height = 8, width = 12)
ggsave(sprintf('%s/individualAY/AY_allATE_boxplot_negative.svg', plot_savepath), height = 8, width = 12)
ggsave(sprintf('%s/individualAY/AY_allATE_boxplot_negative.jpg', plot_savepath), height = 8, width = 12)


# Choosing some POSITIVE AY pairs to look at
some_pair_id = df |> 
  filter(type == 'positive' & EstMethod == 'lmYA' & abs(ATE) > .03) |> 
  distinct(AY_idx, A, Y, ATE) |> 
  arrange(desc(abs(ATE)))
plot_df = df |> 
  filter(AY_idx %in% some_pair_id[1:9, 'AY_idx']) |>
  mutate(AY_name = paste0(A, ":", Y))
ggplot(plot_df,
       aes(x = EstMethod, y = ATE)) +
  geom_hline(aes(yintercept = 0), alpha = 1, linetype = 'solid', color = 'gray26') +
  geom_boxplot(alpha = .3, linewidth = .3,
               outlier.shape = '.', outlier.size = .3,
               outlier.stroke = 0, outlier.alpha = .5,
               fill = 'seagreen3', color = 'seagreen') +
  # coord_cartesian(ylim = c(-1, 2)) +
  scale_y_continuous(limits = quantile(plot_df$ATE, c(0.1, 0.9), na.rm=TRUE)) +
  facet_wrap(vars(AY_name), nrow = 3, scales = 'free_y') +
  labs(x = 'Estimation Method', y = 'ATE Estimate') +
  theme(axis.text.x = element_text(angle = 45, vjust = 1.2, hjust = 1.1, size = 6))
ggsave(sprintf('%s/individualAY/AY_allATE_boxplot_positive.pdf', plot_savepath), height = 8, width = 12)
ggsave(sprintf('%s/individualAY/AY_allATE_boxplot_positive.svg', plot_savepath), height = 8, width = 12)
ggsave(sprintf('%s/individualAY/AY_allATE_boxplot_positive.jpg', plot_savepath), height = 8, width = 12)





# Looking at estimates over different number of NCS
# negative pairs
some_pair_id = df |> 
  filter(type == 'negative' & EstMethod == 'lmYA' & abs(ATE) > .03) |> 
  distinct(AY_idx, A, Y, ATE) |> 
  arrange(desc(abs(ATE)))
plot_df = df |> 
  filter(type == 'negative' & method_type %in% c('naive', '2SLS', 'GMM', 'MEst', 'OCBLin') &
           AY_idx %in% some_pair_id[1:3, 'AY_idx']) |>
  mutate(AY_name = paste0(A, ":", Y))
ggplot(plot_df,
       aes(x = EstMethod_l, y = ATE, fill = method, color = method)) +
  geom_hline(aes(yintercept = 0), alpha = .5, linetype = 'solid', color = 'gray26') +
  geom_boxplot(alpha = .3, linewidth = .3,
               outlier.shape = '.', outlier.size = .3,
               outlier.stroke = 0, outlier.alpha = .5 #,
               # fill = 'darkorchid1', 
               # color = 'darkorchid3'
  ) +
  scale_x_discrete(limits = xorder_l) +
  # coord_cartesian(ylim = c(-1, 2)) +
  # scale_y_continuous(limits = quantile(plot_df$ATE, c(0.1, 0.9), na.rm=TRUE)) +
  coord_cartesian(ylim = c(-8, 2)) +
  facet_wrap(vars(AY_name), nrow = 3, scales = 'free_y') +
  labs(x = 'Estimation Method', y = 'ATE Estimate') +
  theme(axis.text.x = element_text(angle = 45, vjust = 1.2, hjust = 1.1, size = 6))

ggsave(sprintf('%s/individualAY/AY_allATE_boxplot_negative_numNC.pdf', plot_savepath), height = 8, width = 12)
ggsave(sprintf('%s/individualAY/AY_allATE_boxplot_negative_numNC.svg', plot_savepath), height = 8, width = 12)
ggsave(sprintf('%s/individualAY/AY_allATE_boxplot_negative_numNC.jpg', plot_savepath), height = 8, width = 12)

# positive pairs
some_pair_id = df |> 
  filter(type == 'positive' & EstMethod == 'lmYA' & abs(ATE) > .03) |> 
  distinct(AY_idx, A, Y, ATE) |> 
  arrange(desc(abs(ATE)))
plot_df = df |> 
  filter(type == 'positive' & method_type %in% c('naive', '2SLS', 'GMM', 'OCBLin') &
           AY_idx %in% some_pair_id[1:3, 'AY_idx']) |>
  mutate(AY_name = paste0(A, ":", Y))
ggplot(plot_df,
       aes(x = EstMethod_l, y = ATE, fill = method)) +
  geom_hline(aes(yintercept = 0), alpha = 1, linetype = 'solid', color = 'gray26') +
  geom_boxplot(alpha = .3, linewidth = .3,
               outlier.shape = '.', outlier.size = .3,
               outlier.stroke = 0, outlier.alpha = .5 #,
               # fill = 'darkorchid1', 
               # color = 'darkorchid3'
  ) +
  scale_x_discrete(limits = xorder_l) +
  coord_cartesian(ylim = c(-5, 2)) +
  # scale_y_continuous(limits = quantile(plot_df$ATE, c(0.1, 0.9), na.rm=TRUE)) +
  # coord_cartesian(ylim = c(-8.329534,  4.251341 )) +
  facet_wrap(vars(AY_name), nrow = 3, scales = 'free_y') +
  labs(x = 'Estimation Method', y = 'ATE Estimate') +
  theme(axis.text.x = element_text(angle = 45, vjust = 1.2, hjust = 1.1, size = 6))
ggsave(sprintf('%s/individualAY/AY_allATE_boxplot_positive_numNC.pdf', plot_savepath), height = 8, width = 12)
ggsave(sprintf('%s/individualAY/AY_allATE_boxplot_positive_numNC.svg', plot_savepath), height = 8, width = 12)
ggsave(sprintf('%s/individualAY/AY_allATE_boxplot_positive_numNC.jpg', plot_savepath), height = 8, width = 12)


# negative pairs
some_pair_id = df |> 
  filter(type == 'negative' & EstMethod == 'lmYA' & abs(ATE) > .03) |> 
  distinct(AY_idx, A, Y, ATE) |> 
  arrange(desc(abs(ATE)))
plot_df = df |> 
  filter(type == 'negative' & method_type %in% c('naive', '2SLS', 'GMM', 'OCBLin') &
           AY_idx %in% some_pair_id[sample(1:nrow(some_pair_id))[1:9], 'AY_idx']) |> #randomly tak e9
  mutate(AY_name = paste0(A, ":", Y))
ggplot(plot_df,
       aes(x = EstMethod_l, y = ATE, fill = method)) +
  geom_hline(aes(yintercept = 0), alpha = 1, linetype = 'solid', color = 'gray26') +
  geom_boxplot(alpha = .3, linewidth = .3,
               outlier.shape = '.', outlier.size = .3,
               outlier.stroke = 0, outlier.alpha = .5 #,
               # fill = 'darkorchid1', 
               # color = 'darkorchid3'
  ) +
  scale_x_discrete(limits = xorder_l) +
  # coord_cartesian(ylim = c(-1, 2)) +
  # scale_y_continuous(limits = quantile(plot_df$ATE, c(0.1, 0.9), na.rm=TRUE)) +
  # coord_cartesian(ylim = c(-5, 2)) +
  facet_wrap(vars(AY_name), nrow = 3, scales = 'free_y') +
  labs(x = 'Estimation Method', y = 'ATE Estimate') +
  theme(axis.text.x = element_text(angle = 45, vjust = 1.2, hjust = 1.1, size = 6))

ggsave(sprintf('%s/individualAY/AY_allATE_boxplot_negative_numNC_more.pdf', plot_savepath), height = 8, width = 12)
ggsave(sprintf('%s/individualAY/AY_allATE_boxplot_negative_numNC_more.svg', plot_savepath), height = 8, width = 12)
ggsave(sprintf('%s/individualAY/AY_allATE_boxplot_negative_numNC_more.jpg', plot_savepath), height = 8, width = 12)

# positive pairs
some_pair_id = df |> 
  filter(type == 'positive' & EstMethod == 'lmYA' & abs(ATE) > .03) |> 
  distinct(AY_idx, A, Y, ATE) |> 
  arrange(desc(abs(ATE)))
plot_df = df |> 
  filter(type == 'positive' & method_type %in% c('naive', '2SLS', 'GMM', 'OCBLin') &
           AY_idx %in% some_pair_id[sample(1:nrow(some_pair_id))[1:9], 'AY_idx']) |> #randomly tak e9
  mutate(AY_name = paste0(A, ":", Y))
ggplot(plot_df,
       aes(x = EstMethod_l, y = ATE, fill = method)) +
  geom_hline(aes(yintercept = 0), alpha = 1, linetype = 'solid', color = 'gray26') +
  geom_boxplot(alpha = .3, linewidth = .3,
               outlier.shape = '.', outlier.size = .3,
               outlier.stroke = 0, outlier.alpha = .5 #,
               # fill = 'darkorchid1', 
               # color = 'darkorchid3'
  ) +
  scale_x_discrete(limits = xorder_l) +
  # coord_cartesian(ylim = c(-5, 2)) +
  # scale_y_continuous(limits = quantile(plot_df$ATE, c(0.1, 0.9), na.rm=TRUE)) +
  # coord_cartesian(ylim = c(-8.329534,  4.251341 )) +
  facet_wrap(vars(AY_name), ncol=3, nrow = 3, scales = 'free_y') +
  labs(x = 'Estimation Method', y = 'ATE Estimate') +
  theme(axis.text.x = element_text(angle = 45, vjust = 1.2, hjust = 1.1, size = 6))
ggsave(sprintf('%s/individualAY/AY_allATE_boxplot_positive_numNC_more.pdf', plot_savepath), height = 8, width = 12)
ggsave(sprintf('%s/individualAY/AY_allATE_boxplot_positive_numNC_more.svg', plot_savepath), height = 8, width = 12)
ggsave(sprintf('%s/individualAY/AY_allATE_boxplot_positive_numNC_more.jpg', plot_savepath), height = 8, width = 12)


# AY ATE Estimates by NC (over ZW selections): Lines 
# average
# ggplot(df |> # filter(numNC == maxNC | is.na(numNC)) |> 
#          filter(type == 'positive' | type == 'negative') |> 
#          mutate(AY_idx = as.factor(AY_idx)) |>
#          group_by(type, AY_idx, EstMethod_l, method) |>
#          summarize(ATE_median = mean(ATE),
#                    .groups = 'drop'),
#        aes(x = EstMethod_l, y = ATE_median, 
#            color = method, fill = method,
#            group = AY_idx)) +
#   geom_hline(aes(yintercept = 0), alpha = 1, linetype = 'solid', color = 'black', linewidth = .7) +
#   geom_line(aes(color = AY_idx), alpha = .5, linewidth = .3) +
#   # geom_point(aes(color = method), alpha = .3, size = 1) +
#   coord_cartesian(ylim = c(-1, 1)) +
#   scale_x_discrete(limits = xorder_l,
#                    # labels = xorder_l,  
#                    breaks = xorder_l_breaks) +
#   facet_wrap(vars(type), ncol=1, scales = 'free_y',
#              labeller = type_labeller) +
#   labs(title = 'Mean ATE for AY pairs',
#        y = 'Mean ATE', x = 'Estimation Method') +
#   theme(legend.position = 'none',
#         axis.text.x = element_text(angle = 45, vjust = 1.2, hjust = 1.1, size = 10), 
#         panel.grid.major.x = element_line(color = 'gray'))
# ggsave(sprintf('%s/individualAY/AY_meanATE_lineplot_numNC_more.pdf', plot_savepath), height = 8, width =  12)
# ggsave(sprintf('%s/individualAY/AY_meanATE_lineplot_numNC_more.svg', plot_savepath), height = 8, width =  12)
# ggsave(sprintf('%s/individualAY/AY_meanATE_lineplot_numNC_more.jpg', plot_savepath), height = 8, width =  12)



# median
# ggplot(df |> # filter(numNC == maxNC | is.na(numNC)) |> 
#          filter(type == 'positive' | type == 'negative') |> 
#          mutate(AY_idx = as.factor(AY_idx)) |>
#          group_by(type, AY_idx, EstMethod_l, method) |>
#          summarize(ATE_median = median(ATE),
#                    .groups = 'drop'),
#        aes(x = EstMethod_l, y = ATE_median, 
#            color = method, fill = method,
#            group = AY_idx)) +
#   geom_hline(aes(yintercept = 0), alpha = 1, linetype = 'solid', color = 'black', linewidth = .7) +
#   geom_line(aes(color = AY_idx), alpha = .5, linewidth = .3) +
#   # geom_point(aes(color = method), alpha = .3, size = 1) +
#   coord_cartesian(ylim = c(-1, 1)) +
#   scale_x_discrete(limits = xorder_l,
#                    # labels = xorder_l,  
#                    breaks = xorder_l_breaks) +
#   facet_wrap(vars(type), ncol=1, scales = 'free_y',
#              labeller = type_labeller) +
#   labs(title = 'Median ATE for AY pairs',
#        y = 'Median ATE', x = 'Estimation Method') +
#   theme(legend.position = 'none',
#         axis.text.x = element_text(angle = 45, vjust = 1.2, hjust = 1.1, size = 10), 
#         panel.grid.major.x = element_line(color = 'gray'))
# ggsave(sprintf('%s/individualAY/AY_medianATE_lineplot_numNC.pdf', plot_savepath), height = 8, width = 12)
# ggsave(sprintf('%s/individualAY/AY_medianATE_lineplot_numNC.svg', plot_savepath), height = 8, width = 12)
# ggsave(sprintf('%s/individualAY/AY_medianATE_lineplot_numNC.jpg', plot_savepath), height = 8, width = 12)


# =================== Make Plots: ATE/se ======================================
print(sprintf("[%s]        - ATE/se", Sys.time()))
# Save in separate folder
plot_savepath_subfolder = sprintf('%s/ATEoverSE/', plot_savepath)
dir.create(plot_savepath_subfolder, showWarnings = FALSE)


# =================== Make Plots: Ratio (Prox/Naive) ======================================
print(sprintf("[%s]        - Ratio (Prox/Naive)", Sys.time()))
# Save in separate folder
plot_savepath_subfolder = sprintf('%s/numNC/', plot_savepath)
dir.create(plot_savepath_subfolder, showWarnings = FALSE)



# =================== Make Plots: max NCs AY pairs' ATEs ======================================
# print(sprintf("[%s]        - AY pairs' ATEs", Sys.time()))
# 
# 
# maxNC = max(df$numNC, na.rm = TRUE)
# # AY ATE Estimates (average over ZW selections): Boxplot by type
# # mean
# ggplot(df |> filter(numNC == maxNC | is.na(numNC)) |>  
#          group_by(type, AY_idx, EstMethod, method, basis) |> 
#          summarize(ATE = mean(ATE, na.rm = TRUE))) +
#   geom_hline(aes(yintercept = 0), alpha = 1, linetype = 'solid') +
#   geom_boxplot(aes(x = EstMethod, y = (ATE), alpha = basis, fill = method),
#                # fill = 'orange', alpha = .6, 
#                outlier.size = .1, outlier.alpha = .3) +
#   coord_cartesian(ylim = c(-2, 2)) +
#   scale_x_discrete(limits = xorder) +
#   scale_alpha_discrete(guide = 'none') +
#   facet_wrap(vars(type), nrow = 3,
#              labeller = type_labeller) +
#   labs(title = 'Mean ATE for AY pairs (w/ max NC)',
#        y = 'Mean ATE', x = 'Estimation Method') +
#   theme(axis.text.x = element_text(angle = 45, vjust = 1.2, hjust = 1.1, size = 10))
# ggsave(sprintf('%s/AY_meanATE_boxplot_maxNC.pdf', plot_savepath), height = 5, width = 8)
# ggsave(sprintf('%s/AY_meanATE_boxplot_maxNC.svg', plot_savepath), height = 5, width = 8)
# ggsave(sprintf('%s/AY_meanATE_boxplot_maxNC.jpg', plot_savepath), height = 5, width = 8)
# 
# # median
# ggplot(df |> filter(numNC == maxNC | is.na(numNC)) |> 
#          filter(EstMethod %in% c('lmYA', 'OCB2SLSpci2s', 'OCB2SLSReg', 'OCBMEstRWbasis1', 'OCBMEstRWRegbasis1')) |>
#          group_by(type, AY_idx, EstMethod, method, basis) |> 
#          summarize(ATE = median(ATE, na.rm = TRUE))) +
#   geom_hline(aes(yintercept = 0), alpha = 1, linetype = 'solid', color = 'gray') +
#   geom_boxplot(aes(x = EstMethod, y = (ATE), 
#                    alpha = basis, fill = method),
#                # fill = 'orange', alpha = .6, 
#                outlier.size = .1, outlier.alpha = .3) +
#   coord_cartesian(ylim = c(-1.5, .5)) +
#   # coord_cartesian(ylim = c(-10, 10)) +
#   # scale_x_discrete(limits = xorder) +
#   scale_alpha_discrete(guide = 'none') +
#   facet_wrap(vars(type), nrow = 3,
#              labeller = type_labeller) +
#   # labs(title = 'Median ATE for AY pairs (w/ max NC)',
#   #      y = 'Median ATE', x = 'Estimation Method') +
#   labs(title = 'ATE for AY pairs (w/ 20 NC Pairs)',
#        y = 'ATE', x = 'Estimation Method') +
#   theme(axis.text.x = element_text(angle = 45, vjust = 1.1, hjust = 1.2, size = 14), 
#         legend.position = 'none')
# ggsave(sprintf('%s/AY_medianATE_boxplot_maxNC.pdf', plot_savepath), height = 6, width = 8)
# ggsave(sprintf('%s/AY_medianATE_boxplot_maxNC.svg', plot_savepath), height = 6, width = 8)
# ggsave(sprintf('%s/AY_medianATE_boxplot_maxNC.jpg', plot_savepath), height = 6, width = 8)
# 
# 
# 
# 
# 
# 
# 
# 
# # AY ATE Estimates (average over ZW selections): Lines 
# # median
# ggplot(df |> filter(numNC == maxNC | is.na(numNC)) |> 
#          filter(type == 'positive' | type == 'negative') |> 
#          group_by(type, AY_idx, EstMethod, method, basis) |> 
#          mutate(AY_idx = as.factor(AY_idx)) |>
#          summarize(ATE_median = median(ATE),
#                    .groups = 'drop'),
#        aes(x = EstMethod,  y = ATE_median,
#            color = type,
#            fill = type,
#            group = AY_idx)) +
#   geom_hline(aes(yintercept = 0), alpha = 1, linetype = 'solid', color = 'black', linewidth = .7) +
#   geom_line(alpha = .5, linewidth = .3) +
#   coord_cartesian(ylim = c(-1, 1)) +
#   scale_x_discrete(limits = xorder) +
#   facet_wrap(vars(type), ncol=1, scales = 'free_y',
#              labeller = type_labeller) +
#   labs(title = 'Median ATE for AY pairs',
#        y = 'Median ATE', x = 'Estimation Method') +
#   theme(legend.position = 'none',
#         axis.text.x = element_text(angle = 45, vjust = 1.2, hjust = 1.1, size = 10))
# ggsave(sprintf('%s/AY_medianATE_lineplot.pdf', plot_savepath), height = 8, width = 8)
# ggsave(sprintf('%s/AY_medianATE_lineplot.svg', plot_savepath), height = 8, width = 8)
# ggsave(sprintf('%s/AY_medianATE_lineplot.jpg', plot_savepath), height = 8, width = 8)
# 





# # =================== Make Plots: AY pairs' Variances ======================================
# print(sprintf("[%s]        - AY pairs' Variances", Sys.time()))
# 
# 
# 
# # variance/spread over ZW selections
# ggplot(df |> group_by(type, AY_idx, EstMethod) |> 
#          summarize(ATEvar = var(ATE, na.rm = TRUE),
#                    .groups = 'drop'),
#        # summarize(ATEvar = IQR(ATE, na.rm = TRUE)),
#        aes(x = EstMethod, 
#            y = ATEvar |> log(), 
#            group = type, 
#            color = type, 
#            fill  = type)) +
#   geom_point(alpha = .5, size = .3) +
#   geom_smooth(linewidth = 1, method = 'loess', span = .2) +
#   scale_x_discrete(c('lmYA', paste0('CB', 1:30))) +
#   # scale_y_continuous(limits = c(0, 1)) +
#   coord_cartesian(ylim = c(-0, 10)) +
#   scale_color_manual(values = c('darkorchid3',  'seagreen')) +
#   scale_fill_manual(values = c('darkorchid3',  'seagreen')) +
#   labs(title = 'Smoothed Variance of ATE Estimates', 
#        x = 'Estimation Method', 
#        y = 'log(Variance of ATE Estimates)')  +
#   theme(axis.text.x = element_text(angle = 45, vjust = 1.2, hjust = 1.1, size = 10),
#         # legend.position = 'top', 
#         legend.position = c(.9, .02),
#         legend.justification=c(.5, 0))
# ggsave(sprintf('%s/AY_varATE_smoothline.pdf', plot_savepath), height = 8, width = 8)
# ggsave(sprintf('%s/AY_varATE_smoothline.svg', plot_savepath), height = 8, width = 8)
# ggsave(sprintf('%s/AY_varATE_smoothline.jpg', plot_savepath), height = 8, width = 8)
# 
# # IQR/spread over ZW selections
# ggplot(df |> group_by(type, AY_idx, EstMethod) |> 
#          summarize(ATEiqr = IQR(ATE, na.rm = TRUE),
#                    .groups = 'drop'),
#        aes(x = EstMethod, 
#            y = ATEiqr |> log(), 
#            group = type, 
#            color = type, 
#            fill  = type)) +
#   geom_point(alpha = .5, size = .3) +
#   geom_smooth(linewidth = 1, method = 'loess', span = .2) +
#   scale_x_discrete(c('lmYA', paste0('CB', 1:30))) +
#   # coord_cartesian(ylim = c(-0, 10)) +
#   scale_color_manual(values = c('darkorchid3',  'seagreen')) +
#   scale_fill_manual(values = c('darkorchid3',  'seagreen')) +
#   labs(title = 'Smoothed IQR of ATE Estimates', 
#        x = 'Estimation Method', 
#        y = 'log(IQR of ATE Estimates)')  +
#   theme(axis.text.x = element_text(angle = 45, vjust = 1.2, hjust = 1.1, size = 10),
#         # legend.position = 'top', 
#         legend.position = c(.9, .02),
#         legend.justification=c(.5, 0))
# ggsave(sprintf('%s/AY_iqrATE_smoothline.pdf', plot_savepath), height = 8, width = 8)
# ggsave(sprintf('%s/AY_iqrATE_smoothline.svg', plot_savepath), height = 8, width = 8)
# ggsave(sprintf('%s/AY_iqrATE_smoothline.jpg', plot_savepath), height = 8, width = 8)





# =================== Make Plots: Ratio to Biased lmYA Estimate ================
# save ATE ratio in a separate folder
# dir.create(sprintf('%s/ATEratio/', plot_savepath))
# 
# # look at the ratio between new ATE estimate and biased lmYA estimate
# # df_lmYA = df |>
# #   filter(method == 'lmYA') |> 
# #   select(A, Y, type, ATE) |>
# #   distinct() |>
# #   rename(ATElmYA = ATE)
# 
# df_ratio = merge(df, 
#                  df_lmYA |> 
#                    select(A, Y, ATE) |>
#                    rename(ATElmYA = ATE)) |> 
#   mutate(ATEratio = ATE / ATElmYA)
# 
# 
# # AY ATE Ratio: Lines 
# # negative
# ggplot(# df_ratio[1:1000, ],
#   # aes(x = EstMethod_l, y = ATEratio,
#   #     color = ATElmYA, 
#   #     # fill = type,
#   #     group = interaction(AY_idx, ZW_idx))
#   df_ratio |>
#     # filter(abs(ATElmYA) > .5) |> # will filter out in outer call
#     filter(numNC > 11 | EstMethod_l == 'lmYA') |>
#     group_by(type, AY_idx, ATElmYA,  EstMethod_l) |>
#     summarize(ATEratio = median(ATEratio), .groups = 'drop') |> 
#     filter(type == 'negative'),
#   aes(x = EstMethod_l, y = ATEratio,
#       color = log(abs(ATElmYA)), 
#       # fill = type,
#       group = interaction(AY_idx))
# ) +
#   geom_hline(aes(yintercept = 0), alpha = 1, linetype = 'solid', color = 'black', linewidth = .7) +
#   geom_line(alpha = .7, linewidth = .3) +
#   # coord_cartesian(ylim = c(-1, 1.5)) +
#   coord_cartesian(ylim = c(-1.5, 2)) +
#   # scale_x_discrete(limits = xorder_l) +
#   viridis::scale_color_viridis() +
#   facet_wrap(vars(type), ncol=1, scales = 'free_y',
#              labeller = type_labeller) +
#   labs(title = 'ATE Ratio for AY pairs',
#        y = 'ATE Ratio', x = 'Estimation Method') +
#   theme(legend.position = 'right',
#         axis.text.x = element_text(angle = 45, vjust = 1.2, hjust = 1.1, size = 10))
# ggsave(sprintf('%s/ATEratio/AY_ATEratio_lineplot_negative.pdf', plot_savepath), height = 8, width = 8)
# ggsave(sprintf('%s/ATEratio/AY_ATEratio_lineplot_negative.svg', plot_savepath), height = 8, width = 8)
# ggsave(sprintf('%s/ATEratio/AY_ATEratio_lineplot_negative.jpg', plot_savepath), height = 8, width = 8)
# 
# 
# # positive
# ggplot(# df_ratio[1:1000, ],
#   # aes(x = EstMethod_l, y = ATEratio,
#   #     color = ATElmYA, 
#   #     # fill = type,
#   #     group = interaction(AY_idx, ZW_idx))
#   df_ratio |>
#     # filter(abs(ATElmYA) > .5) |> # will filter out in outer call
#     filter(numNC > 11 | EstMethod_l == 'lmYA') |>
#     group_by(type, AY_idx, ATElmYA,  EstMethod_l) |>
#     summarize(ATEratio = median(ATEratio), .groups = 'drop') |> 
#     filter(type == 'positive'),
#   aes(x = EstMethod_l, y = ATEratio,
#       color = log(abs(ATElmYA)), 
#       # fill = type,
#       group = interaction(AY_idx))
# ) +
#   geom_hline(aes(yintercept = 0), alpha = 1, linetype = 'solid', color = 'black', linewidth = .7) +
#   geom_line(alpha = .7, linewidth = .3) +
#   # coord_cartesian(ylim = c(-1, 1.5)) +
#   coord_cartesian(ylim = c(-1.5, 2)) +
#   # scale_x_discrete(limits = xorder_l) +
#   viridis::scale_color_viridis() +
#   facet_wrap(vars(type), ncol=1, scales = 'free_y',
#              labeller = type_labeller) +
#   labs(title = 'ATE Ratio for AY pairs',
#        y = 'ATE Ratio', x = 'Estimation Method') +
#   theme(legend.position = 'right',
#         axis.text.x = element_text(angle = 45, vjust = 1.2, hjust = 1.1, size = 10))
# ggsave(sprintf('%s/ATEratio/AY_ATEratio_lineplot_positive.pdf', plot_savepath), height = 8, width = 8)
# ggsave(sprintf('%s/ATEratio/AY_ATEratio_lineplot_positive.svg', plot_savepath), height = 8, width = 8)
# ggsave(sprintf('%s/ATEratio/AY_ATEratio_lineplot_positive.jpg', plot_savepath), height = 8, width = 8)
# 
# 


# # Histogram of Biased lmYA estimates
# ggplot(df_lmYA,
#        aes(x = ATElmYA, y = after_stat(density))) +
#   geom_histogram(fill = 'orange', alpha = .6, binwidth = .1) +
#   geom_density(color = 'orange') +
#   geom_vline(aes(xintercept = 0), color = 'darkorange') +
#   scale_x_continuous(limits = c(-4, 4)) +
#   facet_wrap(vars(type), nrow = 3, labeller = type_labeller) +
#   # facet_grid(type ~ . , labeller = type_labeller) +
#   labs(title = 'Histogram of lmYA ATE Estimates',  
#        x = 'ATE Estimate') 
# 
# ggplot(df_lmYA,
#        aes(x = ATElmYA |> abs(), y = after_stat(density))) +
#   geom_histogram(fill = 'orange', alpha = .6, binwidth = .05) +
#   geom_density(color = 'orange') +
#   geom_vline(aes(xintercept = 0), color = 'darkorange') +
#   scale_x_continuous(limits = c(-.1, 4)) +
#   facet_wrap(vars(type), nrow = 3, labeller = type_labeller) +
#   # facet_grid(type ~ . , labeller = type_labeller) +
#   labs(title = 'Histogram of lmYA ATE Estimates',  
#        x = 'ATE Estimate') 


# =================== Make Plots: Median AY ATEs by numNC  =======================
# # (will just be the single estimate when only one ZW set per AY pair e.g. SPCA)
# 
# # Median ATE by numNC (boxplot) 
# ggplot(df |> # filter(numNC == maxNC | is.na(numNC)) |> 
#          group_by(type, AY_idx, EstMethod_l, method, basis) |> 
#          summarize(ATE = median(ATE, na.rm = TRUE))) +
#   geom_hline(aes(yintercept = 0), alpha = 1, linetype = 'solid', color = 'gray') +
#   geom_boxplot(aes(x = EstMethod_l, y = (ATE), 
#                    alpha = basis, fill = method),
#                # fill = 'orange', alpha = .6, 
#                outlier.size = .1, outlier.alpha = .3) +
#   coord_cartesian(ylim = c(-1.5, 1)) +
#   # coord_cartesian(ylim = c(-10, 10)) +
#   # scale_x_discrete(limits = xorder) +
#   scale_alpha_discrete(guide = 'none') +
#   facet_wrap(vars(type), nrow = 3,
#              labeller = type_labeller) +
#   labs(title = 'Median ATE for AY pairs',
#        y = 'Median ATE', x = 'Estimation Method') +
#   theme(axis.text.x = element_text(angle = 45, vjust = 1.1, hjust = 1.2, size = 7), 
#         legend.position = 'none')
# 
# ggsave(sprintf('%s/AY_medianATE_boxplot_numNC.pdf', plot_savepath), height = 6, width = 8)
# ggsave(sprintf('%s/AY_medianATE_boxplot_numNC.svg', plot_savepath), height = 6, width = 8)
# ggsave(sprintf('%s/AY_medianATE_boxplot_numNC.jpg', plot_savepath), height = 6, width = 8)
# 
# # Median ATE by numNC (violinplot)
# ggplot(df |> # filter(numNC == maxNC | is.na(numNC)) |> 
#          group_by(type, AY_idx, EstMethod_l, method, basis) |> 
#          summarize(ATE = median(ATE, na.rm = TRUE))) +
#   geom_hline(aes(yintercept = 0), alpha = 1, linetype = 'solid', color = 'gray') +
#   geom_violin(aes(x = EstMethod_l, y = (ATE), 
#                    alpha = basis, fill = method),
#                   color = 'black', scale = 'count') +
#   coord_cartesian(ylim = c(-1.5, 1)) +
#   # coord_cartesian(ylim = c(-10, 10)) +
#   # scale_x_discrete(limits = xorder) +
#   scale_alpha_discrete(guide = 'none') +
#   facet_wrap(vars(type), nrow = 3,
#              labeller = type_labeller) +
#   labs(title = 'Median ATE for AY pairs (w/ max NC)',
#        y = 'Median ATE', x = 'Estimation Method') +
#   theme(axis.text.x = element_text(angle = 45, vjust = 1.1, hjust = 1.2, size = 7), 
#         legend.position = 'none') 
# 
# ggsave(sprintf('%s/AY_medianATE_violinplot_numNC.pdf', plot_savepath), height = 6, width = 8)
# ggsave(sprintf('%s/AY_medianATE_violinplot_numNC.svg', plot_savepath), height = 6, width = 8)
# ggsave(sprintf('%s/AY_medianATE_violinplot_numNC.jpg', plot_savepath), height = 6, width = 8)
# 
# 
# # =================== Make Plots: by NCs ATEs ======================================
# print(sprintf("[%s]        - all ATEs", Sys.time()))
# 
# 
# 
# # median
# ggplot(df |> filter(numNC == maxNC | is.na(numNC)) |> 
#          group_by(type, AY_idx, EstMethod, method, basis) |> 
#          summarize(ATE = median(ATE, na.rm = TRUE))) +
#   geom_hline(aes(yintercept = 0), alpha = 1, linetype = 'solid', color = 'gray') +
#   geom_boxplot(aes(x = EstMethod, y = (ATE), 
#                    alpha = basis, fill = method),
#                # fill = 'orange', alpha = .6, 
#                outlier.size = .1, outlier.alpha = .3) +
#   coord_cartesian(ylim = c(-1.5, 1)) +
#   # coord_cartesian(ylim = c(-10, 10)) +
#   scale_x_discrete(limits = xorder) +
#   scale_alpha_discrete(guide = 'none') +
#   facet_wrap(vars(type), nrow = 3,
#              labeller = type_labeller) +
#   labs(title = 'Median ATE for AY pairs (w/ max NC)',
#        y = 'Median ATE', x = 'Estimation Method') +
#   theme(axis.text.x = element_text(angle = 45, vjust = 1.1, hjust = 1.2, size = 7), 
#         legend.position = 'none')
# ggsave(sprintf('%s/AY_medianATE_boxplot_maxNC.pdf', plot_savepath), height = 6, width = 8)
# ggsave(sprintf('%s/AY_medianATE_boxplot_maxNC.svg', plot_savepath), height = 6, width = 8)
# ggsave(sprintf('%s/AY_medianATE_boxplot_maxNC.jpg', plot_savepath), height = 6, width = 8)
# 
# for(numNC_cur in unique(df$numNC)[!is.na(unique(df$numNC))]) {
#   print(sprintf('%s/numNC=%s/ATE_hist_maxNC.pdf', plot_savepath, numNC_cur))
#   dir.create(sprintf('%s/numNC=%s/', plot_savepath, numNC_cur))
#   
#   # all ATE Estimates: Histogram by type
#   ggplot(df |> filter(numNC == numNC_cur),
#          aes(x = ATE, y = after_stat(density))) +
#     geom_histogram(fill = 'orange', alpha = .6, binwidth = .1) +
#     geom_density(color = 'orange') +
#     geom_vline(aes(xintercept = 0), color = 'darkorange') +
#     scale_x_continuous(limits = c(-4, 4)) +
#     facet_wrap(vars(type), nrow = 3, labeller = type_labeller) +
#     # facet_grid(type ~ . , labeller = type_labeller) +
#     labs(title = sprintf('Histogram of All ATE Estimates w/ numNC = %s', numNC_cur),  
#          x = 'ATE Estimate') 
#   
#   ggsave(sprintf('%s/numNC=%s/ATE_hist_maxNC.pdf', plot_savepath, numNC_cur), height = 5, width = 5)
#   ggsave(sprintf('%s/numNC=%s/ATE_hist_maxNC.svg', plot_savepath, numNC_cur), height = 5, width = 5)
#   ggsave(sprintf('%s/numNC=%s/ATE_hist_maxNC.jpg', plot_savepath, numNC_cur), height = 5, width = 5)
#   
#   # all ATE Estimates: Boxplot by type
#   # median
#   ggplot(df |> filter(numNC == numNC_cur | is.na(numNC)) |> 
#            group_by(type, AY_idx, EstMethod, method, basis) |> 
#            summarize(ATE = median(ATE, na.rm = TRUE))) +
#     geom_hline(aes(yintercept = 0), alpha = 1, linetype = 'solid') +
#     geom_boxplot(aes(x = EstMethod, y = ATE, 
#                      # pattern_fill  = basis,
#                      alpha = basis, 
#                      fill = method),
#                  # fill = 'orange', 
#                  # alpha = .6,
#                  linewidth = .1,
#                  outlier.shape = NA, # don't show outliers
#                  outlier.size = .1, outlier.alpha = .3) +
#     coord_cartesian(ylim = ylim*2) +
#     scale_x_discrete(limits = xorder) +
#     scale_alpha_discrete(guide = 'none') +
#     facet_wrap(vars(type), nrow = 3) +
#     labs(title = sprintf('Median ATE Estimates (w/ numNC=%s)', numNC_cur),
#          x = 'Estimation Method', y = 'ATE Estimate') +
#     theme(axis.text.x = element_text(angle = 45, vjust = 1.2, hjust = 1.3, size = 10))
#   ggsave(sprintf('%s/numNC=%s/medianATE_boxplot.pdf', plot_savepath, numNC_cur), height = 5, width = 8)
#   ggsave(sprintf('%s/numNC=%s/medianATE_boxplot.svg', plot_savepath, numNC_cur), height = 5, width = 8)
#   ggsave(sprintf('%s/numNC=%s/medianATE_boxplot.jpg', plot_savepath, numNC_cur), height = 5, width = 8)
#   
#   
#   # all ATE Estimates: Boxplot by type
#   ggplot(df |> filter(numNC == numNC_cur | is.na(numNC)) |> 
#            group_by(type, AY_idx, EstMethod, method, basis) |> 
#            summarize(ATE = median(ATE, na.rm = TRUE))) +
#     geom_hline(aes(yintercept = 0), alpha = 1, linetype = 'solid') +
#     geom_boxplot(aes(x = EstMethod, y = ATE, 
#                      # pattern_fill  = basis,
#                      alpha = basis, 
#                      fill = method),
#                  # fill = 'orange', 
#                  # alpha = .6,
#                  linewidth = .1,
#                  outlier.shape = NA, # don't show outliers
#                  outlier.size = .1, outlier.alpha = .3) +
#     coord_cartesian(ylim = ylim*2) +
#     scale_x_discrete(limits = xorder[-1]) +
#     scale_alpha_discrete(guide = 'none') +
#     facet_wrap(vars(type), nrow = 3) +
#     labs(title = sprintf('Median ATE Estimates (w/ numNC=%s)', numNC_cur),
#          x = 'Estimation Method', y = 'ATE Estimate') +
#     theme(axis.text.x = element_text(angle = 45, vjust = 1.2, hjust = 1.3, size = 10))
#   ggsave(sprintf('%s/numNC=%s/medianATE_boxplot_nolmYA.pdf', plot_savepath, numNC_cur), height = 5, width = 8)
#   ggsave(sprintf('%s/numNC=%s/medianATE_boxplot_nolmYA.svg', plot_savepath, numNC_cur), height = 5, width = 8)
#   ggsave(sprintf('%s/numNC=%s/medianATE_boxplot_nolmYA.jpg', plot_savepath, numNC_cur), height = 5, width = 8)
#   
# }


# =================== Make Plots: looking at pvalues ======================================
# # folder for p-value plots to be saved in 
# dir.create(sprintf('%s/pvals/', plot_savepath), showWarnings = FALSE)
# 
# # ---- p-values from z-tests on each AY pair means of ZW choices, estimate variance from
# # all means for particular estimation method setting (e.g. OCB2SLSReg with 5 NCs)
# 
# # using estimated standard deviation of the estimates for non-causal AY for each EstMethod_l
# # (Approx with a normal... but could just use the quantiles... but then REALLY need to use separate sample)
# EstMethod_sd = df |> 
#   filter(type == 'negative') |>
#   # filter(EstMethod_l != 'lmYA') |>
#   group_by(AY_idx, EstMethod_l) |> # first find ATEmean for each AY, EstMethod_l 
#   summarize(ATEmean = mean(ATE),
#             .groups = 'drop') |>
#   group_by(EstMethod_l) |> 
#   summarize(EstMethod_l_sd = sd(ATEmean))
# 
# # dim(EstMethod_sd); head(EstMethod_sd)
# 
# # should be faster access to sd than filtering
# EstMethod_sd_list = setNames(EstMethod_sd$EstMethod_l_sd, EstMethod_sd$EstMethod_l)
# mypvalfunc <- function(x, m)  {
#   # m =  'OCB2SLS10'; x = -.02
#   # pnorm(0, mean = 0, sd = .03)
#   pnorm(q=x, 
#         mean = 0, 
#         sd = EstMethod_sd_list[m])
#         # sd = EstMethod_sd |> filter(EstMethod_l == m) |> pull(EstMethod_l_sd))
# }
# 
# df_ATEmean = df |> 
#   group_by(AY_idx, type, numNC, method, EstMethod_l) |>
#   summarize(ATEmean = mean(ATE),
#             .groups = 'drop')
# 
# # dim(df_ATEmean); head(df_ATEmean)
# # mapply(FUN = mypvalfunc, 
# #        x = c(-1, 0, .2),
# #        m = c('OCB2SLS10', 'OCB2SLS10', 'OCB2SLS50'))
# df_ATEmean$zpvals = mapply(FUN = mypvalfunc, 
#        x = df_ATEmean$ATEmean,
#        m = df_ATEmean$EstMethod_l)
# 
# # head(df_ATEmean)
# # df_ATEmean |> filter(is.na(zpvals))
# # ggplot(df_ATEmean, 
# #        aes(x = zpvals)) +
# #   geom_histogram() +
# #   facet_wrap(vars(type), ncol = 1)
# # df_ATEmean$zpvals |> hist()
# 
# 
# # # z pvalues violin plot per EstMethod
# # ggplot(df_ATEmean |> filter(numNC > 0 | is.na(numNC)),
# #        aes(x = EstMethod_l, y = zpvals, 
# #            fill = method, color = method)) +
# #   geom_hline(aes(yintercept = 0.05), alpha = .8, linetype = 'solid', color = 'red') +
# #   geom_jitter(shape=16, position=position_jitter(0.2),
# #               alpha = .6, 
# #               # color = 'darkorange2', 
# #               size = 1) +
# #   geom_violin(alpha = .3, 
# #               color = 'black',
# #               # fill = 'orange',
# #               draw_quantiles = c(.05),
# #               trim = TRUE,
# #               bounds = c(0, 1),
# #               scale = TRUE) +
# #   # geom_boxplot(# fill = 'orange', 
# #   #              # alpha = .6,
# #   #              linewidth = .1,
# #   #              outlier.shape = NA, # don't show outliers
# #   #              outlier.size = .1, outlier.alpha = .3) +
# #   # coord_cartesian(ylim = ylim*2) +
# #   scale_x_discrete(limits = xorder_l,
# #                    breaks = xorder_l_breaks) +
# #   # scale_alpha_discrete(guide = 'none') +
# #   facet_wrap(vars(type), nrow = 3) +
# #   labs(title = 'Violin plot of z-test pvals (sd est neg AY pairs)',
# #        x = 'Estimation Method', y = 'pval') +
# #   theme(axis.text.x = element_text(angle = 45, vjust = 1.2, hjust = 1.3, size = 10))
# # ggsave(sprintf('%s/pvals/zpvals_violinplot.pdf', plot_savepath), height = 10, width = 16)
# # ggsave(sprintf('%s/pvals/zpvals_violinplot.svg', plot_savepath), height = 10, width = 16)
# # ggsave(sprintf('%s/pvals/zpvals_violinplot.jpg', plot_savepath), height = 10, width = 16)
# 
# 
# 
# # Include the lmYA pvalues (pvalues from linear model of Y ~ A)
# df_lmYA$numNC = NA
# df_lmYA$method = 'lmYA0'
# df_lmYA$EstMethod_l = 'lmYA0'
# 
# 
# df_pvals = dplyr::bind_rows(df_ATEmean |> dplyr::rename(pval = zpvals), 
#                             df_lmYA |> dplyr::select(-A, -A_chr, -Y, -Y_chr, -ZW_idx, -ATE))
# head(df_pvals)
# 
# # z pvalues violin plot per EstMethod incl lmYA pvals 
# ggplot(df_pvals |> filter(numNC > 0 | is.na(numNC)),
#        aes(x = EstMethod_l, y = pval, 
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
#   scale_x_discrete(limits = c('lmYA0', xorder_l), #) + #,
#                    breaks = c('lmYA0', xorder_l_breaks)) +
#   # scale_alpha_discrete(guide = 'none') +
#   facet_wrap(vars(type), nrow = 3) +
#   labs(title = 'Violin plot of lmYA pval/z-test pvals (sd est neg AY pairs)',
#        x = 'Estimation Method', y = 'pval') +
#   theme(axis.text.x = element_text(angle = 45, vjust = 1.2, hjust = 1.3, size = 10))
# ggsave(sprintf('%s/pvals/zpvals_violinplot.pdf', plot_savepath), height = 10, width = 20)
# ggsave(sprintf('%s/pvals/zpvals_violinplot.svg', plot_savepath), height = 10, width = 20)
# ggsave(sprintf('%s/pvals/zpvals_violinplot.jpg', plot_savepath), height = 10, width = 20)
# 
# 
# # ---- p-values from t-tests on each AY pair over ZW choices 
# # e.g. for A5 --> Y4, there were 5 sets of ZW choices. 
# # Perform t-test on these 5 estimates, to test against ATE=0
# # No, the estimates for each AY is too correlated, so underestimates the variance
# # simple way to get a quick 'p-value'? 
# # t-test
# myttestfunc <- function(x) {
#   if(length(x) <= 1 | all(is.na(x))) {
#     return(NA)
#   }
#   t.test(x, mu=0)$p.value
# }
# # res = t.test(3:10, mu=0)
# 
# ttest_pval = df |> 
#   # filter(type == 'negative') |>
#   filter(EstMethod_l != 'lmYA') |>
#   group_by(AY_idx, EstMethod_l) |>
#   summarize(pval = myttestfunc(ATE),
#             ATEmean = mean(ATE),
#             ATEmedian = median(ATE),
#             ATEsd = sd(ATE),
#             count = n(),
#             .groups = 'drop')
# 
# 
# # head(ttest_pval)
# 
# df_tpval = merge(df |> select(AY_idx, type, A, A_chr, Y, Y_chr, 
#                              method, method_type, numNC, basis, 
#                              EstMethod, EstMethod_l) |> distinct(),
#                 ttest_pval)
# # add lmYA0 pvalues (from linear model)
# df_tpval = dplyr::bind_rows(df_tpval,
#                             df_lmYA |> dplyr::select(-A, -A_chr, -Y, -Y_chr, -ZW_idx, -ATE))
# # head(df_tpval)
# # dim(df_tpval); dim(df)
# # df |> filter(AY_idx == 2 & EstMethod_l == 'OCB2SLS1')
# # hist(ttest_pval$pval)
# 
# ggplot(df_tpval, 
#        aes(x = EstMethod_l, y = pval,
#            fill = method, color = method)) +
#   geom_hline(aes(yintercept = 0), alpha = 1, linetype = 'solid') +
#   geom_jitter(shape=16, position=position_jitter(0.2),
#               alpha = .6, 
#               # color = 'darkorange2', 
#               size = 1) +
#   geom_violin(alpha = .3, 
#               color = 'black',
#               # fill = 'orange',
#               # draw_quantiles = c(.05),
#               trim = TRUE,
#               # bounds = c(0, 1),
#               scale = TRUE) +
#   # geom_boxplot(aes( 
#   #                  # pattern_fill  = basis,
#   #                  alpha = basis, 
#   #                  fill = method),
#   #              # fill = 'orange', 
#   #              # alpha = .6,
#   #              linewidth = .1,
#   #              # outlier.shape = NA, # don't show outliers
#   #              outlier.size = 1, outlier.alpha = 1) +
#   # coord_cartesian(ylim = ylim*2) +
#   # scale_y_log10() +
#   scale_x_discrete(limits = c('lmYA0', xorder_l), #) + #,
#                    breaks = c('lmYA0', xorder_l_breaks)) +
#   # scale_y_continuous(limits = c(0, 1)) +
#   # scale_alpha_discrete(guide = 'none') +
#   facet_wrap(vars(type), nrow = 3) +
#   labs(title = 'Violin plot of lmYA pval/t-test pvals',
#        x = 'Estimation Method', y = 'pval') +
#   theme(axis.text.x = element_text(angle = 45, vjust = 1.2, hjust = 1.3, size = 10))
# ggsave(sprintf('%s/pvals/tpvals_violinplot.pdf', plot_savepath), height = 10, width = 20)
# ggsave(sprintf('%s/pvals/tpvals_violinplot.svg', plot_savepath), height = 10, width = 20)
# ggsave(sprintf('%s/pvals/tpvals_violinplot.jpg', plot_savepath), height = 10, width = 20)
# 


# numNC_cur = 5
# df_tpval |> filter(method == 'OCB2SLSReg') |> filter(numNC == numNC_cur) |> pull(pval) |> hist()
# df_tpval |> filter(method == 'OCB2SLSReg') |> filter(numNC == numNC_cur) |> pull(ATEmean) |> hist()
# df_tpval |> filter(method == 'OCB2SLSReg') |> filter(numNC == numNC_cur) |> pull(ATEsd) |> hist()
# df_tpval |> filter(method == 'OCB2SLSReg') |> filter(numNC == numNC_cur) |> pull(count) |> hist(breaks = 1:10)




dev.off()


# =================== END ======================================================
print(sprintf("[%s] END", Sys.time()))



# ############################################################
# ggplot(df, 
#        aes(x = EstMethod, 
#            y = abs(ATE), 
#            group = interaction(AY_idx, ZW_idx), 
#            color = type)) +
#   geom_line() +
#   scale_x_discrete(c('lmYA', paste0('CB', 1:15))) +
#   scale_y_continuous(limits = c(-10, 10)) +
#   scale_color_manual(values = c('orange', 'red', 'green'))
# 
# 
# head(df)
# 
# # average over ZW selections
# ggplot(df |> group_by(type, AY_idx, EstMethod) |> 
#          summarize(ATE = mean(ATE, na.rm = TRUE)), 
#        aes(x = EstMethod, 
#            y = (ATE), 
#            group = AY_idx, 
#            fill = type,
#            color = type)) +
#   geom_hline(aes(yintercept = 0), alpha = 1, linetype = 'dashed') +
#   geom_line(alpha = .3) +
#   scale_x_discrete(c('lmYA', paste0('CB', 1:30))) +
#   coord_cartesian(ylim = c(-10, 10)) +
#   scale_color_manual(values = c('orange', 'red', 'green'))
# 
# 
# 
# 
# 
# 
# # average over YA and ZW selections
# ggplot(df |> group_by(type, EstMethod) |> 
#          summarize(ATE_mean= mean(ATE, na.rm = T),
#                    ATE_se  =   sd(ATE, na.rm = T)/sqrt(n()),
#                    .groups = 'drop') |>
#          rename(ATE = ATE_mean),
#        aes(x = EstMethod, y = ATE, group = type, color = type,
#            ymin = ATE - 2*ATE_se, 
#            ymax = ATE + 2*ATE_se)) +
#   geom_hline(aes(yintercept = 0), alpha = 1, linetype = 'dashed') +
#   geom_crossbar(orientation = 'x', alpha = .4, position = 'dodge', width = .5) +
#   # geom_line(linewidth = .5) +
#   scale_x_discrete(c('lmYA', paste0('CB', 1:15))) +
#   coord_cartesian(ylim = c(-8, 6)) +
#   scale_color_manual(values = c('orange', 'red', 'green')) +
#   labs(title = 'Average ATE by EstMethod')
# 
# 
# 
# 
# # average over ZW selections
# ggplot(df |> group_by(type, AY_idx, EstMethod) |> 
#          summarize(ATE = median(ATE, na.rm = TRUE)), 
#        aes(x = EstMethod, 
#            y = (ATE), 
#            group = type, 
#            color = type, 
#            fill  = type)) +
#   geom_hline(aes(yintercept = 0), alpha = 1, linetype = 'dashed') +
#   geom_point(size = .5, alpha = .2) +
#   geom_smooth(linewidth = 1, alpha = .1) +
#   scale_x_discrete(c('lmYA', paste0('CB', 1:15))) +
#   coord_cartesian(ylim = c(-10, 5)) +
#   scale_color_manual(values = c('orange', 'red', 'green')) +
#   scale_fill_manual(values = c('orange', 'red', 'green')) +
#   labs(title = 'Smoothed ATE of Estimators')
# 
# # variance/spread over ZW selections
# ggplot(df |> group_by(type, AY_idx, EstMethod) |> 
#          summarize(ATEvar = var(ATE, na.rm = TRUE)),
#          # summarize(ATEvar = IQR(ATE, na.rm = TRUE)),
#        aes(x = EstMethod, 
#            y = ATEvar |> log(), 
#            group = type, 
#            color = type, 
#            fill  = type)) +
#   geom_point(alpha = .5, size = .5) +
#   geom_smooth(linewidth = 1) +
#   scale_x_discrete(c('lmYA', paste0('CB', 1:15))) +
#   # scale_y_continuous(limits = c(0, 1)) +
#   scale_color_manual(values = c('orange', 'red', 'green')) +
#   scale_fill_manual(values = c('orange', 'red', 'green')) +
#   labs(title = 'Smoothed Variance of Estimators')
# 
# 
# # IQR/spread over ZW selections
# ggplot(df |> group_by(type, AY_idx, EstMethod) |> 
#          # summarize(ATEvar = var(ATE, na.rm = TRUE)),
#          summarize(ATEiqr = IQR(ATE, na.rm = TRUE)),
#        aes(x = EstMethod, 
#            y = ATEiqr,# |> log(), 
#            group = type, 
#            color = type, 
#            fill  = type)) +
#   geom_point(alpha = .5, size = .5) +
#   geom_smooth(linewidth = 1) +
#   scale_x_discrete(c('lmYA', paste0('CB', 1:15))) +
#   # scale_y_continuous(limits = c(0, 1)) +
#   scale_color_manual(values = c('orange', 'red', 'green')) +
#   scale_fill_manual(values = c('orange', 'red', 'green')) +
#   labs(title = 'Smoothed IQR of Estimators')
# 
# 
# # variance over AY and ZW selections
# ggplot(df |> group_by(type, EstMethod) |> 
#          summarize(ATEvar = var(ATE, na.rm = TRUE)), 
#        aes(x = EstMethod, 
#            y = ATEvar |> log(), 
#            group = type, 
#            color = type, 
#            fill  = type)) +
#   geom_line() +
#   # geom_smooth(linewidth = 1, se=F) +
#   scale_x_discrete(c('lmYA', paste0('CB', 1:15))) +
#   # scale_y_continuous(limits = c(-1, 1)) +
#   scale_color_manual(values = c('orange', 'red', 'green')) +
#   scale_fill_manual(values = c('orange', 'red', 'green'))
# 
# 
# 
# df |> filter(EstMethod == 'lmYA' & type == 'negative') |> 
#   distinct(A, Y, EstMethod, ATE) 
# 
# 
# (grna_odm |> ondisc::get_feature_covariates())['PDL1g1', ]
# 
# ggplot(df |> 
#          filter(#EstMethod != 'lmYA' & 
#                 type == 'negative' &
#                 AY_idx <= 12 ),
#        aes(x = EstMethod, y = ATE)) +
#   geom_boxplot() +
#   coord_cartesian(ylim = c(-1, 1)) +
#   facet_wrap(vars(AY_idx), nrow = 3)
# 
# 
# df |> filter(EstMethod != 'lmYA' & type == 'negative') |> 
#   group_by(AY_idx) |> summarize(ATE = mean(ATE))
# 
# 
# # median w [q, 1-q] for every YA_idx 
# ggplot(df |> filter(#EstMethod != 'lmYA' & 
#                       type == 'negative') |> 
#          group_by(AY_idx, EstMethod) |> 
#          mutate(AY_idx = as.factor(AY_idx)) |>
#          summarize(ATE_median = median(ATE),
#                    ATE_IQR    = IQR(ATE, na.rm = TRUE),
#                    ATE_upper  = quantile(ATE, probs = 1 - .15),
#                    ATE_lower  = quantile(ATE, probs =     .15),
#                    .groups = 'drop'),
#        aes(x = EstMethod, 
#            # y = ATE, #|> abs() |> log(), 
#            color = AY_idx,
#            fill = AY_idx,
#            group = AY_idx)) +
#   geom_ribbon(aes(ymin = ATE_lower, ymax = ATE_upper), 
#               alpha = .3,
#               linewidth = 0) +
#   geom_line(aes(y = ATE_median), alpha = .7) +
#   theme(legend.position = 'none')
#   
# 
# # smoothing on estimates for every YA_idx
# ggplot(df |> filter(#EstMethod != 'lmYA' & 
#   type == 'negative') |> 
#     mutate(AY_idx = as.factor(AY_idx)),
#   aes(x = EstMethod, 
#       y = ATE, #|> abs() |> log(),
#       color = AY_idx,
#       fill = AY_idx,
#       group = AY_idx)) +
#   geom_smooth(alpha = .3,
#               linewidth = 0, na.rm = TRUE, span = .1) +
#   theme(legend.position = 'none')
# 
# 
# 
# # median w CI for every YA_idx
# ggplot(df |> filter(#EstMethod != 'lmYA' & 
#   type == 'negative') |> 
#     group_by(AY_idx, EstMethod) |> 
#     mutate(AY_idx = as.factor(AY_idx)) |>
#     summarize(ATE_median = median(ATE),
#               ATE_upper  = sort(ATE, na.last = TRUE)[ceiling(n()*.5 + 1.96*sqrt(n()*.5*(1-.5)))],
#               ATE_lower  = sort(ATE, na.last = TRUE)[  floor(n()*.5 - 1.96*sqrt(n()*.5*(1-.5)))],
#               .groups = 'drop'),
#   aes(x = EstMethod, 
#       # y = ATE, #|> abs() |> log(), 
#       color = AY_idx,
#       fill = AY_idx,
#       group = AY_idx)) +
#   # geom_ribbon(aes(ymin = ATE_lower, ymax = ATE_upper), 
#   #             alpha = .3,
#   #             linewidth = 0) +
#   geom_line(aes(y = ATE_median), alpha = .7) +
#   theme(legend.position = 'none')
# 
# 
# 
# AY_idx_larger = df |> filter(EstMethod == 'lmYA' & abs(ATE) >= .8) |>
#   # filter(AY_idx != 40) |> # remove AY_idx == 40 (outlier)
#   pull(AY_idx) |> unique()
# ggplot(df |> filter(type == 'negative' & AY_idx %in% AY_idx_larger &
#                       EstMethod %in% c('lmYA', paste0('CB', seq(1, 50, by = 5)))) |>
#          mutate(AY_idx = as.factor(AY_idx)),
#        aes(x = EstMethod,  y = ATE,
#            color = AY_idx,
#            fill = AY_idx)) +
#   geom_hline(aes(yintercept = 0), alpha = 1, linetype = 'solid', color = 'black', linewidth = .7) +
#   # geom_line(alpha = .5, linewidth = .3) +
#   geom_boxplot(outlier.shape = NA, alpha = .1, size = .1) +
#   scale_y_continuous(limits = c(-4, 4), breaks = seq(-20, 20, by = 1)) +
#   # facet_wrap(vars(type), ncol=1, scales = 'free_y') +
#   labs(title = 'Median ATE for AY pairs',
#        y = 'Median ATE', x = 'Estimation Method') +
#   theme(legend.position = 'none',
#         axis.text.x = element_text(angle = 45, vjust = 1.2, hjust = 1.1, size = 10))


