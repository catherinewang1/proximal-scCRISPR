


args = commandArgs(trailingOnly = TRUE)
args = c('laptop', 'allPerturbations')


require(assertthat) # for some assert statements

library(dplyr)
library(ggplot2)
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

dir.create(sprintf('%s/spca/cbgenes/%s/pvalues/', save_dir, AYZW_setting_name), recursive = TRUE, showWarnings = FALSE)
cbgenes_setting_name = 'spca/cbgenes'
CB_setting_name = 'simple'

# load chosen AYZW names
AY   = read.csv(sprintf('%s/spca/cbgenes/%s/AY.csv', save_dir, AYZW_setting_name))

# Load lin regressions
ATE_lm = read.csv(sprintf('%s/spca/cbgenes/%s/lmYAU/ATE_lmYAU.csv', save_dir, AYZW_setting_name))

# Load Poisson and Negative binomial regressions
# ATE_glm = read.csv(sprintf('%s/spca/cbgenes/%s/glmYAU/ATE_glmYAU.csv', save_dir, AYZW_setting_name))
ATE_glm = read.csv(sprintf('%s/spca/cbgenes/%s/glmYAU/ATE_glmYAU_loglibsize.csv', save_dir, AYZW_setting_name))

# Load CB methods (pci2s and MEst using max numNC)
ATE_CB = read.csv(sprintf('%s/%s/%s/%s/ATE.csv', save_dir, cbgenes_setting_name, AYZW_setting_name, CB_setting_name))
ATE_CB = ATE_CB |> filter(EstMethod == 'OCB2SLSpci2s' | EstMethod == 'OCBMEstRWbasis1') |>
                   filter(numNC == 20) |> select(-X)
ATE_CB[ATE_CB$EstMethod == 'OCBMEstRWBasis1', 'method_type'] = 'mest'

dim(ATE_CB); head(ATE_CB)
ATE_CB |> pull(method_type)  |> table()

head(ATE_lm)
head(ATE_glm)


ATE = rbind(ATE_lm, ATE_glm, ATE_CB |> select(-EstMethod, -EstMethod_l))

# ATE = read.csv(sprintf('%s/spca/cbgenes/%s/glmYAU/ATE_glmYAU.csv', save_dir, AYZW_setting_name))
df = merge(AY |> mutate(AY_idx = 1:nrow(AY)), 
           ATE, 
           by = 'AY_idx')
head(df); dim(df)
df$method |> table()


#' @param myType (string) eg positive, negative, maybe
#' @param myMethod (string)
get_pval_plotdf <- function(myType, myMethod) {
  # myType = 'negative'; myMethod = 'lmYA'; myMethod = 'OCB2SLSpci2s'
  pvals = df |> filter(type == myType & method == myMethod) |> pull(pval) |> sort()
  data.frame(EstMethod = myMethod,
             pvals = pvals,
             theoretical = ppoints(length(pvals)))
} 
get_pval_plotdf('negative', 'poisYAU') |> pull(pvals) |> hist()
get_pval_plotdf('negative', 'OCB2SLSpci2s') |> pull(pvals) |> hist()
get_pval_plotdf('negative', 'OCBMEstRW') |> pull(pvals) |> hist()


method_list = c('lmYA', 'lmYAU', 
                'poisYA', 'poisYAU',
                'nbYA', 'nbYAU', 
                'OCB2SLSpci2s', 'OCBMEstRW' )
method_type_list = c('linearmodel', 'linearmodel',
                     'poisson', 'poisson',
                     'negativebinomial', 'negativebinomial',
                     'pci2s', 'mest')
method_adjusted = c('unadjusted', 'adjusted',
                    'unadjusted', 'adjusted',
                    'unadjusted', 'adjusted',
                    'adjusted', 'adjusted')


#' @param myType (string) eg positive, negative, maybe
make_pval_plot <- function(myType) {
  
}


# collect df
myType = 'negative'
plot_df = NULL
for(j in 1:length(method_list)) {
  plot_df = rbind(plot_df, 
                  cbind(get_pval_plotdf(myType, method_list[j]), 
                        method_type = method_type_list[j],
                        adjusted    = method_adjusted[j]))
}
# dim(plot_df)



plot_df = plot_df  |>
  mutate(EstMethod = factor(EstMethod, levels = method_list)) |>
  mutate(pvals_trans = -log(pvals, base = 10),
         theoretical_trans = -log(theoretical, base = 10))




# ====== save tests' p-value results


#' get the p-value for myType and myMethod. keep AY_idx info
#' @param myType (string) eg positive, negative, maybe
#' @param myMethod (string)
get_pval_plotdf_AYidx <- function(myType, myMethod) {
  # myType = 'negative'; myMethod = 'lmYA'; myMethod = 'OCB2SLSpci2s'
  pvals = df |> filter(type == myType & method == myMethod) |> select(AY_idx, pval) |> arrange(pval)
  data.frame(AY_idx = pvals$AY_idx,
             type = myType, 
             EstMethod = myMethod,
             pval = pvals$pval,
             theoretical = ppoints(nrow(pvals)))
} 


pvalplot_df = NULL
for(myType in c('negative', 'positive')) {
  for(j in 1:length(method_list)) {
    pvalplot_df = rbind(pvalplot_df, 
                            cbind(get_pval_plotdf_AYidx(myType, method_list[j]), 
                                  method_type = method_type_list[j],
                                  adjusted    = method_adjusted[j]))
  }
}

pvalplot_df = pvalplot_df  |>
  mutate(EstMethod = factor(EstMethod, levels = method_list)) |>
  mutate(pvals_trans = -log(pval, base = 10),
         theoretical_trans = -log(theoretical, base = 10))
write.csv(x = pvalplot_df, 
          file = sprintf('%s/spca/cbgenes/%s/pvalues/pvals_papalexi.csv', save_dir, AYZW_setting_name),
          row.names = FALSE)

# pvalplot_df |> group_by(type, EstMethod) |> summarize(count = n()) # sanity check


# # ====== positive (not used)
# positiveplot_df = NULL
# for(j in 1:length(method_list)) {
#   positiveplot_df = rbind(positiveplot_df, 
#                   cbind(get_pval_plotdf_AYidx('positive', method_list[j]), 
#                         method_type = method_type_list[j],
#                         adjusted    = method_adjusted[j]))
# }
# positiveplot_df = positiveplot_df  |>
#   mutate(EstMethod = factor(EstMethod, levels = method_list)) |>
#   mutate(pvals_trans = -log(pval, base = 10),
#          theoretical_trans = -log(theoretical, base = 10))
# 
# write.csv(x = positiveplot_df, 
#           file = sprintf('%s/spca/cbgenes/%s/pvalues/pvals_%s.csv', save_dir, AYZW_setting_name, 'positive'),
#           row.names = FALSE)
# rm(positiveplot_df)
# ======



# save plot_df for 2 Step Hyp Testing 


EstMethodColors = c("#66C2A5", "#378970",
                    "#FC8D62", "#c34f22",
                    "#8DA0CB", "#51689c",
                    "#e93ba6", "#60ad25")

# qqplot of pvals vs unif
# as points
# p1 = ggplot(plot_df, aes(sample=pvals, group = EstMethod, color = EstMethod, linetype = adjusted)) +
#   geom_abline(slope = 1, intercept = 0, color = 'purple3', linewidth = 1, alpha = .8) +
#   # stat_qq_line(distribution = stats::qunif, color = 'purple', line.p = c(0, 1)) +
#   stat_qq(distribution = stats::qunif, size = 1, alpha = .8) +
#   # scale_color_brewer(palette = "Set2") +
#   scale_color_discrete(type = EstMethodColors) +
#   labs(title = 'QQ Plot of p-values against Unif(0,1)',
#        x = 'theoretical Unif(0,1)', y = 'sample p-values') +
#   theme(legend.position = 'inside', 
#         legend.position.inside = c(.1, .75))
# as line 
p1 = ggplot(plot_df, 
            aes(x = theoretical, y = pvals, group = EstMethod, color = EstMethod, linetype = adjusted)) +
  geom_abline(slope = 1, intercept = 0, 
              # color = 'purple3',
              color = 'black',
              linewidth = 1, alpha = .8) +
  # geom_point(size = 1, alpha = .8) +
  geom_line( linewidth = 1, alpha = .8) +
  # scale_color_brewer(palette = "Set2") +
  scale_color_discrete(type = EstMethodColors) +
  labs(title = 'QQ Plot of p-values against Unif(0,1)',
       x = 'theoretical Unif(0,1)', y = 'sample p-values') +
  theme(# legend.position = 'inside', 
         # legend.position.inside = c(.1, .75)
         legend.position = 'right')

ggsave(plot = p1,
       filename = sprintf('%s/spca/cbgenes/%s/pvalues/%s_pval_qqplot.pdf', save_dir, AYZW_setting_name, myType),
       height = 5, width = 6.5)

# qqplot of transformed pvals vs unif
p2 = ggplot(plot_df, 
            aes(# x = theoretical, y = pvals,
              x = theoretical_trans, y = pvals_trans,
              group = EstMethod, color = EstMethod)) +
  
  # stat_qq_line(distribution = stats::qunif, color = 'purple', line.p = c(0, 1)) +
  # stat_qq(distribution = stats::qunif, size = 1, alpha = .8) +
  geom_abline(slope = 1, intercept = 0, 
              # color = 'purple3',
              color = 'black', 
              linewidth = 1, alpha = .8) +
  geom_point(size = 1) +
  # scale_color_brewer(palette = "Set2") +
  scale_color_discrete(type = EstMethodColors) +
  labs(title = '-log(p-values) vs -log(Unif(0,1))',
       x = 'theoretical -log(Unif(0,1))', y = 'sample -log(pval)') +
  theme(# legend.position = 'inside', 
        #legend.position.inside = c(.1, .75)
        legend.position = 'right')
p2
ggsave(plot = p2,
       filename = sprintf('%s/spca/cbgenes/%s/pvalues/%s_pval_qqplot_trans.pdf', save_dir, AYZW_setting_name, myType),
       height = 5, width = 6.5)

ggsave(plot = p2 + 
         scale_y_continuous(limits = c(0, 25)) +
         scale_x_continuous(limits = c(0, 2.5)),
       filename = sprintf('%s/spca/cbgenes/%s/pvalues/%s_pval_qqplot_trans_cutoff.pdf', save_dir, AYZW_setting_name, myType),
       height = 5, width = 6.5)

# Sample plots but exclude Proximal Methods (eg only naive + covariates)

p1ii = ggplot(plot_df |> filter(method_type != 'pci2s' & method_type != 'mest'), 
            aes(x = theoretical, y = pvals, group = EstMethod, color = EstMethod, linetype = adjusted)) +
  geom_abline(slope = 1, intercept = 0, 
              # color = 'purple3',
              color = 'black', 
              linewidth = 1, alpha = .8) +
  # geom_point(size = 1, alpha = .8) +
  geom_line( linewidth = 1, alpha = .8) +
  # scale_color_brewer(palette = "Set2") +
  scale_color_discrete(type = EstMethodColors) +
  labs(title = 'QQ Plot of p-values against Unif(0,1)',
       x = 'theoretical Unif(0,1)', y = 'sample p-values') +
  theme(# legend.position = 'inside', 
    # legend.position.inside = c(.1, .75)
    legend.position = 'right')

ggsave(plot = p1ii,
       filename = sprintf('%s/spca/cbgenes/%s/pvalues/%s_pval_qqplot_noproximal.pdf', save_dir, AYZW_setting_name, myType),
       height = 5, width = 6)

# qqplot of transformed pvals vs unif
p2ii = ggplot(plot_df |> filter(method_type != 'pci2s' & method_type != 'mest'), 
            aes(# x = theoretical, y = pvals,
              x = theoretical_trans, y = pvals_trans,
              group = EstMethod, color = EstMethod)) +
  geom_point(size = .9) +
  # stat_qq_line(distribution = stats::qunif, color = 'purple', line.p = c(0, 1)) +
  # stat_qq(distribution = stats::qunif, size = 1, alpha = .8) +
  geom_abline(slope = 1, intercept = 0, 
              # color = 'purple3',
              color = 'black',  linewidth = 1, alpha = .8) +
  # scale_color_brewer(palette = "Set2") +
  scale_color_discrete(type = EstMethodColors) +
  labs(title = '-log(p-values) vs -log(Unif(0,1))',
       x = 'theoretical -log(Unif(0,1))', y = 'sample -log(pval)') +
  theme(# legend.position = 'inside', 
    #legend.position.inside = c(.1, .75)
    legend.position = 'right')
p2ii
ggsave(plot = p2ii,
       filename = sprintf('%s/spca/cbgenes/%s/pvalues/%s_pval_qqplot_noproximal_trans.pdf', save_dir, AYZW_setting_name, myType),
       height = 5, width = 6)




ggsave(plot = p2ii + 
              scale_y_continuous(limits = c(0, 25)) +
              scale_x_continuous(limits = c(0, 2.5)),
       filename = sprintf('%s/spca/cbgenes/%s/pvalues/%s_pval_qqplot_noproximal_trans_cutoff.pdf', save_dir, AYZW_setting_name, myType),
       height = 5, width = 6)






