# ------------------------------------------------------------------------------------------------ #
#     Estimate ATE/Log FC
# For only 1/few choices of Proximal....
# - Positive Tests: perform all
# - Discovery Tests: 
#    + Perform initial Linear Reg for ALL/A LOT OF tests
#    + For those with significance (small p-values), perform proximal
# ------------------------------------------------------------------------------------------------ #
args = commandArgs(trailingOnly = TRUE)
# args = c('macbook')


assertthat::assert_that(length(args) > 0, msg="must give arg for specifying device eg 'Rscript <filename>.R ubergenno'")
DEVICE = args[1]
source('../PATHS.R') # load in data_dir and save_dir and CODE_DIR, depending on DEVICE value
assertthat::assert_that(!is.null(data_dir), msg='first arg must be: laptop, desktop, or ubergenno')

library(dplyr)
library(future.apply)
# options(future.globals.maxSize= 850*1024^2) #1st num is MB
options(future.globals.maxSize= 2500*1024^2) #1st num is MB
# plan(multisession, workers = 24)
plan(multisession, workers = 15)
# plan(sequential)

subfoldername = 'excitatoryneurons_alldiscovery'




dir.create(sprintf('%s/AY/', save_dir))
dir.create(sprintf('%s/AY/%s/', save_dir, subfoldername))
dir.create(sprintf('%s/AY/%s/initiallinreg/', save_dir, subfoldername))

# ============================================================================== #
# =================== LOAD ==========================================
# ============================================================================== #


# process to #features x #cells to be consistent
# 
# grna: perturbation-by-cell 0/1 trt assignment
# gene:      gene-by-cell count gene expression 
# gene_norm: gene-by-cell normalized gene expression
# 
# but others may be #cells x #features
# metadata: #cells x 13 variables




# =================== \__Load metadata: carried over from Broad ===========================
# metadata = read.csv(sprintf('%s/processGEO/grna.csv', data_dir)) # for all
metadata0 = read.csv(sprintf('%s/cell_metadata.csv', save_dir), row.name = 'X') # excitatory neuron cells w/ a pert assignment w/ > 80(?) sample size
gene_dev = read.csv(sprintf('%s/gene_dev.csv', save_dir))
# Add TF info, but does not have for GFP (will be NA cols for is_transcription_factor and importance_rank)
# merge(metadata, gene_dev |> select(gene_name, is_transcription_factor, importance_rank), 
#       by.x = 'Perturbation_alias', by.y = 'gene_name') |> dim()
# all(metadata$Perturbation %in% gene_dev$gene_name)
# all(metadata$Perturbation_alias %in% gene_dev$gene_name)
# unique(metadata$Perturbation_alias)[which(!unique(metadata$Perturbation_alias) %in% gene_dev$gene_name)]
# metadata |> filter(Perturbation_alias %in% c('Mll1', 'Myst4'))
# metadata |> filter(Perturbation_alias %in% c('Kmt2a'))
# 'Kmt2a' %in% gene_dev$gene_name
metadata = merge(metadata0, gene_dev |> select(gene_name, is_transcription_factor, importance_rank), 
                 by.x = 'Perturbation_alias', by.y = 'gene_name', all.x = TRUE)
row.names(metadata) = metadata$NAME
metadata = metadata[row.names(metadata0), ] # reorder to same as original


# dim(metadata); dim(metadata0)
# row.names(metadata)[1:5];row.names(metadata0)[1:5]
# all(row.names(metadata) %in% row.names(metadata0))
# all(row.names(metadata0) %in% row.names(metadata))
# all(row.names(metadata0) == row.names(metadata))
rm(metadata0, gene_dev)

# =================== \__Load grna: carried over from Broad ===========================
grna = read.csv(sprintf('%s/grna.csv', save_dir), row.name = 'X') # for cells w a pert assignment 22085 x 39
grna = t(grna) #  30 perturbations x 4744 cells
row.names(grna) = stringr::str_to_title(row.names(grna)) # hopefully these match broad Perturbation names now

# =================== \__Load counts: processed from GEO ===========================
# load in top 4000 genes (count and continuous)
gene = readRDS(sprintf('%s/gene_counts.rds', save_dir))
gene = as.matrix(gene)

# =================== \__Load continuous: normalized counts ==============================================
gene_norm = readRDS(sprintf('%s/gene_norm.rds', save_dir))
gene_norm = as.matrix(gene_norm)

(dim(gene_norm) == dim(gene)) |> all() # same dimensions (all genes)
(row.names(gene_norm) == row.names(gene)) |> all() # same genes and order
(row.names(gene_norm) == row.names(gene)) |> table() # same genes and order
(colnames(gene_norm) == colnames(gene)) |> all()  # same cells and order

# =================== \__Reorder cells in data to match gene/gene_norm order ===========================
grna      = grna[, colnames(gene)]
metadata  = metadata[colnames(gene), ]


(colnames(gene_norm) == colnames(grna)) |> table() # same cells
(colnames(gene_norm) == row.names(metadata)) |> table() # same cells

format(object.size(gene), units = 'Gb') # less when loaded as original sparse matrix form. but then transformed to matrix
format(object.size(gene_norm), units = 'Gb')
format(object.size(grna), units = 'Gb')
format(object.size(metadata), units = 'Gb')



# =================== \__Load gene info (dev, TF, ...) ==============================================

gene_dev_df = read.csv(file = sprintf('%s/gene_dev.csv', save_dir)) # load gene dev info
gene_dev_df = gene_dev_df |> filter(gene_name %in% row.names(gene)) # only choose top 4000 and grna target
gene_dev_df$gene_name_upper = toupper(gene_dev_df$gene_name)
# gene_dev_df |> filter(gene_name_upper == 'PISD') # there is Pisd and PISD apparently??
row.names(gene_dev_df) = gene_dev_df$gene_name
gene_dev_df = gene_dev_df[row.names(gene), ] # select the rows=genes that are top or grna target (the selected gene and gene_norm)


(row.names(gene_norm) %in% gene_dev_df$gene_name) |> all() 
(row.names(gene) %in% gene_dev_df$gene_name) |> all() 

(row.names(gene_norm) == gene_dev_df$gene_name) |> table()
(row.names(gene) == gene_dev_df$gene_name) |> table()

gene_dev_df |> head()
# gene_name gene_idx deviance importance_rank is_grna_target genename_cleaned is_transcription_factor gene_name_upper
# Actb        Actb     9932 579441.3               1          FALSE             actb                   FALSE            ACTB
# Nrgn        Nrgn    17935 242049.9               2          FALSE             nrgn                   FALSE            NRGN
# Malat1    Malat1    27229 232525.0               3          FALSE           malat1                   FALSE          MALAT1
# Tuba1a    Tuba1a    23642 187071.7               4          FALSE           tuba1a                   FALSE          TUBA1A
# Pcp4        Pcp4    24555 179699.0               5          FALSE             pcp4                   FALSE            PCP4
# Calm2      Calm2    25756 136040.4               6          FALSE            calm2                   FALSE           CALM2



# ============================================================================== #
# =================== CHOOSE AY TESTS ==========================================
# ============================================================================== #
# filter for perturbations (in this example dataset, they seem ok already, but we can filter here too)
rowSums(grna) |> summary()
rowSums(grna) |> hist(breaks = 50)
rowSums(grna) |> hist(breaks = 50, xlim = c(0, 500))
GRNA_SAMPLESIZE_MIN = 50 # already filtered out in prev step at 50, but can set higher to filter even more

all_As = row.names(grna)[rowSums(grna) >= GRNA_SAMPLESIZE_MIN]
all_As = setdiff(all_As, 'Gfp')

# filter for genes: larger #cells with non-zero expr and higher rank
rowSums(gene == 0) |> hist(breaks = 50)
GENE_NUMNONZEROCELLS_MIN = 500 # num of cells w nonzero counts
# GENE_IMPORTANCERANK_MAX = 1500 # maximum importance rank of gene, ranked by gene deviance
GENE_IMPORTANCERANK_MAX = 1500 # maximum importance rank of gene, ranked by gene deviance, only for singlegene NC
GENE_NUMRESPONSES = 1000 # number of disc to pull 

# all_Ys = row.names(gene)[(rowSums(gene == 0) >= GENE_NUMNONZEROCELLS_MIN) & 
#                            (gene_dev_df$importance_rank <= GENE_IMPORTANCERANK_MAX)] # make sure order matches

genes_nonzero = row.names(gene)[(rowSums(gene == 0) >= GENE_NUMNONZEROCELLS_MIN)] # genes that meet the nonzero threshold
all_Ys_df = gene_dev_df |>  
            filter(gene_name %in% genes_nonzero) |>  
            filter(!is_grna_target & !is_transcription_factor) |>
            arrange(importance_rank)

validsinglegenes = all_Ys_df |> filter(importance_rank <= GENE_IMPORTANCERANK_MAX) |> pull(gene_name) # valid names for single genes
all_Ys = all_Ys_df |> head(min(GENE_NUMRESPONSES), GENE_NUMRESPONSES) |> pull(gene_name)


# --------  \__'positive' tests  ------ 
# also require 'positive' tests, even if their targeted gene did not pass the qc (some grna targeted genes are not here in the subset...)
perturbation_alias_df = read.csv(sprintf('%s/perturbation_alias.csv', save_dir))
row.names(perturbation_alias_df) = tolower(perturbation_alias_df$Perturbation)

pos_tests_grna = c() # the grna/Perturbation name
pos_tests_alias = c() # the gene alias name (name in gene and gene_norm) 
for(pert_name in all_As) {
  pert_name_alias_lower = perturbation_alias_df[tolower(pert_name), 'Perturbation_alias_lower']
  if(pert_name_alias_lower %in% tolower(row.names(gene))) {
    pos_tests_grna  = c(pos_tests_grna, pert_name)
    pos_tests_alias = c(pos_tests_alias, stringr::str_to_title(pert_name_alias_lower))
  }
}


# pos_tests_vec = row.names(grna)[row.names(grna) %in% row.names(gene)]
pos_tests_df = data.frame(grna = pos_tests_grna, gene = pos_tests_alias, type = 'positive')




# --------  \__form all discovery tests  ------ 
# all discovery tests:
all_disc_df = expand.grid(grna=all_As, gene=all_Ys, stringsAsFactors = FALSE) 
all_disc_df = merge(all_disc_df, perturbation_alias_df |> select(Perturbation, Perturbation_alias),
                     by.x = 'grna', by.y = 'Perturbation')
all_disc_df = all_disc_df |> dplyr::filter(gene != Perturbation_alias)

all(all_disc_df$grna %in% row.names(grna))
all(all_disc_df$gene %in% row.names(gene))
all(all_disc_df$gene_norm %in% row.names(gene_norm))


write.csv(all_disc_df, sprintf('%s/AY/%s/initiallinreg/all_discovery_AY.csv', save_dir, subfoldername))


# --------  \__run initial lin reg fits  ------ 
gfp_idx = which(grna['Gfp', ]==1)
# convert 0/1 matrix into a column
grna_A = row.names(grna)[apply(X = (grna==1), FUN = which, MARGIN = 2)] |> factor()
grna_A = relevel(grna_A, 'Gfp')

libsize_log = log((metadata$nUMI |> as.numeric())+ 1e-20)
# all(row.names(metadata) == colnames(gene))


#' #' Requires in environment:
#' #'   gfp_idx, gene_norm, grna
#' #' @param AY_idx (integer) idx of the test in the df of all discovery pairs
#' fit_lm <- function(AY_idx, AY_df) {
#'   # AY_idx = 2
#' 
#'   cur_grna = AY_df[AY_idx, 'grna']
#'   cur_gene = AY_df[AY_idx, 'gene']
#' 
#'   cur_idx = c(which(grna[cur_grna, ]==1), gfp_idx)
#'   linreg_df = data.frame(Y = gene_norm[cur_gene, cur_idx],
#'                          A =      grna[cur_grna, cur_idx])
#' 
#'   linreg_fit = lm('Y ~ A', data = linreg_df)
#'   linreg_fit_summ = summary(linreg_fit)
#'   linreg_res = (data.frame(linreg_fit_summ$coefficients))['A', ]
#'   colnames(linreg_res) = c('estimate', 'se', 'tstat', 'pval')
#'   row.names(linreg_res) = NULL
#' 
#'   cbind(data.frame(grna = cur_grna, gene = cur_gene), linreg_res)
#' }



#' Requires in environment:
#'   gfp_idx, gene_norm, grna
#' @param AY_idx (integer) idx of the test in the df of all discovery pairs
fit_glms <- function(AY_idx, AY_df) {
  # AY_idx = 2
  
  cur_grna = AY_df[AY_idx, 'grna']
  cur_gene = AY_df[AY_idx, 'gene']
  
  cur_idx = c(which(grna[cur_grna, ]==1), gfp_idx)
  df = data.frame(Y_count =      gene[cur_gene, cur_idx],
                         Y_norm  = gene_norm[cur_gene, cur_idx],
                         A       =      grna[cur_grna, cur_idx], 
                         libsize_log = libsize_log[cur_idx])
  
  # Lin Reg
  linreg_fit = lm('Y_norm ~ A', data = df)
  linreg_fit_summ = summary(linreg_fit)
  linreg_res = (data.frame(linreg_fit_summ$coefficients))[-1, ]
  colnames(linreg_res) = c('estimate', 'se', 'tstat', 'pval')
  row.names(linreg_res) = NULL
  
  # Neg Bin
  nbfit = MASS::glm.nb(formula='Y_count ~ A + offset(libsize_log)', data = df)
  nbfit_summary = summary(nbfit)
  nb_res = (data.frame(nbfit_summary$coefficients))[-1, ]
  colnames(nb_res) = c('estimate', 'se', 'tstat', 'pval')
  row.names(nb_res) = NULL
  
  
  return(rbind(linreg_res |> mutate(method =     'lm', grna=cur_grna, gene=cur_gene, .before='estimate'), 
               nb_res     |> mutate(method = 'negbin', grna=cur_grna, gene=cur_gene, .before='estimate')))
}




#' #' Fit 1 Lin Reg per outcome gene (vs each gene-grna perturbation)
#' #' Requires in environment:
#' #'   grna_A, gene_norm,
#' #' @param cur_gene (character) response gene Y name (needs to be a rowname of both gene and gene_norm)
#' fit_lm_one <- function(cur_gene) {
#'   
#'   linreg_df = data.frame(Y = gene_norm[cur_gene, ],
#'                          A_ =      grna_A)
#'   
#'   linreg_fit = lm('Y ~ A_', data = linreg_df)
#'   linreg_fit_summ = summary(linreg_fit)
#'   linreg_res = (data.frame(linreg_fit_summ$coefficients))[-1, ]
#'   colnames(linreg_res) = c('estimate', 'se', 'tstat', 'pval')
#'   linreg_res = linreg_res |> dplyr::mutate(grna=sapply(row.names(linreg_res), FUN=function(x){substr(x, 3, nchar(x))}), 
#'                                            gene=cur_gene, .before=1)
#'   row.names(linreg_res) = NULL
#'   
#'   return(linreg_res)
#' }






#' Fit 1 Neg Bin GLM AND Lin Reg per outcome gene (vs each gene-grna perturbation)
#' Requires in environment:
#'   grna_A, gene, gene_norm, libsize_log
#' @param cur_gene (character) response gene Y name (needs to be a rowname of both gene and gene_norm)
fit_glms_one <- function(cur_gene) {
  
  df = data.frame(Y_count =      gene[cur_gene, ],
                  Y_norm  = gene_norm[cur_gene, ],
                  A_      = grna_A, 
                  libsize_log = libsize_log)
  
  # Lin Reg
  linreg_fit = lm('Y_norm ~ A_', data = df)
  linreg_fit_summ = summary(linreg_fit)
  linreg_res = (data.frame(linreg_fit_summ$coefficients))[-1, ]
  colnames(linreg_res) = c('estimate', 'se', 'tstat', 'pval')
  linreg_res = linreg_res |> dplyr::mutate(grna=sapply(row.names(linreg_res), FUN=function(x){substr(x, 3, nchar(x))}), 
                                           gene=cur_gene, .before=1)
  row.names(linreg_res) = NULL
  
  # Neg Bin
  nbfit = MASS::glm.nb(formula='Y_count ~ A_ + offset(libsize_log)', data = df)
  nbfit_summary = summary(nbfit)
  nb_res = (data.frame(nbfit_summary$coefficients))[-1, ]
  colnames(nb_res) = c('estimate', 'se', 'tstat', 'pval')
  nb_res = nb_res |> dplyr::mutate(grna=sapply(row.names(nb_res), FUN=function(x){substr(x, 3, nchar(x))}), 
                                   gene=cur_gene, .before=1)
  row.names(nb_res) = NULL
  
  
  return(rbind(linreg_res |> mutate(method = 'lm', .before='estimate'), 
               nb_res |> mutate(method = 'negbin', .before='estimate')))
}


# ================ === === === === === 


# more relaxed filtering and also allows TFs and grna targets
many_Ys = gene_dev_df |>  
  filter(gene_name %in% genes_nonzero &
           (is.na(importance_rank) | importance_rank<=2500)) |> pull(gene_name)
many_Ys = union(many_Ys, perturbation_alias_df$Perturbation_alias)
many_Ys = many_Ys[many_Ys %in% row.names(gene)]



if('all_glm_one.csv' %in% list.files(sprintf('%s/AY/%s/initiallinreg/', save_dir, subfoldername))) {
  print(sprintf("[%s]        - Loading glm (one fit per gene) on all pairs", Sys.time()))
  all_glm_one_df = read.csv(sprintf('%s/AY/%s/initiallinreg/all_glm_one.csv', save_dir, subfoldername))
} else {
  print(sprintf("[%s]        - Fitting glm (one fit per gene) on all pairs", Sys.time()))
  t0 = Sys.time()
  all_glm_one = future.apply::future_mapply(fit_glms_one,
                                        cur_gene = many_Ys, 
                                        future.globals = TRUE,
                                        future.seed = 123456, 
                                        SIMPLIFY = FALSE)
  t1 = Sys.time()
  print(sprintf("[%s]        - %2.2f mins", Sys.time(), difftime(t1, t0, units = "mins")))
  
  all_glm_one_df = NULL
  for(idx in 1:length(all_glm_one)) {   all_glm_one_df = rbind(all_glm_one_df, all_glm_one[[idx]])   }
  write.csv(all_glm_one_df, sprintf('%s/AY/%s/initiallinreg/all_glm_one.csv', save_dir, subfoldername))
  rm(all_glm_one)
}



if('all_glm.csv' %in% list.files(sprintf('%s/AY/%s/initiallinreg/', save_dir, subfoldername))) {
  print(sprintf("[%s]        - Loading glm (one fit per grna-gene) on all pairs", Sys.time()))
  all_glm_df = read.csv(sprintf('%s/AY/%s/initiallinreg/all_glm.csv', save_dir, subfoldername))
} else {
  print(sprintf("[%s]        - Fitting glm (one fit per grna-gene) on all pairs", Sys.time()))
  t0 = Sys.time()
  all_glm = future.apply::future_mapply(fit_glms,
                                        AY_idx = 1:nrow(all_glm_one_df),
                                        MoreArgs = list(AY_df = all_glm_one_df |> select(grna, gene)),
                                        future.globals = TRUE,
                                        future.seed = 123456, 
                                        SIMPLIFY = FALSE)
  t1 = Sys.time()
  print(sprintf("[%s]        - %2.2f mins", Sys.time(), difftime(t1, t0, units = "mins")))
  
  all_glm_df = NULL
  for(idx in 1:length(all_glm)) {   all_glm_df = rbind(all_glm_df, all_glm[[idx]])   }
  write.csv(all_glm_df, sprintf('%s/AY/%s/initiallinreg/all_glm.csv', save_dir, subfoldername))
  rm(all_glm)
}


disc_lm_df = all_glm_df


if(F) {
  

all_glm_compare = merge(all_glm_one_df |> 
               select(method, grna, gene, estimate, pval, se) |>
               rename(estimate_one=estimate, pval_one=pval, se_one=se), 
             all_glm_df |> 
               select(method, grna, gene, estimate, pval, se) |>
               rename(estimate_pairwise=estimate, pval_pairwise=pval, se_pairwise=se), 
             by = c('method', 'grna', 'gene'))

all_glm_compare = merge(perturbation_alias_df |> rename(grna=Perturbation, grna_target=Perturbation_alias) |> select(grna, grna_target),
                       all_glm_compare, 
                       all.y = TRUE)
all_glm_compare = all_glm_compare |> mutate(type_bool = (grna_target == gene), 
                                          type = sapply(type_bool, function(x){if(x){'positive'} else{'discovery'}}))

# add treatment sample size
grna_samplesize = data.frame(grna_samplesize = rowSums(grna)) |> tibble::rownames_to_column('grna')
all_glm_compare = merge(all_glm_compare, grna_samplesize, by = 'grna')




hist(all_glm_compare$estimate_one - all_glm_compare$estimate_pairwise)
hist(all_glm_compare$pval_one - all_glm_compare$pval_pairwise)

# those with large differences seem to be ones with large estimates and large se's (so have large p-values anyway) 
all_glm_compare |> mutate( estimate_diff = estimate_one - estimate_pairwise) |> filter(abs(estimate_diff) > 5) |> arrange(desc(abs(estimate_diff)))
all_glm_compare |> mutate( estimate_diff = estimate_one - estimate_pairwise) |> arrange(desc(abs(estimate_diff)))



# compare one glm vs pairwise glms: estimates
ggplot(all_glm_compare, aes(x = estimate_one, y = estimate_pairwise, color = grna_samplesize)) +
  geom_abline(aes(intercept = 0, slope = 1)) +
  geom_point() + 
  facet_grid(rows = vars(type), cols = vars(method), scales = 'free')


# compare one glm vs pairwise glms: pvals
ggplot(all_glm_compare, aes(x = pval_one, y = pval_pairwise)) +
  geom_abline(aes(intercept = 0, slope = 1)) +
  geom_point() + 
  facet_grid(rows = vars(type), cols = vars(method), scales = 'free')


ggplot(all_glm_compare, aes(x = estimate_one)) +
  geom_histogram() + 
  facet_grid(rows = vars(type), cols = vars(method), scales = 'free')



ggplot(all_glm_compare, aes(x = pval_one)) +
  geom_histogram() + 
  facet_grid(rows = vars(type), cols = vars(method), scales = 'free')




ggplot(all_glm_compare, aes(x = estimate_pairwise)) +
  geom_histogram() + 
  facet_grid(rows = vars(type), cols = vars(method), scales = 'free')



ggplot(all_glm_compare, aes(x = pval_pairwise)) +
  geom_histogram() + 
  facet_grid(rows = vars(type), cols = vars(method), scales = 'free')





row.names(grna)[row.names(grna) %in% many_Ys]
row.names(grna)[!row.names(grna) %in% many_Ys]


all_glm_one_df |> filter(grna == gene)

all_lm$Satb2
all(colnames(grna) == colnames(gene))

perturbation_alias_df |> filter(Perturbation %in% row.names(grna)) 
perturbation_alias_df |> filter(!Perturbation_alias %in% row.names(grna)) 


all_glm_one_df = merge(perturbation_alias_df |> rename(grna=Perturbation, grna_target=Perturbation_alias) |> select(grna, grna_target),
                       all_glm_one_df, 
                  all.y = TRUE)
all_glm_one_df = all_glm_one_df |> mutate(type_bool = (grna_target == gene), 
                    type = sapply(type_bool, function(x){if(x){'positive'} else{'discovery'}}))


all_glm_one_df |>     filter(type == 'positive') |> pull(pval) |> hist(main = 'lm (1 fit per gene) pval of positive tests', breaks = seq(0, 1, length.out = 100))


all_glm_one_df |> mutate(type_bool = (grna_target == gene), 
                    type = sapply(type_bool, function(x){if(x){'positive'} else{'discovery'}})) |> 
  filter(type == 'positive') |> pull(pval) |> hist(main = 'lm (1 fit per gene) pval of positive tests', breaks = seq(0, 1, length.out = 100))

ggplot(all_glm_one_df, aes(x = pval)) + 
  geom_histogram() +
  facet_grid(rows = vars(type), cols = vars(method), scales = 'free')






dim(grna)

# rds matches h5 save
library(h5lite)
h5file <- sprintf('%s/gene.h5', save_dir)
gene_counts_h5 <- h5_read(h5file, "counts/counts") # R reads in row and colnames

colnames(grna)[1:2]
all(colnames(grna) == colnames(gene_counts_h5))
all(gene_counts_h5 == gene)



}








# ================ === === === === === 


# 
# if('all_discovery_linereg.csv' %in% list.files(sprintf('%s/AY/%s/initiallinreg/', save_dir, subfoldername))) {
#   print(sprintf("[%s]        - Loading lm on all discovery pairs", Sys.time()))
#   disc_lm_df = read.csv(sprintf('%s/AY/%s/initiallinreg/all_discovery_linereg.csv', save_dir, subfoldername))
# } else {
#   print(sprintf("[%s]        - Fitting lm on all discovery pairs", Sys.time()))
#   t0 = Sys.time()
#   disc_lm = future.apply::future_mapply(fit_lm,
#                                         AY_idx = 1:nrow(all_disc_df), 
#                                         MoreArgs = list(AY_df = all_disc_df),
#                                         future.globals = TRUE,
#                                         future.seed = 123456, 
#                                         SIMPLIFY = FALSE)
#   t1 = Sys.time()
#   print(sprintf("[%s]        - %2.2f mins", Sys.time(), difftime(t1, t0, units = "mins")))
#   
#   disc_lm_df = NULL
#   for(idx in 1:length(disc_lm)) {   disc_lm_df = rbind(disc_lm_df, disc_lm[[idx]])   }
#   write.csv(disc_lm_df, sprintf('%s/AY/%s/initiallinreg/all_discovery_linereg.csv', save_dir, subfoldername))
# }
# 
# if('positive_linereg.csv' %in% list.files(sprintf('%s/AY/%s/initiallinreg/', save_dir, subfoldername))) {
#   print(sprintf("[%s]        - Loading lm on all positive pairs", Sys.time()))
#   pos_lm_df = read.csv(sprintf('%s/AY/%s/initiallinreg/positive_linereg.csv', save_dir, subfoldername))
# } else {
#   print(sprintf("[%s]        - Fitting lm on all positive pairs", Sys.time()))
#   t0 = Sys.time()
#   pos_lm = future.apply::future_mapply(fit_lm,
#                                        AY_idx = 1:nrow(pos_tests_df), 
#                                        MoreArgs = list(AY_df = pos_tests_df),
#                                        future.globals = TRUE,
#                                        future.seed = 123456, 
#                                        SIMPLIFY = FALSE)
#   t1 = Sys.time()
#   print(sprintf("[%s]        - %2.2f mins", Sys.time(), difftime(t1, t0, units = "mins")))
#   
#   pos_lm_df = NULL
#   for(idx in 1:length(pos_lm)) {   pos_lm_df = rbind(pos_lm_df, pos_lm[[idx]])   }
#   write.csv(pos_lm_df, sprintf('%s/AY/%s/initiallinreg/positive_linereg.csv', save_dir, subfoldername))
# }
# 



# --------  \__choose discovery tests  ------ 
# NUM_DISC_TESTS = 2000
# NUM_DISC_TESTS = nrow(disc_lm_df)
disc_tests_df = expand.grid(grna=all_As, gene=all_Ys, stringsAsFactors = FALSE) |> dplyr::filter(grna != gene) |> mutate(type='discovery')
# disc_tests_df = disc_tests_df[sample(nrow(disc_tests_df), size = NUM_DISC_TESTS), ] |> mutate(type = 'discovery')

# # take the top?? or take the top per grna?
# take_top_disc = 'overall'
# # take_top_disc = 'bygrna'
# if(take_top_disc == 'overall') {
#   disc_tests_df = disc_lm_df |> arrange(pval) |> head(NUM_DISC_TESTS)
# 
#   # disc_tests_df$grna |> table()
# } else if (take_top_disc == 'bygrna') {
#   num_grna = length(unique(disc_tests_df$grna))
#   disc_tests_df = disc_lm_df |> group_by(grna) |> arrange(pval) |> slice_head(n=ceiling(NUM_DISC_TESTS / num_grna))
# 
#   # disc_tests_df$grna |> table()
#   rm(num_grna)
# }
# disc_tests_df = ungroup(disc_tests_df) |> select(grna, gene) |> mutate(type = 'discovery')
# 
# disc_tests_df = disc_lm_df |> select(grna, gene) |> mutate(type = 'discovery')


# --------  \__assemble AY tests together and save  ------ 
AY = rbind(pos_tests_df, disc_tests_df)


write.csv(x = AY, file = sprintf('%s/AY/%s/AY.csv', save_dir, subfoldername), row.names = FALSE) # csv 



rm(GRNA_SAMPLESIZE_MIN, GENE_NUMNONZEROCELLS_MIN, GENE_IMPORTANCERANK_MAX, NUM_DISC_TESTS)
rm(all_As, all_Ys, disc_tests_df, pos_tests_df, gene_dev_df) # prob keep all_Ys to have a collection of genes to use as NCs
rm(pos_tests_grna, pos_tests_alias)

# all(validsinglegenes %in% row.names(gene))
# all(validsinglegenes %in% row.names(gene_norm))
# all(AY$gene %in% row.names(gene))
# AY$gene[!AY$gene %in% row.names(gene)] # make sure positive tests' genes names are modified if they have alias (no "Mll1"  "Myst4" "Scn2a")


# ============================================================================== #
# =================== ESTIMATE EFFECTS (parallel) ================================
# ============================================================================== #
dir.create(sprintf('%s/AY/%s/intermediateATEs/', save_dir, subfoldername))




NCs_pca = readRDS(sprintf('%s/pca/NCloadings.rds', save_dir)) |> as.matrix() # PCA loadings
# numNCs = c(1, 3, 5, 10, 15)
# numNCs = c(1, 5) # test w smaller set of numNCs
numNCs = c(1, 3, 5, 8, 10)
libsize_log = log((metadata$nUMI |> as.numeric())+ 1e-20)
# control_cell_idx = which(colSums(grna) == 0) # cells with GFP- DIFFERENT!
control_cell_idx = which(grna['Gfp', ] == 1) # cells with GFP
# control_cell_idx2 = which(metadata$Perturbation == 'GFP') # cells with GFP  table(control_cell_idx == control_cell_idx2)


# double check cells are in the same order
(colnames(gene) == colnames(grna)) |> all()
(colnames(gene) == colnames(gene_norm)) |> all()
(colnames(gene) == metadata$NAME) |> all()
(colnames(gene) == row.names(NCs_pca)) |> all()



# remove all non-relevant objects in environment
rm(list = setdiff(objects(),
           c('grna', 'gene', 'gene_norm', 'AY',
             'libsize_log', 'control_cell_idx', 'NCs_pca', 'validsinglegenes', 'numNCs',
             'DEVICE', 'save_dir', 'util_dir', 'data_dir', 'subfoldername')))



#' estimate ATEs/logFoldChanges/Coefficients using loaded in vars in environment
#' Specifically: 
#'    grna, gene, gene_norm, AY, validsinglegenes
#'    libsize_log, control_cell_idx, NCs, numNCs
#'    
#' res is dataframe with cols:
#' - grna, gene, type 
#' - ... etc ...
#' - tstat: test statistic either the t value (lm) or z value (glm, proximal)
#' @param AY_idx (integer) index of the AY test 
#' 
estimate_effects <- function(AY_idx, save_intermediateATEs=TRUE) {
  res = NULL
  A_name = AY[AY_idx, 'grna']
  Y_name = AY[AY_idx, 'gene']
  
  
  
  # ===== assemble df =====
  grna_cell_idx = which(grna[A_name, ] == 1)
  AY_data_idx = c(control_cell_idx, grna_cell_idx)
  
  df = data.frame(# A2 = c(rep(0, length(control_cell_idx)), rep(1, length(grna_cell_idx))), # another way to get 0/1
    A           = grna[A_name, AY_data_idx],
    Y           = gene[Y_name, AY_data_idx], 
    Y_norm      = gene_norm[Y_name, AY_data_idx],
    libsize_log = libsize_log[AY_data_idx])
  
  
  

  
  
  # ===== estimate =====
  # ===== \__Linear Model without any covariates =====
  t0 = Sys.time()
  lmfit = lm(formula='Y_norm ~ A', data = df)
  lmfit_summary = summary(lmfit)
  t1 = Sys.time()
  
  res = bind_rows(res, 
                  data.frame(
                    NC_type     = NA,
                    method      = 'lm',
                    method_type = 'glm',
                    numNC       = NA,
                    ATE = lmfit_summary$coefficients['A', 'Estimate'],
                    se  = lmfit_summary$coefficients['A', 'Std. Error'],
                    tstat= lmfit_summary$coefficients['A', 't value'],
                    pval= lmfit_summary$coefficients['A', 'Pr(>|t|)'],
                    time_sec = difftime(t1, t0, units = 'secs') |> as.numeric()))
  rm(lmfit, lmfit_summary, t1, t0)
  
  
  # ===== \__ Negative Binomial Fit without any covariates ====
  t0 = Sys.time()
  nbfit = MASS::glm.nb(formula='Y ~ A + offset(libsize_log)', data = df)
  nbfit_summary = summary(nbfit)
  t1 = Sys.time()
  
  res = bind_rows(res, 
                  data.frame(
                    NC_type     = NA,
                    method      = 'negbin',
                    method_type = 'glm',
                    numNC       = NA,
                    ATE = nbfit_summary$coefficients['A', 'Estimate'],
                    se  = nbfit_summary$coefficients['A', 'Std. Error'],
                    tstat= nbfit_summary$coefficients['A', 'z value'],
                    pval= nbfit_summary$coefficients['A', 'Pr(>|z|)'],
                    time_sec = difftime(t1, t0, units = 'secs') |> as.numeric()))
  rm(nbfit, nbfit_summary, t1, t0)
  

  # ===== \__ Proximal PCA ====
  # Proximal with PCA loadings
  dfZ = NCs_pca[AY_data_idx, seq(from = 2, to = min(2*max(numNCs), ncol(NCs_pca)), by = 2)] # evens which((1:ncol(NCs)) %% 2 == 0)
  dfW = NCs_pca[AY_data_idx, seq(from = 1, to = min(2*max(numNCs), ncol(NCs_pca)), by = 2)] # odds
  colnames(dfZ) = paste0('Z', 1:ncol(dfZ))
  colnames(dfW) = paste0('W', 1:ncol(dfW))
  for(numNCs_ in numNCs)  {
    # ===== \__\__Proximal Linear ====
    t0 = Sys.time()
    pci2s_res = tryCatch({  # give pci2s results if works
      pci2s::p2sls.lm(
        Y = df$Y_norm, 
        A = df$A, 
        W = dfZ[,1:numNCs_] |> as.data.frame(), 
        Z = dfW[,1:numNCs_] |> as.data.frame(), 
        variance = TRUE)
    },
    error = function(cond) { # give NA results if errored
      return(list(summary_second_stage=
                    matrix(c(NA, NA, NA, NA), nrow=1, ncol=4, 
                           dimnames = list(c('A'), c('Estimate', 'Std. Error', 'z value', 'Pr(>|z|)')))))
    })
    t1 = Sys.time()
    
    res = bind_rows(res, 
                    data.frame(
                      NC_type     = 'PCA',
                      method      = 'proximallinear',
                      method_type = 'proximal',
                      numNC       = numNCs_,
                      ATE  = pci2s_res$summary_second_stage['A', 'Estimate'],
                      se   = pci2s_res$summary_second_stage['A', 'Std. Error'],
                      tstat= pci2s_res$summary_second_stage['A', 'z value'],
                      pval = pci2s_res$summary_second_stage['A', 'Pr(>|z|)'],
                      time_sec = difftime(t1, t0, units = 'secs') |> as.numeric()))
    
    rm(pci2s_res, t0, t1)
    
    # ===== \__\__Proximal Negative Binomial ====
    if(F){   # do not run
    t0 = Sys.time()
    pci2s_res = tryCatch({  # give pci2s results if works
      pci2s::p2sls.negbin(
        Y = df$Y, 
        A = df$A, 
        W = dfZ[,1:numNCs_] |> as.data.frame(), 
        Z = dfW[,1:numNCs_] |> as.data.frame(), 
        offset = df$libsize_log,
        nco_type = rep("linear", numNCs_),
        # nco_args = lapply(X = 1:num_NCs, FUN = function(x){list(init=NA, offset=log(df_all$library_size))}), # No offsets for continuous versions
        variance = TRUE,
        verbose = FALSE)
    },
    error = function(cond) { # give NA results if errored
      return(list(summary_second_stage=
                    matrix(c(NA, NA, NA, NA), nrow=1, ncol=4, 
                           dimnames = list(c('A'), c('Estimate', 'Std. Error', 'z value', 'Pr(>|z|)')))))
    })
    t1 = Sys.time()
    
    res = bind_rows(res, 
                    data.frame(
                      NC_type     = 'PCA',
                      method      = 'proximalnegbin',
                      method_type = 'proximal',
                      numNC       = numNCs_,
                      ATE  = pci2s_res$summary_second_stage['A', 'Estimate'],
                      se   = pci2s_res$summary_second_stage['A', 'Std. Error'],
                      tstat= pci2s_res$summary_second_stage['A', 'z value'],
                      pval = pci2s_res$summary_second_stage['A', 'Pr(>|z|)'],
                      time_sec = difftime(t1, t0, units = 'secs') |> as.numeric()))
    
    rm(pci2s_res, t0, t1)
    }
    
  }
  rm(dfZ, dfW, numNCs_)
  
  
  
  # ===== \__ Proximal singlegenes ====
  RUN_PROXIMAL_SINGLEGENES = T
  if(RUN_PROXIMAL_SINGLEGENES) {
    # ===== \__\__Proximal Linear ====
    # Proximal with singlegenes
    NC_names = setdiff(validsinglegenes, c(A_name, Y_name)) # set of genes to use as NCs (not outcome and not grna target [assume grna name is the targeted gene name])
    NC_names = sample(NC_names, 2*max(numNCs))
    # Z_names = NC_names[1:max(numNCs)]
    # W_names = NC_names[(max(numNCs)+1):(2*max(numNCs))]
    
    NCs_singlegene = gene_norm[NC_names, AY_data_idx] |> t() |> as.data.frame()
    dfZ = NCs_singlegene[ , seq(from = 2, to = min(2*max(numNCs), ncol(NCs_singlegene)), by = 2)] # evens which((1:ncol(NCs)) %% 2 == 0)
    dfW = NCs_singlegene[ , seq(from = 1, to = min(2*max(numNCs), ncol(NCs_singlegene)), by = 2)] # odds
    colnames(dfZ) = paste0('Z', 1:ncol(dfZ))
    colnames(dfW) = paste0('W', 1:ncol(dfW))
    for(numNCs_ in numNCs)  {
      t0 = Sys.time()
      pci2s_res = tryCatch({  # give pci2s results if works
        pci2s::p2sls.lm(
          Y = df$Y_norm, 
          A = df$A, 
          W = dfZ[,1:numNCs_], 
          Z = dfW[,1:numNCs_], 
          variance = TRUE)
      },
      error = function(cond) { # give NA results if errored
        return(list(summary_second_stage=
                      matrix(c(NA, NA, NA, NA), nrow=1, ncol=4, 
                             dimnames = list(c('A'), c('Estimate', 'Std. Error', 'z value', 'Pr(>|z|)')))))
      })
      t1 = Sys.time()
      
      res = bind_rows(res, 
                      data.frame(
                        NC_type     = 'singlegene',
                        method      = 'proximallinear',
                        method_type = 'proximal',
                        numNC       = numNCs_,
                        ATE  = pci2s_res$summary_second_stage['A', 'Estimate'],
                        se   = pci2s_res$summary_second_stage['A', 'Std. Error'],
                        tstat= pci2s_res$summary_second_stage['A', 'z value'],
                        pval = pci2s_res$summary_second_stage['A', 'Pr(>|z|)'],
                        time_sec = difftime(t1, t0, units = 'secs') |> as.numeric()))
      
      rm(pci2s_res, t0, t1)
      
      # ===== \__\__Proximal Negative Binomial ====
      if(F){ # do not run
      t0 = Sys.time()
      pci2s_res = tryCatch({  # give pci2s results if works
        pci2s::p2sls.negbin(
          Y = df$Y, 
          A = df$A, 
          W = dfZ[,1:numNCs_] |> as.data.frame(), 
          Z = dfW[,1:numNCs_] |> as.data.frame(), 
          offset = df$libsize_log,
          nco_type = rep("linear", numNCs_),
          # nco_args = lapply(X = 1:num_NCs, FUN = function(x){list(init=NA, offset=log(df_all$library_size))}), # No offsets for continuous versions
          variance = TRUE,
          verbose = FALSE)
      },
      error = function(cond) { # give NA results if errored
        return(list(summary_second_stage=
                      matrix(c(NA, NA, NA, NA), nrow=1, ncol=4, 
                             dimnames = list(c('A'), c('Estimate', 'Std. Error', 'z value', 'Pr(>|z|)')))))
      })
      t1 = Sys.time()
      
      res = bind_rows(res, 
                      data.frame(
                        NC_type     = 'singlegene',
                        method      = 'proximalnegbin',
                        method_type = 'proximal',
                        numNC       = numNCs_,
                        ATE  = pci2s_res$summary_second_stage['A', 'Estimate'],
                        se   = pci2s_res$summary_second_stage['A', 'Std. Error'],
                        tstat= pci2s_res$summary_second_stage['A', 'z value'],
                        pval = pci2s_res$summary_second_stage['A', 'Pr(>|z|)'],
                        time_sec = difftime(t1, t0, units = 'secs') |> as.numeric()))
      
      rm(pci2s_res, t0, t1)
      }
    }
    rm(dfZ, dfW, numNCs_)
  }

  
  # ===== return ==== 
  # n_trtmnt = length(grna_cell_idx) # sample size of treatments (# of cells receiving this perturbation)
  res = res |> mutate(AY_idx = AY_idx, grna = A_name, gene = Y_name, type = AY[AY_idx, 'type'], n_trtmnt = length(grna_cell_idx), .before = 1)
  
  if(save_intermediateATEs) {
    write.csv(res, sprintf('%s/AY/%s/intermediateATEs/ATE%s.csv', save_dir, subfoldername, AY_idx))
  }
  return(res)
}

estimate_effects_errorhandling <- function(AY_idx) {
  res_df = tryCatch({estimate_effects(AY_idx=AY_idx)},
                    error = function(cond) {
                      # message(sprintf('Error est ATE w pci2s with %s', 
                      #                 AY_idx)) 
                      return(NULL)
                    })
  # if errored, return NULL
  if(is.null(res_df)) {
    return(NULL)
  } else {
    return(res_df)
  }
}




# test a run
# estimate_effects(4)
# estimate_effects_errorhandling(4)


# =================== \__Get ATEs (parallel) ==========================================================
print(sprintf("[%s]    - Get Estimates (parallel)", Sys.time()))

# whichROWS = 1:50
whichROWS = 1:nrow(AY)

t0 = Sys.time()
effects_par = future.apply::future_mapply(estimate_effects_errorhandling,
                                      AY_idx = whichROWS, 
                                      future.globals = TRUE,
                                      future.seed = 123456, 
                                      SIMPLIFY = FALSE)
t1 = Sys.time()
print(sprintf("[%s]        - %2.2f mins", Sys.time(), difftime(t1, t0, units = "mins")))


saveRDS(effects_par, file = sprintf('%s/AY/%s/effects_par.rds', save_dir, subfoldername)) # save this parallel res as rds



# =================== \__Combine ATEs (into one df) ===================================================
print(sprintf("[%s]    - Combine ATEs", Sys.time()))


effects_df = NULL
for(whichROWS_idx in 1:length(whichROWS)) { # whichROWS_idx is always 1, 2, ...
  cur_effects_par = effects_par[[whichROWS_idx]]
  
  if(!is.null(cur_effects_par)) {
    # AY test is whichROWS[whichROWS_idx]  = AY_idx!!! 
    effects_df = dplyr::bind_rows(effects_df, cur_effects_par)

  }
  rm(cur_effects_par)
}
write.csv(x = effects_df, file = sprintf('%s/AY/%s/effects.csv', save_dir, subfoldername), row.names = FALSE)






# ==== TRASH ==== TRASH ==== TRASH ==== TRASH ==== TRASH ==== TRASH ==== TRASH 
# ============================= TRASH ========================================
# ==== TRASH ==== TRASH ==== TRASH ==== TRASH ==== TRASH ==== TRASH ==== TRASH 



# effects_df |> select(AY_idx, type) |> distinct() |> pull(type) |> table()
# effects_df |> filter(method == 'lm' & type == 'positive') |> pull(pval) |> hist(breaks = 100)
# effects_df |> filter(method == 'negbin' & type == 'positive') |> pull(pval) |> hist(breaks = 100)
# 
# effects_df |> filter(method == 'proximallinear' & type == 'positive' & numNC == 5) |> pull(pval) |> hist(breaks = 100)
# effects_df |> filter(method == 'proximallinear' & type == 'discovery' & numNC == 5) |> pull(pval) |> hist(breaks = 100)
# effects_df |> filter(method == 'negbin' & type == 'positive') |> pull(pval) |> hist(breaks = 100)



# ==== TRASH ==== TRASH ==== TRASH ==== TRASH ==== TRASH ==== TRASH ==== TRASH 
# ============================= TRASH ========================================
# ==== TRASH ==== TRASH ==== TRASH ==== TRASH ==== TRASH ==== TRASH ==== TRASH 




# ==== TRASH ==== TRASH ==== TRASH ==== TRASH ==== TRASH ==== TRASH ==== TRASH 
# ============================= TRASH ========================================
# ==== TRASH ==== TRASH ==== TRASH ==== TRASH ==== TRASH ==== TRASH ==== TRASH 

if(F) {
  
  # ============================================================================== #
  # =================== Some Estimates (seq) =====================================
  # ============================================================================== #
  
  NCs = readRDS(sprintf('%s/pca/NCloadings.rds', save_dir))  # PCA loadings
  numNCs = c(1, 3, 5, 10, 15, 20)
  
  libsize_log = log((metadata$nUMI |> as.numeric())+ 1e-20)
  
  
  control_cell_idx = which(colSums(grna) == 0) # cells with GFP
  # control_cell_idx2 = which(metadata$Perturbation == 'GFP') # cells with GFP  table(control_cell_idx == control_cell_idx2)
  
  res_all = NULL
  # whichROWS = 1:nrow(AY)
  whichROWS = 1:5
  for(AY_idx in whichROWS) {
    
    res = NULL
    A_name = AY[AY_idx, 'grna']
    Y_name = AY[AY_idx, 'gene']
    
    # assemble df
    grna_cell_idx = which(grna[A_name, ] == 1)
    AY_data_idx = c(control_cell_idx, grna_cell_idx)
    
    df = data.frame(# A2 = c(rep(0, length(control_cell_idx)), rep(1, length(grna_cell_idx))), # another way to get 0/1
      A           = grna[A_name, AY_data_idx],
      Y           = gene[Y_name, AY_data_idx], 
      Y_norm      = gene_norm[Y_name, AY_data_idx],
      libsize_log = libsize_log[AY_data_idx])
    
    # Negative Binomial Fit without any covariates
    formula = 'Y ~ A + offset(libsize_log)'
    t0 = Sys.time()
    nbfit = MASS::glm.nb(formula=formula, data = df)
    nbfit_summary = summary(nbfit)
    t1 = Sys.time()
    
    res = bind_rows(res, 
                    data.frame(
                      NC_type     = NA,
                      method      = 'negbin',
                      method_type = 'glm',
                      numNC       = NA,
                      ATE = nbfit_summary$coefficients['A', 'Estimate'],
                      se  = nbfit_summary$coefficients['A', 'Std. Error'],
                      zval= nbfit_summary$coefficients['A', 'z value'],
                      pval= nbfit_summary$coefficients['A', 'Pr(>|z|)'],
                      time_sec = difftime(t1, t0, units = 'secs') |> as.numeric()))
    
    
    # Proximal with PCA loadings
    dfZ = NCs[AY_data_idx, seq(from = 2, to = min(2*max(numNCs), ncol(NCs)), by = 2)] # evens which((1:ncol(NCs)) %% 2 == 0)
    dfW = NCs[AY_data_idx, seq(from = 1, to = min(2*max(numNCs), ncol(NCs)), by = 2)] # odds
    colnames(dfZ) = paste0('Z', 1:ncol(dfZ))
    colnames(dfW) = paste0('W', 1:ncol(dfW))
    
    for(numNCs_ in numNCs)  {
      t0 = Sys.time()
      pci2s_res = pci2s::p2sls.lm(
        Y = df$Y_norm, 
        A = df$A, 
        W = dfZ[,1:numNCs_], 
        Z = dfW[,1:numNCs_], 
        variance = TRUE)
      t1 = Sys.time()
      
      res = bind_rows(res, 
                      data.frame(
                        NC_type     = 'PCA',
                        method      = 'proximal',
                        method_type = 'proximal',
                        numNC       = numNCs_,
                        ATE = pci2s_res$summary_second_stage['A', 'Estimate'],
                        se  = pci2s_res$summary_second_stage['A', 'Std. Error'],
                        zval= pci2s_res$summary_second_stage['A', 'z value'],
                        pval= pci2s_res$summary_second_stage['A', 'Pr(>|z|)'],
                        time_sec = difftime(t1, t0, units = 'secs') |> as.numeric()))
      
      rm(pci2s_res, t0, t1)
    }
    rm(dfZ, dfW)
    
    res_all = bind_rows(res_all, 
                        res |> mutate(AY_idx = AY_idx, grna = A_name, gene = Y_name, type = AY[AY_idx, 'type'], .before = 1))
    
  }
  
  # ============================================================================== #
  # =================== Negative Binomial GLM (for a lot) ========================
  # ============================================================================== #
  
  # sometimes the names are bad for the glm formulas (e.g. AB-C). Remove non letter and number characters 
  clean_name <- function(x) {  gsub(pattern='[^a-zA-Z0-9]', replacement='', x=x)  }
  
  
  
  # Neg Bin fits without any covariates
  libsize = data.frame(libsize_log=log(metadata$nUMI |> as.numeric() + 1e-20))
  # libsize_log = log((metadata$nUMI |> as.numeric())+ 1e-20)
  
  
  
  
  # which_genes = colnames(Y) # all genes
  N_genes = 1000   # or randomly choose some
  set.seed(12345)
  
  perturbations_names = unique(metadata$Perturbation) |> sort()
  pert_targetgene_in_subsetdata = intersect(perturbations_names, colnames(Y)) # always test these bc they are pert targets
  which_genes = union(pert_targetgene_in_subsetdata, 
                      sample(setdiff(colnames(Y), perturbations_names),  replace = F, size = N_genes - length(pert_targetgene_in_subsetdata)))
  
  # fit and save res: 1 glm per gene (vs individual pert-gene fits)
  res = NULL
  for(Y_name in which_genes) {
    
    df = cbind(Y[, Y_name, drop=FALSE] |> `colnames<-`(paste0('gene_', clean_name(Y_name))), 
               libsize,
               A
    )
    
    
    # formula = paste0('gene_', Y_name, ' ~ offset(libsize_log) + ', paste0(colnames(A), collapse = ' + '))
    formula = paste0('gene_', clean_name(Y_name), 
                     ' ~ offset(libsize_log) + ', 
                     paste0(sapply(colnames(A), clean_name), collapse = ' + '))
    nbfit = MASS::glm.nb(formula=formula, data = df)
    
    res = rbind(res, 
                summary(nbfit)$coefficients[-1, ] |> 
                  `colnames<-`(c('estimate', 'se', 'zvalue', 'pvalue')) |> 
                  data.frame() |> 
                  dplyr::mutate(grna=colnames(A), gene = Y_name, .after="estimate") )
    
    
    # summary(nbfit)$coefficients[colnames(A), ] |> 
    #   `colnames<-`(c('estimate', 'se', 'zvalue', 'pvalue')) |> 
    #   data.frame() |> 
    #   tibble::rownames_to_column("grna") |> 
    #   dplyr::mutate(gene = Y_name, .after="grna")
  }
  
  
  
  
  
}

