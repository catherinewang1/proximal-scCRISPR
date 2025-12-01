# ---------------------------------------------------------------------------- #
#           Perform Oracle Estimates- using measured confounders               #
# Oracle Estimators:                                                           #
# - simple lm(Y ~ A + X)                                                       #
# - SCEPTRE (using raw counts)                                                 #
# SCEPTRE is the ideal estimator, but SCEPTRE uses raw counts and assumes a    #
# non-identity link function (for NB or Pois?). So, we should really compare   #
# these linear proximal methods                                                #
# ---------------------------------------------------------------------------- #
args = commandArgs(trailingOnly = TRUE)
args = c('laptop', '2000', 'E_PCA')

require(assertthat) # for some assert statements
library(Matrix)
library(dplyr)
library(tibble)
library(ggplot2)    # plotting
library(cowplot)


library(future.apply)
options(future.globals.maxSize= 850*1024^2) #1st num is MB
plan(multisession, workers = 8)
# plan(sequential)

theme_set(theme_cowplot() +
            theme(plot.title = element_text(hjust = .5),
                  plot.subtitle = element_text(hjust = .5),
                  strip.background = element_rect(color = 'black', fill = 'white')))

DEVICE = args[1]
source('../PATHS.R') # load in data_dir and save_dir and CODE_DIR, depending on DEVICE value

assertthat::assert_that(!is.null(data_dir), msg='first arg must be: laptop, desktop, or ubergenno')

assertthat::assert_that(length(args) > 1, msg="must give arg for specifying num imp genes 'Rscript <filename>.R ubergenno 4000'")
NUM_IMPORTANTGENES = as.integer(args[2]) # should be max importance from setting... TODO: change to remove this input. extract this value from saved setting

assertthat::assert_that(length(args) > 2, msg="must give arg for specifying chosen AYZW name 'Rscript <filename>.R ubergenno C'")
AYZW_setting_name = args[3]

plot_savepath = sprintf('%s/oracle/%s/', save_dir, AYZW_setting_name)
dir.create(plot_savepath, recursive = TRUE)


# =================== Start ====================================================
print(sprintf("[%s] START: Oracles (continuous Y)", Sys.time()))


source(sprintf('%s/CBEstAll.R', util_dir)) # for functions to estimate here

# load chosen AYZW names
AY   = read.csv(sprintf('%s/cbgenes/%s/AY.csv', save_dir, AYZW_setting_name))
AYZW = readRDS(sprintf('%s/cbgenes/%s/AYZW.rds', save_dir, AYZW_setting_name))

# load gene importance info
imp_gene_names = readRDS(sprintf('%s/important_genes_name.rds', save_dir))
imp_gene_idx   = readRDS(sprintf('%s/important_genes_idx.rds',  save_dir))
imp_gene = data.frame(gene     = imp_gene_names,
                      gene_idx = imp_gene_idx,
                      gene_imp_rank = 1:length(imp_gene_names))


# create gene and grna ondisc managers
gene_odm <- ondisc::read_odm(odm_fp      = paste0(data_dir, "/papalexi-2021/processed/gene/expression_matrix.odm"),
                             metadata_fp = paste0(data_dir, "/papalexi-2021/processed/gene/metadata.rds"))
grna_odm <- ondisc::read_odm(odm_fp      = paste0(data_dir, "/papalexi-2021/processed/grna_assignment/assignment_matrix.odm"),
                             metadata_fp = paste0(data_dir, "/papalexi-2021/processed/grna_assignment/metadata.rds"))

# load grna assignments (load all into memory)
grna = grna_odm[[,1:ncol(grna_odm)]] # |> as.matrix() # ~110 x 20729 = #grnas x #cells
grna_rownames = grna_odm |> ondisc::get_feature_covariates() |> rownames()

# load normalized gene exp
h5file      = paste0(save_dir, "/gene.h5"); print(h5file)
reading_hd5file  = rhdf5::H5Fopen(name = h5file)
readin_gene_norm = reading_hd5file&'gene_norm'
gene_norm = readin_gene_norm[1:NUM_IMPORTANTGENES, ] # dim = 4000 x 20729 = #important x #cells
rhdf5::h5closeAll()
invisible(gc(verbose=FALSE))



# =================== Define File Specific Functions ======================================
print(sprintf("[%s]    - Define File Specific Functions", Sys.time()))

# All the NT grna idx
NT_names = grna_odm |> ondisc::get_feature_covariates() |>
  filter(target_type == 'non-targeting') |> rownames()
NT_idx = which(apply(X = grna_odm[[NT_names, ]], MARGIN = 2, FUN = sum) > 0)



cell_covariates = gene_odm |> ondisc::get_cell_covariates()
# dim(cell_covariates); head(cell_covariates)


#' make function that gets the lmYAX dataframe (trtmt A, resp Y, confdrs X)
get_lmYAXdf_make <- function(AY, AYZW, 
                           gene_norm, grna_rownames, grna, cell_covariates,
                           NT_idx, imp_gene_names) {
  get_importance_rank = get_importance_rank_make(imp_gene_names) 
  
  get_lmYAX <- function(AY_idx) {
    # Get treatment vector
    AY_row = AY[AY_idx, ]
    A_grna_idx = which(grna_rownames == AY_row$A)
    A_idx = which(as.logical(grna[A_grna_idx, 1:ncol(grna)]))
    # subset cells of 'treated' (w A grna) and 'control' (NT gran)
    A = grna[A_grna_idx, c(A_idx, NT_idx)]
    
    # Get response vector
    Y = gene_norm[get_importance_rank(AY_row$Y), c(A_idx, NT_idx)]
    
    # Get covariates/confounders 
    X = cell_covariates[c(A_idx, NT_idx), ]
    
    df = data.frame(Y, A, X)
    
    return(df)
  }
  
  return(get_lmYAX)
  
}

#' make function that gets the lmYAX dataframe (trtmt A, resp Y, confdrs X)
get_lmYAXdf_counts_make <- function(AY, AYZW, 
                             gene_odm, grna_rownames, grna, cell_covariates,
                             NT_idx, imp_gene_names) {
  
  get_importance_rank = get_importance_rank_make(imp_gene_names) 
  
  get_lmYAX_counts <- function(AY_idx) {
    # Get treatment vector
    AY_row = AY[AY_idx, ]
    A_grna_idx = which(grna_rownames == AY_row$A)
    A_idx = which(as.logical(grna[A_grna_idx, 1:ncol(grna)]))
    # subset cells of 'treated' (w A grna) and 'control' (NT gran)
    A = grna[A_grna_idx, c(A_idx, NT_idx)]
    
    
    # Get response vector
    Y = gene_odm[[AY_row$Y, c(A_idx, NT_idx)]] |> as.vector()
    
    
    # Get covariates/confounders 
    X = cell_covariates[c(A_idx, NT_idx), ]
    
    df = data.frame(Y, A, X)
    
    return(df)
  }
  
  return(get_lmYAX_counts)
}

#' Get resulting ATE and p-value estimate using the confounders
#' (includes unadjusted, biased lmYA, and confounder adjusted lmYAX)
#' 
getlmYAX_lmresult <- function(df) {
  df$A = as.integer(df$A)
  lmYA  = summary(lm(Y ~ A, df))$coefficients
  lmYAX = summary(lm(Y ~ ., df))$coefficients
  
  
  # lmYA: extract ATE and pval if present
  if('A' %in% rownames(lmYA)) {
    lmYAdf = data.frame(lmYA_ATE  =  lmYA['A', 'Estimate'],
                        lmYA_pval =  lmYA['A', 'Pr(>|t|)'])
  } else {
    lmYAdf = data.frame(lmYA_ATE=NA,
                        lmYA_pval=NA)
  }
  
  
  # lmYAX: extract ATE and pval if present
  if('A' %in% rownames(lmYAX)) {
    lmYAXdf = data.frame(lmYAX_ATE  =  lmYAX['A', 'Estimate'],
                         lmYAX_pval =  lmYAX['A', 'Pr(>|t|)'])
  } else {
    lmYAXdf = data.frame(lmYAX_ATE=NA,
                         lmYAX_pval=NA)
  }
  
  return(cbind(lmYAdf, lmYAXdf))
}


# library(npcausal)

#' get ATEs using 
getlmYAX_npcausalresult <- function(df) {
  
  # df = my_getlmYAXdf(3)
  npres = npcausal::ate(y=df$Y,
                        a=as.integer(df$A),
                        # x = df |> dplyr::select(n_umis, p_mito),
                        x=df |> dplyr::select(-Y, -A), 
                        nsplits = 1)
  
  data.frame(npcausal_ATE = npres$res[3, 'est'],
             npcausal_pval= npres$res[3, 'pval'])
  
}






my_getlmYAXdf = get_lmYAXdf_make(AY=AY, AYZW=AYZW, 
                                 gene_norm=gene_norm, 
                                 grna_rownames=grna_rownames, 
                                 grna=grna, 
                                 cell_covariates=cell_covariates,
                                 NT_idx=NT_idx, 
                                 imp_gene_names=imp_gene_names)
my_getlmYAXcountsdf = get_lmYAXdf_counts_make(AY=AY, AYZW=AYZW, 
                                 gene_odm=gene_odm, 
                                 grna_rownames=grna_rownames, 
                                 grna=grna, 
                                 cell_covariates=cell_covariates,
                                 NT_idx=NT_idx, 
                                 imp_gene_names=imp_gene_names)




ATEs_lm = NULL
ATEs_npcausal = NULL
for(AY_idx in 1:nrow(AY)) {
  dfi = my_getlmYAXdf(AY_idx = AY_idx)
  
  ATEs_lm = rbind(ATEs_lm, 
                  cbind(data.frame(AY_idx=AY_idx),
                        getlmYAX_lmresult(dfi)))
  
  ATEs_npcausal = rbind(ATEs_npcausal, 
                        suppressWarnings({
                          cbind(data.frame(AY_idx = AY_idx), getlmYAX_npcausalresult(dfi))
                        }))
}


ATEs = merge(ATEs_lm, ATEs_npcausal,
             by = 'AY_idx')




AY_oracle = merge(AY |> dplyr::mutate(AY_idx = 1:nrow(AY)), ATEs, 
                  by = 'AY_idx')

AY_oracle |> dplyr::filter(type == 'positive') |> pull(lmYA_pval) |> hist()
AY_oracle |> dplyr::filter(type == 'positive') |> pull(lmYAX_pval) |> hist()

AY_oracle |> dplyr::filter(type == 'negative') |> pull(lmYA_pval) |> hist()
AY_oracle |> dplyr::filter(type == 'negative') |> pull(lmYAX_pval) |> hist()



AY_oracle |> dplyr::filter(type == 'negative') |> pull(lmYA_ATE) |> hist()
AY_oracle |> dplyr::filter(type == 'negative') |> pull(lmYAX_ATE) |> hist()

AY_oracle_tall = AY_oracle |> reshape2::melt(id.vars = c('type', 'A', 'A_chr', 'Y', 'Y_chr'), 
                            measure.vars = c('lmYA_pval', 'lmYAX_pval'),
                            variable.name = 'method', 
                            value.name = 'pval')
head(AY_oracle_tall)
AY_oracle_tall$method |> table()
dim(AY_oracle)




ggplot(AY_oracle_tall) +
  geom_jitter(aes(x = pval, y = method, 
                 color = type))

AY_oracle |> dplyr::filter(type == 'negative') |> mutate(ratio = abs(lmYAX_ATE/lmYA_ATE)) |> pull(ratio) |> hist(breaks = 0:25)

AY_oracle |> dplyr::filter(lmYAX_pval < .05 & lmYAX_ATE > .05) |>
  dplyr::arrange(lmYAX_pval)


#                              Y    A n_nonzero n_umis  lane bio_rep phase     p_mito
# l1_ATCATCTCATTGCGGC -2.9587450 TRUE      5150  24715 Lane1   rep_1   G2M 0.07274934
# l1_CGATGGCAGCTGTTCA  0.3500470 TRUE      3332   9816 Lane1   rep_1    G1 0.05093725
head(dfi)


ggplot(dfi) +
  # geom_point(aes(x = p_mito, y = Y,
  #            color = phase), alpha = .5) +
  geom_boxplot(aes(x = lane, y = p_mito))


lm(Y ~ A + n_nonzero + n_umis * lane + p_mito * lane, dfi) |> summary() |> coefficients()



write.csv(x = AY_oracle, file = sprintf('%sAY_oracle.csv', plot_savepath))
oracle_pvals = AY_oracle |> reshape2::melt(id.vars = c('type', 'A', 'A_chr', 'Y', 'Y_chr'), 
                                    measure.vars = c('lmYA_pval', 'lmYAX_pval', 'npcausal_pval'),
                                    variable.name = 'method', 
                                    value.name = 'pval')

oracle_ATEs = AY_oracle |> reshape2::melt(id.vars = c('type', 'A', 'A_chr', 'Y', 'Y_chr'), 
                                   measure.vars = c('lmYA_ATE', 'lmYAX_ATE', 'npcausal_ATE'),
                                   variable.name = 'method', 
                                   value.name = 'ATE')
# Nicer Labels for type
type_labeller = labeller(type = 
                           c("negative" = "Non-Causal A-Y Pairs",
                             "maybe"    = "Candidate A-Y Pairs",
                             "positive" = "Causal A-Y Pairs"),
                          method = 
                           c("lmYA_pval" = "lmYA",
                             "lmYAX_pval"= "lmYAX",
                             "npcausal_pval"="npcausal",
                             "lmYA_ATE" = "lmYA",
                             "lmYAX_ATE"= "lmYAX",
                             "npcausal_ATE"="npcausal"))


ggplot(oracle_pvals,
       aes(x = pval)) +
  geom_histogram(binwidth = .02, alpha = .8, fill = 'orange') +
  geom_vline(aes(xintercept = .05), color = 'black', linetype = 'dashed') +
  scale_x_continuous(limits = c(-.01, 1.01)) +
  facet_grid(rows = vars(method), cols = vars(type), labeller = type_labeller)
ggsave(sprintf('%splots/pval_histogram.pdf', plot_savepath), height = 5, width = 10)


ggplot(oracle_ATEs,
       aes(x = ATE)) +
  geom_histogram(binwidth = .05, alpha = .8, fill = 'orange') +
  geom_vline(aes(xintercept = .0), color = 'black', linetype = 'dashed') +
  facet_grid(rows = vars(method), cols = vars(type), labeller = type_labeller)
ggsave(sprintf('%splots/ATE_histogram.pdf', plot_savepath), height = 5, width = 10)






