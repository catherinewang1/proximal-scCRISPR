# -------------------------------------------------------------------------------- #
#                   Do CB using genes as NCE/NCO                                   #
# 3.1 - add chromosome number information                                          #
# 3.2 - choose A,Y,Z,W combinations                                                #
# 3.2 - CB Effect Estimate                                                         #
# Requires: prev saved normalized gene expression (HDF5)                           #
#           prev saved chromosome information                                      #
# Ouputs: (nothing) but saves                                                      #
#         CBGENE_AYZW in the rds file                                              #
#                         "<save_dir>/cbgenes/<AYZW_setting_name>/CBGENE_AYZW.rds" #
# -------------------------------------------------------------------------------- #
args = commandArgs(trailingOnly = TRUE)
# args = c('laptop', 'A1')

require(assertthat) # for some assert statements
library(dplyr)
library(tibble)
library(ggplot2)    # plotting
library(cowplot)

theme_set(theme_cowplot() +
          theme(plot.title = element_text(hjust = .5),
                plot.subtitle = element_text(hjust = .5)))

assertthat::assert_that(length(args) > 0, msg="must give arg for specifying device eg 'Rscript <filename>.R ubergenno'")


DEVICE = args[1]
source('../PATHS.R') # load in data_dir and save_dir and CODE_DIR, depending on DEVICE value
# DEVICE = args[1]
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
# assertthat::assert_that(!is.null(data_dir), msg='first arg must be: laptop, desktop, or ubergenno')


assertthat::assert_that(length(args) > 1, msg="must give arg for specifying chosen AYZW name 'Rscript <filename>.R ubergenno C'")
AYZW_setting_name = (args[2])

# SETTINGS HERE
# setting = list(seed = 942346,
#                NUM_A             = NA, # NA if all As
#                NUM_Y_PER_A_NEG   = 5,
#                NUM_Y_PER_A_MAYBE = 0,
#                MAX_Y_IMPORTANCE  = 2000, # limit how 'unimportant' a gene can be
#                # NUM_NCENCO_pairs  = 100,  # number of NCE/NCO pairs (dimU, length of ZWs)
#                NUM_NCE           = NA,  # number of NCE per AY test (prev #NCE/NCO equal) (NA=all avail)
#                NUM_NCO           = NA,  # number of NCO per AY test 
#                NUM_NCENCO_per_AY = 5,   # number of NCE/NCO sets per AY test
#                NUM_AY_POS        = NA    # number of known causal/positive AY tests
# )
# set.seed(setting$seed)

# change setting list to test for ALL AY pairs. Too many to test for all possible,
# so choose all A's, all A's on the same chromosome, and 
setting = list(seed = 942346,
               NUM_A             = 5, # NA if all As (probably SHOULD)
               NUM_Y_PER_A_NEG   = 5,  # NA if all (probably should NOT)
               NUM_Y_PER_A_MAYBE = 5, # NA if all on same chromosome (probably SHOULD)
               MAX_Y_IMPORTANCE  = 2000, # limit how 'unimportant' a response gene can be
               # NUM_NCENCO_pairs  = 100,  # number of NCE/NCO pairs (dimU, length of ZWs)
               # NUM_NCE           = NA,  # number of NCE per AY test (prev #NCE/NCO equal) (NA=all avail)
               # NUM_NCO           = NA,  # number of NCO per AY test 
               # NUM_NCENCO_per_AY = 1,   # number of NCE/NCO sets per AY test
               NUM_AY_POS        = 5    # number of known causal/positive AY tests, NA if all As (probably SHOULD)
)

set.seed(setting$seed)

# =================== Start ====================================================
print(sprintf("[%s] START: CB Choose AYZW", Sys.time()))


# =================== Set up saving dir + save setting ======================================
print(sprintf("[%s]        - Set up saving dir + save setting", Sys.time()))
dir.create(sprintf('%s/spca/cbgenes/', save_dir), showWarnings = FALSE)
dir.create(sprintf('%s/spca/cbgenes/%s', save_dir, AYZW_setting_name), showWarnings = FALSE)

capture.output(print(setting), file = sprintf('%s/spca/cbgenes/%s/AYZW_setting.txt', save_dir, AYZW_setting_name))
saveRDS(setting,
        sprintf('%s/spca/cbgenes/%s/AYZW_setting.rds', save_dir, AYZW_setting_name))

# =================== load ======================================
print(sprintf("[%s]        - load", Sys.time()))

# gene normalized info
gene_odm <- ondisc::read_odm(odm_fp      = paste0(data_dir, "/papalexi-2021/processed/gene/expression_matrix.odm"),
                             metadata_fp = paste0(data_dir, "/papalexi-2021/processed/gene/metadata.rds"))
grna_odm <- ondisc::read_odm(odm_fp      = paste0(data_dir, "/papalexi-2021/processed/grna_assignment/assignment_matrix.odm"),
                             metadata_fp = paste0(data_dir, "/papalexi-2021/processed/grna_assignment/metadata.rds"))


chr_factors = c(as.character(1:22), 'A', 'Y', NA)
# load gene chr info
gene_chr = read.csv(sprintf('%s/chromosome/gene_chromosome.csv', save_dir)) |> # 4000 x 7
           mutate(chromosome_name = 
                      factor(chromosome_name, ordered = TRUE,
                             levels = chr_factors))
# load grna chr info
grna_chr = read.csv(sprintf('%s/chromosome/grna_chromosome.csv', save_dir)) |> # 110 x 5
           mutate(target_chromosome_name = 
                      factor(target_chromosome_name, ordered = TRUE,
                             levels = chr_factors))

# load gene info (top XX genes, removing TF and Targets)
gene_TF_Target = read.csv(sprintf('%s/important_genes_name_TF_Target.csv', save_dir)) # top important genes noting the TF and gRNA targets

# # check that the gene_TF_Target info and gene_chr match up
# #   - importance_rank = importance_rank
# #   - gene = wikigene_name
# #   - idx = gene_idx
# a = gene_TF_Target |> select(importance_rank, gene, idx) 
# b = gene_chr |> select(importance_rank, wikigene_name, chromosome_name, gene_idx) |> arrange(importance_rank)
# 
# all(a$importance_rank == b$importance_rank)
# all(a$gene == b$wikigene_name)
# all(a$idx == b$gene_idx)
# 
# 
# c = merge(a, b, by = 'importance_rank')
# head(c)
# rm(a,b,c)
# gene_TF_Target |> select(importance_rank, gene, TF, grna_target) |> head()

# new gene meta information with: importance_rank, gene, gene_idx, chromosome_name, TF, grna_target 
# gene_metainfo = merge(gene_TF_Target |> rename(gene_idx = idx),
#                       gene_chr |> select(importance_rank, chromosome_name),
#                       by = 'importance_rank')


gene_metainfo = merge(gene_chr,
                      gene_TF_Target |> select(importance_rank, gene, TF, grna_target),
                      by.x = c('wikigene_name', 'importance_rank'),
                      by.y = c('gene',          'importance_rank'))


# =================== Make Barplots of #genes per chromosome info ======================================
print(sprintf("[%s]    - Make Barplots of #genes per chromosome info", Sys.time()))

# Barplot: genes per chromosome
# gene_chr |> group_by(chromosome_name) |> summarize(count = n())
ggplot(gene_chr) +
  geom_bar(aes(x = chromosome_name), fill = 'aquamarine3') +
  labs(title = 'Number of Genes per Chromosome',
       subtitle = sprintf('(of Top %d)', nrow(gene_chr)),
       x = 'Chromosome', y = 'Count') +
  scale_y_continuous(expand = c(0, 0, 0, 0)) +
  scale_x_discrete(breaks = c(chr_factors), drop = F)

ggsave(filename = sprintf('%s/chromosome/barplot_genes_per_chromosome.svg', save_dir),
       height = 4, width = 6)


# Barplot: gRNA target genes per chromosome
ggplot(grna_chr |> filter(target != 'non-targeting')) +
  geom_bar(aes(x = target_chromosome_name), fill = 'aquamarine3') +
  labs(title = 'Number of gRNA Target Genes per Chromosome',
       subtitle = sprintf('(Targeting gRNA)'),
       x = 'Chromosome', y = 'Count') +
  scale_y_continuous(expand = c(0, 0, 0, 0)) +
  scale_x_discrete(breaks = c(chr_factors), drop = F)

ggsave(filename = sprintf('%s/chromosome/barplot_grnas_per_chromosome.svg', save_dir),
       height = 4, width = 6)


# =================== choose some AY  ======================================
print(sprintf("[%s]    - Choose some AY ", Sys.time()))

# =================== choose some AY (Negative/Maybe)  ======================================
print(sprintf("[%s]        - Negative/Maybe", Sys.time()))

# num grna small enough here papalexi-2021 to do full comb? (cannot/shouldn't do in gasperini-2019)
# # 101 x 3731
# AY = merge(grna_chr |> filter(target != 'non-targeting') |>
#                        select(A = grna, A_chr = target_chromosome_name),
#            gene_chr |> filter(!is.na(chromosome_name)) |> 
#                        select(Y = wikigene_name, Y_chr = chromosome_name))

# messy but whatever
AY = NULL
all_A = grna_chr |> filter(target != 'non-targeting') |>
           filter((!is.na(target_chromosome_name)) & 
                    (n_nonzero >= 200)) |>
           select(A = grna, A_chr = target_chromosome_name)
all_Y = gene_metainfo |> filter(!is.na(chromosome_name)) |> 
           filter(importance_rank <= setting$MAX_Y_IMPORTANCE) |>
           select(Y = wikigene_name, Y_chr = chromosome_name)

# define all_A_idx, idx for A to get Y values for
if(is.na(setting$NUM_A)) {
  all_A_idx = 1:nrow(all_A)
} else if(setting$NUM_A > nrow(all_A)) {
  print(sprintf("[%s]            - bad setting$NUM_A value (> possible As)", Sys.time()))
} else {
  all_A_idx = sample(1:nrow(all_A), setting$NUM_A, replace = FALSE)
}


for(i in all_A_idx) {
  A_chr_i   = all_A[i, 'A_chr'] # A target chr
  maybe_idx = which(all_Y$Y_chr == A_chr_i)
  neg_idx   = which(all_Y$Y_chr != A_chr_i)
  
  # if NA, choose ALL 
  if(is.na(setting$NUM_Y_PER_A_MAYBE)) {
    maybe_Y_i = all_Y[maybe_idx, ]
  } else {
    maybe_Y_i = all_Y[sample(x = maybe_idx,  # rand choose maybe Y's
                             size = min(length(maybe_idx), setting$NUM_Y_PER_A_MAYBE), 
                             replace = FALSE), ]
  }
  
  # if NA, choose ALL
  if(is.na(setting$NUM_Y_PER_A_NEG)) {
    neg_Y_i   = all_Y[neg_idx, ]
  } else {
    neg_Y_i   = all_Y[sample(x = neg_idx,  # rand choose negative Y's
                             size = min(length(neg_idx), setting$NUM_Y_PER_A_NEG), 
                             replace = FALSE), ]
  }
  
  AY = rbind(AY,
             dplyr::bind_cols(all_A[i, ], 
                              rbind(maybe_Y_i, neg_Y_i)))
}
AY = AY |> mutate(type = case_when(A_chr != Y_chr ~ 'negative', 
                                             TRUE ~ 'maybe')) |>
     select(type, A, A_chr, Y, Y_chr)


# =================== choose some AY (Positive)  ======================================
print(sprintf("[%s]        - Positive", Sys.time()))

# 101 perturbations, but some AY positive might be filtered out bc: not in top XX important
AY_pos = grna_chr |> filter(target != 'non-targeting') |> 
         mutate(type = 'positive') |> 
         select(type, A = grna, Y = target, 
                A_chr = target_chromosome_name, Y_chr = target_chromosome_name)

AY_pos = AY_pos |> filter(Y %in% all_Y$Y) # only allow 'important' and abundant Y
if(!is.na(setting$NUM_AY_POS)) {
  AY_pos = AY_pos[sample(1:setting$NUM_AY_POS), ] # randomly sample NUM_AY_POS pos tests
}


# # c('PDL1', 'PDL2', 'CD86') %in% (gene_odm |> ondisc::get_feature_covariates() |> rownames())
# # could add pos protein effects... but new ones don't have genes measured  PDL1, PDL2
# AY_pos = rbind(AY_pos,
#                grna_chr |> filter(!is.na(known_protein_effect)) |> 
#                  mutate(type = 'positive') |>
#                  select(type, A = grna, Y = known_protein_effect)) |> distinct()


# =================== Combine chosen AY   ======================================
print(sprintf("[%s]        - Combine", Sys.time()))

AY = rbind(AY, AY_pos)
# AY |> group_by(type) |> summarize(count = n())



# =================== DO NOT NEED TO CHOOSE NCE/NCO WHEN USING SPARSE PCA ===========
# # =================== choose some ZW (NCE/NCO) ======================================
# print(sprintf("[%s]    - Choose some ZW (NCE/NCO)", Sys.time()))
# 
# # create list. indexed by
# # $A_name
# #   $Y_name
# #     $ZW_idx
# #       $'AY_idx'  = integer
# #       $'Z_names' = vector of str
# #       $'W_names' = vector of str
# #       $'A_chr'   = str (1:22, A,Y) 
# #       $'Y_chr'   = str (1:22, A,Y)
# #       $'Z_chrs'   = vector of str (1:22, A,Y) 
# #       $'W_chrs'   = vector of str (1:22, A,Y) 
# 
# 
# 
# AYZW = list()
# 
# for(i in 1:nrow(AY)) {
#   
#   A_name = AY[i, 'A']
#   Y_name = AY[i, 'Y']
#   A_chr  = AY[i, 'A_chr']
#   Y_chr  = AY[i, 'Y_chr']
#   if(i %% 20 == 0) {
#     print(sprintf("[%s]           i=%d: A=%s Y=%s", Sys.time(), i, A_name, Y_name))
#   }
#   
#   
#   for(j in 1:setting$NUM_NCENCO_per_AY) {
#     ZW_chr = setdiff(chr_factors, c(A_chr, Y_chr)) # chr for ZW (no AY chrs)
#     Z_chrs = ZW_chr[runif(length(ZW_chr)) < .5]    # chr for Z  (1/2 each)
#     W_chrs = setdiff(ZW_chr, Z_chrs)               # chr for W  (take remaining)
#     
# 
#     # allow setting$NUM_NCE/O=NA option
#     if(is.na(setting$NUM_NCE)) { # no sampling, take all
#       Z = all_Y |> filter(Y_chr %in% Z_chrs)
#     } else {
#       Z = all_Y |> filter(Y_chr %in% Z_chrs) |> slice_sample(n = setting$NUM_NCE)
#     }
#     
#     if(is.na(setting$NUM_NCO)) {
#       W = all_Y |> filter(Y_chr %in% W_chrs)
#     } else {
#       W = all_Y |> filter(Y_chr %in% W_chrs) |> slice_sample(n = setting$NUM_NCO)
#     }
#     
#     
# 
#     # A_name$Y_name$ZW_idx$...
#     AYZW[[A_name]][[Y_name]][[j]] = list()
#     AYZW[[A_name]][[Y_name]][[j]][['AY_idx']] = i
#     
#     AYZW[[A_name]][[Y_name]][[j]][['A_chr']] = A_chr
#     AYZW[[A_name]][[Y_name]][[j]][['Y_chr']] = Y_chr
#     
#     AYZW[[A_name]][[Y_name]][[j]][['Z_names']] = Z$Y
#     AYZW[[A_name]][[Y_name]][[j]][['W_names']] = W$Y
#     
#     AYZW[[A_name]][[Y_name]][[j]][['Z_chrs']]  = Z$Y_chr
#     AYZW[[A_name]][[Y_name]][[j]][['W_chrs']]  = W$Y_chr
#   }
# }
# 
# 
# # AYZW
# # object.size(AYZW)
# # names(AYZW)
# # names(AYZW$PDL1g1)
# # names(AYZW$ATF2g4)
# # names(AYZW$ATF2g4$SERPINE2)
# # 
# # AYZW$ATF2g4$SERPINE2[[1]]
# # length(AYZW$ATF2g4$SERPINE2)


# =================== Saving AYZW names ======================================
print(sprintf("[%s]    - Saving AYZW names", Sys.time()))
print(sprintf("[%s]        - nrow(AY) = %s", Sys.time(), nrow(AY)))
write.csv(AY, sprintf('%s/spca/cbgenes/%s/AY.csv', save_dir, AYZW_setting_name), row.names = FALSE)
# saveRDS(AYZW, sprintf('%s/spca/cbgenes/%s/AYZW.rds', save_dir, AYZW_setting_name))


# =================== END ======================================================
print(sprintf("[%s] END", Sys.time()))

