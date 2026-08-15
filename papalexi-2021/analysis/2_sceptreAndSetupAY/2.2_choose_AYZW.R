# -------------------------------------------------------------------------------- #
#                   Choose AY tests (and NCE/NCOs)                                 #
# Requires: prev saved normalized gene expression (HDF5)                           #
#           prev saved chromosome information                                      #
# Ouputs: (nothing) but saves                                                      #
#         CBGENE_AYZW in the rds file                                              #
#                         "<save_dir>/cbgenes/<AYZW_setting_name>/CBGENE_AYZW.rds" #
# -------------------------------------------------------------------------------- #
args = commandArgs(trailingOnly = TRUE)
# args = c('macbook', 'A1')
# args = c('laptop', 'C1')


suppressPackageStartupMessages(require(assertthat)) # for some assert statements
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(ggplot2))    # plotting
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(dplyr))

theme_set(theme_cowplot() +
            theme(plot.title = element_text(hjust = .5),
                  plot.subtitle = element_text(hjust = .5)))

assertthat::assert_that(length(args) > 0, msg="must give arg for specifying device eg 'Rscript <filename>.R ubergenno'")


DEVICE = args[1]
source('../PATHS.R') # load in data_dir and save_dir and CODE_DIR, depending on DEVICE value

assertthat::assert_that(length(args) > 1, msg="must give arg for specifying chosen AYZW name 'Rscript <filename>.R ubergenno C'")
AYZW_setting_name = (args[2])





# AY test settings in separate AY_settings.R file and load in `settings` object
source(sprintf('%s/AY/AY_settings.R', save_dir))
setting = settings[[AYZW_setting_name]]; rm(settings)
assertthat::assert_that(!is.null(setting), msg="bad setting name (e.g. A, A1, B, ...). make sure included in AY_settings.R")
set.seed(setting$seed)
# =================================================================================================#
# =================== START =======================================================================
# =================================================================================================#
print(sprintf("[%s] START: Choose AY and ZW", Sys.time()))


# =================== Set up saving dir + save setting ======================================#
print(sprintf("[%s]        - Set up saving dir + save setting", Sys.time()))
dir.create(sprintf('%s/AY/', save_dir), showWarnings = FALSE)
dir.create(sprintf('%s/AY/%s', save_dir, AYZW_setting_name), showWarnings = FALSE)

# also save settings in folder to make sure (in case settings .r file changes later)
capture.output(print(setting), file = sprintf('%s/AY/%s/AYZW_setting.txt', save_dir, AYZW_setting_name))
saveRDS(setting,
        sprintf('%s/AY/%s/AYZW_setting.rds', save_dir, AYZW_setting_name))

# =================================================================================================#
# =================== Load ========================================================================
# =================================================================================================#
print(sprintf("[%s]        - Load", Sys.time()))

# =================== * SCEPTRE saves ========================================================================
SCEPTRE_savepath = sprintf('%s/sceptre/', save_dir)
sceptre_results_power       = readRDS(sprintf('%s/results_run_power_check.rds',        SCEPTRE_savepath))
sceptre_results_discovery   = readRDS(sprintf('%s/results_run_discovery_analysis.rds', SCEPTRE_savepath))
sceptre_results_calibration = readRDS(sprintf('%s/results_run_calibration_check.rds',  SCEPTRE_savepath))

# =================== * count and grna data ========================================================================
gene_odm <- ondisc::read_odm(odm_fp      = paste0(data_dir, "/papalexi-2021/processed/gene/expression_matrix.odm"),
                             metadata_fp = paste0(data_dir, "/papalexi-2021/processed/gene/metadata.rds"))
grna_odm <- ondisc::read_odm(odm_fp      = paste0(data_dir, "/papalexi-2021/processed/grna_assignment/assignment_matrix.odm"),
                             metadata_fp = paste0(data_dir, "/papalexi-2021/processed/grna_assignment/metadata.rds"))

# =================== * chromosome ========================================================================
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
# gene_TF_Target = read.csv(sprintf('%s/important_genes_name_TF_Target.csv', save_dir)) # top important genes noting the TF and gRNA targets

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


# gene_metainfo = merge(gene_chr,
#                       gene_TF_Target |> select(importance_rank, gene, TF, grna_target),
#                       by.x = c('wikigene_name', 'importance_rank'),
#                       by.y = c('gene',          'importance_rank'))

# TF already removed prev, don't need to add this info (loaded in gene_chr and grna_chr start w/ TF removed already)
gene_metainfo = gene_chr
# add n_nonzero_counts to gene info (#cells ith non-zero counts) to threshold later
gene_metainfo$n_nonzero_cell = rowSums(as.matrix(gene_odm[[gene_chr[, 'wikigene_name'], ]] > 0))
# ggplot(gene_metainfo, aes(x=n_nonzero_cell)) + geom_histogram(binwidth = 100)

# =================== Make Barplots of #genes per chromosome info==================================#
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




# =================================================================================================#
# =================== Choose AY ===================================================================
# =================================================================================================#
print(sprintf("[%s]    - Choose AY ", Sys.time()))


AY = NULL

# I am not sure this cleaned dataset has chromosome info for NT perturbations
# since we don't really care about chromosomes anyway, I will just ignore for now. 
NT_A = grna_chr |> filter(target == 'non-targeting') |>
  filter(# (!is.na(target_chromosome_name)) & # ignore chr for now, bc some chr info not found
    (n_nonzero >= setting$PERTURBATION_N_NONZERO_CELLS_MIN)) |>
  select(A = grna, A_chr = target_chromosome_name)

Targ_A = grna_chr |> filter(target != 'non-targeting') |>
  filter((!is.na(target_chromosome_name)) & 
           (n_nonzero >= setting$PERTURBATION_N_NONZERO_CELLS_MIN)) |>
  select(A = grna, A_chr = target_chromosome_name)

all_Y = gene_metainfo |> filter(!is.na(chromosome_name)) |> # keep chr here for now... can remove if chr doesn't matter anyway
  filter(importance_rank <= setting$MAX_Y_IMPORTANCE & 
           n_nonzero_cell >= setting$GENE_N_NONZERO_CELLS_MIN ) |>
  select(Y = wikigene_name, Y_chr = chromosome_name)

# =================== * Choose Null (Negative) AY ===================================================
if(!is.na(setting$NUM_NULL) & setting$NUM_NULL > 0) {
  # filter A and Y for more QC
  possible_sceptre_calibration = sceptre_results_calibration |> filter(grna_id%in%NT_A$A & response_id%in%all_Y$Y)
  
  # choose tests
  if(nrow(possible_sceptre_calibration) < setting$NUM_NULL) {
    print('  When choosing AY Nulls, not enough SCEPTRE calibration tests that 
             passed additional QC (MAX_Y_IMPORTANCE, GENE_N_NONZERO_CELLS_MIN,is.na(chromosome_name).
             Try including more SCEPTRE tests (in 2.1...., change MAX_NUM_DISCOVERY_PAIRS) or lowering requirements. ')
    
    AY_null = possible_sceptre_calibration
  } else {
    AY_null = possible_sceptre_calibration[sample(1:nrow(possible_sceptre_calibration), size=setting$NUM_NULL, replace=FALSE), ] 
  }
  
  AY_null = AY_null |> 
    dplyr::mutate(A = grna_id, Y = response_id) |> 
    dplyr::select(A, Y)
  
  # add chr info
  AY_null = merge(merge(AY_null, NT_A, by = 'A', all.x = TRUE, all.y = FALSE), 
                  all_Y, by = 'Y', all.x = TRUE, all.y = FALSE) |> 
    dplyr::mutate(type = 'negative') |>
    dplyr::select(type, A, A_chr, Y, Y_chr)
  
  # add to overall AY df
  AY = rbind(AY, AY_null)
  rm(possible_sceptre_calibration, AY_null)
}

# =================== * Choose Unknown (Maybe) AY ===================================================
if(!is.na(setting$NUM_MAYBE) & setting$NUM_MAYBE > 0) {
  # filter A and Y for more QC
  possible_sceptre_discovery = sceptre_results_discovery |> filter(grna_id%in%Targ_A$A & response_id%in%all_Y$Y)
  
  # choose tests
  if(nrow(possible_sceptre_discovery) < setting$NUM_MAYBE) {
    print('  When choosing AY Maybe, not enough SCEPTRE discovery tests that 
             passed additional QC (MAX_Y_IMPORTANCE, GENE_N_NONZERO_CELLS_MIN,is.na(chromosome_name).
             Try including more SCEPTRE tests (in 2.1...., change MAX_NUM_DISCOVERY_PAIRS) or lowering requirements. ')
    
    AY_maybe = possible_sceptre_discovery
  } else {
    AY_maybe = possible_sceptre_discovery[sample(1:nrow(possible_sceptre_discovery), size=setting$NUM_MAYBE, replace=FALSE), ] 
  }
  
  AY_maybe = AY_maybe |> 
    dplyr::mutate(A = grna_id, Y = response_id) |> 
    dplyr::select(A, Y)
  
  # add chr info
  AY_maybe = merge(merge(AY_maybe, NT_A, by = 'A', all.x = TRUE, all.y = FALSE), 
                  all_Y, by = 'Y', all.x = TRUE, all.y = FALSE) |> 
    dplyr::mutate(type = 'maybe') |>
    dplyr::select(type, A, A_chr, Y, Y_chr)
  
  # add to overall AY df
  AY = rbind(AY, AY_maybe)
  rm(possible_sceptre_discovery, AY_maybe)
}




# =================== * Choose Alternative (Positive) AY ============================================
print(sprintf("[%s]        - Positive", Sys.time()))


# lower requirements for positive tests, but still need to require that the response is in the top XX
max_imp_gene = nrow(gene_chr) # only normalized these top XX important genes, so we can only use these
min_n_nonzero = .5*setting$GENE_N_NONZERO_CELLS_MIN # allow for a smaller threshold for these positive tests

# 101 perturbations, but some AY positive might be filtered out bc: not in top XX important
AY_pos = grna_chr |> 
  dplyr::filter(target != 'non-targeting') |> 
  dplyr::mutate(type = 'positive') |> 
  dplyr::select(type, A = grna, Y = target, 
                A_chr = target_chromosome_name, 
                Y_chr = target_chromosome_name)


# print(dim(AY_pos))
# print(head(AY_pos))

Y_pos = gene_metainfo |> 
  filter(importance_rank <= max_imp_gene &     # needs to be in (larger) top XXX s.t. this gene is normalized (this line is useless)
           n_nonzero_cell >= min_n_nonzero) |> # min # cells w nonzero counts for this gene 
  pull(wikigene_name)

AY_pos = AY_pos |> filter(Y %in% Y_pos) # only allow 'important' and abundant Y



AY_pos = merge(
  sceptre_results_power |> 
          dplyr::mutate(A = grna_id, Y = response_id) |> 
          dplyr::select(A, Y),
  AY_pos, by = c('A', 'Y'))

# print(dim(AY_pos))
# print(head(AY_pos))


if(!is.na(setting$NUM_ALT)) {
  AY_pos = AY_pos[sample(1:nrow(AY_pos), size=setting$NUM_ALT), ] # randomly sample NUM_ALT pos tests
}

AY = rbind(AY, AY_pos)

rm(AY_pos, Y_pos, max_imp_gene, min_n_nonzero)


AY = data.frame(AY)



# =================== Choose NCE/NCO (Z/W) single genes ===========================================
# =================================================================================================#

print(sprintf("[%s]    - Choose some ZW (NCE/NCO)", Sys.time()))

# create list. indexed by
# $A_name
#   $Y_name
#     $ZW_idx
#       $'AY_idx'  = integer
#       $'Z_names' = vector of str
#       $'W_names' = vector of str
#       $'A_chr'   = str (1:22, A,Y) 
#       $'Y_chr'   = str (1:22, A,Y)
#       $'Z_chrs'   = vector of str (1:22, A,Y) 
#       $'W_chrs'   = vector of str (1:22, A,Y) 





AYZW = list()

for(i in 1:nrow(AY)) {
  
  A_name = AY[i, 'A']
  Y_name = AY[i, 'Y']
  A_chr  = AY[i, 'A_chr']
  Y_chr  = AY[i, 'Y_chr']
  if(i %% 20 == 0) {
    print(sprintf("[%s]           i=%d: A=%s Y=%s", Sys.time(), i, A_name, Y_name))
  }
  
  
  # # add to AYZW list structure if not present, but we shouldn't need to do this?? it should automatically add?
  # if(!(A_name %in% names(AYZW))) {
  #   AYZW[[A_name]] = list()
  # }
  # if(!(Y_name %in% names(AYZW[[A_name]]))) {
  #   AYZW[[A_name]][[Y_name]] = list()
  # }
  
  for(j in 1:setting$NUM_NCENCO_per_AY) {
    ZW_chr = setdiff(chr_factors, 
                     c( if(is.na(A_chr)) NULL else A_chr, 
                        if(is.na(Y_chr)) NULL else Y_chr)) # chr for ZW (no AY chrs)
    Z_chrs = ZW_chr[runif(length(ZW_chr)) < .5]    # chr for Z  (1/2 each)
    W_chrs = setdiff(ZW_chr, Z_chrs)               # chr for W  (take remaining)
    
    
    # allow setting$NUM_NCE/O=NA option
    if(!is.null(setting$NUM_NCENCO_pairs) && !is.na(setting$NUM_NCENCO_pairs)) {  
      # take equal number of ZW (no longer used...keep in case we do again in future)
      # will supersede individually set NUM_NCE and NUM_NCO
      Z = all_Y |> filter(Y_chr %in% Z_chrs) |> slice_sample(n = setting$NUM_NCENCO_pairs)
      W = all_Y |> filter(Y_chr %in% W_chrs) |> slice_sample(n = setting$NUM_NCENCO_pairs)
    } else { #NUM_NCENCO_pairs not specified, try NUM_NCE and NUM_NCO
      if(is.na(setting$NUM_NCE)) { # no sampling, take all
        Z = all_Y |> filter(Y_chr %in% Z_chrs)
      } else {
        Z = all_Y |> filter(Y_chr %in% Z_chrs) |> slice_sample(n = setting$NUM_NCE)
      }
      
      if(is.na(setting$NUM_NCO)) {
        W = all_Y |> filter(Y_chr %in% W_chrs)
      } else {
        W = all_Y |> filter(Y_chr %in% W_chrs) |> slice_sample(n = setting$NUM_NCO)
      }
    }
    
    
    
    # A_name$Y_name$ZW_idx$...
    AYZW[[A_name]][[Y_name]][[j]] = list()
    AYZW[[A_name]][[Y_name]][[j]][['AY_idx']] = i
    
    AYZW[[A_name]][[Y_name]][[j]][['A_chr']] = A_chr
    AYZW[[A_name]][[Y_name]][[j]][['Y_chr']] = Y_chr
    
    AYZW[[A_name]][[Y_name]][[j]][['Z_names']] = Z$Y
    AYZW[[A_name]][[Y_name]][[j]][['W_names']] = W$Y
    
    AYZW[[A_name]][[Y_name]][[j]][['Z_chrs']]  = Z$Y_chr
    AYZW[[A_name]][[Y_name]][[j]][['W_chrs']]  = W$Y_chr
  }
}


# AYZW
# object.size(AYZW)
# format(object.size(AYZW), units = 'MB')
# names(AYZW)
# names(AYZW$PDL1g1)
# names(AYZW$ATF2g4)
# names(AYZW$ATF2g4$SERPINE2)
# 
# AYZW$ATF2g4$SERPINE2[[1]]
# length(AYZW$ATF2g4$SERPINE2)


# =================================================================================================#
# =================== Save  =======================================================================
# =================================================================================================#
print(sprintf("[%s]    - Save AY and chosen NCE/NCO", Sys.time()))
print(sprintf("[%s]        - nrow(AY) = %s", Sys.time(), nrow(AY)))
write.csv(AY, sprintf('%s/AY/%s/AY.csv', save_dir, AYZW_setting_name), row.names = FALSE)
saveRDS(AYZW, sprintf('%s/AY/%s/AYZW.rds', save_dir, AYZW_setting_name))


# =================== END =========================================================================
print(sprintf("[%s] END", Sys.time()))




# === > TRASH === === === -----
if(F) {
  # =================================================================================================#
  # =================== Choose AY (OLD) =============================================================
  # =================================================================================================#
  print(sprintf("[%s]    - Choose AY ", Sys.time()))
  
  
  AY = NULL
  
  
  
  
  # I am not sure this cleaned dataset has chromosome info for NT perturbations
  # since we don't really care about chromosomes anyway, I will just ignore for now. 
  NT_A = grna_chr |> filter(target == 'non-targeting') |>
    filter(# (!is.na(target_chromosome_name)) & # ignore chr for now, bc some chr info not found
      (n_nonzero >= setting$PERTURBATION_N_NONZERO_CELLS_MIN)) |>
    select(A = grna, A_chr = target_chromosome_name)
  
  Targ_A = grna_chr |> filter(target != 'non-targeting') |>
    filter((!is.na(target_chromosome_name)) & 
             (n_nonzero >= setting$PERTURBATION_N_NONZERO_CELLS_MIN)) |>
    select(A = grna, A_chr = target_chromosome_name)
  
  all_Y = gene_metainfo |> filter(!is.na(chromosome_name)) |> # keep chr here for now... can remove if chr doesn't matter anyway
    filter(importance_rank <= setting$MAX_Y_IMPORTANCE & 
             n_nonzero_cell >= setting$GENE_N_NONZERO_CELLS_MIN ) |>
    select(Y = wikigene_name, Y_chr = chromosome_name)
  
  
  # =================== * Choose Null (Negative) AY ===================================================
  if(!is.na(setting$NUM_NULL) & setting$NUM_NULL > 0) {
    # total number of possible tests is nrow(NT_A) x nrow(all_Y)
    AY_null_idx = sample.int(n = nrow(NT_A)*nrow(all_Y), size = setting$NUM_NULL, replace = FALSE)
    AY_null = data.frame(A = rep(NA, times = setting$NUM_NULL), 
                         Y = rep(NA, times = setting$NUM_NULL))
    for(i in 1:setting$NUM_NULL) {
      AY_null[i, 'A'] =  NT_A[ceiling(AY_null_idx[i] /  nrow(all_Y)), 'A']
      AY_null[i, 'Y'] = all_Y[1 +    (AY_null_idx[i] %% nrow(all_Y)), 'Y']
    }
    
    # AY_null$A |> table() |> barplot() # should be around even
    # AY_null$Y |> table() |> barplot()
    # head(AY_null); dim(AY_null)
    
    
    AY_null = merge(merge(AY_null, NT_A, all.x = TRUE, all.y = FALSE), 
                    all_Y, all.x = TRUE, all.y = FALSE) |> 
      dplyr::mutate(type = 'negative') |>
      dplyr::select(type, A, A_chr, Y, Y_chr)
    
    
    AY = rbind(AY, AY_null)
    
    rm(AY_null_idx, i, AY_null)
  }
  
  
  # =================== * Choose Unknown (maybe) AY ===================================================
  if(!is.na(setting$NUM_MAYBE) & setting$NUM_MAYBE > 0) {
    
    # total number of possible tests is nrow(Targ_A) x nrow(all_Y)
    AY_maybe_idx = sample.int(n = nrow(Targ_A)*nrow(all_Y), size = setting$NUM_MAYBE, replace = FALSE)
    AY_maybe = data.frame(A = rep(NA, times = setting$NUM_MAYBE), 
                          Y = rep(NA, times = setting$NUM_MAYBE))
    for(i in 1:setting$NUM_MAYBE) {
      AY_maybe[i, 'A'] = Targ_A[ceiling(AY_maybe_idx[i] /  nrow(all_Y)), 'A']
      AY_maybe[i, 'Y'] =  all_Y[1 +    (AY_maybe_idx[i] %% nrow(all_Y)), 'Y']
    }
    
    # AY_maybe$A |> table() |> barplot() # should be around even
    # AY_maybe$Y |> table() |> barplot()
    # head(AY_maybe); dim(AY_maybe)
    
    
    AY_maybe = merge(merge(AY_maybe, Targ_A, all.x = TRUE, all.y = FALSE), 
                     all_Y, all.x = TRUE, all.y = FALSE) |> 
      dplyr::mutate(type = 'maybe') |>
      dplyr::select(type, A, A_chr, Y, Y_chr)
    
    AY = rbind(AY, AY_maybe)
    rm(AY_maybe_idx, i, AY_maybe)
  }
  
  
  # =================== * Choose Alternative (Positive) AY ============================================
  print(sprintf("[%s]        - Positive", Sys.time()))
  
  # 101 perturbations, but some AY positive might be filtered out bc: not in top XX important
  AY_pos = grna_chr |> 
    dplyr::filter(target != 'non-targeting') |> 
    dplyr::mutate(type = 'positive') |> 
    dplyr::select(type, A = grna, Y = target, 
                  A_chr = target_chromosome_name, 
                  Y_chr = target_chromosome_name)
  
  max_imp_gene = nrow(gene_chr) # only normalized these top XX important genes, so we can only use these
  min_n_nonzero = .5*setting$GENE_N_NONZERO_CELLS_MIN # allow for a smaller threshold for these positive tests
  
  # print(dim(AY_pos))
  # print(head(AY_pos))
  
  Y_pos = gene_metainfo |> 
    filter(importance_rank <= max_imp_gene &     # needs to be in (larger) top XXX s.t. this gene is normalized (this line is useless)
             n_nonzero_cell >= min_n_nonzero) |> # min # cells w nonzero counts for this gene 
    pull(wikigene_name)
  
  AY_pos = AY_pos |> filter(Y %in% Y_pos) # only allow 'important' and abundant Y
  
  
  # print(dim(AY_pos))
  # print(head(AY_pos))
  
  
  if(!is.na(setting$NUM_ALT)) {
    AY_pos = AY_pos[sample(1:nrow(AY_pos), size=setting$NUM_ALT), ] # randomly sample NUM_ALT pos tests
  }
  
  AY = rbind(AY, AY_pos)
  
  rm(AY_pos, Y_pos, max_imp_gene, min_n_nonzero)
  
  # # c('PDL1', 'PDL2', 'CD86') %in% (gene_odm |> ondisc::get_feature_covariates() |> rownames())
  # # could add pos protein effects... but new ones don't have genes measured  PDL1, PDL2
  # AY_pos = rbind(AY_pos,
  #                grna_chr |> filter(!is.na(known_protein_effect)) |> 
  #                  mutate(type = 'positive') |>
  #                  select(type, A = grna, Y = known_protein_effect)) |> distinct()
  
  
  # =================================================================================================#
}

