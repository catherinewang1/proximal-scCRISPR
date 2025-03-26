# -------------------------------------------------------------------------------- #
#                   Do CB using `super`genes as NCE/NCO                            #
# 7.1 - normalize genes                                                            #
# 7.2 - construct supergenes and supergene count expression                        #
# 7.3 - choose A,Y,Z,W combinations                                                #
# 7.4 - CB Effect Estimate                                                         #
# Requires: prev saved normalized gene expression (HDF5)                           #
#           prev saved chromosome information                                      #
# Ouputs: (nothing) but saves                                                      #
#         CBGENE_AYZW in the rds file                                              #
#                      "<save_dir>/supergene/<AYZW_setting_name>/AYZW.rds"         #
# -------------------------------------------------------------------------------- #
args = commandArgs(trailingOnly = TRUE)
args = c('laptop', 'Atest')



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
setting = list(seed = 942346,
               NUM_A             = NA, # NA if all As
               NUM_Y_PER_A_NEG   = 5,
               NUM_Y_PER_A_MAYBE = 5,
               MAX_Y_IMPORTANCE  = 2000, # limit how 'unimportant' a gene can be
               # NUM_NCENCO_pairs  = 100,  # number of NCE/NCO pairs (dimU, length of ZWs)
               NUM_NCE           = 10,  # number of NCE per AY test (prev #NCE/NCO equal) (NA=all avail)
               NUM_NCO           = 10,  # number of NCO per AY test 
               NUM_NCENCO_per_AY = 1,   # number of NCE/NCO sets per AY test
               NUM_AY_POS        = NA    # number of known causal/positive AY tests
)
set.seed(setting$seed)


# =================== Start ====================================================
print(sprintf("[%s] START: CB Choose AYZW", Sys.time()))


# =================== Set up saving dir + save setting ======================================
print(sprintf("[%s]        - Set up saving dir + save setting", Sys.time()))
dir.create(sprintf('%s/supergene/', save_dir), showWarnings = FALSE)
dir.create(sprintf('%s/supergene/%s', save_dir, AYZW_setting_name), showWarnings = FALSE)

capture.output(print(setting), file = sprintf('%s/supergene/%s/AYZW_setting.txt', save_dir, AYZW_setting_name))
saveRDS(setting,
        sprintf('%s/supergene/%s/AYZW_setting.rds', save_dir, AYZW_setting_name))

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


# load grna chr info and Add TF or not TF info (for A's)
grna_chr = read.csv(sprintf('%s/chromosome/grna_chromosome.csv', save_dir)) |> # 110 x 5
  mutate(target_chromosome_name = 
           factor(target_chromosome_name, ordered = TRUE,
                  levels = chr_factors))
# TF info
tf_raw = readxl::read_xlsx(path = paste0(data_dir, "/../extra/transcriptionfactorlist.xlsx"),  # given from Kathryn over slack (suppl of a paper?)
                           sheet = 2) |> suppressMessages() # suppress messages on how they renamed columns
tf = tf_raw[-1 , c('...2', 'Is TF?')]
colnames(tf) = c('gene_name', 'TF')

grna_chr = merge(grna_chr, tf,
                 by.x = 'target',
                 by.y = 'gene_name')



# load supergene info
supergene_membership = read.csv(sprintf('%s/supergene/supergene_membership.csv', save_dir))
sgene_chr = merge(supergene_membership,
                  gene_chr |> select(wikigene_name, chromosome_name),
                  by.x = 'gene_name',
                  by.y = 'wikigene_name',
                  all.x=TRUE, all.y=FALSE)
rm(supergene_membership)
# chrs per supergene
ggplot(sgene_chr |> 
         distinct(supergene_membership, chromosome_name) |> 
         group_by(supergene_membership) |> 
         summarize(nchrs = n()),
       aes(x = supergene_membership, y = nchrs)) +
  geom_col() +
  labs(title = 'Number of Distinct Chrs for each Supergene',
       x = 'Supergene', y = 'Num Chrs')

# chrs per supergene
ggplot(sgene_chr |> 
         distinct(supergene_membership, chromosome_name) |> 
         group_by(supergene_membership) |> 
         summarize(nchrs = n()),
       aes(x = nchrs, y = after_stat(density))) +
  geom_histogram() +
  labs(title = 'Histogram of Distinct Chrs for each Supergene',
       x = 'Num Chrs', y = 'count')


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
all_Y = gene_chr |> filter(!is.na(chromosome_name)) |> 
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
  
  maybe_Y_i = all_Y[sample(x = maybe_idx,  # rand choose maybe Y's
                           size = min(length(maybe_idx), setting$NUM_Y_PER_A_MAYBE), 
                           replace = FALSE), ]
  
  neg_Y_i   = all_Y[sample(x = neg_idx,  # rand choose negative Y's
                           size = min(length(neg_idx), setting$NUM_Y_PER_A_NEG), 
                           replace = FALSE), ]
  AY = rbind(AY,
             dplyr::bind_cols(all_A[i, ], 
                              rbind(maybe_Y_i, neg_Y_i)))
}
AY = AY |> mutate(type = case_when(A_chr != Y_chr ~ 'negative', 
                                   TRUE ~ 'maybe')) |>
  select(type, A, A_chr, Y, Y_chr)


# =================== choose some AY (Positive)  ======================================
print(sprintf("[%s]        - Positive", Sys.time()))

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




# TODO: !!! This has to change, instead of individual genes, choose supergenes, and ensure no same chrs
# =================== choose some ZW (NCE/NCO) ======================================
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
  
  
  for(j in 1:setting$NUM_NCENCO_per_AY) {
    # ======== get available sgenes for NCE and NCO
    # can set a threshold (e.g. <=c% can be on a violating chrs)? (for now, strict? <= 0) (or .1 for testing)
    chr_viol_thresh = 0.1
    # Find supergenes without (A) chrs (NCOs)
    W_sgenes_nonviol = sgene_chr |> 
                      group_by(supergene_membership) |> 
                      summarize(percent_violating_chr = mean(chromosome_name %in% c(A_chr))) |>
                      filter(percent_violating_chr <= chr_viol_thresh) |>
                      arrange(percent_violating_chr, supergene_membership)
    
    # Find supergenes without (Y) chrs (NCEs)
    Z_sgenes_nonviol = sgene_chr |> 
      group_by(supergene_membership) |> 
      summarize(percent_violating_chr = mean(chromosome_name %in% c(Y_chr))) |>
      filter(percent_violating_chr <= chr_viol_thresh) |>
      arrange(percent_violating_chr, supergene_membership)
    
    # ======== Assign sgenes to NCE or NCO
    
    # ineff but ok for now, idealy set data structure w fast operations
    # number of NC(O or E) to have
    if(is.na(setting$NUM_NCE)) {
      num_nce = nrow(Z_sgenes_nonviol) # max avail NCs (probably stop before)
    } else {
      num_nce = setting$NUM_NCE        # chosen # of NCs
    }
    if(is.na(setting$NUM_NCO)) {
      num_nco = nrow(W_sgenes_nonviol) # max avail NCs (probably stop before)
    } else {
      num_nco = setting$NUM_NCO        # chosen # of NCs
    }
    
    # pool of sgenes for NCE/NCO
    W_sgenes_remaining = W_sgenes_nonviol$supergene_membership
    Z_sgenes_remaining = Z_sgenes_nonviol$supergene_membership
    
    
    # # debugging
    # num_nce = 4
    # num_nco = 4
    # W_sgenes_remaining = c(1, 2, 3, 4, 5, 6, 7, 9)
    # Z_sgenes_remaining = c(1, 3, 8, 4, 7, 6, 10)
    
    # # don't do this, this can result in suboptimal choices, start with nothing
    # #    assign sgenes only ok for either NCE or NCO
    # W_sgenes = setdiff(W_sgenes_remaining, Z_sgenes_remaining)
    # Z_sgenes = setdiff(Z_sgenes_remaining, W_sgenes_remaining)
    # if(length(W_sgenes) > num_nco) { W_sgenes = W_sgenes[1:num_nco]}
    # if(length(Z_sgenes) > num_nce) { Z_sgenes = Z_sgenes[1:num_nce]}
    # 
    # # update remaining available sgenes
    # W_sgenes_remaining = setdiff(W_sgenes_remaining, W_sgenes)
    # Z_sgenes_remaining = setdiff(Z_sgenes_remaining, Z_sgenes)
    W_sgenes = c()
    Z_sgenes = c()
    
    
    # # Stopping condition: got enough NCs OR Not enough available
    # # continuing cond: not enough NCs AND some still avail
    # # Still need NCO
    # (length(W_sgenes) < num_nco) & (length(W_sgenes_remaining) > 0)
    # # Or
    # 
    # # Still need NCO
    # (length(Z_sgenes) < num_nce) & (length(Z_sgenes_remaining) > 0)
    # 
    
    while(((length(W_sgenes) < num_nco) & (length(W_sgenes_remaining) > 0))
          |
          ((length(Z_sgenes) < num_nce) & (length(Z_sgenes_remaining) > 0))  ) {
      chosen_sgene = NA
      chosen_NCO  = NA
      # if both need more sgenes
      if((length(W_sgenes) < num_nco) & (length(W_sgenes_remaining) > 0)
         &
         (length(Z_sgenes) < num_nce) & (length(Z_sgenes_remaining) > 0)) {
        
        # give to the one with fewer NCs so far, prioritizing NCO
        if(length(W_sgenes) <= length(Z_sgenes)) {
          chosen_NCO = TRUE
        } else {
          chosen_NCO = FALSE
        }
        
      } else if((length(W_sgenes) < num_nco) & (length(W_sgenes_remaining) > 0)) {
        # only NCO can get it
        chosen_NCO = TRUE
      } else {
        # only NCE can get it
        chosen_NCO = FALSE
      }
      
      # give next best to W
      if(chosen_NCO) {
        chosen_sgene = W_sgenes_remaining[1]
        W_sgenes = c(W_sgenes, chosen_sgene)
      } else { # give next best to W
        chosen_sgene = Z_sgenes_remaining[1]
        Z_sgenes = c(Z_sgenes, chosen_sgene)
      }
      
      # remove from remaining
      W_sgenes_remaining = setdiff(W_sgenes_remaining, chosen_sgene)
      Z_sgenes_remaining = setdiff(Z_sgenes_remaining, chosen_sgene)
    }
    
    
    # double check
    
    # no intersecting between Z and W
    assertthat::assert_that(length(intersect(W_sgenes, Z_sgenes)) == 0, msg = 'checking: no intersection between Z and W')
    # chosen W_sgenes are from original list
    assertthat::assert_that(all(W_sgenes %in% W_sgenes_nonviol$supergene_membership), msg = 'checking: chosen W_sgenes are from original list')
    # chosen Z_sgenes are from original list
    assertthat::assert_that(all(Z_sgenes %in% Z_sgenes_nonviol$supergene_membership), msg = 'checking: chosen Z_sgenes are from original list')

    
    # # allow setting$NUM_NCE/O=NA option
    # if(is.na(setting$NUM_NCE)) { # no sampling, take all
    #   Z = all_Y |> filter(Y_chr %in% Z_chrs)
    # } else {
    #   Z = all_Y |> filter(Y_chr %in% Z_chrs) |> slice_sample(n = setting$NUM_NCE)
    # }
    # 
    # if(is.na(setting$NUM_NCO)) {
    #   W = all_Y |> filter(Y_chr %in% W_chrs)
    # } else {
    #   W = all_Y |> filter(Y_chr %in% W_chrs) |> slice_sample(n = setting$NUM_NCO)
    # }
    
    
    
    # A_name$Y_name$ZW_idx$...
    AYZW[[A_name]][[Y_name]][[j]] = list()
    AYZW[[A_name]][[Y_name]][[j]][['AY_idx']] = i
    
    AYZW[[A_name]][[Y_name]][[j]][['A_chr']] = A_chr
    AYZW[[A_name]][[Y_name]][[j]][['Y_chr']] = Y_chr
    
    AYZW[[A_name]][[Y_name]][[j]][['Z_sgenes']] = Z_sgenes
    AYZW[[A_name]][[Y_name]][[j]][['W_sgenes']] = W_sgenes
    
    # AYZW[[A_name]][[Y_name]][[j]][['Z_chrs']]  = Z$Y_chr
    # AYZW[[A_name]][[Y_name]][[j]][['W_chrs']]  = W$Y_chr
  }
}


# AYZW
# object.size(AYZW)
# names(AYZW)
# names(AYZW$PDL1g1)
# names(AYZW$ATF2g4)
# names(AYZW$ATF2g4$SERPINE2)
# 
# AYZW$ATF2g4$SERPINE2[[1]]
# length(AYZW$ATF2g4$SERPINE2)




# =================== Saving AYZW names ======================================
print(sprintf("[%s]    - Saving AYZW names", Sys.time()))
print(sprintf("[%s]        - nrow(AY) = %s", Sys.time(), nrow(AY)))
write.csv(AY, sprintf('%s/supergene/%s/AY.csv', save_dir, AYZW_setting_name), row.names = FALSE)
saveRDS(AYZW, sprintf('%s/supergene/%s/AYZW.rds', save_dir, AYZW_setting_name))


# =================== END ======================================================
print(sprintf("[%s] END", Sys.time()))
