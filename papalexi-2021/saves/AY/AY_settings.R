
# Settings for selecting perturbation (A) - gene (Y) tests
# 
# Creates a named list `settings` where each of the settings names are like A, A1, B, C, ... etc...
#   Indexed by the setting like settings$A or settings[['A']] which provides another list with
#   parameters for the analysis setting like specifying the number of AY tests, some data filtering thresholds, etc...
# 

settings = list(
  # Settings for 'A'
  "A" = list(seed = 1345678,
             # AY test settings
             NUM_ALT           = NA,   # num of Alt tests: NA if all As (probably SHOULD)
             NUM_NULL          = 1000, # num of Null tests (Non-Targeting to Y)
             NUM_MAYBE         = 200,  # num of Unknown tests: targeting perturbation A to randomly chosen gene Y (mostly null though)
             MAX_Y_IMPORTANCE  = 1500, # limit how 'unimportant' a response gene can be
             NUM_NCE           = 50,   # number of individual gene NCE (Z) per AY test 
             NUM_NCO           = 50,   # number of individual gene NCO (W) per AY test 
             NUM_NCENCO_per_AY = 1,    # number of NCE/NCO sets per AY test 
             # quality control settings- 
             PERTURBATION_N_NONZERO_CELLS_MIN = 100, # grna_chr$n_nonzero min (but allow alt tests)
             GENE_N_NONZERO_CELLS_MIN = 2000         # gene_metainfo$n_nonzero_cell min
  ),
  # Settings for 'A1'- small set to test code
  "A1" = list(seed = 1345678,
              # AY test settings
              NUM_ALT           = 10,   # num of Alt tests: NA if all As (probably SHOULD)
              NUM_NULL          = 20,  # num of Null tests (Non-Targeting to Y)
              NUM_MAYBE         = 10,    # num of Unknown tests: targeting perturbation A to randomly chosen gene Y (mostly null though)
              MAX_Y_IMPORTANCE  = 1500, # limit how 'unimportant' a response gene can be
              NUM_NCE           = 50,   # number of individual gene NCE (Z) per AY test 
              NUM_NCO           = 50,   # number of individual gene NCO (W) per AY test 
              NUM_NCENCO_per_AY = 1,    # number of NCE/NCO sets per AY test 
              # quality control settings- 
              PERTURBATION_N_NONZERO_CELLS_MIN = 100, # grna_chr$n_nonzero min (but allow alt tests)
              GENE_N_NONZERO_CELLS_MIN = 2000         # gene_metainfo$n_nonzero_cell min
  ),
  # Settings for 'B': Equal num of NULL and Maybe (this can compare NT-A and Targ-A tests)
  "B" = list(seed = 1345678,
             # AY test settings
             NUM_ALT           = NA,   # num of Alt tests: NA if all As (probably SHOULD)
             NUM_NULL          = 500, # num of Null tests (Non-Targeting to Y)
             NUM_MAYBE         = 500,  # num of Unknown tests: targeting perturbation A to randomly chosen gene Y (mostly null though)
             MAX_Y_IMPORTANCE  = 1500, # limit how 'unimportant' a response gene can be
             NUM_NCE           = 50,   # number of individual gene NCE (Z) per AY test 
             NUM_NCO           = 50,   # number of individual gene NCO (W) per AY test 
             NUM_NCENCO_per_AY = 1,    # number of NCE/NCO sets per AY test 
             # quality control settings- 
             PERTURBATION_N_NONZERO_CELLS_MIN = 100, # grna_chr$n_nonzero min (but allow alt tests)
             GENE_N_NONZERO_CELLS_MIN = 2000         # gene_metainfo$n_nonzero_cell min
  ),
  "TEST" = list(TESTINGFORMAT=NA)
)





# settings = list(
#     # Settings for 'A'
#     "A" = list(seed = 1345678,
#                # AY test settings
#                NUM_ALT           = NA,   # num of Alt tests: NA if all As (probably SHOULD)
#                NUM_NULL          = 1000, # num of Null tests (Non-Targeting to Y)
#                NUM_MAYBE         = 100,    # num of Unknown tests: targeting perturbation A to randomly chosen gene Y (mostly null though)
#                MAX_Y_IMPORTANCE  = 1500, # limit how 'unimportant' a response gene can be
#                NUM_NCE           = 50,   # number of individual gene NCE (Z) per AY test 
#                NUM_NCO           = 50,   # number of individual gene NCO (W) per AY test 
#                NUM_NCENCO_per_AY = 1,    # number of NCE/NCO sets per AY test 
#                # quality control settings- 
#                PERTURBATION_N_NONZERO_CELLS_MIN = 100, # grna_chr$n_nonzero min (but allow alt tests)
#                GENE_N_NONZERO_CELLS_MIN = 2000         # gene_metainfo$n_nonzero_cell min
#     ),
#     # Settings for 'A1'- small set to test code
#     "A1" = list(seed = 1345678,
#                # AY test settings
#                NUM_ALT           = 10,   # num of Alt tests: NA if all As (probably SHOULD)
#                NUM_NULL          = 50,  # num of Null tests (Non-Targeting to Y)
#                NUM_MAYBE         = 10,    # num of Unknown tests: targeting perturbation A to randomly chosen gene Y (mostly null though)
#                MAX_Y_IMPORTANCE  = 1500, # limit how 'unimportant' a response gene can be
#                NUM_NCE           = 50,   # number of individual gene NCE (Z) per AY test 
#                NUM_NCO           = 50,   # number of individual gene NCO (W) per AY test 
#                NUM_NCENCO_per_AY = 1,    # number of NCE/NCO sets per AY test 
#                # quality control settings- 
#                PERTURBATION_N_NONZERO_CELLS_MIN = 100, # grna_chr$n_nonzero min (but allow alt tests)
#                GENE_N_NONZERO_CELLS_MIN = 2000         # gene_metainfo$n_nonzero_cell min
#     ),
#     "TEST" = list(TESTINGFORMAT=NA)
# )