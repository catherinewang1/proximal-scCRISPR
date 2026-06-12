
# Settings for selecting perturbation (A) - gene (Y) tesets

settings = list(
    # Settings for 'A'
    "A" = list(seed = 1345678,
               NUM_A             = NA, # NA if all As (probably SHOULD)
               NUM_Y_PER_A_NEG   = 20,  # NA if all (probably should NOT)
               NUM_Y_PER_A_MAYBE = 0, # NA if all on same chromosome (probably SHOULD)
               MAX_Y_IMPORTANCE  = 1500, # limit how 'unimportant' a response gene can be
               NUM_NCENCO_pairs  = 100,  # number of NCE/NCO pairs (dimU, length of ZWs)
               NUM_NCE           = NA,  # number of NCE per AY test (prev #NCE/NCO equal) (NA=all avail)
               NUM_NCO           = NA,  # number of NCO per AY test 
               NUM_NCENCO_per_AY = 1,   # number of NCE/NCO sets per AY test 
               NUM_AY_POS        = NA,    # number of known causal/positive AY tests, NA if all As (probably SHOULD)
               PERTURBATION_N_NONZERO_CELLS_MIN = 100, # grna_chr$n_nonzero min (allo pos tests)
               GENE_N_NONZERO_CELLS_MIN = 2000  # gene_metainfo$n_nonzero_cell min
    ),
    # Settings for 'A1'- small set to test code
    "A1" = list(seed = 1345678,
               NUM_A             = 5, # NA if all As (probably SHOULD)
               NUM_Y_PER_A_NEG   = 5,  # NA if all (probably should NOT)
               NUM_Y_PER_A_MAYBE = 0, # NA if all on same chromosome (probably SHOULD)
               MAX_Y_IMPORTANCE  = 2000, # limit how 'unimportant' a response gene can be
               NUM_NCENCO_pairs  = 10,  # number of NCE/NCO pairs (dimU, length of ZWs)
               NUM_NCE           = NA,  # number of NCE per AY test (prev #NCE/NCO equal) (NA=all avail)
               NUM_NCO           = NA,  # number of NCO per AY test 
               NUM_NCENCO_per_AY = 1,   # number of NCE/NCO sets per AY test
               NUM_AY_POS        = 5,    # number of known causal/positive AY tests, NA if all As (probably SHOULD)
               PERTURBATION_N_NONZERO_CELLS_MIN = 100, # grna_chr$n_nonzero min
               GENE_N_NONZERO_CELLS_MIN = 2000 # gene_metainfo$n_nonzero_cell min
    ),
    # Settings for 'B'- more AY tests, require grnas have larger sample size
    "B" = list(seed = 54321,
               NUM_A             = NA, # NA if all As (probably SHOULD)
               NUM_Y_PER_A_NEG   = 50,  # NA if all (probably should NOT)
               NUM_Y_PER_A_MAYBE = 0, # NA if all on same chromosome (probably SHOULD)
               MAX_Y_IMPORTANCE  = 2000, # limit how 'unimportant' a response gene can be
               NUM_NCENCO_pairs  = 100,  # number of NCE/NCO pairs (dimU, length of ZWs)
               NUM_NCE           = NA,  # number of NCE per AY test (prev #NCE/NCO equal) (NA=all avail)
               NUM_NCO           = NA,  # number of NCO per AY test 
               NUM_NCENCO_per_AY = 1,   # number of NCE/NCO sets per AY test
               NUM_AY_POS        = NA,    # number of known causal/positive AY tests, NA if all As (probably SHOULD)
               PERTURBATION_N_NONZERO_CELLS_MIN = 100, # grna_chr$n_nonzero min (allow pos tests)
               GENE_N_NONZERO_CELLS_MIN = 2000  # gene_metainfo$n_nonzero_cell min
    ),
    # Settings for 'C'- force a matrix: choose NUM_A grnas and  perturbations
    "C" = list(seed = 54321,
               NUM_A             = NA,   # NA if all As (probably SHOULD)
               FORCE_AY_MATRIX   = TRUE, # force chosen A and Y to be all cross comb of set of A and Y
               NUM_Y_PER_A_NEG   = 40,  # NA if all (probably should NOT)
               NUM_Y_PER_A_MAYBE = 0,   # NA if all on same chromosome (probably SHOULD)
               MAX_Y_IMPORTANCE  = 1500, # limit how 'unimportant' a response gene can be
               NUM_NCENCO_pairs  = 100,  # number of NCE/NCO pairs (dimU, length of ZWs)
               NUM_NCE           = NA,  # number of NCE per AY test (prev #NCE/NCO equal) (NA=all avail)
               NUM_NCO           = NA,  # number of NCO per AY test 
               NUM_NCENCO_per_AY = 1,   # number of NCE/NCO sets per AY test
               NUM_AY_POS        = NA,    # number of known causal/positive AY tests, NA if all As (probably SHOULD)
               PERTURBATION_N_NONZERO_CELLS_MIN = 100, # grna_chr$n_nonzero min (allow pos tests)
               GENE_N_NONZERO_CELLS_MIN = 2000  # gene_metainfo$n_nonzero_cell min
    ),
    # Settings for 'C'- force a matrix: choose NUM_A grnas and  perturbations, smaller for testing
    "C1" = list(seed = 54321,
               NUM_A             = 5,   # NA if all As (probably SHOULD)
               FORCE_AY_MATRIX   = TRUE, # force chosen A and Y to be all cross comb of set of A and Y
               NUM_Y_PER_A_NEG   = 5,  # NA if all (probably should NOT)
               NUM_Y_PER_A_MAYBE = 0,   # NA if all on same chromosome (probably SHOULD)
               MAX_Y_IMPORTANCE  = 1500, # limit how 'unimportant' a response gene can be
               NUM_NCENCO_pairs  = 100,  # number of NCE/NCO pairs (dimU, length of ZWs)
               NUM_NCE           = NA,  # number of NCE per AY test (prev #NCE/NCO equal) (NA=all avail)
               NUM_NCO           = NA,  # number of NCO per AY test 
               NUM_NCENCO_per_AY = 1,   # number of NCE/NCO sets per AY test
               NUM_AY_POS        = NA,    # number of known causal/positive AY tests, NA if all As (probably SHOULD)
               PERTURBATION_N_NONZERO_CELLS_MIN = 100, # grna_chr$n_nonzero min (allow pos tests)
               GENE_N_NONZERO_CELLS_MIN = 2000  # gene_metainfo$n_nonzero_cell min
    ),
    # Settings for 'TEST'
    "TEST" = list(TESTINGFORMAT=NA)

)