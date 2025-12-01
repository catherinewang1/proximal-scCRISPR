# R script for loading in PATHS to 
# data directory
# save directory
# etc...
# Requires DEVICE argument that is one of: 'laptop', 'desktop', 'ubergenno'

# location of papalexi-2021 folder
data_dir = switch(DEVICE,
               'laptop'='C:/Users/Cathe/Documents/School/Genetic Hypothesis Testing using Negative Controls/genData/papalexi', 
               'ubergenno'='/raid6/Catherine/papalexi')

# location of intermediate save files/plots are written
save_dir = switch(DEVICE,
               'laptop'='C:/Users/Cathe/Documents/School/Genetic Hypothesis Testing using Negative Controls/proximal-scCRISPR-github/proximal-scCRISPR/papalexi-2021/saves', 
               'ubergenno'='/home/catheri2/proximal-scCRISPR-github/papalexi-2021/saves')

# really should be named util_dir
util_dir = switch(DEVICE,
                  'laptop'='C:/Users/Cathe/Documents/School/Genetic Hypothesis Testing using Negative Controls/proximal-scCRISPR-github/proximal-scCRISPR/utils', 
                  'ubergenno'='/home/catheri2/proximal-scCRISPR-github/utils')


assertthat::assert_that(!is.null(data_dir), msg='DEVICE must be: laptop or ubergenno')
