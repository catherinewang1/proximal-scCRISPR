# R script for loading in PATHS to 
# data directory
# save directory
# etc...
# Requires DEVICE argument that is one of: 'laptop', 'desktop', 'ubergenno'

# location of papalexi-2021 folder
data_dir = switch(DEVICE,
                  'laptop'='xxxxxxx',
                  'macbook'='/Users/catherinewang/Documents/School/genData/jin/',
                  'ubergenno'='xxxxxxx')

# location of intermediate save files/plots are written
save_dir = switch(DEVICE,
                  'laptop'='C:/Users/Cathe/Documents/School/Genetic Hypothesis Testing using Negative Controls/proximal-scCRISPR-github/proximal-scCRISPR/jin-2020/saves', 
                  'macbook'='/Users/catherinewang/Documents/School/proximal-scCRISPR-github/jin-2020/saves',
                  'ubergenno'='/home/catheri2/proximal-scCRISPR-github/jin-2020/saves')

# really should be named util_dir
util_dir = switch(DEVICE,
                  'laptop'='C:/Users/Cathe/Documents/School/Genetic Hypothesis Testing using Negative Controls/proximal-scCRISPR-github/proximal-scCRISPR/utils', 
                  'macbook'='/Users/catherinewang/Documents/School/proximal-scCRISPR-github/utils',
                  'ubergenno'='/home/catheri2/proximal-scCRISPR-github/utils')


assertthat::assert_that(!is.null(data_dir), msg='DEVICE must be: laptop or macbook or ubergenno')
