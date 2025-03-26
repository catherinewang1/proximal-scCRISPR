# R script for loading in PATHS to 
# data directory
# save directory
# etc...
# Requires DEVICE argument that is one of: 'laptop', 'desktop', 'ubergenno'

# location of papalexi-2021 folder
data_dir = switch(DEVICE,
               'laptop'='C:/Users/Cathe/Documents/School/Genetic Hypothesis Testing using Negative Controls/genData/papalexi', 
               'desktop'='C:/Users/Catherine W/Documents/Research/genData/papalexi<UPDATE!>', 
               'ubergenno'='/raid6/Catherine/papalexi<UPDATE!>')
# location of intermediate save files/plots are written
save_dir = switch(DEVICE,
               'laptop'='C:/Users/Cathe/Documents/School/Genetic Hypothesis Testing using Negative Controls/proximal-scCRISPR-github/proximal-scCRISPR/papalexi-2021/saves', 
               'desktop'='C:/Users/Catherine W/Documents/Research/DoubleBridge/saves/papalexi_saves<UPDATE!>', 
               'ubergenno'='/raid6/Catherine/papalexi/papalexi_saves<UPDATE!>')
# # location of utils code
# CODE_DIR = switch(DEVICE,
#                   'laptop'='C:/Users/Cathe/Documents/School/Genetic Hypothesis Testing using Negative Controls/Proximal/utils/dataProcessing', 
#                   'desktop'='C:/Users/Catherine W/Documents/Research/DoubleBridge/code/utils', 
#                   'ubergenno'='/raid6/home/catheri2/DoubleBridge/code/utils')

# really should be named util_dir
util_dir = switch(DEVICE,
                  'laptop'='C:/Users/Cathe/Documents/School/Genetic Hypothesis Testing using Negative Controls/proximal-scCRISPR-github/proximal-scCRISPR/utils', 
                  'desktop'='C:/Users/Catherine W/Documents/Research/DoubleBridge/code/utils', 
                  'ubergenno'='/raid6/home/catheri2/DoubleBridge/code/utils')


assertthat::assert_that(!is.null(data_dir), msg='first arg must be: laptop, desktop, or ubergenno')
