# R script for loading in PATHS to 
# data directory
# save directory
# etc...
# Requires DEVICE argument that is one of: 'laptop', 'desktop', 'ubergenno'

# location of papalexi-2021 folder
data_dir = switch(DEVICE,
               'laptop'='C:/Users/Cathe/Documents/School/Genetic Hypothesis Testing using Negative Controls/genData/papalexi', 
               'ubergenno'='<UPDATE!>')

# location of intermediate save files/plots are written
save_dir = switch(DEVICE,
               'laptop'='C:/Users/Cathe/Documents/School/Genetic Hypothesis Testing using Negative Controls/Proximal2025/papalexi-2021/saves', 
               'ubergenno'='<UPDATE!>')

# really should be named util_dir
util_dir = switch(DEVICE,
                  'laptop'='C:/Users/Cathe/Documents/School/Genetic Hypothesis Testing using Negative Controls/Proximal2025/utils', 
                  'ubergenno'='<UPDATE!>')


assertthat::assert_that(!is.null(data_dir), msg='DEVICE must be: laptop or ubergenno')
