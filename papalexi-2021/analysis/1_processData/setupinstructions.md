


# Setup `ondisc` (and good for installing packages from source generally)



SETUP ON NEW DEVICE
ideal: change code to suit new functions from ondisc 
trying to use new functions
It doesn't seem to load in the format that the older odm version is saved in 
gene_odm <- ondisc::initialize_odm_from_backing_file(odm_file = paste0(data_dir, "/papalexi-2021/processed/gene/expression_matrix.odm"))

So the options are t try to clean the data to fit with the new format
 - (1) either CellRanger format or R matrix format
 - OR (2)  install the older package ondisc 1.1 (this was very easy on windows and linux. this took forever on mac)

## (1) Use new version of `ondisc` properly

Trying (1) loading in prev saved data using newer functions or seeing how hard it might be to change
library(Seurat)
library(SeuratData)
SeuratData::InstallData(ds = "thp1.eccite")
SeuratData::AvailableData() |> row.names()
SeuratData::InstalledData()
eccite_obj <- SeuratData::LoadData(ds = "thp1.eccite")

## (2) install older package version of ondisc 1.1

 - probably all my issues came from me trying to install from source (using the tar.gz file) 
   and trying to just specify the version from CRAN/github might've just fixed everything



macOS: https://mac.r-project.org/tools/


Need xcode-select and gfortran
  - install gfortran from the one given by R! not just newest on homebrew
Install these other packages using homebrew https://github.com/Rdatatable/data.table/wiki/Installation
  - brew install libomp
  - brew install llvm

Need to add to path https://stackoverflow.com/questions/69639782/installing-gfortran-on-macbook-with-apple-m1-chip-for-use-in-r
  - if ~/.R/Makevars file does not exist, create it. ~ is the home dir for username/    https://github.com/Rdatatable/data.table/issues/4437
  - whatever this is: i'm crying, this took so long. Needed to add to LDFLAGS the path to openssl installed. 


```
FC = /opt/gfortran/bin/gfortran
F77 = /opt/gfortran/bin/gfortran
FLIBS = -L/opt/gfortran/bin/gfortran -lgfortran -lquadmath -lm

# Newly installed Homebrew is located in
#  - /opt/homebrew for ARM Macs (M1 and its successors)
#  - /usr/local for Intel Macs
HOMEBREW_LOC=/opt/homebrew
# If you downloaded llvm manually above, replace with your chosen NEW_PATH/clang
LLVM_LOC = $(HOMEBREW_LOC)/Cellar/llvm/22.1.8
CC=$(LLVM_LOC)/bin/clang -fopenmp
CXX=$(LLVM_LOC)/bin/clang++ -fopenmp
# -O3 should be faster than -O2 (default) level optimisation ..
CFLAGS=-g -O3 -Wall -pedantic -std=gnu99 -mtune=native -pipe
CXXFLAGS=-g -O3 -Wall -pedantic -std=c++14 -mtune=native -pipe
LDFLAGS=-L$(HOMEBREW_LOC)/Cellar/gettext/1.0/lib -L$(HOMEBREW_LOC)/Cellar/openssl@3/3.6.3/lib -L$(LLVM_LOC)/lib -Wl,-rpath,$(LLVM_LOC)/lib
CPPFLAGS=-I$(HOMEBREW_LOC)/Cellar/gettext/1.0 -I$(LLVM_LOC)/include -I/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/include
```


code to install from source
```
install.packages(paste0(data_dir, '/ondisc-1.1.tar.gz'), repos = NULL, type ='source')
devtools::install(paste0(data_dir, '/ondisc-1.1/'))
install.packages("/Users/catherinewang/Documents/School/genData/papalexi/ondisc-1.1.tar.gz", repos = NULL, type ='source')

install.packages('/Users/catherinewang/Downloads/ondisc-1.1.tar.gz', repos=NULL, type='source')
install.packages('/Users/catherinewang/Downloads/quantreg_6.1.tar.gz', repos=NULL, type='source') testing random other package to install from source
```




# Install `pci2s` 



https://github.com/KenLi93/pci2s


Install using the github repo: `remotes::install_github("KenLi93/pci2s")`











