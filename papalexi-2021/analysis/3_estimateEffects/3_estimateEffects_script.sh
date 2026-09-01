# nice -10 Rscript 3.1_papalexi_continuous.R ubergenno A1 PCA-SPCA8.0-SPCA34.5-WGCNA-singlegene
# nice -10 Rscript 3.1_papalexi_continuous.R ubergenno A PCA-SPCA8.0-SPCA34.5-WGCNA-singlegene
# nice -10 Rscript 3.1_papalexi_continuous.R ubergenno B PCA-SPCA8.0-SPCA34.5-WGCNA-singlegene
# nice -10 Rscript 3.1_papalexi_continuous.R ubergenno C1 PCA-SPCA8.0-SPCA34.5-WGCNA-singlegene
# nice -10 Rscript 3.1_papalexi_continuous.R ubergenno C PCA-SPCA8.0-SPCA34.5-WGCNA-singlegene

# nice -10 Rscript 3.2_papalexi_count.R ubergenno A1 simpleCount-proximalNegBinCountCount-proximalNegBinSinglegeneSinglegene-proximalNegBinPCAPCA
# nice -10 Rscript 3.2_papalexi_count.R ubergenno A simpleCount-proximalNegBinCountCount-proximalNegBinSinglegeneSinglegene-proximalNegBinPCAPCA
# nice -10 Rscript 3.2_papalexi_count.R ubergenno B simpleCount-proximalNegBinCountCount-proximalNegBinSinglegeneSinglegene-proximalNegBinPCAPCA
# nice -10 Rscript 3.2_papalexi_countGLM.R ubergenno C1
# nice -10 Rscript 3.2_papalexi_countGLM.R ubergenno C

# Rscript 3.2_papalexi_countGLM.R laptop A1
# Rscript 3.2_papalexi_countGLM.R laptop A
# Rscript 3.2_papalexi_countGLM.R laptop B
# Rscript 3.2_papalexi_countGLM.R laptop C1
# Rscript 3.2_papalexi_countGLM.R laptop C



# Perform analysis setting by setting (so that we have a whole setting's analysis done first)
# nice -10 Rscript 3.1_papalexi_continuous.R ubergenno A1 PCA-SPCA8.0-SPCA34.5-WGCNA-singlegene
# nice -10 Rscript 3.2_papalexi_count.R      ubergenno A1 simpleCount-proximalNegBinCountCount-proximalNegBinSinglegeneSinglegene-proximalNegBinPCAPCA

# nice -10 Rscript 3.1_papalexi_continuous.R ubergenno A PCA-SPCA8.0-SPCA34.5-WGCNA-singlegene
# nice -10 Rscript 3.2_papalexi_count.R      ubergenno A simpleCount-proximalNegBinCountCount-proximalNegBinSinglegeneSinglegene-proximalNegBinPCAPCA


# nice -10 Rscript 3.1_papalexi_continuous.R ubergenno B PCA-SPCA8.0-SPCA34.5-WGCNA-singlegene
# nice -10 Rscript 3.2_papalexi_count.R      ubergenno B simpleCount-proximalNegBinCountCount-proximalNegBinSinglegeneSinglegene-proximalNegBinPCAPCA



# ================================
# Standard Analyses first
# ================================
# nice -10 Rscript 3.1_papalexi_continuous.R macbook A1 simpleLinear
# nice -10 Rscript 3.2_papalexi_count.R      macbook A1 simpleCount

# nice -10 Rscript 3.1_papalexi_continuous.R macbook A simpleLinear
# nice -10 Rscript 3.2_papalexi_count.R      macbook A simpleCount

# nice -10 Rscript 3.1_papalexi_continuous.R macbook B simpleLinear
# nice -10 Rscript 3.2_papalexi_count.R      macbook B simpleCount




# ==============================================
# Small subset of proximal continuous methods 
# ==============================================
# nice -10 Rscript 3.1_papalexi_continuous.R macbook A1 PCA-singlegene
# nice -10 Rscript 3.1_papalexi_continuous.R macbook A  PCA-singlegene


# nice -10 Rscript 3.2_papalexi_count.R macbook A1 simpleCount-proximalNegBinSinglegeneSinglegene-proximalNegBinPCAPCA
# nice -10 Rscript 3.2_papalexi_count.R macbook A  simpleCount-proximalNegBinSinglegeneSinglegene-proximalNegBinPCAPCA


# ==============================================
# Limited Methods Comparison: ubergenno
# (continuous: LinReg, PCA, singlegene
#  count: Pois, NegBin, singlegene, PCA)
# ==============================================
# A1
nice -10 Rscript 3.1_papalexi_continuous.R macbook A1 PCA-singlegene
nice -10 Rscript 3.2_papalexi_count.R macbook      A1 simpleCount-proximalNegBinSinglegeneSinglegene-proximalNegBinPCAPCA

# A
nice -10 Rscript 3.1_papalexi_continuous.R macbook A PCA-singlegene
nice -10 Rscript 3.2_papalexi_count.R macbook      A simpleCount-proximalNegBinSinglegeneSinglegene-proximalNegBinPCAPCA


# B
nice -10 Rscript 3.1_papalexi_continuous.R macbook B PCA-singlegene
nice -10 Rscript 3.2_papalexi_count.R macbook      B simpleCount-proximalNegBinSinglegeneSinglegene-proximalNegBinPCAPCA










