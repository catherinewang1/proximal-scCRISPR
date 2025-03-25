# ==================================================================================================
# =================== Choose Top XX nonTF genes ====================================================
# ==================================================================================================
# print(sprintf("[%s]    - Choose Top %d nonTF genes", Sys.time(), NUM_IMPORTANT_GENES))
myPrintColor(sprintf("[%s]    - Choose Top %d nonTF genes", Sys.time(), NUM_IMPORTANT_GENES))


# =================== load list of genes that are Transcription Factors ===========================
# print(sprintf("[%s]          loading list of TF genes", Sys.time()))
myPrintColor(sprintf("[%s]          loading list of TF genes", Sys.time()))
# read in xlsx sheet
tf_raw = readxl::read_xlsx(path = paste0(data_dir, "/../extra/transcriptionfactorlist.xlsx"),  # given from Kathryn over slack (suppl of a paper?)
                           sheet = 2) |> suppressMessages() # suppress messages on how they renamed columns
# clean up a bit
#      2nd col ('...2') is gene name, 4th col is TF indicator 'Is TF?'
#      First row is not data (2 rows of colnames)
tf = tf_raw[-1 , c('...2', 'Is TF?')]
colnames(tf) = c('gene_name', 'TF')
# table(tf$TF) # only levels are No, Yes
#        No  Yes 
#       1126 1639 
TF_names = tf |> dplyr::filter(TF == 'Yes') |> dplyr::pull(gene_name)