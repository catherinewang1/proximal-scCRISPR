


require(ondisc)
require(dplyr)
require(ggplot2)







#' Construct a dataframe showing which gRNA each cell received
#' 
#' @param grna_odm (ondisc odm) ondisc package's manager for interacting with data
#'            e.g. grna_odm <- read_odm(odm_fp      = "assignment_matrix.odm",
#'                                      metadata_fp = "metadata.rds")
#' @return dataframe with each row indicating cell index, gRNA received, features of gRNA received
#'   colnames: cell_idx, gRNA_idx, gRNA, target_type
get_grna_assignment <- function(grna_odm) {
  grna <- grna_odm[[,1:ncol(grna_odm)]] # entire grna ds (eg 110 grnas x 20729 cells)
  
  # make df of which grna each cell got
  grna_key = data.frame(gRNA_idx = 1:nrow(grna_odm), 
                        gRNA = grna_odm |> get_feature_covariates() |> rownames(),
                        target_type = grna_odm |> get_feature_covariates() |> select(target_type))
  grna_assignment_0 = data.frame(cell_idx = 1:ncol(grna_odm), 
                                 gRNA_idx = apply(grna, MARGIN = 2, FUN = which))
  
  grna_assignment = merge(grna_assignment_0, grna_key, by = 'gRNA_idx') |> 
    arrange(cell_idx) |> 
    select(cell_idx, gRNA_idx, gRNA, target_type)
  
  grna_assignment
}

#' Construct 2 dataframes:
#' allcells : concatenating info about grna assignment, 2 dimensional dim reduction, cell covariates
#' gRNAcells: subsetting df_allcells to those receiving specified gRNA and NT. Also add info on 
#'         is_targeting (T/F) whether received gRNA is targeting or not
#'         is_exposure (1/0) whether a cell is 'exposed' (whether a separate NT cell should be matched
#'                     to this cell). Randomly, some NT cells will be chosen to be added to all gRNA cells
#'         gRNA_label (string) simplify gRNA column. All 'NTgX' values will be changed to 'NT'
#' @param grna_assignment (dataframe) dataframe depicting which gRNA each cell received
#' @param cell_covariates (dataframe) features of cells (eg potential confounding)
#' @param gRNA_names (vector of strings) gRNA chosen to focus on/make plots of
#' @param NT_names (vector of strings) names of NTgXX 
#' @param NT_exposed_size (int) number of NT cells to be chosen to be 'exposed' (number of NT cells
#'                              to be put into set of cells that will receive a match) if not given,
#'                              choose the size of the largest gRNA group
#' @param seed (int) random seed for choosing some NT cells for exposure 
get_dimred <- function(grna_assignment,
                       df_dimred,
                       cell_covariates,
                       gRNA_names,
                       NT_names,
                       NT_exposed_size=NULL,
                       seed = 12345) {
  # # debugging
  # num_imp_genes = 4000              # number of important genes (for saving name)
  # dimred_type = 'pca'               # type of dim reduction (pca, umap, pca_umap)
  # grna_assignment = grna_assignment # dataframe cell x feat describing gRNA assignment for each cell
  # df_dimred = as.data.frame(papalexi_pca$x[,1:2]) # cell x 2 dim reduction
  # cell_covariates = gene_odm |> get_cell_covariates() # covariates for cells (potential confounders)
  # gRNA_names = chosen_grna_T
  # NT_names = all_grna_NT
  # seed = 12345
  
  colnames(df_dimred) = c('dim1', 'dim2')
  
  df_allcells = cbind(grna_assignment, df_dimred, cell_covariates) 
  
  df = df_allcells |> 
    filter(gRNA %in% c(gRNA_names, NT_names)) |> 
    mutate(is_targeting = (target_type == 'gene')) |> 
    mutate(is_exposure = is_targeting)
  
  # choose some NT to be 1 in exposure col (used in matching as response)
  # every Targeting cell + these chosen NT cells, will receive a matched other NT cell
  set.seed(seed) 
  if(is.null(NT_exposed_size)) {
    NT_exposed_size = df |> filter(gRNA %in% gRNA_names) |> 
      group_by(gRNA) |> summarise(count = n()) |> 
      pull(count) |> max()
  }
  NT_exposure_idx = sample(which(!df$is_targeting), 
                           size = NT_exposed_size, 
                           replace = F)
  df$is_exposure[NT_exposure_idx] = 1 # take some NT to be exposure
  
  # make col of gRNA assignments, combining NTgX values to NT
  df = df |> mutate(gRNA_label = gRNA)
  df[df$target_type == 'non-targeting', 'gRNA_label'] = 'NT'
  
  return(list(allcells=df_allcells, 
              gRNAcells=df))
}

#' Plot dimension reductions for the specified input parameters
#' (first 5 inputs are for plotting names and saved plot names)
#' @param save_dir (string) file directory of saves
#' @param num_imp_genes (int) number of important genes chosen 
#' @param dimred_type (string) type of dimension reduction performed (pca, umap, pca_umap)
#' @param gRNA_names (vector of strings) gRNA chosen to focus on/make plots of
#' @param NT_names (vector of strings) names of NTgXX 
#' @param gRNAcells (dataframe) second return of get_dimred. cells with given gRNA or NT 
#'                              dim reduction values
#' @param allcells (dataframe) first return of get_dimred. all cells
plot_dimred <- function(save_dir, 
                        num_imp_genes,
                        dimred_type,
                        gRNA_names,
                        NT_names,
                        gRNAcells,
                        allcells) {
  
  x_label = paste0(dimred_type, '1')
  y_label = paste0(dimred_type, '2')
  
  # create new folder to save imgs in
  # dir.create(sprintf('%s/papalexi_saves/%s/%s/', 
  #                    save_dir, dimred_type, num_imp_genes), 
  #            showWarnings = FALSE)
  dir.create(sprintf('%s/%s/%s/', 
                     save_dir, dimred_type, num_imp_genes), 
             showWarnings = FALSE, recursive = TRUE)
  
  # Plot single gRNA vs NT
  for(i in 1:length(gRNA_names)) {
    gRNA_i = gRNA_names[i]
    # gRNA i vs NT
    ggplot(gRNAcells |> filter(gRNA %in% c(gRNA_i, all_grna_NT)), 
           aes(x = dim1, y = dim2, color = gRNA_label)) +
      geom_point(alpha = .4, size = .75) +
      labs(title = paste0(gRNA_i, ' vs NT'), color = 'gRNA',
           x = x_label, y = y_label) +
      theme_bw() +
      theme(plot.title = element_text(hjust = .5),
            legend.position = c(.85, .84),
            legend.box.background = element_rect(colour = "black"))
    
    
    ggsave(filename=sprintf('%s/%s/%s/%s_%s_gRNA_%s.png', 
                            save_dir, dimred_type, num_imp_genes, dimred_type, num_imp_genes, gRNA_i),
           height = 5, width = 6, dpi = 320)
  }
  # Plot all gRNA vs NT
  ggplot(gRNAcells, 
         aes(x = dim1, y = dim2, color = gRNA_label)) +
    geom_point(alpha = .3, size = .6) +
    labs(title = paste0('some gRNA', ' vs NT'), color = 'gRNA',
         x = x_label, y = y_label) +
    theme_bw() +
    theme(plot.title = element_text(hjust = .5),
          legend.position = c(.85, .84),
          legend.box.background = element_rect(colour = "black"))
  
  
  ggsave(filename=sprintf('%s/%s/%s/%s_%s_gRNA_all.png', 
                          save_dir, dimred_type, num_imp_genes, dimred_type, num_imp_genes),
         height = 5, width = 6, dpi = 320)
  
  # Plot by potential confounders
  for(color_name in c('n_nonzero', 'n_umis', 'lane', 'bio_rep', 'phase', 'p_mito')) {
    ggplot(allcells,
           aes(x = dim1, y = dim2, 
               color=eval(parse(text = color_name)))) +
      geom_point(alpha=.2, size = .75) +
      labs(title = sprintf('%s - %s', num_imp_genes,color_name),
           x = x_label, y = y_label,
           color = color_name) +
      theme_bw() +
      theme(plot.title = element_text(hjust = .5),
            legend.position = c(.9, .84),
            legend.box.background = element_rect(color = "black"),
            legend.key.size = unit(.45, 'cm'))
    
    ggsave(filename=sprintf('%s/%s/%s/%s_%s_cell_%s.png', 
                            save_dir, dimred_type, num_imp_genes, dimred_type, num_imp_genes, color_name),
           height = 5, width = 6, dpi = 320)
    
  }
}



#' with same other inputs, make different df's with different dim reduction types
get_dimred_bydimredmethod <- function(grna_assignment,
                                      cell_covariates,
                                      gRNA_names,
                                      NT_names,
                                      NT_exposed_size=NULL,
                                      seed=12345) {
  inner_func <- function(df_dimred) {
    get_dimred(grna_assignment = grna_assignment, 
               df_dimred = df_dimred, 
               cell_covariates = cell_covariates,
               gRNA_names = chosen_grna_T,
               NT_names = all_grna_NT,
               NT_exposed_size=NT_exposed_size,
               seed=seed) 
  }
  
  return(inner_func)
}


plot_dimred_bydimredmethod <- function(save_dir, 
                                       num_imp_genes, 
                                       gRNA_names,
                                       NT_names) {
  
  inner_func <- function(dimred_type, gRNAcells, allcells) {
    plot_dimred(save_dir=save_dir, 
                num_imp_genes=num_imp_genes, 
                dimred_type = dimred_type, 
                gRNA_names=gRNA_names,
                NT_names=NT_names, 
                gRNAcells = gRNAcells, 
                allcells = allcells)
  }
  
  return(inner_func)
}



