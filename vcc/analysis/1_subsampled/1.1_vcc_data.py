
from matplotlib import axis
import scanpy as sc
import numpy as np


# Paths for vcc dataset
PATHS = {'data_dir':  '/Users/catherinewang/Documents/School/genData/vcc/',
         'save_dir': '/Users/catherinewang/Documents/School/proximal-scCRISPR-github/vcc/saves',
         'util_dir': '/Users/catherinewang/Documents/School/proximal-scCRISPR-github/utils'}
PATHS['save_dir'] = PATHS['save_dir'] + '_subsample'



h5ad_file  = PATHS['data_dir'] + '/VCC-subsampled.h5ad' 

adata_subset = sc.read_h5ad(PATHS['data_dir'] + 'VCC-subsampled.h5ad')
adata        = sc.read_h5ad(PATHS['data_dir'] + 'VCC.h5ad')


adata.X   # <- gene expr
adata.obs # <- cell data, but no nUMIs/library size
adata_subset.X   # <- gene expr
adata_subset.obs # <- cell data, but no nUMIs/library size?



jin2020 = sc.read_h5ad(PATHS['data_dir'] + '../jin/causarrayexampledata/perturbseq-exneu.h5ad')
jin2020.obs

adata.obs.head

adata.layers
adata.n_obs
adata.n_vars
adata.obs
adata.obs_keys
adata.obsm
adata.obsm_keys
adata.raw
adata.var
adata.var_keys


adata.X.shape # (491046, 17516)

temp = adata.X[1:10, :]
temp.sum(axis=1)
np.array(temp).sum(axis = 2)
np.matrix(temp)
np.array(temp).view()



np.ndarray()


temp[1, 4]

adata_subset.X.shape

adata_subset.var_names
adata_subset.obs_names

adata_subset.obs_names
type(adata_subset.raw)


# choose the colls that appear in the subsample
adata.X.shape
adata.X[adata_subset.obs_names[1:5], 1:3]
adata_subset.obs_names[1:4]

subsampled_cells = np.where(np.isin(adata.obs_names, adata_subset.obs_names))[0]

X_subsampledcells_allgenes = adata.X[subsampled_cells, :]


X_subsampledcells_allgenes.shape


X_subsampledcells_allgenes[0:4, :].sum(axis=1)[:, 0]

adata_subset.X
adata_subset.obs
adata_subset
adata_subset.layers
adata_subset.n_obs
adata_subset.n_vars
adata_subset.obs
adata_subset.obs_keys
adata_subset.obsm
adata_subset.obsm_keys
adata_subset.raw
adata_subset.var
adata_subset.var_keys




adata_subset.X[0:4, :].sum(axis=1)[:, 0]





# # R
# # load h5ad
# h5ad_file  = sprintf('%s/VCC-subsampled.h5ad', data_dir)
# vcc_seurat = anndataR::read_h5ad(h5ad_file, as = "Seurat")

# gene_continuous = vcc_seurat[['RNA']]['X']
# cell_metadata   = vcc_seurat[[]]


# Sys.setenv(UV_OFFLINE=1)
# vcc_full = anndata::read_h5ad(sprintf('%s/VCC.h5ad', data_dir))