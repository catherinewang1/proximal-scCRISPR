"""
Concatenate the 2 datasets with different perturbations together. 
Filter for the speicified celltypes (neuron type and also neuron subgroup using the neuron_type and predicted_group columns respectively)

"""


import os

# Parameters for creating this subset of data
PARAMS = {
  'TAG': 'L45Glut',
  'subfolder': '3_L45glutcells',
  'RUN_EVEN_IF_SAVED'   : True,
  'seed':12345,
  'GRNA_MIN_CELLS'      : 20, # filter for grnas with at least this number of cells (keep low here, set higher later)
  'GENE_MIN_NONZERO_TRTCELLS' : 3000, # filter for genes with at least this number of TREATED cells (not Non_target)
  'neuron_type' : "Glut",                     # subset for this neuron type
  'neuron_subgroup' : "005 L4-5 IT CTX Glut", # subset for this neuron subtype predicted_group=='005 L4-5 IT CTX Glut’
  'NUM_SUBSAMPLE_CONTROL_CELLS': 10000, 
  'NUM_FAKE_GRNA': 10,
  'FAKE_GRNA_SIZE': 10, # number of nt grnas to constrruct fake grnas
  'FAKE_GRNA_MIN_NUMCELLS': 100   # minimum size of fake grna
}


# Paths for original jin-2025 dataset
PATHS = {'data_dir':  '/Users/catherinewang/Documents/School/genData/jin-2025/',
         'save_dir': '/Users/catherinewang/Documents/School/proximal-scCRISPR-github/jin-2025/saves',
         'util_dir': '/Users/catherinewang/Documents/School/proximal-scCRISPR-github/utils'}
PATHS['save_dir'] = PATHS['save_dir'] + '/' + PARAMS['subfolder']


os.makedirs(PATHS['save_dir'], exist_ok=True) # create dir for saved subsetted data

h5adfilename = 'jin2025_' + PARAMS['TAG'] + '.h5ad'

# if PARAMS['RUN_EVEN_IF_SAVED'] or not h5adfilename in os.listdir(PATHS['save_dir']):
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt


# load and concatenate 2 datasets
jin2025_10pert = sc.read_h5ad(PATHS['data_dir'] + 'subset_adata_SCARF_10celltypes_10perturbations.h5ad')
jin2025_48pert = sc.read_h5ad(PATHS['data_dir'] + 'subset_adata_SCARF_10celltypes_48perturbations.h5ad')
jin2025_10pert.obs_names = jin2025_10pert.obs_names + ':10pert' # change cellnames to specify data source
jin2025_48pert.obs_names = jin2025_48pert.obs_names + ':48pert'

jin2025 = sc.concat({'10pert': jin2025_10pert, '48pert': jin2025_48pert}, axis='obs') # combine into h5ad format



# =================================================================================================================
# Subset for celltype and then do initial filtering for grna size and gene non-zero size
# new jin2025_celltype anndata object
# =================================================================================================================

# subset for celltype and celltype subgroup
if PARAMS['neuron_type']!=None and PARAMS['neuron_subgroup']== None: 
    # choose only neuron type
    subtype_cells = jin2025.obs.query(f"neuron_type=={PARAMS['neuron_type']}").index
elif PARAMS['neuron_type']==None and PARAMS['neuron_subgroup']!= None:  
    # choose on only celltype subgroup
    subtype_cells = jin2025.obs.query(f" predicted_group=={PARAMS['neuron_subgroup']}").index
elif PARAMS['neuron_type']!=None and PARAMS['neuron_subgroup']!= None:   
    # choose on both
    subtype_cells = jin2025.obs.query(f"neuron_type=='{PARAMS['neuron_type']}' and predicted_group=='{PARAMS['neuron_subgroup']}'").index
else: 
    # both are None, use all cells (probably not good! this is too many)
    subtype_cells = jin2025.obs.index



# subtype_cells
# subset for glut cells
jin2025_celltype = jin2025[subtype_cells, ]


# Filter grnas initially based on sample size
grna_counts = jin2025_celltype.obs['gene_target'].value_counts()
grna_subset = grna_counts[grna_counts>PARAMS['GRNA_MIN_CELLS']].index.tolist()

grna_subset_cells = jin2025_celltype.obs[jin2025_celltype.obs['gene_target'].isin(grna_subset)].index
jin2025_celltype = jin2025_celltype[grna_subset_cells, ] # UPDATE DATA


# Filter genes initially based on sample size on treated cells
trt_cellnames = jin2025_celltype.obs.loc[jin2025_celltype.obs['gene_target'] != 'Non_target'].index
gene_counts_trt = jin2025_celltype[trt_cellnames, ].X.count_nonzero(0)
gene_subset = jin2025_celltype.var_names[gene_counts_trt >= PARAMS['GENE_MIN_NONZERO_TRTCELLS']]
gene_subset = np.unique(np.append(gene_subset, np.setdiff1d(grna_subset, 'Non_target'))) # add grna targets if filtered out here

# # look at dist of nonzero gene counts
# plt.hist(gene_counts_trt, bins=30, color='skyblue', edgecolor='black')
# plt.show()


# for overall
# gene_counts = pd.Series((jin2025_celltype.X!=0).sum(0).tolist()[0]) # slow
# gene_counts = jin2025_celltype.X.count_nonzero(0)



jin2025_celltype = jin2025_celltype[:, gene_subset] # UPDATE DATA

# =================================================================================================================
# For the 399 NT grnas, reserve some for overall control cells and some for fake perturbations
# For fake perturbations, combine many, PARAMS['FAKE_GRNA_SIZE'], perturbations together
# and subsample for sizes similar to actual treated (but filter for at least PARAMS['FAKE_GRNA_MIN_NUMCELLS'])
# In the h5ad subset, create a new gene_target_w_fakeNT column that adds the new fake NT groupings 
# =================================================================================================================

NT_names = jin2025_celltype.obs.query("gene_target == 'Non_target'")['guide_call'].unique().to_numpy()

rng = np.random.default_rng(PARAMS['seed'])
rng.shuffle(NT_names)

# trt sizes for the fake NT
trt_sizes = grna_counts[(grna_counts>= PARAMS['FAKE_GRNA_MIN_NUMCELLS']) & (grna_counts.index != 'Non_target')].tolist()


fakeNT_cellnames = np.array([])
fake_grna_map = {} 

for i in np.arange(PARAMS['NUM_FAKE_GRNA']): 
    # choose NT combo for this fake NT
    cur_NT_names = NT_names[(i*PARAMS['FAKE_GRNA_SIZE']) + np.arange(PARAMS['FAKE_GRNA_SIZE'])]
    cur_fakeNT_name = '&'.join(cur_NT_names)

    # subset cells for this fake NT
    cur_fakeNT_samplesize = rng.choice(trt_sizes, 1)[0]
    cur_fakeNTcellnames = jin2025_celltype.obs[jin2025_celltype.obs['guide_call'].isin(cur_NT_names)].index.tolist()
    cur_fakeNTcellnames = rng.choice(a=cur_fakeNTcellnames,   size = min(cur_fakeNT_samplesize, len(cur_fakeNTcellnames)), replace=False) # subsample if too many (should usually subsample)

    fakeNT_cellnames = np.append(fakeNT_cellnames, cur_fakeNTcellnames) # add the chosen fake NT cells

    for nt_ in cur_NT_names:
        fake_grna_map[nt_] = cur_fakeNT_name


NT_controls = np.setdiff1d(NT_names, list(fake_grna_map.keys())) # rest of the NT names


# create new col that includes fake NT
gene_target_w_fakeNT = [jin2025_celltype.obs['gene_target'].iloc[i, ] if (fake_grna_map.get(jin2025_celltype.obs['guide_call'].iloc[i, ])==None) 
                        else fake_grna_map.get(jin2025_celltype.obs['guide_call'].iloc[i, ]) 
                        for i in np.arange(jin2025_celltype.obs.shape[0])]

jin2025_celltype.obs['gene_target_w_fakeNT'] = gene_target_w_fakeNT


# subsample rest of controls
trt_cells   = jin2025_celltype.obs.query("gene_target_w_fakeNT != 'Non_target'").index.tolist()
cntrl_cells = jin2025_celltype.obs[jin2025_celltype.obs['guide_call'].isin(NT_controls)].index
cntrl_cells = rng.choice(cntrl_cells, size=min(PARAMS['NUM_SUBSAMPLE_CONTROL_CELLS'], len(cntrl_cells)), replace=False).tolist()
subset_cells = trt_cells + cntrl_cells

# subset
jin2025_subset = jin2025_celltype[subset_cells, ]

# save
jin2025_subset.write_h5ad(
    f'{PATHS['save_dir']}/jin2025_{PARAMS['TAG']}.h5ad',
    # PATHS['save_dir'] + '/jin2025_.h5ad',
    compression='lzf'
)




