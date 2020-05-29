"""
# Preprocess scripts
Bring all datasets in a consistent format. Every dataset
needs to be saved in `results/data_processed` as a `h5ad` file.

## Input
* count data as downloaded from the public database
* sample/cell-level metadata wherever available

## Output
* a `scanpy` `AnnData` object for each dataset, containing raw (UMI) counts and the following metadata as `.obs`.

## Metadata
Metadata needs to have consistent column names. Columns in **bold** are mandatory.

* **samples**: unique identifier of each sample ('batch')
* **patient**: unique identifier of a patient
* **origin**, origin of the biopsy. controlled vocabs: `tumor_primary`, `normal_adjacent`, `tumor_edge`, `blood_peripheral`, `lymph_node`.
* **replicate**, same patient and origin, but different biopsy
* **platform**, experimental platform. controlled vocabs: `10x_3p_v1`, `10x_3p_v2`, `smartseq2`, `indrop_v2`, `10x_5p`
* **tumor_type**, tissue of origin. use TCGA cancer identifiers, such as `BRCA`.
* **dataset**, the dataset identifier. Used later to visualize batches.


A consistencyt check for the columns and their contents is implemented in
`lib/scio.py:check_obs`.
"""

import scanpy.api as sc
import anndata
import pandas as pd
import os
import sys
sys.path.append("lib")
from scio import concatenate, check_obs, check_var

DATASET = "vanderburg_01"
OBS_DATA = "tables/{}_samples.csv".format(DATASET)
OUTPUT_DIR = "./"

obs = pd.read_csv(OBS_DATA)


dataset_samples = obs["samples"].values


filenames = ["data/cellranger/{}_GEX/outs/raw_feature_bc_matrix.h5".format(sample[1:])
             for sample in dataset_samples]

adatas = [sc.read_10x_h5(filename, genome="GRCh38") for filename in filenames]

adatas2 = []
for adata, sample in zip(adatas, dataset_samples):
    duplicated = adata.var_names.duplicated()
    print("Removing {} gene symbols because they are duplicated".format(sum(duplicated)))
    adata = adata[:, ~duplicated].copy()
    adata.obs['samples'] = sample
    adatas2.append(adata)

adata = concatenate(adatas2, merge_var_cols=["gene_ids"])
adata.obs = adata.obs.join(obs.set_index("samples"), on="samples", how="left")

adata.obs["dataset"] = DATASET

check_obs(adata)
check_var(adata)

adata.write(os.path.join(OUTPUT_DIR, "adata.h5ad"), compression="lzf")


