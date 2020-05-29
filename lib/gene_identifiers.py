from anndata import AnnData
import scipy as sp
import numpy as np
import pandas as pd

MART = pd.read_csv("tables/biomart.tsv", sep="\t")
MART = MART.rename({
    "Gene stable ID": "ensg",
    "Gene name": "gene_name",
    "HGNC symbol": "gene_symbol"}, axis="columns")

def map_to_ensembl(adata):
    """
    map gene names to ensembl using biomart.
    """
    var_df = pd.DataFrame().assign(var_names = adata.var_names)
    var_df_ensg = var_df.join(MART[["ensg", "gene_name"]].set_index("gene_name"), on="var_names", how="left")

    var_df_ensg_nona = var_df_ensg.dropna()
    var_df_ensg_nona_nodups = var_df_ensg_nona.drop_duplicates(subset=["var_names"], keep=False)

    print("[remapping] dropped {} genes because no matching ENSG was found.".format(
                len(var_df_ensg) - len(var_df_ensg_nona)))
    print("[remapping] dropped {} genes because they mapped to multiple ENSG.".format(
                len(set(var_df_ensg_nona.loc[var_df_ensg_nona.duplicated("var_names"), "var_names"].values))))

    adata = adata[:, var_df_ensg_nona_nodups.index.values]

    adata.var["gene_symbols"] = adata.var_names
    adata.var_names = var_df_ensg_nona_nodups["ensg"]

    return adata


def collapse_gene_symbols(adata, fun=np.max):
    """
    make gene symbols unique. Aggregate multiple
    gene symbols of the same name by `fun`.

    Will densify the matrix, although it could
    also be implemented in a sparse way.
    """
    # work on numpy arrays. This should also work on sparse matrices.
    new_var_names = np.unique(adata.var_names)
    X2 = np.zeros((adata.shape[0], len(new_var_names)))
    for i, sym in enumerate(new_var_names):
        X2[:, i] = fun(adata.X[:, adata.var_names == sym], axis=0)

    new_var = pd.DataFrame().assign(gene_symbols = new_var_names).set_index("gene_symbols", inplace=False)
    return AnnData(X2, obs = adata.obs, var=new_var)






