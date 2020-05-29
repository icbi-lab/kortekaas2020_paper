"""
Patched versions of SCENIC code copied over
from their pySCENIC repository. 
"""
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
from math import ceil, floor, sqrt
from scipy.spatial.distance import jensenshannon

def plot_rss(rss, cell_type, top_n=5, max_n=None, ax=None, fontdict=None):
    if ax is None:
        _, ax = plt.subplots(1, 1, figsize=(4, 4))
    if max_n is None:
        max_n = rss.shape[1]
    data = rss.T[cell_type].sort_values(ascending=False)[0:max_n]
    ax.plot(np.arange(len(data)), data, '.')
    ax.set_ylim([floor(data.min()*100.0)/100.0, ceil(data.max()*100.0)/100.0])
    ax.set_ylabel('RSS')
    ax.set_xlabel('Regulon')
    ax.set_title(cell_type)
    ax.set_xticklabels([])

    font = {
        'color':  'red',
        'weight': 'normal',
        'size': 6,
    }
    if fontdict is not None:
      font.update(fontdict)

    for idx, (regulon_name, rss_val) in enumerate(zip(data[0:top_n].index, data[0:top_n].values)):
        ax.plot([idx, idx], [rss_val, rss_val], 'r.')
        ax.text(idx+(max_n/25), rss_val, regulon_name, fontdict=font, horizontalalignment='left', verticalalignment='center')


def regulon_specificity_scores(auc_mtx, cell_type_series):
    """
    Calculates the Regulon Specificty Scores (RSS). [doi: 10.1016/j.celrep.2018.10.045]
    :param auc_mtx: The dataframe with the AUC values for all cells and regulons (n_cells x n_regulons).
    :param cell_type_series: A pandas Series object with cell identifiers as index and cell type labels as values.
    :return: A pandas dataframe with the RSS values (cell type x regulon).
    """

    cell_types = list(cell_type_series.unique())
    n_types = len(cell_types)
    regulons = list(auc_mtx.columns)
    n_regulons = len(regulons)
    rss_values = np.empty(shape=(n_types, n_regulons), dtype=np.float)

    def rss(aucs, labels):
        # jensenshannon function provides distance which is the sqrt of the JS divergence.
        return 1.0 - jensenshannon(aucs/aucs.sum(), labels/labels.sum())

    for cidx, regulon_name in enumerate(regulons):
        for ridx, cell_type in enumerate(cell_types):
            rss_values[ridx, cidx] = rss(auc_mtx[regulon_name], (cell_type_series == cell_type).astype(int))

    return pd.DataFrame(data=rss_values, index=cell_types, columns=regulons)
