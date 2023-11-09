"""run with 'sc2' conda environment """

# import statements
import logging
import json
import numpy as np
import scanpy as sc
import pandas as pd
from anndata import AnnData
from typing import Union
import scipy
import scipy.stats as stats
from statsmodels.stats.multitest import multipletests
import time

# function definitions
def adata_split(
        adata: AnnData,
        sample_id: str,
        sample_high: Union[str, list],
        sample_low: Union[str, list]
):
    """
    split main AnnData object into separate AnnData objects for high- and low-scoring sample(s)
    """
    adata_high = None
    adata_low = None
    if type(sample_high) is str:
        adata_high = adata[adata.obs[sample_id] == sample_high]
    elif type(sample_high) is list:
        adata_high = adata[adata.obs[sample_id].isin(sample_high)]
    else:
        logging.warning("sample_high must be type 'str' or 'list'")

    if type(sample_low) is str:
        adata_low = adata[adata.obs[sample_id] == sample_low]
    elif type(sample_low) is list:
        adata_low = adata[adata.obs[sample_id].isin(sample_low)]
    else:
        logging.warning("sample_high must be type 'str' or 'list'")

    return adata_high, adata_low

def setup(config_file: str):
    """

    """
    config = json.load(config_file)
    full_adata_filename = config['files']['full_adata_filename']
    out_dir = config['files']['out_dir']
    sample_high = config['samples']['high']
    sample_low = config['samples']['low']
    sample_id = config['parameters']['sample_id']
    celltypes = config['parameters']['celltypes_compare']
    mode = config['parameters']['mode']
    cell_threshold = int(config['parameters']['cell_threshold'])
    adata = sc.read_h5ad(full_adata_filename)
    adata_high, adata_low = adata_split(adata=adata, sample_id=sample_id, sample_high=sample_high, sample_low=sample_low)

    return adata_high, adata_low, out_dir, celltypes, mode, cell_threshold

def cell_count_checks(
        adata_high: AnnData,
        adata_low: AnnData
):
    """

    """
    pass

def run_mw_test(elem, vecA, vecB):
    """
    run mann-whitney u test
    """
    mw_test = stats.mannwhitneyu(vecA, vecB, alternative='two-sided')
    t = mw_test[0]
    p = mw_test[1]
    direction = None
    meanA = np.mean(vecA)
    meanB = np.mean(vecB)

    if meanA == 0 and meanB == 0:
        direction = 'skip'
        logging.info('mean failed for ' + elem + ' -- mean expression is zero for both sample groups')
    elif meanA > meanB:
        direction = 'low'
    elif meanA < meanB:
        direction = 'high'
    elif meanA == meanB:
        logging.info('mean failed for ' + elem + ' -- means are non-zero but identical')
    else:
        logging.warning('mean failed for ' + elem + ' -- check data types')  # different message

    return t, p, direction

def calc_scores(
        adata_high: AnnData,
        adata_low: AnnData,
        celltypes: list,
        mode: str,
        cell_threshold: int = 10
):
    """

    """
    celltype_dfs = []
    for celltype in celltypes:
        print('calculating scores for celltype: ' + celltype)
        adata_high_cut = adata_high[adata_high.obs['celltype'] == celltype]
        adata_low_cut = adata_low[adata_low.obs['celltype'] == celltype]

        ## this belongs in cell count function
        #logging.info('low sample group has ' + str(adata_low_cut.shape[0]) + ' cells of type: ' + celltype)
        #logging.info('high sample group has ' + str(adata_high_cut.shape[0]) + ' cells of type: ' + celltype)

        if adata_low_cut.shape[0] < cell_threshold or adata_high_cut.shape[0] < cell_threshold:
            logging.warning(
                'insufficient cells in these samples to analyze celltype: ' + celltype + '; skip cell type or pick another sample group')
        else:
            if mode == 'tf':
                res_add = []
                ## optimize as with gene section below
                tfs = adata_low_cut.obsm['X_regulonsAUC'].columns
                low_df = adata_low_cut.obsm['X_regulonsAUC']
                high_df = adata_high_cut.obsm['X_regulonsAUC']
                for tf_n in range(len(tfs)):
                    tf_name = tfs[tf_n]
                    vecA = low_df[tf_name]
                    vecB = high_df[tf_name]
                    t, p, direction = run_mw_test(tf_name, vecA, vecB)
                    res = (tf_name[:-3], t, p, direction)  # SETUP FOR TF COMPARISON
                    res_add.append(res)
                df = pd.DataFrame(res_add, columns=['tf', 'test_statistic', 'p_value', 'direction'])
                new_col = [celltype] * df.shape[0]
                df['celltype'] = new_col
                celltype_dfs.append(df)

            elif mode == 'gene':
                res_add = []
                g_names = adata_low_cut.var.index
                low_x = adata_low_cut.X.toarray()  # '.toarray()' added to fix array format issue
                high_x = adata_high_cut.X.toarray()
                for g in range(adata_low_cut.shape[1]):
                    g_name = g_names[g]
                    vecA = low_x[:, g]
                    vecB = high_x[:, g]
                    t, p, direction = run_mw_test(g_name, vecA, vecB)
                    res = (g_name, t, p, direction)
                    res_add.append(res)
                df = pd.DataFrame(res_add, columns=['gene', 'test_statistic', 'p_value', 'direction'])
                new_col = [celltype] * df.shape[0]
                df['celltype'] = new_col
                celltype_dfs.append(df)

    else:
        logging.warning("mode parameter must be either 'tf' or 'gene'")

    results_df = pd.concat(celltype_dfs)

    return results_df

def pval_correct(
        results: pd.DataFrame,
        method: str = 'fdr_bh',
        alpha: float = 0.05,
        sort: bool = False,
        only_sig: bool = True
):
    """

    """
    correction_input = results['p_value'].to_numpy()
    reject, pvals_corrected, alphacSidak, alphacBonf = multipletests(correction_input, alpha=alpha, method=method)
    results['adjusted_p_value'] = pvals_corrected
    if only_sig:
        results['reject'] = reject
    if sort:
        results.sort_values(by='adjusted_p_value', inplace=True)

    return results

def analyze(config_file: str):
    """

    """
    adata_high, adata_low, out_dir, celltypes, mode, cell_threshold = setup(config_file)
    #cell_count_checks(adata_high, adata_low)
    results = calc_scores(adata_high=adata_high, adata_low=adata_low, celltypes=celltypes, mode=mode, cell_threshold=cell_threshold)
    results = pval_correct(results)
    results.sort_values(by='adjusted_p_value', inplace=True)
    print(results)
    results.to_csv(out_dir + 'combined_bh.csv')

# main function definition
def main():
    tic = time.perf_counter()


    toc = time.perf_counter()
    timer = toc - tic
    #display_results(ligand_scores=scores, runtime=timer)
    #save_results()

# run main function
if __name__ == "__main__":
    main()

