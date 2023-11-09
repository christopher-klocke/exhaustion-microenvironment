# import statements
import sys
import numpy as np 
import pandas as pd 
import loompy as lp
import scanpy as sc
import anndata as ad
import time
from multiprocessing import Pool, cpu_count
import tqdm
from arboreto.utils import load_tf_names
from arboreto.algo import genie3, grnboost2, _prepare_input
from arboreto.core import SGBM_KWARGS, RF_KWARGS, EARLY_STOP_WINDOW_LENGTH
from arboreto.core import to_tf_matrix, target_gene_indices, infer_partial_network
from pathlib import PurePath


print('arboreto libraries imported')
p_num_workers = int(sys.argv[1]) ## number of workers to use -- rws has 8 cores, use 6 max
p_output_file = sys.argv[2] ## location of the output file
p_tfs_fname = sys.argv[3] ## location of the transcription factor list (mouse, human, etc.)
p_expression_matrix_fname = sys.argv[4] ## location of the input loom expression data file
p_method = 'grnboost2'
#p_cell_id_attribute = 'obs_names' ############ these need to match the way the loom file was written from \ 
#p_gene_attribute = 'var_names' ################## anndata object or will cause a problem 
p_cell_id_attribute = "CellID"
p_gene_attribute = "Gene"

# run the arboreto function calls
# use this option to run GRNBoost2 -- use gradient boost regressor to calculate co-expression modules
method_params = ['GBM', SGBM_KWARGS]  # regressor_type  # regressor_kwargs

# use this option to run GENIE3 -- use random forest regressor to calculate co-expression modules
#method_params = ['RF', RF_KWARGS]  # regressor_type  # regressor_kwargs

# from pyscenic.cli.utils:
ATTRIBUTE_NAME_CELL_IDENTIFIER = "CellID" ## this sets a default value that is not used (keyword argument supplied)
ATTRIBUTE_NAME_GENE = "Gene" ## this sets a default value that is not used (keyword argument supplied)
ATTRIBUTE_NAME_REGULONS_AUC = "RegulonsAUC"
ATTRIBUTE_NAME_REGULONS = "Regulons"
ATTRIBUTE_NAME_METADATA = "MetaData"

def load_exp_matrix_as_loom(
    fname,
    return_sparse=False,
    attribute_name_cell_id: str = ATTRIBUTE_NAME_CELL_IDENTIFIER,
    attribute_name_gene: str = ATTRIBUTE_NAME_GENE,
) -> pd.DataFrame:
    """
    Load expression matrix from loom file.

    :param fname: The name of the loom file to load.
    :return: A 2-dimensional dataframe (rows = cells x columns = genes).
    """
    if return_sparse:
        with lp.connect(fname, mode='r', validate=False) as ds:
            ex_mtx = ds.layers[''].sparse().T.tocsc()
            genes = pd.Series(ds.ra[attribute_name_gene])
            cells = ds.ca[attribute_name_cell_id]
        return ex_mtx, genes, cells

    else:
        with lp.connect(fname, mode='r', validate=False) as ds:
            # The orientation of the loom file is always:
            #   - Columns represent cells or aggregates of cells
            # 	- Rows represent genes
            return pd.DataFrame(
                data=ds[:, :], index=ds.ra[attribute_name_gene], columns=ds.ca[attribute_name_cell_id]
            ).T


def suffixes_to_separator(extension):
    if '.csv' in extension:
        return ','
    if '.tsv' in extension:
        return '\t'


# MODIFYING THIS FUNCTION TO TAKE IN anndata OBJECT DIRECTLY -- LATER ##########################
def load_exp_matrix(
    fname: str,
    #anndata_object, 
    transpose: bool = False,
    return_sparse: bool = False,
    attribute_name_cell_id: str = ATTRIBUTE_NAME_CELL_IDENTIFIER,
    attribute_name_gene: str = ATTRIBUTE_NAME_GENE,
) -> pd.DataFrame:
    """
    Load expression matrix from disk.

    Supported file formats are CSV, TSV and LOOM.

    :param fname: The name of the file that contains the expression matrix.
    :param transpose: Is the expression matrix stored as (rows = genes x columns = cells)?
    :param return_sparse: Returns a sparse matrix when loading from loom
    :return: A 2-dimensional dataframe (rows = cells x columns = genes).
    """

    extension = PurePath(fname).suffixes
    #if is_valid_suffix(extension, 'grn'):
    if True:  ## USING THIS AS A STOPGAP FOR NOW, FIX LATER ################################ 
        
        if '.loom' in extension:
            return load_exp_matrix_as_loom(fname, return_sparse, attribute_name_cell_id, attribute_name_gene)
        elif '.h5ad' in extension:
            from anndata import read_h5ad

            adata = read_h5ad(filename=fname, backed='r')
            #adata = anndata_object

            if return_sparse:
                # expr, gene, cell:
                return adata.X.value, adata.var_names.values, adata.obs_names.values
            else:
                return pd.DataFrame(
                    adata.X.value.todense(), index=adata.obs_names.values, columns=adata.var_names.values
                )

        else:
            df = pd.read_csv(fname, sep=suffixes_to_separator(extension), header=0, index_col=0)
            return df.T if transpose else df
    else:
        raise ValueError("Unknown file format \"{}\".".format(fname))

# TROUBLESHOOT THIS
# I manually changed this to True, make sure this is correct and doesn't
# change anything negatively
#p_sparse=False
p_sparse=True
p_seed = 1
p_transpose1 = 'yes' ## this may cause problems, troubleshoot it
print('arboreto params defined')


def run_infer_partial_network(target_gene_index):
    target_gene_name = gene_names[target_gene_index]
    target_gene_expression = ex_matrix[:, target_gene_index]

    n = infer_partial_network(
        regressor_type=method_params[0],
        regressor_kwargs=method_params[1],
        tf_matrix=tf_matrix,
        tf_matrix_gene_names=tf_matrix_gene_names,
        target_gene_name=target_gene_name,
        target_gene_expression=target_gene_expression,
        include_meta=False,
        early_stop_window_length=EARLY_STOP_WINDOW_LENGTH,
        seed=p_seed,
    )
    return n

if __name__ == '__main__':
    start_time = time.time()
    ex_matrix = load_exp_matrix( ############## CHANGING THIS TO FEED IN ADATA OBJECT DIRECTLY
        p_expression_matrix_fname, (p_transpose1 == 'yes'), p_sparse, p_cell_id_attribute, p_gene_attribute
        #adata, (p_transpose1 == 'yes'), p_sparse, p_cell_id_attribute, p_gene_attribute
    )
    if p_sparse:
        gene_names = ex_matrix[1]
        ex_matrix = ex_matrix[0]
    else:
        gene_names = ex_matrix.columns
    end_time = time.time()
    print(
        f'Loaded expression matrix of {ex_matrix.shape[0]} cells and {ex_matrix.shape[1]} genes in {end_time - start_time} seconds...',
        file=sys.stdout, ############## CHECK this sys.stdout part...
    )
    tf_names = load_tf_names(p_tfs_fname)
    print(f'Loaded {len(tf_names)} TFs...', file=sys.stdout)
    ex_matrix, gene_names, tf_names = _prepare_input(ex_matrix, gene_names, tf_names)
    tf_matrix, tf_matrix_gene_names = to_tf_matrix(ex_matrix, gene_names, tf_names)
    print(f'starting {p_method} using {p_num_workers} processes...', file=sys.stdout)
    start_time = time.time()
    with Pool(p_num_workers) as p:
        adjs = list(
            tqdm.tqdm(
                p.imap(run_infer_partial_network, target_gene_indices(gene_names, target_genes='all'), chunksize=1),
                total=len(gene_names),
            )
        )
    adj = pd.concat(adjs).sort_values(by='importance', ascending=False)
    end_time = time.time()
    print(f'Done in {end_time - start_time} seconds.', file=sys.stdout)
    extension = PurePath(p_output_file).suffixes
    adj.to_csv(p_output_file, index=False, sep=suffixes_to_separator(extension))

print('End of arboreto. Starting ctx...')
