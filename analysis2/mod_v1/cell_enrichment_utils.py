import os
import scanpy as sc
import pandas as pd 
import numpy as np
from anndata import AnnData
from typing import Dict, Optional, Union
from itertools import chain
import warnings
import logging
from scipy.stats import spearmanr
import matplotlib.pyplot as plt 
from statsmodels.stats.multitest import multipletests
import pathway_analyzer as pa


def create_cell_gmt(
    adata: AnnData, outfile : str, sample_id_name: str = 'patient'): 
    """ This function takes an AnnData object and an output filename as input. 
    It creates a dictionary from the AnnData object, with sample IDs as keys and
    a list of their corresponding cell IDs as values. It then creates a '.gmt-like'
    file (with the .gmt suffix), writing each sample ID and its corresponding cell
    IDs to its own line, separated by tabs. 

    Params
    ------
    adata
        Anndata object
    outfile
        Location to write output .gmt file 
    sample_id_name
        sample ID column label in adata.obs (typically 'patient' or 'batch')

    Returns
    -------
    No object returned. Output written to 'outfile' in '.gmt' format. 

    Example
    --------
    >>> create_cell_gmt() ## UPDATE
    """

    if outfile[-4:]=='.gmt': 
        if os.path.exists(outfile):
            os.remove(outfile)

        sample_dict = {}
        cell_id = adata.obs.index
        sample_id = adata.obs[sample_id_name]
        
        for i in range(len(cell_id)): 
            if sample_id[i] in sample_dict.keys(): 
                sample_dict[sample_id[i]].append(cell_id[i])
            else: 
                sample_dict[sample_id[i]] = [cell_id[i]]

        with open(outfile, 'a') as write_out: 
            for i in sample_dict.keys(): 
                add_line = i + '\t' + '\t'.join(j for j in sample_dict[i]) + '\n'
                write_out.write(add_line)

    else: 
        logging.error('Output file must use ".gmt" suffix')


def create_cell_adata(
    adata: AnnData, 
    mode: str, 
    pt_label: str = 'monocle3_pseudotime', 
    gene_set_label: str = 'X_regulonsAUC', 
    celltype_name: str = 'celltype'): 
    """Returns modified adata that can be used as input into 'GSEAPyWrapper()'
    for cell enrichment analysis. 

    Params
    ------
    adata
        Anndata object -- in 'pseudotime' mode, contains pseudotime values for
        each cell in 'adata.obs' column; in 'gene_sets' mode, contains AUCell 
        activity scores for each cell in Pandas dataframe stored in adata.obsm
    mode
        Mode to run in -- either 'pseudotime' or 'gene_sets'
    pt_label
        Pseudotime column label in adata.obs
    gene_set_label
        Label of Pandas dataframe in adata.obsm containing AUCell activity
        scores
    celltype_name
        fill in 

    Returns
    -------
    Anndata object -- not a normal one though -- describe ## 

    """

    if mode == 'pseudotime': 
        pt = adata.obs[pt_label]
        adata_mod = sc.AnnData(pd.DataFrame(pt)).T
     
    elif mode == 'gene_sets': 
        df = adata.obsm[gene_set_label]
        adata_mod = sc.AnnData(df.T)

        assert celltype_name in adata.obs.keys(), "cell type field in input AnnData object is missing or labeled incorrectly"

        ## add 'celltypes' from adata.obs ##########
        adata_mod.var[celltype_name] = adata.obs[celltype_name]

    else: 
        print('Not an accepted mode. Set "mode" to "pseudotime" or "gene_sets"')

    return adata_mod


class CSEAPyWrapper(pa.PathwayAnalyzer):
    """
    adapted from 'GSEAPyWrapper(PathwayAnalyzer)' in pathway_analyzer (from Reactome code) 

    A wrapper of CSEAPy -- Cell Set Enrichment Analysis in Python 

    adding this here to avoid modifying original 'pathway_analyzer' script -- just import 
    needed functions / classes and add new functionality 

    needed to add this to pass through non-default arguments for 'min_size' and 'max_size' 
    """
    def __init__(self):
        super().__init__()
        self.adata_key = pa.SSGSEA_KEY

        ## add min, max arguments -- also key? 

    def csea(self,
               adata: Union[str, AnnData],
               reactome_gmt: Union[str, dict],
               data_key: str,
               weigted_score_type: float = 0.25,
               n_processes: int = 1
               #min_size: int = 1,
               #max_size: int = 100000
               ) -> Optional[pd.DataFrame]:
        """
        Perform ssgsea analysis
        :param adata: an AnnData object or a hd5 file
        :param reactome_gmt:
        :param data_key
        :param weigted_score_type
        :param n_processes
        :return:
        """
        if data_key is not None:
            self.adata_key = data_key
        adata, adata_df = pa._load_data(adata)
        genesets_dict = pa._load_reactome_gmt(reactome_gmt)

        ss = pa.ReactomeSSGSEA(adata_df,
                            genesets_dict,
                            weighted_score_type=weigted_score_type,
                            scale=False,
                            processes=n_processes)

        ss.min_size=1 ## to fix filtering issue -- set manually 
        ss.max_size=100000
        ss.run()
        # ss = ssgsea_runner.ssgsea(adata_df,
        #                genesets_dict,
        #                weighted_score_type=weigted_score_type,
        #                scale=False,
        #                processes=n_processes,
        #                outdir=None)  # Don't want to keep the output directory
        # Keep the raw enrichment scores
        if isinstance(adata, AnnData):
            adata.obsm[self.adata_key] = pd.DataFrame(ss.resultsOnSamples).T  # Need to revert it back
            self.adata = adata
        else:
            return pd.DataFrame(ss.resultsOnSamples)


def add_celltypes(
    adata: AnnData, 
    celltype_dict: dict, 
    celltype_name: str = 'celltype',
    clusters_name: str = 'leiden',
    category_convert: bool = False
): 
    """
    
    given adata with leiden clusters but no celltype labels, and 
    dict mapping celltypes to leiden cluster IDs
    
    add celltypes to adata in adata.obs['celltype']

    celltype_dict ## dict with celltype string as key and list of leiden as value 

    """

    ## build lookup dictionary from input dictionary 

    values_list_structured = list(celltype_dict.values())
    values_list_flat = [str(i) for i in list(map(int, chain.from_iterable(values_list_structured)))]

    assert len(set(values_list_flat)) == len(values_list_flat), "Ensure no cluster ID is assigned to more than one cell type"

    lookup_dict = {}

    for cluster_id in values_list_flat: 
        for celltype in celltype_dict.keys(): 
            if cluster_id in celltype_dict[celltype]: 
                lookup_dict[cluster_id] = celltype
                #break ####### check this -- think about best way to break out of loop once found -- avoid unnecessary search ############
        if cluster_id not in lookup_dict.keys(): 
            warnings.warn("WARNING: cluster ID was not found in cell type dictionary")
        #assert cluster_id in lookup_dict.keys(), "cluster ID was not found in cell type dictionary"

    adata.obs[celltype_name] = [lookup_dict[i] for i in list(adata.obs[clusters_name])]

    if category_convert: 
        adata.obs[celltype_name] = adata.obs[celltype_name].astype("category")

    return adata 


def create_enrichment_dict(
    adata_mod_pt: str, 
    gene_sets_adata_list: str
) -> dict: 
    """Creates enrichment dictionary. 

    Params
    ------
    adata_mod_pt
        Anndata object -- describe ##

    gene_sets_adata_list
        describe ##

    Returns
    -------
    dictionary
    """
    adata_mod_pseudotime = sc.read_h5ad(adata_mod_pt)

    cell_enrichment_dict = {
        'pseudotime': adata_mod_pseudotime.obsm[pa.SSGSEA_KEY], 
        'gene_sets': {}
        }

    with open(gene_sets_adata_list, 'r') as input1: 
        for line in input1: 
            line = line.rstrip() 
            adata_celltype_mod = sc.read_h5ad(line)
            line = line.split('/')
            cellname = line[-1][25:-5] #######
            cell_enrichment_dict['gene_sets'][cellname] = adata_celltype_mod.obsm[pa.SSGSEA_KEY]

    return cell_enrichment_dict


def cell_enrichment_correlation(
    enrich_dict: dict, 
    permute: bool = False, 
    n_perm: int = 1000
): 
    """Run correlation 

    input is output dictionary from 'create_enrichment_dict()' function 

    convert to object method? 

    pulling from 'run_score_correlation()' from 'utils2.py' 

    Params
    ------
    enrich_dict

    permute

    n_perm


    Returns
    -------
    Correlation dictionary -- used as input into 
    """

    #enrich_dict['pseudotime'] ## pandas df with sample IDs and sample cell enrichment cd8 exhaustion scores

    #enrich_dict['gene_sets'] ## dictionary of pandas dfs keyed by cell type -- pandas df with sample IDs and sample cell enrichment aucell activity scores 

    correlation_dict = {}

    for cell_type in enrich_dict['gene_sets'].keys(): 

        act_df = enrich_dict['gene_sets'][cell_type]

        pt_df = enrich_dict['pseudotime']

        if len(pt_df.columns) != len(act_df.index): 

            print("Sample IDs do not match")

            ## remove any sample IDs not present in pseudotime ranking 
            samples_remove_act = [i for i in list(act_df.columns) if i not in pt_df.columns] ## can do at higher level -- running step repeatedly #######
                
            print("Samples removed from activity score dataframe: " + str(samples_remove_act)) ###

            for i in samples_remove_act: 

                act_df = act_df.drop(i, axis=1)

            samples_remove_pt = [i for i in list(pt_df.columns) if i not in act_df.columns]

            print("Samples removed from pseudotime dataframe for cell type " + cell_type + ": " + str(samples_remove_pt))

            for i in samples_remove_pt: 

                pt_df = pt_df.drop(i, axis=1)

        #assert list(pt_df.columns) == list(act_df.columns), "Check matching of sample IDs" ######

        for i in list(act_df.index): ## iterate through gene sets

            rho, p = spearmanr(pt_df.T.to_numpy(), act_df.loc[i, :].to_numpy()) ## calculate Spearman Rank correlation and corresponding p-value

            if permute: 
                print('running permutations (this part may take a while)...')
                rho_values = []
                for n in range(n_perm): 
                    shuffled = act_df.loc[i, :].to_numpy()
                    np.random.shuffle(shuffled)
                    rho_x, p_x = spearmanr(pt_df.T.to_numpy(), shuffled)
                    rho_values.append(rho_x)

                rho_values = np.asarray(rho_values)
                if rho > 0: 
                    p_perm = sum(rho_values >= rho) / n_perm

                elif rho < 0: 
                    p_perm = sum(rho_values <= rho) / n_perm
                
                else: 
                    logging.warning('rho is exactly zero')
                    p_perm = 0.5 ######## FIX #########

                mean_rho = np.mean(rho_values)

                ## add results to dictionary
                if cell_type in correlation_dict.keys(): 
                    correlation_dict[cell_type][i] = (rho, p, p_perm, mean_rho) ## add to results dictionary
                else: 
                    correlation_dict[cell_type] = {}
                    correlation_dict[cell_type][i] = (rho, p, p_perm, mean_rho) ## add to results dictionary

                print('permutation complete -- celltype: ' + cell_type + ' gene set: ' + i)

            else: 

                ## add results to dictionary
                if cell_type in correlation_dict.keys(): 
                    correlation_dict[cell_type][i] = (rho, p) ## add to results dictionary
                else: 
                    correlation_dict[cell_type] = {}
                    correlation_dict[cell_type][i] = (rho, p) ## add to results dictionary

    return correlation_dict


def add_fdr_perm(
    res_dict: Optional[dict] = None, 
    res_df: Optional[pd.DataFrame] = None, 
    df_in: bool = False, 
    alpha: float = 0.05
) -> pd.DataFrame: 

    """Creates enrichment dictionary. 

    can pass in dictionary (if running all at once) or df (read from file if 
    run beforehand)
        
    correlation_dict[cell_type][i] = (rho, p, p_perm, mean_rho)

    Params
    ------
    res_dict

    res_df

    df_in

    alpha

    Returns
    -------
    dataframe
    """

    ## create dataframe 

    if df_in: ## inputting df directly -- read from csv file 
        if res_df is None: 
            logging.error('Must supply dataframe if "df_in == True"')
            ## raise error here? 
        else: 
            df = res_df

    else: 

        col_names = ['celltype', 'gene_set', 'rho', 'p_value', 'p_perm', 'mean_rho']
        df = pd.DataFrame(columns=col_names)

        for celltype in res_dict.keys(): 
            for gene_set in res_dict[celltype].keys(): 
                res = res_dict[celltype][gene_set]
                #df = df.append({'celltype': celltype, 'gene_set': gene_set, 'rho': res[0], 'p_value': res[1]}, ignore_index=True)

                df_add = pd.DataFrame([[celltype, gene_set, res[0], res[1], res[2], res[3]]], columns=['celltype', 'gene_set', 'rho', 'p_value', 'p_perm', 'mean_rho'])
                df = pd.concat([df, df_add])

    ## run p-value correction on dataframe column, create adjusted p-value column 

    correction_input = df['p_perm'].to_numpy()

    reject, pvals_corrected, alphacSidak, alphacBonf = multipletests(correction_input, alpha=alpha, method='fdr_bh')

    df['adjusted_p_value'] = pvals_corrected
    df['reject'] = reject

    return df


def export_results_csv(
    out_dir: str, 
    scores_dict: Optional[dict] = None, 
    scores_df: Optional[pd.DataFrame] = None, 
    fdr: bool = False, 
    permute: bool = False
): 
    """ Export results to .csv files. 

    Input file is output dictionary from 'cell_enrichment_correlation()' 

    out_dir is string desribing path to a directory -- output .csv files 
    written here 

    dict[celltype][gene_set] = (rho, p)

    correlation_dict[cell_type][i] = (rho, p, p_perm, mean_rho)

    ['celltype', 'gene_set', 'rho', 'p_value', 'p_perm', 'mean_rho', 'adjusted_p_value', 'reject']
    """

    if fdr and permute: 
        if scores_df is None: 
            logging.error('scores_df argument not supplied')
        else: 
            outfile = out_dir + 'results_permute_fdr.csv'
            with open(outfile, 'w') as output_write: 
                output_write.write('celltype,gene_set,rho,p_value,p_perm,mean_rho,adjusted_p_value,reject' + '\n')
                for i in range(scores_df.shape[0]): 
                    row = scores_df.iloc[i]
                    line = row['celltype'] + ',' + row['gene_set'] + ',' + str(row['rho']) + ',' + str(row['p_value']) + ',' + str(row['p_perm']) + ',' + str(row['mean_rho']) + ',' + str(row['adjusted_p_value']) + ',' + str(row['reject'])
                    output_write.write(line + '\n')

    elif fdr: 
        pass ## add functionality 

    elif permute: 
        if scores_dict is None: 
            logging.error('scores_dict argument not supplied')
        else: 
            outfile = out_dir + 'results_permute.csv'
            with open(outfile, 'w') as output_write: 
                output_write.write('celltype,gene_set,rho,p_value,p_perm,mean_rho' + '\n')
                for celltype in scores_dict.keys(): 
                    for gene_set in scores_dict[celltype].keys(): 
                        res = scores_dict[celltype][gene_set]
                        output_write.write(celltype + ',' + gene_set + ',' + str(res[0]) + ',' + str(res[1]) + ',' + str(res[2]) + ',' + str(res[3]) + '\n')

    else: 
        if scores_dict is None: 
            logging.error('scores_dict argument not supplied')
        else: 
            for celltype in scores_dict.keys(): 
                outfile = out_dir + celltype + '.csv'
                with open(outfile, 'w') as output_write: 
                    output_write.write('gene set,rho,p-value' + '\n')
                    for gene_set in scores_dict[celltype].keys(): 
                        res = scores_dict[celltype][gene_set]
                        output_write.write(gene_set + ',' + str(res[0]) + ',' + str(res[1]) + '\n')


def csea_wrapper(
    adata_pt_filename: str, 
    adata_act_filename: str, 
    results_dir: str, 
    celltype_dict: dict, 
    remove_inf: bool = True, 
    pt_name: str = 'monocle3_pseudotime', 
    inf_name: str = 'inf', 
    set_n_perm: int = 1000, 
    alpha: float = 0.05, 
    remove_pt_clusters : Optional[list] = None, 
    sample_id: str = 'patient',
):
    """ wrapper script to replace all of the function calls in 
    'cell_enrichment_li1.py' 
    FINISH DESCRIPTION HERE

    takes as input strings to two Anndata objects -- one for pseudotime and 
    one for AUCell scores

    Params
    ------
    adata_pt_filename
        Anndata object containing pseudotime values calculated with Monocle3. 
        Created with (list script) 
    adata_act_filename
        Anndata object containing...
    results_dir
        Output directory -- intermediate files and output files stored here
    remove_inf

    pt_name

    inf_name

    set_n_perm

    alpha

    remove_pt_clusters

    sample_id

    Returns
    -------
    Returns dictionary with...

    Example
    --------
    >>> import cell_enrichment_utils as ceu
    >>> wang_path = '/Users/klockec/Documents/data/wang_data/'
    >>> wang_act_filename = 'wang_dd_viz_ready_regs_cut.h5ad'
    >>> wang_pt_input_filename = 'wang_t_monocle1.h5ad'
    >>> results_dir = '/Users/klockec/Documents/data/analysis_files/p3/mod1/wang_run2/'
    >>> cell_types_input = {'macrophage': ['4', '8'], 'plasma': ['7'], 
    'B': ['1'], 'CD8_T': ['0', '2', '3', '5', '6']}
    >>> ceu.csea_wrapper(adata_pt_filename=(wang_path + wang_pt_input_filename), adata_act_filename=(wang_path + wang_filename), results_dir=results_dir, celltype_dict=cell_types_input, sample_id='batch', set_n_perm=10, alpha=0.2)
    """
    logging.info('loading the data...')

    if os.path.isfile(adata_pt_filename): 
        adata_pt = sc.read_h5ad(adata_pt_filename) ## anndata object with monocle pseudotime results
    else: 
        logging.error('valid pseudotime anndata not provided')

    if os.path.isfile(adata_act_filename): 
        adata_act = sc.read_h5ad(adata_act_filename) ## anndata object with AUCell activity scores 
    else: 
        logging.error('valid gene set anndata not provided')

    ## remove 'inf' pseudotime values 
    if remove_inf: 
        adata_pt = adata_pt[adata_pt.obs[pt_name] != float(inf_name)]

    if remove_pt_clusters is not None: ## remove CD4 clusters
        for i in remove_pt_clusters: 
            adata_pt = adata_pt[adata_pt.obs['leiden'] != i] ## more efficient / pythonic way to do this -- list comparison -- not if/else, use list comprehension or R-like vector comparison 

    logging.info('running in pseudotime mode')

    pt_cell_gmt_filename = results_dir + 'pt_cell_gmt.gmt'
    create_cell_gmt(adata=adata_pt, outfile=pt_cell_gmt_filename, sample_id_name=sample_id)
    adata_mod_pseudotime = create_cell_adata(adata=adata_pt, mode='pseudotime')

    cell_enrichment_wrapper_pt = CSEAPyWrapper()
    cell_enrichment_wrapper_pt.csea(adata_mod_pseudotime, reactome_gmt=pt_cell_gmt_filename, data_key=pa.SSGSEA_KEY, weigted_score_type=0, n_processes=1) ## typo in 'weighted' #####
    adata_mod_pseudotime.write_h5ad(results_dir + 'cell_enrichment_pt.h5ad') ## write to pseudotime output .h5ad file 

    logging.info('pseudotime section complete')
    logging.info('running in gene set mode')

    add_celltypes(adata=adata_act, celltype_dict=celltype_dict)

    ## need to iterate through, selecting one cell type at a time
    ## need cell-type-specific cell set .gmt files 

    adata_mod_regulons = create_cell_adata(adata=adata_act, mode='gene_sets')
    cell_types_list = list(set(list(adata_mod_regulons.var['celltype'])))

    celltypes_filename = results_dir + 'celltype_files.txt'
    if os.path.exists(celltypes_filename):
        os.remove(celltypes_filename)

    for celltype in cell_types_list: 
        adata_celltype_original = adata_act[adata_act.obs['celltype'] == celltype]
        celltype_gmt_filename = results_dir + celltype + '_celltype_gmt.gmt'
        logging.info(celltype_gmt_filename) ##########
        create_cell_gmt(adata=adata_celltype_original, outfile=celltype_gmt_filename, sample_id_name=sample_id)
        adata_celltype_mod = adata_mod_regulons[:, adata_mod_regulons.var['celltype'] == celltype]
        cell_enrichment_wrapper_regulons = CSEAPyWrapper()
        cell_enrichment_wrapper_regulons.csea(adata_celltype_mod, reactome_gmt=celltype_gmt_filename, data_key=pa.SSGSEA_KEY, weigted_score_type=0, n_processes=1) ## typo in 'weighted' #####
        outfile = results_dir + 'cell_enrichment_activity_' + celltype + '.h5ad' ##### 
        adata_celltype_mod.write_h5ad(outfile)

        with open(results_dir + 'celltype_files.txt', 'a') as output1: 
            output1.write(outfile + '\n')

    logging.info('gene set section complete')

    ## create enrichment dictionary 
    cell_enrichment_dict = create_enrichment_dict(
        adata_mod_pt=(results_dir + 'cell_enrichment_pt.h5ad'), 
        gene_sets_adata_list=(results_dir + 'celltype_files.txt'))

    ## remove PBMC samples
    samples_remove_pt = [i for i in list(cell_enrichment_dict['pseudotime'].columns) if i[-4:] == 'PBMC']

    for i in samples_remove_pt: 

        cell_enrichment_dict['pseudotime'] = cell_enrichment_dict['pseudotime'].drop(i, axis=1)

    logging.info('calculating correlations...')

    results_dict = cell_enrichment_correlation(
        enrich_dict = cell_enrichment_dict, permute=True, 
        n_perm=set_n_perm)

    logging.info('writing in-process results to file...')

    export_results_csv(scores_dict=results_dict, out_dir=results_dir, 
        permute=True)

    outfile = results_dir + 'results_permute.csv'

    df = pd.read_csv(outfile) 

    results_df = add_fdr_perm(res_df=df, df_in=True, alpha=alpha)

    logging.info('writing final results to file...')
    
    export_results_csv(scores_df=results_df, out_dir=results_dir, 
        permute=True, fdr=True)

    logging.info('CSEA complete. ')




























def get_significant_results(scores_dict, alpha=0.05):
    """
    print most significant results -- all with p-values below chosen threshold 
    
    """
    for celltype in scores_dict.keys(): 
        g_dict = scores_dict[celltype]
        for r in g_dict.keys(): 
            p = g_dict[r][1]
            if p < alpha: ## is p value below chosen alpha? 
                print('cell type: ' + celltype + ' gene set: ' + str(r) + ' p-value: ' + str(p))

def plot_scores(
    enrich_dict: dict, 
    results_dict: dict, 
    celltype, 
    gene_set, 
    save=False, 
    savefile: Optional[str] = None
): 
    """ plotting function for significant results 

    updated version

    enrich_dict['pseudotime'] -- pandas df with pseudotime values 

    enrich_dict['gene_sets'][cell_type] -- pandas df with activity values 

    correlation_dict[cell_type][i] = (rho, p) ## add to results dictionary

    """

    rho = results_dict[celltype][gene_set][0]
    p = results_dict[celltype][gene_set][1]

    logging.info('cell type: ' + celltype + '\n' + 'gene set: ' + gene_set + '\n' + 'rho: ' + str(rho) + '\n' + 'p-value: ' + str(p))

    act_df = enrich_dict['gene_sets'][celltype]
    pt_df = enrich_dict['pseudotime']

    if len(pt_df.columns) != len(act_df.index): 

        print("Sample IDs do not match")

        ## remove any sample IDs not present in pseudotime ranking 
        samples_remove_act = [i for i in list(act_df.columns) if i not in pt_df.columns] ## can do at higher level -- running step repeatedly #######
                
        logging.warning("Samples removed from activity score dataframe: " + str(samples_remove_act)) ###

        for i in samples_remove_act: 
            act_df = act_df.drop(i, axis=1)

        samples_remove_pt = [i for i in list(pt_df.columns) if i not in act_df.columns]

        logging.warning("Samples removed from pseudotime dataframe for cell type " + celltype + ": " + str(samples_remove_pt))

        for i in samples_remove_pt: 
            pt_df = pt_df.drop(i, axis=1)

    plt.clf()
    x = np.array(pt_df)
    y = np.array(act_df.loc[gene_set])

    if save: 
        plt.scatter(x, y) ## need to add s -- counts -- point size 
        plt.xlabel('sample-level exhaustion score')
        plt.ylabel('gene set activity')
        title_text = 'cell type: ' + celltype + ' gene set: ' + gene_set + ' rho: ' + str(rho) + ' p-value: ' + str(p)
        plt.title(title_text)
        plt.savefig(savefile)

    else: 
        plt.scatter(x, y) ## need to add s -- counts -- point size 
        plt.xlabel('sample-level exhaustion score')
        plt.ylabel('gene set activity')
        title_text = 'cell type: ' + celltype + ' gene set: ' + gene_set + ' rho: ' + str(rho) + ' p-value: ' + str(p)
        plt.title(title_text)
        plt.show(block=True)

def add_fdr_correction(res_dict): 
    """ Add FDR-adjusted p-values to results dictionary -- no -- return df -- update 

    input is output dictionary from 'cell_enrichment_correlation()'

    output is same dictionary, with adjusted p values added -- no -- return df

    can be used as input into 'export_results.csv()' 
    
    



    pull out all p-values from dictionary -- in labeled and organized way

    run correction

    put back in 



    can use pandas df -- column to track celltype, column for regulon, column for p values

    can run one correction on p value vector, add adjusted p values as another column 

    then subset out by celltype column (index df), put back into dict by match? 

    not the most efficient, but it will work -- can optimize later -- runtime will be short anyway 

    """

    ## create dataframe 

    col_names = ['celltype', 'gene_set', 'rho', 'p_value']
    df = pd.DataFrame(columns=col_names)

    for celltype in res_dict.keys(): 
        for gene_set in res_dict[celltype].keys(): 
            res = res_dict[celltype][gene_set]
            #df = df.append({'celltype': celltype, 'gene_set': gene_set, 'rho': res[0], 'p_value': res[1]}, ignore_index=True)

            df_add = pd.DataFrame([[celltype, gene_set, res[0], res[1]]], columns=['celltype', 'gene_set', 'rho', 'p_value'])
            df = pd.concat([df, df_add])

    ## run p-value correction on dataframe column, create adjusted p-value column 

    correction_input = df['p_value'].to_numpy()

    reject, pvals_corrected, alphacSidak, alphacBonf = multipletests(correction_input, alpha=0.05, method='fdr_bh') ## CHECK METHOD CHOICE ('fdr_bh' vs. 'fdr_by')

    df['adjusted_p_value'] = pvals_corrected
    df['reject'] = reject

    ## data from df back into dictionary 

    ## 'reject' is boolean vector, true for those than can be rejected 

    ## 'pvals_corrected' is adjusted p-values 

    #multipletests(pvals, alpha=0.05, method='hs', is_sorted=False, returnsorted=False)

    #return res_dict

    return df 

def export_results_csv2(
    scores_df: pd.DataFrame, 
    out_dir: str, 
    fdr: bool = False
): 
    """ Export results to .csv files. 

    Input file is output dictionary from 'cell_enrichment_correlation()' 

    out_dir is string desribing path to a directory -- output .csv files 
    written here 



    dict[celltype][gene_set] = (rho, p)
 
    """

    if fdr: 
        outfile = out_dir + 'results.csv'
        with open(outfile, 'w') as output_write: 
            output_write.write('celltype,gene_set,rho,p_value,adjusted_p_value,reject' + '\n')
            for i in range(scores_df.shape[0]): 
                row = scores_df.iloc[i]
                output_write.write(row['celltype'] + ',' + row['gene_set'] + ',' + str(row['rho']) + ',' + str(row['p_value']) + ',' + str(row['adjusted_p_value']) + ',' + str(row['reject']) + '\n')

    else: 
        pass 


def get_gene_set_results(
    scores_dict: dict, 
    gene_set
):
    """get results for one TF / pathway / gene set for all cell type clusters"""
    for celltype in scores_dict.keys(): 
        g_dict = scores_dict[celltype]
        try: 
            res = g_dict[gene_set]
            print('cell type: ' + celltype + ' gene set: ' + gene_set + ' p-value: ' + str(res[1]) + ' rho: ' + str(res[0]))
        except KeyError: 
            print('key ' + gene_set + ' not in ' + celltype + ' object')

def main(): 
    pass

if __name__ == "__main__": 
    main()

