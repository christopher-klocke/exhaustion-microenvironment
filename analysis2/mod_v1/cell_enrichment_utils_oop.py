"""adapted from 'cell_enrichment_utils.py'"""

import os
from statistics import correlation
from unittest import result
import scanpy as sc
import pandas as pd 
import numpy as np
from anndata import AnnData
from typing import Dict, Optional, Union
from itertools import chain
import logging
from scipy.stats import spearmanr
import matplotlib.pyplot as plt 
from statsmodels.stats.multitest import multipletests
import random
import pathway_analyzer as pa
import cell_enrichment_utils as ceu


def _load_adata(adata_filename: str) -> AnnData: 
    """load AnnData object from .h5ad file """ 
    if os.path.isfile(adata_filename): 
        if adata_filename[-5:] == '.h5ad': 
            logging.info('loading the data...')
            adata = sc.read_h5ad(adata_filename)
            return adata
        else: 
            logging.error('input file must be in ".h5ad" format')
    else: 
        logging.error('filename for input .h5ad file not provided')


def _load_cell_gmt(): 
    """load pre-existing cell gmt from file"""
    pass


class DataObject(): 
    """docstring"""
    def __init__(
        self, 
        adata: AnnData,
        out_dir: str,
        sample_id: str = 'patient'
    ): 
        self.adata = adata
        self.out_dir = out_dir
        self.sample_id = sample_id

    def create_cell_gmt(
        self, 
        outfile: str, 
        write_dict: bool = True, 
        return_dict: bool = False, 
        permute: bool = False
    ): 
        if outfile[-4:]=='.gmt': 
            if os.path.exists(outfile):
                os.remove(outfile)

            sample_dict = {}
            cell_id = self.adata.obs.index
            sample_id = self.adata.obs[self.sample_id]
            if permute: 
                np.random.shuffle(sample_id)                
            for i in range(len(cell_id)): 
                if sample_id[i] in sample_dict.keys(): 
                    sample_dict[sample_id[i]].append(cell_id[i])
                else: 
                    sample_dict[sample_id[i]] = [cell_id[i]]

            if write_dict: 
                with open(outfile, 'a') as write_out: 
                    for i in sample_dict.keys(): 
                        add_line = i + '\t' + '\t'.join(j for j in sample_dict[i]) + '\n'
                        write_out.write(add_line)

            if return_dict: 
                return sample_dict

        else: 
            logging.error('Output file must use ".gmt" suffix')

    def create_cell_adata(
        self, 
        mode: str, 
        pt_label: str = 'monocle3_pseudotime', 
        gene_set_label: str = 'X_regulonsAUC', 
        celltype_name: str = 'celltype'
    ): 

        if mode == 'pseudotime': 
            pt = self.adata.obs[pt_label]
            self.adata_mod = sc.AnnData(pd.DataFrame(pt)).T ## address warning
            #self.adata_mod = sc.AnnData(X=pd.DataFrame(pt), dtype=pd.DataFrame(pt).dtypes).T ## address warning ## THIS WILL CAUSE ERROR, FIX DIFFERENTLY
        
        elif mode == 'gene_sets': 
            if self.adata.shape[0] == 0: ## can remove this check -- redundant -- addressed further upstream #########
                logging.warning('no matching cells for gene set')
            else: 
                df = self.adata.obsm[gene_set_label]
                self.adata_mod = sc.AnnData(df.T)
                #self.adata_mod = sc.AnnData(X=df.T, dtype=df.T.dtypes) ## address warning ## THIS WILL CAUSE ERROR, FIX DIFFERENTLY
                if celltype_name in self.adata.obs.keys(): 
                    self.adata_mod.var[celltype_name] = self.adata.obs[celltype_name]
                else: 
                    logging.warning('cell type field in input AnnData object is missing or labeled incorrectly')

        else: 
            print('Not an accepted mode. Set "mode" to "pseudotime" or "gene_sets"')


class PseudotimeData(DataObject): 
    """docstring"""
    def run(self, permute: bool = False, pt_cell_gmt: Optional[str] = None): 
        logging.info('running in pseudotime mode')
        if pt_cell_gmt == None: 
            pt_cell_gmt = self.out_dir + 'pt_cell_gmt.gmt'
        self.create_cell_gmt(outfile=pt_cell_gmt, permute=permute)
        self.create_cell_adata(mode='pseudotime') ## this step is being re-run -- save this and speed up permutations ## 
        enrichment_wrapper = ceu.CSEAPyWrapper()
        enrichment_wrapper.csea(self.adata_mod, reactome_gmt=pt_cell_gmt, 
                                data_key=pa.SSGSEA_KEY, weigted_score_type=0, 
                                n_processes=1) ## typo in 'weighted'
        #self.adata_mod.write_h5ad(self.out_dir + 'cell_enrichment_pt.h5ad')
        logging.info('pseudotime section complete')


class ActivityData(DataObject): 
    """docstring"""
    def __init__(
        self, 
        adata: AnnData, 
        out_dir: str, 
        celltype: str, 
        sample_id: str = 'patient'
    ):
        super().__init__(adata, out_dir, sample_id)
        self.celltype = celltype


class CellTypesActivity(): 
    """class to hold 'ActivityData' objects for each cell type"""
    def __init__(
        self, 
        adata: AnnData, 
        celltype_dict: dict, 
        out_dir: str,
        sample_id: str = 'patient'
    ):
        self.adata = adata
        self.celltype_dict = celltype_dict
        self.out_dir = out_dir
        self.sample_id = sample_id

    def add_celltypes(self): 
        """given adata with leiden clusters but no celltype labels, and 
        dict mapping celltypes to leiden cluster IDs
        
        add celltypes to adata in adata.obs['celltype']

        celltype_dict ## dict with celltype string as key and list of leiden as value 
        """

        ## build lookup dictionary from input dictionary 
        values_list_structured = list(self.celltype_dict.values())
        values_list_flat = [str(i) for i in list(map(int, chain.from_iterable(values_list_structured)))]

        if len(set(values_list_flat)) == len(values_list_flat): 
            lookup_dict = {}
            for cluster_id in values_list_flat: 
                for celltype in self.celltype_dict.keys(): 
                    if cluster_id in self.celltype_dict[celltype]: 
                        lookup_dict[cluster_id] = celltype
                        #break ####### check this -- think about best way to break out of loop once found -- avoid unnecessary search ############
                if cluster_id not in lookup_dict.keys(): 
                    logging.warning("WARNING: cluster ID was not found in cell type dictionary")

            self.adata.obs['celltype'] = [lookup_dict[i] for i in list(self.adata.obs['leiden'])]

        else: 
            logging.warning('Ensure no cluster ID is assigned to more than one cell type')

    def run(self, permute: bool = False, act_cell_gmt_prefix: Optional[str] = None): 
        logging.info('running in gene set mode')
        self.add_celltypes()
        self.celltype_objects = []
        for celltype in self.celltype_dict.keys(): 
            adata_cut = self.adata[self.adata.obs['celltype'] == celltype]
            if adata_cut.shape[0] == 0: 
                logging.warning('no cells of type: ' + str(celltype) + ' present')
            else: 
                activity_object = ActivityData(adata=adata_cut, 
                                                celltype=celltype,
                                                out_dir=self.out_dir, 
                                                sample_id=self.sample_id)
                if act_cell_gmt_prefix == None: 
                    act_cell_gmt = self.out_dir + 'act_' + celltype + '_cell_gmt.gmt'
                else: 
                    act_cell_gmt = self.out_dir + act_cell_gmt_prefix + celltype + '_cell_gmt.gmt'
                activity_object.create_cell_gmt(outfile=act_cell_gmt, permute=permute)
                activity_object.create_cell_adata(mode='gene_sets')
                enrichment_wrapper = ceu.CSEAPyWrapper()
                enrichment_wrapper.csea(activity_object.adata_mod, 
                                        reactome_gmt=act_cell_gmt, 
                                        data_key=pa.SSGSEA_KEY, 
                                        weigted_score_type=0, 
                                        n_processes=1) ## typo in 'weighted'
                self.celltype_objects.append(activity_object)
        logging.info('gene set section complete')


"""delete modified anndatas and text file afterward -- keep data on permutations though (log?) 

add option to not add cell types -- if this column already exists in act anndata

should write to files and read from them to save ram, then delete all afterward? 
or just keep in ram? design with options to run either way? not too hard to add

also set this up so that CSEA can be run once and stored, and then correlation() can load those 
results from file. or can just run csea without correlation. etc. keep it modular 

"""


class RunCSEA(): 
    """docstring"""
    def __init__(
        self, 
        adata_pseudotime: str, 
        adata_activity: str, 
        celltypes_dict: dict,
        output_dir: str,
        sample_id: str = 'patient', 
        fdr: float = 0.05, 
        pt_data: PseudotimeData = None, 
        act_data: CellTypesActivity = None, 
        set_num_permutations: int = 1000
    ): 
        self.adata_pt = _load_adata(adata_pseudotime)
        self.adata_act = _load_adata(adata_activity)
        self.celltypes_dict = celltypes_dict
        self.out_dir = output_dir
        self.sample_id = sample_id
        self.fdr = fdr
        self.pt_data = pt_data
        self.act_data = act_data
        self.n_perm = set_num_permutations

    def csea(self,
             just_pt: bool = False
    ):
        logging.info('running CSEA...')
        self.pt_data = PseudotimeData(adata=self.adata_pt, out_dir=self.out_dir, sample_id=self.sample_id)
        self.pt_data.run()
        if not just_pt:
            self.act_data = CellTypesActivity(adata=self.adata_act,
                                        celltype_dict=self.celltypes_dict,
                                        out_dir=self.out_dir,
                                        sample_id=self.sample_id)
            self.act_data.run()
        logging.info('CSEA complete.')

    def correlation(self, 
                    permute: bool = False, 
                    pt_data: Optional[PseudotimeData] = None, 
                    act_data: Optional[CellTypesActivity] = None, 
                    export_enrich_values: bool = True, 
                    n_perm: Optional[int] = None
    ): 
        """
        
        """
        if permute: 
            if pt_data == None or act_data == None: 
                logging.error('"pt_data" and "act_data" must be provided as arguments in permutation mode.')
        if not permute: 
            pt_data = self.pt_data
            act_data = self.act_data
        if self.pt_data == None or self.act_data == None: 
            logging.error('Must run RunCSEA.csea() before running correlation step.')
        else: 
            logging.info('calculating correlations...')
            correlation_results = {}
            pt_df = pt_data.adata_mod.obsm[pa.SSGSEA_KEY]
            if export_enrich_values: 
                if permute: 
                    pt_df.to_csv(self.out_dir + 'pt_perm' + str(n_perm) + '_enrich_scores.csv')
                else: 
                    pt_df.to_csv(self.out_dir + 'pt_main_enrich_scores.csv')
            for object in act_data.celltype_objects: 
                act_df = object.adata_mod.obsm[pa.SSGSEA_KEY]
                if export_enrich_values: 
                    if permute: 
                        act_df.to_csv(self.out_dir + object.celltype + '_perm' + str(n_perm) + '_act_enrich_scores.csv')
                    else: 
                        act_df.to_csv(self.out_dir + object.celltype + '_main_act_enrich_scores.csv')
                if len(pt_df.columns) != len(act_df.index): 
                    logging.warning('sample IDs do not match')
                else: 
                    ## remove any sample IDs not present in both objects
                    samples_remove_act = [i for i in list(act_df.columns) if i not in pt_df.columns] ## can do at higher level -- running step repeatedly #######
                    logging.info("Samples removed from activity score dataframe: " + str(samples_remove_act))
                    for i in samples_remove_act: 
                        act_df = act_df.drop(i, axis=1)
                    samples_remove_pt = [i for i in list(pt_df.columns) if i not in act_df.columns]
                    logging.info("Samples removed from pseudotime dataframe for cell type " + object.celltype + ": " + str(samples_remove_pt))
                    for i in samples_remove_pt: 
                        pt_df = pt_df.drop(i, axis=1)
                correlation_results[object.celltype] = []
                for gene_set in list(act_df.index): ## iterate through gene sets
                    pt_vect = pt_df.T.to_numpy()
                    act_vect = act_df.loc[gene_set, :].to_numpy()
                    if pt_vect.size != act_vect.size:
                        shared_samples = set(pt_df.columns) & set(act_df.columns)                        
                        if len(shared_samples) > 0: ## if any shared samples
                            pt_df_intersect = pt_df.loc[:, list(shared_samples)]
                            act_df_intersect = act_df.loc[:, list(shared_samples)]
                            pt_vect = pt_df_intersect.T.to_numpy() ## repeating lines from above -- rewrite with fewer lines 
                            act_vect = act_df_intersect.loc[gene_set, :].to_numpy() ## 
                            logging.info('used ' + str(len(shared_samples)) + ' for celltype: ' + object.celltype + '; gene set: ' + gene_set)
                        else: 
                            logging.warning('gene set skipped -- no shared samples between pt and act dfs')

                    rho, p = spearmanr(pt_vect, act_vect) ## calculate Spearman Rank correlation and corresponding p-value
                    correlation_results[object.celltype].append((gene_set, rho, p))

        if permute: 
            return correlation_results
        else: 
            self.correlation_results = correlation_results
        logging.info('correlation calculation complete.')

    def permute(self, 
                random_state: int = 1, 
                intermittent_write: bool = True, 
                interval: int = 50, 
                export_enrich_values: bool = True
    ): 
        logging.info('running permutations -- this part may take a while...')

        ## there will be plenty of opportunity for optimization here after the first implementation 
        
        self.permute_raw = {}
        #for celltype in self.celltypes_dict.keys(): 
            #self.permute_raw[celltype] = pd.DataFrame()

        ## for each permutation run, add vector of rho values for all tfs, to pd df for each celltype 
            
        for n in range(self.n_perm): 
            logging.info('running permutation number: ' + str(n))
            print('running permutation number: ' + str(n))
            random.seed(random_state + n) ## check this -- may not be necessary ## 
            pt_data = PseudotimeData(adata=self.adata_pt, 
                                    out_dir=self.out_dir, 
                                    sample_id=self.sample_id)
            act_data = CellTypesActivity(adata=self.adata_act, 
                                        celltype_dict=self.celltypes_dict, 
                                        out_dir=self.out_dir, 
                                        sample_id=self.sample_id)
            pt_gmt = self.out_dir + 'pt_cell_gmt_permute' + str(n) + '.gmt'
            pt_data.run(permute=True, pt_cell_gmt=pt_gmt)
            act_gmt_prefix = 'permute_' + str(n) + '_'
            act_data.run(permute=True, act_cell_gmt_prefix=act_gmt_prefix)
            correlation_results = self.correlation(permute=True, 
                                                   pt_data=pt_data, 
                                                   act_data=act_data, 
                                                   export_enrich_values=export_enrich_values, ## PEP8 fix 
                                                   n_perm=n)

            for celltype in correlation_results: 
                gene_sets = [x[0] for x in correlation_results[celltype]]
                rho_values = [x[1] for x in correlation_results[celltype]]
                df = pd.DataFrame(rho_values).T
                df.columns=gene_sets
                if celltype in self.permute_raw.keys(): 
                    if len(set(self.permute_raw['macrophage'].columns) & set(df.columns)) > 0: ##### 
                        self.permute_raw[celltype] = pd.merge(left=self.permute_raw[celltype], right=df, how='outer') ## CHECK THIS
                    else: 
                        logging.warning('no matching columns to merge on') ###############
                else: 
                    self.permute_raw[celltype] = df
                
            if intermittent_write: 
                if n % interval == 0: 
                    logging.info('CHECKPOINT -- write permutation p values in case run fails later')
                    self.export_permutation_runs(mid_perm=n)

        logging.info('permutations complete.')

    def p_values(self): 
        """docstring"""
        self.permute_results = {}
        for celltype in self.correlation_results.keys():
            results_list = self.correlation_results[celltype]
            gene_sets = [x[0] for x in results_list]
            rho_values = [x[1] for x in results_list]
            results_df = pd.DataFrame(rho_values)
            results_df.columns = ['rho']
            results_df['p_perm'] = np.nan
            results_df['mean_permuted_rho'] = np.nan
            results_df.index = gene_sets
            for gene_set in results_df.index: 
                real_rho = results_df['rho'].loc[gene_set]
                perm_rhos = self.permute_raw[celltype][gene_set]
                n_perm = self.permute_raw[celltype].shape[0] ## make sure this will always work
                mean_rho = np.mean(perm_rhos)
                if real_rho > 0: 
                    p_perm = sum(perm_rhos >= real_rho) / n_perm
                elif real_rho < 0: 
                    p_perm = sum(perm_rhos <= real_rho) / n_perm
                else: 
                    logging.warning('rho is exactly zero') ## check this 
                    p_perm = 0.5 ## check this 
                results_df['p_perm'].loc[gene_set] = p_perm
                results_df['mean_permuted_rho'].loc[gene_set] = mean_rho

            self.permute_results[celltype] = results_df
            for celltype in self.permute_results.keys(): 
                self.permute_results[celltype].to_csv(self.out_dir + celltype + '_results.csv')

        logging.info('calculation of p-values via permutation test complete.')

    def export_permutation_runs(self, mid_perm: Optional[int] = None): 
        """export rho values for all permutation runs"""
        logging.info('writing permutation rho values to file...')
        for celltype in self.permute_raw.keys(): 
            if mid_perm != None: 
                self.permute_raw[celltype].to_csv(self.out_dir + celltype + '_intermittent_' + str(mid_perm) + '_permutation_rhos.csv')
            else: 
                self.permute_raw[celltype].to_csv(self.out_dir + celltype + '_permutation_rhos.csv')

    def export_enrichment_values(self): 
        """export GSEA enrichment score dataframes for pseudotime and TF activity to output .csv files"""
        pass

    def export_sample_pt_scores(self):
        """export sample-level pseudotime scores to file"""
        pt_df = self.pt_data.adata_mod.obsm[pa.SSGSEA_KEY]
        pt_df.to_csv(self.out_dir + 'pt_main_enrich_scores.csv')

