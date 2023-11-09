"""generate violin plots

save-as of 'violins1.py'
"""

## import statements
import scanpy as sc
import os
import seaborn as sns
import matplotlib.pyplot as plt
from anndata import AnnData
from typing import Optional

## file paths
dir = '/Users/klockec/OneDrive - Oregon Health & Science University/data_backup_rws07890/data_january2023_backup/'  # prefix updated for local
#li_output_filename = 'data/li_data/li_t_monocle1.h5ad'  # which to plot? check
li_output_filename = 'data/li_data/li_t_monocle2.h5ad'  # which to plot? check
#yost_output_filename = 'data/yost_data/yost_t_monocle1.h5ad'
#yost_output_filename = 'data/yost_data/yost_t_monocle1_v2.h5ad'
yost_output_filename = 'data/yost_data/yost_t_monocle1_v2_cut.h5ad'  # only keep cells with pt values
wang_output_filename = 'data/wang_data/wang_t_monocle1.h5ad'
#plots_store = '/Users/klockec/OneDrive - Oregon Health & Science University/data/analysis_files/p3/manuscript_figs/yost_t/'
plots_store = '/Users/klockec/OneDrive - Oregon Health & Science University/data/analysis_files/p3/manuscript_figs/violins/'

## setup
random1 = 1
os.chdir(plots_store)

## load the data
adata_li = sc.read_h5ad(dir + li_output_filename)
adata_yost = sc.read_h5ad(dir + yost_output_filename)
adata_wang = sc.read_h5ad(dir + wang_output_filename)

## subset to just bcc
adata_yost_bcc = adata_yost[adata_yost.obs['tumor_type'] == 'bcc']

## function to generate violin plots
def violin_plot(
        adata: AnnData,
        sample_id: str = 'patient',
        pt_name: str = 'monocle3_pseudotime',
        show: bool = True,
        save: Optional[str] = None,
        order: Optional[list] = None,
        width: Optional[float] = None
):
    sns.set_theme(style="whitegrid")
    ax = ''  # scope correctly
    if order is not None and width is None:
        ax = sns.violinplot(x=sample_id, y=pt_name, data=adata.obs, order=order)
    elif order is not None and width is not None:
        ax = sns.violinplot(x=sample_id, y=pt_name, data=adata.obs, order=order, width=width)
    else:
        ax = sns.violinplot(x=sample_id, y=pt_name, data=adata.obs)
    ax.set(xlabel='patient ID', ylabel='CD8 T cell pseudotime')
    plt.title('CD8 T cell Exhaustion by Sample')

    if show:
        plt.show(block=True)
    if save is not None:
        plt.savefig(save)

## generate violin plots
#violin_plot(adata=adata_li, order=['p27_PBMC', 'p17_PBMC', 'p13_PBMC', 'p17', 'p1', 'p11', 'p24', 'p16', 'p13', 'p15', 'p21', 'p25', 'p3', 'p26', 'p18', 'p12', 'p19', 'p23', 'p27'])
#violin_plot(adata=adata_li, order=['p27_PBMC', 'p17_PBMC', 'p13_PBMC', 'p17', 'p1', 'p11', 'p24', 'p16', 'p13', 'p15', 'p21', 'p25', 'p3', 'p26', 'p18', 'p12', 'p19', 'p23', 'p27'], width=1.2)
violin_plot(adata=adata_li, order=['p17', 'p1', 'p11', 'p24', 'p16', 'p13', 'p15', 'p21', 'p25', 'p3', 'p26', 'p18', 'p12', 'p19', 'p23', 'p27'], width=1.2)
#violin_plot(adata=adata_li, save='_li_pt_violin.pdf')
violin_plot(adata=adata_wang, sample_id='batch', order=['2', '3', '5', '1', '0', '4'])
#violin_plot(adata=adata_wang, sample_id='batch', order=['2', '3', '5', '1', '0', '4'], save='_wang_pt_violin.pdf')
#violin_plot(adata=adata_yost, sample_id='tumor_type')
violin_plot(adata=adata_yost, sample_id='sample', order=['su003', 'su007', 'su004', 'su002', 'su012', 'su001', 'su010', 'su005', 'su009', 'su008', 'su006', 'su013', 'su011', 'su014', ])
#violin_plot(adata=adata_yost_bcc, sample_id='sample')

"""
find the sample score gmt files 

plot prelim violins, then iterate: 

    order by sample score 
    
    annotate with sample score (manually for now if needed) 
    
    
'/Users/klockec/OneDrive - Oregon Health & Science University/data/analysis_files/p3/csea_results/li_run1/pt_enrich_main.csv' 

'/Users/klockec/OneDrive - Oregon Health & Science University/data/analysis_files/p3/mod1/wang_run8_csea_oop_pt_main_enrich_scores.csv' 

"""
