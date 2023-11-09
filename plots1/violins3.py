"""generate violin plots

save-as of 'violins2.py'

scores (for order) from:

in:

'/Users/klockec/OneDrive - Oregon Health & Science University/data/analysis_files/p3/mod1/pt_optimize/'

li_run4_csea_oop/
yost_run4_csea_oop/
wang_run4_csea_oop/
"""

# import statements
import scanpy as sc
import os
import seaborn as sns
import matplotlib.pyplot as plt
from anndata import AnnData
from typing import Optional

## file paths
li_pt_filename = '/Users/klockec/OneDrive - Oregon Health & Science University/data_backup_rws07890/data_january2023_backup/data/li_data/li_t_monocle_run4.h5ad'  # from 'p3/analysis2/monocle3_r/run4_collect_li.py'
yost_pt_filename = '/Users/klockec/OneDrive - Oregon Health & Science University/data_backup_rws07890/data_january2023_backup/data/yost_data/yost_t_monocle_run4.h5ad'  # from 'p3/analysis2/monocle3_r/run4_collect_yost.py'
wang_pt_filename = '/Users/klockec/OneDrive - Oregon Health & Science University/data_backup_rws07890/data_january2023_backup/data/wang_data/wang_t_monocle_run4.h5ad'  # from 'p3/analysis2/monocle3_r/run4_collect_wang.py'
plots_store = '/Users/klockec/OneDrive - Oregon Health & Science University/data/analysis_files/p3/manuscript_figs/violins2/'

# setup
random1 = 1
os.chdir(plots_store)

# load the data
adata_li = sc.read_h5ad(li_pt_filename)
adata_yost = sc.read_h5ad(yost_pt_filename)
adata_wang = sc.read_h5ad(wang_pt_filename)

# remove 'inf' values (from other partitions -- monocle3 clustering)
adata_li = adata_li[adata_li.obs['monocle3_pseudotime'] != float('inf')]
adata_yost = adata_yost[adata_yost.obs['monocle3_pseudotime'] != float('inf')]
adata_wang = adata_wang[adata_wang.obs['monocle3_pseudotime'] != float('inf')]

# subset to just bcc
adata_yost_bcc = adata_yost[adata_yost.obs['tumor_type'] == 'bcc']

# function to generate violin plots
def violin_plot(
        adata: AnnData,
        sample_id: str = 'patient',
        pt_name: str = 'monocle3_pseudotime',
        fig_size: tuple = (10, 6),
        font_scale: float = 1,
        show: bool = True,
        save: Optional[str] = None,
        order: Optional[list] = None,
        width: Optional[float] = None,
        scale: Optional[str] = None
):
    plt.figure(figsize=fig_size)
    sns.set_theme(style="whitegrid", font_scale=font_scale)
    #sns.set_style("whitegrid")
    #sns.set(font_scale=font_scale)
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

# generate violin plots
font_scale = 1.6
violin_plot(adata=adata_li, order=['p1', 'p24', 'p13', 'p16', 'p17', 'p15', 'p26', 'p11', 'p21', 'p18', 'p12', 'p19', 'p3', 'p23', 'p25', 'p27'], save='li_pt_violin.pdf', show=False, fig_size=(16,6), font_scale=font_scale)
#violin_plot(adata=adata_li, order=['p1', 'p24', 'p13', 'p16', 'p17', 'p15', 'p26', 'p11', 'p21', 'p18', 'p12', 'p19', 'p3', 'p23', 'p25', 'p27'])
#plt.close()
violin_plot(adata=adata_yost, sample_id='sample', order=['su002', 'su007', 'su009', 'su012', 'su004', 'su001', 'su003', 'su010', 'su005', 'su006', 'su011', 'su008', 'su014', 'su013'], save='yost_pt_violin.pdf', show=False, fig_size=(14,6), font_scale=font_scale)
#violin_plot(adata=adata_yost, sample_id='sample', order=['su002', 'su007', 'su009', 'su012', 'su004', 'su001', 'su003', 'su010', 'su005', 'su006', 'su011', 'su008', 'su014', 'su013'])
#plt.close()
violin_plot(adata=adata_wang, sample_id='batch', order=['2', '3', '0', '1', '4', '5'], save='wang_pt_violin.pdf', show=False, fig_size=(10,6), font_scale=font_scale)
#violin_plot(adata=adata_wang, sample_id='batch', order=['2', '3', '0', '1', '4', '5'])
#plt.close()
