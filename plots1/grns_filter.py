"""run with 'sc2' conda env """

# import statements
import scanpy as sc
import pandas as pd
import networkx as nx
from typing import Optional
from anndata import AnnData
from scipy.stats import spearmanr

# setup
li_adata_filename = "/Users/klockec/OneDrive - Oregon Health & Science University/data_backup_rws07890/data_january2023_backup/data/li_data/li_dd_viz_ready_regs_updated.h5ad"  # from 'config_li_tf_v4.json'
yost_adata_filename = "/Users/klockec/OneDrive - Oregon Health & Science University/data_backup_rws07890/data_january2023_backup/data/yost_data/yost_bcc_full_processed_dd_regs_mod2.h5ad"  # from 'config_yost_tf_v4.json'
wang_adata_filename = "/Users/klockec/OneDrive - Oregon Health & Science University/data_backup_rws07890/data_january2023_backup/data/wang_data/wang_dd_viz_ready_regs_cut_updated.h5ad"  # from 'config_wang_tf_v4.json'
dataset1_results = "/Users/klockec/OneDrive - Oregon Health & Science University/data/analysis_files/p3/mod1/li_parallel_tfs_run4/combined_bh.csv"
dataset1_links = "/Users/klockec/OneDrive - Oregon Health & Science University/data_backup_rws07890/data_january2023_backup/data/analysis_files/p3/pyscenic_runs/p3/li_run/ctx_output1.gmt"
dataset2_results = "/Users/klockec/OneDrive - Oregon Health & Science University/data/analysis_files/p3/mod1/yost_parallel_tfs_run4/combined_bh.csv"
dataset2_links = "/Users/klockec/OneDrive - Oregon Health & Science University/data_backup_rws07890/data_january2023_backup/data/analysis_files/p3/pyscenic_runs/p3/yost_bcc_run1/ctx_output1.gmt"
dataset3_results = "/Users/klockec/OneDrive - Oregon Health & Science University/data/analysis_files/p3/mod1/wang_parallel_tfs_run4/combined_bh.csv"
dataset3_links = "/Users/klockec/OneDrive - Oregon Health & Science University/data_backup_rws07890/data_january2023_backup/data/analysis_files/p3/pyscenic_runs/p3/wang_run/ctx_output1.gmt"
dataset1_link_scores = "/Users/klockec/OneDrive - Oregon Health & Science University/data_backup_rws07890/data_january2023_backup/data/analysis_files/p3/pyscenic_runs/p3/li_run/grnboost_out.tsv"
dataset2_link_scores = "/Users/klockec/OneDrive - Oregon Health & Science University/data_backup_rws07890/data_january2023_backup/data/analysis_files/p3/pyscenic_runs/p3/yost_bcc_run1/grnboost_out.tsv"
dataset3_link_scores = "/Users/klockec/OneDrive - Oregon Health & Science University/data_backup_rws07890/data_january2023_backup/data/analysis_files/p3/pyscenic_runs/p3/wang_run/grnboost_out.tsv"
out_dir1 = "/Users/klockec/OneDrive - Oregon Health & Science University/data/analysis_files/p3/mod1/li_grn_graphs/"
out_dir2 = "/Users/klockec/OneDrive - Oregon Health & Science University/data/analysis_files/p3/mod1/yost_grn_graphs/"
out_dir3 = "/Users/klockec/OneDrive - Oregon Health & Science University/data/analysis_files/p3/mod1/wang_grn_graphs/"

# function definitions
def gmt_to_dict(
        gmt_filename: str
):
    """

    """
    gmt_dict = {}
    with open(gmt_filename, 'r') as read_file:
        for line in read_file:
            line = line.rstrip()
            line = line.split('\t')
            tf = line[0][:-3]
            targets = line[2:]
            gmt_dict[tf] = targets

    return gmt_dict

def links_cut(
        tf_list: list,
        links: dict,
        self_links: bool = False
):
    """
    built from 'find_grn_links()' in 'overlap.py' -- look there to add in consensus GRN functionality
    """
    links_keep = []
    for tf in tf_list:
        if tf in links.keys():
            targets = links[tf]
            for target in targets:
                if target in tf_list:
                    if self_links:
                        links_keep.append((tf, target))
                    else:
                        if tf != target:
                            links_keep.append((tf, target))

    return links_keep

def save_graph(
        tf_df: list,
        links: dict,
        celltype: str,
        dir: str,
        sig_threshold: float = 0.05,
        write_graph: bool = True,
        print_n_tfs: bool = False,
        outfile: Optional[str] = None
):
    """
    create networkx graph object from graph df and save as graphml file
    """
    celltype_tfs = tf_df[tf_df['celltype'] == celltype]
    directional_tfs = celltype_tfs[celltype_tfs['direction'] == dir]
    tfs_sig = list(directional_tfs[directional_tfs['adjusted_p_value'] < sig_threshold]['tf'])
    if print_n_tfs:
        print('tfs: ' + str(len(tfs_sig)))
    links1 = links_cut(tf_list=tfs_sig, links=links)
    if write_graph:
        G = nx.from_edgelist(links1)
        nx.write_graphml(G=G, path=outfile)

    return links1

def save_graph_simple(
        edges: list,
        outfile: str
):
    """
    function to save graph from edge list
    """
    G = nx.from_edgelist(edges)
    nx.write_graphml(G=G, path=outfile)

def tf_correlation(
        tf1: str,
        tf2: str,
        df: pd.DataFrame
):
    """
    calculate the correlation between AUCell activity scores for two transcription factors in a chosen subset of cells
    """
    vec1 = df[[(tf1 + '(+)')]]
    vec2 = df[[(tf2 + '(+)')]]
    corr = spearmanr(vec1, vec2)
    p_val = corr[1]

    return p_val

def prune_edges(
        links: list,
        adata: AnnData,
        p_threshold: float = 0.5
):
    """
    choose subset of GRN edges -- use correlation b/t TF AUCell activities
    """
    links_keep = []
    df = adata.obsm['X_regulonsAUC']
    for pair in links:
        p_val = tf_correlation(tf1=pair[0], tf2=pair[1], df=df)
        if p_val > p_threshold:
            links_keep.append(pair)

    return links_keep

def scores_to_dict(
        scores_df: pd.DataFrame
):
    """scores from df to dict for easy lookup """
    scores_dict = {}
    for nrow in range(scores_df.shape[0]):
        row = scores_df.iloc[nrow]
        if row[0] not in scores_dict.keys():
            scores_dict[row[0]] = {}
        scores_dict[row[0]][row[1]] = row[2]

    return scores_dict

def prune_edges_alt(
        links: list,
        scores: pd.DataFrame,
        score_threshold: float = 1.0
):
    """
    choose subset of GRN edges -- use scores from GRNBoost
    """
    links_keep = []

    scores_cut = scores[scores['importance'] > score_threshold]
    tfs = set(scores_cut['TF'])
    targets = set(scores_cut['target'])
    for pair in links:
        if pair[0] in tfs:
            if pair[1] in targets:
                links_keep.append(pair)

    return links_keep

# load the data
adata_li = sc.read_h5ad(li_adata_filename)
adata_yost = sc.read_h5ad(yost_adata_filename)
adata_wang = sc.read_h5ad(wang_adata_filename)
results1 = pd.read_csv(dataset1_results)
links1 = gmt_to_dict(dataset1_links)
link_scores1 = pd.read_csv(dataset1_link_scores, sep='\t')
results2 = pd.read_csv(dataset2_results)
links2 = gmt_to_dict(dataset2_links)
link_scores2 = pd.read_csv(dataset2_link_scores, sep='\t')
results3 = pd.read_csv(dataset3_results)
links3 = gmt_to_dict(dataset3_links)
link_scores3 = pd.read_csv(dataset3_link_scores, sep='\t')

# save graphs
# ready -- li mac high, li cd8 low, yost cd8 low

#li_cd8_high = save_graph(tf_df=results1, links=links1, celltype='CD8_T_cells', dir='high', write_graph=False, sig_threshold=0.0000000000000000000000000001, print_n_tfs=True)
#li_cd8_high_pruned = prune_edges_alt(links=li_cd8_high, scores=link_scores1, score_threshold=60)
#print('edges: ' + str(len(li_cd8_high_pruned)))
#save_graph_simple(edges = li_cd8_high_pruned, outfile=(out_dir1 + 'cd8_high_grn.graphml'))

#li_cd8_low = save_graph(tf_df=results1, links=links1, celltype='CD8_T_cells', dir='low', write_graph=False, sig_threshold=0.00000000000000000000000000000001, print_n_tfs=True)
#li_cd8_low_pruned = prune_edges_alt(links=li_cd8_low, scores=link_scores1, score_threshold=10)
#print('edges: ' + str(len(li_cd8_low_pruned)))
#save_graph_simple(edges = li_cd8_low_pruned, outfile=(out_dir1 + 'cd8_low_grn2.graphml'))

#li_mac_high = save_graph(tf_df=results1, links=links1, celltype='macrophage', dir='high', write_graph=False, sig_threshold=0.0000001, print_n_tfs=True)
#li_mac_high_pruned = prune_edges_alt(links=li_mac_high, scores=link_scores1, score_threshold=50)
#print('edges: ' + str(len(li_mac_high_pruned)))
#save_graph_simple(edges=li_mac_high_pruned, outfile=(out_dir1 + 'macrophage_high_grn2.graphml'))

#li_mac_low = save_graph(tf_df=results1, links=links1, celltype='macrophage', dir='low', write_graph=False, sig_threshold=0.05, print_n_tfs=True)
#li_mac_low_pruned = prune_edges_alt(links=li_mac_low, scores=link_scores1, score_threshold=50)
#print('edges: ' + str(len(li_mac_low_pruned)))
#save_graph_simple(edges=li_mac_low_pruned, outfile=(out_dir1 + 'macrophage_low_grn.graphml'))

#yost_cd8_low = save_graph(tf_df=results2, links=links2, celltype='CD8_T_cells', dir='low', write_graph=False, sig_threshold=0.0000000000000000000000000000000001, print_n_tfs=True)
#yost_cd8_low_pruned = prune_edges_alt(links=yost_cd8_low, scores=link_scores2, score_threshold=100)
#print('edges: ' + str(len(yost_cd8_low_pruned)))
#save_graph_simple(edges=yost_cd8_low_pruned, outfile=(out_dir2 + 'cd8_low_grn2.graphml'))

yost_mac_high = save_graph(tf_df=results2, links=links2, celltype='macrophage', dir='high', write_graph=False, sig_threshold=0.0001, print_n_tfs=True)
yost_mac_high_pruned = prune_edges_alt(links=yost_mac_high, scores=link_scores2, score_threshold=110)
print('edges: ' + str(len(yost_mac_high_pruned)))
save_graph_simple(edges=yost_mac_high_pruned, outfile=(out_dir2 + 'macrophage_high_grn2.graphml'))

#yost_mac_high = save_graph(tf_df=results2, links=links2, celltype='macrophage', dir='high', outfile=(out_dir2 + 'macrophage_high_grn.graphml'))


exit()

# generate graphml files
li_cd8_high = save_graph(tf_df=results1, links=links1, celltype='CD8_T_cells', dir='high', outfile=(out_dir1 + 'cd8_high_grn.graphml'))
li_cd8_low = save_graph(tf_df=results1, links=links1, celltype='CD8_T_cells', dir='low', outfile=(out_dir1 + 'cd8_low_grn.graphml'))
li_mac_high = save_graph(tf_df=results1, links=links1, celltype='macrophage', dir='high', outfile=(out_dir1 + 'macrophage_high_grn.graphml'))
li_mac_low = save_graph(tf_df=results1, links=links1, celltype='macrophage', dir='low', outfile=(out_dir1 + 'macrophage_low_grn.graphml'))
yost_cd8_high = save_graph(tf_df=results2, links=links2, celltype='CD8_T_cells', dir='high', outfile=(out_dir2 + 'cd8_high_grn.graphml'))
yost_cd8_low = save_graph(tf_df=results2, links=links2, celltype='CD8_T_cells', dir='low', outfile=(out_dir2 + 'cd8_low_grn.graphml'))
yost_mac_high = save_graph(tf_df=results2, links=links2, celltype='macrophage', dir='high', outfile=(out_dir2 + 'macrophage_high_grn.graphml'))
yost_mac_low = save_graph(tf_df=results2, links=links2, celltype='macrophage', dir='low', outfile=(out_dir2 + 'macrophage_low_grn.graphml'))
wang_cd8_high = save_graph(tf_df=results3, links=links3, celltype='CD8_T_cells', dir='high', outfile=(out_dir3 + 'cd8_high_grn.graphml'))
wang_cd8_low = save_graph(tf_df=results3, links=links3, celltype='CD8_T_cells', dir='low', outfile=(out_dir3 + 'cd8_low_grn.graphml'))
wang_mac_high = save_graph(tf_df=results3, links=links3, celltype='macrophage', dir='high', outfile=(out_dir3 + 'macrophage_high_grn.graphml'))
wang_mac_low = save_graph(tf_df=results3, links=links3, celltype='macrophage', dir='low', outfile=(out_dir3 + 'macrophage_low_grn.graphml'))