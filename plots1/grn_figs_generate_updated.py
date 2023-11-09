"""run with 'network' conda env """
# import statements
import sys
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt

#sys.path.insert(0, '/Users/klockec/Documents/code/github_vscode/immune-tme-inference/analysis/p3/analysis2/mod_v1')
#from overlap import gmt_to_dict

# workaround for now:
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

# function definitions
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
        outfile: str,
        sig_threshold: float = 0.05
):
    """
    create networkx graph object from graph df and save as graphml file
    """
    celltype_tfs = tf_df[tf_df['celltype'] == celltype]
    directional_tfs = celltype_tfs[celltype_tfs['direction'] == dir]
    tfs_sig = list(directional_tfs[directional_tfs['adjusted_p_value'] < sig_threshold]['tf'])
    links1 = links_cut(tf_list=tfs_sig, links=links)
    G = nx.from_edgelist(links1)
    nx.write_graphml(G=G, path=outfile)

    return links1

def graph_from_final_edgelists(
        list1: list,
        list2: list,
        outfile: str,
        exclude_second: bool = False
):
    """
    abbreviated version of above function -- for intersection across datasets -- generate graphml file

    if exclude_second is true, instead of looking for intersection between lists, will look for items in first list
    and not in second list
    """
    shared = []
    if exclude_second:
        shared = [x for x in list1 if x not in list2]
    else:
        shared = list(set(list1) & set(list2))
    G = nx.from_edgelist(shared)
    nx.write_graphml(G=G, path=outfile)

    return shared

# setup
dataset1_results = "/Users/klockec/OneDrive - Oregon Health & Science University/data/analysis_files/p3/mod1/li_parallel_tfs_run4/combined_bh.csv"
dataset1_links = "/Users/klockec/OneDrive - Oregon Health & Science University/data_backup_rws07890/data_january2023_backup/data/analysis_files/p3/pyscenic_runs/p3/li_run/ctx_output1.gmt"
dataset2_results = "/Users/klockec/OneDrive - Oregon Health & Science University/data/analysis_files/p3/mod1/yost_parallel_tfs_run4/combined_bh.csv"
dataset2_links = "/Users/klockec/OneDrive - Oregon Health & Science University/data_backup_rws07890/data_january2023_backup/data/analysis_files/p3/pyscenic_runs/p3/yost_bcc_run1/ctx_output1.gmt"
dataset3_results = "/Users/klockec/OneDrive - Oregon Health & Science University/data/analysis_files/p3/mod1/wang_parallel_tfs_run4/combined_bh.csv"
dataset3_links = "/Users/klockec/OneDrive - Oregon Health & Science University/data_backup_rws07890/data_january2023_backup/data/analysis_files/p3/pyscenic_runs/p3/wang_run/ctx_output1.gmt"
out_dir1 = "/Users/klockec/OneDrive - Oregon Health & Science University/data/analysis_files/p3/mod1/li_grn_graphs/"
out_dir2 = "/Users/klockec/OneDrive - Oregon Health & Science University/data/analysis_files/p3/mod1/yost_grn_graphs/"
out_dir3 = "/Users/klockec/OneDrive - Oregon Health & Science University/data/analysis_files/p3/mod1/wang_grn_graphs/"
out_dir_shared = "/Users/klockec/OneDrive - Oregon Health & Science University/data/analysis_files/p3/mod1/tumor_tumor_grn_graphs/"
out_dir_tumor_viral = "/Users/klockec/OneDrive - Oregon Health & Science University/data/analysis_files/p3/mod1/tumor_viral_grn_graphs/"
out_dir_tumor_not_viral = "/Users/klockec/OneDrive - Oregon Health & Science University/data/analysis_files/p3/mod1/tumor_not_viral_grn_graphs/"
out_dir_viral_not_tumor = "/Users/klockec/OneDrive - Oregon Health & Science University/data/analysis_files/p3/mod1/viral_not_tumor_grn_graphs/"

# load the data
results1 = pd.read_csv(dataset1_results)
links1 = gmt_to_dict(dataset1_links)
results2 = pd.read_csv(dataset2_results)
links2 = gmt_to_dict(dataset2_links)
results3 = pd.read_csv(dataset3_results)
links3 = gmt_to_dict(dataset3_links)

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

# for subsetting -- add later if needed
"""
# cut down tfs by in- or out-degree 
#https://stackoverflow.com/questions/22391433/count-the-frequency-that-a-value-occurs-in-a-dataframe-column
out_deg = [i[0] for i in links]
in_deg = [i[1] for i in links]
df_degs = pd.DataFrame(out_deg, columns=['out_degree'])
df_degs['in_degree'] = in_deg
df_out_deg = pd.DataFrame(df_degs['out_degree'].value_counts())
#df_in_deg = df_degs['in_degree'].value_counts()
df_out_deg.columns = ['out_degree']
df_out_deg['tf'] = df_out_deg.index
out_deg_min = 2
df_out_deg_cut = df_out_deg[df_out_deg['out_degree'] >= out_deg_min]
#links2 = [i for i in links if (i[0] in list(df_out_deg_cut['tf']) and (i[1] in list(df_out_deg_cut['tf'])))]  # cut again  #######
links2 = [i for i in links if i[0] in list(df_out_deg_cut['tf'])]  # cut again
print(str(len(links2)))
"""

# find overlapping graphs across datasets -- tumor / tumor
shared_cd8_high = graph_from_final_edgelists(list1=li_cd8_high, list2=yost_cd8_high, outfile=(out_dir_shared + 'cd8_high_grn.graphml'))
shared_cd8_low = graph_from_final_edgelists(list1=li_cd8_low, list2=yost_cd8_low, outfile=(out_dir_shared + 'cd8_low_grn.graphml'))
shared_macrophage_high = graph_from_final_edgelists(list1=li_mac_high, list2=yost_mac_high, outfile=(out_dir_shared + 'macrophage_high_grn.graphml'))
shared_macrophage_low = graph_from_final_edgelists(list1=li_mac_low, list2=yost_mac_low, outfile=(out_dir_shared + 'macrophage_low_grn.graphml'))

# shared across all three -- tumor / viral
tumor_viral_cd8_high = graph_from_final_edgelists(list1=shared_cd8_high, list2=wang_cd8_high, outfile=(out_dir_tumor_viral + 'cd8_high_grn.graphml'))
tumor_viral_cd8_low = graph_from_final_edgelists(list1=shared_cd8_low, list2=wang_cd8_low, outfile=(out_dir_tumor_viral + 'cd8_low_grn.graphml'))
tumor_viral_macrophage_high = graph_from_final_edgelists(list1=shared_macrophage_high, list2=wang_mac_high, outfile=(out_dir_tumor_viral + 'macrophage_high_grn.graphml'))
tumor_viral_macrophage_low = graph_from_final_edgelists(list1=shared_macrophage_low, list2=wang_mac_low, outfile=(out_dir_tumor_viral + 'macrophage_low_grn.graphml'))

# tumor but not viral (tumor-specific)
tumor_not_viral_cd8_high = graph_from_final_edgelists(list1=shared_cd8_high, list2=wang_cd8_high, outfile=(out_dir_tumor_not_viral + 'cd8_high_grn.graphml'), exclude_second=True)
tumor_not_viral_cd8_low = graph_from_final_edgelists(list1=shared_cd8_low, list2=wang_cd8_low, outfile=(out_dir_tumor_not_viral + 'cd8_low_grn.graphml'), exclude_second=True)
tumor_not_viral_macrophage_high = graph_from_final_edgelists(list1=shared_macrophage_high, list2=wang_mac_high, outfile=(out_dir_tumor_not_viral + 'macrophage_high_grn.graphml'), exclude_second=True)
tumor_not_viral_macrophage_low = graph_from_final_edgelists(list1=shared_macrophage_low, list2=wang_mac_low, outfile=(out_dir_tumor_not_viral + 'macrophage_low_grn.graphml'), exclude_second=True)

# viral but not tumor (viral-specific)
"""
to make this comparison correctly, need to ensure link is found in neither tumor dataset 
"""
#viral_not_tumor_cd8_high = graph_from_final_edgelists(list1=shared_cd8_high, list2=wang_cd8_high, outfile=(out_dir_tumor_not_viral + 'cd8_high_grn.graphml'), exclude_second=True)
#viral_not_tumor_cd8_low = graph_from_final_edgelists(list1=shared_cd8_low, list2=wang_cd8_low, outfile=(out_dir_tumor_not_viral + 'cd8_low_grn.graphml'), exclude_second=True)
#viral_not_tumor_macrophage_high = graph_from_final_edgelists(list1=shared_macrophage_high, list2=wang_mac_high, outfile=(out_dir_tumor_not_viral + 'macrophage_high_grn.graphml'), exclude_second=True)
#viral_not_tumor_macrophage_low = graph_from_final_edgelists(list1=shared_macrophage_low, list2=wang_mac_low, outfile=(out_dir_tumor_not_viral + 'macrophage_low_grn.graphml'), exclude_second=True)

