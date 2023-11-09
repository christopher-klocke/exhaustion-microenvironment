""" run with 'sc2' conda environment"""
#  import statements
import logging
import json
import pandas as pd
import scipy.stats as stats
from statsmodels.stats.multitest import multipletests
import time
import networkx as nx
import matplotlib.pyplot as plt


#  function definitions
def gmt_to_dict(
        gmt_filename: str
):
    gmt_dict = {}
    with open(gmt_filename, 'r') as read_file:
        for line in read_file:
            line = line.rstrip()
            line = line.split('\t')
            tf = line[0][:-3]
            targets = line[2:]
            gmt_dict[tf] = targets

    return gmt_dict

def setup(config_file: str):
    config = json.load(config_file)
    results1 = pd.read_csv(config['files']['dataset1_results'])
    results2 = pd.read_csv(config['files']['dataset2_results'])
    links1 = gmt_to_dict(config['files']['dataset1_links'])
    links2 = gmt_to_dict(config['files']['dataset2_links'])
    out_dir = config['files']['out_dir']
    celltypes = config['parameters']['celltypes_compare']
    tf_num = int(config['parameters']['tf_num'])

    return results1, results2, links1, links2, out_dir, celltypes, tf_num

def overlap_stage1(
        res1: pd.DataFrame,
        res2: pd.DataFrame,
        celltypes: list,
        tf_background: int,
        direction: str = 'high',
        alpha: float = 0.05
):
    res1 = res1[['celltype', 'tf', 'direction', 'adjusted_p_value']]  # select necessary columns
    res2 = res2[['celltype', 'tf', 'direction', 'adjusted_p_value']]
    res1 = res1[res1['adjusted_p_value'] < alpha]  # only significant results
    res2 = res2[res2['adjusted_p_value'] < alpha]
    if direction != 'high' and direction != 'low':
        logging.warning("direction must be set to 'high' or 'low'")
    res1 = res1[res1['direction'] == direction]  # select only high or low results
    res2 = res2[res2['direction'] == direction]
    merge_cols = ['celltype', 'tf']
    overlap_tfs = pd.merge(res1, res2, how='inner', left_on=merge_cols,
                           right_on=merge_cols)  # https://stackoverflow.com/questions/41815079/pandas-merge-join-two-data-frames-on-multiple-columns

    celltype_dfs = []
    for celltype in celltypes:
        overlap_n = overlap_tfs[overlap_tfs['celltype'] == celltype]
        sig1 = res1[res1['celltype'] == celltype]
        sig2 = res2[res2['celltype'] == celltype]
        a = overlap_n.shape[0]
        b = sig1.shape[0] - a
        c = sig2.shape[0] - a
        d = tf_background - (a + b + c)
        compare = [[a, b], [c, d]]
        results = stats.fisher_exact(compare, alternative='greater')
        celltype_results = [direction, results[1], results[0], a, b, c, d]
        df = pd.DataFrame(celltype_results).T
        df.columns = ['direction', 'p_value', 'statistic', 'shared_tfs', 'just_A_tfs', 'just_B_tfs', 'tfs_neither']
        new_col = [celltype]
        df['celltype'] = new_col
        df = df[['celltype', 'direction', 'p_value', 'statistic', 'shared_tfs', 'just_A_tfs', 'just_B_tfs', 'tfs_neither']]  # fix column order
        celltype_dfs.append(df)

    results_df = pd.concat(celltype_dfs)

    return results_df, overlap_tfs

def overlap_stage2(
        links1: dict,
        links2: dict,
        celltypes: list,
        s1_res: pd.DataFrame,
        tf_target_background: int = 410583  # from dorothea, per R dorothea package; from: 'dorothea_hs %>% filter(mor == 1)'
):
    """
    calculate the number of shared and non-shared TF-target relationships between the two datasets

    perform overlap analysis calculation with Fisher exact test
    """
    shared_links = 0
    links_A = 0
    links_B = 0
    tfs_A = [i for i in links1.keys() if i not in links2.keys()]
    tfs_B = [i for i in links2.keys() if i not in links1.keys()]
    tfs_shared = [i for i in links1.keys() if i in links2.keys()]
    for tf in tfs_A:
        links_A += len(links1[tf])
    for tf in tfs_B:
        links_B += len(links2[tf])
    for tf in tfs_shared:
        targets1 = links1[tf]
        targets2 = links2[tf]
        links_A += len([i for i in targets1 if i not in targets2])
        links_B += len([i for i in targets2 if i not in targets1])
        shared_links += len([i for i in targets1 if i in targets2])

    tf_target_background = tf_target_background - (shared_links + links_A + links_B)
    compare = [[shared_links, links_A], [links_B, tf_target_background]]
    results = stats.fisher_exact(compare, alternative='greater')
    links_results = [results[1], results[0], shared_links, links_A, links_B, tf_target_background]
    df = pd.DataFrame(links_results).T
    df.columns = ['p_value', 'statistic', 'shared_links', 'just_A_links', 'just_B_links', 'links_neither']

    return df

def pval_correct(
        results: pd.DataFrame,
        method: str = 'fdr_bh',
        alpha: float = 0.05
):
    correction_input = results['p_value'].to_numpy()
    reject, pvals_corrected, alphacSidak, alphacBonf = multipletests(correction_input, alpha=alpha, method=method)
    results['adjusted_p_value'] = pvals_corrected
    results['reject'] = reject

    return results

def find_grn_links(
        tf_list: list,
        links1: dict,
        links2: dict,
        self_links: bool = False
):
    """
    -takes a list of significant shared tfs as input
    -also takes 2 tf-target link dictionaries as input
    -returns any tf-target links within list of transcription factors as list of 2-tuples
    """
    shared_links = []
    for tf in tf_list:
        if tf in links1.keys() and tf in links2.keys():
            targets1 = links1[tf]
            targets2 = links2[tf]
            shared = [x for x in targets1 if x in targets2]
            shared_sig = [x for x in shared if x in tf_list]
            if len(shared_sig) > 0:
                for target in shared_sig:
                    if self_links:
                        shared_links.append((tf, target))
                    else:
                        if tf != target:
                            shared_links.append((tf, target))

    return shared_links

def grn_viz(
        nodes: list,
        edges: list,
        show: bool = True
):
    #G = nx.from_edgelist(edges)
    #G = nx.DiGraph(G)  ##
    G = nx.DiGraph()
    G.add_edges_from(edges)
    seed = 13648  # Seed random number generators for reproducibility
    font_size = 4.5
    node_color = "blue"
    node_size = 350
    spacing = 0.7
    pos = nx.spring_layout(G, seed=seed, k=spacing)
    #pos = nx.shell_layout(G)
    #pos = nx.spectral_layout(G)
    #cmap = plt.cm.plasma
    cmap = plt.cm.Blues
    node_sizes = [node_size] * len(G)
    #nodes = nx.draw_networkx_nodes(G, pos, node_size=node_sizes, node_color="indigo")
    nodes = nx.draw_networkx_nodes(G, pos, node_color=node_color, node_size=node_sizes, alpha=0.6, label=True)
    edges = nx.draw_networkx_edges(
        G,
        pos,
        #node_size=node_sizes,
        arrowstyle="->",
        arrowsize=10,
        #edge_color=edge_colors,
        edge_cmap=cmap,
        width=2,
    )
    nx.draw_networkx_labels(G, pos, font_size=font_size)
    #nx.draw_networkx(G, width=weights*1)
    #nx.draw_networkx(G)
    if show:
        plt.show()

def overlap_analysis(config_file: str):
    results1, results2, links1, links2, out_dir, celltypes, tf_num = setup(config_file)

    #  stage 1
    s1_high, tfs_high = overlap_stage1(res1=results1, res2=results2, celltypes=celltypes, tf_background=tf_num, direction='high')
    s1_low, tfs_low = overlap_stage1(res1=results1, res2=results2, celltypes=celltypes, tf_background=tf_num, direction='low')
    s1 = pd.concat([s1_high, s1_low])
    s1 = pval_correct(s1, method='bonferroni')
    s1.sort_values(by='celltype', inplace=True)
    s1.sort_values(by='adjusted_p_value', inplace=True)
    tfs_high.to_csv(out_dir + 'stage1_tfs_high.csv')
    tfs_low.to_csv(out_dir + 'stage1_tfs_low.csv')
    print(s1)
    s1.to_csv(out_dir + 'stage1_results.csv')
    #  stage 2
    s2 = overlap_stage2(links1=links1, links2=links2, celltypes=celltypes, s1_res=s1)
    print(s2)
    s2.to_csv(out_dir + 'stage2_results.csv')
    #mac_high_tfs = list(tfs_high[tfs_high['celltype'] == 'macrophage']['tf'])
    #print(mac_high_tfs)
    #mac_high_links = find_grn_links(tf_list=mac_high_tfs, links1=links1, links2=links2)
    #print(mac_high_links)
    #mac_high_nodes_cut = ['STAT3', 'FOXP3', 'NFKB1', 'NFKB2', 'IRF1', 'FOSL2', 'CEBPB', 'ETV3']
    #mac_high_links_cut = [i for i in mac_high_links if (i[0] in mac_high_nodes_cut) and i[1] in mac_high_nodes_cut]
    #grn_viz(nodes=mac_high_tfs, edges=mac_high_links)
    #grn_viz(nodes=mac_high_tfs, edges=mac_high_links_cut)
    #mac_low_tfs = list(tfs_low[tfs_low['celltype'] == 'macrophage']['tf'])
    #print(mac_low_tfs)
    #mac_low_links = find_grn_links(tf_list=mac_low_tfs, links1=links1, links2=links2)
    #print(mac_low_links)
    #mac_low_nodes_cut = ['EBF1', 'LYL1', 'ATF1']
    #mac_low_links_cut = [i for i in mac_low_links if (i[0] in mac_low_nodes_cut) and i[1] in mac_low_nodes_cut]
    #grn_viz(nodes=mac_low_tfs, edges=mac_low_links_cut)

#  main function definition
def main():
    tic = time.perf_counter()
    toc = time.perf_counter()
    timer = toc - tic
    # display_results(ligand_scores=scores, runtime=timer)
    # save_results()

#  run main function
if __name__ == "__main__":
    main()
