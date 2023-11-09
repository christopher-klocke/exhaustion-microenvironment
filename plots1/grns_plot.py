"""run with 'network' conda env """
# import statements
import sys
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt

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
        tf_df: pd.DataFrame,
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

def dict_to_set(in_dict: dict):
    """
    generate set including all genes listed in keys and values (list) of TF/target dict
    """
    out_list = []
    for tf in in_dict.keys():
        out_list.append(tf)
        for target in in_dict[tf]:
            out_list.append(target)
    out_set = set(out_list)

    return out_set

def dict_aggregate(dicts: list):
    agg_dict = {}
    for d in dicts:
        for k in d.keys():
            if k in agg_dict.keys():
                for v in d[k]:
                    agg_dict[k].append(v)
            else:
                agg_dict[k] = d[k]

    for k in agg_dict.keys():
        agg_dict[k] = set(agg_dict[k])

    return agg_dict

def df_cut1(
        df: pd.DataFrame,
        celltype: str,
        dir: str,
        significant: bool = True
):
    """
    subset pandas dataframe by celltype and direction
    """
    df = df[df["celltype"] == celltype]
    df = df[df["direction"] == dir]
    if significant:
        df = df[df["reject"] == True]

    return df

def cut_by_scores(
        links: list,
        scores_df: pd.DataFrame,
        score_threshold: float = 10,
        edge_weights: bool = False,
        as_df: bool = False
):
    """
    cut links by scores of correlation relationships from GRNBoost
    """
    scores_df = scores_df[scores_df['importance'] > score_threshold]  # subset grnboost scores df by score threshold
    links_updated = []  # initialize output
    links_df = pd.DataFrame(links, columns=['tf', 'target'])
    unique_tfs = set(links_df['tf'])
    for tf in unique_tfs:
        score_df_cut = scores_df[scores_df['TF'] == tf]  # subset full grnboost df to regulator TF of interest
        links_df_cut = links_df[links_df['tf'] == tf]
        for target in links_df_cut['target']:
            if target in list(score_df_cut['target']):
                if edge_weights:
                    weight = float(score_df_cut.loc[score_df_cut['target'] == target, 'importance'])
                    links_updated.append((tf, target, weight))
                else:
                    links_updated.append((tf, target))
    if as_df:
        links_updated = pd.DataFrame(links_updated, columns=['TF', 'target', 'importance'])

    return links_updated

def links_agg(
        t1: pd.DataFrame,
        t2: pd.DataFrame,
        v1: pd.DataFrame,
        threshold: float
):
    """
    already have importance scores for each TF/target pair for GRN edges; need to add one of three categories (tumor, viral, both)

    do this with joins and anti joins -- see code from: 'gene_pathway_overlap.py'

    input is list of pandas dataframes, each of which is output from 'cut_by_scores()'
    """
    join_on = ['TF', 'target']
    tumor_both = t1.merge(t2, how='inner', left_on=join_on, right_on=join_on)
    tumor_only = tumor_both.merge(v1, on=join_on, how='left', indicator=True)
    tumor_only = tumor_only[tumor_only['_merge'] == 'left_only']
    tumor_only = tumor_only.drop('_merge', axis=1)
    tumor_only['importance'] = tumor_only[['importance_x', 'importance_y']].max(axis=1)
    viral_only = v1.merge(tumor_both, on=join_on, how='left', indicator=True)
    viral_only = viral_only[viral_only['_merge'] == 'left_only']
    viral_only = viral_only.drop('_merge', axis=1)
    shared = tumor_both.merge(v1, how='inner', left_on=join_on, right_on=join_on)
    shared['importance'] = shared[['importance_x', 'importance_y']].max(axis=1)
    tumor_only = tumor_only[['TF', 'target', 'importance']]
    viral_only = viral_only[['TF', 'target', 'importance']]
    shared = shared[['TF', 'target', 'importance']]
    tumor_only['category'] = ['tumor_only'] * tumor_only.shape[0]
    viral_only['category'] = ['viral_only'] * viral_only.shape[0]
    shared['category'] = ['shared'] * shared.shape[0]
    dfs = [tumor_only, viral_only, shared]
    dfs_concat = pd.concat(dfs, ignore_index=True)
    dfs_concat['name'] = dfs_concat['TF'] + ' (-) ' + dfs_concat['target']
    dfs_concat = dfs_concat[dfs_concat['importance'] >= threshold]

    return dfs_concat

def main():
    pass

if __name__ == "__main__":
    main()

