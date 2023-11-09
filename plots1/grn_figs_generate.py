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

# setup
#dataset1_results = "/Users/klockec/OneDrive - Oregon Health & Science University/data/analysis_files/p3/mod1/li_parallel_tfs_run3/combined_bh.csv"
dataset1_results = "/Users/klockec/OneDrive - Oregon Health & Science University/data/analysis_files/p3/mod1/li_parallel_tfs_run4/combined_bh.csv"
dataset1_links = "/Users/klockec/OneDrive - Oregon Health & Science University/data_backup_rws07890/data_january2023_backup/data/analysis_files/p3/pyscenic_runs/p3/li_run/ctx_output1.gmt"




mac_high_graph_path = "/Users/klockec/OneDrive - Oregon Health & Science University/data/analysis_files/p3/mod1/li_grn_graphs/mac_high_grn1.graphml"
mac_high_graph_path_cut = "/Users/klockec/OneDrive - Oregon Health & Science University/data/analysis_files/p3/mod1/li_grn_graphs/mac_high_grn2.graphml"
mac_high_graph_path_cut2 = "/Users/klockec/OneDrive - Oregon Health & Science University/data/analysis_files/p3/mod1/li_grn_graphs/mac_high_grn_v2.graphml"
mac_high_graph_path_cut3 = "/Users/klockec/OneDrive - Oregon Health & Science University/data/analysis_files/p3/mod1/li_grn_graphs/mac_high_grn_v3.graphml"
mac_high_graph_path_cut4 = "/Users/klockec/OneDrive - Oregon Health & Science University/data/analysis_files/p3/mod1/li_grn_graphs/mac_high_grn_v4.graphml"
mac_high_graph_path_cut5 = "/Users/klockec/OneDrive - Oregon Health & Science University/data/analysis_files/p3/mod1/li_grn_graphs/mac_high_grn_v5.graphml"
mac_low_graph_path_cut1 = "/Users/klockec/OneDrive - Oregon Health & Science University/data/analysis_files/p3/mod1/li_grn_graphs/mac_low_grn_v1.graphml"
mac_high_graph_path_cut6 = "/Users/klockec/OneDrive - Oregon Health & Science University/data/analysis_files/p3/mod1/li_grn_graphs/mac_high_grn_v6.graphml"

# load the data
results1 = pd.read_csv(dataset1_results)
links1 = gmt_to_dict(dataset1_links)

#sig_threshold = 0.05
#sig_threshold = 0.0000001
sig_threshold = 0.0000000001
#sig_threshold = 0.000000001


mac_tfs = results1[results1['celltype'] == 'macrophage']
mac_high_tfs = mac_tfs[mac_tfs['direction'] == 'high']
mac_high_tfs_sig = list(mac_high_tfs[mac_high_tfs['adjusted_p_value'] < sig_threshold]['tf'])


mac_low_tfs = mac_tfs[mac_tfs['direction'] == 'low']
mac_low_tfs_sig = list(mac_low_tfs[mac_low_tfs['adjusted_p_value'] < sig_threshold]['tf'])




tfs_subset = mac_high_tfs_sig  # UPDATE HERE
#tfs_subset = mac_low_tfs_sig


links = links_cut(tf_list=tfs_subset, links=links1)

print(str(len(tfs_subset)))
print(str(len(links)))

exit()

#print(links)


#https://stackoverflow.com/questions/22391433/count-the-frequency-that-a-value-occurs-in-a-dataframe-column
out_deg = [i[0] for i in links]
in_deg = [i[1] for i in links]
df_degs = pd.DataFrame(out_deg, columns=['out_degree'])
df_degs['in_degree'] = in_deg

df_out_deg = pd.DataFrame(df_degs['out_degree'].value_counts())
#df_in_deg = df_degs['in_degree'].value_counts()
df_out_deg.columns = ['out_degree']
df_out_deg['tf'] = df_out_deg.index

print(df_out_deg)

out_deg_min = 2
df_out_deg_cut = df_out_deg[df_out_deg['out_degree'] >= out_deg_min]

#links2 = [i for i in links if (i[0] in list(df_out_deg_cut['tf']) and (i[1] in list(df_out_deg_cut['tf'])))]  # cut again  #######
links2 = [i for i in links if i[0] in list(df_out_deg_cut['tf'])]  # cut again


print(str(len(links2)))



#print(mac_high_tfs_sig)
#print(links)


G = nx.from_edgelist(links2)  ##########
#G = nx.DiGraph()
#G.add_edges_from(links)

#nx.write_graphml(G=G, path=mac_high_graph_path)  ###########
#nx.write_graphml(G=G, path=mac_high_graph_path_cut6)  ###########

