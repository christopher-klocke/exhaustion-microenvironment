"""
conda create --name plot_graph

conda activate plot_graph

conda install -c anaconda pandas

conda install -c conda-forge igraph

conda install -c anaconda networkx

conda install -c conda-forge matplotlib



conda remove networkx

pip install networkx
"""

## import statements
import pandas as pd

## setup
li_regs_gmt = '/Users/klockec/OneDrive - Oregon Health & Science University/data_backup_rws07890/data_september2022_backup/pyscenic_runs/p3/li_run/ctx_output1.gmt'
tfs_check_file = '/Users/klockec/OneDrive - Oregon Health & Science University/data/analysis_files/p3/manuscript_figs/grn_reconstruction/sig_tf_overlap_li_wang_cd8_v1.csv'
graph_outfile = '/Users/klockec/OneDrive - Oregon Health & Science University/data/analysis_files/p3/manuscript_figs/grn_reconstruction/cd8_graph1.csv'
graph_plotted_file = '/Users/klockec/OneDrive - Oregon Health & Science University/data/analysis_files/p3/manuscript_figs/grn_reconstruction/cd8_graph1_plot_v2.png'

## load the data
tfs_check = list(pd.read_csv(tfs_check_file)['tf'])
regs_dict = {}

## process the data
with open(li_regs_gmt, 'r') as infile:
    for line in infile:
        line = line.rstrip()
        line = line.split('\t')
        name = line[0]
        desc = line[1]
        targets = line[2:]
        tf = desc.split(',')[0][3:]
        if tf in tfs_check:
            #targets_keep = [i for i in targets if i in tfs_check]  # KEEP AUTOREGULATION
            targets_keep = [i for i in targets if i in tfs_check and i != tf]  # REMOVE AUTOREGULATION
            if len(targets_keep) != 0:
                regs_dict[tf] = targets_keep

"""
with open(graph_outfile, 'w') as write_out:
    for tf in regs_dict.keys():
        for target in regs_dict[tf]:
            line = tf + ',' + target + '\n'
            write_out.write(line)
"""

import networkx as nx
import matplotlib.pyplot as plt

import numpy as np
import matplotlib.colors as mcolors
import matplotlib.cm as cm

G = nx.from_dict_of_lists(regs_dict)

## for threshold graph
from networkx.algorithms.threshold import find_threshold_graph
G = find_threshold_graph(G)

pos = nx.spring_layout(G, k=2.5)
#pos = nx.circular_layout(G)
#pos = nx.spiral_layout(G)

#nx.draw(G, pos=pos, node_color='r', edge_color='b')
labels = nx.draw_networkx_labels(G, pos=pos, font_size=6)
edges = nx.draw_networkx_edges(G, pos, edge_color='green', arrows=True, arrowstyle='-|>')

deg_centrality = nx.degree_centrality(G)

cent = np.fromiter(deg_centrality.values(), float)
sizes = cent / np.max(cent) * 200
normalize = mcolors.Normalize(vmin=cent.min(), vmax=cent.max())
colormap = cm.viridis

scalarmappaple = cm.ScalarMappable(norm=normalize, cmap=colormap)
scalarmappaple.set_array(cent)

plt.colorbar(scalarmappaple)
nx.draw(G, pos, node_size=sizes, node_color=sizes, cmap=colormap)

plt.show()
plt.savefig(graph_plotted_file)
