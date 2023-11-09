""" run with 'plotting' conda env

create venn diagrams

https://www.geeksforgeeks.org/how-to-create-and-customize-venn-diagrams-in-python/
"""

# import statements
import os
import pandas as pd
from matplotlib_venn import venn2
from matplotlib import pyplot as plt

# setup
df_filename = '/Users/klockec/OneDrive - Oregon Health & Science University/data/analysis_files/p3/mod1/tumor_viral_overlap_li_wang_run3/stage1_results.csv'
plots_store = '/Users/klockec/OneDrive - Oregon Health & Science University/data/analysis_files/p3/manuscript_figs/venns_tumor_viral/'
os.chdir(plots_store)

# load the data
df = pd.read_csv(df_filename)

# create visualizations
for res in range(df.shape[0]):
    res_row = df.iloc[res]
    shared = res_row['shared_tfs']
    left = res_row['just_A_tfs']
    right = res_row['just_B_tfs']
    direction = res_row['direction']
    celltype = res_row['celltype']
    adj_p_value = res_row['adjusted_p_value']
    plt.figure(figsize=(2.2, 2.2))
    venn2(subsets=(left, right, shared), set_labels=('melanoma', 'HIV'))
    #plt.title(celltype + ', ' + direction + ': ' + str(adj_p_value))
    #plt.title(str(round(adj_p_value, 2)))
    #plt.show()
    plt.savefig('venn_' + celltype + '_' + direction + '.pdf')
    plt.close()

