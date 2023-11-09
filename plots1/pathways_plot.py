"""run with 'sc2' conda env

input files from 'gene_pathway_enrichment1.py'
"""
# import statements
import json
import pandas as pd
import math
import matplotlib.pyplot as plt
import numpy as np
import os
from textwrap import wrap

# function definitions
def get_children(pathway: dict):
    children = []
    if 'children' in pathway.keys():
        children = pathway['children']
        return children

def get_all_children(pathway: dict):
    """
    have growing queue of all children

    have growing list of pathways

    iterate through growing children queue until exhausted

    when operating on a given child dict, add pathway name to pathway list, get children, add children to chidren queue

    do this until child queue exhausted

    return pathway list
    """
    all_pathways = []
    all_pathways.append(pathway['name'])
    child_list = get_children(pathway)
    while len(child_list) > 0:
        child = child_list.pop(0)
        all_pathways.append(child['name'])
        c2 = get_children(child)
        if c2 is not None:
            for c in c2:
                child_list.append(c)

    return all_pathways

def add_path_cat(
        df: pd.DataFrame,
        #pathway_hierarchy_filename: str = '/Users/klockec/OneDrive - Oregon Health & Science University/data/analysis_files/p3/reactome_files/PathwayHierarchy.json'
        pathway_hierarchy_filename: str = '/Users/klockec/OneDrive - Oregon Health & Science University/data/analysis_files/p3/reactome_files/PathwayHierarchy_copy.json'
):
    """
    add pathway category -- for 'pathway' in input df, add another feature to tell which top-level pathway
    this pathway is a child of (or equal to, if pathway is top-level)

    need to check if query pathway is in total pathway set, otherwise while loop will never terminate
    """
    with open(pathway_hierarchy_filename) as json_file:
        data = json.load(json_file)
        top_dict = {}
        all_pathways = []
        for top_path in data:
            sub_pathways = get_all_children(top_path)
            for i in sub_pathways:
                all_pathways.append(i)
            top_dict[top_path['name']] = sub_pathways
        lookup = list(df['pathway'])
        top_level_pathways = []
        for l in lookup:
            if l in all_pathways:
                search = True
                while search:
                    for top in top_dict.keys():
                        if l in top_dict[top]:
                            top_level_pathways.append(top)
                            search = False
                            break
            else:
                print("pathway: " + l + " not in hierarchy file and will be labeled as 'unknown'")
                top_level_pathways.append('unknown')
        df['top_level_pathway'] = top_level_pathways

        return df

def process_df(
        df1: pd.DataFrame,
        df2: pd.DataFrame,
        df3: pd.DataFrame,
        threshold1: float,
        threshold2: float,
        pathway_categories: bool = True,
        sorted: bool = True
):
    """

    """
    df1['adjusted_p_value_t1'] = df1['adjusted_p_value']
    df2['adjusted_p_value_t2'] = df2['adjusted_p_value']
    df3['adjusted_p_value_v1'] = df3['adjusted_p_value']
    df = df1.merge(df2, on='pathway').merge(df3, on='pathway')
    df = df[['pathway', 'adjusted_p_value_t1', 'adjusted_p_value_t2', 'adjusted_p_value_v1']]
    df = df[( (df['adjusted_p_value_t1'] <= threshold1) | (df['adjusted_p_value_t2'] <= threshold1) ) | (df['adjusted_p_value_v1'] <= threshold1)]
    df['neg_logp_t1'] = [math.log(i, 10) * -1 for i in list(df['adjusted_p_value_t1'])]
    df['neg_logp_t2'] = [math.log(i, 10) * -1 for i in list(df['adjusted_p_value_t2'])]
    df['neg_logp_v1'] = [math.log(i, 10) * -1 for i in list(df['adjusted_p_value_v1'])]
    category = []
    pathway_level = []
    for row_num in range(df.shape[0]):
        t1p = df['adjusted_p_value_t1'].iloc[row_num]
        t2p = df['adjusted_p_value_t2'].iloc[row_num]
        v1p = df['adjusted_p_value_v1'].iloc[row_num]
        t1l = df['neg_logp_t1'].iloc[row_num]
        t2l = df['neg_logp_t2'].iloc[row_num]
        v1l = df['neg_logp_v1'].iloc[row_num]
        if t1p < threshold2 and t2p < threshold2 and v1p >= threshold2:
            category.append('tumor_specific')
            pathway_level.append(np.mean([t1l, t2l]))
        elif t1p >= threshold2 and t2p >= threshold2 and v1p < threshold2:
            category.append('viral_specific')
            pathway_level.append(v1l)
        elif t1p < threshold2 and t2p < threshold2 and v1p < threshold2:
            category.append('shared')
            pathway_level.append(np.mean([t1l, t2l, v1l]))
        else:
            category.append('skip')
            pathway_level.append(0)
    df['category'] = category
    df['pathway_level'] = pathway_level
    df = df[df['category'] != 'skip']
    if pathway_categories:
        df = add_path_cat(df)
    if sorted:
        df = df.sort_values(by='top_level_pathway')

    return df

def plot_pathways(
        df: pd.DataFrame,
        subset: str,
        xlim: float = 10,
        cutoff: float = 1,
        text_offset: float = 0.5,
        fig_width: float = 10,
        fig_height: float = 4,
        wrap_text_cutoff: int = 40,
        cmap: dict = {'Immune System': 'green', 'Cellular responses to stimuli': 'yellow', 'Disease': 'red', 'Cell-Cell communication': 'blue', 'Autophagy': 'purple', 'Signal Transduction': 'black'},
        show: bool = True,
        save: str = None,
        legend_location: str = 'upper right',
        title: bool = False,
        legend: bool = False
):
    df_subset = df[df['category'] == subset]
    df_subset = df_subset[df_subset['pathway_level'] >= cutoff]
    pathway = (df_subset['pathway'])
    pathway_enrich = {
        'tumor -- melanoma' : (df_subset['neg_logp_t1'], 'forestgreen'),
        'tumor -- BCC' : (df_subset['neg_logp_t2'], 'palegreen'),
        'viral -- HIV' : (df_subset['neg_logp_v1'], 'lightsalmon'),
    }
    x = np.arange(len(pathway))  # the label locations
    width = 0.25  # the width of the bars
    multiplier = 0
    fig, ax = plt.subplots(layout='constrained')
    for attribute, measurement in pathway_enrich.items():
        offset = width * multiplier
        rects = ax.barh(x + offset, measurement[0], width, label=attribute, color=measurement[1])
        #ax.bar_label(rects, padding=3)
        multiplier += 1

    # Add some text for labels, title and custom x-axis tick labels, etc.
    ax.set_xlabel('-log(adjusted p-value)')
    if title:
        ax.set_title('Pathway enrichment scores')
    ax.set_yticks(x + width, pathway)

    if legend:
        ax.legend(loc=legend_location)
    ax.set_xlim(0, xlim)
    fig.set_figwidth(fig_width)
    fig.set_figheight(fig_height)
    plt.subplots_adjust(left=text_offset)

    tick_labels = ax.get_yticklabels()
    for i, label in enumerate(tick_labels):
        color = df_subset.iloc[i]['top_level_pathway']
        wrapped_text = '\n'.join(wrap(label.get_text(), wrap_text_cutoff))
        tick_labels[i].set_text(wrapped_text)
        tick_labels[i].set_color(cmap[color])

    ax.set_yticklabels(tick_labels)

    red_line_position = 1.3  # Adjust this value as needed
    ax.axvline(x=red_line_position, color='red', label='Red Line', linewidth=1)

    if save is not None:
        plt.savefig(save)

    if show:
        plt.show(block=True)

# setup
results_load_dir = "/Users/klockec/OneDrive - Oregon Health & Science University/data/analysis_files/p3/mod1/deg_gsea/"
plots_store = "/Users/klockec/OneDrive - Oregon Health & Science University/data/analysis_files/p3/mod1/deg_gsea/pathway_plots/"
li_mac_high_filename = 'li_deg_gsea_bh_macrophage_high_all.csv'
yost_mac_high_filename = 'yost_deg_gsea_bh_macrophage_high_all.csv'
wang_mac_high_filename = 'wang_deg_gsea_bh_macrophage_high_all.csv'
li_cd8_high_filename = 'li_deg_gsea_bh_CD8_T_cells_high_all.csv'
yost_cd8_high_filename = 'yost_deg_gsea_bh_CD8_T_cells_high_all.csv'
wang_cd8_high_filename = 'wang_deg_gsea_bh_CD8_T_cells_high_all.csv'
li_cd8_low_filename = 'li_deg_gsea_bh_CD8_T_cells_low_all.csv'
yost_cd8_low_filename = 'yost_deg_gsea_bh_CD8_T_cells_low_all.csv'
wang_cd8_low_filename = 'wang_deg_gsea_bh_CD8_T_cells_low_all.csv'

os.chdir(plots_store)  # set the working directory (to save the plots)

# load the data
li_mac_high = pd.read_csv(results_load_dir + li_mac_high_filename)
yost_mac_high = pd.read_csv(results_load_dir + yost_mac_high_filename)
wang_mac_high = pd.read_csv(results_load_dir + wang_mac_high_filename)
li_cd8_high = pd.read_csv(results_load_dir + li_cd8_high_filename)
yost_cd8_high = pd.read_csv(results_load_dir + yost_cd8_high_filename)
wang_cd8_high = pd.read_csv(results_load_dir + wang_cd8_high_filename)
li_cd8_low = pd.read_csv(results_load_dir + li_cd8_low_filename)
yost_cd8_low = pd.read_csv(results_load_dir + yost_cd8_low_filename)
wang_cd8_low = pd.read_csv(results_load_dir + wang_cd8_low_filename)

custom_cmap1 = {'Immune System': 'green',
               'Cellular responses to stimuli': 'yellow',
               'Disease': 'red',
               'Cell-Cell communication': 'blue',
               'Autophagy': 'purple',
               'Signal Transduction': 'black',
               'Cell Cycle': 'orange',
               'DNA Replication': 'magenta',
               'Metabolism of RNA': 'chocolate',
               'Metabolism of proteins': 'darkslategray',
               'Gene expression (Transcription)': 'lavender',
               'Metabolism': 'olive',
               'DNA Repair': 'sienna',
               'Vesicle-mediated transport': 'turquoise',
               'Developmental Biology': 'yellowgreen',
               'Programmed Cell Death': 'plum',
               'unknown': 'salmon'
               }

other = 'black'
custom_cmap2 = {'Immune System': 'dodgerblue',
               'Cellular responses to stimuli': 'darkkhaki',
               'Disease': other,
               'Cell-Cell communication': 'navy',
               'Autophagy': other,
               'Signal Transduction': 'purple',
               'Cell Cycle': 'red',
               'DNA Replication': other,
               'Metabolism of RNA': other,
               'Metabolism of proteins': other,
               'Gene expression (Transcription)': other,
               'Metabolism': other,
               'DNA Repair': other,
               'Vesicle-mediated transport': other,
               'Developmental Biology': 'gainsboro',
               'Programmed Cell Death': other,
               'unknown': other
               }

custom_cmap = custom_cmap2

df_mac_high = process_df(df1=li_mac_high, df2=yost_mac_high, df3=wang_mac_high, threshold1=0.05, threshold2=0.05)

#threshold_new = 0.03
#threshold_new = 0.04  # used for viral and shared -- macrophage high
#threshold_tumor_specific = 0.05  # used for tumor-specific (updated version)
#df = df_mac_high
#df = df[( (df['adjusted_p_value_t1'] <= threshold_new) | (df['adjusted_p_value_t2'] <= threshold_new) ) | (df['adjusted_p_value_v1'] <= threshold_new)]
#remove = ['Uptake and function of diphtheria toxin', 'Chaperone Mediated Autophagy', 'PERK regulates gene expression', 'Unfolded Protein Response (UPR)']
#df = df[-df['pathway'].isin(remove, )]
#df_mac_high_mod2 = df_mac_high[df_mac_high['pathway'] != 'Disease']  # remove 'Disease' pathway

#plot_pathways(df=df_mac_high, subset='tumor_specific', xlim=10, text_offset=0.55, fig_width=10, fig_height=4, cmap=custom_cmap, wrap_text_cutoff=40, save='mac_high_tumor_specific_pathways_v3.pdf')  # plot macrophage high, tumor-specific pathways
#plot_pathways(df=df, subset='viral_specific', xlim=6, legend_location='lower right', cmap=custom_cmap, wrap_text_cutoff=45, show=True, save='mac_high_viral_specific_pathways_v4.pdf')  # plot macrophage high, viral-specific pathways
#plot_pathways(df=df_mac_high_mod2, subset='shared', xlim=26, cmap=custom_cmap, wrap_text_cutoff=100, save='mac_high_shared_pathways_v2.pdf')  # plot macrophage high, shared pathways

# unabridged figure for supplemental -- macrophage high
#plot_pathways(df=df_mac_high, subset='tumor_specific', xlim=10, text_offset=0.55, fig_width=10, fig_height=4, cmap=custom_cmap, wrap_text_cutoff=40, save='mac_high_tumor_specific_pathways_uncut.pdf')  # plot macrophage high, tumor-specific pathways
#plot_pathways(df=df_mac_high, subset='viral_specific', xlim=6, legend_location='lower right', cmap=custom_cmap, wrap_text_cutoff=45, fig_height=6.5, save='mac_high_viral_specific_pathways_uncut_v2.pdf')  # plot macrophage high, viral-specific pathways
#plot_pathways(df=df_mac_high, subset='shared', xlim=26, cmap=custom_cmap, wrap_text_cutoff=100, save='mac_high_shared_pathways_uncut.pdf')  # plot macrophage high, shared pathways

#df_cd8_high = process_df(df1=li_cd8_high, df2=yost_cd8_high, df3=wang_cd8_high, threshold1=0.05, threshold2=0.05)
#plot_pathways(df=df_cd8_high, subset='tumor_specific', xlim=80, cutoff=20, cmap=custom_cmap, fig_height=7, save='cd8_high_tumor_specific_pathways_v3.pdf')  # plot cd8 high, tumor-specific pathways
#plot_pathways(df=df_cd8_high, subset='viral_specific', xlim=7, cutoff=2, cmap=custom_cmap, fig_height=9.5, wrap_text_cutoff=60, save='cd8_high_viral_specific_pathways_v2.pdf')  # plot cd8 high, viral-specific pathways
#plot_pathways(df=df_cd8_high, subset='shared', xlim=42, cutoff=2, cmap=custom_cmap, fig_height=8.5, wrap_text_cutoff=50, save='cd8_high_shared_pathways_v3.pdf')  # plot cd8 high, shared pathways

df_cd8_low = process_df(df1=li_cd8_low, df2=yost_cd8_low, df3=wang_cd8_low, threshold1=0.05, threshold2=0.05)
#plot_pathways(df=df_cd8_low, subset='tumor_specific', xlim=15, cmap=custom_cmap, save='cd8_low_tumor_specific_pathways.pdf')  # plot cd8 low, tumor-specific pathways
#plot_pathways(df=df_cd8_low, subset='viral_specific', xlim=30, cutoff=10, cmap=custom_cmap, fig_height=22, wrap_text_cutoff=70, save='cd8_low_viral_specific_pathways_v4.pdf')  # plot cd8 low, viral-specific pathways
#plot_pathways(df=df_cd8_low, subset='shared', xlim=58, cutoff=20, cmap=custom_cmap, fig_height=10, wrap_text_cutoff=60, save='cd8_low_shared_pathways_v3.pdf')  # plot cd8 low, shared pathways

# cd8 low viral-specific figure is still a mess; try cutting it in half
#df_cd8_low_part1 = df_cd8_low
#df_cd8_low_part2 = df_cd8_low

"""
['Autophagy' 'Cell Cycle' 'Cell-Cell communication'
 'Cellular responses to stimuli' 'Chromatin organization'
 'Circadian Clock' 'DNA Repair' 'DNA Replication' 'Developmental Biology'
 'Disease' 'Gene expression (Transcription)' 'Hemostasis' 'Immune System'
 'Metabolism' 'Metabolism of RNA' 'Metabolism of proteins'
 'Programmed Cell Death' 'Protein localization' 'Signal Transduction'
 'Transport of small molecules' 'Vesicle-mediated transport' 'unknown']
"""
df = df_cd8_low[df_cd8_low['category'] == 'viral_specific']
#print(df['top_level_pathway'].unique())
#df1 = df[df['top_level_pathway'] == 'Cell Cycle']
#plot_pathways(df=df1, subset='viral_specific', xlim=15, cutoff=10, cmap=custom_cmap, fig_height=10, wrap_text_cutoff=70, save='cd8_low_viral_specific_pathways_cell_cycle.pdf')
#pathways2 = ['Immune System', 'Autophagy', 'Cell-Cell communication', 'Cellular responses to stimuli', 'Chromatin organization', 'Circadian Clock', 'Circadian Clock', 'DNA Repair', 'DNA Replication', 'Developmental Biology', 'Disease', 'Gene expression (Transcription)', 'Hemostasis', 'Metabolism of RNA', 'Metabolism of proteins']
#df2 = df[df['top_level_pathway'].isin(pathways2)]
#plot_pathways(df=df2, subset='viral_specific', xlim=30, cutoff=10, cmap=custom_cmap, fig_height=8, fig_width=12, wrap_text_cutoff=75, save='cd8_low_viral_specific_pathways_pathways2.pdf')
pathways3 = ['Metabolism', 'Programmed Cell Death', 'Protein localization', 'Signal Transduction', 'Transport of small molecules', 'Vesicle-mediated transport', 'unknown']
df3 = df[df['top_level_pathway'].isin(pathways3)]
plot_pathways(df=df3, subset='viral_specific', xlim=30, cutoff=10, cmap=custom_cmap, fig_height=8, fig_width=12, wrap_text_cutoff=70, save='cd8_low_viral_specific_pathways_pathways3.pdf')
