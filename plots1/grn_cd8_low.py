"""run with 'network' conda env """
# import statements
import sys
import pandas as pd
import networkx as nx
import grns_plot as gp

# setup
dataset1_results = "/Users/klockec/OneDrive - Oregon Health & Science University/data/analysis_files/p3/mod1/li_parallel_tfs_run4/combined_bh.csv"
dataset1_links = "/Users/klockec/OneDrive - Oregon Health & Science University/data_backup_rws07890/data_january2023_backup/data/analysis_files/p3/pyscenic_runs/p3/li_run/ctx_output1.gmt"
dataset2_results = "/Users/klockec/OneDrive - Oregon Health & Science University/data/analysis_files/p3/mod1/yost_parallel_tfs_run4/combined_bh.csv"
dataset2_links = "/Users/klockec/OneDrive - Oregon Health & Science University/data_backup_rws07890/data_january2023_backup/data/analysis_files/p3/pyscenic_runs/p3/yost_bcc_run1/ctx_output1.gmt"
dataset3_results = "/Users/klockec/OneDrive - Oregon Health & Science University/data/analysis_files/p3/mod1/wang_parallel_tfs_run4/combined_bh.csv"
dataset3_links = "/Users/klockec/OneDrive - Oregon Health & Science University/data_backup_rws07890/data_january2023_backup/data/analysis_files/p3/pyscenic_runs/p3/wang_run/ctx_output1.gmt"
out_dir1 = "/Users/klockec/OneDrive - Oregon Health & Science University/data/analysis_files/p3/mod1/li_grn_graphs2/"
out_dir2 = "/Users/klockec/OneDrive - Oregon Health & Science University/data/analysis_files/p3/mod1/yost_grn_graphs2/"
out_dir3 = "/Users/klockec/OneDrive - Oregon Health & Science University/data/analysis_files/p3/mod1/wang_grn_graphs2/"
out_dir_shared = "/Users/klockec/OneDrive - Oregon Health & Science University/data/analysis_files/p3/mod1/tumor_tumor_grn_graphs/"
out_dir_tumor_viral = "/Users/klockec/OneDrive - Oregon Health & Science University/data/analysis_files/p3/mod1/tumor_viral_grn_graphs/"
out_dir_tumor_not_viral = "/Users/klockec/OneDrive - Oregon Health & Science University/data/analysis_files/p3/mod1/tumor_not_viral_grn_graphs/"
out_dir_viral_not_tumor = "/Users/klockec/OneDrive - Oregon Health & Science University/data/analysis_files/p3/mod1/viral_not_tumor_grn_graphs/"

# load the data
results1 = pd.read_csv(dataset1_results)
links1 = gp.gmt_to_dict(dataset1_links)
results2 = pd.read_csv(dataset2_results)
links2 = gp.gmt_to_dict(dataset2_links)
results3 = pd.read_csv(dataset3_results)
links3 = gp.gmt_to_dict(dataset3_links)
pick_celltype = 'CD8_T_cells'
pick_dir = 'low'
out_dir_agg = '/Users/klockec/OneDrive - Oregon Health & Science University/data/analysis_files/p3/mod1/grn_agg_fig/'
res1 = gp.df_cut1(df=results1, celltype=pick_celltype, dir=pick_dir)
res2 = gp.df_cut1(df=results2, celltype=pick_celltype, dir=pick_dir)
res3 = gp.df_cut1(df=results3, celltype=pick_celltype, dir=pick_dir)
tfs1 = list(res1["tf"])
tfs2 = list(res2["tf"])
tfs3 = list(res3["tf"])
all_tfs = []
all_tfs += tfs1
all_tfs += tfs2
all_tfs += tfs3
all_tfs = set(all_tfs)
df_rows = []
for tf in all_tfs:
    if tf in tfs1 and tf in tfs2 and tf in tfs3:  # significant in all 3 datasets
        df_rows.append((tf, 'tumor_viral_both'))
    elif tf in tfs1 and tf in tfs2 and tf not in tfs3:  # significant in both tumor but not viral
        df_rows.append((tf, 'tumor_specific'))
    elif (tf in tfs1 or tf in tfs2) and tf in tfs3:  # significant in one tumor and viral
        df_rows.append((tf, 'tumor_viral_split'))
    elif tf not in tfs1 and tf not in tfs2 and tf in tfs3:  #significant in just viral
        df_rows.append((tf, 'viral_specific'))

df = pd.DataFrame(df_rows, columns=['tf', 'annotation'])
df = df[df["annotation"] != 'tumor_viral_split']  # remove 'tumor_viral_split' for now

# remove 'tumor_viral_split'
tfs1 = [i for i in tfs1 if i in list(df["tf"])]
tfs2 = [i for i in tfs2 if i in list(df["tf"])]
tfs3 = [i for i in tfs3 if i in list(df["tf"])]

# cut links down by significant TFs
links1 = gp.links_cut(tf_list=tfs1, links=links1, self_links=True)
links2 = gp.links_cut(tf_list=tfs2, links=links2, self_links=True)
links3 = gp.links_cut(tf_list=tfs3, links=links3, self_links=True)

lower_score_threshold = 1
higher_score_threshold = 30
file_suffix = 'mod1'

li_links_scores_filename = "/Users/klockec/OneDrive - Oregon Health & Science University/data_backup_rws07890/data_january2023_backup/data/analysis_files/p3/pyscenic_runs/p3/li_run/grnboost_out.tsv"
yost_links_scores_filename = "/Users/klockec/OneDrive - Oregon Health & Science University/data_backup_rws07890/data_january2023_backup/data/analysis_files/p3/pyscenic_runs/p3/yost_bcc_run1/grnboost_out.tsv"
wang_links_scores_filename = "/Users/klockec/OneDrive - Oregon Health & Science University/data_backup_rws07890/data_january2023_backup/data/analysis_files/p3/pyscenic_runs/p3/wang_run/grnboost_out.tsv"
li_links_scores = pd.read_csv(li_links_scores_filename, sep='\t')
yost_links_scores = pd.read_csv(yost_links_scores_filename, sep='\t')
wang_links_scores = pd.read_csv(wang_links_scores_filename, sep='\t')
links1_update = gp.cut_by_scores(links=links1, scores_df=li_links_scores, score_threshold=lower_score_threshold, edge_weights=True, as_df=True)
links2_update = gp.cut_by_scores(links=links2, scores_df=yost_links_scores, score_threshold=lower_score_threshold, edge_weights=True, as_df=True)
links3_update = gp.cut_by_scores(links=links3, scores_df=wang_links_scores, score_threshold=lower_score_threshold, edge_weights=True, as_df=True)
datasets1 = {'t1': links1_update, 't2': links2_update, 'v1': links3_update, 'threshold': higher_score_threshold}
links_agg1 = gp.links_agg(**datasets1)
links_agg1.to_csv((out_dir_agg + 'cd8_low_annotations_agg_' + file_suffix + '.csv'), index=False)
df.to_csv((out_dir_agg + "cd8_low_tf_annotation_" + file_suffix + ".csv"), index=False)  # node annotation file -- save

# save updated .graphml object
links_agg1_tograph = links_agg1[['TF', 'target', 'importance']]
links_agg1_tograph.columns = ['source', 'target', 'importance']
G = nx.from_pandas_edgelist(links_agg1_tograph, 'source', 'target', create_using=nx.DiGraph())
nx.write_graphml(G=G, path=(out_dir_agg + "cd8_low_agg_" + file_suffix + ".graphml"))

