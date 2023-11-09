""" run with 'sc2' conda env

create plots comparing samples from patients with responder / non-resonder status

"""
## import statements
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt

## setup
scores_filename = '/Users/klockec/OneDrive - Oregon Health & Science University/data/analysis_files/p3/mod1/yost_r_nr_scores_updated.csv'
#mac_high_tfs_shared_filename = '/Users/klockec/OneDrive - Oregon Health & Science University/data/analysis_files/p3/mod1/tumor_tumor_overlap_run5/mac_high_tfs_shared.txt'
#outfile1 = '/Users/klockec/OneDrive - Oregon Health & Science University/data/analysis_files/p3/manuscript_figs/r_nr/yost_scatter2.pdf'
outfile3 = '/Users/klockec/OneDrive - Oregon Health & Science University/data/analysis_files/p3/manuscript_figs/r_nr/yost_scatter3.pdf'

# load the data
df = pd.read_csv(scores_filename)

#param_set = [200, 80, 3, 4.5, 0.8, 100]
param_set = [200, 80, 2, 5, 0.8, 100]

## visualize the data
mpl.rcParams['figure.dpi'] = param_set[0]
#point_size = 60
point_size = param_set[1]
fig = plt.figure(figsize = (param_set[2], param_set[3]))
sns.set_theme(style="whitegrid", palette="muted", font_scale=param_set[4])
ax = sns.swarmplot(data=df, x="response", y="sample score", sizes=(point_size, point_size))  # Draw a categorical scatterplot to show each observation
# hue="tumor type",
ax.set(xlabel="response to immunotherapy")
ax.set(ylabel="sample-level exhaustion score")
plt.tight_layout(pad=param_set[5])

plt.show()
#ax.figure.savefig(outfile1)
#ax.figure.savefig(outfile3)
ax.figure.savefig(outfile3, bbox_inches="tight")

#fig.savefig('/Users/klockec/OneDrive - Oregon Health & Science University/data/analysis_files/p3/manuscript_figs/r_nr/yost_scatter1.png')


#pandas.DataFrame.plot. df.plot(kind='box', subplots=True, layout=(7, 2), figsize=(7, 10)); plt.tight_layout()


exit()

#  PART TWO -- MAC HIGH TF PANEL VS. R / NR
mac_high_tfs = pd.read_csv(mac_high_tfs_shared_filename, header=None)
mac_high_tfs = list(mac_high_tfs.iloc[:, 0])
adata_yost_filename = "/Users/klockec/OneDrive - Oregon Health & Science University/data_backup_rws07890/data_january2023_backup/data/yost_data/yost_bcc_full_processed_dd_regs_mod2.h5ad"
adata = sc.read_h5ad(adata_yost_filename)
adata_macs = adata[adata.obs['celltype'] == 'macrophage']


sc.tl.score_genes(adata=adata_macs, gene_list=mac_high_tfs)

samples = []
scores = []
for sample in set(adata_macs.obs['sample']):
    adata_sample = adata_macs[adata_macs.obs['sample'] == sample]
    mean_score = np.mean(adata_sample.obs['score'])
    samples.append(sample)
    scores.append(mean_score)

df2 = pd.DataFrame(samples, columns=['sample ID'])
df2['scores'] = scores

df_full = pd.merge(df, df2, how='outer', left_on='sample ID', right_on='sample ID')



print(df_full)



## visualize the data
mpl.rcParams['figure.dpi'] = 200
#point_size = 60
point_size = 200
fig = plt.figure(figsize = (6, 7))
sns.set_theme(style="whitegrid", palette="muted")
ax = sns.swarmplot(data=df_full, x="response", y="scores", sizes=(point_size, point_size))  # Draw a categorical scatterplot to show each observation
# hue="tumor type",
ax.set(xlabel="response to immunotherapy")
ax.set(ylabel="macrophage exhaustion-associated TFs")
plt.tight_layout(pad=100)

plt.show()
#ax.figure.savefig('/Users/klockec/OneDrive - Oregon Health & Science University/data/analysis_files/p3/manuscript_figs/r_nr/yost_scatter1.png')
#fig.savefig('/Users/klockec/OneDrive - Oregon Health & Science University/data/analysis_files/p3/manuscript_figs/r_nr/yost_scatter1.png')
