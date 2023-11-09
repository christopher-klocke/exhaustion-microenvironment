"""generate fake violin plots for overview figure """

# import statements
import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

## file paths
plots_store = '/Users/klockec/OneDrive - Oregon Health & Science University/data/analysis_files/p3/manuscript_figs/violins2/'
fig_name = 'synthetic_violin1.pdf'

# setup
random1 = 1
fig_size = (10, 10)
font_scale = 1
os.chdir(plots_store)

# generate synthetic data
list1 = [1, 2, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 6, 6, 7, 8, 15, 20, 52, 52, 53, 58]
list2 = [1, 2, 3, 3, 3, 4, 5, 5, 5, 6, 7, 10, 12, 14, 16, 50, 50, 51, 51, 51, 51, 52, 53, 53, 55, 56, 56, 62]
list3 = [1, 2, 3, 4, 5, 6, 50, 50, 51, 51, 51, 53, 53, 53, 53, 55, 55, 55, 55, 55, 56, 56, 56, 56, 56, 58, 60, 65]
df = pd.DataFrame({'sample1': list1, 'sample2': list2, 'sample3': list3})

# generate visualization
plt.figure(figsize=fig_size)
sns.set_theme(style="white", font_scale=font_scale)
ax = sns.violinplot(data=df)
#plt.show(block=True)
plt.savefig(fig_name)
plt.close()
