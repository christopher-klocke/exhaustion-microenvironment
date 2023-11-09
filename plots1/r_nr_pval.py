""" run with 'sc2' conda env

run statistical test on the relationship between sample-level exhaustion score and patient responder / non-resonder status

"""
# import statements
import pandas as pd
import scipy.stats as stats

# setup
scores_filename = '/Users/klockec/OneDrive - Oregon Health & Science University/data/analysis_files/p3/mod1/yost_r_nr_scores_updated.csv'

# load the data
df = pd.read_csv(scores_filename)
df.columns = ['sample_id', 'sample_score', 'response', 'tumor_type']
df = df[df['tumor_type'] == 'bcc']

# ruh the test
responder_scores = df[df['response'] == 'Yes ']['sample_score']
non_responder_scores = df[df['response'] == 'No ']['sample_score']
pval = stats.mannwhitneyu(responder_scores, non_responder_scores, alternative='two-sided')
#pval = stats.ttest_ind(responder_scores, non_responder_scores, alternative='two-sided')
print(pval)
