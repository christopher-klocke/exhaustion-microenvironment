"""
Classes and functions defined here are used by analysis scripts, e.g.: 
    li_analyze_oop.py

"""




## import statements

import scanpy as sc ########
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt 
import seaborn as sns #######@
from scipy.stats import spearmanr, kstest
from sklearn.cluster import k_means
from sklearn.decomposition import PCA





"""
take anndata object with pseudotime values as input and return: 


-construct numpy array

-perform k-means clustering on pseudotime values 


"""
def cluster_pt(adata, n_pt, n_sample, n_clusters, seed): 

    adata.obs[n_pt].to_numpy()

    adata.obs.index.to_numpy()

    a = adata.obs.index.to_numpy()
    b = adata.obs[n_pt].to_numpy()

    np.vstack((a,b)).T

    c = adata.obs[n_sample].to_numpy()

    np.vstack((a,b, c)).T

    ## yes, this works 
    ## yields array with 3 columns -- cell ID, pseudotime, and patient ID


    d = np.vstack((a,b, c)).T

    d[:,1]

    ## or just use this: 

    b

    ## keep it simple 



    k_means(X=b, n_clusters=n_clusters)


    ## ValueError: Expected 2D array, got 1D array instead:
    ## Reshape your data either using array.reshape(-1, 1) if your data has a single feature or array.reshape(1, -1) if it contains a single sample.

    b.reshape(-1, 1)

    k_means(X=b.reshape(-1, 1), n_clusters=n_clusters)

    ## that problem is fixed -- now fix 'Inf' problem 
    ## ValueError: Input contains NaN, infinity or a value too large for dtype('float64')

    ## sum([i == float('inf') for i in list(b)]) ## 31 'inf' values in li cd8 monocle results -- remove them 

    b[b != float('inf')] 

    len(b) ## 19997
    len(b[b != float('inf')]) ## 19966
    ## 31 'inf' values removed 

    ## need to use same boolean vector to index cell ID array -- remove those as well 

    a[b != float('inf')] 
    c[b != float('inf')] 

    ## this worked 
    ## try the clustering again 

    b_cut = b[b != float('inf')] 

    k_means(X=b_cut.reshape(-1, 1), n_clusters=n_clusters, random_state=seed)

    ## this worked 



    


def cluster_activities(adata, type, n_clusters, seed): 


    ## now create array for activity levels from other anndata object 
    ## pandas df is: 
    adata.obsm[type]

    ## as numpy array: 
    adata.obsm[type].to_numpy() ## need to add in cell IDs if needed too 



    aucell = adata.obsm[type].to_numpy()

    aucell[:,0] ## where 0 is first regulon -- iterate through 

    x=0

    k_means(X=aucell[:,x].reshape(-1, 1), n_clusters=n_clusters, random_state=seed)


## select 'k' for k-means clustering via histogram -- use 2 for now to set up script -- 



















## perform clustering
k_means(X=?, n_clusters=10, random_state=seed) ## optimize param 

## get cluster proportions
    
## run PCA
pca = PCA(n_components=15) ## optimize param
pca.fit(X) 
pca1 = pca.transform(X) ## ???

## plot by sample 

## Kolmogorov-Smirnov test
kstest() 

## calculate score 

## adata.obs['monocle3_pseudotime']



















## class definitions

"""
Instances of the CellType class can be initialized for each cell type in the 
dataset. These CellType object instances hold the cell type name ('type') and 
a list of clusters ('clusters') corresponding to the cell type in the dataset.

"""
class CellType(): 
    def __init__(self, type, clusters):
        self.type = type
        self.clusters = clusters

"""
Instances of the Sample class can be initialized for each sample in the 
dataset. These Sample object instances hold the sample id ('id'), the number 
of cells from this sample in the dataset ('count'), and the sample score 
('sample_score'), calculated during the sample-level ordering step and 
corresponding to the proportion of CD8 T cells (or other cell type of 
interest) with pseudotime values that fall above a chosen threshold -- 
calculated with FINISH DESCRIPTION 
"""
class Sample(): 
    def __init__(self, id, count): 
        self.id = id
        self.count = count
        self.sample_score = '' ## FIX? #################

"""
Instances of the GeneSetActivity class can be initialized for each cell type /
 sample combination in the dataset. These GeneSetActivity objects hold the 
 cell type ('celltype'), the sample ('sample'), the number of cells described 
by this combination, and a 'scores' dictionary. This dictionary has keys 
describing gene sets of interest (e.g. pyscenic regulons, Reactome pathways) 
and values (float) giving the average activity score (calculated with AUCell) 
of each gene set for this group of cells. 

The 'add_scores()' class method adds gene set ids and their average activity 
score values to the 'scores' attribute -- a dictionary. It takes as input an 
anndata object ('adata'), 'score_df' -- denoting the gene set pd.df name (e.g.
 'X_regulonsAUC')

By default, the 'add_scores()' method removes the last 3 characters from gene 
set name strings (e.g. '(+)' for pyscenic regulons). To retain full gene set 
name (e.g. for Reactome pathways), set 'full_name' to True. 
"""
class GeneSetActivity(): 
    def __init__(self, celltype, sample, count): 
        self.celltype = celltype
        self.sample = sample
        self.count = count
        self.scores = {}
    def add_scores(self, adata, score_df, full_name=False):
        df = adata.obsm[score_df]
        for tf in df.columns:
            q = df[[tf]]
            ## calculate average AUCell activity score for cells in this group
            score = sum(q[tf]) / len(q)
            if full_name: 
                self.scores[tf] = score
            else:
                self.scores[tf[:-3]] = score

## function definitions

"""


"""
def lists_create(adata, sample_name, gene_set, types, full_name=False): 
    types_list = []
    for elem in types: 
        types_list.append( CellType(type=elem[0], clusters=elem[1]) )
    samples_list = []
    activity_objects_list = []
    for x in set(list(adata.obs[sample_name])): ## FOR EACH PATIENT ID
        a = adata[adata.obs[sample_name] == x] ## anndata subsetted by patient
        samples_list.append( Sample(id=x, count=a.X.shape[0]) )
        for y in types_list: ## FOR EACH CELL TYPE
            b = y.clusters
            ## iterate through obs in a with list comprehension, for each, boolean retain if cluster ID in list
            c = a[[x in b for x in a.obs['leiden']]] ## subset anndata by cell type 
            if c.X.shape[0] != 0:
                sample_object = [n for n in samples_list if n.id == x][0] ## REWORK -- there has to be a better way to do this ############ ALSO may not even work -- check the reference / variable scope 
                g = GeneSetActivity(celltype=y, sample=sample_object, count=c.X.shape[0])
                g.add_scores(adata=c, score_df=gene_set, full_name=full_name)
                activity_objects_list.append(g)
    return types_list, samples_list, activity_objects_list

"""
THIS SHOULD BE TURNED INTO AN OBJECT METHOD


Create a dictionary that lists sample ids as keys and sample-level trajectory 
scores as values. Sample-level trajectory scores are just the proportion of 
cells from a sample that are above the threshold. 

input: sample_id_list, adata, sample_name, threshold

output: sample_id_list objects will be modified -- no data structure returned

"""
def add_sample_scores(sample_ids, samples, adata, sample_name, threshold): 
    sample_score_dict = {}
    for sample_id in sample_ids: 
        ## create patient-specific anndata objects -- subsetted by patient id 
        ## anndata object, subsetted by given sample; selecting dpt scores from adata.obs
        ## result is list of scores 
        sample_values = adata[adata.obs[sample_name] == sample_id].obs['dpt_pseudotime']
        sample_score = sum(sample_values > threshold) / len(sample_values) ## calculate proportion of cells in sample with dpt value above threshold -- this is sample score 
        sample_score_dict[sample_id] = sample_score
    ## first -- add 'sample_score_dict' sample scores to sample objects -- 
    for id in sample_score_dict.keys(): ## iterate through sample_score_dict keys (sample ids)
        for x in samples: 
            if id == x.id: ## find matching sample object
                x.sample_score = sample_score_dict[id] ## add score as attribute to sample object

"""
This function takes as input a list of CellType objects (each describing a 
cell type in the dataset), a list of Sample objects (each describing a sample 
in the dataset), and a list of GeneSetActivity objects (each describing
a population of cells specific to a particular cell type / sample 
combination). The Sample objects must already have sample_score values 
assigned. The function will calculate Spearman / rank-based correlation scores
 and corresponding p-values for the sample-level ordering values and the 
 various gene sets described in the GeneSetActivity.scores dictionaries. 

Arguments: 

celltypes

samples

activities

Output: 

A dictionary of results -- keys are cell types, values are dictionaries in 
which keys are gene set ids and values are tuples containing rho (correlation 
score) and corresponding p-value

"""
def run_score_correlation(celltypes, samples, activities): 
    final_scores_dict = {}
    for celltype in celltypes: ## iterate through cell type objects
        current_type = celltype.type
        correlate_tuples_list = []
        final_scores_dict[current_type] = {}
        for r in activities: ## search activity object list
            if r.celltype.type == current_type: ## for those matching current cell type
                for s in samples: 
                    if r.sample.id == s.id: 
                        if s.sample_score != '': ## not the best way to do this -- non-matching samples (e.g. removed PBMC controls) should not be getting through -- fix upstream code ########## dev 
                            correlate_tuples_list.append((s.sample_score, r.scores)) ## tuple is sample score, dictionary -- ensure strings to floats, etc. ############ dev -- LOOKS FINE FOR NOW
        if correlate_tuples_list != []: ## FIX BUG -- ERROR IF AN EMPTY LIST RETURNED -- BYPASS THIS ##########
            for j in correlate_tuples_list[0][1].keys(): ## iterate through gene sets / regulons
                sample_scores = [] ## update code to make more efficient -- create this list once instead of every time for each gene set ############################# dev 
                gene_set_scores = []
                for k in correlate_tuples_list: ## iterate through tuples list
                    sample_scores.append(k[0])
                    gene_set_scores.append(k[1][j])
                corr_dict = {}
                corr_dict['sample_scores'] = sample_scores
                corr_dict['gene_set_scores'] = gene_set_scores
                df = pd.DataFrame(corr_dict) ## create DataFrame
                rho, p = spearmanr(df['sample_scores'], df['gene_set_scores']) ## calculate Spearman Rank correlation and corresponding p-value
                final_scores_dict[current_type][j] = (rho, p) ## add to results dictionary
    return final_scores_dict

def get_significant_results(scores_dict, alpha=0.05):
    for celltype in scores_dict.keys(): 
        g_dict = scores_dict[celltype]
        for r in g_dict.keys(): 
            p = g_dict[r][1]
            if p < alpha: ## is p value below chosen alpha? 
                print('cell type: ' + celltype + ' gene set: ' + r + ' p-value: ' + str(p))

"""
get results for one TF / pathway / gene set for all cell type clusters

"""
def get_gene_set_results(scores_dict, gene_set):
    for celltype in scores_dict.keys(): 
        g_dict = scores_dict[celltype]
        try: 
            res = g_dict[gene_set]
            print('cell type: ' + celltype + ' gene set: ' + gene_set + ' p-value: ' + str(res[1]) + ' rho: ' + str(res[0]))
        except KeyError: 
            print('key ' + gene_set + ' not in ' + celltype + ' object')

"""
plotting function for significant results 

"""
def plot_scores(results, celltype, gene_set, samples, activities, save=False, savefile=''): 
    res = results[celltype][gene_set]
    print('cell type: ' + celltype + '\n' + 'gene set: ' + gene_set + '\n' + 'rho: ' + str(res[0]) + '\n' + 'p-value: ' + str(res[1]))

    tuples_list = []
    for r in activities: 
        if r.celltype.type == celltype: 
            for s in samples: 
                if r.sample.id == s.id: 
                    if s.sample_score != '': ## FIX -- SEE similar section in 'run_score_correlation' 
                        tuples_list.append((s.sample_score, r.scores[gene_set]))
    sample_scores = []
    gene_set_scores = []
    for k in tuples_list: 
        sample_scores.append(k[0])
        gene_set_scores.append(k[1])

    ########### 

    if save: 

        x = np.array(sample_scores)
        y = np.array(gene_set_scores)
        plt.scatter(x, y)
        plt.savefig(savefile)
        plt.show(block=True)

    else: 
        x = np.array(sample_scores)
        y = np.array(gene_set_scores)
        plt.scatter(x, y)
        plt.show(block=True)












"""
plotting function for significant results 

updated version

"""
def plot_scores_v2(results, celltype, gene_set, samples, activities, save=False, savefile=''): 
    res = results[celltype][gene_set]
    print('cell type: ' + celltype + '\n' + 'gene set: ' + gene_set + '\n' + 'rho: ' + str(res[0]) + '\n' + 'p-value: ' + str(res[1]))

    tuples_list = []
    for r in activities: 
        if r.celltype.type == celltype: 
            for s in samples: 
                if r.sample.id == s.id: 
                    if s.sample_score != '': ## FIX -- SEE similar section in 'run_score_correlation' 
                        tuples_list.append((s.sample_score, r.scores[gene_set]))
    sample_scores = []
    gene_set_scores = []
    for k in tuples_list: 
        sample_scores.append(k[0])
        gene_set_scores.append(k[1])

    ########### 

    if save: 

        plt.clf()

        x = np.array(sample_scores)
        y = np.array(gene_set_scores)
        plt.scatter(x, y) ## need to add s -- counts -- point size 

        plt.xlabel('sample-level exhaustion score')
        plt.ylabel('gene set activity')

        # giving title to the plot
        title_text = 'cell type: ' + celltype + ' gene set: ' + gene_set + ' rho: ' + str(round(res[0], 2)) + ' p-value: ' + str(round(res[1], 5))
        plt.title(title_text)


        """
        # Fit linear regression via least squares with numpy.polyfit
        # It returns an slope (b) and intercept (a)
        # deg=1 means linear fit (i.e. polynomial of degree 1)
        b, a = np.polyfit(x, y, deg=1)

        # Create sequence of 100 numbers from 0 to 100 
        xseq = np.linspace(0, 10, num=100)

        # Plot regression line
        #plt.plot(xseq, a + b * xseq, color="k", lw=2.5)
        """



        """
        # use regplot
        sns.regplot(x = "sepal_length",
                    y = "petal_length", 
                    ci = None,
                    data = df)


        pd.DataFrame()

        """



        plt.savefig(savefile)
        #plt.show(block=True) ################

    else: 

        plt.clf()

        x = np.array(sample_scores)
        y = np.array(gene_set_scores)
        plt.scatter(x, y) ## need to add s -- counts -- point size 

        plt.xlabel('sample-level exhaustion score')
        plt.ylabel('gene set activity')

        # giving title to the plot
        title_text = 'cell type: ' + celltype + ' gene set: ' + gene_set + ' rho: ' + str(round(res[0], 2)) + ' p-value: ' + str(round(res[1], 5))
        plt.title(title_text)

        plt.show(block=True)