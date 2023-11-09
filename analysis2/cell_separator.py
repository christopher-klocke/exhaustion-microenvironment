"""

Separate cell types that are mixed in a single cluster

"""

## import statements
import matplotlib.pyplot as plt

## function definitions

"""
This function takes as input an anndata object, a list of cluster ids, and 
a dictionary with keys 'higher' and 'lower' corresponding to lists of genes 
expected to be more highly or lowly expressed in one cell type of interest vs.
 the other. 


"""
def separate_score(adata, clusters, genes_dict, save=False, savefile=''): 
    cells = adata[[x in clusters for x in adata.obs['leiden']]]

    up = genes_dict['higher']
    down = genes_dict['lower']

    scores = []

    total = cells.shape[0] ##############
    counter = 0 ###############

    for cell in cells.obs.index: ## iterate through cells in subsetted anndata
        score = 0
        for gene in up: ## iterate through up genes
            score += float(cells[cell, gene].X.todense())
        
        for gene in down: 
            score -= float(cells[cell, gene].X.todense())

        scores.append((cell, score))

        counter += 1 #################
        print(str(total - counter)) ###############

    if save: 
        with open(savefile, 'w') as outfile: 
            for elem in scores: 
                outfile.write(elem[0] + ',' + str(elem[1]) + '\n')

    return scores
            

## plot histogram
def separate_plot_histograms(scores_list=False, file=False, threshold=False): 

    scores = []

    if scores_list != False: 
        scores=scores_list
 
    if file != False: 
        with open(file, 'r') as infile: 
            for line in infile: 
                line=line.rstrip()
                line = line.split(sep=',')
                scores.append((line[0], line[1]))

    just_scores = [float(x[1]) for x in scores]

    fig, axs = plt.subplots(1, 1, sharey=True, tight_layout=True)
    #plt.axvline(x=threshold, color='red') ## plot threshold line 
    axs.hist(just_scores)
    plt.show(block=True)

## separate cells -- two lists of cell ids ## NOTE: two thresholds, ignore ambiguous middle? 
def separate_id_lists(): 
    pass 

