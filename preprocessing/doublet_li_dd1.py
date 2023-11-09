"""

run DoubletDetection method on Li dataset 

https://github.com/JonathanShor/DoubletDetection#running-doubletdetection

run with 'doubletdetection' conda env (instructions in 'roadmap')




raw_counts is a scRNA-seq count matrix (cells by genes), and is array-like

labels is a 1-dimensional numpy ndarray with the value 1 representing a detected doublet, 0 a singlet, and np.nan an ambiguous cell

scores is a 1-dimensional numpy ndarray representing a score for how likely a cell is to be a doublet. The score is used to create the labels


need to correct by sample 

"""

## import statements 

import scanpy as sc
import doubletdetection as dd

## setup 

li_dir = '/Users/klockec/Documents/data/li_data/'

li_adata_filename = 'li_viz_ready.h5ad'


li_dd1_filename = 'li_dd1_save.csv'

## load the data

adata = sc.read_h5ad(li_dir + li_adata_filename)




## main function definition
def main(): 

    ## process the data

    doublet_dict = {}

    for p in adata.obs['patient']: 
        adata_x = adata[adata.obs['patient'] == p]
        raw_counts = adata_x.X # raw_counts is a cells by genes count matrix

        cell_ids = adata_x.obs.index

        clf = dd.BoostClassifier(clustering_algorithm='leiden')
        labels = clf.fit(raw_counts).predict()
        scores = clf.doublet_score() # higher means more likely to be doublet

        doublet_dict[p] = (cell_ids, labels, scores)


    """
    save the scores and labels for each sample -- can use a dictionary

    keys are sample ids, values are 2-tuples of arrays -- labels and scores 

    then can look at individually or aggregate and visualize scores on UMAP plot -- pull UMAP from processed data -- load in another anndata file (from later in analysis) 



    need to add array of cell IDs as well -- make it a 3-tuple 

    """

    print(doublet_dict) #################

    ## PLOT RESULTS 

    ### plot UMAP of (processed) adata annotated with 



    """
    
    need to be able to save and re-load -- 

    save as pandas df as attribute of anndata object? or just as independent df? 

    dataframe -- score per 
    
    
    """

    with open((li_dir + li_dd1_filename), 'w') as output1: 
        for p in doublet_dict: 
            cell_ids = doublet_dict[p][0]
            labels = doublet_dict[p][1]
            scores = doublet_dict[p][2]

            counter = len(cell_ids)
            for i in range(counter): 
                output1.write(cell_ids[i] + ',' + p + ',' + str(labels[i]) + ',' str(scores[i]) + '\n')


## run main function 
if __name__ == "__main__": 
    main()
