"""

This script will be used to compare gene sets (regulons, pathways) to look for shared genes and calculate percentage overlap. 

First application -- look at Li POU6F1 and Wang POU5F1 regulon gene sets -- expression of each in macrophages correlated 
with CD8 T cell sample-level exhaustion trajectories. 

"""

## import statements


## file paths 

li_dir = '/Users/klockec/Documents/data/analysis_files/p3/pyscenic_runs/p3/li_run/'

wang_dir = '/Users/klockec/Documents/data/analysis_files/p3/pyscenic_runs/p3/wang_run/'

gmt_filename = 'ctx_output1.gmt'

## load the data

li_dict = {}

wang_dict = {}

with open((li_dir + gmt_filename), 'r') as li_gmt: 
    for line in li_gmt: 
        line = line.rstrip() 
        line = line.split(sep='\t')
        tf = line[0][:-3]
        genes = line[2:]
        li_dict[tf] = genes

with open((wang_dir + gmt_filename), 'r') as wang_gmt: 
    for line in wang_gmt: 
        line = line.rstrip() 
        line = line.split(sep='\t')
        tf = line[0][:-3]
        genes = line[2:]
        wang_dict[tf] = genes

## process the data 

"""

have two dictionaries -- now, pull gene lists of interest and compare 

"""

