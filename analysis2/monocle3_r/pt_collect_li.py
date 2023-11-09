"""
load in results from monocle3 (run in R)

save-as of 'run2_collect_li.py'

from 'li_pt_eval4.R'
"""

## import statements
import scanpy as sc

## setup
input_filename = '/Users/klockec/OneDrive - Oregon Health & Science University/data_backup_rws07890/data_january2023_backup/data/li_data/li_T_refined.h5ad'
pt_filename = "/Users/klockec/OneDrive - Oregon Health & Science University/data/analysis_files/p3/mod1/pt_optimize/li_runA/pseudotime.csv"
output_filename = '/Users/klockec/OneDrive - Oregon Health & Science University/data_backup_rws07890/data_january2023_backup/data/li_data/li_t_monocle_run4.h5ad'

## function definitions
def add_pt_values(adata, pt_filename, remove_last: bool = False):
    pt_dict = {}
    with open(pt_filename, 'r') as input1: 
        for line in input1: 
            line = line.rstrip()
            if line[0:2] != '""': 
                line = line.split(',')
                cellID = line[0][1:-1]
                pt_value = float(line[1])
                pt_dict[cellID] = pt_value
    pseudotime_list = []
    for i in list(adata.obs.index):
        if remove_last:
            cell_id = '-'.join(i.split('-')[:-1])  # https://stackoverflow.com/questions/44778/how-would-you-make-a-comma-separated-string-from-a-list-of-strings
            pseudotime_list.append(pt_dict[cell_id])
        else:
            pseudotime_list.append(pt_dict[i])
    adata.obs['monocle3_pseudotime'] = pseudotime_list

    return adata
    
## main function definition
def main():
    adata = sc.read_h5ad(input_filename)
    #adata = add_pt_values(adata=adata, pt_filename=pt_filename)
    adata = add_pt_values(adata=adata, pt_filename=pt_filename, remove_last=True)
    adata.write_h5ad(output_filename)

## run main function
if __name__ == "__main__": 
    main()
