"""

look at T cell highly variable gene sets

"""

## import statements 

import scanpy as sc

## setup 

random1 = 1

## file paths 

li_path = '/Users/klockec/Documents/data/li_data/' ## file path for rws07890

li_just_T_filename = 'li_T_viz_ready.h5ad'

wang_dir = '/Users/klockec/Documents/data/wang_data/'

wang_just_T_filename = 'wang_T_viz_ready.h5ad'

## load the data

adata_li = sc.read_h5ad(li_path + li_just_T_filename)

adata_wang = sc.read_h5ad(wang_dir + wang_just_T_filename)

## process the data

li_list = list(adata_li[:, adata_li.var['highly_variable']].var.index)
wang_list = list(adata_wang[:, adata_wang.var['highly_variable']].var.index)

print('Li T cell highly variable genes: ')
print(str(len(li_list)))
print(li_list)

print('Wang T cell highly variable genes: ')
print(str(len(wang_list)))
print(wang_list)

shared = []
just_li = []
just_wang = []

full_list = []

for x in li_list: 
    full_list.append(x)

for x in wang_list: 
    full_list.append(x)

full_set = list(set(full_list))

for x in full_set: 
    if x in li_list: 
        if x in wang_list: 
            shared.append(x) 
        else: 
            just_li.append(x)
    else:
        just_wang.append(x)

print('shared genes: ')
print(shared)
print('just Li: ')
print(just_li)
print('just wang: ')
print(just_wang)

"""


"""

