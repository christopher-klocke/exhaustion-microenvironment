"""-run with 'sc2' conda env
#nohup python3 csea_oop_wang_updated.py >csea_oop_wang4_output.log 2>&1 &
"""
# set up logging
import logging
logging.basicConfig(filename='csea_oop_wang4.log', filemode='w', level=logging.INFO,
                    format='%(asctime)s - %(levelname)s - %(message)s')
logging.info('run start')

# import statements
import cell_enrichment_utils_oop as ceu_oop


# setup
pt_input_filename = '/Users/klockec/OneDrive - Oregon Health & Science University/data_backup_rws07890/data_january2023_backup/data/wang_data/wang_t_monocle_run4.h5ad'  # from 'p3/analysis2/monocle3_r/run4_collect_wang.py'
act_input_filename = ''  # not needed
results_dir = '/Users/klockec/OneDrive - Oregon Health & Science University/data/analysis_files/p3/mod1/wang_run4_csea_oop/'

#random_state = 1
cell_types_input = {'CD8_T': ['0', '2', '3', '5', '6']}  # don't need to update -- will not be used here

csea_args = {'adata_pseudotime': pt_input_filename,
            'adata_activity': act_input_filename,
            'celltypes_dict': cell_types_input, 
            'output_dir': results_dir, 
            'sample_id': 'batch',
            'set_num_permutations': 1000, 
            'fdr': 0.05
            }

wang_csea = ceu_oop.RunCSEA(**csea_args)
wang_csea.csea(just_pt=True)
wang_csea.export_sample_pt_scores()
