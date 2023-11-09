""" run with 'sc2' conda env 

caffeinate -i python3 csea_oop_yost_updated.py >csea_oop_yost4_output.log 2>&1 &
#nohup python3 csea_oop_yost_updated.py >csea_oop_yost4_output.log 2>&1 &

save-as of 'csea_oop_li_updated.py'
"""

import logging 

logging.basicConfig(filename='csea_oop_yost4.log', filemode='w', level=logging.INFO,
                    format='%(asctime)s - %(levelname)s - %(message)s')
logging.info('run start')

import cell_enrichment_utils_oop as ceu_oop

pt_input_filename = '/Users/klockec/OneDrive - Oregon Health & Science University/data_backup_rws07890/data_january2023_backup/data/yost_data/yost_t_monocle_run4.h5ad'  # from 'p3/analysis2/monocle3_r/run4_collect_yost.py'
act_input_filename = ''  # not needed
results_dir = '/Users/klockec/OneDrive - Oregon Health & Science University/data/analysis_files/p3/mod1/yost_run4_csea_oop/'

#random_state = 1
cell_types_input = {'CD8_T': ['0', '2', '3', '5', '6']}  # don't need to update -- will not be used here

csea_args = {'adata_pseudotime': pt_input_filename,
            'adata_activity': act_input_filename,
            'celltypes_dict': cell_types_input, 
            'output_dir': results_dir, 
            'sample_id': 'sample',
            'set_num_permutations': 1000, 
            'fdr': 0.05
            }

yost_csea = ceu_oop.RunCSEA(**csea_args)
yost_csea.csea(just_pt=True)
yost_csea.export_sample_pt_scores()
