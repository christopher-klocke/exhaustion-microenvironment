""" run with 'sc2' conda env 

caffeinate -i python3 csea_oop_li_updated.py >csea_oop_li4_output.log 2>&1 &
#nohup python3 csea_oop_li_updated.py >csea_oop_li4_output.log 2>&1 &

save-as of 'csea_oop_li1.py'
"""

import logging 

logging.basicConfig(filename='csea_oop_li4.log', filemode='w', level=logging.INFO,
                    format='%(asctime)s - %(levelname)s - %(message)s')
logging.info('run start')

import cell_enrichment_utils_oop as ceu_oop

#li_path = '/Users/klockec/Documents/data/li_data/' ## file path for rws07890
li_pt_input_filename = '/Users/klockec/OneDrive - Oregon Health & Science University/data_backup_rws07890/data_january2023_backup/data/li_data/li_t_monocle_run4.h5ad'  # from 'p3/analysis2/monocle3_r/run4_collect_li.py'
li_act_input_filename = '/Users/klockec/OneDrive - Oregon Health & Science University/data_backup_rws07890/data_january2023_backup/data/li_data/li_dd_viz_ready_regs.h5ad'
results_dir = '/Users/klockec/OneDrive - Oregon Health & Science University/data/analysis_files/p3/mod1/li_run4_csea_oop/'

#random_state = 1
cell_types_input = {'macrophage': ['3', '11'], 'NK': ['1'], 'B cells': ['5'], 'plasma': ['8'], 'ILC???': ['12'], 'CD8 T cells': ['0', '2', '4', '6', '9', '10'], '???': ['7']}

csea_args = {'adata_pseudotime': li_pt_input_filename,
            'adata_activity': li_act_input_filename,
            'celltypes_dict': cell_types_input, 
            'output_dir': results_dir, 
            'sample_id': 'patient',
            'set_num_permutations': 1000, 
            'fdr': 0.05
            }

li_csea = ceu_oop.RunCSEA(**csea_args)
#li_csea.csea()
li_csea.csea(just_pt=True)
#li_csea.correlation()
#li_csea.permute()
#li_csea.p_values()
#li_csea.export_permutation_runs()
li_csea.export_sample_pt_scores()
