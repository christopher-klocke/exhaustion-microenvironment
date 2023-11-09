# bash script to perform pre-processing for viral dataset
# run with 'sc' conda environment (described in 'sc_conda_env.ipynb')
python select_protein_coding_viral.py
python process1_wang.py
python process_wang_no6.py
