# bash script to perform pre-processing
# run with 'sc' conda environment (described in 'sc_conda_env.ipynb')
python process1_yost.py
python process1_li.py
python select_protein_coding.py
python process2_yost.py
python process2_li.py
python metadata_filtering_yost.py
