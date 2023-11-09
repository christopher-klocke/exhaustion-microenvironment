# bash script to run data download process for viral dataset
# run with 'data_download' conda environment active (see docstring of 'data_download_yost.py' script for further details)
python data_download_wang_A.py
tar -xvf "/Users/klockec/Documents/data/wang_data/download/wang.csv.tar" -C "/Users/klockec/Documents/data/wang_data/download"
ls "/Users/klockec/Documents/data/wang_data/download/" > "/Users/klockec/Documents/data/wang_data/download/wang_zipped_filenames.txt"
python data_download_wang_B.py
