# bash script to run data download process
# run with 'data_download' conda environment active (see docstring of 'data_download_yost.py' script for further details)
python data_download_yost.py
python data_download_li_A.py
tar -xvkf "/Users/klockec/Documents/data/li_data/download/GSE123139_RAW.tar" -C "/Users/klockec/Documents/data/li_data/download"
ls "/Users/klockec/Documents/data/li_data/download/" > "/Users/klockec/Documents/data/li_data/download/zipped_sample_files.txt"
python data_download_li_B.py
ls "/Users/klockec/Documents/data/li_data/download/" > "/Users/klockec/Documents/data/li_data/download/sample_files_unzipped.txt"
python data_download_li_C.py
python data_download_li_add_metadata.py
