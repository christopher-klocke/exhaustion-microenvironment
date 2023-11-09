"""-Use this script to download the data from Gene Expression Omnibus (GEO).
-Use the following to create the conda environment "data_download" to run this script:
    conda create --name data_download
    conda activate data_download
    conda install -c conda-forge urllib3
"""
# import statements
from download_unzip import download_and_unzip


# download the data
dir_name = '/Users/klockec/Documents/data/li_data/download/' ## file path for rws07890
links_and_files = [('https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE123139&format=file', 'GSE123139_RAW.tar')]
for elem in links_and_files: 
    download_url = elem[0]
    unzipped_filename = dir_name + elem[1]
    download_and_unzip(url=download_url, unzipped=unzipped_filename, file_zipped=False)
    