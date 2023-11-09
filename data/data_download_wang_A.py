""" from 'data_download_viral.ipynb' """
# import statements
from download_unzip import download_and_unzip


# file paths
wang_dir = '/Users/klockec/Documents/data/wang_data/download/'

# download the data
download_and_unzip(url='https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE157829&format=file', unzipped=(wang_dir + 'wang.csv.tar'), file_zipped=False)
