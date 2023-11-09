"""-Use this script to download the data from Gene Expression Omnibus (GEO).
-Use the following to create the conda environment "data_download" to run this script:
    conda create --name data_download
    conda activate data_download
    conda install -c conda-forge urllib3
    conda install -c anaconda numpy
    conda install -c anaconda pandas
    conda install -c conda-forge scanpy python-igraph leidenalg
"""
# import statements
from download_unzip import download_and_unzip


# download the data
dir_name = '/Users/klockec/Documents/data/yost_data/download/' # file path for rws07890
links_and_files = [('https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE123813&format=file&file=GSE123813%5Fbcc%5Fall%5Fmetadata%2Etxt%2Egz', 'GSE123813_bcc_all_metadata.txt.gz', 'GSE123813_bcc_all_metadata.txt'), 
    ('https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE123813&format=file&file=GSE123813%5Fbcc%5FscRNA%5Fcounts%2Etxt%2Egz', 'GSE123813_bcc_scRNA_counts.txt.gz', 'GSE123813_bcc_scRNA_counts.txt'), 
    ('https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE123813&format=file&file=GSE123813%5Fbcc%5Ftcell%5Fmetadata%2Etxt%2Egz', 'GSE123813_bcc_tcell_metadata.txt.gz', 'GSE123813_bcc_tcell_metadata.txt'), 
    ('https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE123813&format=file&file=GSE123813%5Fbcc%5Ftcr%2Etxt%2Egz', 'GSE123813_bcc_tcr.txt.gz', 'GSE123813_bcc_tcr.txt'), 
    ('https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE123813&format=file&file=GSE123813%5Fscc%5Fmetadata%2Etxt%2Egz', 'GSE123813_scc_metadata.txt.gz', 'GSE123813_scc_metadata.txt'), 
    ('https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE123813&format=file&file=GSE123813%5Fscc%5FscRNA%5Fcounts%2Etxt%2Egz', 'GSE123813_scc_scRNA_counts.txt.gz', 'GSE123813_scc_scRNA_counts.txt'), 
    ('https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE123813&format=file&file=GSE123813%5Fscc%5Ftcr%2Etxt%2Egz', 'GSE123813_scc_tcr.txt.gz', 'GSE123813_scc_tcr.txt')
]

for elem in links_and_files: 
    download_url = elem[0]
    zipped_filename = dir_name + elem[1]
    unzipped_filename = dir_name + elem[2]
    download_and_unzip(url=download_url, zipped=zipped_filename, unzipped=unzipped_filename)
