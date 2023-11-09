"""continuation of Li dataset download steps """
# import statements
from download_unzip import unzip

# process the data
dir_name = '/Users/klockec/Documents/data/li_data/download/'
zipped_filenames_path = '/Users/klockec/Documents/data/li_data/download/zipped_sample_files.txt'
zipped_filenames_list = []

with open(zipped_filenames_path, 'r') as input1: 
    for line in input1: 
        line = line.rstrip()
        if line[0:3] == 'GSM': 
            zipped_filenames_list.append(line)

for elem in zipped_filenames_list: 
    unzip(zipped=(dir_name + elem))
