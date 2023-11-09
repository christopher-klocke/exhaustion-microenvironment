"""

functions for downloading and unzipping data files -- used by 'data_download_*.py' scripts 

"""

## import statements 

import urllib3
import gzip
import shutil

## function definitions

def download_and_unzip(url, unzipped, zipped='', file_zipped=True): 
    http = urllib3.PoolManager()
    get_object = http.request('GET', url, preload_content=False)

    if file_zipped:
        filename = zipped
    else:
        filename = unzipped

    with open(filename, 'wb') as out: 
        while True: 
            data = get_object.read()
            if not data: 
                break
            out.write(data)

    get_object.release_conn()

    if file_zipped:

        with gzip.open(zipped, 'rb') as f_in:
            with open(unzipped, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)

def unzip(zipped): ## function argument is path to the zipped file
    unzipped=zipped[:-3] ## this assumes zipped file suffix is '.gz', modify if different
    with gzip.open(zipped, 'rb') as f_in: 
        with open(unzipped, 'wb') as f_out: 
            shutil.copyfileobj(f_in, f_out)

