import numpy as np
import pandas as pd
import os
import requests
from zipfile import ZipFile
from io import BytesIO
from multiprocessing import Pool
from functools import partial

dfMetadata = pd.read_csv('data/expression/AMBA_metadata.csv')

dfTemp = dfMetadata.loc[:50].copy()


def download_zip(experimentid, outdir):
    
    url = 'http://api.brain-map.org/grid_data/download/{}'.format(experimentid)
    print('Experiment {} pid {}'.format(experimentid, os.getpid()))
    response = requests.get(url)
    sourceZip = ZipFile(BytesIO(response.content))
    sourceZip.extract('energy.raw', path=outdir)
    sourceZip.close()
    print('Experiment {} pid {}'.format(experimentid, os.getpid()))
    os.rename(outdir+'energy.raw', outdir+str(experimentid)+'.raw')

    
experimentids = [dfTemp.loc[i, 'experiment_id'] for i in range(0, dfTemp.shape[0])]

pool = Pool(4)
download_func = partial(download_zip, outdir = 'data/expression/coronal/')
results = pool.map(download_func, experimentids)
pool.close()
pool.join()