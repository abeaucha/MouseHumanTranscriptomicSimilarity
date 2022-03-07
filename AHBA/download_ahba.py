# download_ahba.py 
#
# Download AHBA ontology and microarray data from the web API
# 
# Antoine Beauchamp
# Created: February 10th, 2022
# Edited: March 5th, 2022

# Packages -------------------------------------------------------

import os
import argparse
import abagen
import requests
import json

# Functions -------------------------------------------------------

def parse_args():
    
    parser = argparse.ArgumentParser()
    
    parser.add_argument(
        '--outdir',
        type = str,
        default = 'data/',
        help = 'Directory in which to download the data'
    )
    
    args = vars(parser.parse_args())
    
    return args



def main():

    #Get command line arguments
    args = parse_args()
    outdir = args['outdir']
    
    #Format the path properly
    outdir = os.path.join(outdir, '')
    
    if not os.path.exists(outdir):
        os.mkdir(outdir)
        
    #Download AHBA hierarchical ontology as JSON
    print('Downloading AHBA hierarchical ontology...')
    hierarchy_url = "http://api.brain-map.org/api/v2/structure_graph_download/10.json"
    hierarchy_url_get = requests.get(hierarchy_url)
    hierarchy_dict = json.loads(hierarchy_url_get.text)
    with open(outdir+"AHBA_hierarchy_definitions.json", "w") as outfile: 
        json.dump(hierarchy_dict, outfile)
    
    #Include 'microarray' sub-directory. abagen will download to ~/ otherwise
    outdir = os.path.join(outdir, 'microarray', '')
    if not os.path.exists(outdir):
        os.mkdir(outdir)

    #Download microarray data
    print('Downloading AHBA microarray data...')
    files = abagen.fetch_microarray(donors = 'all', data_dir = outdir)

    return
    
    
if __name__ == '__main__':
    main()