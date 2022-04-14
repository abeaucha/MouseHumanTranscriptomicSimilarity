# ----------------------------------------------------------------------------
# download_ahba.py 
# Author: Antoine Beauchamp
# Created: February 10th, 2022

"""
Download data from the Allen Human Brain Atlas

Description
-----------
This script downloads the microarray data sets from all six donors in the
AHBA. It also downloads the hierarchical ontology from the AHBA.
"""

# Packages -------------------------------------------------------------------

import os
import argparse
import abagen
import requests
import json

# Command line arguments -----------------------------------------------------

def parse_args():
   
    """Parse command line arguments"""

    parser = argparse.ArgumentParser(
                 formatter_class = argparse.ArgumentDefaultsHelpFormatter
             )
    
    parser.add_argument(
        '--outdir',
        type = str,
        default = 'data/',
        help = "Directory in which to download the data"
    )
    
    parser.add_argument(
        '--verbose',
        type = str,
        default = 'true',
        choices = ['true', 'false'],
        help = "Verbosity"
    )
    
    args = vars(parser.parse_args())
    
    return args

# Main -----------------------------------------------------------------------

def main():

    #Get command line arguments
    args = parse_args()
    outdir = args['outdir']
    verbose = True if args['verbose'] == 'true' else False
    
    #Format the path properly
    outdir = os.path.join(outdir, '')
    
    if not os.path.exists(outdir):
        os.mkdir(outdir)
        
    #Download AHBA hierarchical ontology as JSON
    if verbose:
        print("Downloading AHBA hierarchical ontology...")
    hierarchy_url = ('http://api.brain-map.org/api/v2/'
                     'structure_graph_download/10.json')
    hierarchy_url_get = requests.get(hierarchy_url)
    hierarchy_dict = json.loads(hierarchy_url_get.text)
    with open(outdir+'AHBA_hierarchy_definitions.json', 'w') as outfile: 
        json.dump(hierarchy_dict, outfile)
    
    #Include 'microarray' sub-directory. abagen will download to ~/ otherwise
    outdir = os.path.join(outdir, 'microarray', '')
    if not os.path.exists(outdir):
        os.mkdir(outdir)

    #Download microarray data
    if verbose:
        print("Downloading AHBA microarray data...")
    files = abagen.fetch_microarray(donors = 'all', data_dir = outdir)

    return
    
if __name__ == '__main__':
    main()
