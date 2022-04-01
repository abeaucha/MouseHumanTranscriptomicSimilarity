# -----------------------------------------------------------------------------
# download_AMBA.py
# Authors: Antoine Beauchamp, Yohan Yee

"""
Download in-situ hybridization data sets from the Allen Mouse Brain Atlas
"""

# Packages --------------------------------------------------------------------

import os
import argparse
import requests
import numpy                as np
import pandas               as pd
import multiprocessing      as mp
from pyminc.volumes.factory import volumeFromFile, volumeFromDescription
from zipfile                import ZipFile
from functools              import partial
from tqdm                   import tqdm

# Functions -------------------------------------------------------------------

def parse_args():

    """Parse command line arguments"""

    parser = argparse.ArgumentParser(
                 formatter_class = argparse.ArgumentDefaultsHelpFormatter
             )

    parser.add_argument(
        '--dataset',
        type = str,
        default = 'coronal',
        choices = ['coronal', 'sagittal'],
        help = 'AMBA dataset to import.'
    )

    parser.add_argument(
        '--outdir',
        type = str,
        default = 'expression/',
        help = 'Path to directory in which to download the data.'
    )
    
    parser.add_argument(
        '--metadata',
        type = str,
        default = 'metadata.csv',
        help = 'File in --outdir containing AMBA metadata.'
    )
    
    parser.add_argument(
        '--parallel',
        type = str,
        default = 'false',
        choices = ['true', 'false'],
        help = 'Option to run in parallel.'
    )
    
    parser.add_argument(
        '--nproc',
        type = int,
        default = mp.cpu_count(),
        help = 'Number of CPUs to use in parallel.'
    )
    
    parser.add_argument(
        '--verbose',
        type = str,
        default = 'true',
        choices = ['true', 'false'],
        help = 'Verbosity.'
    )

    args = vars(parser.parse_args())

    return args


def fetch_metadata(dataset = 'coronal', outdir='./', 
                   outfile = 'metadata.csv', verbose = True):

    """
    Download metadata for in-situ hybridization data sets
    
    Arguments
    ---------
    dataset: str, optional
        String indicating whether to download the coronal or sagittal
        (default 'coronal')
    outdir: str, optional
        Directory in which to download the metadata.
        (default './')
    outfile: str, optional
        Name of csv file in which to save the metadata. 
        (default 'metadata.csv')
        
    Returns
    -------
    None 
    """        

    abi_query_metadata = ("http://api.brain-map.org/api/v2/data/SectionDataSet/"
                          "query.csv?criteria="
                          "[failed$eqfalse],"
                          "plane_of_section[name$eq{}],"
                          "products[abbreviation$eqMouse],"
                          "treatments[name$eqISH],"
                          "genes"
                          "&tabular="
                          "data_sets.id+as+experiment_id,"
                          "data_sets.section_thickness,"
                          "data_sets.specimen_id,"
                          "plane_of_sections.name+as+plane,"
                          "genes.acronym+as+gene,"
                          "genes.name+as+gene_name,"
                          "genes.chromosome_id,"
                          "genes.entrez_id,"
                          "genes.genomic_reference_update_id,"
                          "genes.homologene_id,"
                          "genes.organism_id"
                          "&start_row=0"
                          "&num_rows=all".format(dataset))

    (pd.read_csv(abi_query_metadata)
        .drop_duplicates()
        .to_csv(outdir+outfile, index=False))

    if verbose:
        print('Metadata downloaded at: {}'.format(outdir+outfile))
    
    return 


def fetch_expression(experiment_id, outdir = './tmp/'):

    """
    Download the expression energy for an ISH experiment

    Description
    -----------
    This function downloads the gridded expression energy for a single
    in-situ hybridization experiment in the Allen Mouse Brain Atlas.
    The downloaded data is in RAW format and is stored in a file named
    according to the experiment ID. 

    Arguments
    ---------
    experiment_id: int
        The ID of the in-situ hybridization experiment to download.
    outdir: str, optional
        Directory in which to download the expression energy file.
        (default './tmp/')   

    Returns
    -------
    outfile: str
        The name of the file containing the expression energy.
    success: int
        Integer indicating whether the download was successful.
    """
    
    if not os.path.exists(outdir):
        os.mkdir(outdir)

    tmpdir = outdir+str(experiment_id)+'/' 
    os.mkdir(tmpdir)
        
    abi_query_expr = ('http://api.brain-map.org/grid_data/download/{}'
                      .format(experiment_id))
    amba_request = requests.get(abi_query_expr)

    tmpfile = tmpdir+str(experiment_id)+'.zip'
    with open(tmpfile, 'wb') as file:
        file.write(amba_request.content)

    outfile = outdir+str(experiment_id)+'.raw'
    with ZipFile(tmpfile, 'r') as file:
        try:
            file.extract('energy.raw', path = tmpdir)
            os.rename(tmpdir+'energy.raw', outfile)
            success = 1
        except KeyError as err:
            print('Error for experiment {}: {}. Ignoring.'
                  .format(experiment_id, err))
            success = 0
            
    os.remove(tmpfile)
    os.rmdir(tmpdir)
 
    return outfile, success


def rawtominc_wrapper(infile, outfile = None, keep_raw = False):
    
    """
    Convert a RAW file to MINC format

    Arguments
    ---------
    infile: str
    outfile: str, optional
        Name of output MINC file. If None, will use the same name as
        the RAW file.
        (default None)
    keep_raw: bool, optional
        Option to keep RAW file. If False, the input file will be
        deleted.
        (default False)

    Returns
    -------
    outfile: str
        Name of output MINC file.
    success: int
        Integer indicating whether conversion was a success.
    """
    
    if outfile is None:
        outfile = infile.replace('.raw', '.mnc')
    
    try: 
        rawtominc = ('cat {} '
                     '| rawtominc {} '
                     '-signed -float '
                     '-ounsigned '
                     '-oshort '
                     '-xstep 0.2 '
                     '-ystep 0.2 '
                     '-zstep 0.2 '
                     '-clobber 58 41 67'.format(infile, outfile))
        success = 1
    except: 
        success = 0
    
    os.system(rawtominc)
    
    if keep_raw is not True:
        os.remove(infile)
        
    return outfile, success


def transform_space(infile, outfile = None, voxel_orientation = 'RAS', 
                    world_space = 'MICe', expansion_factor = 1.0, 
                    volume_type = None, data_type = None, labels = False):

    """ 
    Transform the coordinate space of a MINC file

    Author: Yohan Yee

    Arguments
    ---------
    infile: str
        Name of the MINC file to transform.
    outfile: str, optional
        Name of the output MINC file. If None, the input file will be
        overwritten. 
        (default None)
    voxel_orientation: str, optional
        (default 'RAS')
    world_space: str, optional
        (default 'MICe')
    expansion_factor: float, optional
        (default 1.0)
    volume_type: 
        (default None)
    data_type: 
        (default None)
    labels: bool, optional
        (default None)

    Returns
    -------
    outfile: str
        Name of the transformed file.
    """
    
    def reorient_to_standard(dat):
        dat = np.rot90(dat, k=1, axes=(0, 2))
        dat = np.rot90(dat, k=1, axes=(0, 1))

        shape = dat.shape
        dat = np.ravel(dat)
        dat = np.reshape(dat, shape)

        return(dat)

    def do_nothing(dat):
        return(dat)

    # %% Coordinate definitions

    # Centers are listed as x,y,z; reverse these when writing out
    # Centers listed in um in CCFv3 coordinates
    centers_RAS = {"MICe"   :   [5700, 7900, 2700],
                   "CCFv3"  :   [0, 13200, 8000]}
    centers_PIR = {"MICe"   :   [5300, 5300, 5700],
                   "CCFv3"  :   [0, 0, 0]}

    # Direction cosines
    direction_cosines_RAS = {"MICe" :[[1, 0, 0], [0, 1, 0], [0, 0, 1]],
                             "CCFv3":[[0, 0, 1], [-1, 0, 0], [0, -1, 0]]}

    direction_cosines_PIR = {"MICe" :[[0, -1, 0], [0, 0, -1], [1, 0, 0]],
                             "CCFv3":[[1, 0, 0], [0, 1, 0], [0, 0, 1]]}

    # Map arguments to functions/dicts/values
    map_voxel_orientations = {"RAS":reorient_to_standard,
                              "PIR":do_nothing}

    map_centers = {"RAS":centers_RAS,
                   "PIR":centers_PIR}

    map_dir_cosines = {"RAS":direction_cosines_RAS,
                       "PIR":direction_cosines_PIR}

    size_10 = 1320*800*1140
    size_25 = 528*320*456
    size_50 = 264*160*228
    size_100 = 132*80*114
    size_200 = 58*41*67

    map_resolutions = {size_10: 10,
                       size_25: 25,
                       size_50: 50,
                       size_100: 100,
                       size_200: 200}
    
    vol = volumeFromFile(infile)

    res = map_resolutions[vol.data.size]

    # Voxel orientation
    if voxel_orientation in map_voxel_orientations:
        new_data = map_voxel_orientations[voxel_orientation](vol.data)
    else:
        print("Invalid voxel orientation")
        sys.exit(1)

    # World coordinate system
    centers = [expansion_factor*c/(1000) 
               for c in map_centers[voxel_orientation][world_space]]
    steps = [expansion_factor*res/1000] * 3
    xdc = map_dir_cosines[voxel_orientation][world_space][0]
    ydc = map_dir_cosines[voxel_orientation][world_space][1]
    zdc = map_dir_cosines[voxel_orientation][world_space][2]

    # Types
    vtype = vol.volumeType if volume_type is None else volume_type
    dtype = vol.dtype if data_type is None else data_type
    labels = vol.labels if labels is None else labels
    
    if outfile is None:
        outfile = infile
        tmpfile = infile.replace('.mnc', '')+'_tmp.mnc'
    else:
        tmpfile = outfile

    outvol = volumeFromDescription(outputFilename=tmpfile,
                                   dimnames=["zspace", "yspace", "xspace"],
                                   sizes=new_data.shape,
                                   starts=[-c for c in reversed(centers)],
                                   steps=[s for s in reversed(steps)],
                                   x_dir_cosines=xdc,
                                   y_dir_cosines=ydc,
                                   z_dir_cosines=zdc,
                                   volumeType=vtype,
                                   dtype=dtype,
                                   labels=labels)

    outvol.data = new_data
    outvol.writeFile()
    outvol.closeVolume()
    
    if outfile == infile:
        os.rename(tmpfile, infile)
        
    return outfile
    

def download_data(experiment, outdir):
    
    """
    Download and transform the data from an ISH experiment

    Arguments
    ---------
    experiment: list
        List containing ISH experiment information. Element 0 must
        contain the experiment ID as `int`. Element 1 must contain
        the gene acronym as `str`.
    outdir: str
        Directory in which to download the ISH image.

    Returns
    -------
    None
    """
    
    experiment_id = experiment[0]
    gene = experiment[1]
    
    rawfile, success = fetch_expression(experiment_id, outdir = outdir)

    if success == 1:
        mincfile, success = rawtominc_wrapper(infile = rawfile)
    
        outfile = transform_space(infile = mincfile, voxel_orientation = 'RAS',
                                  world_space = 'MICe', expansion_factor = 1.0)
    
        os.rename(outfile, outdir+'{}_{}.mnc'.format(gene, experiment_id))
    
    return
    
    
# Main -----------------------------------------------------------------------    

def main():

    #Get command line arguments
    args = parse_args()
    dataset = args['dataset']
    outdir = args['outdir']
    metadata = args['metadata']
    parallel = True if args['parallel'] == 'true' else False
    verbose = True if args['verbose'] == 'true' else False
   
    #Format the path properly
    outdir = os.path.join(outdir, '')

    #If outdir does not exist, create it
    if os.path.exists(outdir) == False:
        if verbose:
            print('Output directory {} not found. Creating it...'.format(outdir))
        os.mkdir(outdir)
        
    #If AMBA metadata file not found, download it from the web
    if os.path.isfile(outdir+metadata) == False:
        if verbose:
            print('Metadata file {} not found in {}. Fetching from API...'
                  .format(metadata, outdir))
        fetch_metadata(dataset = dataset,
                       outdir = outdir,
                       outfile = metadata,
                       verbose = verbose)
        
    #Import AMBA metadata
    dfMetadata = pd.read_csv(outdir+metadata, index_col = None)
    
    #Extract experiment IDs and gene names from metadata
    nrows = dfMetadata.shape[0]
    experiments = [(dfMetadata.loc[i, 'experiment_id'], dfMetadata.loc[i, 'gene']) 
                   for i in range(0, nrows)]
    
    #Create output sub-directory based on data set specified
    outdir = os.path.join(outdir, dataset, '')
    if os.path.exists(outdir) == False:
        os.mkdir(outdir)

    #Partial version of function for iteration
    download_data_partial = partial(download_data, outdir = outdir)
    
#     size = 1.0 if dataset == 'coronal' else 3.0

#     response_flag = 0
#     while response_flag != 1:
#         response = input(('Download will take up approximately {}GB of space.'
#                           'Proceed? (y/n) '.format(size)))
#         if response == 'y':
#             response_flag = 1
#         elif response == 'n':
#             print('Aborting.')
#             quit()
#         else:
#             print('Input not recognized (y/n): {}'.format(response))
    
    if verbose:
        print('Downloading {} AMBA dataset to: {}'.format(dataset, outdir))
        
    if parallel:
        
        nproc = args['nproc']
        pool = mp.Pool(nproc)
        
        if verbose:
            print('Running in parallel on {} CPUs...'.format(nproc))
        
        #Download data in parallel. Show progress bar.
        results_tqdm = []
        for result in tqdm(pool.imap(download_data_partial, experiments),
                           total = len(experiments)):
            results_tqdm.append(result)
            
        pool.close()
        pool.join()
        
    else:
        
        experiments_tqdm = tqdm(experiments)
        results = list(map(download_data_partial, experiments_tqdm))
        
    return
    
if __name__ == '__main__':
    main()