import os
import argparse
import numpy as np
import pandas as pd
import requests
from pyminc.volumes.factory import *
from zipfile import ZipFile
from multiprocessing import Pool
from functools import partial
from itertools import starmap


def parse_args():

    parser = argparse.ArgumentParser()

    parser.add_argument(
        '--dataset',
        type = str,
        default = 'coronal',
        choices = ['coronal', 'sagittal'],
        help = 'AMBA dataset to import'
    )

    parser.add_argument(
        '--datadir',
        type = str,
        default = 'data/expression/',
        help = 'Path to directory in which to download the data'
    )
    
    parser.add_argument(
        '--metadata',
        type = str,
        default = 'AMBA_metadata.csv',
        help = 'File containing AMBA metadata'
    )
    
    parser.add_argument(
        '--parallel',
        type = str,
        default = 'false',
        choices = ['true', 'false'],
        help = 'Run in parallel'
    )
    
    parser.add_argument(
        '--nproc',
        type = int,
        default = 1
        help = 'Number of CPUs to use in parallel'
    )

    args = vars(parser.parse_args())

    return args


def fetch_metadata(dataset = 'coronal', outdir='./', outfile = 'AMBA_metadata.csv'):

    """ """        

    abi_query_metadata = "http://api.brain-map.org/api/v2/data/SectionDataSet/query.csv?"+\
"criteria=[failed$eqfalse],plane_of_section[name$eq{}],products[abbreviation$eqMouse],treatments[name$eqISH],genes&".format(dataset)+\
"tabular=data_sets.id+as+experiment_id,data_sets.section_thickness,data_sets.specimen_id,"+\
"plane_of_sections.name+as+plane,"+\
"genes.acronym+as+gene,genes.name+as+gene_name,genes.chromosome_id,genes.entrez_id,genes.genomic_reference_update_id,genes.homologene_id,genes.organism_id&"+\
"start_row=0&num_rows=all"

    pd.read_csv(abi_query_metadata).to_csv(outdir+outfile, index=False)

    print('Metadata downloaded at: {}'.format(outdir+outfile))
    
    return 


def fetch_expression(experiment_id, outdir = './tmp/'):

    """ """
    
    if not os.path.exists(outdir):
        os.mkdir(outdir)

    tmpdir = outdir+str(experiment_id)+'/' 
    os.mkdir(tmpdir)
        
    abi_query_expr = 'http://api.brain-map.org/grid_data/download/{}'.format(experiment_id)
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
            print('Error for experiment {}: {}'.format(experiment_id, err))
            success = 0
            
    os.remove(tmpfile)
    os.rmdir(tmpdir)
 
    return outfile, success


def rawtominc_wrapper(infile, outfile = None, keep_raw = False):
    
    """ """
    
    if outfile is None:
        outfile = infile.replace('.raw', '.mnc')
    
    try: 
        rawtominc = 'cat {} | rawtominc {} -signed -float -ounsigned -oshort -xstep 0.2 -ystep 0.2 -zstep 0.2 -clobber 58 41 67'.format(infile, outfile)
        success = 1
    except: 
        success = 0
    
    os.system(rawtominc)
    
    if keep_raw is not True:
        os.remove(infile)
        
    return outfile, success


def transform_space(infile, outfile = None, voxel_orientation = 'RAS', world_space = 'MICe', expansion_factor = 1.0, volume_type = None, data_type = None, labels = False):

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
    direction_cosines_RAS = {"MICe"     :   [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
                             "CCFv3"    :   [[0, 0, 1], [-1, 0, 0], [0, -1, 0]]}

    direction_cosines_PIR = {"MICe"     :   [[0, -1, 0], [0, 0, -1], [1, 0, 0]],
                             "CCFv3"    :   [[1, 0, 0], [0, 1, 0], [0, 0, 1]]}

    # Map arguments to functions/dicts/values
    map_voxel_orientations = {"RAS" :   reorient_to_standard,
                              "PIR":    do_nothing}

    map_centers = {"RAS"    :   centers_RAS,
                   "PIR"    :   centers_PIR}

    map_dir_cosines = {"RAS"    :   direction_cosines_RAS,
                       "PIR"    :   direction_cosines_PIR}

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

    
def download_data(experiment_id, gene, outdir):
    
    rawfile, success = fetch_expression(experiment_id, outdir = outdir)

    mincfile, success = rawtominc_wrapper(infile = rawfile)
    
    outfile = transform_space(infile = mincfile, voxel_orientation = 'RAS', world_space = 'MICe', expansion_factor = 1.0)
    
    os.rename(outfile, outdir+'{}_{}.mnc'.format(gene, experiment_id))
    
    return


    
def main():

    args = parse_args()
    dataset = args['dataset']
    datadir = args['datadir']
    metadata = args['metadata']
    parallel = if args['parallel'] == True
    
    if os.path.isfile(datadir+metadata) == False:
        print('AMBA metadata file {} not found in {}. Fetching from API.'.format(metadata,datadir))
        fetch_metadata(dataset = dataset,
                       outdir = datadir,
                       outfile = metadata)
        
    dfMetadata = pd.read_csv(datadir+metadata, index_col=None)
    
    experiment_ids = [dfMetadata.loc[i, 'experiment_id'] for i in range(0, dfMetadata.shape[0])]
    genes = [dfMetadata.loc[i, 'gene'] for i in range(0, dfMetadata.shape[0])]
    
    experiments = [(dfMetadata.loc[i, 'experiment_id'], dfMetadata.loc[i, 'gene']) for i in range(0, dfMetadata.shape[0])]
    
    outdir = datadir+dataset+'/'
    download_data_partial = partial(download_data, outdir = outdir)

    nproc = 4
    parallel = True
    if parallel:
        pool = Pool(nproc)
        results = pool.starmap(download_data_partial, experiments[:20])
        pool.close()
        pool.join()
    else:
        results = list(starmap(download_data_partial, experiments[:5]))
    
    
#     experiment_id = experiment_ids[0]
#     gene = genes[0]
    
#     rawfile, success = fetch_expression(experiment_id, outdir = outdir)

#     mincfile, success = rawtominc_wrapper(infile = rawfile)
    
#     outfile = transform_space(infile = mincfile, voxel_orientation = 'RAS', world_space = 'MICe', expansion_factor = 1.0)
    
#     os.rename(outfile, outdir+'{}_{}.mnc'.format(gene, experiment_id))

#     nproc = 4
#     pool = Pool(nproc)
#     results = pool.starmap(download_data_partial, zip(experiment_ids[:20], genes[:20]))
#     pool.close()
#     pool.join()
 
    


if __name__ == '__main__':
    main()