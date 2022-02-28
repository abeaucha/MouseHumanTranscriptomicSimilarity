import os
import argparse
import numpy as np
import pandas as pd
import requests
from pyminc.volumes.factory import *
from zipfile import ZipFile


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

    args = vars(parser.parse_args())

    return args



def fetch_metadata(dataset = 'coronal', outdir='./', outfile = None):

    """ """

    if outfile is None:
        outfile = 'AMBA_metadata_{}.csv'.format(dataset)
        outfile = outdir+outfile
    else:
        outfile = outdir+outfile
            

    if os.path.isfile(outfile):
        metadata = pd.read_csv(outfile, index_col = None)
    else:
        abi_query_metadata = "http://api.brain-map.org/api/v2/data/SectionDataSet/query.csv?"+\
"criteria=[failed$eqfalse],plane_of_section[name$eq{}],products[abbreviation$eqMouse],treatments[name$eqISH],genes&".format(dataset)+\
"tabular=data_sets.id+as+experiment_id,data_sets.section_thickness,data_sets.specimen_id,"+\
"plane_of_sections.name+as+plane,"+\
"genes.acronym+as+gene,genes.name+as+gene_name,genes.chromosome_id,genes.entrez_id,genes.genomic_reference_update_id,genes.homologene_id,genes.organism_id&"+\
"start_row=0&num_rows=all"

        metadata = pd.read_csv(abi_query_metadata)
        metadata.to_csv(outfile, index=False)

    return metadata, outfile



def fetch_expression(experiment_id, outdir = './tmp/'):

    """ """

    if not os.path.exists(outdir):
        os.mkdir(outdir)

    abi_query_expr = 'http://api.brain-map.org/grid_data/download/{}'.format(experiment_id)

    amba_request = requests.get(abi_query_expr)

    tmpfile = outdir+str(experiment_id)+'.zip'
    with open(tmpfile, 'wb') as file:
        file.write(amba_request.content)

    with ZipFile(tmpfile, 'r') as file:
        try:
            file.extract('energy.raw', path = outdir)
            os.rename(outdir+'energy.raw', outdir+str(experiment_id)+'.raw')
            success = 1
        except KeyError as err:
            print('Error for experiment {}: {}'.format(experiment_id, err))
            success = 0
            
    os.remove(tmpfile)
 
    return success


def transform_space(infile, outfile, voxel_orientation = 'RAS', world_space = 'MICe', expansion_factor = 1.0, volume_type = None, data_type = None, labels = False):

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

    outvol = volumeFromDescription(outputFilename=outfile,
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


    
def main():

    args = parse_args()

    dataset = args['dataset']
    datadir = args['datadir']

    dfMetadata, metadatafile = fetch_metadata(dataset = dataset, 
                                              outdir = datadir)
    
    print(metadatafile)
    
    dfMetadata['success'] = 0

    dfTemp = dfMetadata.loc[:5]
    for index, row in dfTemp.iterrows():

        experiment_id = row['experiment_id']

        outdir = datadir+'{}/'.format(dataset)
        success = fetch_expression(experiment_id = experiment_id, 
                                   outdir = outdir)

        
        if bool(success):
            
            dfMetadata.loc[index,'success'] = 1
        
            infile = outdir+'{}.raw'.format(experiment_id)
            outfile = outdir+'{}_tmp.mnc'.format(experiment_id)

            cmd = 'cat {} | rawtominc {} -signed -float -ounsigned -oshort -xstep 0.2 -ystep 0.2 -zstep 0.2 -clobber 58 41 67'.format(infile, outfile)

            os.system(cmd)
            os.remove(infile)

            infile = outfile
            outfile = outdir+'{}_tmp2.mnc'.format(experiment_id)
            transform_space(infile = infile,
                            outfile = outfile,
                            voxel_orientation = 'RAS',
                            world_space = 'MICe',
                            expansion_factor = 1.0)
            os.remove(infile)
    
    
            gene_id = row['gene']
    
            infile = outfile
            outfile = outdir+'{}_{}.mnc'.format(gene_id, experiment_id)
            os.rename(infile, outfile)

    dfMetadata.to_csv(metadatafile)


if __name__ == '__main__':
    main()