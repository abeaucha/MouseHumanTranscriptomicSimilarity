# ----------------------------------------------------------------------------
# build_sample_matrix.py 
# Author: Antoine Beauchamp
# Created: June 22nd, 2022

"""
Build a gene-by-sample expression matrix.

Description
-----------
This is a script to process the Allen Human Brain Atlas microarray gene 
expression data and build a gene-by-sample expression matrix. 
The processing steps are implemented using the abagen package.
"""


# Packages -------------------------------------------------------------------

import argparse
import os
import abagen
import pandas as pd


# Command line arguments -----------------------------------------------------

def parse_args():
    
    """Parse command line arguments"""
    
    parser = argparse.ArgumentParser(
        formatter_class = argparse.ArgumentDefaultsHelpFormatter
    )
    
    parser.add_argument(
        '--datadir',
        type = str,
        default = 'data/',
        help = ("Directory containing the AHBA data sets. If the data are not 
                 found, they will be downloaded.")
    )
    
    parser.add_argument(
        '--outdir',
        type = str,
        default = 'data/',
        help = "Directory in which to save the sample expression matrix."
    )
    
    parser.add_argument(
        '--donorsfile',
        type = str,
        default = 'data/donors.csv',
        help = ("Path to .csv file containing donor naming conventions.")
    )
    
    parser.add_argument(
        '--ibf-threshold',
        type = float,
        default = 0.5,
        help = ("Threshold for intensity-based filtering.")
    )
    
    parser.add_argument(
        '--probe-selection',
        type = str,
        default = 'diff_stability',
        help = ("Method used to subset multiple probes.")
    )
    
    parser.add_argument(
        '--donor-probes',
        type = str,
        default = 'aggregate',
        help = ("")
    )
    
    parser.add_argument(
        '--sim-threshold',
        type = float,
        help = ("")
    )
    
    parser.add_argument(
        '--sample-norm',
        type = str,
        default = 'srs',
        help = ("")
    )
    
    parser.add_argument(
        '--gene_norm',
        type = str,
        default = 'srs',
        help = ("")
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


# Functions ------------------------------------------------------------------

def get_sample_metadata(data_dir, donors):

    """
    Get sample metadata
    
    Arguments
    ---------
    data_dir: str
        Path to directory containing AHBA data sets.
    donors: str
        Path to .csv file containing donor naming conventions.
    
    Returns
    -------
    df_samples: pandas.core.frame.DataFrame
        Dataframe containing sample metadata.
    """
    
    donors = pd.read_csv(donors)
    sample_data = []
    for i, row in donors.iterrows():
        
        donor_dir = 'normalized_microarray_{}'.format(row['donorFileID'])
        samplefile = os.path.join(data_dir, donor_dir, 'SampleAnnot.csv')
        df_sample = pd.read_csv(samplefile)
    
        sample_ids = (df_sample
                      .apply(lambda x: '%s-%s-%s' % (x['structure_id'], 
                                                     x['slab_num'], 
                                                     x['well_id']), 
                             axis = 1))
        df_sample['SampleID'] = sample_ids
    
        df_sample['Donor'] = row['donorID']
    
        sample_data.append(df_sample)
        
    df_samples = pd.concat(sample_data, 
                           axis = 0, 
                           ignore_index = True)
    
    return df_samples
    

# Main -----------------------------------------------------------------------

def main():
    
    #Parse command line arguments
    args = parse_args()
    datadir = args['datadir']
    outdir = args['outdir']
    donorsfile = args['donorsfile']
    ibf_threshold = args['ibf_threshold']
    probe_selection = args['probe_selection']
    donor_probes = args['donor_probes']
    sim_threshold = args['sim_threshold']
    sample_norm = args['sample_norm']
    gene_norm = args['gene_norm']
    verbose = True if args['verbose'] == 'true' else False
    verbose_int = 1 if verbose else 0
    
    #Format the paths properly
    outdir = os.path.join(outdir, '')
    datadir = os.path.join(datadir, '')
    
    #Create directories if needed
    if not os.path.exists(outdir):
        os.makedirs(outdir)
        
    if not os.path.exists(datadir):
        os.makedirs(datadir)
    
    if verbose:
        print("Creating sample expression matrix...")
        
    #Build sample expression matrix
    expression, coords = abagen.get_samples_in_mask(mask=None, 
                                                    data_dir = datadir, 
                                                    donors = 'all',
                                                    ibf_threshold = ibf_threshold,
                                                    probe_selection = probe_selection,
                                                    donor_probes = donor_probes,
                                                    sim_threshold = sim_threshold,
                                                    sample_norm = sample_norm,
                                                    gene_norm = gene_norm,
                                                    verbose = verbose_int)
    
    if verbose:
        print("Getting sample metadata...")
    
    #Build sample metadata data frame
    df_samples = get_sample_metadata(data_dir = datadir,
                                     donors = donorsfile)
    
    if verbose:
        print("Matching sample expression and metadata...")
    
    #Match ordering of expression and metadata
    expression = (expression
                  .merge(df_samples.loc[:,['well_id', 'SampleID']],
                         on = 'well_id', how = 'left')
                  .drop(labels = 'well_id', axis = 1))
    
    df_samples = (df_samples
                  .set_index('SampleID')
                  .reindex(index = expression['SampleID'])
                  .reset_index(level = 'SampleID'))
    
    #Transpose expression data
    expression = expression.set_index('SampleID')
    expression.index.name = None
    expression = expression.transpose().sort_index()
    
    if verbose:
        print("Writing sample expression matrix to file...")
    
    #Write expression matrix to file
    exprfile = 'HumanExpressionMatrix_samples_pipeline_abagen.csv'
    exprfile = os.path.join(outdir, exprfile)
    expression.to_csv(exprfile, index_label = 'Gene')
    
    if verbose:
        print("Writing sample metadata to file...")
    
    #Write metadata to file
    samplefile = 'SampleInformation_pipeline_abagen.csv'
    samplefile = os.path.join(outdir, samplefile)
    df_samples.to_csv(samplefile, index = False)

    return
    
    
if __name__=='__main__':
    main()
