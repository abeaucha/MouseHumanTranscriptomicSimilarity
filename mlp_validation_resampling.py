# ----------------------------------------------------------------------------
# mlp_validation_resampling.py
# Author: Antoine Beauchamp
# Created: February 3rd, 2021

"""


Description
-----------

"""

# Packages -------------------------------------------------------------------

import sys
import argparse
import os
import random
import pandas                 as pd
import numpy                  as np
from   datatable              import fread
from   itertools              import product

from sklearn.impute           import SimpleImputer
from sklearn.preprocessing    import StandardScaler, FunctionTransformer
from sklearn.pipeline         import Pipeline
from sklearn.metrics          import accuracy_score, confusion_matrix

from skorch                   import NeuralNetClassifier
from skorch.toy               import make_classifier
from skorch.helper            import DataFrameTransformer
from skorch.callbacks         import LRScheduler, EpochScoring
from skorch.dataset           import Dataset

import torch.optim
from torch                    import manual_seed
from torch.optim.lr_scheduler import OneCycleLR
from torch.cuda               import is_available
from torch.utils.data         import Subset


# Command line arguments ------------------------------------------------------------------

def parse_args():
    
    parser = argparse.ArgumentParser(
        formatter_class = argparse.ArgumentDefaultsHelpFormatter
    )
    
    parser.add_argument(
        '--datadir',
        type = str,
        default = 'data/',
        help = "Directory containing input data."
    )
    
    parser.add_argument(
        '--outdir',
        type = str,
        default = 'data/MLP_validation/',
        help = "Directory in which to export performance data."
    )
    
    parser.add_argument(
        '--outfile',
        type = str,
        help = "Name of CSV file in which to export performance data."
    )
    
    parser.add_argument(
        '--labels',
        type = str,
        default = 'region5',
        choices = ['region5', 
                   'region11',
                   'region28',
                   'region46',
                   'region67',
                   'region134'],
        help = "Class of mouse labels on which to train the network."
    )
    
    parser.add_argument(
        '--nsamples',
        type = int,
        default = 1,
        help = ("Number of times to sample training and validation data sets.")
    )
    
    parser.add_argument(
        '--nunits',
        nargs = '*',
        type = int,
        default = [200],
        help = "List containing the number of hidden units to tune over."
    )

    parser.add_argument(
        '--nlayers',
        nargs = '*',
        type = int,
        default = [3],
        help = "List containing the number of hidden layers to tune over."
    )

    parser.add_argument(
        '--dropout',
        nargs = '*',
        type = float,
        default = [0.0],
        help = "List containing dropout rates to tune over."
    )

    parser.add_argument(
        '--L2',
        nargs = '*',
        type = float,
        default = [0.0],
        help = "List containing weight decay values to tune over."
    )

    parser.add_argument(
        '--nepochs',
        nargs = '*',
        type = int,
        default = [200],
        help = "List containing the number of epochs to train over."
    )

    parser.add_argument(
        '--learningrate',
        nargs = '*',
        type = float,
        default = [1e-5],
        help = "List containing maximum learning rates to tune over."
    )
    
    parser.add_argument(
        '--totalsteps',
        nargs = '*',
        type = int,
        help = "List containing the total number of optimization steps to tune over."
    )
    
    parser.add_argument(
        '--optimizer',
        nargs = '*',
        type = str,
        default = ['SGD'],
        help = "List containing torch.optim algorithms to tune over."
    )

    parser.add_argument(
        '--seed',
        type = int,
        help = ("Random seed for network training.")
    )
    
    args = vars(parser.parse_args())

    return args
    
    
# Functions ------------------------------------------------------------------
    
def buildTrainValidationSets(coronal, sagittal, seed = None):
    
    """ 
    Generate training and validation data sets.
    
    Description
    -----------
    
    Arguments
    ---------
    coronal:
    sagittal:
    seed:
    
    Returns
    -------
    train:
    validation:
    """
    
    #Get genes in the coronal data set (includes duplicates)
    genes_coronal = coronal.columns.str.replace('\.\.\.[0-9]+', '', regex = True)
    
    #Get unique genes, present in sagittal and coronal data
    genes_unique = np.unique(genes_coronal)
    
    #Initialize data frames for training and validation sets
    train = pd.DataFrame(np.empty((coronal.shape[0], len(genes_unique)), dtype = 'float'),
                         columns = genes_unique)
    
    validation = pd.DataFrame(np.empty_like(train),
                              columns = genes_unique)
    
    #Initialize random number generator
    rng = np.random.default_rng(seed = seed)
    
    #Iterate over unique genes
    for i, gene in enumerate(genes_unique):
        
        #If the gene has a replicated coronal experiment, use that
        if np.sum(genes_coronal == gene) > 1:
            
            #Extract replicated experiments
            df_choices = coronal.loc[:, genes_coronal == gene]
            
            #Randomly choose one of the experiments for the training set
            choices = np.arange(0, df_choices.shape[1])
            choice_train = rng.choice(choices, 1)
            
            #Randomly choose one of the remaining experiments for the validation set
            choice_val = rng.choice(choices[choices != choice_train], 1)
            
            #Assign data to training and validation sets
            train.iloc[:,i] = df_choices.iloc[:, choice_train]
            validation.iloc[:, i] = df_choices.iloc[:, choice_val]
            
        #If the gene is unique in the coronal set, choose between coronal and sagittal
        else:
            
            #Random binary choice
            choice_train = rng.choice([0,1],1)
            
            #Assign coronal data to training and sagittal to validation, or vice versa
            if choice_train[0] == 0:
                train.iloc[:,i] = coronal.loc[:, coronal.columns == gene]
                validation.iloc[:,i] = sagittal.loc[:, sagittal.columns == gene]
            else:
                train.iloc[:,i] = sagittal.loc[:, sagittal.columns == gene]
                validation.iloc[:,i] = coronal.loc[:, coronal.columns == gene]
                
    return train, validation


def train_val_split(dataset, y):
    
    """
    """
    
    train_ind = [i for i in range(int(len(dataset)/2))]
    val_ind = [i for i in range(int(len(dataset)/2),len(dataset))] 
    dataset_train = Subset(dataset, train_ind)
    dataset_val = Subset(dataset, val_ind)
    return dataset_train, dataset_val


def calculate_accuracy(net, X, y):
    
    """ 
    Compute the prediction accuracy.
    
    Arguments
    ---------
    net:
        The neural network classifier.
    X: numpy.ndarray
        The input data tensor.
    y: numpy.ndarray
        The true labels.

    Returns
    -------
    acc: numpy.float64
        The accuracy score.
    """
    
    y_pred = net.predict(X)
    acc = accuracy_score(y, y_pred)
    return acc


def get_training_history(net):

    """ 
    Obtain training history data.
    
    Arguments
    ---------
    net:
        The neural network classifier.
        
    Returns
    -------
    df_epochs: pandas.core.frame.DataFrame
        Data frame containing history data.
    """

    for i in range(len(net.__dict__['history_'])):
        epoch_dict = net.__dict__['history_'][i]
        if i == 0:
            df_epochs = (pd.DataFrame(epoch_dict)
                         .drop(columns = 'batches')
                         .drop_duplicates())
        else:
            df_epochs_tmp = (pd.DataFrame(epoch_dict)
                             .drop(columns = 'batches')
                             .drop_duplicates())
            df_epochs = pd.concat([df_epochs, 
                                   df_epochs_tmp], 
                                  axis = 0)
    return df_epochs


# Main ------------------------------------------------------------------------

def main():
    
    #Load command line arguments
    args = parse_args()
    
    datadir = args['datadir']
    outdir = args['outdir']
    
    datadir = os.path.join(datadir, '')
    outdir = os.path.join(outdir, '')
    
    if os.path.exists(outdir) == False:
        print('Output directory {} not found. Creating it...'.format(outdir))
        os.makedirs(outdir)
   
    
    
    # Importing data -------------------------------------------------------------------

    print("Importing data...")
    
    #Input files
    file_ref = ('MouseExpressionMatrix_'
                'voxel_coronal_maskcoronal_'
                'log2_grouped_imputed_labelled.csv')
    file_sagittal = ('MouseExpressionMatrix_'
                     'voxel_sagittal_masksagittal_'
                     'log2_grouped_imputed_labelled.csv')
    file_coronal = ('MouseExpressionMatrix_'
                    'voxel_coronal_masksagittal_'
                    'log2_imputed_labelled.csv')

    file_ref = os.path.join(datadir, file_ref)
    file_sagittal = os.path.join(datadir, file_sagittal)
    file_coronal = os.path.join(datadir, file_coronal)

    df_ref = (fread(file_ref, header = True)
                  .to_pandas())
    df_sagittal = (fread(file_sagittal, header = True)
                  .to_pandas())
    df_coronal = (fread(file_coronal, header = True)
                  .to_pandas())


    # Processing -------------------------------------------------------------------------

    print("Cleaning and preparing data...")

    #Identify and remove label columns
    ind_coronal_labels = df_coronal.columns.str.match('Region')
    df_input_coronal = df_coronal.loc[:, ~ind_coronal_labels]

    ind_sagittal_labels = df_sagittal.columns.str.match('Region')
    df_input_sagittal = df_sagittal.loc[:, ~ind_sagittal_labels]

    #Obtain reference gene set
    ind_ref_labels = df_ref.columns.str.contains('Region')
    genes_ref = df_ref.columns[~ind_ref_labels].str.replace('\.\.\.[0-9]+', '', regex = True)
    genes_ref = set(genes_ref)

    #Extract genes in coronal and sagittal data sets
    genes_coronal = df_input_coronal.columns.str.replace('\.\.\.[0-9]+', '', regex = True)
    genes_sagittal = df_input_sagittal.columns.str.replace('\.\.\.[0-9]+', '', regex = True)

    genes_coronal_set = set(genes_coronal)
    genes_sagittal_set = set(genes_sagittal)

    #Compute gene set intersection
    genes = set.intersection(genes_ref, genes_coronal_set, genes_sagittal_set)

    #Filter sagittal and coronal data sets for genes in the intersection
    ind_sagittal_genes = genes_sagittal.isin(genes)
    df_input_sagittal = df_input_sagittal.loc[:, ind_sagittal_genes]

    ind_coronal_genes = genes_coronal.isin(genes)
    df_input_coronal = df_input_coronal.loc[:, ind_coronal_genes]

    #Extract labels
    labelcol = args['labels'].title()
    df_labels = df_coronal[[labelcol]].copy()

    #Convert labels to category
    df_labels.loc[:,labelcol] = df_labels.loc[:,labelcol].astype('category')

    
    
    # ---------

    print("Setting up hyperparameter grid...")
    
    #Define a dictionary containing the grid values
    dict_grid = {'sample':[i for i in range(1, args['nsamples']+1)],
                 'hidden_units':args['nunits'],
                'hidden_layers':args['nlayers'],
                'dropout':args['dropout'],
                'weight_decay':args['L2'],
                'max_epochs':args['nepochs'],
                'total_steps':(args['nepochs'] if args['totalsteps'] is None else args['totalsteps']),
                'learning_rate':args['learningrate'],
                'optimizer':args['optimizer'],
                'seed':[args['seed']]}

    #Expand the dictionary grid into a data frame containing all combinations
    df_params = pd.DataFrame([row for row in product(*dict_grid.values())], 
                             columns = dict_grid.keys())
    
    #Iterate over unique samples
    for sample in np.unique(df_params['sample']):

        print("On sample {}".format(sample))

        print("Generating training and validation sets...")

        #For the given sample, build the training and validation sets
        df_training, df_validation = buildTrainValidationSets(df_input_coronal,
                                                              df_input_sagittal,
                                                              seed = sample)    


        print("Preprocessing data...")

        #Initialize classes for imputing, scaling, centering and transposing
        scale = StandardScaler()
        center = StandardScaler(with_std = False)
        transpose = FunctionTransformer(np.transpose)
        
        processing_pipeline = Pipeline([('transpose1', transpose),
                                        ('scale', scale),
                                        ('transpose2', transpose),
                                        ('center', center)])
        
        #Fit the pipeline to the training data and transform
        X_train = processing_pipeline.fit_transform(df_training.to_numpy())
        X_val = processing_pipeline.fit_transform(df_validation.to_numpy())
        
        X = np.concatenate((X_train, X_val), axis = 0)
        X = X.astype(np.float32)
        
        #Transform labels into data preferred by the network
        y = DataFrameTransformer().fit_transform(df_labels)[labelcol]
        y = np.concatenate((y, y))
        
        #For the given training/validation sample, extract hyperparameters to iterate over
        df_params_sample = df_params[df_params['sample'] == sample].copy()
        df_params_sample['parameter_set'] = [i+1 for i in range(df_params_sample.shape[0])]

        #Iterate over hyperparameter combinations
        for index, row in df_params_sample.iterrows():

            #Extract hyperparameters
            parameter_set = int(row['parameter_set'])
            hidden_units = int(row['hidden_units'])
            hidden_layers = int(row['hidden_layers'])
            dropout = row['dropout']
            weight_decay = row['weight_decay']
            max_epochs = row['max_epochs']
            total_steps = row['total_steps']
            learning_rate = row['learning_rate']
            optimizer = row['optimizer']
            seed = row['seed']
            
            print(('\nParameter set {}\n'
                   '  Labels: {}\n'
                   '  Hidden units: {}\n' 
                   '  Hidden layers: {}\n' 
                   '  Dropout: {}\n'
                   '  L2: {}\n'
                   '  Max epochs: {}\n'
                   '  Total steps: {}\n'
                   '  Learning rate: {}\n'
                   '  Optimizer: {}\n'.format(parameter_set,
                                              args['labels'].title(), 
                                            hidden_units, 
                                            hidden_layers, 
                                            dropout, 
                                            weight_decay,
                                            max_epochs,
                                            total_steps,
                                            learning_rate,
                                            optimizer)))
            
            optimizer = getattr(sys.modules['torch.optim'], 
                                optimizer)
            
            if seed is not None:
                np.random.seed(seed)
                manual_seed(seed)
                random.seed(seed)
            
            #Generate classifier module with specified architecture
            net_module = make_classifier(input_units = X.shape[1],
                                        output_units = len(np.unique(y)),
                                        hidden_units = hidden_units,
                                        num_hidden = hidden_layers,
                                        dropout = dropout)
            
            if is_available() == True:
                print("GPU available. Training network using GPU ...")
                device = 'cuda'
            else:
                print("GPU unavailable. Training network using CPU ...")
                device = 'cpu'
                
            net_callbacks = [('lr_scheduler',
                                          LRScheduler(policy=OneCycleLR,
                                                      total_steps=total_steps, 
                                                      cycle_momentum=False,   
                                                      max_lr=learning_rate)),
                                        EpochScoring(calculate_accuracy, 
                                                     use_caching = False,
                                                    lower_is_better = False,
                                                    on_train = True,
                                                    name = 'train_acc')] 
            
            net = NeuralNetClassifier(net_module,
                                      train_split = train_val_split,
                                      optimizer = optimizer,
                                      optimizer__weight_decay = weight_decay, 
                                      max_epochs = max_epochs,
                                      callbacks = net_callbacks,
                                      device = device)

            #Fit the network to the training data
            net.fit(X, y)
    
            #Get training history
            df_epochs_iter = get_training_history(net)
            df_epochs_iter['parameter_set'] = parameter_set
    
             #Combine history with parameters
            df_performance_iter = pd.merge(df_params_sample, 
                                           df_epochs_iter, 
                                           on = 'parameter_set')
            
            #Concatenate current iteration with previous
            if parameter_set == 1:
                df_performance_sample = df_performance_iter
            else:
                df_performance_sample = pd.concat([df_performance_sample,
                                                   df_performance_iter],
                                                  axis = 0)
            
        #Concatenate current sample with previous
        if sample == 1:
            df_performance = df_performance_sample
        else:
            df_performance = pd.concat([df_performance, 
                                        df_performance_sample], 
                                       axis = 0)

    #Write to file
    outfile = args['outfile']
    if outfile is None:
        outfile = 'MLP_validation_resampling_{}.csv'.format(args['labels'])
    df_performance.to_csv(os.path.join(outdir, outfile), index=False)
    
    
if __name__ == '__main__':
    main()
