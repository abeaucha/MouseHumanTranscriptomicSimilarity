# ----------------------------------------------------------------------------
# mlp_validation_naive.py
# Author: Antoine Beauchamp

"""


Description
-----------

"""

# Packages -------------------------------------------------------------------

import sys
import argparse
import os
import random
import pandas                   as pd
import numpy                    as np
from   datatable                import fread
from   itertools                import product

import torch.optim
from   torch                    import manual_seed
from   torch.optim.lr_scheduler import OneCycleLR
from   torch.cuda               import is_available

from   skorch                   import NeuralNetClassifier
from   skorch.toy               import make_classifier
from   skorch.callbacks         import LRScheduler, EpochScoring
from   skorch.helper            import DataFrameTransformer
from   skorch.dataset           import ValidSplit

from   sklearn.metrics          import accuracy_score


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
        '--nfolds',
        type = int,
        default = 5,
        help = ("Number of folds to use for cross-validation.")
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
        nargs = '*',
        type = int,
        help = ("List containing random seeds to tune over.")
    )
    
    args = vars(parser.parse_args())
    
    return args
    
    
# Functions ------------------------------------------------------------------
    
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
    
    
    # Import data -------------------------------------------------------------

    print("Importing data...")
    
    #Training data
    infile = ('MouseExpressionMatrix_'
              'voxel_coronal_maskcoronal_'
              'log2_grouped_imputed_labelled_scaled.csv')
    infile = os.path.join(datadir, infile)

    dfExprVoxel = fread(infile, header = True).to_pandas()
    

    # Process data ------------------------------------------------------------

    print("Preparing data for learning...")

    #Identify which columns contain labels
    indLabels = dfExprVoxel.columns.str.match('Region')

    #Extract matrix of gene expression values
    dfInput = dfExprVoxel.loc[:,~indLabels]
    
    labelcol = args['labels'].capitalize()
    
    print("Using labels: {}".format(labelcol))
    
    #Create a new data frame containing intermediate labels
    dfLabels = dfExprVoxel[[labelcol]].copy()

    #Convert labels to category type
    dfLabels.loc[:,labelcol] = dfLabels.loc[:,labelcol].astype('category')

    # Create an instance of the transformer
    dftx = DataFrameTransformer()

    # Fit and transform the input and label data frames
    X_temp = dftx.fit_transform(dfInput)
    y_temp = dftx.fit_transform(dfLabels)

    # Extract the arrays from the dictionaries.
    X = X_temp['X']
    y = y_temp[labelcol]
    
    
    # Initialize the network --------------------------------------------------

    print("Setting up hyperparameter grid...")
    
    #Define a dictionary containing the grid values
    dict_grid = {'nfolds':[args['nfolds']],
                 'hidden_units':args['nunits'],
                'hidden_layers':args['nlayers'],
                'dropout':args['dropout'],
                'weight_decay':args['L2'],
                'max_epochs':args['nepochs'],
                'total_steps':(args['nepochs'] if args['totalsteps'] is None else args['totalsteps']),
                'learning_rate':args['learningrate'],
                'optimizer':args['optimizer'],
                'seed':[args['seed']] if args['seed'] is None else args['seed']}
    
    df_params = pd.DataFrame([row for row in product(*dict_grid.values())], 
                             columns = dict_grid.keys())
    df_params['parameter_set'] = [i+1 for i in range(df_params.shape[0])]
    
    for index, row in df_params.iterrows():
    
        #Extract hyperparameters from the current set
        parameter_set = int(row['parameter_set'])
        nfolds = int(row['nfolds'])
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
               '  Optimizer: {}\n'
               '  Seed: {}\n'.format(parameter_set,
                                     args['labels'].title(), 
                                     hidden_units, 
                                     hidden_layers, 
                                     dropout, 
                                     weight_decay,
                                     max_epochs,
                                     total_steps,
                                     learning_rate,
                                     optimizer,
                                     seed)))
        
        optimizer = getattr(sys.modules['torch.optim'], 
                            optimizer)
            
        if seed is not None:
            np.random.seed(seed)
            manual_seed(seed)
            random.seed(seed)
            
        #Initialize the classifier module
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

        #Initialize the network optimization
        net_callbacks = [('lr_scheduler', LRScheduler(policy=OneCycleLR,
                                                      total_steps=total_steps,
                                                      cycle_momentum=False,  
                                                      max_lr=learning_rate)),
                         EpochScoring(calculate_accuracy, 
                                      use_caching = False,
                                      lower_is_better = False,
                                      on_train = True,
                                      name = 'train_acc')]
            
        net = NeuralNetClassifier(net_module,
                                  train_split = ValidSplit(nfolds),
                                  optimizer = optimizer,
                                  optimizer__weight_decay = weight_decay,
                                  max_epochs = max_epochs,
                                  callbacks = net_callbacks,
                                  device = device)
    
        #Train the network
        net.fit(X, y)
        
        #Get training history
        df_epochs = get_training_history(net)
        df_epochs['parameter_set'] = parameter_set
        
        #Combine history with parameters
        df_performance_iter = pd.merge(df_params,
                                       df_epochs,
                                       on = 'parameter_set')
        
        #Concatenate current performance data with previous
        if parameter_set == 1:
            df_performance = df_performance_iter
        else:
            df_performance = pd.concat([df_performance,
                                        df_performance_iter],
                                       axis = 0)

    #Write to file
    outfile = args['outfile']
    if outfile is None:
        outfile = 'MLP_validation_naive_{}.csv'.format(args['labels'])
    df_performance.to_csv(os.path.join(outdir, outfile), index = False)
        
        
if __name__ == '__main__':
    main()
