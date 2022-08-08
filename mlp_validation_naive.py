# ----------------------------------------------------------------------------
# train_multilayer_perceptron.py 
# Author: Antoine Beauchamp
# Created: December 15th, 2020

"""
Train a multi-layer perceptron

Description
-----------
This script trains a multi-layer perceptron neural network to classify mouse
voxels into a specified number of neuroanatomical atlas regions. 

Once the network is trained, the final layer is removed to obtain the
transformation from the input space to the latent space defined by the last
hidden layer. 

The modified network architecture is used to transform aggregated regional
mouse and human gene expression matrices from the input space into the latent
space. An option also exists to transform voxel- and sample-wise expression
matrices as well.
"""

# Packages -------------------------------------------------------------------

import sys

import pandas                 as pd
import numpy                  as np
import random
import argparse
import os
from datatable                import fread
from itertools                import product

import torch
import torch.nn.functional    as F
from torch                    import nn
from torch.optim              import SGD, AdamW
from torch.optim.lr_scheduler import OneCycleLR
from torch.cuda               import is_available

from skorch                   import NeuralNetClassifier
from skorch.callbacks         import LRScheduler, EpochScoring
from skorch.helper            import DataFrameTransformer

from sklearn.metrics          import accuracy_score, confusion_matrix


# Functions ------------------------------------------------------------------

def parse_args():
    
    parser = argparse.ArgumentParser(
        formatter_class = argparse.ArgumentDefaultsHelpFormatter
    )
    
    parser.add_argument(
        "--datadir",
        type = str,
        default = 'data/',
        help = "Directory containing input data."
    )
    
    parser.add_argument(
        "--outdir",
        type = str,
        default = 'data/MLP_outcomes/',
        help = "Directory in which to write neural net outcomes."
    )
    
    parser.add_argument(
        "--labels",
        type = str,
        default = 'region5',
        choices = ['region5', 
                   'region11',
                   'region28',
                   'region46',
                   'region67',
                   'region134'],
        help = "Class of mouse labels on which to train."
    )
    
    parser.add_argument(
        "--nunits",
        nargs = "*",
        type = int,
        default = [500],
        help = "List containing the number of hidden units to tune over."
    )
    
    parser.add_argument(
        "--L2",
        nargs = "*",
        type = float,
        default = [1e-6],
        help = "List containing weight decay values to tune over."
    )
    
    parser.add_argument(
        "--nepochs",
        type = int,
        default = 200,
        help = "Number of epochs to train over."
    )
    
    parser.add_argument(
        "--learningrate",
        type = float,
        default = 1e-5,
        help = "Learning rate during training."
    )
    
    parser.add_argument(
        '--totalsteps',
        type = int,
        help = "Number of steps to use in optimizer."
    )
    
    parser.add_argument(
        '--optimizer',
        type = str,
        default = 'SGD',
        choices = ['SGD', 'AdamW'],
        help = ("Neural network optimizer.")
    )
    
    parser.add_argument(
        "--confusionmatrix",
        type = str,
        default = 'false',
        choices = ['true', 'false'],
        help = ("Option to compute confusion matrix from training set "
                "predictions.")
    )
    
    parser.add_argument(
        '--seed',
        type = int,
        help = ("Random seed")
    )
    
    args = vars(parser.parse_args())
    
    return args
    
    
def calculate_accuracy(net, X, y):
    y_pred = net.predict(X)
    acc = accuracy_score(y, y_pred)
    return acc

#Define network architecture
class ClassifierModule(nn.Module):
        
        def __init__(
            self,
            input_units,
            output_units,
            hidden_units,
            apply_output_layer = True #Flag to apply output layer
        ):
            super(ClassifierModule, self).__init__()

            self.apply_output_layer = apply_output_layer

            self.hidden1 = nn.Linear(input_units, hidden_units)
            self.hidden2 = nn.Linear(hidden_units, hidden_units)
            self.hidden3 = nn.Linear(hidden_units, hidden_units)
            self.output = nn.Linear(hidden_units, output_units)

        def forward(self, X, **kwargs):
            X = F.relu(self.hidden1(X))
            X = F.relu(self.hidden2(X))
            X = F.relu(self.hidden3(X))

            #If flag is True, apply output layer
            if self.apply_output_layer is True:
                X = F.softmax(self.output(X), dim = -1)

            return X

    
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
        os.mkdir(outdir)
    
    # Import data -------------------------------------------------------------

    #Set up files for import
    #Mouse voxelwise data to train over
    file_voxel = ("MouseExpressionMatrix_"
                  "voxel_coronal_maskcoronal_"
                  "log2_grouped_imputed_labelled_scaled.csv")
    filepath_voxel = os.path.join(datadir, file_voxel)
    
    print("Importing data...")

    #Import data
    dfExprVoxel = (fread(filepath_voxel, header = True)
                   .to_pandas())

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

    print("Initializing neural network...")
    
    #Get network parameters from command line args
    dict_grid = {'hidden_units':args['nunits'],
                 'weight_decay':args['L2']}
    
    df_params = pd.DataFrame([row for row in product(*dict_grid.values())], 
                             columns = dict_grid.keys())
    
    max_epochs = args['nepochs']
    total_steps = args['totalsteps']
    learning_rate = args['learningrate']
    optimizer = args['optimizer']
    seed = args['seed']
    
    if total_steps is None:
        total_steps = max_epochs
        
    df_params['max_epochs'] = max_epochs
    df_params['totalsteps'] = total_steps
    df_params['learningrate'] = learning_rate
    df_params['optimizer'] = optimizer
        
    if optimizer == 'AdamW':
        optimizer = AdamW
    elif optimizer == 'SGD':
        optimizer = SGD
    else:
        raise ValueError
    
    if is_available() == True:
        print("GPU available. Training network using GPU...")
        device = 'cuda'
    else:
        print("GPU unavailable. Training network using CPU...")
        device = 'cpu'
        
    for i, row in df_params.iterrows():

        hidden_units = int(row['hidden_units'])
        weight_decay = row['weight_decay']
        
        if seed is not None:
            np.random.seed(seed)
            torch.manual_seed(seed)
            random.seed(seed)
        
        net_module = ClassifierModule(input_units = X.shape[1],
                                      output_units = len(np.unique(y)),
                                      hidden_units = hidden_units)

    #Create the classifier
        net = NeuralNetClassifier(net_module,
                                  optimizer = optimizer,
                                  optimizer__weight_decay = weight_decay,
                                  max_epochs = max_epochs,
                                  callbacks = [('lr_scheduler',
                                                LRScheduler(policy=OneCycleLR,
                                                            total_steps=total_steps,
                                                            cycle_momentum=False,  
                                                            max_lr=learning_rate)),
                                               EpochScoring(calculate_accuracy, 
                                                            use_caching = False,
                                                            lower_is_better = False,
                                                            on_train = True,
                                                            name = 'train_acc')],
                                  device = device)

    
    # Train the network ------------------------------------------------------
    
    #Fit the network
        net.fit(X, y)
    
        for i in range(max_epochs):
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
    
        
        df_performance = pd.concat([df_params, 
                                    df_epochs],
                                   axis = 1)
        
    outfile = 'MLP_validation_naive_{}.csv'.format(args['labels'])
    df_performance.to_csv(os.path.join(outdir, outfile), index = False)
        
    
#     #Match the dummy variable labels to the region names
#     dfLabels['y'] = y
#     dfLabelsUnique = dfLabels.sort_values('y').drop_duplicates()
    
    
    # Compute training confusion matrix --------------------------------------
    
#     #Switch to compute confusion matrix
#     if args['confusionmatrix'] == 'true':
    
#         print("Computing confusion matrix from training set...")
    
#         #Compute confusion matrix and store as data frame
#         dfConfusionMat = pd.DataFrame(confusion_matrix(y, y_pred))
    
#         #Assign region names to the confusion matrix columns
#         dfConfusionMat.columns = dfLabelsUnique[labelcol].astype('str')
    
#         #Assign region names to the confusion matrix rows
#         dfConfusionMat['TrueLabels'] = (dfLabelsUnique[labelcol]
#                                         .astype('str')
#                                         .reset_index(drop = True))
        
#         #File to save confusion matrix
#         fileConfMat = "MLP_ConfusionMatrix_Training"+\
#         "_"+args['labels'].capitalize()+\
#         "_Layers3"+\
#         "_Units"+str(args['nunits'])+\
#         "_L2"+str(args['L2'])+".csv"
    
#         #Write confusion matrix to file
#         dfConfusionMat.to_csv(os.path.join(outdir, fileConfMat),
#                               index = False)
        
        
if __name__ == "__main__":
    main()
