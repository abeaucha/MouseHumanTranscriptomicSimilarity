# MouseHumanTranscriptomicSimilarity

Antoine Beauchamp

This repository contains all of the data, code, and text necessary to fully re-create the paper "Whole-brain comparison of rodent and human brains using spatial transcriptomics" by A. Beauchamp et al. The paper is currently available as a [preprint on BioRxiv](https://doi.org/10.1101/2022.03.18.484766).

## Project overview

We explore the idea of using the brain-wide spatial expression patterns of mouse-human homologous genes to create a common space in which to make direct and quantitative comparisons between the brains of mice and humans. To build this common space, we use data sets from the Allen Institute for Brain Science, which are readily accessible from the web.

## Requirements

## Usage

All aspects of this project, including downloading and analyzing the data, generating figures, and rendering the final manuscript, were done in a programmatic manner that is automatically re-createable. 

The manuscript can be generated in its entirety from scratch using the following steps:

1. Source the `pipeline_main` shell script. This pipeline is responsible for downloading and processing the spatial transcriptomics data, as well as training the perceptron neural network. 
2. Source the `compile_manuscript` shell script in the manuscript sub-directory of interest (e.g. manuscript/submissions/biorxiv/). This script runs all analyses for the paper, generates all associated figures, and renders the final PDF. 

