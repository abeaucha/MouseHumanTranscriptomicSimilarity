# MouseHumanTranscriptomicSimilarity

Antoine Beauchamp

This repository contains all of the data, code, and text necessary to fully re-create the eLife article:

Beauchamp A, Yee Y, Darwin BC, Mars RB, Lerch JP. 2022. _Whole-brain comparison of rodent and human brains using spatial transcriptomics_. eLife e79418. 

The article is [openly accessible via the eLife website](https://elifesciences.org/articles/79418).


## Project overview

We explore the idea of using the brain-wide spatial expression patterns of mouse-human homologous genes to create a common space in which to make direct and quantitative comparisons between the brains of mice and humans. To build this common space, we use data sets from the Allen Institute for Brain Science, which are readily accessible from the web.

## Requirements

## Usage

All aspects of this project, including downloading and analyzing the data, generating figures, and rendering the final manuscript, were done in a programmatic manner that is automatically re-createable. 

The manuscript can be generated in its entirety from scratch using the following steps:

1. Source the `pipeline_main` shell script. This pipeline is responsible for downloading and processing the spatial transcriptomics data, as well as training the perceptron neural network. 
2. Source the `compile_manuscript` shell script in the manuscript sub-directory of interest (e.g. manuscript/submissions/biorxiv/). This script runs all analyses for the paper, generates all associated figures, and renders the final PDF. 

