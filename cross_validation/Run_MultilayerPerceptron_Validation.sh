#!/bin/bash -l

#PBS -N MLP_Region67_Validation
#PBS -q gpu
#PBS -l nodes=1:ppn=4:gpus=1
#PBS -l mem=32G
#PBS -l walltime=48:00:00

cd $PBS_O_WORKDIR

module load pytorch/1.5.0-conda3.7-GPU

source PaperDescriptiveEnvHPF/bin/activate

python3 Model_MultilayerPerceptron_Validation_CoronalSagittalSampling.py --labels region67 --nsamples 1 --nunits 200 --nlayers 3 --dropout 0 --L2 1e-6 --confusionmatrix true

#mv Data/ValidationExperiments/MLP_Validation_CoronalSagittalSampling_Region67_ImputeKNN_temp.csv Data/ValidationExperiments/MLP_Validation_CoronalSagittalSampling_Region67_ImputeKNN_temp_05.csv