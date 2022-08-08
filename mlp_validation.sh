#!/bin/bash

source activate_venv

datadir=data/
outdir=data/MLP_validation/

nunits=200
L2=1e-6
nepochs=200
learningrate=1e-5
optimizer='AdamW'

python3 mlp_validation_naive.py \
	 --datadir $datadir \
	 --outdir $outdir \
	 --labels region67 \
	 --nunits $nunits \
	 --L2 $L2 \
	 --nepochs $nepochs \
	 --learningrate $learningrate \
     --optimizer $optimizer \
	 --seed 1
     
mv ${outdir}MLP_validation_naive_region67.csv \
${outdir}MLP_validation_naive_region67_original.csv
     
python3 mlp_validation_resampling.py \
    --datadir ${datadir} \
    --outdir ${outdir} \
    --labels region67 \
    --nsamples 10 \
    --nunits $nunits \
    --nlayers 3 \
    --dropout 0 \
    --L2 $L2 \
    --nepochs $nepochs \
    --optimizer $optimizer \
    --learningrate $learningrate \
    --confusionmatrix false \
    --seed 1
    
mv ${outdir}MLP_Validation_CoronalSagittalSampling_Region67.csv \
${outdir}MLP_Validation_CoronalSagittalSampling_Region67_original.csv
     
nunits=200
L2=0
nepochs=50
learningrate=0.01
optimizer='SGD'     
     
python3 mlp_validation_naive.py \
	 --datadir $datadir \
	 --outdir $outdir \
	 --labels region67 \
	 --nunits $nunits \
	 --L2 $L2 \
	 --nepochs $nepochs \
	 --learningrate $learningrate \
     --optimizer $optimizer \
	 --seed 1
     
mv ${outdir}MLP_validation_naive_region67.csv \
${outdir}MLP_validation_naive_region67_updated.csv
     
python3 mlp_validation_resampling.py \
    --datadir ${datadir} \
    --outdir ${outdir} \
    --labels region67 \
    --nsamples 10 \
    --nunits $nunits \
    --nlayers 3 \
    --dropout 0 \
    --L2 $L2 \
    --nepochs $nepochs \
    --optimizer $optimizer \
    --learningrate $learningrate \
    --confusionmatrix false \
    --seed 1
        
mv ${outdir}MLP_Validation_CoronalSagittalSampling_Region67.csv \
${outdir}MLP_Validation_CoronalSagittalSampling_Region67_updated.csv

deactivate