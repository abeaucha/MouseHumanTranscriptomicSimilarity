#!/bin/bash

source activate_venv.sh

datadir=data/
outdir=data/MLP_validation/

method='resampling'
# method='naive'

nsamples=5
nlayers=3
nunits='200 500 1000'
L2='0 1e-3 1e-6'
dropout=0
nepochs=200
learningrate='1e-1 1e-2 1e-3 1e-4 1e-5'
optimizer='SGD AdamW' 
seed=1

python3 mlp_validation_resampling.py \
    --datadir ${datadir} \
    --outdir ${outdir} \
    --labels region67 \
    --nsamples $nsamples \
    --nunits $nunits \
    --nlayers $nlayers \
    --dropout $dropout \
    --L2 $L2 \
    --nepochs $nepochs \
    --optimizer $optimizer \
    --learningrate $learningrate \
    --seed $seed
    
# python3 mlp_validation_naive.py \
# 	 --datadir $datadir \
# 	 --outdir $outdir \
# 	 --labels region67 \
# 	 --nunits $nunits \
# 	 --L2 $L2 \
# 	 --nepochs $nepochs \
# 	 --learningrate $learningrate \
#      --optimizer $optimizer \
# 	 --seed $seed
     
# nunits=200
# L2=0
# nepochs=50
# learningrate=0.01
# optimizer='SGD'     
     
# python3 mlp_validation_resampling.py \
#     --datadir ${datadir} \
#     --outdir ${outdir} \
#     --labels region67 \
#     --nsamples 10 \
#     --nunits $nunits \
#     --nlayers 3 \
#     --dropout 0 \
#     --L2 $L2 \
#     --nepochs $nepochs \
#     --optimizer $optimizer \
#     --learningrate $learningrate \
#     --confusionmatrix false \
#     --seed 1
        
deactivate