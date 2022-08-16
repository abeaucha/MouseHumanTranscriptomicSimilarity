#!/bin/bash

source activate_venv.sh

datadir=data/
outdir=data/MLP_validation/

method='resampling'
# method='naive'

nfolds=5
nlayers=3
nunits='200 500 1000'
L2='0 1e-6 1e-3'
dropout=0
nepochs=200
learningrate='1e-1 1e-2 1e-3 1e-4 1e-5'
optimizer='SGD AdamW' 
seed=1

if [ $method = 'resampling' ];
then
    python3 mlp_validation_resampling.py \
        --datadir ${datadir} \
        --outdir ${outdir} \
        --labels region67 \
        --nsamples $nfolds \
        --nunits $nunits \
        --nlayers $nlayers \
        --L2 $L2 \
        --dropout $dropout \
        --nepochs $nepochs \
        --learningrate $learningrate \
        --optimizer $optimizer \
        --seed $seed
elif [ $method = 'naive' ];
then 
    python3 mlp_validation_naive.py \
    --datadir $datadir \
	 --outdir $outdir \
	 --labels region67 \
     --nfolds $nfolds \
	 --nunits $nunits \
     --nlayers $nlayers \
	 --L2 $L2 \
     --dropout $dropout \
	 --nepochs $nepochs \
	 --learningrate $learningrate \
     --optimizer $optimizer \
     --seed $seed
fi

#Original parameters
# nunits=200
# nlayers=3
# dropout=0
# L2=1e-6
# nepochs=200
# learningrate=1e-5
# optimizer='AdamW'

#Updated parameters
# nunits=200
# nlayers=3
# dropout=0
# L2=0
# nepochs=50
# learningrate=0.01
# optimizer='SGD'     
     
deactivate