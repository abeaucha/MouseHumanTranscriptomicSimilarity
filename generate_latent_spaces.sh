#!/bin/bash -l

# ----------------------------------------------------------------------------
# generate_latent_spaces
# Author: Antoine Beauchamp
#
# This script generates 500 gene expression latent spaces by repeatedly 
# training the multi-layer perceptron neural network and using the network
# architecture to transform the input space into the latent space.

source activate_venv.sh

datadir=data/
# datadir=data/isocortex/

outdir=data/MLP_outcomes/
# outdir=data/MLP_outcomes_isocortex/

if [ ! -d "$outdir" ]; then
	mkdir "$outdir"
fi

nunits=200
L2=0.0

for i in {1..500}
do
 echo "Iteration $i"

 python3 train_multilayer_perceptron.py \
	 --datadir $datadir \
	 --outdir $outdir \
	 --labels region67 \
	 --mousedata region67 \
	 --humandata region88 \
	 --nunits $nunits \
	 --L2 $L2 \
	 --nepochs 25 \
	 --totalsteps 50 \
	 --learningrate 0.01 \
	 --confusionmatrix false \
	 --voxeltransform true \
	 --seed $i

#  mv ${outdir}MLP_Region67_Layers3_Units${nunits}_L2${L2}_IntegratedGradients.csv \
#      ${outdir}MLP_Region67_Layers3_Units${nunits}_L2${L2}_IntegratedGradients_$i.csv

 mv ${outdir}MLP_Region67_Layers3_Units${nunits}_L2${L2}_HumanProb_Region88.csv \
	 ${outdir}MLP_Region67_Layers3_Units${nunits}_L2${L2}_HumanProb_Region88_$i.csv

 mv ${outdir}MLP_Region67_Layers3_Units${nunits}_L2${L2}_HumanTx_Region88.csv \
	 ${outdir}MLP_Region67_Layers3_Units${nunits}_L2${L2}_HumanTx_Region88_$i.csv

 mv ${outdir}MLP_Region67_Layers3_Units${nunits}_L2${L2}_MouseProb_Region67.csv \
	 ${outdir}MLP_Region67_Layers3_Units${nunits}_L2${L2}_MouseProb_Region67_$i.csv

 mv ${outdir}MLP_Region67_Layers3_Units${nunits}_L2${L2}_MouseTx_Region67.csv \
	 ${outdir}MLP_Region67_Layers3_Units${nunits}_L2${L2}_MouseTx_Region67_$i.csv

mv ${outdir}MLP_Region67_Layers3_Units200_L2${L2}_MouseVoxelTx.csv \
	${outdir}MLP_Region67_Layers3_Units200_L2${L2}_MouseVoxelTx_$i.csv

mv ${outdir}MLP_Region67_Layers3_Units200_L2${L2}_HumanVoxelTx.csv \
	${outdir}MLP_Region67_Layers3_Units200_L2${L2}_HumanVoxelTx_$i.csv

done

deactivate
