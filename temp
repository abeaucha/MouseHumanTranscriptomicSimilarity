#!/bin/bash -l

for i in {251..500}
do
 echo "Iteration $i"
 python3 Model_MultilayerPerceptron.py --labels region67 --mousedata region67 --humandata region88 --nunits 200 --confusionmatrix true --voxeltransform true
 mv Data/MLP_ConfusionMatrix_Training_Region67_Layers3_Units200_L21e-06.csv Data/LatentSpaceStability/MLP_ConfusionMatrix_Training_Region67_Layers3_Units200_L21e-06_$i.csv
 mv Data/MLP_Region67_Layers3_Units200_L21e-06_HumanProb_Region88.csv Data/LatentSpaceStability/MLP_Region67_Layers3_Units200_L21e-06_HumanProb_Region88_$i.csv 
 mv Data/MLP_Region67_Layers3_Units200_L21e-06_HumanTx_Region88.csv Data/LatentSpaceStability/MLP_Region67_Layers3_Units200_L21e-06_HumanTx_Region88_$i.csv 
 mv Data/MLP_Region67_Layers3_Units200_L21e-06_MouseProb_Region67.csv Data/LatentSpaceStability/MLP_Region67_Layers3_Units200_L21e-06_MouseProb_Region67_$i.csv
 mv Data/MLP_Region67_Layers3_Units200_L21e-06_MouseTx_Region67.csv Data/LatentSpaceStability/MLP_Region67_Layers3_Units200_L21e-06_MouseTx_Region67_$i.csv
mv Data/MLP_Region67_Layers3_Units200_L21e-06_MouseVoxelTx.csv Data/LatentSpaceStability/MLP_Region67_Layers3_Units200_L21e-06_MouseVoxelTx_$i.csv
mv Data/MLP_Region67_Layers3_Units200_L21e-06_HumanVoxelTx.csv Data/LatentSpaceStability/MLP_Region67_Layers3_Units200_L21e-06_HumanVoxelTx_$i.csv
done

