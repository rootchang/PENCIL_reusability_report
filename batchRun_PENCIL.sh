#!/bin/bash

# Create arrays with your changing parameters
classweight=(1.5 1.6 1.7 1.8 1.9 1.95 1.96 1.97 1.98 1.99 2.01 2.02 2.03 2.04 2.05 2.06 2.07 2.08 2.09 2.1 2.2 2.3 2.4 2.5)

for idx in ${!classweight[@]}; do
  python 02.Run_PENCIL_csv.py /data/Lab_ruppin/tiangen/CancerProject/02.PENCIL/02.Input/GSE120575/exp_data_mvg2000.csv /data/Lab_ruppin/tiangen/CancerProject/02.PENCIL/02.Input/GSE120575/label_info.csv /data/Lab_ruppin/tiangen/CancerProject/02.PENCIL/02.Input/GSE120575/embedding-umap.csv GSE120575_Tissue_CD8T ResponseInfo multi-classification ${classweight[idx]}
done
