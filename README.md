# PENCIL_usage
This is a reproducibility report of PENCIL (https://doi.org/10.1038/s42256-023-00656-y)

This is an ongoing project. To do list:

```
1. UMAP of patients, R/NR binarized patient labels, PENCIL annotation (R/NR/Rej.)
2. Marker gene analysis
3. Enriched pathway analysis
4. Known signature enrichment (e.g. CD8T exh., cyt., TRM., ...)
4. Use the learned biological information to predict on new datasets
5. Benchmark PENCIL with other relevant tools

```


## 1. How to run PENCIL on Biowulf

(1.1) One should use GPU to accelarate the speed of PENCIL, otherwise it would be extremely slow

```{bash}
sinteractive --mem=128g --cpus-per-task=32 --gres=gpu:k80:2 --time=24:00:00
```

(1.2) I use python/3.8 on the Biowulf to run PENCIL. I did not check if python/3.9 also works or not.

```{bash}
module load python/3.8
```

(1.3) I encountered error during running PENCIL initially, then I changed the version of the numpy package in python to 1.21, it worked out.

```{bash}
pip install numpy==1.21
```



## 2. Analysis for individual datasets

(2.1) I have made a pipeline `02.Run_PENCIL.py` for running PENCIL with the .h5ad input files. The basic usage is as follows:

```{bash}
python 02.Run_PENCIL.py p1 p2 p3 p4 p5 p6
```

Table 1. Explanation for the parameters used by `02.Run_PENCIL.py`.

| Parameter | Description | Default |
|----------|----------|----------|
| p1   | Input .h5ad file | required, no default |
| p2   | phenotype name, stored in seurat object meta.data (r) or adata.obs (python) | required, no default |
| p3   | experiment id, PENCIL will make a new folder in the results directory use this name | required, no default |
| p4   | PENCIL mode, can be either 'multi-classification' or 'regression' | required, no default |
| p5   | To use HVGs (1) or all genes (0) for PENCIL prediction | required, no default |
| p6   | To interrupt the program (for one to adjust the variable 'values_to_remove' in line 32 in 02.Run_PENCIL.py to clean the phenotypic values) before running PENCIL (1) or not (0) | required, no default |

(2.2) I have prepared .h5ad files as input for PENCIL for all the datasets. The location of these datasets is `/data/Lab_ruppin/tiangen/CancerProject/02.PENCIL/02.Input`.

(2.2.1) To run PENCIL for the GSE120575 melanoma dataset using tissue CD8T/CD4T/B cells that presented in the original PENCIL paper, run the following code:

```{bash}

# tissue CD8T All cells
python 02.Run_PENCIL_csv.py /data/Lab_ruppin/tiangen/CancerProject/02.PENCIL/03.Results/GSE120575/seu_Tissue_CD8T.csv /data/Lab_ruppin/tiangen/CancerProject/02.PENCIL/03.Results/GSE120575/cell_label_info.csv /data/Lab_ruppin/tiangen/CancerProject/02.PENCIL/03.Results/GSE120575/seu_Tissue_CD8T_embedding_umap.csv GSE120575_Tissue_CD8T ResponseInfo multi-classification 2 /data/Lab_ruppin/tiangen/CancerProject/02.PENCIL/03.Results/GSE120575/seu_Tissue_CD8T_Test.csv

# tissue CD8T Train cells
python 02.Run_PENCIL_csv.py /data/Lab_ruppin/tiangen/CancerProject/02.PENCIL/03.Results/GSE120575/seu_Tissue_CD8T_Train.csv /data/Lab_ruppin/tiangen/CancerProject/02.PENCIL/03.Results/GSE120575/cell_label_info_Train.csv /data/Lab_ruppin/tiangen/CancerProject/02.PENCIL/03.Results/GSE120575/seu_Tissue_CD8T_Train_embedding_umap.csv GSE120575_Tissue_CD8T_Train ResponseInfo multi-classification 2 /data/Lab_ruppin/tiangen/CancerProject/02.PENCIL/03.Results/GSE120575/seu_Tissue_CD8T_Test.csv

```

To run PENCIL on the wrongly processed GSE120575 melanoma dataset by the PENCIL authors, using the following code:

```

# tissue CD8T cells
python 02.Run_PENCIL.py /data/Lab_ruppin/tiangen/CancerProject/02.PENCIL/02.Input/GSE120575/seu_Tissue_CD8T_PENCIL_wrong.h5ad ResponseInfo GSE120575_Tissue_CD8T_PENCIL_wrong multi-classification 1 0

# OR
python 02.Run_PENCIL_csv.py /data/Lab_ruppin/tiangen/CancerProject/02.PENCIL/02.Input/GSE120575/exp_data_mvg2000.csv /data/Lab_ruppin/tiangen/CancerProject/02.PENCIL/02.Input/GSE120575/label_info.csv /data/Lab_ruppin/tiangen/CancerProject/02.PENCIL/02.Input/GSE120575/embedding-umap.csv GSE120575_Tissue_CD8T ResponseInfo multi-classification

exp_data_mvg2000.csv
exp_data_mvg2000_unscale.csv
exp_data_mvg2000_log2TPM.csv

exp_data_mvg2000_log2TPM_self.csv
exp_data_mvg2000_scale_self.csv

```



(2.2.2) To run PENCIL for the GSE200996 HNSCC dataset using tissue/PBMC CD8T/CD4T/B cells, run the following code:

```{bash}
# tissue CD8T cells
python 02.Run_PENCIL.py /data/Lab_ruppin/tiangen/CancerProject/02.PENCIL/02.Input/GSE200996/seu_Tissue_CD8T.h5ad Any_response GSE200996_Tissue_CD8T multi-classification 1 0

# tissue CD4T cells
python 02.Run_PENCIL.py /data/Lab_ruppin/tiangen/CancerProject/02.PENCIL/02.Input/GSE200996/seu_Tissue_CD4T.h5ad Any_response GSE200996_Tissue_CD4T multi-classification 1 0

# tissue B cells
python 02.Run_PENCIL.py /data/Lab_ruppin/tiangen/CancerProject/02.PENCIL/02.Input/GSE200996/seu_Tissue_B.h5ad Any_response GSE200996_Tissue_B multi-classification 1 0

# PBMC CD8T cells
python 02.Run_PENCIL.py /data/Lab_ruppin/tiangen/CancerProject/02.PENCIL/02.Input/GSE200996/seu_PBMC_CD8T.h5ad Any_response GSE200996_PBMC_CD8T multi-classification 1 0

# PBMC CD4T cells
python 02.Run_PENCIL.py /data/Lab_ruppin/tiangen/CancerProject/02.PENCIL/02.Input/GSE200996/seu_PBMC_CD4T.h5ad Any_response GSE200996_PBMC_CD4T multi-classification 1 0

# PBMC B cells
python 02.Run_PENCIL.py /data/Lab_ruppin/tiangen/CancerProject/02.PENCIL/02.Input/GSE200996/seu_PBMC_B.h5ad Any_response GSE200996_PBMC_B multi-classification 1 0

```

(2.2.3) To run PENCIL for the GSE200996 HNSCC dataset using PBMC CD8T cells, run the following code:

```{bash}
python 02.Run_PENCIL.py /data/Lab_ruppin/tiangen/CancerProject/02.PENCIL/02.Input/GSE200996/seu_PBMC_CD8T.h5ad Any_response GSE200996_PBMC_CD8T multi-classification 1 0
```
