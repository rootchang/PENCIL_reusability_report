import sys
from pencil import *
from matplotlib import pyplot as plt
import scanpy as sc
import os
os.environ['CUDA_VISIBLE_DEVICES'] = '0,1' #select a gpu id, e.g. '0', '0,1', use all gpus by setting to ''. otherwise set to '-1'.
import warnings
warnings.filterwarnings('ignore') #to ignore warnings on the output.
import time


start_time  = time.time()
h5ad_fn = sys.argv[1] # '../02.Input/GSE200996/seu_Tissue_Tcell.h5ad'
phenotype = sys.argv[2] # 'RECIST_response' # 'Volumetric_response'
data_name = sys.argv[3] # 'GSE200996'
mode = sys.argv[4] # 'multi-classification'  'regression'
use_HVG = int(sys.argv[5]) # 1: use high-variable-genes with scaled data. 0: use all genes with log1p expression
interrupt = int(sys.argv[6]) # 1: Interrupt to check the raw phenotypic labels. 0: no interruption.


print('************** Step 0: loading data ...')
# Note that adata.X stores the scaled data of HVGs
#           adata.raw.X stores the log1p normalized expression of all genes
adata = sc.read_h5ad(h5ad_fn)
print(adata)


print('************** Step 1: preparing PENCIL input parameters ...')
labels_raw = adata.obs[phenotype]
print('Raw phenotypic labels: ', set(labels_raw))
if interrupt:
    raise Exception('Interrupt to check the raw phenotypic labels')
# Specify the phenotypic values to be removed
values_to_remove = ['not measurable']
# Subset the AnnData object
adata = adata[~adata.obs[phenotype].isin(values_to_remove)]#.copy()
labels_raw = adata.obs[phenotype]
print('Clean labels: ', set(labels_raw))
labels_raw = pd.Categorical(labels_raw)
if use_HVG:
    data, labels = adata.X.copy(), labels_raw.codes
else:
    data = adata.raw.X.copy().toarray()
    labels = labels_raw.codes
class_names = list(labels_raw.categories)
emd = adata.obsm['X_umap']
print('data.shape: ', data.shape)
print('labels.shape: ', labels.shape)
print('class_names: ', class_names)
print('emd.shape: ', emd.shape)


print('************** Step 2: running PENCIL ...')
# class_weights = [2.0, 1.0]
pencil = Pencil(mode, select_genes=True, seed=1234, data_name=data_name, expr_id=phenotype, mlflow_record=False,
                model_types=['linear', 'non-linear'])
pred, confidence = pencil.fit_transform(
    data, labels,
    test=True,
    # c=None, # default is None
    shuffle_rate=1/4, # 1/3, # default is 1/3
    lambda_L1=1e-5, # 1e-4, # default is 1e-5
    lambda_L2=1e-3, # default is 1e-3
    lr=0.01, # default is 0.01
    # epochs=500, # default is 500
    # pre_train_epochs=500, # default = 500 for classification
    class_weights=None, # class_weights, # default is None
    class_names=class_names, # default is None
    emd=emd,
    plot_show=True
    )
w = pencil.gene_weights(plot=True)
plt.close()
# pred_new, confidence_new = pencil.transform(data_new) # apply trained PENCIL on new data

print('All done! Time used: %.2f s'%(time.time() - start_time))