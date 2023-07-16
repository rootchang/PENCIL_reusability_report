import sys
from pencil import *
from matplotlib import pyplot as plt
import scanpy as sc
import os
os.environ['CUDA_VISIBLE_DEVICES'] = '0,1' #select a gpu id, e.g. '0', '0,1', use all gpus by setting to ''. otherwise set to '-1'.
import warnings
warnings.filterwarnings('ignore') #to ignore warnings on the output.
import time
from scipy.sparse import issparse
from collections import Counter

start_time  = time.time()
h5ad_fn = sys.argv[1] # '../02.Input/GSE200996/seu_Tissue_Tcell.h5ad'
phenotype = sys.argv[2] # 'RECIST_response' # 'Volumetric_response' # "Any_response"
data_name = sys.argv[3] # 'GSE200996'
mode = sys.argv[4] # 'multi-classification'  'regression'
use_HVG = int(sys.argv[5]) # 1: use high-variable-genes with scaled data. 0: use all genes with log1p expression (?)
use_scaledData = int(sys.argv[6]) # 1: use scaled data (for HVG only). 0: use normalized but unscaled data (for either HVGs or all genes).
interrupt = int(sys.argv[7]) # 1: Interrupt to check the feasibility of data (e.g., remove unwanted phenotypic labels and determine the class weights). 0: no interruption.


print('************** Step 0: loading data ...')
# Note that adata.X stores the scaled data of HVGs
#           adata.raw.X stores the log1p normalized expression of all genes
adata = sc.read_h5ad(h5ad_fn)
print(adata)

print('************** Step 1: preparing PENCIL input parameters ...')
# Specify the phenotypic values to be removed
values_to_remove = ['NA'] # ['not measurable']
# Subset the AnnData object
adata = adata[~adata.obs[phenotype].isin(values_to_remove)]#.copy()
labels_raw = adata.obs[phenotype]
if mode == 'multi-classification':
    labels_raw = pd.Categorical(labels_raw)
    class_names = list(labels_raw.categories)
    labels = labels_raw.codes
    #label_dict = dict(zip(class_names, range(len(class_names))))
    #labels = list(map(lambda s: label_dict[s], labels_raw))
    #labels = np.array(labels)
else:
    labels = np.array(labels_raw.values, dtype=float)
    #print('class_labels: ', labels)
if use_HVG:
    if use_scaledData:
        if issparse(adata.X):
            data = adata.X.copy().todense()
        else:
            data = adata.X.copy()
    else:
        #high_var_genes = adata.raw[:, adata.var['vst.variable']].var_names
        if issparse(adata.raw.X):
            data = adata.raw[:, adata.var['vst.variable']].X.copy().todense()
        else:
            data = adata.raw[:, adata.var['vst.variable']].X.copy()
else:
    if issparse(adata.raw.X):
        data = adata.raw.X.copy().todense()
    else:
        data = adata.raw.X.copy()
class_weights = [1.0, 2.0] # [1.0, 2.0] # None
emd = adata.obsm['X_umap']
print('data.shape: ', data.shape)
if use_HVG:
    print('Genes (first 5): ', adata.var['vst.variable'].index.tolist()[0:5])
else:
    print('Genes (first 5): ', list(adata.var_names[0:5]))
print('Cells (first 5): ', list(adata.obs_names[0:5]))
print('Expression data (first 5*5):\n', data[0:5,0:5])
if mode == 'multi-classification':
    print('class_names: ', class_names)
print('Unique labels and occurrence: ', sorted(list(Counter(labels).items()))) # sorted(list(Counter(labels).items()))
print('class_weights: ', class_weights)
print('labels.shape: ', labels.shape)
print('emd.shape: ', emd.shape)

if interrupt:
    raise Exception('Interrupt to check the feasibility of input data.')

print('************** Step 2: running PENCIL ...')
if mode == 'multi-classification':
    pencil = Pencil(mode, select_genes=True, seed=1234, data_name=data_name, expr_id=phenotype, model_types=['linear', 'non-linear'])
    pred, confidence = pencil.fit_transform(
        data, labels,
        test=True,
        # c=None, # default is None
        shuffle_rate=1/3, # 1/4, # default is 1/3
        lambda_L1=1e-4, # 1e-5, # default is 1e-5
        lambda_L2=1e-3, # 1e-3, default is 1e-3
        lr=0.01, # default is 0.01
        epochs=500, # default is 500
        # pre_train_epochs=500, # default = 500 for classification
        class_weights=class_weights, # class_weights, # default is None
        class_names=class_names, # default is None
        emd=emd,
        plot_show=True
        )
else:
    pencil = Pencil(mode, select_genes=True, seed=1234, data_name=data_name, expr_id=phenotype, mlflow_record=True,
                    dropouts=[0.0, 0.0])  # `select_genes` can also be set to False, if True, pencil will output a weight vector for all 10 PCs.
    with mlflow.start_run():
        pred, confidence = pencil.fit_transform(
            data, labels,
            test=True,
            shuffle_rate=1 / 5,
            lambda_L1=1e-5,
            lambda_L2=0.0,
            lr=0.01,
            epochs=2000,
            class_weights=None,
            emd=emd,
            plot_show=True
        )

w = pencil.gene_weights(plot=True)
plt.close()
# pred_new, confidence_new = pencil.transform(data_new) # apply trained PENCIL on new data

print('All done! Time used: %.2f s'%(time.time() - start_time))