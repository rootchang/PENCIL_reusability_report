
#########################################################################################################
# This version of PENCIL takes .csv expression file and .csv cell label file as input
#########################################################################################################

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
exp_fn = sys.argv[1] # '../02.Input/GSE200996/seu_Tissue_CD8T_PENCIL_wrong.csv'
anno_fn = sys.argv[2] # '../02.Input/GSE200996/cellLabel_ICBresponse.csv'
embedding_fn = sys.argv[3] # '../02.Input/GSE200996/seu_Tissue_CD8T_embedding_umap.csv'
data_name = sys.argv[4] # 'GSE200996' 'GSE120575_Tissue_CD8T'
phenotype = sys.argv[5] # 'ResponseInfo'
mode = sys.argv[6] # 'multi-classification'  'regression'
class_weights = [float(eval(x)) for x in sys.argv[7].split(',')]
exp_newdata_fn_list = sys.argv[8].split(',') # new test datasets, separates by ','

print('************** Step 0: loading data ...')
exp_df=pd.read_csv(exp_fn, sep=',',index_col=0)
data=exp_df.values.T
anno_df = pd.read_csv(anno_fn, sep=',', index_col=0, dtype=str)
labels_raw = anno_df.values.flatten()

print('************** Step 1: preparing PENCIL input parameters ...')
print('Cell labels: ', set(labels_raw))
if mode == 'multi-classification':
    labels_raw = pd.Categorical(labels_raw)
    class_names = list(labels_raw.categories)
    labels = labels_raw.codes
else:
    labels = np.array(labels_raw.values, dtype=float)

#class_weights = [1.0, class_weight] # [1.0, 2.0] # None

print('Expression data (first 5*5):\n', data[0:5,0:5])
if mode == 'multi-classification':
    print('class_names: ', class_names)
print('Unique labels and occurrence: ', sorted(list(Counter(labels).items()))) # sorted(list(Counter(labels).items()))
print('class_weights: ', class_weights)
print('labels.shape: ', labels.shape)


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
        anno_file=anno_fn,
        embedding_file=embedding_fn,
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
            anno_file=anno_fn,
            embedding_file=embedding_fn,
            plot_show=True
        )

w = pencil.gene_weights(plot=True)
plt.close()

print('************** Step 3: predicting on new data ...')
for i in range(len(exp_newdata_fn_list)):
    print('    Working on test set %d ...'%(i+1))
    exp_newdata_fn = exp_newdata_fn_list[i]
    exp_new_df=pd.read_csv(exp_newdata_fn, sep=',',index_col=0)
    data_new=exp_new_df.values.T
    pred_new, confidence_new = pencil.transform(data_new) # apply trained PENCIL on new data
    save_fn = './results/'+data_name+'/py/'+phenotype+'/predicted_labels_test'+str(i+1)+'.csv'
    df = pd.DataFrame({'predicted_label': pred_new, 'confidence': confidence_new})
    df.loc[df['confidence'] < 0, 'predicted_label'] = 'Rejected'
    df.to_csv(save_fn, index=False)

print('All done! Time used: %.2f s'%(time.time() - start_time))