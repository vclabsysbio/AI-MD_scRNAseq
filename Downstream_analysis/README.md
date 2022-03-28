# Downstream analysis

A feature-barcode matrix derived from Cellranger will be used for scRNA-seq analysis and visualization. Scanpy is a toolkit for analysing single-cell gene expression data (Wolf *et al.*, 2018). There are several tutorials for analyzing scRNA-seq data, including clustering, trajectory inference, and integrating datasets, as well as a variety of data types, including PBMC and neuron. Scanpy's analysis pipeline, algorithm, and methods are similar to Seurat, an R-based package (Stuart *et al.*, 2019).
    
_All commands were modified from [Scanpy](https://scanpy-tutorials.readthedocs.io/en/latest/pbmc3k.html) and [RAPIDS scRNA-seq](https://github.com/clara-parabricks/rapids-single-cell-examples)_

- [Tools](#Tools)
- [Data analysis pipeline](#Data-analysis-pipeline)
- [Commands](#Commands)
- [CPU and GPU performance testing](#CPU-and-GPU-performance-testing)



## Tools
**Python** 
- version 3.7

**Packages**
- scanpy==1.8.2
- anndata==0.8.0
- umap==0.5.1
- numpy==1.21.2
- scipy==1.6.0
- pandas==1.3.3
- scikit-learn==0.24.2
- statsmodels==0.13.0
- python-igraph==0.9.9
- pynndescent==0.5.4-0.5.6
- igraph==0.9.7
- leidenalg==0.8.9

## Data analysis pipeline
- **Download GitHub** (GPU only)
- **Install packages** (ICT-HPC server only)
- **Import requirements**
- **Input data & load data**
- **Prepare Data** (GPU only)
- **Preprocessing**
- **Normalization & Scaling the data**
- **Select Most Variable Genes**
- **Regress out confounding factors (number of counts, mitochondrial gene expression)**
- **Perform linear dimensional reduction**
- **Clustering**
- **Run non-linear dimensional reduction (UMAP)**
- **Finding marker genes  & Differential expression analysis**
- **Cell type identification**
- **Visualization (dot plot & violin plot)** (optional)
- **Save file**


## Commands
### Download GitHub (GPU only)
```
git clone https://github.com/clara-parabricks/rapids-single-cell-examples.git
```

### Install packages
``` {python}
import sys
sys.executable
!{sys.executable} -m pip install wget
!{sys.executable} -m pip install numpy
!{sys.executable} -m pip install scipy
!{sys.executable} -m pip install pandas
!{sys.executable} -m pip install scikit-learn
!{sys.executable} -m pip install statsmodels
!{sys.executable} -m pip install leidenalg
!{sys.executable} -m pip install anndata
!{sys.executable} -m pip install scanpy
!{sys.executable} -m pip install scirpy
```

### Import requirements
**CPU version**
``` {python}
import time
import os, wget

import numpy as np
import pandas as pd
import scanpy as sc

import warnings
warnings.filterwarnings('ignore', 'Expected ')
warnings.simplefilter('ignore')
```

**GPU version**
``` {python}
import scanpy as sc
import anndata
#
import numpy as np
import pandas as pd
#
import time
import os, wget

import cudf
import cupy as cp

from cuml.decomposition import PCA
from cuml.manifold import TSNE
from cuml.cluster import KMeans
from cuml.preprocessing import StandardScaler

import rapids_scanpy_funcs

import warnings
warnings.filterwarnings('ignore', 'Expected ')
warnings.simplefilter('ignore')

import rmm

rmm.reinitialize(
    managed_memory=True, # Allows oversubscription
    pool_allocator=False, # default is False
    devices=0, # GPU device IDs to register. By default registers only GPU 0.
)

cp.cuda.set_allocator(rmm.rmm_cupy_allocator)
```

### Input data & load data
**CPU and GPU version**
``` {python}
adata = sc.read_10x_mtx(
    './data/',  # the directory with the `.mtx` file
    var_names='gene_symbols',
    cache=True) 
adata.var_names_make_unique()
```

### Prepare Data (GPU only)
``` {python}
genes = cudf.Series(adata.var_names)
sparse_gpu_array = cp.sparse.csr_matrix(adata.X)
```

### Preprocessing
**CPU version**
``` {python}
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)

adata.var['mt'] = adata.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

adata = adata[adata.obs.n_genes_by_counts < 3000, :]
adata = adata[adata.obs.pct_counts_mt < 20, :]
```

**GPU version**
``` {python}
sparse_gpu_array = rapids_scanpy_funcs.filter_cells(sparse_gpu_array, min_genes=200, max_genes=3000)
sparse_gpu_array, genes = rapids_scanpy_funcs.filter_genes(sparse_gpu_array, genes, min_cells=3)
```


### Normalization & Scaling the data
**CPU version**
``` {python}
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
```

**GPU version**
``` {python}
sparse_gpu_array = rapids_scanpy_funcs.normalize_total(sparse_gpu_array, target_sum=1e4)
sparse_gpu_array = sparse_gpu_array.log1p()
```

### Select Most Variable Genes
**CPU version**
``` {python}
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)

adata = adata[:, adata.var.highly_variable]
```

**GPU version**
``` {python}
markers = ['CD8A', 'CD8B', 'IL7R', 'CCR7', 'CST7', 'KLRB1', 
                'S100A4', 'CD14', 'LYZ', 'LGALS3', 'S100A8', 
                'MS4A1','CD79A', #, 'GNLY', 'NKG7'
                'FCGR3A', 'MS4A7', 'PPBP', 'FCER1A', 'CST3'] # Marker genes for visualization

tmp_norm = sparse_gpu_array.tocsc()
marker_genes_raw = {
    ("%s_raw" % marker): tmp_norm[:, genes[genes == marker].index[0]].todense().ravel()
    for marker in markers
}

del tmp_norm

hvg = rapids_scanpy_funcs.highly_variable_genes(sparse_gpu_array, genes, n_top_genes=5000)

sparse_gpu_array = sparse_gpu_array[:, hvg]
genes = genes[hvg].reset_index(drop=True)
sparse_gpu_array.shape
```

### Regress out confounding factors (number of counts, mitochondrial gene expression)
**CPU version**
``` {python}
sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
```

**GPU version**
``` {python}
mito_genes = genes.str.startswith("MT-")
n_counts = sparse_gpu_array.sum(axis=1)
percent_mito = (sparse_gpu_array[:,mito_genes].sum(axis=1) / n_counts).ravel()

n_counts = cp.array(n_counts).ravel()
percent_mito = cp.array(percent_mito).ravel()

sparse_gpu_array = rapids_scanpy_funcs.regress_out(sparse_gpu_array, n_counts, percent_mito)
```

### Perform linear dimensional reduction
**CPU version**
``` {python}
sc.tl.pca(adata, svd_solver='arpack')
```

**GPU version**
``` {python}
adata = anndata.AnnData(sparse_gpu_array.get())
adata.var_names = genes.to_pandas()

for name, data in marker_genes_raw.items():
    adata.obs[name] = data.get()

adata.obsm["X_pca"] = PCA(n_components=50, output_type="numpy").fit_transform(adata.X)
```

### Clustering
**CPU version**
``` {python}
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=15)
sc.tl.leiden(adata, resolution=0.3)
```

**GPU version**
``` {python}
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=15, method='rapids')
adata.obs['leiden'] = rapids_scanpy_funcs.leiden(adata, resolution=0.3)
```

### Run non-linear dimensional reduction (UMAP)
**CPU version**
``` {python}
sc.tl.umap(adata)
```

**GPU version**
``` {python}
sc.tl.umap(adata, min_dist=0.5, spread=1.0, method='rapids')
```

### Finding marker genes  & Differential expression analysis
**CPU version**
``` {python}
sc.tl.rank_genes_groups(adata, groupby="leiden", n_genes=20, groups='all', reference='rest', method='wilcoxon')
```

**GPU version**
``` {python}
cluster_labels = cudf.Series.from_categorical(adata.obs["leiden"].cat)
genes = cudf.Series(genes)

scores, names, reference = rapids_scanpy_funcs.rank_genes_groups(
    sparse_gpu_array, 
    cluster_labels, 
    genes, 
    n_genes=20, groups='all', reference='rest')
```

### Cell type identification
**CPU and GPU version**
``` {python}
new_cluster_names = [
    'CD8+ T', 'Naive CD4+ T',
    'NK', 'Memory CD4+ T',
    'CD14+ Mono', 'B',
    'FCGR3A+ Mono', 'Platelets',
    'DC']
adata.obs['leiden'] =adata.obs['leiden'].cat.rename_categories(new_cluster_names)
sc.pl.umap(adata, color='leiden', legend_loc='on data', title='', frameon=False, save='.pdf')
```
![Alt text](https://github.com/vclabsysbio/AI-MD_scRNAseq/blob/main/Downstream_analysis/Figures/UMAP_plot_CPU_vs_GPU.jpg?raw=true "UMAP")

### Visualization (dot plot & violin plot) (optional)
**CPU and GPU version**
``` {python}
#Dot plot
sc.pl.dotplot(adata, marker_genes, groupby='leiden')
#Violin plot
sc.pl.stacked_violin(adata, marker_genes, groupby='leiden', rotation=90)
```

### Save file
**CPU and GPU version**
``` {python}
results_file = './write/filename.h5ad'
```
``` {python}
adata.write(results_file)
```

## CPU and GPU performance testing
### Comparison


| Analysis steps                | CPU of 246 server        | CPU of ICT-HPC server          | GPU of ICT-HPC server          |
|----------------------|--------------|---------------|--------------|
| Input data & load data     | 682 ms  | 582 ms  | 442 ms  | 
| Prepare Data | -   | -   | 3.83 s    |
| Preprocessing | 385 ms   | 337 ms   | 91.7 ms    |
| Normalization & Scaling | 140 ms | 135 ms | 4.2 ms |
| Select Most Variable Genes | 778 ms | 1.14 s | 500 ms | 
| Regress out  | 14.9 s | 12.8 s | 2min 39s |
| PCA | 805 ms | 4min 15s | 786 ms |
| Clustering | 3.88 s | 9.36 s | 18.2 s | 
| UMAP | 7.88 s | 1min 34s | 357 ms |
| Finding marker genes & DE analysis | 2.65 s | 1.88 s | 2.31 s |
| Cell type identification | 762 ms | 659 ms | 689 ms |
| Visualization | 2.28 s | 2.83 s | 3.12 s |
| Save file | 312 ms | 210 ms | 210 ms |

![Alt text](https://github.com/vclabsysbio/AI-MD_scRNAseq/blob/main/Downstream_analysis/Figures/Comparison_Scanpy_on_CPU_vs_GPU.jpg?raw=true "Comparison")

**246 server**
- anaconda3
- Jupyter notebook

**ICT-HPC server**
- image: kubeflow/ngc-rapidsai:21.10-cuda11.0-runtime-ubuntu20.04
- Jupyterlab


## References
Stuart, T., Butler, A., Hoffman, P., Hafemeister, C., Papalexi, E., Mauck III, W. M., ... & Satija, R. 
(2019). *Comprehensive integration of single-cell data. Cell, 177*(7), 1888-1902.

Wolf, F. A., Angerer, P., & Theis, F. J. (2018). SCANPY: large-scale single-cell gene expression 
data analysis. *Genome biology, 19*(1), 1 
