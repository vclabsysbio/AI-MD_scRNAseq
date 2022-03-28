# Downstream analysis
- [**Tools**](#Tools)
- [**Data analysis pipeline**](#Data-analysis-pipeline)
- [**Commands**](#Commands)
- [**CPU and GPU performance testing**](#CPU-and-GPU-performance-testing)

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
### Input data & load data
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

### Input data & load data
``` {python}
adata = sc.read_10x_mtx(
    './data/',  # the directory with the `.mtx` file
    var_names='gene_symbols',
    cache=True) 
adata.var_names_make_unique()
```
### Prepare Data (GPU only)
``` {python}

```

``` {python}

```



## CPU and GPU performance testing
### Comparison


| tools                | CPU of 246 server        | CPU of ICT-HPC server          | GPU of ICT-HPC server          |
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

![Alt text](https://github.com/vclabsysbio/AI-MD_scRNAseq/blob/main/Downstream_analysis/Picture1.jpg?raw=true "Comparison")

**246 server**
- anaconda3
- Jupyter notebook

**ICT-HPC server**
- image: kubeflow/ngc-rapidsai:21.10-cuda11.0-runtime-ubuntu20.04
- Jupyterlab

