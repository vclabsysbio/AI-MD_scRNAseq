# Demutiplexing
- [Tools](#Tools)
- [Data analysis pipeline](#Data-analysis-pipeline)
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

### Create 
```

```

## CPU and GPU performance testing
### Comparison



**246 server**
- anaconda3
- Jupyter notebook

**ICT-HPC server**
- image: kubeflow/ngc-rapidsai:21.10-cuda11.0-runtime-ubuntu20.04
- Jupyterlab

