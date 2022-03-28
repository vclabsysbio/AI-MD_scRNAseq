# AI-MD single-cell RNA-seq

## Description
   This GitHub respiratory includes tools and commands for analyzing single-cell RNA-seq (scRNA-seq) on a general or HPC server using CPU and GPU. The main analysis processes are as follows: 1. Mapping sequencing data with a human genome reference for counting gene expression and immune receptor (TCR/BCR), 2. Processing genotype data and demultiplexing pooled scRNA-seq, and 3. Data analysis and visualization.


- [Requirements](#Requirements)
- [Running steps](#Running-steps)
- [CPU specifications](#CPU-specifications)
- [GPU specifications](#GPU-specifications)

## Requirements
### Tools
- CrossMap (v0.6.1)
- vcftools (v0.1.16)
- bcftools (v1.11)
- samtools (v1.10)
- popscle (include _demuxlet_ v2) [GitHub](https://github.com/statgen/popscle)
- Cellranger (v6.1.2)
- Python (v3.7)
- Scanpy (v1.8.2) (Python package)

### Datasets
- GEX
   - Library kit - Chromium Next GEM Single Cell 5p RNA library v2
   - Total Read1 + Read2 = 2,372,137,918 reads
- TCR
   - Library kit - Chromium Next GEM Single Cell 5p RNA library v2
   - Total Read1 + Read2 = 141,348,498 reads
- BCR
   - Library kit - Chromium Next GEM Single Cell 5p RNA library v2
   - Total Read1 + Read2 = 62,366,650 reads
- GSA
   - Platform - Illumina Infinium SNP Genotyping Array (GSA MG v2)
   - Total SNPs 766,221 SNPs

### References
- hg19 To Hg38 over chain (December 31, 2013) [Download](http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz)
- Human Genome GRCh38 Reference (January 15, 2014) [Download](https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.fa.gz)
- GEX Reference - Human GRCh38 2020-A (July 7, 2020) [Download](https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2020-A.tar.gz)
- VDJ Reference - Human GRCh38 5.0.0 (November 19, 2020) [Download](https://cf.10xgenomics.com/supp/cell-vdj/refdata-cellranger-vdj-GRCh38-alts-ensembl-5.0.0.tar.gz)

## Running steps
1. [**Mapping (Using Cellranger)**](https://github.com/vclabsysbio/AI-MD_scRNAseq/tree/main/cellranger)
2. [**Demultiplexing**](https://github.com/vclabsysbio/AI-MD_scRNAseq/tree/main/popscle)
   1. Liftover (Using CrossMap)
   2. VCF preprocessing
   3. VCF filtering
   4. Demutiplexing (Using popscle)
3. [**Downstream analysis**](https://github.com/vclabsysbio/AI-MD_scRNAseq/tree/main/Downstream_analysis)
   1. Preprocessing
   2. Normalization & Scaling the data
   3. Select Most Variable Genes
   4. Regress out confounding factors
   5. Perform linear dimensional reduction
   6. Clustering
   7. Run non-linear dimensional reduction (UMAP)
   8. Finding marker genes & Differential expression analysis
   9. Cell type identification
   10. Visualization

## CPU specifications
### ICBS server specifications
- **IP address:** 10.90.202.246
- **Vendor id:** GenuineIntel
- **Model name:** Intel Core Processor (Broadwell, IBRS)
- **Operating system:** Ubuntu 16.04.6 LTS (xenial)
- **CPUs:** 39
- **Memory (RAM):** 480GB

### ICT-HPC server specifications 
- **IP address:** 10.134.1.9
- **Vendor id:** AuthenticAMD
- **Model name:** AMD EPYC 7742 64-Core Processor
- **Operating System:** Ubuntu 20.04.2 LTS (focal)
- **CPUs:** 256
- **Memory (RAM):** 2TB

## GPU specifications
### ICT-HPC server specifications
- **IP address:** 10.134.1.9
- **Vendor id:** AuthenticAMD
- **Model name:** AMD EPYC 7742 64-Core Processor
- **Operating System:** Ubuntu 20.04.2 LTS (focal)
- **CPUs:** 256
- **Memory (RAM):** 2TB
-  **GPUs:** NVIDIA DGX A100 8 cards
