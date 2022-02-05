# AI-MD single-cell RNA-seq

- [Requirements](#Requirements)
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

### Datasets

1. GEX
   - Library kit - Chromium Next GEM Single Cell 5p RNA library v2
   - Total Read1 + Read2 = 2,372,137,918 reads
2. TCR
   - Library kit - Chromium Next GEM Single Cell 5p RNA library v2
   - Total Read1 + Read2 = 141,348,498 reads
3. BCR
   - Library kit - Chromium Next GEM Single Cell 5p RNA library v2
   - Total Read1 + Read2 = 62,366,650 reads
4. GSA
   - Platform - Illumina Infinium SNP Genotyping Array (GSA MG v2)
   - Total SNPs 766,221 SNPs

### References
- hg19 To Hg38 over chain (December 31, 2013) [Download](http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz)
- Human Genome GRCh38 Reference (January 15, 2014) [Download](https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.fa.gz)
- GEX Reference - Human GRCh38 2020-A (July 7, 2020) [Download](https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2020-A.tar.gz)
- VDJ Reference - Human GRCh38 5.0.0 (November 19, 2020) [Download](https://cf.10xgenomics.com/supp/cell-vdj/refdata-cellranger-vdj-GRCh38-alts-ensembl-5.0.0.tar.gz)

## CPU specifications

### ICBS server Specifications
- **IP address:** 10.90.202.246
- **Vendor id:** GenuineIntel
- **Model name:** Intel Core Processor (Broadwell, IBRS)
- **Operating system:** Ubuntu 16.04.6 LTS (xenial)
- **CPUs:** 39
- **Memory (RAM):** 480GB

### ICT-HPC server Specifications 
- **IP address:** 10.134.1.9
- **Vendor id:** AuthenticAMD
- **Model name:** AMD EPYC 7742 64-Core Processor
- **Operating System:** Ubuntu 20.04.2 LTS (focal)
- **CPUs:** 256
- **Memory (RAM):** 2TB

## GPU specifications

### ICT-HPC server Specifications
- **IP address:** 10.134.1.9
- **Vendor id:** AuthenticAMD
- **Model name:** AMD EPYC 7742 64-Core Processor
- **Operating System:** Ubuntu 20.04.2 LTS (focal)
- **CPUs:** 256
- **Memory (RAM):** 2TB
