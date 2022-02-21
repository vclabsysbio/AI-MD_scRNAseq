# Mapping
- [**Tools**](#Tools)
- [**Commands**](#Commands)
  - Create Cellranger singularity file
  - GEX mapping
  - VDJ mapping
    - TCR
    - BCR
- [**CPU performance testing**](#CPU-performance-testing)
  - ICBS server
  - ICT-HPC server
- [**GPU performance testing**](#GPU-performance-testing)
  - ICT-HPC server

## Tools
- Cellranger (v6.1.2)

## Commands
### Create Cellranger singularity file
```
vim Dockerfile
sudo docker build -t cellranger612 .
sudo docker tag cellranger:latest tipsarin/cellranger612:latest
sudo docker push tipsarin/cellranger612:latest
```

```
singularity pull cellranger612.sif docker://docker.io/tipsarin/cellranger612
```

### GEX mapping

```
cellranger count --id=$id_name \
                 --transcriptome=$ref_gex \
                 --fastqs=$fastq_path \
                 --sample=$sample_name \
                 --expect-cells=$expect_cells \
                 --localcores=8 \
                 --localmem=64 \
                 --chemistry=SC5P-R2
```

### VDJ mapping
- **TCR mapping**

**Inputs**
- Fastq path - ex. `./HN00144497_10X_RawData_Outs/CVc_B2rxn2_GEX/HFTM7CCX2`
  - Fastq files
  - `CVc_B2rxn2_TCR_S3_L004_I1_001.fastq.gz`
  - `CVc_B2rxn2_TCR_S3_L004_I2_001.fastq.gz`
  - `CVc_B2rxn2_TCR_S3_L004_R1_001.fastq.gz`
  - `CVc_B2rxn2_TCR_S3_L004_R2_001.fastq.gz`
- VDJ reference path - ex. `./refdata-cellranger-vdj-GRCh38-alts-ensembl-5.0.0`
```
cellranger vdj --id=$id_name \
               --reference=$ref_vdj \
               --fastqs=$fastq_path \
               --sample=$sample_name \
               --localcores=8 \
               --localmem=64 \
               --chain=TR
```
**Outputs**
- Run summary HTML - ex. `web_summary.html`
- Run summary CSV - ex. `metrics_summary.csv`
- Clonotype info - ex. `clonotypes.csv`
- Filtered contig sequences FASTA - ex. `filtered_contig.fasta`
- Filtered contig sequences FASTQ - ex. `filtered_contig.fastq`
- Filtered contigs (CSV) - ex. `filtered_contig_annotations.csv`
- All-contig FASTA - ex. `all_contig.fasta`
- All-contig FASTA index - ex. `all_contig.fasta.fai`
- All-contig FASTQ - ex. `all_contig.fastq`
- Read-contig alignments - ex. `all_contig.bam`
- Read-contig alignment index - ex. `all_contig.bam.bai`
- All contig annotations (JSON) - ex. `all_contig_annotations.json`
- All contig annotations (BED) - ex. `all_contig_annotations.bed`
- All contig annotations (CSV) - ex. `all_contig_annotations.csv`
- Barcodes that are declared to be targetted cells - ex. `cell_barcodes.json`
- Clonotype consensus FASTA - ex. `consensus.fasta`
- Clonotype consensus FASTA index - ex. `consensus.fasta.fai`
- Contig-consensus alignments - ex. `consensus.bam`
- Contig-consensus alignment index - ex. `consensus.bam.bai`
- Clonotype consensus annotations (CSV) - ex. `consensus_annotations.csv`
- Concatenated reference sequences - ex. `concat_ref.fasta`
- Concatenated reference index - ex. `concat_ref.fasta.fai`
- Contig-reference alignments - ex. `concat_ref.bam`
- Contig-reference alignment index - ex. `concat_ref.bam.bai`
- Loupe V(D)J Browser file - ex. `vloupe.vloupe`



- **BCR mapping**
- 
**Inputs**
- Fastq path - ex. `./HN00144497_10X_RawData_Outs/CVc_B2rxn2_GEX/HFTM7CCX2`
  - Fastq files
  - `CVc_B2rxn2_BCR_S3_L004_I1_001.fastq.gz`
  - `CVc_B2rxn2_BCR_S3_L004_I2_001.fastq.gz`
  - `CVc_B2rxn2_BCR_S3_L004_R1_001.fastq.gz`
  - `CVc_B2rxn2_BCR_S3_L004_R2_001.fastq.gz`
- VDJ reference path - ex. `./refdata-cellranger-vdj-GRCh38-alts-ensembl-5.0.0`
```
cellranger vdj --id=$id_name \
               --reference=$ref_vdj \
               --fastqs=$fastq_path \
               --sample=$sample_name \
               --localcores=8 \
               --localmem=64 \
               --chain=IG
```
**Outputs**
- Run summary HTML - ex. `web_summary.html`
- Run summary CSV - ex. `metrics_summary.csv`
- Clonotype info - ex. `clonotypes.csv`
- Filtered contig sequences FASTA - ex. `filtered_contig.fasta`
- Filtered contig sequences FASTQ - ex. `filtered_contig.fastq`
- Filtered contigs (CSV) - ex. `filtered_contig_annotations.csv`
- All-contig FASTA - ex. `all_contig.fasta`
- All-contig FASTA index - ex. `all_contig.fasta.fai`
- All-contig FASTQ - ex. `all_contig.fastq`
- Read-contig alignments - ex. `all_contig.bam`
- Read-contig alignment index - ex. `all_contig.bam.bai`
- All contig annotations (JSON) - ex. `all_contig_annotations.json`
- All contig annotations (BED) - ex. `all_contig_annotations.bed`
- All contig annotations (CSV) - ex. `all_contig_annotations.csv`
- Barcodes that are declared to be targetted cells - ex. `cell_barcodes.json`
- Clonotype consensus FASTA - ex. `consensus.fasta`
- Clonotype consensus FASTA index - ex. `consensus.fasta.fai`
- Contig-consensus alignments - ex. `consensus.bam`
- Contig-consensus alignment index - ex. `consensus.bam.bai`
- Clonotype consensus annotations (CSV) - ex. `consensus_annotations.csv`
- Concatenated reference sequences - ex. `concat_ref.fasta`
- Concatenated reference index - ex. `concat_ref.fasta.fai`
- Contig-reference alignments - ex. `concat_ref.bam`
- Contig-reference alignment index - ex. `concat_ref.bam.bai`
- Loupe V(D)J Browser file - ex. `vloupe.vloupe`

## CPU performance testing
### ICBS server

Settings
- CPUs: 8
- Memory: 64GB


| tools                | real         | user          | sys          |
|----------------------|--------------|---------------|--------------|
| Cellranger count     | 1233m38.423s | 6615m45.040s  | 248m10.156s  | 
| (TCR) Cellranger vdj | 142m44.837s  | 479m20.684s   | 27m56.532s   |
| (BCR) Cellranger vdj | 84m51.065s   | 193m31.028s   | 18m19.740s   |



GEX reference: refdata-gex-GRCh38-2020-A

VDJ reference: refdata-cellranger-vdj-GRCh38-alts-ensembl-5.0.0

### ICT-HPC server

Settings
- CPUs: 8
- Memory: 64GB

| tools                | real         | user          | sys          |
|----------------------|--------------|---------------|--------------|
| Cellranger count     | 485m42.891s  | 2622m40.515s  | 127m33.903s  | 
| (TCR) Cellranger vdj | 55m32.961s   | 315m19.836s   | 9m13.907s    |
| (BCR) Cellranger vdj | 26m11.690s   | 148m12.192s   | 3m42.816s    |

## GPU performance testing
