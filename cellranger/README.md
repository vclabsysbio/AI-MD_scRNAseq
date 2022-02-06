# Mapping
- [**Tools**](#Tools)
- [**Commands**](#Commands)
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
```
cellranger vdj --id=$id_name \
               --reference=$ref_vdj \
               --fastqs=$fastq_path \
               --sample=$sample_name \
               --localcores=8 \
               --localmem=64 \
               --chain=TR
```

- **BCR mapping**

```
cellranger vdj --id=$id_name \
               --reference=$ref_vdj \
               --fastqs=$fastq_path \
               --sample=$sample_name \
               --localcores=8 \
               --localmem=64 \
               --chain=IG
```

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
