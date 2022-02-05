# Mapping

## CPU performance testing
### ICBS server

Settings
- CPUs: 8
- Memory: 64GB


| tools                | real         | user          | sys          |
|----------------------|--------------|---------------|--------------|
| Cellranger count     | 1233m38.423s | 6615m45.040s  | 248m10.156s  | 
| (TCR) Cellranger vdj |
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
