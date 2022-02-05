# Demutiplexing
- [Tools](#Tools)
- [Commands](#Commands)
- [CPU performance testing](#CPU-performance-testing)
- [GPU performance testing](#GPU-performance-testing)

## Tools
- CrossMap (v0.6.1)
- vcftools (v0.1.16)
- bcftools (v1.11)
- samtools (v1.10)
- popscle (include _demuxlet_ v2) [GitHub](https://github.com/statgen/popscle)

## Commands
- **CrossMap** (optional)
- **VCF proprocessing** (optional)
- **VCF filtering** (required)
- **Demultiplexing**

### CrossMap
**Inputs**
- Liftover file - ex. `./hg19ToHg38.over.chain.gz`
- Genotyping data vcf filename - ex. `HB00004766`
- Human reference genome fasta file - ex. `./hg38.fa`
- Outputs filename - ex. `HB00004766.out.hg38`
```
CrossMap.py vcf $LIFTOVER_FILE ./${VCF_OLD_FILENAME}.vcf $HUMAN_GENOME_REF_FILE ${VCF_FILENAME}.vcf
```
**Outputs**
- Map vcf file - ex. `HB00004766.out.hg38.vcf`
- Unmap vcf file - ex. `HB00004766.out.hg38.vcf.unmap`

### VCF filtering
**Inputs**
- Map vcf file - ex. `HB00004766.out.hg38.vcf`
```
vcftools --vcf ./6${VCF_FILENAME}_subsetsamples.vcf \
         --remove-indels \
         --max-missing-count 0 \
         --maf 0.01 \
         --max-maf 0.99 \
         --max-alleles 2 \
         --not-chr chrM \
         --not-chr chrX \
         --not-chr chrY \
         --remove-filtered-geno-all \
         --recode \
         --recode-INFO-all \
         --out 7${VCF_FILENAME}_excluded
```
**Outputs**
- vcf file - ex. `HB00004766.out.hg38_excluded.recode.vcf`

### VCF proprocessing
**Inputs**
```
(optional)
```
**Outputs**

### Demultiplexing
**Inputs**
- BAM file - ex. `./possorted_genome_bam.bam`
- vcf file - ex. `./HB00004766.out.hg38_excluded.recode.vcf`
- Barcode file - ex. `./barcodes.tsv`
- Sample list file - ex. `./sample_name_CVc_B2.txt`
- Output filename - ex. `CVc_B2rxn2`
```
popscle demuxlet --sam $BAM_FILE \
	         --tag-group CB \
	         --tag-UMI UB \
	         --vcf ./${NEW_VCF_FILENAME}_excluded.vcf \
	         --field GT \
	         --out $DEMUXLET_OUTPUT \
	         --group-list $BARCODE_FILE \
	         --sm-list $SAMPLE_LIST_FILE 
```
**Outputs**
- Demultiplexed file - ex. `CVc_B2rxn2.best`

## CPU performance testing
### ICBS server
Settings
- CPUs: xx
- Memory: xxxGB

| tools                | real         | user          | sys          |
|----------------------|--------------|---------------|--------------|
| CrossMap             | 1m6.421s     | 0m33.508s     | 0m15.204s    | 
| VCF preprocessing    | 0m29.314s    | 0m23.636s     | 0m6.472s     | 
| VCF filtering        | 0m7.860s     | 0m7.592s      | 0m0.236s     |
| demuxlet             | 1666m37.889s | 1583m21.216s  | 27m49.324s   |


### ICT-HPC server
Settings
- CPUs: 64
- Memory: 800GB

| tools                | real         | user          | sys          |
|----------------------|--------------|---------------|--------------|
| CrossMap             | 0m21.475s    | 0m17.708s     | 0m5.220s     | 
| VCF preprocessing    | 0m20.619s    | 0m16.475s     | 0m3.030s     | 
| VCF filtering        | 0m5.754s     | 0m5.492s      | 0m0.236s     |
| demuxlet             | 133m46.388s  | 132m36.092s   | 0m29.843s    |

## GPU performance testing
