cd /path/ && \

export VCF_OLD_FILENAME=filename
export LIFTOVER_FILE=/path/filename
export HUMAN_GENOME_REF_FILE=/path/filename
export VCF_FILENAME=filename && \
export NEW_VCF_FILENAME=filename && \
export BAM_FILE=/path/filename && \
export REF_FASTA_FILE=/path/filename && \
export SAMPLE_LIST_FILE=/path/filename && \
export BARCODE_FILE=/path/filename && \
export DEMUXLET_OUTPUT=/path/filename && \

echo "Liftover hg19 to GRCh38" && \
time CrossMap.py vcf $LIFTOVER_FILE ./${VCF_OLD_FILENAME}.vcf $HUMAN_GENOME_REF_FILE ${VCF_FILENAME}.vcf && \

echo "Files for preprocessiong" && \
echo "VCF file = $VCF_FILENAME" && \
echo "Bam file (Output from cellranger) = $BAM_FILE" && \
echo "Sample list = $SAMPLE_LIST_FILE" && \
cat $SAMPLE_LIST_FILE && \
echo "Fasta file from Human genome reference = $REF_FASTA_FILE" && \

echo "Files for demultiplexing" && \
echo "Bam file (Output from cellranger) = $BAM_FILE" && \
echo "Sample list = $SAMPLE_LIST_FILE" && \
cat $SAMPLE_LIST_FILE && \
echo "Filtered barcode file (Output from cellranger) = $BARCODE_FILE" && \
echo "Output file = $DEMUXLET_OUTPUT" && \

echo "Step 1.1: print Chromosome name in VCF file" && \
awk '{print $1}' ./${VCF_FILENAME}.vcf | uniq && \

echo "Step 1.2: print Chromosome name in BAM file" && \
samtools idxstats $BAM_FILE | cut -f 1 | head -26 && \

echo "Step 2: Add 'chr' prefix in VCF header" && \
bcftools view --header-only ./${VCF_FILENAME}.vcf | sed 's/##contig=<ID=/##contig=<ID=chr/' > ./2.1${VCF_FILENAME}_header.vcf && \

echo "Step 3: Add 'chr' prefix in VCF body" && \
bcftools view --no-header ./${VCF_FILENAME}.vcf | sed 's/^/chr/' | sed 's/^chrMT/chrM/' > ./2.2${VCF_FILENAME}_body.vcf && \

echo "Step 4: Sort chromosome name (e.g. 1, 11, 12) and edit in VCF header" && \
head -5 ./2.1${VCF_FILENAME}_header.vcf > ./2.1-1${VCF_FILENAME}_header_top.vcf && \

sed -n '6,326p' ./2.1${VCF_FILENAME}_header.vcf > ./2.1-2${VCF_FILENAME}_header_center.vcf && \
sort ./2.1-2${VCF_FILENAME}_header_center.vcf > ./2.1-2${VCF_FILENAME}_header_center_sorted.vcf && \

sed -n '455p' ./2.1${VCF_FILENAME}_header.vcf > ./2.1-2${VCF_FILENAME}_header_center_chrX.vcf && \
sed -n '459p' ./2.1${VCF_FILENAME}_header.vcf > ./2.1-2${VCF_FILENAME}_header_center_chrY.vcf && \
cat ./2.1-2${VCF_FILENAME}_header_center_chrX.vcf ./2.1-2${VCF_FILENAME}_header_center_chrY.vcf > ./2.1-2${VCF_FILENAME}_header_center_chrX_Y.vcf && \

sed -n '327p' ./2.1${VCF_FILENAME}_header.vcf > ./2.1-2${VCF_FILENAME}_header_center_chrM.vcf && \
sed -n '461,466p' ./2.1${VCF_FILENAME}_header.vcf > ./2.1-3${VCF_FILENAME}_header_bottom.vcf && \

cat ./2.1-1${VCF_FILENAME}_header_top.vcf ./2.1-2${VCF_FILENAME}_header_center_sorted.vcf ./2.1-3${VCF_FILENAME}_header_bottom.vcf > ./3.1${VCF_FILENAME}_header_edited.vcf && \

echo "Step 5: Sort chromosome name (e.g. 1, 11, 12) and edit in VCF body" && \

vcf-sort ./2.2${VCF_FILENAME}_body.vcf > ./3.2${VCF_FILENAME}_body_sorted.vcf && \

echo "Step 6: Merge VCF edited files (VCF header and VCF body)" && \
cat ./3.1${VCF_FILENAME}_header_edited.vcf ./3.2${VCF_FILENAME}_body_sorted.vcf > ./4${VCF_FILENAME}_chrname_edited_sorted.vcf && \

rm ./2.1${VCF_FILENAME}_header.vcf && \
rm ./2.1-1${VCF_FILENAME}_header_top.vcf && \
rm ./2.1-2${VCF_FILENAME}_header_center.vcf && \
rm ./2.1-2${VCF_FILENAME}_header_center_sorted.vcf &&\
rm ./2.1-2${VCF_FILENAME}_header_center_chrX_Y.vcf && \
rm ./2.1-2${VCF_FILENAME}_header_center_chrM.vcf && \
rm ./2.1-3${VCF_FILENAME}_header_bottom.vcf && \

rm ./2.2${VCF_FILENAME}_body.vcf && \

rm ./3.1${VCF_FILENAME}_header_edited.vcf && \
rm ./3.2${VCF_FILENAME}_body_sorted.vcf && \

echo "Step 7: Subset chromosomes" && \
bcftools reheader ./4${VCF_FILENAME}_chrname_edited_sorted.vcf -f $REF_FASTA_FILE > ./5${VCF_FILENAME}_subsetchr.vcf && \

echo "Check your file:Remove some chromosomes (e.g.chr14_GL000009v2_random)" && \
vcftools --vcf ./5${VCF_FILENAME}_subsetchr.vcf --chr chr1 --chr chr10 --chr chr11 --chr chr12 --chr chr13 --chr chr14 --chr chr15 --chr chr16 --chr chr17 --chr chr18 --chr chr19 --chr chr2 --chr chr20 \
    --chr chr21 --chr chr22 --chr chr3 --chr chr4 --chr chr5 --chr chr6 --chr chr7 --chr chr8 --chr chr9 --chr chrM --chr chrX --chr chrY --recode --recode-INFO-all --out 5${VCF_FILENAME}_subsetchr_body && \
mv ./5${VCF_FILENAME}_subsetchr_body.recode.vcf ./5${VCF_FILENAME}_subsetchr.vcf && \

echo "Step 8: Subset Samples" && \
bcftools view -S $SAMPLE_LIST_FILE -o ./6${VCF_FILENAME}_subsetsamples.vcf ./5${VCF_FILENAME}_subsetchr.vcf && \

rm ./4${VCF_FILENAME}_chrname_edited_sorted.vcf && \
rm ./5${VCF_FILENAME}_subsetchr.vcf && \

echo "Remove indel,  0/0 1/1 in all samples"
time vcftools --vcf ./6${VCF_FILENAME}_subsetsamples.vcf \
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
         --out 7${VCF_FILENAME}_rm_indel
mv ./7${VCF_FILENAME}_rm_indel.recode.vcf ./${NEW_VCF_FILENAME}_excluded.vcf && \

rm ./6${VCF_FILENAME}_subsetsamples.vcf && \

echo "VCF Preprocessing Finish !!" && \

echo "dumuxlet" && \
time popscle demuxlet --sam $BAM_FILE \
	--tag-group CB \
	--tag-UMI UB \
	--vcf ./${NEW_VCF_FILENAME}_excluded.vcf \
	--field GT \
	--out $DEMUXLET_OUTPUT \
	--group-list $BARCODE_FILE \
	--sm-list $SAMPLE_LIST_FILE 
#	--alpha 0.5
#	--min-snp 50 \

echo "!!!! Demultiplexing Finish !!!!"
echo "Please check your output and log file"
