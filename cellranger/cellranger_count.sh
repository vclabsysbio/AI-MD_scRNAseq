cd /output-path && \

ref_gex="/path/refdata-gex-GRCh38-2020-A" && \
fastq_path="/path-of-fastq-files/" && \
id_name="idname" && \
sample_name="prefix-of-fastq-file" && \
expect_cells="22000" && \

time cellranger count --id=$id_name \
--transcriptome=$ref_gex \
--fastqs=$fastq_path \
--sample=$sample_name \
--expect-cells=$expect_cells \
--localcores=8 \
--localmem=64 \
--chemistry=SC5P-R2
