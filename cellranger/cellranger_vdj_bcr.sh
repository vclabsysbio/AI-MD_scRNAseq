cd /output-path && \

ref_vdj="/path/refdata-cellranger-vdj-GRCh38-alts-ensembl-5.0.0" && \
fastq_path="/path-of-BCR-fastq-files/" && \
id_name="idname" && \
sample_name="prefix-of-BCR-fastq-files" && \

time cellranger vdj --id=$id_name \
--reference=$ref_vdj \
--fastqs=$fastq_path \
--sample=$sample_name \
--localcores=8 \
--localmem=64 \
--chain=IG
