cd /home/sarintip.ngu/Covid_vac/GEX/cellranger_outputs

ref_gex="/home/sarintip.ngu/References/refdata-gex-GRCh38-2020-A"

fastq_path="/home/sarintip.ngu/Covid_vac/GEX/GEX_CV2_1/H37NYDSX3"

id_name="CovidVac_B2_rxn1_612"

sample_name="GEX_CV2_1"

time cellranger count --id=$id_name \
--transcriptome=$ref_gex \
--fastqs=$fastq_path \
--sample=$sample_name \
--expect-cells=20000 \
--localcores=8 \
--localmem=64 \
--chemistry=SC5P-R2 \
--r1-length=26
