# Setting up the Sun Grid Environment
#$ -S /bin/bash
#$ -V
#$ -N wgs_vep
#$ -j y
#$ -o projects/wgs/qsub_out/$JOB_NAME.o$JOB_ID
#$ -l mem_free=5G
#$ -m ae
#$ -M jonathan.chung@phd.einstein.yu.edu

echo "==============================="
echo "starting on        : $(date)"
echo "Running on node    : $(hostname)"
echo "Current job ID     : $JOB_ID"
echo "Current job name   : $JOB_NAME"
echo "==============================="

perl ~/programs/ensembl-tools-release-75/scripts/variant_effect_predictor/variant_effect_predictor.pl \
    -i ~/projects/wgs/data/wgs_phased_indel.vcf \
    --cache \
    --no_progress \
    --buffer_size 20000 \
    --offline \
    --sift b \
    --polyphen b \
    --symbol \
    --numbers \
    --biotype \
    --total_length \
    -o ~/projects/wgs/results/vep/vep.wgs_phased_indel.vcf \
    --vcf \
    --fields Consequence,Codons,Amino_acids,Gene,SYMBOL,Feature,EXON,PolyPhen,SIFT,Protein_position,BIOTYPE

