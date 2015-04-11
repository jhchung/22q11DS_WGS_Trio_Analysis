# Setting up the Sun Grid Environment
#$ -S /bin/bash
#$ -V
#$ -N wgs_annovar
#$ -j y
#$ -o projects/wgs/qsub_out/$JOB_NAME.o$JOB_ID
#$ -l mem_free=40G
#$ -m ae
#$ -M jonathan.chung@phd.einstein.yu.edu

echo "==============================="
echo "starting on        : $(date)"
echo "Running on node    : $(hostname)"
echo "Current job ID     : $JOB_ID"
echo "Current job name   : $JOB_NAME"
echo "==============================="

snakemake -s ~/projects/wgs/src/annovar_wgs_samples/annovar_annotation.snakefile \
    checkpoint1 \
    --forceall
