# Setting up the Sun Grid Environment
#$ -S /bin/bash
#$ -V
#$ -N wgs_annovar
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

 vcftools \
 --vcf ~/projects/wgs/data/wgs_phased_snp_pass.recode.vcf \
 --chr chr22 \
 --recode \
 --recode-INFO-all \
 --out ~/projects/wgs/data/for_deyou/wgs_phased_snp_chr22
 
 

