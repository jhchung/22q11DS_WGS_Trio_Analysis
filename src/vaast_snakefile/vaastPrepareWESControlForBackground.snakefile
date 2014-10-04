# Analysis of whole genome sequencing data using VAAST

workdir: "/home/jchung/projects/wgs/"
vaast_sort_gff = "/home/jchung/programs/VAAST_Code_2.0.1/bin/vaast_tools/vaast_sort_gff"
vaast_converter = "/home/jchung/programs/VAAST_Code_2.0.1/bin/vaast_tools/vaast_converter"

refseqFeatures  = "/home/jchung/programs/VAAST_Code_2.0.1/data/hg19/Features/refGene_hg19.gff3"
refseqFeaturesSort = "/home/jchung/programs/VAAST_Code_2.0.1/data/hg19/Features/refGene_hg19.sorted.gff3"
hg19Fasta = "/home/jchung/programs/VAAST_Code_2.0.1/data/hg19/Fasta/vaast_hsap_chrs_hg19.fa"

wesFamFile = "/home/jchung/projects/WES/data/skat_test/wes186_skat_TOF4.fam"
wesIdForWgsSampleFile = "/home/jchung/projects/wgs/data/wesIdForWgsSamples.txt"

def removeValuesFromList(theList, val):
   return [value for value in theList if value != val]

import csv
caseSamples = []
controlSamples = []

with open(wesFamFile, "r") as wesFam:
    filereader = csv.reader(wesFam, delimiter = " ")
    for row in filereader:
        if row[5] == "1":
            controlSamples.append(row[1])
        elif row[5] == "2":
            caseSamples.append(row[1])

wesIdForWgsSample = []
with open(wesIdForWgsSampleFile, "r") as idFile:
    idFileReader = csv.reader(idFile, delimiter = "\t")
    for row in idFileReader:
        controlSamples = removeValuesFromList(controlSamples, row[1])
        caseSamples = removeValuesFromList(caseSamples, row[1])
    
rule all:
    input: "results/vaastOutput/wesBackground/VST/wesControls.cdr"
    
rule runVSTOnWesControls:
    input: expand("results/vaastOutput/wesBackground/VAT/refGene.{CONTROLSAMPLES}.vat.gvf", CONTROLSAMPLES = controlSamples)
    output: "results/vaastOutput/wesBackground/VST/wesControls.cdr"
    shell: """
    VST \
    -o 'U(0..100)' \
    -b hg19 \
    -chunk 100000000 \
    {input} \
    > {output} """

rule runVAT:
    input: "results/vaastOutput/wesBackground/GVF/{CONTROLSAMPLES}.gvf",
           refseqFeaturesSort,
           hg19Fasta
    output: "results/vaastOutput/wesBackground/VAT/refGene.{CONTROLSAMPLES}.vat.gvf"
    shell: """
    VAT \
    --features {input[1]} \
    --build hg19 \
    -c 100000000 \
    --fasta {input[2]} \
    {input[0]} \
    > {output}"""

rule sortGff:
    input: refseqFeatures
    output: refseqFeaturesSort
    shell: """
    {vaast_sort_gff} \
    {input} > {output}"""
    
rule convertVCF2GVF:
    input: "data/forVAAST/wes186_snp_pass_autosomal.recode.vcf"
    params: outputDir = "results/vaastOutput/wesBackground/GVF/"
    output: expand("results/vaastOutput/wesBackground/GVF/{CONTROLSAMPLES}.gvf", CONTROLSAMPLES = controlSamples)
    shell: """
    {vaast_converter} \
    --build hg19 \
    --path {params.outputDir} \
    {input}
    """
    
rule filterWgsFile:
    input: "/home/jchung/dat/whole_exome_sequencing/morrow_rsng_r186.batch_1.all.polymorphic.filtered.snps.vcf",
           "data/forVAAST/refseq_gene_hg19_20131111_wes_format.bed"
    params: outputPrefix = "data/forVAAST/wes186_snp_pass_autosomal"
    output: "data/forVAAST/wes186_snp_pass_autosomal.recode.vcf"
    shell: """
    {vcftools} \
    --remove-filtered-all \
    --vcf {input[0]} \
    --not-chr X \
    --not-chr Y \
    --out {params.outputPrefix} \
    --recode \
    --recode-INFO-all
    """
    
rule formatUCSCBedFile:
    input: "data/forVAAST/refseq_gene_hg19_20131111.bed"
    output: "data/forVAAST/refseq_gene_hg19_20131111_wes_format.bed"
    shell: """
    Rscript src/formatUcscBed.R {input} {output}
    """
    
