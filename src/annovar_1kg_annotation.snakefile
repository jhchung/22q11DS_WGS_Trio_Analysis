import re

workdir: "/home/jchung/projects/wgs/"
annovar_dir = '/home/jchung/programs/annovar/annovar2013Aug23/'
vcftools = "/apps1/vcftools/vcftools_0.1.11/cpp/vcftools"
vcfconcat = "/apps1/vcftools/vcftools_0.1.11/perl/vcf-concat"
chromosome = "1 2 3 4 5 6 7 8 10 11 12 14 15 16 17 18 20 22 X".split()

rule all:
    input:
        "/cork/jchung/wgs/results/variant_counts/1kg/tbx1_variants.1kg_eur.all_counts.txt"
    
rule combine_counts:
    input:
        dynamic("/cork/jchung/wgs/results/variant_counts/1kg/tbx1_variants.{SAMPLEID}.counts.txt")
    output:
        "/cork/jchung/wgs/results/variant_counts/1kg/tbx1_variants.1kg_eur.all_counts.txt"
    shell: """
    Rscript src/combine_variant_counts.R \
    {output} \
    {input} 
    """
    
rule count_1kg_variants:
    input:
        "/cork/jchung/wgs/results/annovar/1kg_EUR/tbx1_variants.{SAMPLEID}.hg19_multianno.txt",
        "data/tbx1_pathway_genes_2014_03_27.txt"
    output:
        "/cork/jchung/wgs/results/variant_counts/1kg/tbx1_variants.{SAMPLEID}.counts.txt"
    shell: """
    perl src/perl/count_var_per_gene.pl \
    {input[1]} \
    {input[0]} \
    {output}
    """
    
rule annotate_sample_variants:
    input:
        "/cork/jchung/wgs/data/1kg_vcf_filtered_avinput/tbx1_variants.{SAMPLEID}.avinput"
    params:
        output_prefix = "/cork/jchung/wgs/results/annovar/1kg_EUR/tbx1_variants.{SAMPLEID}"
    output: 
        "/cork/jchung/wgs/results/annovar/1kg_EUR/tbx1_variants.{SAMPLEID}.hg19_multianno.txt"
    shell: """
    {annovar_dir}/table_annovar.pl \
    {input} \
    {annovar_dir}/humandb/ \
    --protocol refGene \
    --operation g \
    --buildver hg19 \
    --outfile {params.output_prefix} \
    --otherinfo
    """
    
# Use `--format vcf4` to extract only variants with alternate alleles
rule convert_vcf_to_avinput:
    input:
        "/cork/jchung/wgs/data/1kg_vcf_filtered/tbx1_variants_allchr.vcf"
    params:
        output = "/cork/jchung/wgs/data/1kg_vcf_filtered_avinput/tbx1_variants"
    output:
        dynamic("/cork/jchung/wgs/data/1kg_vcf_filtered_avinput/tbx1_variants.{SAMPLEID}.avinput")
    shell: """
    {annovar_dir}/convert2annovar.pl \
    --format vcf4 {input} \
    --outfile {params.output} \
    --allsample \
    --include \
    --comment
    """
    
rule combine_vcf_files:
    input:
        expand("/cork/jchung/wgs/data/1kg_vcf_filtered/tbx1_variants_chr{CHROM}.recode.vcf", CHROM = chromosome)
    output:
        "/cork/jchung/wgs/data/1kg_vcf_filtered/tbx1_variants_allchr.vcf"
    shell: """
    perl {vcfconcat} \
    {input} \
    > {output}
    """
    
rule extract_variants_from_vcf:
    input:
        "/cork/jchung/wgs/data/1kg_vcf/ALL.chr{CHROM}.phase1_release_v2.20101123.snps_indels_svs.vcf",
        "/cork/jchung/wgs/results/annovar/1kg/tbx1_pathway_variant_positions.txt",
        "data/1kg_EUR_subjects.txt"
    params:
        output = "/cork/jchung/wgs/data/1kg_vcf_filtered/tbx1_variants_chr{CHROM}"
    output:
        "/cork/jchung/wgs/data/1kg_vcf_filtered/tbx1_variants_chr{CHROM}.recode.vcf"
    shell: """
    {vcftools} \
    --vcf {input[0]} \
    --out {params.output} \
    --positions {input[1]} \
    --keep {input[2]} \
    --recode \
    --recode-INFO-all 
    """
    
rule get_variant_positions:
    input:
        "/cork/jchung/wgs/results/annovar/1kg/1kg_sites.hg19_multianno.genic.tbx1.txt",
    output:
        "/cork/jchung/wgs/results/annovar/1kg/tbx1_pathway_variant_positions.txt"
    shell: """
    awk '{{print $1"\t"$2}}' {input} > {output}
    """
    
rule extract_tbx1_genes:
    input:
        "/cork/jchung/wgs/results/annovar/1kg/1kg_sites.hg19_multianno.genic.txt",
        "data/tbx1_pathway_genes_2014_03_27.txt"
    output:
        "/cork/jchung/wgs/results/annovar/1kg/1kg_sites.hg19_multianno.genic.tbx1.txt"
    run:
        variant_file = open(input[0], "r")
        out_file = open(output[0], "w")
        
        with open(input[1], "r") as f:
            genes = f.read().splitlines()
            
        for variant_annotation in variant_file:
            if any(gene in variant_annotation for gene in genes):
                out_file.write(variant_annotation)
    
rule all_lod100:
    input: "/cork/jchung/wgs/results/annovar/1kg/1kg_sites.hg19_multianno.txt"
    output: "/cork/jchung/wgs/results/annovar/1kg/1kg_sites.hg19_multianno.lod100.txt"
    run:
        in_file = open(input[0], "r")
        out_file = open(output[0], "w")
        lod_index = None
        
        for line in in_file:
            if line.startswith("Chr\tStart"):
                split_line = line.split("\t")
                lod_index = split_line.index("phastConsElements46way")
            else:
                split_line = line.split("\t")
                lod_column = split_line[lod_index]
                
                if not lod_column == "":
                    lod_score = lod_column.split(";")[1]
                    lod_score = lod_score.split("=")[2]
                    if int(lod_score) >= 100:
                        out_file.write(line)
        in_file.close()
        out_file.close()
    
rule genic_lod100:
    input: "/cork/jchung/wgs/results/annovar/1kg/1kg_sites.hg19_multianno.genic.txt"
    output: "/cork/jchung/wgs/results/annovar/1kg/1kg_sites.hg19_multianno.genic.lod100.txt"
    run:
        in_file = open(input[0], "r")
        out_file = open(output[0], "w")
        lod_index = None
        
        for line in in_file:
            if line.startswith("Chr\tStart"):
                split_line = line.split("\t")
                lod_index = split_line.index("phastConsElements46way")
            else:
                split_line = line.split("\t")
                lod_column = split_line[lod_index]
                
                if not lod_column == "":
                    lod_score = lod_column.split(";")[1]
                    lod_score = lod_score.split("=")[2]
                    if int(lod_score) >= 100:
                        out_file.write(line)
        in_file.close()
        out_file.close()
            
rule extract_genic:
    input: "/cork/jchung/wgs/results/annovar/1kg/1kg_sites.hg19_multianno.txt"
    output: "/cork/jchung/wgs/results/annovar/1kg/1kg_sites.hg19_multianno.genic.txt"
    run:
        in_file = open(input[0], "r")
        out_file = open(output[0], "w")
        
        re_pattern = re.compile("exonic|splicing|UTR5|URT3|intronic")
        for line in in_file:
            if line.startswith("Chr\tStart"):
                out_file.write(line)
                
            elif not line.startswith("#"):
                if re_pattern.search(line):
                    out_file.write(line)
        
        in_file.close()
        out_file.close()
        
rule table_annotation:
    input: "/cork/jchung/wgs/data/for_annovar/1kg/1kg_sites.avinput"
    params: output_prefix = "/cork/jchung/wgs/results/annovar/1kg/1kg_sites"
    output: "/cork/jchung/wgs/results/annovar/1kg/1kg_sites.hg19_multianno.txt"
    shell: """
    {annovar_dir}/table_annovar.pl \
    {input} \
    {annovar_dir}/humandb/ \
    --protocol refGene,phastConsElements46way \
    --operation g,r \
    --buildver hg19 \
    --outfile {params.output_prefix} \
    --otherinfo
    """
    
# Use `--format vcf4old` to keep all variants regardless of genotype
rule convertToAnnovar:
    input: "/cork/jchung/wgs/data/1kg_vcf/ALL.wgs.phase1_release_v2.20101123.snps_indels_sv.sites.vcf"
    output: "/cork/jchung/wgs/data/for_annovar/1kg/1kg_sites.avinput"
    message: "Converting {input} into format for annovar. Use '--format vcf4old` to keep all variants"
    shell: """
    {annovar_dir}/convert2annovar.pl \
    -format vcf4old {input} \
    -outfile {output} \
    -includeinfo \
    -filter PASS \
    -allallele \
    -comment
    """

