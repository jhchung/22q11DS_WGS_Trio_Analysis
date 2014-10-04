import re

workdir: "/home/jchung/projects/wgs/"
variantFiles = "wgs.k1_only.maf1.deleterious.variant.txt wgs.k2_only.maf1.deleterious.variant.txt wgs.maf1.deleterious.variant.txt".split()
annovar_dir = '/home/jchung/programs/annovar/annovar2013Aug23/'
vcfFiles = "wgs_phased_indel.vcf wgs_phased_snp.vcf".split()
vcftools = "/apps1/vcftools/vcftools_0.1.11/cpp/vcftools"
sampleid = "BM1452.001 BM1452.100 BM1452.200 BM1453.001 BM1453.100 BM1453.200".split()
var_type = "snp indel".split()

def extract_conservation(input, output, conservation_column, score_cutoff):
    in_file = open(input, "r")
    out_file = open(output, "w")
    lod_index = None
    
    for line in in_file:
        if line.startswith("Chr\tStart"):
            split_line = line.split("\t")
            lod_index = split_line.index(conservation_column)
        else:
            split_line = line.split("\t")
            lod_column = split_line[lod_index]
            
            if not lod_column == "":
                lod_score = lod_column.split(";")[1]
                lod_score = lod_score.split("=")[2]
                if int(lod_score) >= score_cutoff:
                    out_file.write(line)
    in_file.close()
    out_file.close()

def extract_novel_variants(input, output):
    in_file = open(input, "r")
    out_file = open(output, "w")
    
    for line in in_file:
        split_line = line.split("\t")
        
        if line.startswith("Chr\tStart"):
            esp_index = split_line.index("esp6500si_all")
            thousand_genomes_index = split_line.index("1000g2012apr_all")
            snp138_index = split_line.index("snp138")
            out_file.write(line)
        elif not line.startswith("#"):
            esp_annot = split_line[esp_index] == ""
            thousand_genomes_annot = split_line[thousand_genomes_index] == ""
            snp138_annot = split_line[snp138_index] == ""
            
            if all([esp_annot, thousand_genomes_annot, snp138_annot]):
                out_file.write(line)
    in_file.close()
    out_file.close()

def extract_genic_variants(input, output):
    in_file = open(input, "r")
    out_file = open(output, "w")
    
    re_pattern = re.compile("exonic|splicing|UTR5|URT3|intronic")
    for line in in_file:
        if line.startswith("Chr\tStart"):
            out_file.write(line)
            
        elif not line.startswith("#"):
            if re_pattern.search(line):
                out_file.write(line)
    
    in_file.close()
    out_file.close()

def extract_deleterious_variants(input, output):
    in_file = open(input, "r")
    out_file = open(output, "w")
    
    re_frameshift = re.compile("frameshift insertion | frameshift deletion | frameshift block substitution")
    re_stopalter = re.compile("stopgain | stoploss")
    re_nonsyn = re.compile("nonsynonymous")
    re_splice = re.compile("splicing")
    
    # Parse nonsynonymous SNVs
    for line in in_file:
        split_line = line.split("\t")
        
        if line.startswith("Chr\tStart"):
            sift_pred_index = split_line.index("LJB23_SIFT_pred")
            polyphen2_hdiv_pred_index = split_line.index("LJB23_Polyphen2_HDIV_pred")
            lrt_pred_index = split_line.index("LJB23_LRT_pred")
            mutationtaster_pred_index = split_line.index("LJB23_MutationTaster_pred")
            out_file.write(line)
        elif re_nonsyn.search(line):
            # Test for deleteriousness
            sift = split_line[sift_pred_index] == "T"
            polyphen2 = (split_line[polyphen2_hdiv_pred_index] == "D") or (split_line[polyphen2_hdiv_pred_index] == "P")
            lrt = split_line[lrt_pred_index] == "D"
            mutationtaster = (split_line[mutationtaster_pred_index] == "A") or (split_line[mutationtaster_pred_index] == "D")
            
            if (sum([sift, polyphen2, lrt, mutationtaster]) > 1):
                out_file.write(line)
        elif re_frameshift.search(line):
            out_file.write(line)
        elif re_stopalter.search(line):
            out_file.write(line)
        elif re_splice.search(line):
            out_file.write(line)
            
    in_file.close()
    out_file.close()

def extract_rare_variants(input, output, maf_cutoff):
    in_file = open(input, "r")
    out_file = open(output, "w")
    
    for line in in_file:
        split_line = line.split("\t")
        
        if line.startswith("Chr\tStart"):
            allelefreq_index = split_line.index("1000g2012apr_all")
            out_file.write(line)
        elif not line.startswith("#"):
            if split_line[allelefreq_index] == "":
                out_file.write(line)
            elif float(split_line[allelefreq_index]) < float(maf_cutoff):
                out_file.write(line)
    
    in_file.close()
    out_file.close()

def split_indel_and_snv(input, output_prefix):
    in_file = open(input, "r")
    snv_out_file = open(output_prefix + ".snv.txt", "w")
    indel_out_file = open(output_prefix + ".indel.txt", "w")

    for line in in_file:
        line_split = line.split("\t")
        if line.startswith("Chr\tStart"):
            snv_out_file.write(line)
            indel_out_file.write(line)
            ref_index = line_split.index("Ref")
            alt_index = line_split.index("Alt")
        elif line_split[ref_index] == "-" or line_split[alt_index] == "-" or len(line_split[ref_index]) > 1 or len(line_split[alt_index]) > 1:
            indel_out_file.write(line)
        else:
            snv_out_file.write(line)
    in_file.close()
    snv_out_file.close()
    indel_out_file.close()

rule all:
    input:
        "results/variant_counts/wgs_phased_all_variants.hg19_multianno.genic.lod100.tbx1_counts.txt",
        "results/variant_counts/wgs_phased_all_variants.hg19_multianno.lod100.tbx1_counts.txt",
        "results/novel_variants/novel_deleterious.snv.txt",
        "results/novel_variants/novel_deleterious.indel.txt",
        "results/rare_variants/rare_genic_deleterious.snv.txt",
        "results/rare_variants/rare_genic_deleterious.indel.txt"

rule split_indel_snv_novel:
    input: "results/novel_variants/novel_deleterious_variants.txt"
    params: out_prefix = "results/novel_variants/novel_deleterious"
    output:
        "results/novel_variants/novel_deleterious.snv.txt",
        "results/novel_variants/novel_deleterious.indel.txt"
    run:
        split_indel_and_snv(input[0], params.out_prefix)

rule split_indel_snv_rare:
    input: "results/rare_variants/rare_genic_deleterious_variants.txt"
    params: out_prefix = "results/rare_variants/rare_genic_deleterious"
    output: "results/rare_variants/rare_genic_deleterious.snv.txt",
        "results/rare_variants/rare_genic_deleterious.indel.txt"
    run:
        split_indel_and_snv(input[0], params.out_prefix)

rule extract_rare_deleterious:
    input:
        "results/rare_variants/rare_genic_variants.txt"
    output:
        "results/rare_variants/rare_genic_deleterious_variants.txt"
    run:
        extract_deleterious_variants(input[0], output[0])

rule extract_rare_genic_variants:
    input: "results/rare_variants/rare_variants.txt"
    output: "results/rare_variants/rare_genic_variants.txt"
    run:
        extract_genic_variants(input[0], output[0])

rule extract_rare_variants:
    input: "results/annovar_output/multiannotation/wgs_phased_all_variants.hg19_multianno.txt"
    params: maf_cutoff = "0.01"
    output: "results/rare_variants/rare_variants.txt"
    message: "Extracting variants with MAF < 0.01 in 1000 genomes"
    run:
        extract_rare_variants(input[0], output[0], params.maf_cutoff)

rule extract_novel_deleterious:
    input:
        "results/novel_variants/novel_genic_variants.txt"
    output:
        "results/novel_variants/novel_deleterious_variants.txt"
    run:
        extract_deleterious_variants(input[0], output[0])

rule extract_novel_genic:
    input:
        "results/novel_variants/novel_variants.txt"
    output:
        "results/novel_variants/novel_genic_variants.txt"
    run:
        extract_genic_variants(input[0], output[0])

rule extract_novel_variants:
    input: "results/annovar_output/multiannotation/wgs_phased_all_variants.hg19_multianno.txt"
    output: "results/novel_variants/novel_variants.txt"
    run:
        extract_novel_variants(input[0], output[0])

rule count_variants_tbx1_pathway_genic_lod100:
    input:
        "results/annovar_output/genic_lod100/wgs_phased_all_variants.hg19_multianno.genic.lod100.txt",
        "data/tbx1_pathway_genes_2014_03_27.txt"
    output:
        "results/variant_counts/wgs_phased_all_variants.hg19_multianno.genic.lod100.tbx1_counts.txt"
    shell: """
    perl src/perl/count_var_per_gene.pl \
    {input[1]} \
    {input[0]} \
    {output}
    """

rule count_variants_tbx1_pathway_all_lod100:
    input:
        "results/annovar_output/all_lod100/wgs_phased_all_variants.hg19_multianno.lod100.txt",
        "data/tbx1_pathway_genes_2014_03_27.txt"
    output:
        "results/variant_counts/wgs_phased_all_variants.hg19_multianno.lod100.tbx1_counts.txt"
    shell: """
    perl src/perl/count_var_per_gene.pl \
    {input[1]} \
    {input[0]} \
    {output}
    """

rule all_lod100:
    input: "results/annovar_output/multiannotation/wgs_phased_all_variants.hg19_multianno.txt"
    output: "results/annovar_output/all_lod100/wgs_phased_all_variants.hg19_multianno.lod100.txt"
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
    input: "results/annovar_output/genic/wgs_phased_all_variants.hg19_multianno.genic.txt"
    params:
        conservation_column = "phastConsElements46way",
        lod_cutoff = "100"
    output: "results/annovar_output/genic_lod100/wgs_phased_all_variants.hg19_multianno.genic.lod100.txt"
    run:
        extract_conservation(input[0], 
                             output[0], 
                             params.conservation_column, 
                             int(params.lod_cutoff))


rule extract_genic:
    input: "results/annovar_output/multiannotation/wgs_phased_all_variants.hg19_multianno.txt"
    output: "results/annovar_output/genic/wgs_phased_all_variants.hg19_multianno.genic.txt"
    run:
        extract_genic_variants(input[0], output[0])

rule table_annotation_rare:
    input: "data/for_annovar/wgs_phased_all_variants.avinput"
    params: output_prefix = "results/annovar_output/multiannotation_rare/wgs_phased_all_variants"
    output: "results/annovar_output/multiannotation_rare/wgs_phased_all_variants.hg19_multianno.txt"
    shell: """
    {annovar_dir}/table_annovar.pl \
    {input} \
    {annovar_dir}/humandb/ \
    --protocol refGene,phastConsElements46way,genomicSuperDups,esp6500si_all,1000g2012apr_all,snp138NonFlagged,ljb23_all \
    --operation g,r,r,f,f,f,f \
    --buildver hg19 \
    --outfile {params.output_prefix} \
    --otherinfo
    """


    

########################
# Variants for individual subjects
########################

rule test2:
    input: 
        expand("results/annovar_output/variant_count_summary/wgs_phased_{VAR_TYPE}.sample.{SAMPLEID}.summary.txt", VAR_TYPE = var_type, SAMPLEID = sampleid)

rule summarize_indiv_variant_counts:
    input: 
        multianno = "results/annovar_output/indiv_multiannotation/wgs_phased_{VAR_TYPE}.sample.{SAMPLEID}.hg19_multianno.txt",
        dbsnp = "results/annovar_output/indiv_multiannotation/wgs_phased_{VAR_TYPE}.sample.{SAMPLEID}.hg19_snp138_dropped",
        onekg = "results/annovar_output/indiv_multiannotation/wgs_phased_{VAR_TYPE}.sample.{SAMPLEID}.hg19_ALL.sites.2012_04_dropped",
        exonic = "results/annovar_output/indiv_multiannotation/wgs_phased_{VAR_TYPE}.sample.{SAMPLEID}.refGene.exonic_variant_function"
    output: "results/annovar_output/variant_count_summary/wgs_phased_{VAR_TYPE}.sample.{SAMPLEID}.summary.txt"
    shell: """
        # total variants
        grep -v -P "Chr\tStart" {input.multianno} | grep -v -E "^#" | awk 'END {{print "total\t", NR}}' > {output}
        
        # Count of variants found in dbSNP
        awk 'END {{print "In dbSNP138\t", NR}}' {input.dbsnp} >> {output}
        
        # Count of variants found in 1kg
        awk 'END {{print "In 1kg (2012/10)\t", NR}}' {input.onekg} >> {output}
        
        # Count exonic variants
        awk 'END {{print "Exonic\t", NR}}' {input.exonic} >> {output}
        
        # Count nonsynonymous variants
        grep "nonsynonymous" {input.exonic} | awk 'END {{print "Nonsynonymous\t", NR}}' >> {output}
        
        # Count stop gain/loss
        grep "stop" {input.multianno} | awk 'END {{print "Stop gain/loss\t", NR}}' >> {output}
        
        # Count splice altering
        grep "splicing" {input.multianno} | awk 'END {{print "Splice altering\t", NR}}' >> {output}
        
        # Count frameshift
        grep -w "frameshift" {input.multianno} | awk 'END {{print "Frameshift\t", NR}}' >> {output}
    """
    
# rule extract_subject_specific_annotations:
    # input: 
        # "results/indiv_sample_variants/wgs_phased_{VAR_TYPE}.sample.{SAMPLEID}.variants.txt",
        # "results/annovar_output/multiannotation/wgs_phased_{VAR_TYPE}.hg19_multianno.txt"
    # output:
        # "results/annovar_output/indiv_multiannotation/wgs_phased_{VAR_TYPE}.sample.{SAMPLEID}.hg19_multianno.txt"
    # run:
        # table_anno_file = open(input[1], "r")
        # out_file = open(output[0], "w")
        
        # with open(input[0], "r") as var_file:
            # variants = var_file.read()
        
        # for line in table_anno_file:
            # test_string = line.split("\t")[0:3]
            # test_string = "\t".join(test_string)
            
            # if line.startswith("Chr\tStart") or line.startswith("#"):
                # out_file.write(line)
            # elif any(variant in line for variant in variants):
                # out_file.write(line)

# rule get_subject_specific_variants:
    # input: "results/indiv_avinput/wgs_phased_{VAR_TYPE}.sample.{SAMPLEID}.avinput"
    # output: "results/indiv_sample_variants/wgs_phased_{VAR_TYPE}.sample.{SAMPLEID}.variants.txt"
    # shell: """
        # grep -v -P "Chr\tStart" {input} | \
        # grep -v -E "^#" | \
        # awk '{{print $1,"\t",$2,"\t",$3,"\t",$4,"\t",$5}}' \
        # > {output}
    # """

rule table_annotation_indiv:
    input: "results/indiv_avinput/wgs_phased_{VAR_TYPE}.sample.{SAMPLEID}.avinput"
    params: output_prefix = "results/annovar_output/indiv_multiannotation/wgs_phased_{VAR_TYPE}.sample.{SAMPLEID}"
    output: 
        "results/annovar_output/indiv_multiannotation/wgs_phased_{VAR_TYPE}.sample.{SAMPLEID}.hg19_multianno.txt",
        "results/annovar_output/indiv_multiannotation/wgs_phased_{VAR_TYPE}.sample.{SAMPLEID}.hg19_snp138_dropped",
        "results/annovar_output/indiv_multiannotation/wgs_phased_{VAR_TYPE}.sample.{SAMPLEID}.hg19_ALL.sites.2012_04_dropped",
        "results/annovar_output/indiv_multiannotation/wgs_phased_{VAR_TYPE}.sample.{SAMPLEID}.refGene.exonic_variant_function"
    shell: """
    {annovar_dir}/table_annovar.pl \
    {input} \
    {annovar_dir}/humandb/ \
    --protocol refGene,phastConsElements46way,genomicSuperDups,esp6500si_all,1000g2012apr_all,snp138,ljb23_all \
    --operation g,r,r,f,f,f,f \
    --buildver hg19 \
    --outfile {params.output_prefix} \
    --otherinfo
    """

rule convert_to_annovar_indiv:
    input: "data/indiv_samples/wgs_phased_{VAR_TYPE}.sample.{SAMPLEID}.recode.vcf"
    output: "results/indiv_avinput/wgs_phased_{VAR_TYPE}.sample.{SAMPLEID}.avinput"
    shell: """
        {annovar_dir}/convert2annovar.pl \
        -format vcf4 {input} \
        -outfile {output} \
        -includeinfo \
        -filter PASS \
        -comment
    """
    
rule split_vcf_file:
    input: "data/wgs_phased_{VAR_TYPE}.vcf"
    params: 
        subject = "{SAMPLEID}",
        output_prefix = "data/indiv_samples/wgs_phased_{VAR_TYPE}.sample.{SAMPLEID}"
    output: "data/indiv_samples/wgs_phased_{VAR_TYPE}.sample.{SAMPLEID}.recode.vcf"
    shell: """
        {vcftools} \
        --vcf {input} \
        --indv {params.subject} \
        --out {params.output_prefix} \
        --recode \
        --recode-INFO-all \
        
    """

########################
# All sample table annotations
########################

rule test:
    input: expand("results/annovar_output/variant_count_summary/wgs_phased_{VAR_TYPE}.summary.txt", VAR_TYPE = var_type)
    
rule summarize_variant_counts:
    input: 
        multianno = "results/annovar_output/multiannotation/wgs_phased_{VAR_TYPE}.hg19_multianno.txt",
        dbsnp = "results/annovar_output/multiannotation/wgs_phased_{VAR_TYPE}.hg19_snp138_dropped",
        onekg = "results/annovar_output/multiannotation/wgs_phased_{VAR_TYPE}.hg19_ALL.sites.2012_04_dropped",
        exonic = "results/annovar_output/multiannotation/wgs_phased_{VAR_TYPE}.refGene.exonic_variant_function"
    output: "results/annovar_output/variant_count_summary/wgs_phased_{VAR_TYPE}.summary.txt"
    shell: """
        # total variants
        grep -v -P "Chr\tStart" {input.multianno} | grep -v -E "^#" | awk 'END {{print "total\t", NR}}' > {output}
        
        # Count of variants found in dbSNP
        awk 'END {{print "In dbSNP138\t", NR}}' {input.dbsnp} >> {output}
        
        # Count of variants found in 1kg
        awk 'END {{print "In 1kg (2012/10)\t", NR}}' {input.onekg} >> {output}
        
        # Count exonic variants
        awk 'END {{print "Exonic\t", NR}}' {input.exonic} >> {output}
        
        # Count nonsynonymous variants
        grep "nonsynonymous" {input.exonic} | awk 'END {{print "Nonsynonymous\t", NR}}' >> {output}
        
        # Count stop gain/loss
        grep "stop" {input.multianno} | awk 'END {{print "Stop gain/loss\t", NR}}' >> {output}
        
        # Count splice altering
        grep "splicing" {input.multianno} | awk 'END {{print "Splice altering\t", NR}}' >> {output}
        
        # Count frameshift
        grep -w "frameshift" {input.multianno} | awk 'END {{print "Frameshift\t", NR}}' >> {output}
    """

rule table_annotation_novel:
    input: "data/for_annovar/wgs_phased_{VAR_TYPE}.avinput"
    params: output_prefix = "results/annovar_output/multiannotation/wgs_phased_{VAR_TYPE}"
    output:
        "results/annovar_output/multiannotation/wgs_phased_{VAR_TYPE}.hg19_multianno.txt",
        "results/annovar_output/multiannotation/wgs_phased_{VAR_TYPE}.hg19_snp138_dropped",
        "results/annovar_output/multiannotation/wgs_phased_{VAR_TYPE}.hg19_ALL.sites.2012_04_dropped",
        "results/annovar_output/multiannotation/wgs_phased_{VAR_TYPE}.refGene.exonic_variant_function"
    output: "results/annovar_output/variant_count_summary/wgs_phased_{VAR_TYPE}.summary.txt"
    shell: """
    {annovar_dir}/table_annovar.pl \
    {input} \
    {annovar_dir}/humandb/ \
    --protocol refGene,phastConsElements46way,genomicSuperDups,esp6500si_all,1000g2012apr_all,snp138,ljb23_all \
    --operation g,r,r,f,f,f,f \
    --buildver hg19 \
    --outfile {params.output_prefix} \
    --otherinfo
    """

rule convertToAnnovar_allsample:
    input: "data/wgs_phased_{VAR_TYPE}.vcf"
    output: "data/for_annovar/wgs_phased_{VAR_TYPE}.avinput"
    message: "Converting {input} into format for annovar. Use '--format vcf4old` to keep all variants"
    shell:"""
        {annovar_dir}/convert2annovar.pl \
        -format vcf4old {input} \
        -outfile {output} \
        -includeinfo \
        -filter PASS \
        -allallele \
        -comment
    """
