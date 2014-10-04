import re
from snakemake.utils import R
from subprocess import call
import os

workdir: "/home/jchung/projects/wgs/"
# Path to required programs
annovar_dir = '/home/jchung/programs/annovar/annovar2013Aug23/'
vcftools = "/apps1/vcftools/vcftools_0.1.11/cpp/vcftools"
vcftools_vcfisec = "/apps1/vcftools/vcftools_0.1.11/perl/vcf-isec"

# Wildcards
variantFiles = ["wgs.k1_only.maf1.deleterious.variant.txt",
                "wgs.k2_only.maf1.deleterious.variant.txt",
                "wgs.maf1.deleterious.variant.txt"]
vcfFiles     = ["wgs_phased_indel.vcf", "wgs_phased_snp.vcf"]
sampleid     = ["BM1452.001", "BM1452.100", "BM1452.200", "BM1453.001", 
                "BM1453.100", "BM1453.200"]
probandid    = ["BM1452.001", "BM1453.001"]
var_type     = ["snp", "indel"]
var_freq     = ["maf_0.01", "maf_0.05", "novel", "all"]



################################################################################
# Functions
################################################################################
def remove_superdup(input, output):
    """Remove variants that are found in genomic superdup regions.
    
    Args:
        input: Path to the input file.
        output: Path to the output file.
    Returns:
        Creates an output file without variants in genomic superdups.
    """
    in_file = open(input, "r")
    out_file = open(output, "w")

    for line in in_file:
        if line.startswith("Chr\tStart"):
            split_line = line.split("\t")
            superdup_index = split_line.index("genomicSuperDups")
            out_file.write(line)
        elif line.startswith("#"):
            out_file.write(line)
        else:
            split_line = line.split("\t")
            if split_line[superdup_index] == "":
                out_file.write(line)
    in_file.close()
    out_file.close()

def extract_conserved_elements(input, output, conservation_column, score_cutoff):
    in_file = open(input, "r")
    out_file = open(output, "w")
    score_index = None
    
    for line in in_file:
        if line.startswith("Chr\tStart"):
            split_line = line.split("\t")
            score_index = split_line.index(conservation_column)
            out_file.write(line)
        elif line.startswith("#"):
            out_file.write(line)
        else:
            split_line = line.split("\t")
            score_column = split_line[score_index]
            
            if not score_column == "":
                score = score_column.split(";")[1]
                score = score.split("=")[2]
                if int(score) >= int(score_cutoff):
                    out_file.write(line)
    in_file.close()
    out_file.close()

def parse_coding_vs_noncoding(input, output_exonic, output_noncoding):
    """Parse exonic and noncoding variants from ANNOVAR results.
    
    Split input file into exonic and noncoding variants.
    
    Exonic variants are defined by annovar as:
        
        * exonic
        * splicing
        * UTR5
        * UTR3
    
    Non-coding variants are defined by annovar as:
        
        * ncRNA
        * intronic
        * upstream
        * downstream
        * intergenic
    
    Args:
        input: Path to input file.
        output: Path to output file.
    Returns:
        Creates an output file with only noncoding variants.
    """
    in_file = open(input, "r")
    exonic_out_file = open(output_exonic, "w")
    noncoding_out_file = open(output_noncoding, "w")
    
    re_pattern = re.compile("\sexonic\s|\ssplicing\s|\sUTR5\s|\sURT3\s")
    for line in in_file:
        if line.startswith("Chr\tStart"):
            exonic_out_file.write(line)
            noncoding_out_file.write(line)
        elif line.startswith("#"):
            exonic_out_file.write(line)
            noncoding_out_file.write(line)
        elif not line.startswith("#"):
            if re_pattern.search(line):
                exonic_out_file.write(line)
            else:
                noncoding_out_file.write(line)
    in_file.close()
    exonic_out_file.close()
    noncoding_out_file.close()

def extract_nonsynonymous_variants(input, output):
    """Extract nonsynonymous variants from the input file.
    
    Search variant annotation for:
    
        * splicing
        * nonsynonymous SNV
        * stopgain
        * stoploss
        * frameshift insertion
        * frameshift deletion
        * frameshift block substitution
    
    Args:
        input: Path to input file
        output: Path to output file
    Returns:
        Creates an output file with nonsynonymous variants
    """
    in_file = open(input, "r")
    out_file = open(output, "w")
    
    re_nonsynonymous = re.compile(
        "|".join(["\snonsynonymous SNV\s", "\sstopgain\s", "\sstoploss\s", 
                    "\ssplicing\s", "\sframeshift insertion\s", 
                    "\sframeshift deletion\s", 
                    "\sframeshift block substitution\s"])
    )
    
    for line in in_file:
        if line.startswith("Chr\tStart"):
            out_file.write(line)
        elif line.startswith("#"):
            out_file.write(line)
        elif re_nonsynonymous.search(line):
            out_file.write(line)
    in_file.close()
    out_file.close()

def extract_deleterious_variants(input, output, deleterious_score_cutoff):
    """Extract deleterious variants.
    
    Extract variants that are predicted to be deleterious.
    
    Nonsynonymous variants are selected if they are predicted to be deleterious
    by at least the number of methods defined by `deleterious_score_cutoff'.
    
    Frameshift, stop altering and splicing altering variants are automatically
    predicted to be deleterious.
    
    Args:
        input: Path to input file.
        output: Path to output file.
        deleterious_score_cutoff: Integer indicating the number of methods 
            needed to define a nonsynonymous variant as deleterious.
    Return:
        Creates an output file with deleterious variants.
    """
    deleterious_score_cutoff = int(deleterious_score_cutoff)
    in_file = open(input, "r")
    out_file = open(output, "w")
    
    re_frameshift = re.compile("\sframeshift insertion\s|\sframeshift deletion\s|\sframeshift block substitution\s")
    re_stopalter = re.compile("\sstopgain\s|\sstoploss\s")
    re_nonsyn = re.compile("nonsynonymous")
    re_splice = re.compile("splicing")
    
    
    for line in in_file:
        split_line = line.split("\t")
        
        if line.startswith("Chr\tStart"):
            sift_pred_index = split_line.index("LJB23_SIFT_pred")
            polyphen2_hdiv_pred_index = split_line.index("LJB23_Polyphen2_HDIV_pred")
            lrt_pred_index = split_line.index("LJB23_LRT_pred")
            mutationtaster_pred_index = split_line.index("LJB23_MutationTaster_pred")
            out_file.write(line)
        elif line.startswith("#"):
            out_file.write(line)
        elif re_nonsyn.search(line):
            # Test for deleteriousness of nonsynonymous variants
            sift = split_line[sift_pred_index] == "D"
            polyphen2 = (split_line[polyphen2_hdiv_pred_index] == "D") or (split_line[polyphen2_hdiv_pred_index] == "P")
            lrt = split_line[lrt_pred_index] == "D"
            mutationtaster = (split_line[mutationtaster_pred_index] == "A") or (split_line[mutationtaster_pred_index] == "D")
            
            deleterious_score = sum([sift, polyphen2, lrt, mutationtaster])
            if (deleterious_score >= deleterious_score_cutoff):
                out_file.write(line)
        elif re_frameshift.search(line):
            out_file.write(line)
        elif re_stopalter.search(line):
            out_file.write(line)
        elif re_splice.search(line):
            out_file.write(line)
            
    in_file.close()
    out_file.close()

def extract_novel_variants(
    input, 
    output, 
    dbsnp_database = "snp138", 
    esp_database = "esp6500si_all", 
    thousand_genomes_database = "1000g2012apr_all"):
    """Remove variants from population databases.
    
    Args:
        input: Path to input file.
        output: Path to output file.
        dbsnp_database: Name of the dbsnp_database to filter on.
        esp_database: Name of the exome sequencing project databse to filter on.
        thousand_genomes_database: Name of the 1000 Genomes database to filter on.
    Returns:
        Creates an output file with all 
    """
    in_file = open(input, "r")
    out_file = open(output, "w")
    
    for line in in_file:
        split_line = line.split("\t")
        
        if line.startswith("Chr\tStart"):
            dbsnp_index = split_line.index(dbsnp_database)
            esp_index = split_line.index(esp_database)
            thousand_genomes_index = split_line.index(thousand_genomes_database)
            out_file.write(line)
        elif line.startswith("#"):
            out_file.write(line)
        elif not line.startswith("#"):
            is_dbsnp_novel = split_line[dbsnp_index] == ""
            is_esp_novel = split_line[esp_index] == ""
            is_1kg_novel = split_line[thousand_genomes_index] == ""
            
            if all([is_dbsnp_novel, is_esp_novel, is_1kg_novel]):
                out_file.write(line)
    in_file.close()
    out_file.close()

def extract_rare_variants(input, output, maf_cutoff):
    """Extract rare variants based on 1000 Genomes variant frequencies.
    
    Args:
        input: Path to input file.
        output: Path to output file.
        maf_cutoff: Minor allele frequency cutoff
    """
    in_file = open(input, "r")
    out_file = open(output, "w")
    
    for line in in_file:
        split_line = line.split("\t")
        
        if line.startswith("Chr\tStart"):
            allelefreq_index = split_line.index("1000g2012apr_all")
            out_file.write(line)
        elif line.startswith("#"):
            out_file.write(line)
        elif not line.startswith("#"):
            if split_line[allelefreq_index] == "":
                out_file.write(line)
            elif float(split_line[allelefreq_index]) < float(maf_cutoff):
                out_file.write(line)
    
    in_file.close()
    out_file.close()

def check_alt_in_genotype(genotype):
    """Determine if the genotype has an alternative allele.
    
    Args:
        genotype: genotype in VCF format.
    Returns:
        boolean indicating whether an alternative allele is present. None if the
        genotype is missing.
    """
    genotype = genotype.split(":")[0]
    if "." in genotype:
        return None
    elif "1" in genotype:
        return True
    else:
        return False

def check_alt_allele_status(genotype):
    """Get the allele status of a genotype.
    
    Args:
        genotype: genotype in VCF format.
    Returns:
        string indicating the type allele status. `heterozygous' if one minor
        allele. `homozygous' if two minor alleles. None if missing.
    """
    genotype = genotype.split(":")[0]
    if "." in genotype:
        return None
    elif genotype.count("1") == 1:
        return "heterozygous"
    elif genotype.count("1") == 2:
        return "homozygous"

def check_inherited_variant(proband_genotype, father_genotype, mother_genotype, 
    inherited_from = "father"):
    """Check which parent a variant is inherited from.
    
    Args:
        proband_genotype: genotype of proband in VCF format.
        father_genotype: genotype of father in VCF format.
        mother_genotype: genotype of mother in VCF format.
        inherited_from: string indicating which parent to check.
    Returns:
        boolean indicating if inheritance makes sense.
    """
    proband_alt = check_alt_in_genotype(proband_genotype)
    father_alt = check_alt_in_genotype(father_genotype)
    mother_alt = check_alt_in_genotype(mother_genotype)
    if check_alt_allele_status(proband_genotype) == "homozygous":
        # If proband is homozygous alternate, one allele should be inherited
        # from each parent. Check to make sure this is the case.
        meet_criteria = all([father_alt, mother_alt])
    else:
        if inherited_from.lower() == "father":
            # If proband has an alternate allele, check if it found in the
            # father and not the mother.
            meet_criteria = all([proband_alt, father_alt, not mother_alt])
        elif inherited_from.lower() == "mother":
            # If proband has an alternate allele, check if it found in the
            # mother and not the father.
            meet_criteria = all([proband_alt, not father_alt, mother_alt])
    return meet_criteria

def extract_inherited_variant(input, output,proband_id, father_id, mother_id, 
    inherited_from):
    """
    """
    in_file = open(input, "r")
    out_file = open(output, "w")
    
    for line in in_file:
        if line.startswith("##"):
            out_file.write(line)
        elif line.startswith("#CHROM"):
            header = line.split()
            proband_index = header.index(proband_id)
            father_index = header.index(father_id)
            mother_index = header.index(mother_id)
            out_file.write(line)
        elif not line.startswith("#"):
            split_line = line.split()
            proband_genotype = split_line[proband_index]
            father_genotype = split_line[father_index]
            mother_genotype = split_line[mother_index]
            if "." not in proband_genotype:
                if check_inherited_variant(proband_genotype, father_genotype, mother_genotype, inherited_from):
                    out_file.write(line)
    out_file.close()
    in_file.close()

def check_alt_allele_ratio(
    genotype, allele_counts,hom_ref_cutoff = 0.15, het_low_cutoff = 0.3, 
    het_high_cutoff = 0.7, hom_alt_cutoff = 0.85):
    """Calculate alternate allele ratio and determine if it falls within QC
    measures.
    
    Args:
        
    """
    allele_counts = allele_counts.split(",")
    ref_count = float(allele_counts[0])
    alt_count = float(allele_counts[1])
    total_allele_count = ref_count + alt_count
    
    if total_allele_count == 0:
        return False
    else:
        num_alt_alleles = genotype.count("1")
        alt_allele_ratio = alt_count / (total_allele_count)
        
        if num_alt_alleles == 0:
            if alt_allele_ratio > 0.15:
                return False
            else:
                return True
        elif num_alt_alleles == 1:
            if (alt_allele_ratio < 0.3) or (alt_allele_ratio > 0.7):
                return False
            else:
                return True
        elif num_alt_alleles == 2:
            if alt_allele_ratio < 0.85:
                return False
            else:
                return True

def check_genotype_quality(quality_score):
    if int(quality_score) < 20:
        return False
    else:
        return True

def check_read_depth(read_depth):
    if int(read_depth) < 15:
        return False
    else:
        return True

def genotype_quality_control(genotype_info):
    """Filter genotypes.
    
    Filter genotypes based on:
        * Alternate allele ratio
        * Read depth
        * Genotype quality
    
    Args:
        genotype_info: Genotype to test.
    Returns:
        boolean indicating whether the genotype passes all quality control tests.
    """
    genotype_split = genotype_info.split(":")
    if "." in genotype_split[0]:
        return None
    else:
        pass_alt_allele_ratio = check_alt_allele_ratio(genotype_split[0], genotype_split[1])
        pass_depth = check_read_depth(genotype_split[2])
        pass_genotype_quality = check_genotype_quality(genotype_split[3])
        pass_qc = all([pass_alt_allele_ratio, pass_depth, pass_genotype_quality])
        return pass_qc

def extract_denovo_variants(input, output, proband_id, father_id, mother_id):
    in_file = open(input, "r")
    out_file = open(output, "w")
    
    for line in in_file:
        if line.startswith("##"):
            out_file.write(line)
        elif line.startswith("#CHROM"):
            header = line.split()
            proband_index = header.index(proband_id)
            father_index = header.index(father_id)
            mother_index = header.index(mother_id)
            out_file.write(line)
        elif not line.startswith("#"):
            split_line = line.split()
            proband_genotype = split_line[proband_index]
            father_genotype = split_line[father_index]
            mother_genotype = split_line[mother_index]
            
            # print(proband_genotype)
            if "." not in proband_genotype:
                proband_alt = check_alt_in_genotype(proband_genotype)
                father_alt = check_alt_in_genotype(father_genotype)
                mother_alt = check_alt_in_genotype(mother_genotype)
                
                if not any([proband_alt == None, father_alt == None, mother_alt == None]):
                    # Check genotype quality
                    proband_qc = genotype_quality_control(proband_genotype)
                    father_qc = genotype_quality_control(father_genotype)
                    mother_qc = genotype_quality_control(mother_genotype)
                
                    # print(proband_genotype)
                    # print(proband_alt_allele_ratio)
                    if all([proband_alt,
                        father_alt == False,
                        mother_alt == False,
                        proband_qc,
                        father_qc,
                        mother_qc]):
                        
                        out_file.write(line)
    out_file.close()
    in_file.close()

def summarize_exonic_functions(input, output):
    """Summarize exonic functions by counting up the number of variants with a 
    given function.
    
    Args:
        input: Path to input file.
        output: Path to output file.
    """
    exonic_function_column = "ExonicFunc.refGene"

    all_exonic_function = []
    
    with open(input) as in_file:
        for line in in_file:
            if line.startswith("Chr\tStart"):
                header = line.split("\t")
                exonic_function_index = header.index(exonic_function_column)
            else:
                var_annotation = line.split("\t")
                var_exonic_function = var_annotation[exonic_function_index]
                all_exonic_function.append(var_exonic_function)
    
    frameshift_insertion = all_exonic_function.count("frameshift insertion")
    frameshift_deletion = all_exonic_function.count("frameshift deletion")
    frameshift_block_sub = all_exonic_function.count("frameshift block substitution")
    stopgain = all_exonic_function.count("stopgain")
    stoploss = all_exonic_function.count("stoploss")
    nonframeshift_insertion = all_exonic_function.count("nonframeshift insertion")
    nonframeshift_deletion = all_exonic_function.count("nonframeshift deletion")
    nonframeshift_block_sub = all_exonic_function.count("nonframeshift block substitution")
    nonsynonymous_snv = all_exonic_function.count("nonsynonymous SNV")
    synonymous_snv = all_exonic_function.count("synonymous SNV")
    unknown = all_exonic_function.count("unknown")
    
    with open(output, "w") as out_file:
        out_file.write("Exonic function\tCount\n")
        out_file.write("Frameshift insertion\t{0}\n".format(frameshift_insertion))
        out_file.write("Frameshift deletion\t{0}\n".format(frameshift_deletion))
        out_file.write("Frameshift block substitution\t{0}\n".format(frameshift_block_sub))
        out_file.write("Stopgain\t{0}\n".format(stopgain))
        out_file.write("Stoploss\t{0}\n".format(stoploss))
        out_file.write("Nonframeshift insertion\t{0}\n".format(nonframeshift_insertion))
        out_file.write("Nonframeshift deletion\t{0}\n".format(nonframeshift_deletion))
        out_file.write("Nonframeshift block substitution\t{0}\n".format(nonframeshift_block_sub))
        out_file.write("Nonsynonymous SNV\t{0}\n".format(nonsynonymous_snv))
        out_file.write("Synonymous SNV\t{0}\n".format(synonymous_snv))
        out_file.write("Unknown\t{0}\n".format(unknown))

def summarize_variant_function(input, output):
    """Count variants by variant function.
    
    Args:
        input: Path to input file.
        output: Path to output file.
    """
    variant_function_column = "Func.refGene"
    all_variant_function = []
    with open(input) as in_file:
        for line in in_file:
            if line.startswith("Chr\tStart"):
                header = line.split("\t")
                variant_function_column = header.index(variant_function_column)
            else:
                var_annotation = line.split("\t")
                var_function = var_annotation[variant_function_column]
                all_variant_function.append(var_function)
    
    exonic = all_variant_function.count("exonic")
    splicing = all_variant_function.count("splicing")
    ncrna_exonic = all_variant_function.count("ncRNA_exonic")
    ncrna_intronic = all_variant_function.count("ncRNA_intronic")
    utr5 = all_variant_function.count("UTR5")
    utr3 = all_variant_function.count("UTR3")
    intronic = all_variant_function.count("intronic")
    upstream = all_variant_function.count("upstream")
    downstream = all_variant_function.count("downstream")
    intergenic = all_variant_function.count("intergenic")
    
    with open(output, "w") as out_file:
        out_file.write("Variant function\tCount\n")
        out_file.write("Exonic\t{0}\n".format(exonic))
        out_file.write("splicing\t{0}\n".format(splicing))
        out_file.write("ncRNA exonic\t{0}\n".format(ncrna_exonic))
        out_file.write("ncRNA intronic\t{0}\n".format(ncrna_intronic))
        out_file.write("5'-UTR\t{0}\n".format(utr5))
        out_file.write("3'-UTR\t{0}\n".format(utr3))
        out_file.write("Intronic\t{0}\n".format(intronic))
        out_file.write("Upstream\t{0}\n".format(upstream))
        out_file.write("Downstream\t{0}\n".format(downstream))
        out_file.write("Intergenic\t{0}\n".format(intergenic))

################################################################################
# Begin Rules
################################################################################
rule all:
    input:
        "results/variant_counts/wgs_phased_all_variants.hg19_multianno.genic.lod100.tbx1_counts.txt",
        "results/variant_counts/wgs_phased_all_variants.hg19_multianno.lod100.tbx1_counts.txt",
        "results/novel_variants/novel_deleterious.snv.txt",
        "results/novel_variants/novel_deleterious.indel.txt",
        "results/rare_variants/rare_genic_deleterious.snv.txt",
        "results/rare_variants/rare_genic_deleterious.indel.txt"

rule split_indel_snv_novel:
    input:
        "results/novel_variants/novel_deleterious_variants.txt"
    params:
        out_prefix = "results/novel_variants/novel_deleterious"
    output:
        "results/novel_variants/novel_deleterious.snv.txt",
        "results/novel_variants/novel_deleterious.indel.txt"
    run:
        split_indel_and_snv(input[0], params.out_prefix)

rule split_indel_snv_rare:
    input: 
        "results/rare_variants/rare_genic_deleterious_variants.txt"
    params:
        out_prefix = "results/rare_variants/rare_genic_deleterious"
    output:
        "results/rare_variants/rare_genic_deleterious.snv.txt",
        "results/rare_variants/rare_genic_deleterious.indel.txt"
    run:
        split_indel_and_snv(input[0], params.out_prefix)

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

rule checkpoint1:
    input:
        expand("results/annovar/output/counts/wgs_phased_{VAR_TYPE}.sample.{SAMPLEID}.remaining22q11.filtering_summary.txt", 
               VAR_TYPE = var_type, SAMPLEID = probandid),
        expand("results/annovar/output/counts/remaining_22q11/wgs_phased_{VAR_TYPE}.sample.{SAMPLEID}.remaining22q11.exonic.nonsyn.deleterious.var_function.txt",
               VAR_TYPE = var_type, SAMPLEID = probandid),
        expand("results/annovar/output/counts/remaining_22q11/wgs_phased_{VAR_TYPE}.sample.{SAMPLEID}.remaining22q11.exonic.nonsyn.deleterious.exonic_function.txt", VAR_TYPE = var_type, SAMPLEID = probandid),
        expand("results/annovar/output/counts/wgs_phased_{VAR_TYPE}.sample.{SAMPLEID}.freq.{FREQ}.filtering_summary.txt", FREQ = var_freq, VAR_TYPE = var_type, SAMPLEID = sampleid),
        expand("results/annovar/output/counts/wgs_phased_{VAR_TYPE}.sample.{SAMPLEID}.freq.{FREQ}.nosuperdup.exonic.deleterious.var_function.txt", FREQ = var_freq, VAR_TYPE = var_type, SAMPLEID = sampleid),
        expand("results/annovar/output/counts/wgs_phased_{VAR_TYPE}.sample.{SAMPLEID}.freq.{FREQ}.nosuperdup.exonic.deleterious.exonic_function.txt", FREQ = var_freq, VAR_TYPE = var_type, SAMPLEID = sampleid),
        "results/annovar/output/counts/wgs_phased_snp.all_samples.pop_counts.txt",
        "results/annovar/output/counts/wgs_phased_indel.all_samples.pop_counts.txt"
################
# 
################

rule denovo_novel_rare_all:
    input:
        "results/annovar_output/counts/denovo.wgs_phased_snp.all_samples.summary.txt",
        "results/annovar_output/counts/denovo.wgs_phased_indel.all_samples.summary.txt",
        expand("results/annovar_output/variant_filtering/novel/snp138.wgs_phased_{VAR_TYPE}.sample.{SAMPLEID}.nosuperdup.novel.genic.deleterious.counts.txt", VAR_TYPE = var_type, SAMPLEID = probandid),
        expand("results/annovar_output/variant_filtering/rare/snp138.wgs_phased_{VAR_TYPE}.sample.{SAMPLEID}.nosuperdup.rare.genic.deleterious.counts.txt", VAR_TYPE = var_type, SAMPLEID = probandid),
        expand("results/annovar_output/variant_filtering/remaining_22q11/wgs_phased_{VAR_TYPE}.sample.{SAMPLEID}.remaining22q11.genic.deleterious.txt", VAR_TYPE = var_type, SAMPLEID = probandid),
        expand("results/annovar_output/variant_filtering/remaining_22q11/wgs_phased_{VAR_TYPE}.sample.{SAMPLEID}.remaining22q11.summary.txt", VAR_TYPE = var_type, SAMPLEID = probandid)

################################################################################
# de novo variant extraction
################################################################################
rule denovo_all:
    input:
        "results/annovar_output/counts/denovo.wgs_phased_snp.all_samples.summary.txt",
        "results/annovar_output/counts/denovo.wgs_phased_indel.all_samples.summary.txt"

rule combine_variant_counts_denovo:
    input:
        expand("results/annovar_output/counts/wgs_phased_{VAR_TYPE}.sample.{SAMPLE}.denovo.txt", 
               VAR_TYPE = var_type, SAMPLE = ["BM1452.001", "BM1453.001"])
    output:
        "results/annovar_output/counts/denovo.wgs_phased_snp.all_samples.summary.txt",
        "results/annovar_output/counts/denovo.wgs_phased_indel.all_samples.summary.txt"
    shell:"""
        Rscript src/annovar_wgs_samples/merge_variant_count_summaries.R \
        {output[0]} \
        {output[1]} \
        {input}
    """

rule summarize_denovo_variant_counts:
    input: 
        multianno = "results/annovar_output/multiannotation/denovo/wgs_phased_{VAR_TYPE}.sample.{SAMPLE}.denovo.hg19_multianno.txt",
        dbsnp = "results/annovar_output/multiannotation/denovo/wgs_phased_{VAR_TYPE}.sample.{SAMPLE}.denovo.hg19_snp138_dropped",
        onekg = "results/annovar_output/multiannotation/denovo/wgs_phased_{VAR_TYPE}.sample.{SAMPLE}.denovo.hg19_ALL.sites.2012_04_dropped",
        exonic = "results/annovar_output/multiannotation/denovo/wgs_phased_{VAR_TYPE}.sample.{SAMPLE}.denovo.refGene.exonic_variant_function",
        variant_function = "results/annovar_output/multiannotation/denovo/wgs_phased_{VAR_TYPE}.sample.{SAMPLE}.denovo.refGene.variant_function"
    output: 
        "results/annovar_output/counts/wgs_phased_{VAR_TYPE}.sample.{SAMPLE}.denovo.txt"
    shell: """
        grep -v -P "Chr\tStart" {input.multianno} \
            | grep -v -E "^#" \
            | awk 'END {{print "total\t", NR}}' \
            > {output}
        awk 'END {{print "In dbSNP138\t", NR}}' {input.dbsnp} \
            >> {output}
        awk 'END {{print "In 1kg (2012/10)\t", NR}}' {input.onekg} \
            >> {output}
        awk 'END {{print "Exonic\t", NR}}' {input.exonic} \
            >> {output}
        grep "nonsynonymous" {input.exonic} \
            | awk 'END {{print "Nonsynonymous\t", NR}}' \
            >> {output}
        grep "stop" {input.multianno} \
            | awk 'END {{print "Stop gain/loss\t", NR}}' \
            >> {output}
        grep "splicing" {input.multianno} \
            | awk 'END {{print "Splice altering\t", NR}}' \
            >> {output}
        grep -w "frameshift" {input.multianno} \
            | awk 'END {{print "Frameshift\t", NR}}' \
            >> {output}
        grep -w "intronic" {input.variant_function} \
            | awk 'END {{print "Intronic\t", NR}}' \
            >> {output}
        grep -w "upstream" {input.variant_function} \
            | awk 'END {{print "Upstream\t", NR}}' \
            >> {output}
        grep -w "downstream" {input.variant_function} \
            | awk 'END {{print "Downstream\t", NR}}' \
            >> {output}
        grep -w "intergenic" {input.variant_function} \
            | awk 'END {{print "Intergenic\t", NR}}' \
            >> {output}
    """

rule table_annotation_denovo:
    input: 
        "data/for_annovar/wgs_phased_{VAR_TYPE}.sample.{SAMPLE}.denovo.recode.avinput"
    params: 
        output_prefix = "results/annovar_output/multiannotation/denovo/wgs_phased_{VAR_TYPE}.sample.{SAMPLE}.denovo"
    output: 
        "results/annovar_output/multiannotation/denovo/wgs_phased_{VAR_TYPE}.sample.{SAMPLE}.denovo.hg19_multianno.txt",
        "results/annovar_output/multiannotation/denovo/wgs_phased_{VAR_TYPE}.sample.{SAMPLE}.denovo.hg19_snp138_dropped",
        "results/annovar_output/multiannotation/denovo/wgs_phased_{VAR_TYPE}.sample.{SAMPLE}.denovo.hg19_ALL.sites.2012_04_dropped",
        "results/annovar_output/multiannotation/denovo/wgs_phased_{VAR_TYPE}.sample.{SAMPLE}.denovo.refGene.exonic_variant_function",
        "results/annovar_output/multiannotation/denovo/wgs_phased_{VAR_TYPE}.sample.{SAMPLE}.denovo.refGene.variant_function"
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

rule convert_denovo_to_avinput:
    input:
        "data/indiv_samples/denovo/wgs_phased_{VAR_TYPE}.sample.{SAMPLE}.denovo.recode.vcf"
    output:
        "data/for_annovar/wgs_phased_{VAR_TYPE}.sample.{SAMPLE}.denovo.recode.avinput"
    shell: """
        {annovar_dir}/convert2annovar.pl \
        -format vcf4 {input} \
        -outfile {output} \
        -includeinfo \
        -filter PASS \
        -comment
    """

rule split_denovo_vcf_file:
    input: "data/indiv_samples/denovo/wgs_phased_{VAR_TYPE}.sample.{SAMPLE}.denovo.vcf"
    params: 
        subject = "{SAMPLE}",
        output_prefix = "data/indiv_samples/denovo/wgs_phased_{VAR_TYPE}.sample.{SAMPLE}.denovo"
    output: "data/indiv_samples/denovo/wgs_phased_{VAR_TYPE}.sample.{SAMPLE}.denovo.recode.vcf"
    shell: """
        {vcftools} \
        --vcf {input} \
        --indv {params.subject} \
        --out {params.output_prefix} \
        --recode \
        --recode-INFO-all
    """

rule extract_denovo_variants_bm1453:
    input:
        "data/wgs_phased_{VAR_TYPE}.vcf"
    params:
        proband = "BM1453.001",
        father = "BM1453.200",
        mother = "BM1453.100"
    output:
        "data/indiv_samples/denovo/wgs_phased_{VAR_TYPE}.sample.BM1453.001.denovo.vcf"
    run:
        extract_denovo_variants(input[0], 
                                output[0], 
                                params.proband, 
                                params.father, 
                                params.mother)

rule extract_denovo_variants_bm1452:
    input:
        "data/wgs_phased_{VAR_TYPE}.vcf"
    params:
        proband = "BM1452.001",
        father = "BM1452.200",
        mother = "BM1452.100"
    output:
        "data/indiv_samples/denovo/wgs_phased_{VAR_TYPE}.sample.BM1452.001.denovo.vcf"
    run:
        extract_denovo_variants(input[0], 
                                output[0], 
                                params.proband, 
                                params.father, 
                                params.mother)

################################################################################
# Rare inherited extraction
################################################################################
rule rare_inherited_all:
    input: expand("results/annovar_output/variant_filtering/rare_inherited/inherited.qc_wgs_phased_{VAR_TYPE}.sample.{SAMPLEID}.nosuperdup.rare.genic.deleterious.txt", VAR_TYPE = var_type, SAMPLEID = "BM1452.001")

rule rare_inherited_extract_deleterious:
    input:
        "results/annovar_output/variant_filtering/rare_inherited/inherited.qc_wgs_phased_{VAR_TYPE}.sample.{SAMPLEID}.nosuperdup.rare.genic.txt"
    params: deleterious_score_cutoff = "2"
    output:
        "results/annovar_output/variant_filtering/rare_inherited/inherited.qc_wgs_phased_{VAR_TYPE}.sample.{SAMPLEID}.nosuperdup.rare.genic.deleterious.txt"
    run:
        extract_deleterious_variants(input[0], output[0], params.deleterious_score_cutoff)

rule rare_inherited_extract_genic:
    input: "results/annovar_output/variant_filtering/rare_inherited/inherited.qc_wgs_phased_{VAR_TYPE}.sample.{SAMPLEID}.nosuperdup.rare.txt"
    output: "results/annovar_output/variant_filtering/rare_inherited/inherited.qc_wgs_phased_{VAR_TYPE}.sample.{SAMPLEID}.nosuperdup.rare.genic.txt"
    run:
        extract_exonic_variants(input[0], output[0])

rule rare_inherited_extract_population:
    input:
        "results/annovar_output/variant_filtering/rare_inherited/inherited.qc_wgs_phased_{VAR_TYPE}.sample.{SAMPLEID}.nosuperdup.txt"
    params: maf_cutoff = "0.05"
    output:
        "results/annovar_output/variant_filtering/rare_inherited/inherited.qc_wgs_phased_{VAR_TYPE}.sample.{SAMPLEID}.nosuperdup.rare.txt"
    run:
        extract_rare_variants(input[0], output[0], params.maf_cutoff)

rule rare_inherited_remove_superdup:
    input: "results/annovar_output/inherited_multiannotation/inherited.qc_wgs_phased_{VAR_TYPE}.sample.{SAMPLEID}.hg19_multianno.txt"
    output: "results/annovar_output/variant_filtering/rare_inherited/inherited.qc_wgs_phased_{VAR_TYPE}.sample.{SAMPLEID}.nosuperdup.txt"
    run:
        remove_superdup(input[0], output[0])

rule rare_inherited_table_annotation_indiv:
    input: "results/indiv_avinput/inherited.qc_wgs_phased_{VAR_TYPE}.sample.{SAMPLEID}.avinput"
    params: output_prefix = "results/annovar_output/inherited_multiannotation/inherited.qc_wgs_phased_{VAR_TYPE}.sample.{SAMPLEID}"
    output: 
        "results/annovar_output/inherited_multiannotation/inherited.qc_wgs_phased_{VAR_TYPE}.sample.{SAMPLEID}.hg19_multianno.txt",
        "results/annovar_output/inherited_multiannotation/inherited.qc_wgs_phased_{VAR_TYPE}.sample.{SAMPLEID}.hg19_snp138_dropped",
        "results/annovar_output/inherited_multiannotation/inherited.qc_wgs_phased_{VAR_TYPE}.sample.{SAMPLEID}.hg19_ALL.sites.2012_04_dropped",
        "results/annovar_output/inherited_multiannotation/inherited.qc_wgs_phased_{VAR_TYPE}.sample.{SAMPLEID}.refGene.exonic_variant_function"
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

rule rare_inherited_convert_to_annovar_indiv:
    input: "data/indiv_samples/inherited.qc_wgs_phased_{VAR_TYPE}.sample.{SAMPLEID}.recode.vcf"
    output: "results/indiv_avinput/inherited.qc_wgs_phased_{VAR_TYPE}.sample.{SAMPLEID}.avinput"
    shell: """
        {annovar_dir}/convert2annovar.pl \
        -format vcf4 {input} \
        -outfile {output} \
        -includeinfo \
        -filter PASS \
        -comment
    """

rule rare_inherited_qc_genotype_indiv:
    input: "data/indiv_samples/inherited.wgs_phased_{VAR_TYPE}.sample.{SAMPLEID}.recode.vcf"
    params: sample = "{SAMPLEID}"
    output: "data/indiv_samples/inherited.qc_wgs_phased_{VAR_TYPE}.sample.{SAMPLEID}.recode.vcf"
    run:
        in_file = open(input[0], "r")
        out_file = open(output[0], "w")

        for line in in_file:
            if line.startswith("##"):
                out_file.write(line)
            elif line.startswith("#CHROM"):
                line_split = line.split()
                genotype_index = line_split.index(params.sample)
                out_file.write(line)
            else:
                line_split = line.split()
                genotype = line_split[genotype_index]
                pass_qc = genotype_quality_control(genotype)

                if pass_qc:
                    out_file.write(line)

rule rare_inherited_split_vcf_file:
    input: "data/for_annovar/wgs_phased_{VAR_TYPE}_filtered_inherited.vcf"
    params: 
        subject = "{SAMPLEID}",
        output_prefix = "data/indiv_samples/inherited.wgs_phased_{VAR_TYPE}.sample.{SAMPLEID}"
    output: "data/indiv_samples/inherited.wgs_phased_{VAR_TYPE}.sample.{SAMPLEID}.recode.vcf"
    shell: """
        {vcftools} \
        --vcf {input} \
        --indv {params.subject} \
        --out {params.output_prefix} \
        --recode \
        --recode-INFO-all
    """

rule rare_inherited_extract_inherited:
    input: "data/for_annovar/wgs_phased_{VAR_TYPE}_filtered.vcf"
    params:
        proband_id = "BM1452.001",
        father_id = "BM1452.200",
        mother_id = "BM1452.100",
        inherited_from = "father"
    output: "data/for_annovar/wgs_phased_{VAR_TYPE}_filtered_inherited.vcf"
    run:
        extract_inherited_variant(input[0], output[0], params.proband_id, params.father_id, params.mother_id, params.inherited_from)

################################################################################
# Remaining allele of chr22q11
################################################################################
rule check_22q11_genes:
    input:
        expand("results/annovar/output/counts/wgs_phased_{VAR_TYPE}.sample.{SAMPLEID}.remaining22q11.filtering_summary.txt", VAR_TYPE = var_type, SAMPLEID = probandid),
        expand("results/annovar/output/counts/remaining_22q11/wgs_phased_{VAR_TYPE}.sample.{SAMPLEID}.remaining22q11.exonic.nonsyn.deleterious.var_function.txt", VAR_TYPE = var_type, SAMPLEID = probandid),
        expand("results/annovar/output/counts/remaining_22q11/wgs_phased_{VAR_TYPE}.sample.{SAMPLEID}.remaining22q11.exonic.nonsyn.deleterious.exonic_function.txt", VAR_TYPE = var_type, SAMPLEID = probandid)

rule chr22q11_summarize_variant_function:
    input:
        "results/annovar/output/variant_filtering/remaining_22q11/wgs_phased_{VAR_TYPE}.sample.{SAMPLEID}.remaining22q11.exonic.nonsyn.deleterious.txt"
    output:
        "results/annovar/output/counts/remaining_22q11/wgs_phased_{VAR_TYPE}.sample.{SAMPLEID}.remaining22q11.exonic.nonsyn.deleterious.var_function.txt"
    run:
        summarize_variant_function(input[0], output[0])

rule chr22q11_summarize_exonic_function:
    input:
        "results/annovar/output/variant_filtering/remaining_22q11/wgs_phased_{VAR_TYPE}.sample.{SAMPLEID}.remaining22q11.exonic.nonsyn.deleterious.txt"
    output:
        "results/annovar/output/counts/remaining_22q11/wgs_phased_{VAR_TYPE}.sample.{SAMPLEID}.remaining22q11.exonic.nonsyn.deleterious.exonic_function.txt"
    run:
        summarize_exonic_functions(input[0], output[0])

rule chr22q11_summarize_filtering_steps:
    input:
        multianno_file = os.path.join(
            "results/annovar/output/multiannotation/indiv_samples/{SAMPLEID}",
            "wgs_phased_{VAR_TYPE}.sample.{SAMPLEID}.hg19_multianno.txt"
        ),
        chr22_file = os.path.join(
            "results/annovar/output/variant_filtering/remaining_22q11",
            "wgs_phased_{VAR_TYPE}.sample.{SAMPLEID}.remaining22q11.txt"
        ),
        noncoding_file = os.path.join(
            "results/annovar/output/variant_filtering/remaining_22q11",
            ("wgs_phased_{VAR_TYPE}.sample.{SAMPLEID}."
             "remaining22q11.noncoding.txt")
        ),
        exonic_file = os.path.join(
            "results/annovar/output/variant_filtering/remaining_22q11",
            ("wgs_phased_{VAR_TYPE}.sample.{SAMPLEID}."
             "remaining22q11.exonic.txt")
        ),
        nonsynonymous_file = os.path.join(
            "results/annovar/output/variant_filtering/remaining_22q11",
            ("wgs_phased_{VAR_TYPE}.sample.{SAMPLEID}."
             "remaining22q11.exonic.nonsyn.txt")
        ),
        deleterious_file = os.path.join(
            "results/annovar/output/variant_filtering/remaining_22q11",
            ("wgs_phased_{VAR_TYPE}.sample.{SAMPLEID}." 
             "remaining22q11.exonic.nonsyn.deleterious.txt")
        )
    output:
        os.path.join(
            "results/annovar/output/counts",
            ("wgs_phased_{VAR_TYPE}.sample.{SAMPLEID}."
             "remaining22q11.filtering_summary.txt")
        )
    shell: """
    echo -e "Filtering step\tVariants" > {output}
    grep -v -P "Chr\tStart" {input.multianno_file} \
        | grep -v -E "^#" \
        | awk 'END {{print "22q11.2 deleted region\t", NR}}' \
        >> {output}
    # total variants
    grep -v -P "Chr\tStart" {input.chr22_file} \
        | grep -v -E "^#" \
        | awk 'END {{print "22q11.2 deleted region\t", NR}}' \
        >> {output}
    # non-coding variants
    grep -v -P "Chr\tStart" {input.noncoding_file} \
        | grep -v -E "^#" \
        | awk 'END {{print "Non-coding\t", NR}}' \
        >> {output}
    # Exonic variants
    grep -v -P "Chr\tStart" {input.exonic_file} \
        | grep -v -E "^#" \
        | awk 'END {{print "Exonic\t", NR}}' \
        >> {output}
    # Non-synonymous variants
    grep -v -P "Chr\tStart" {input.nonsynonymous_file} \
        | grep -v -E "^#" \
        | awk 'END {{print "Non-synonymous\t", NR}}' \
        >> {output}
    # Deleterious variants
    grep -v -P "Chr\tStart" {input.deleterious_file} \
        | grep -v -E "^#" \
        | awk 'END {{print "Deleterious\t", NR}}' \
        >> {output}
    """

rule chr22q11_extract_deleterious:
    input: 
        "results/annovar/output/variant_filtering/remaining_22q11/wgs_phased_{VAR_TYPE}.sample.{SAMPLEID}.remaining22q11.exonic.nonsyn.txt"
    params: 
        deleterious_cutoff = "2"
    output: 
        "results/annovar/output/variant_filtering/remaining_22q11/wgs_phased_{VAR_TYPE}.sample.{SAMPLEID}.remaining22q11.exonic.nonsyn.deleterious.txt"
    run:
        extract_deleterious_variants(input[0], output[0], 
                                     params.deleterious_cutoff)

rule chr22q11_extract_nonsynonymous:
    input: 
        "results/annovar/output/variant_filtering/remaining_22q11/wgs_phased_{VAR_TYPE}.sample.{SAMPLEID}.remaining22q11.exonic.txt"
    output: 
        "results/annovar/output/variant_filtering/remaining_22q11/wgs_phased_{VAR_TYPE}.sample.{SAMPLEID}.remaining22q11.exonic.nonsyn.txt"
    run:
        extract_nonsynonymous_variants(input[0], output[0])

rule chr22q11_parse_coding_noncoding:
    input:
        "results/annovar/output/variant_filtering/remaining_22q11/wgs_phased_{VAR_TYPE}.sample.{SAMPLEID}.remaining22q11.txt"
    output:
        exonic = "results/annovar/output/variant_filtering/remaining_22q11/wgs_phased_{VAR_TYPE}.sample.{SAMPLEID}.remaining22q11.exonic.txt",
        noncoding = "results/annovar/output/variant_filtering/remaining_22q11/wgs_phased_{VAR_TYPE}.sample.{SAMPLEID}.remaining22q11.noncoding.txt"
    message:"""
    Parse exonic and noncoding variants.
    
        Input: {input}
        Output: {output}
    """
    run:
        parse_coding_vs_noncoding(input[0], output.exonic, output.noncoding)

rule extract_22q11_genes:
    input:
        "results/annovar/output/multiannotation/indiv_samples/{SAMPLEID}/wgs_phased_{VAR_TYPE}.sample.{SAMPLEID}.hg19_multianno.txt"
    params:
        chromosome = "chr22",
        start = "18656000",
        stop = "21792000"
    output: "results/annovar/output/variant_filtering/remaining_22q11/wgs_phased_{VAR_TYPE}.sample.{SAMPLEID}.remaining22q11.txt"
    message:"""
    Extract all variants on the remaining allele of chromosome 22q11.2.
    
        Chromosome: {params.chromosome}
        Start position: {params.start}
        Stop position: {params.stop}
        Input: {input}
        Ouptut: {output}
    """
    run:
        variants_file = open(input[0], "r")
        out_file = open(output[0], "w")
        
        for line in variants_file:
            if line.startswith("Chr\tStart"):
                out_file.write(line)
            elif line.startswith("#"):
                out_file.write(line)
            else:
                line_split = line.split()
                if all([line_split[0] == params.chromosome, 
                        int(line_split[1]) >= int(params.start), 
                        int(line_split[2]) <= int(params.stop)]):
                    out_file.write(line)

        variants_file.close()
        out_file.close()

################################################################################
# Variant filtering using all variants
################################################################################
rule check_indiv_filtering:
    input:
        expand("results/annovar/output/counts/wgs_phased_{VAR_TYPE}.sample.{SAMPLEID}.freq.{FREQ}.filtering_summary.txt", FREQ = var_freq, VAR_TYPE = var_type, SAMPLEID = sampleid),
        expand("results/annovar/output/counts/wgs_phased_{VAR_TYPE}.sample.{SAMPLEID}.freq.{FREQ}.nosuperdup.exonic.deleterious.var_function.txt", FREQ = var_freq, VAR_TYPE = var_type, SAMPLEID = sampleid),
        expand("results/annovar/output/counts/wgs_phased_{VAR_TYPE}.sample.{SAMPLEID}.freq.{FREQ}.nosuperdup.exonic.deleterious.exonic_function.txt", FREQ = var_freq, VAR_TYPE = var_type, SAMPLEID = sampleid)

rule indiv_summarize_variant_function:
    input:
        "results/annovar/output/variant_filtering/{FREQ}/wgs_phased_{VAR_TYPE}.sample.{SAMPLEID}.freq.{FREQ}.nosuperdup.exonic.deleterious.txt"
    output:
        "results/annovar/output/counts/wgs_phased_{VAR_TYPE}.sample.{SAMPLEID}.freq.{FREQ}.nosuperdup.exonic.deleterious.var_function.txt"
    run:
        summarize_variant_function(input[0], output[0])

rule indiv_summarize_exonic_function:
    input:
        "results/annovar/output/variant_filtering/{FREQ}/wgs_phased_{VAR_TYPE}.sample.{SAMPLEID}.freq.{FREQ}.nosuperdup.exonic.deleterious.txt"
    output:
        "results/annovar/output/counts/wgs_phased_{VAR_TYPE}.sample.{SAMPLEID}.freq.{FREQ}.nosuperdup.exonic.deleterious.exonic_function.txt"
    run:
        summarize_exonic_functions(input[0], output[0])

rule summarize_filtering_steps_indiv:
    input:
        multianno_file = "results/annovar/output/multiannotation/indiv_samples/{SAMPLEID}/wgs_phased_{VAR_TYPE}.sample.{SAMPLEID}.hg19_multianno.txt",
        freq_file = "results/annovar/output/variant_filtering/{FREQ}/wgs_phased_{VAR_TYPE}.sample.{SAMPLEID}.freq.{FREQ}.txt",
        superdup_file = "results/annovar/output/variant_filtering/{FREQ}/wgs_phased_{VAR_TYPE}.sample.{SAMPLEID}.freq.{FREQ}.nosuperdup.txt",
        noncoding_file = "results/annovar/output/variant_filtering/{FREQ}/wgs_phased_{VAR_TYPE}.sample.{SAMPLEID}.freq.{FREQ}.nosuperdup.noncoding.txt",
        exonic_file = "results/annovar/output/variant_filtering/{FREQ}/wgs_phased_{VAR_TYPE}.sample.{SAMPLEID}.freq.{FREQ}.nosuperdup.exonic.txt",
        nonsynonymous_file = "results/annovar/output/variant_filtering/{FREQ}/wgs_phased_{VAR_TYPE}.sample.{SAMPLEID}.freq.{FREQ}.nosuperdup.exonic.nonsyn.txt",
        deleterious_file = "results/annovar/output/variant_filtering/{FREQ}/wgs_phased_{VAR_TYPE}.sample.{SAMPLEID}.freq.{FREQ}.nosuperdup.exonic.deleterious.txt"
    params:
        frequency = "{FREQ}"
    output:
        "results/annovar/output/counts/wgs_phased_{VAR_TYPE}.sample.{SAMPLEID}.freq.{FREQ}.filtering_summary.txt"
    shell: """
    echo -e "Filtering step\tVariants" > {output}
    grep -v -P "Chr\tStart" {input.multianno_file} | grep -v -E "^#" | awk 'END {{print "All variants\t", NR}}' >> {output}
    # total variants
    grep -v -P "Chr\tStart" {input.freq_file} | grep -v -E "^#" | awk 'END {{print "Allele frequency: {params.frequency}\t", NR}}' >> {output}
    # Remove superdup variants
    grep -v -P "Chr\tStart" {input.superdup_file} | grep -v -E "^#" | awk 'END {{print "Excluding genomic superdups\t", NR}}' >> {output}
    # non-coding variants
    grep -v -P "Chr\tStart" {input.noncoding_file} | grep -v -E "^#" | awk 'END {{print "Non-coding variants\t", NR}}' >> {output}
    # Exonic variants
    grep -v -P "Chr\tStart" {input.exonic_file} | grep -v -E "^#" | awk 'END {{print "Exonic variants\t", NR}}' >> {output}
    # Non-synonymous variants
    grep -v -P "Chr\tStart" {input.nonsynonymous_file} | grep -v -E "^#" | awk 'END {{print "Non-synonymous variants\t", NR}}' >> {output}
    # Deleterious variants
    grep -v -P "Chr\tStart" {input.deleterious_file} | grep -v -E "^#" | awk 'END {{print "Deleterious variants\t", NR}}' >> {output}
    """

rule extract_deleterious_indiv:
    input:
        "results/annovar/output/variant_filtering/{FREQ}/wgs_phased_{VAR_TYPE}.sample.{SAMPLEID}.freq.{FREQ}.nosuperdup.exonic.nonsyn.txt"
    params:
        deleterious_count_cutoff = "2"
    output:
        "results/annovar/output/variant_filtering/{FREQ}/wgs_phased_{VAR_TYPE}.sample.{SAMPLEID}.freq.{FREQ}.nosuperdup.exonic.deleterious.txt"
    message: """
    Extract deleterious exonic variants.
    
    Input: {input}
    Output: {output}
    """
    run:
        extract_deleterious_variants(input[0], output[0], params.deleterious_count_cutoff)

rule extract_nonsynonymous_indiv:
    input:
        "results/annovar/output/variant_filtering/{FREQ}/wgs_phased_{VAR_TYPE}.sample.{SAMPLEID}.freq.{FREQ}.nosuperdup.exonic.txt"
    output:
        "results/annovar/output/variant_filtering/{FREQ}/wgs_phased_{VAR_TYPE}.sample.{SAMPLEID}.freq.{FREQ}.nosuperdup.exonic.nonsyn.txt"
    message:"""
    Extract non-synonymous variants.
    
    Input: {input}
    Output: {output}
    """
    run:
        extract_nonsynonymous_variants(input[0], output[0])

rule parse_coding_noncoding_indiv:
    input:
        "results/annovar/output/variant_filtering/{FREQ}/wgs_phased_{VAR_TYPE}.sample.{SAMPLEID}.freq.{FREQ}.nosuperdup.txt"
    output:
        exonic = "results/annovar/output/variant_filtering/{FREQ}/wgs_phased_{VAR_TYPE}.sample.{SAMPLEID}.freq.{FREQ}.nosuperdup.exonic.txt",
        noncoding = "results/annovar/output/variant_filtering/{FREQ}/wgs_phased_{VAR_TYPE}.sample.{SAMPLEID}.freq.{FREQ}.nosuperdup.noncoding.txt"
    message:"""
    Split input file into coding and noncoding variants.
    
    Input: {input}
    Output: {output}
    """
    run:
        parse_coding_vs_noncoding(input[0], output.exonic, output.noncoding)

rule remove_superdup_indiv:
    input: "results/annovar/output/variant_filtering/{FREQ}/wgs_phased_{VAR_TYPE}.sample.{SAMPLEID}.freq.{FREQ}.txt"
    output: "results/annovar/output/variant_filtering/{FREQ}/wgs_phased_{VAR_TYPE}.sample.{SAMPLEID}.freq.{FREQ}.nosuperdup.txt"
    message:"""
    Removing variants found in genomic superdups.
    
    Input: {input}
    Output: {output}
    """
    run:
        remove_superdup(input[0], output[0])

rule filter_on_variant_frequency:
    input:
        "results/annovar/output/multiannotation/indiv_samples/{SAMPLEID}/wgs_phased_{VAR_TYPE}.sample.{SAMPLEID}.hg19_multianno.txt"
    params:
        frequency = "{FREQ}"
    output:
        "results/annovar/output/variant_filtering/{FREQ}/wgs_phased_{VAR_TYPE}.sample.{SAMPLEID}.freq.{FREQ}.txt"
    run:
        if params.frequency == "all":
            shell("cp {input} {output}")
        elif params.frequency == "novel":
            extract_novel_variants(input[0], output[0])
        else:
            maf = float(params.frequency[4:])
            extract_rare_variants(input[0], output[0], maf)

################################################################################
# Perform ANNOVAR on variants from individual subjects
################################################################################

rule combine_population_counts:
    input:
        expand("results/annovar/output/counts/wgs_phased_{VAR_TYPE}.sample.{SAMPLEID}.pop_counts.txt", VAR_TYPE = var_type, SAMPLEID = sampleid),
        expand("results/annovar/output/counts/wgs_phased_{VAR_TYPE}.all_sample.pop_counts.txt", VAR_TYPE = var_type)
    output:
        snp_output = "results/annovar/output/counts/wgs_phased_snp.all_samples.pop_counts.txt",
        indel_ouptut = "results/annovar/output/counts/wgs_phased_indel.all_samples.pop_counts.txt"
    message:"""
    Combine variant counts from individuals into a single table.
    """
    shell:"""
        Rscript src/annovar_wgs_samples/merge_variant_count_summaries.R \
        {output[0]} \
        {output[1]} \
        {input}
    """

rule indiv_population_counts:
    input: 
        multianno = "results/annovar/output/multiannotation/indiv_samples/{SAMPLEID}/wgs_phased_{VAR_TYPE}.sample.{SAMPLEID}.hg19_multianno.txt",
        dbsnp     = "results/annovar/output/multiannotation/indiv_samples/{SAMPLEID}/wgs_phased_{VAR_TYPE}.sample.{SAMPLEID}.hg19_snp138_dropped",
        onekg     = "results/annovar/output/multiannotation/indiv_samples/{SAMPLEID}/wgs_phased_{VAR_TYPE}.sample.{SAMPLEID}.hg19_ALL.sites.2012_04_dropped",
        esp       = "results/annovar/output/multiannotation/indiv_samples/{SAMPLEID}/wgs_phased_{VAR_TYPE}.sample.{SAMPLEID}.hg19_esp6500si_all_dropped"
    output: 
        "results/annovar/output/counts/wgs_phased_{VAR_TYPE}.sample.{SAMPLEID}.pop_counts.txt"
    message:"""
    Summarize counts of the different annotations.
    
    Get the counts of variants that have annotations for:
        
        * All variants
        * Found in dbSNP138
        * Found in 1kg 
        * Found in ESP6500
    
    Input: {input}
    Output: {output}
    """
    shell: """
        # total variants
        grep -v -P "Chr\tStart" {input.multianno} | grep -v -E "^#" | awk 'END {{print "total\t", NR}}' > {output}
        awk 'END {{print "In dbSNP138\t", NR}}' {input.dbsnp} >> {output}
        awk 'END {{print "In 1kg (2012/10)\t", NR}}' {input.onekg} >> {output}
        awk 'END {{print "In ESP6500\t", NR}}' {input.esp} >> {output}
    """

rule table_annotation_indiv:
    input: "results/annovar/input/qc_wgs_phased_{VAR_TYPE}.sample.{SAMPLEID}.avinput"
    params: output_prefix = "results/annovar/output/multiannotation/indiv_samples/{SAMPLEID}/wgs_phased_{VAR_TYPE}.sample.{SAMPLEID}"
    output: 
        "results/annovar/output/multiannotation/indiv_samples/{SAMPLEID}/wgs_phased_{VAR_TYPE}.sample.{SAMPLEID}.hg19_multianno.txt",
        "results/annovar/output/multiannotation/indiv_samples/{SAMPLEID}/wgs_phased_{VAR_TYPE}.sample.{SAMPLEID}.hg19_snp138_dropped",
        "results/annovar/output/multiannotation/indiv_samples/{SAMPLEID}/wgs_phased_{VAR_TYPE}.sample.{SAMPLEID}.hg19_ALL.sites.2012_04_dropped",
        "results/annovar/output/multiannotation/indiv_samples/{SAMPLEID}/wgs_phased_{VAR_TYPE}.sample.{SAMPLEID}.hg19_esp6500si_all_dropped"
    message: """
    Perform table annotation using annovar.
    
    Annotate variants using:
    
        * refGene
        * phastConsElements46way
        * genomicSuperDups
        * esp6500si_all
        * 1000g2012apr_all
        * snp138
        * ljb23_all
    
    Input: {input}
    Output: {output}
    """
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

################################################################################
# All sample CADD table annotations
################################################################################
rule filter_cadd:
    input: "results/annovar_output/cadd_annotation/wgs_phased_snp.cadd.hg19_cadd_dropped"
    params: scaled_cadd_cutoff = "20"
    output: "results/annovar_output/cadd_annotation/wgs_phased_snp.cadd.hg19_cadd_dropped.filtered"
    run:
        in_file = open(input[0], "r")
        out_file = open(output[0], "w")

        for line in in_file:
            split_line = line.split("\t")
            cadd_score = split_line[1].split(",")
            scaled_cadd = cadd_score[1]
            if float(scaled_cadd) >= int(params.scaled_cadd_cutoff):
                out_file.write(line)
        in_file.close()
        out_file.close()

rule annotate_novel_cadd:
    input: "results/indiv_avinput/qc_wgs_phased_snp.sample.BM1452.001.avinput"
    params: output_prefix = "results/annovar_output/cadd_annotation/wgs_phased_snp.cadd"
    output: "results/annovar_output/cadd_annotation/wgs_phased_snp.cadd.hg19_cadd_dropped"
    shell: """
        $HOME/programs/annovar/annovar2014Jul14/annotate_variation.pl \
        {input} \
        $HOME/programs/annovar/annovar2014Jul14/humandb/ \
        -filter \
        -dbtype cadd \
        -buildver hg19 \
        -out {params.output_prefix} \
        -otherinfo
    """

################################################################################
# Prepare vcf files for individual samples
################################################################################
rule convert_to_annovar_indiv:
    input: "results/annovar/input/qc_wgs_phased_{VAR_TYPE}.sample.{SAMPLEID}.recode.vcf"
    output: "results/annovar/input/qc_wgs_phased_{VAR_TYPE}.sample.{SAMPLEID}.avinput"
    message:"""
    Convert VCF file to avinput
    
    Input: {input}
    Output: {output}
    """
    shell: """
        {annovar_dir}/convert2annovar.pl \
        -format vcf4 {input} \
        -outfile {output} \
        -includeinfo \
        -filter PASS \
        -comment
    """

rule qc_genotype_indiv:
    input: 
        "results/annovar/input/wgs_phased_{VAR_TYPE}.sample.{SAMPLEID}.recode.vcf"
    params: 
        sample = "{SAMPLEID}"
    output: 
        "results/annovar/input/qc_wgs_phased_{VAR_TYPE}.sample.{SAMPLEID}.recode.vcf"
    message:"""
    Perform quality control on individual genotypes, checking for:
    
        * Alternate allele ratio
        * Read depth
        * Genotype quality
    
    Sample ID: {params.sample}
    Input: {input}
    Output: {output}
    """
    run:
        in_file = open(input[0], "r")
        out_file = open(output[0], "w")

        for line in in_file:
            if line.startswith("##"):
                out_file.write(line)
            elif line.startswith("#CHROM"):
                line_split = line.split()
                genotype_index = line_split.index(params.sample)
                out_file.write(line)
            else:
                line_split = line.split()
                genotype = line_split[genotype_index]
                pass_qc = genotype_quality_control(genotype)

                if pass_qc:
                    out_file.write(line)

rule split_vcf_file:
    input: "data/wgs_phased_{VAR_TYPE}.vcf"
    params: 
        subject = "{SAMPLEID}",
        output_prefix = "results/annovar/input/wgs_phased_{VAR_TYPE}.sample.{SAMPLEID}"
    output: "results/annovar/input/wgs_phased_{VAR_TYPE}.sample.{SAMPLEID}.recode.vcf"
    message:"""
    Split multi-sample VCF file into individual VCF files: 
    
    Sample ID: {params.subject}
    Input: {input}
    output: {output}
    """
    shell: """
        {vcftools} \
        --vcf {input} \
        --indv {params.subject} \
        --out {params.output_prefix} \
        --recode \
        --recode-INFO-all
    """

################################################################################
# All sample table annotations
################################################################################
rule test_all_sample_annotation:
    input:
        expand(
            os.path.join("results/annovar/output/counts",
                         "wgs_phased_{VAR_TYPE}.all_sample.pop_counts.txt"), 
            VAR_TYPE = var_type
        )

rule all_sample_population_counts:
    input: 
        multianno = os.path.join(
            "results/annovar/output/multiannotation/all_samples",
            "wgs_phased_{VAR_TYPE}.hg19_multianno.txt"),
        dbsnp = os.path.join(
            "results/annovar/output/multiannotation/all_samples",
            "wgs_phased_{VAR_TYPE}.hg19_snp138_dropped"),
        onekg = os.path.join(
            "results/annovar/output/multiannotation/all_samples",
            "wgs_phased_{VAR_TYPE}.hg19_ALL.sites.2012_04_dropped"),
        esp = os.path.join(
            "results/annovar/output/multiannotation/all_samples",
            "wgs_phased_{VAR_TYPE}.hg19_esp6500si_all_dropped")
    output: 
        os.path.join("results/annovar/output/counts",
                     "wgs_phased_{VAR_TYPE}.all_sample.pop_counts.txt")
    message:"""
Summarize counts of the different annotations.

Get the counts of variants that were previously reported:
    
    * All variants
    * Found in dbSNP138
    * Found in 1kg 
    * Found in ESP6500
    
    Input: {input}
    Output: {output}
    """
    shell: """
    echo -e "Category\tCount" \
        > {output}
    grep -v -P "Chr\tStart" {input.multianno} \
        | grep -v -E "^#" \
        | awk 'END {{print "Total variants\t", NR}}' \
        >> {output}
    awk 'END {{print "In dbSNP138\t", NR}}' {input.dbsnp} \
        >> {output}
    awk 'END {{print "In 1kg (2012/10)\t", NR}}' {input.onekg} \
        >> {output}
    awk 'END {{print "In ESP6500\t", NR}}' {input.esp} \
        >> {output}
    """

rule table_annotation_all_samples:
    input: 
        "results/annovar/input/allsample/wgs_phased_{VAR_TYPE}.avinput"
    params: 
        output_prefix = os.path.join(
            "results/annovar/output/multiannotation/all_samples",
            "wgs_phased_{VAR_TYPE}"),
        protocol = ",".join(["refGene", "phastConsElements46way", 
                             "genomicSuperDups", "esp6500si_all", 
                             "1000g2012apr_all", "snp138", "ljb23_all"]),
        operation = "g,r,r,f,f,f,f",
        build = "hg19"
    output:
        os.path.join("results/annovar/output/multiannotation/all_samples/",
                     "wgs_phased_{VAR_TYPE}.hg19_multianno.txt"),
        os.path.join("results/annovar/output/multiannotation/all_samples/",
                     "wgs_phased_{VAR_TYPE}.hg19_snp138_dropped"),
        os.path.join("results/annovar/output/multiannotation/all_samples/",
                     "wgs_phased_{VAR_TYPE}.hg19_ALL.sites.2012_04_dropped"),
        os.path.join("results/annovar/output/multiannotation/all_samples/",
                     "wgs_phased_{VAR_TYPE}.hg19_esp6500si_all_dropped")
    message: """
    Perform table annotation using annovar.
    
    Annotate variants using:
        
        * refGene
        * phastConsElements46way
        * genomicSuperDups
        * esp6500si_all
        * 1000g2012apr_all
        * snp138
        * ljb23_all

    Input: {input}
    Output: {output}
    """
    shell: """
    {annovar_dir}/table_annovar.pl \
    {input} \
    {annovar_dir}/humandb/ \
    --protocol {params.protocol} \
    --operation {params.operation} \
    --buildver {params.build} \
    --outfile {params.output_prefix} \
    --otherinfo
    """

rule convert_allsample_to_avinput:
    input: "results/annovar/input/allsample/wgs_phased_{VAR_TYPE}_filtered.vcf"
    output: "results/annovar/input/allsample/wgs_phased_{VAR_TYPE}.avinput"
    message: """
    Converting {input} into format for annovar. 
    
    convert2annovar options:
        
        '--format vcf4old` to keep all variants
    
    Input: {input}
    Output: {output}
    """
    shell:"""
        {annovar_dir}/convert2annovar.pl \
        --format vcf4old {input} \
        --outfile {output} \
        --includeinfo \
        --filter PASS \
        --allallele \
        --comment
    """

rule filter_vcf_file:
    input: 
        "data/wgs_phased_{VAR_TYPE}.vcf"
    output: 
        "results/annovar/input/allsample/wgs_phased_{VAR_TYPE}_filtered.vcf"
    message:"""
    Perform quality control on individual genotypes, checking for:
    
        * Alternate allele ratio
        * Read depth
        * Genotype quality
        * Alternate allele and non-missing genotype
    
    Input: {input}
    Output: {output}
    """
    run:
        # open input and output files
        in_file = open(input[0], "r")
        out_file = open(output[0], "w")
        
        for line in in_file:
            # Extract header information
            if line.startswith("##"):
                out_file.write(line)
            elif line.startswith("#CHROM\tPOS"):
                header = line.split()
                genotype_start_index = header.index("FORMAT") + 1
                genotype_end_index = len(header)
                out_file.write(line)
            else:
                # Check each genotype to see if it passes QC and that subject
                # has an alternate allele.
                line_split = line.split()
                test_genotype = []
                for genotype_index in range(genotype_start_index, 
                                            genotype_end_index):
                    pass_qc = genotype_quality_control(
                        line_split[genotype_index]
                    )
                    check_alt_allele = check_alt_in_genotype(
                        line_split[genotype_index]
                    )
                    test_genotype.append(all([pass_qc, check_alt_allele]))
                    # test_genotype.append(pass_qc)
                if any(test_genotype):
                    # Output variant if any samples pass QC and have an 
                    # alternate allele.
                    out_file.write(line)
        in_file.close()
        out_file.close()
