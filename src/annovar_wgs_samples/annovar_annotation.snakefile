import re
from snakemake.utils import R
from subprocess import call
import os

workdir: "/home/jchung/projects/wgs/"
# Path to required programs
annovar_dir = '/home/jchung/programs/annovar/annovar2014Jul14/'
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
    
    Creates an output file without variants in genomic superdups.
    
    Args:
        input: Path to the input file.
        output: Path to the output file.
    Returns:
        None
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

def extract_conserved_elements(input, output, conservation_column, 
                               score_cutoff):
    """Extract variants that have conservation scores >= score_cutoff.
    
    Creates an output file with variants that have a conservation score
        >= score cutoff.
    
    Args:
        input: Path to input file.
        output: Path to output file.
        conservation_column: Name of the column containing conservation scores.
        score_cutoff: Minimum conservation score.
    Returns:
        None
    """
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
    
    Creates an output file with only noncoding variants.
    
    Args:
        input: Path to input file.
        output: Path to output file.
    Returns:
        None
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
    
    Creates an output file with nonsynonymous variants.
    
    Args:
        input: Path to input file
        output: Path to output file
    Returns:
        None
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
    """Extract deleterious variants from input file.
    
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
    Returns:
        None
    """
    deleterious_score_cutoff = int(deleterious_score_cutoff)
    in_file = open(input, "r")
    out_file = open(output, "w")
    
    re_frameshift = re.compile(
        "|".join(["\sframeshift insertion\s", "\sframeshift deletion\s",
                  "\sframeshift block substitution\s"])
    )
    re_stopalter = re.compile("\sstopgain\s|\sstoploss\s")
    re_nonsyn = re.compile("nonsynonymous")
    re_splice = re.compile("splicing")
    
    
    for line in in_file:
        split_line = line.split("\t")
        
        if line.startswith("Chr\tStart"):
            sift_pred_index = split_line.index("LJB23_SIFT_pred")
            polyphen2_hdiv_pred_index = split_line.index(
                "LJB23_Polyphen2_HDIV_pred"
            )
            lrt_pred_index = split_line.index("LJB23_LRT_pred")
            mutationtaster_pred_index = split_line.index(
                "LJB23_MutationTaster_pred"
            )
            out_file.write(line)
        elif line.startswith("#"):
            out_file.write(line)
        elif re_nonsyn.search(line):
            # Test for deleteriousness of nonsynonymous variants
            sift = split_line[sift_pred_index] == "D"
            polyphen2 = (split_line[polyphen2_hdiv_pred_index] == "D" or 
                         split_line[polyphen2_hdiv_pred_index] == "P")
            lrt = split_line[lrt_pred_index] == "D"
            mutationtaster = (split_line[mutationtaster_pred_index] == "A" or
                              split_line[mutationtaster_pred_index] == "D")
            
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

def extract_novel_variants(input, output, dbsnp_database = "generic", 
                           esp_database = "esp6500si_all", 
                           thousand_genomes_database = "1000g2012apr_all"):
    """Remove variants from population databases.
    
    Creates an output file with all variants in dbSNP, ESP, and 1KG removed.
    
    Args:
        input: Path to input file.
        output: Path to output file.
        dbsnp_database: Name of column with dbSNP information.
        esp_database: Name of column with exome sequencing project information.
        thousand_genomes_database: Name of column with 1000 Genomes information.
    Returns:
        None
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
    Returns:
        None
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
        bool: True if genotype contains alternate allele, False if homozygous
            reference, None if missing.
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
        String indicating the type allele status. `heterozygous' if one minor
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
        bool: True if inheritance passes sanity check, False otherwise.
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

def extract_inherited_variant(input, output, proband_id, father_id, mother_id, 
    inherited_from):
    """Extract variants that are inherited from specified parent parent.
    
    Args:
        input: Path to input file.
        ouptut: Path to output file.
        proband_id: Character string indicating proband ID.
        father_id: Character string indicating father ID.
        mother_id: Character string indicating mother ID.
        inherited_from: Character string one of ["mother", "father"] indicating
            which parent to check for inheritence.
    Return:
        None
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
                if check_inherited_variant(proband_genotype, father_genotype, 
                                           mother_genotype, inherited_from):
                    out_file.write(line)
    out_file.close()
    in_file.close()

def check_alt_allele_ratio(genotype, allele_counts, hom_ref_cutoff = 0.15,
                           het_low_cutoff = 0.3, het_high_cutoff = 0.7, 
                           hom_alt_cutoff = 0.85):
    """Check if alternate allele ratio is within QC bounds.
    
    Calculate alternate allele ratio and determine if it falls within QC
    measures.
    
    Args:
        genotype: Subject genotype.
        allele_counts: Genotype allele counts.
        hom_ref_cutoff: Maximum alternate allele ratio for homozygous reference
            genotype.
        het_low_cutoff: Minimum alternate allele ratio for heterozygous 
            genotype.
        het_high_cutoff: Maximum alternate allele ratio for heterozygous
            genotype.
        hom_alt_cutoff: Minimum alternate allele ratio for homozygous alternate
            genotype.
    Returns:
        bool: True if genotype passes QC, False otherwise.
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

def check_genotype_quality(quality_score, score_cutoff = 20):
    """Check if genotype quality is >= cutoff.
    
    Args:
        quality_score: Genotype quality score.
        score_cutoff: Minimum genotype quality score.
    Returns:
        bool: True if quality score passes QC, False otherwise.
    """
    if int(quality_score) < score_cutoff:
        return False
    else:
        return True

def check_read_depth(read_depth, depth_cutoff = 15):
    """Check if genotype read depth is >= cutoff.
    
    Args:
        read_depth: Genotype quality score.
        depth_cutoff: Minimum genotype read depth.
    Returns:
        bool: True if read depth passes QC, False otherwise.
    """
    if int(read_depth) < depth_cutoff:
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
        bool: True if genotype passes all QC tests, False otherwise.
    """
    genotype_split = genotype_info.split(":")
    if "." in genotype_split[0]:
        return None
    else:
        pass_alt_allele_ratio = check_alt_allele_ratio(genotype_split[0], 
                                                       genotype_split[1])
        pass_depth = check_read_depth(genotype_split[2])
        pass_genotype_quality = check_genotype_quality(genotype_split[3])
        pass_qc = all([pass_alt_allele_ratio, pass_depth, 
                       pass_genotype_quality])
        return pass_qc

def extract_denovo_variants(input, output, proband_id, father_id, mother_id):
    """Extract *de novo* variants.
    
    Extracts *de novo* variants from input file and saves them into the output.
    
    Args:
        input: Path to input file.
        output: Path to output file.
        proband_id: String describing proband ID.
        father_id: String describing father ID.
        mother_id: String describing mother ID.
    Returns:
        None
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
            
            # print(proband_genotype)
            if "." not in proband_genotype:
                proband_alt = check_alt_in_genotype(proband_genotype)
                father_alt = check_alt_in_genotype(father_genotype)
                mother_alt = check_alt_in_genotype(mother_genotype)
                
                if not any([proband_alt == None, 
                            father_alt == None, 
                            mother_alt == None]):
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
            elif not line.startswith("#"):
                var_annotation = line.split("\t")
                var_exonic_function = var_annotation[exonic_function_index]
                all_exonic_function.append(var_exonic_function)
    
    frameshift_insertion = all_exonic_function.count("frameshift insertion")
    frameshift_deletion = all_exonic_function.count("frameshift deletion")
    frameshift_block_sub = all_exonic_function.count(
        "frameshift block substitution"
    )
    stopgain = all_exonic_function.count("stopgain")
    stoploss = all_exonic_function.count("stoploss")
    nonframeshift_insertion = all_exonic_function.count(
        "nonframeshift insertion"
    )
    nonframeshift_deletion = all_exonic_function.count("nonframeshift deletion")
    nonframeshift_block_sub = all_exonic_function.count(
        "nonframeshift block substitution"
    )
    nonsynonymous_snv = all_exonic_function.count("nonsynonymous SNV")
    synonymous_snv = all_exonic_function.count("synonymous SNV")
    unknown = all_exonic_function.count("unknown")
    
    with open(output, "w") as out_file:
        out_file.write("Exonic function\tCount\n")
        out_file.write(
            "Frameshift insertion\t{0}\n".format(frameshift_insertion)
        )
        out_file.write("Frameshift deletion\t{0}\n".format(frameshift_deletion))
        out_file.write(
            "Frameshift block substitution\t{0}\n".format(frameshift_block_sub)
        )
        out_file.write("Stopgain\t{0}\n".format(stopgain))
        out_file.write("Stoploss\t{0}\n".format(stoploss))
        out_file.write(
            "Nonframeshift insertion\t{0}\n".format(nonframeshift_insertion)
        )
        out_file.write(
            "Nonframeshift deletion\t{0}\n".format(nonframeshift_deletion)
        )
        out_file.write(
            "Nonframeshift block substitution\t{0}\n".format(
                nonframeshift_block_sub
            )
        )
        out_file.write("Nonsynonymous SNV\t{0}\n".format(nonsynonymous_snv))
        out_file.write("Synonymous SNV\t{0}\n".format(synonymous_snv))
        out_file.write("Unknown\t{0}\n".format(unknown))

def summarize_variant_function(input, output):
    """Count variants by variant function.
    
    Args:
        input: Path to input file.
        output: Path to output file.
    Returns:
        None
    """
    variant_function_column = "Func.refGene"
    all_variant_function = []
    with open(input) as in_file:
        for line in in_file:
            if line.startswith("Chr\tStart"):
                header = line.split("\t")
                variant_function_column = header.index(variant_function_column)
            elif not line.startswith("#"):
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

def count_population_database(input, output, dbsnp_name = "snp138",
                              thousand_genomes_name = "1000g2012apr_all",
                              esp_name = "esp6500si_all",
                              custom_dbsnp_name = "generic"):
    """Get a count of population database annotations.
    
    Args:
        input: Path to multianno ANNOVAR file.
        output: Path to output file
    Returns:
        None
    """
    
    dbsnp = []
    custom_dbsnp = []
    thousand_genomes = []
    esp = []
    in_db = 0
    with open(input) as in_file:
        for line in in_file:
            if line.startswith("Chr\tStart"):
                header = line.split("\t")
                dbsnp_index = header.index(dbsnp_name)
                custom_dbsnp_index = header.index(custom_dbsnp_name)
                thousand_index = header.index(thousand_genomes_name)
                esp_index = header.index(esp_name)
            elif not line.startswith("#"):
                var_annotation = line.split("\t")
                dbsnp.append(var_annotation[dbsnp_index])
                thousand_genomes.append(var_annotation[thousand_index])
                esp.append(var_annotation[esp_index])
                if not all([var_annotation[dbsnp_index] == "",
                            var_annotation[thousand_index] == "",
                            var_annotation[esp_index] == ""]):
                    in_db = in_db + 1
    total = len(dbsnp)
    in_dbsnp = total - dbsnp.count("")
    in_custom_dbsnp = total - custom_dbsnp.count("")
    in_thousand_genomes = total - thousand_genomes.count("")
    in_esp = total - esp.count("")
    novel = total - in_db
    
    with open(output, "w") as out_file:
        out_file.write("Database\tCount\n")
        out_file.write("All variants\t{0}\n".format(total))
        out_file.write("In any DB\t{0}\n".format(in_db))
        out_file.write("In dbSNP\t{0}\n".format(in_dbsnp))
        out_file.write("In custom dbsnp138\t{0}\n".format(in_custom_dbsnp))
        out_file.write("In 1KG\t{0}\n".format(in_thousand_genomes))
        out_file.write("In ESP6500\t{0}\n".format(in_esp))
        out_file.write("Novel\t{0}\n".format(novel))

def summarize_filtering_results(multianno_file, var_filter_file, superdup_file,
                                noncoding_file, exonic_file, nonsynonymous_file,
                                deleterious_file, output, var_filter_type):
    """Create table summarizing filtering results.
    
    Args:
        multianno_file: Path to multianno ANNOVAR file.
        var_filter_file: Path to variants after initial filtering.
        superdup_file: Path to file after removeing variants in superdup 
            regions.
        noncoding_file: Path to file containing non-coding variants.
        exonic_file: Path to file containing exonic variants.
        nonsynonymous_file: Path to file with nonsynonymous variants.
        deleterious_file: Path to file with deleterious variants.
        output: Path to output file.
        var_filter_type: Type of variant. One of ``["snv", "indel"]``.
    Returns:
        None
    """
    multianno_count = count_vars_in_file(multianno_file)
    var_filter_count = count_vars_in_file(var_filter_file)
    
    noncoding_count = count_vars_in_file(noncoding_file)
    exonic_count = count_vars_in_file(exonic_file)
    nonsynonymous_count = count_vars_in_file(nonsynonymous_file)
    deleterious_count = count_vars_in_file(deleterious_file)
    
    # For 22q11 region, I don't filter on superdups.
    if superdup_file == "None":
        superdup_count = "NA"
    else:
        superdup_count = count_vars_in_file(superdup_file)
    
    with open(output, "w") as out_file:
        out_file.write("Filtering step\tCount\n")
        out_file.write("All variants\t{0}\n".format(multianno_count))
        out_file.write("{0}\t{1}\n".format(var_filter_type, var_filter_count))
        out_file.write("Genomic SuperDups\t{0}\n".format(superdup_count))
        out_file.write("Non-coding variants\t{0}\n".format(noncoding_count))
        out_file.write("Exonic variants\t{0}\n".format(exonic_count))
        out_file.write("Non-synonymous variants\t{0}\n".format(
            nonsynonymous_count
        ))
        out_file.write("Deleterious variants\t{0}\n".format(deleterious_count))

def count_vars_in_file(input):
    """Count variants in input file.
    
    Args:
        input: Path to input file.
    Returns:
        None
    """
    lines = 0
    with open(input) as in_file:
        for line in in_file:
            if (line.startswith("#") or 
                line.startswith("Chr\tStart")):
                # print("Skipping header")
                pass
            else:
                lines = lines + 1
    return(lines)

################################################################################
# Begin Rules
################################################################################
rule all:
    input:
        os.path.join(
            "results/variant_counts",
            ("wgs_phased_all_variants.hg19_multianno."
             "genic.lod100.tbx1_counts.txt")
         ),
        os.path.join(
            "results/variant_counts",
            "wgs_phased_all_variants.hg19_multianno.lod100.tbx1_counts.txt"
        ),
        "results/novel_variants/novel_deleterious.snv.txt",
        "results/novel_variants/novel_deleterious.indel.txt",
        "results/rare_variants/rare_genic_deleterious.snv.txt",
        "results/rare_variants/rare_genic_deleterious.indel.txt",
        expand(os.path.join("results/annovar/output/refgene/all_samples/",
                     "wgs_phased_{VAR_TYPE}.refgene_hgvs.variant_function"),
              VAR_TYPE = var_type),
        expand(os.path.join("results/annovar/output/refgene/all_samples/",
                     "wgs_phased_{VAR_TYPE}.refgene_hgvs.exonic_variant_function"),
              VAR_TYPE = var_type)

rule check_hgvs_annotation:
    input:
        expand(os.path.join("results/annovar/output/refgene/all_samples/",
                     "wgs_phased_{VAR_TYPE}.refgene_hgvs.variant_function"),
              VAR_TYPE = var_type),
        expand(os.path.join("results/annovar/output/refgene/all_samples/",
                     "wgs_phased_{VAR_TYPE}.refgene_hgvs.exonic_variant_function"),
              VAR_TYPE = var_type)

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
        os.path.join(
            "results/annovar/output/genic_lod100",
            "wgs_phased_all_variants.hg19_multianno.genic.lod100.txt"
        ),
        "data/tbx1_pathway_genes_2014_03_27.txt"
    output:
        os.path.join(
            "results/variant_counts",
            "wgs_phased_all_variants.hg19_multianno.genic.lod100.tbx1_counts.txt"
        )
    shell: """
    perl src/perl/count_var_per_gene.pl \
    {input[1]} \
    {input[0]} \
    {output}
    """

rule count_variants_tbx1_pathway_all_lod100:
    input:
        os.path.join(
            "results/annovar/output/all_lod100",
            "wgs_phased_all_variants.hg19_multianno.lod100.txt"
        ),
        "data/tbx1_pathway_genes_2014_03_27.txt"
    output:
        os.path.join(
            "results/variant_counts",
            "wgs_phased_all_variants.hg19_multianno.lod100.tbx1_counts.txt"
        )
    shell: """
    perl src/perl/count_var_per_gene.pl \
    {input[1]} \
    {input[0]} \
    {output}
    """

rule checkpoint1:
    input:
        expand(
            os.path.join("results/annovar/output/counts",
                        ("wgs_phased_{VAR_TYPE}.sample.{SAMPLEID}."
                         "remaining22q11.filtering_summary.txt")), 
            VAR_TYPE = var_type, 
            SAMPLEID = probandid
        ),
        expand(
            os.path.join(
                "results/annovar/output/counts",
                ("wgs_phased_{VAR_TYPE}.sample.{SAMPLEID}."
                 "remaining22q11.exonic.nonsyn.deleterious.var_function.txt")
            ), 
            VAR_TYPE = var_type, 
            SAMPLEID = probandid
        ),
        expand(
            os.path.join(
                "results/annovar/output/counts",
                ("wgs_phased_{VAR_TYPE}.sample.{SAMPLEID}.remaining22q11."
                 "exonic.nonsyn.deleterious.exonic_function.txt")
            ),
            VAR_TYPE = var_type, 
            SAMPLEID = probandid
        ),
        expand(
            os.path.join(
                "results/annovar/output/counts",
                ("wgs_phased_{VAR_TYPE}.sample.{SAMPLEID}.freq.{FREQ}."
                 "filtering_summary.txt")
            ), 
            FREQ = var_freq, 
            VAR_TYPE = var_type, 
            SAMPLEID = sampleid
        ),
        expand(
            os.path.join(
                "results/annovar/output/counts",
                ("wgs_phased_{VAR_TYPE}.sample.{SAMPLEID}.freq.{FREQ}."
                 "nosuperdup.exonic.deleterious.var_function.txt")
            ), 
            FREQ = var_freq, 
            VAR_TYPE = var_type, 
            SAMPLEID = sampleid
        ),
        expand(
            os.path.join(
                "results/annovar/output/counts",
                ("wgs_phased_{VAR_TYPE}.sample.{SAMPLEID}.freq.{FREQ}."
                 "nosuperdup.exonic.deleterious.exonic_function.txt")
            ),
            FREQ = var_freq, 
            VAR_TYPE = var_type, 
            SAMPLEID = sampleid
        ),
        expand(
            os.path.join(
                "results/annovar/output/variant_filtering/{FREQ}",
                ("wgs_phased_{VAR_TYPE}.sample.{SAMPLEID}.freq.{FREQ}."
                 "nosuperdup.exonic.deleterious.txt")
            ), 
            FREQ = var_freq, 
            VAR_TYPE = var_type, 
            SAMPLEID = sampleid
        ),
        expand(
            os.path.join(
                "results/annovar/output/counts",
                "wgs_phased_{VAR_TYPE}.sample.{SAMPLEID}.pop_counts.txt"
            ),
            VAR_TYPE = var_type,
            SAMPLEID = sampleid
        ),
        expand(
            os.path.join(
                "results/annovar/output/counts/",
                "wgs_phased_{VAR_TYPE}.all_sample.pop_counts.txt"
            ),
            VAR_TYPE = var_type
        ),
        expand(
            os.path.join(
                "results/annovar/output/counts",
                ("wgs_phased_{VAR_TYPE}.all_sample.freq.{FREQ}."
                 "nosuperdup.exonic.deleterious.var_function.txt")
            ),
            VAR_TYPE = var_type,
            FREQ = var_freq
        ),
        expand(
            os.path.join(
                "results/annovar/output/counts",
                ("wgs_phased_{VAR_TYPE}.all_sample.freq.{FREQ}."
                 "nosuperdup.exonic.deleterious.exonic_function.txt")
            ),
            VAR_TYPE = var_type,
            FREQ = var_freq
        ),
        expand(
            os.path.join(
                "results/annovar/output/counts",
                ("wgs_phased_{VAR_TYPE}.all_sample.freq.{FREQ}."
                 "filtering_summary.txt")
            ),
            VAR_TYPE = var_type,
            FREQ = var_freq
        ),
        expand(
            os.path.join(
                "results/annovar/output/counts",
                "wgs_phased_{VAR_TYPE}.sample.{SAMPLE}.denovo.var_function.txt"
            ),
            VAR_TYPE = var_type,
            SAMPLE = probandid
        ),
        expand(
            os.path.join(
                "results/annovar/output/counts",
                ("wgs_phased_{VAR_TYPE}.sample.{SAMPLE}.denovo."
                 "exonic_function.txt")
            ),
            VAR_TYPE = var_type,
            SAMPLE = probandid
        ),
        expand(
            os.path.join(
                "results/annovar/output/counts",
                "wgs_phased_{VAR_TYPE}.sample.{SAMPLE}.denovo_pop_counts.txt"
            ),
            VAR_TYPE = var_type,
            SAMPLE = probandid
        ),
        expand(
            os.path.join(
                "results/annovar/output/counts",
                ("wgs_phased_{VAR_TYPE}.sample.{SAMPLE}.denovo."
                 "filtering_summary.txt")
            ),
            VAR_TYPE = var_type,
            SAMPLE = probandid
        )

################################################################################
# de novo variant extraction
################################################################################
rule denovo_all:
    input:
        expand(
            os.path.join(
                "results/annovar/output/counts",
                "wgs_phased_{VAR_TYPE}.sample.{SAMPLE}.denovo.var_function.txt"
            ),
            VAR_TYPE = var_type,
            SAMPLE = probandid
        ),
        expand(
            os.path.join(
                "results/annovar/output/counts",
                ("wgs_phased_{VAR_TYPE}.sample.{SAMPLE}.denovo."
                 "exonic_function.txt")
            ),
            VAR_TYPE = var_type,
            SAMPLE = probandid
        ),
        expand(
            os.path.join(
                "results/annovar/output/counts",
                "wgs_phased_{VAR_TYPE}.sample.{SAMPLE}.denovo_pop_counts.txt"
            ),
            VAR_TYPE = var_type,
            SAMPLE = probandid
        ),
        expand(
            os.path.join(
                "results/annovar/output/counts",
                ("wgs_phased_{VAR_TYPE}.sample.{SAMPLE}.denovo."
                 "filtering_summary.txt")
            ),
            VAR_TYPE = var_type,
            SAMPLE = probandid
        )

rule denovo_summarize_variant_function:
    input:
        os.path.join(
            "results/annovar/output/multiannotation/denovo/",
            "wgs_phased_{VAR_TYPE}.sample.{SAMPLE}.denovo.hg19_multianno.txt"
        )
    output:
        os.path.join(
            "results/annovar/output/counts",
            "wgs_phased_{VAR_TYPE}.sample.{SAMPLE}.denovo.var_function.txt"
        )
    run:
        summarize_variant_function(input[0], output[0])

rule denovo_summarize_exonic_function:
    input:
        os.path.join(
            "results/annovar/output/multiannotation/denovo/",
            "wgs_phased_{VAR_TYPE}.sample.{SAMPLE}.denovo.hg19_multianno.txt"
        )
    output:
        os.path.join(
            "results/annovar/output/counts",
            "wgs_phased_{VAR_TYPE}.sample.{SAMPLE}.denovo.exonic_function.txt"
        )
    run:
        summarize_exonic_functions(input[0], output[0])

rule denovo_population_counts:
    input: 
        os.path.join(
            "results/annovar/output/multiannotation/denovo/",
            "wgs_phased_{VAR_TYPE}.sample.{SAMPLE}.denovo.hg19_multianno.txt"
        )
    output:
        os.path.join(
            "results/annovar/output/counts",
            "wgs_phased_{VAR_TYPE}.sample.{SAMPLE}.denovo_pop_counts.txt"
        )
    message:"""
    Summarize counts of the different annotations.
    
    Get the counts of variants that have annotations for:
        
        * All variants
        * Found in dbSNP138
        * Found in custom dbsnp138
        * Found in 1kg 
        * Found in ESP6500
        * Found in any DB
        * Novel variants
    
    Input: {input}
    Output: {output}
    """
    run:
        count_population_database(input[0], output[0])

rule denovo_summarize_filtering_steps:
    input:
        multianno_file = os.path.join(
            "results/annovar/output/multiannotation/indiv_samples/{SAMPLE}",
            "wgs_phased_{VAR_TYPE}.sample.{SAMPLE}.hg19_multianno.txt"
        ),
        var_filter_file = os.path.join(
            "results/annovar/output/multiannotation/denovo",
            "wgs_phased_{VAR_TYPE}.sample.{SAMPLE}.denovo.hg19_multianno.txt"
        ),
        superdup_file = os.path.join(
            "results/annovar/output/variant_filtering/denovo",
            "wgs_phased_{VAR_TYPE}.sample.{SAMPLE}.denovo.nosuperdup.txt"
        ),
        noncoding_file = os.path.join(
            "results/annovar/output/variant_filtering/denovo",
            ("wgs_phased_{VAR_TYPE}.sample.{SAMPLE}.denovo.nosuperdup."
             "noncoding.txt")
        ),
        exonic_file = os.path.join(
            "results/annovar/output/variant_filtering/denovo",
            ("wgs_phased_{VAR_TYPE}.sample.{SAMPLE}.denovo.nosuperdup."
             "exonic.txt")
        ),
        nonsynonymous_file = os.path.join(
            "results/annovar/output/variant_filtering/denovo",
            ("wgs_phased_{VAR_TYPE}.sample.{SAMPLE}.denovo.nosuperdup."
             "exonic.nonsyn.txt")
        ),
        deleterious_file = os.path.join(
            "results/annovar/output/variant_filtering/denovo",
            ("wgs_phased_{VAR_TYPE}.sample.{SAMPLE}.denovo.nosuperdup."
             "exonic.deleterious.txt")
        )
    params:
        var_filter_type = "Allele frequency: de novo"
    output:
        os.path.join(
            "results/annovar/output/counts",
            ("wgs_phased_{VAR_TYPE}.sample.{SAMPLE}.denovo."
             "filtering_summary.txt")
        )
    run:
        summarize_filtering_results(
            multianno_file = input.multianno_file,
            var_filter_file = input.var_filter_file,
            superdup_file = input.superdup_file,
            noncoding_file = input.noncoding_file,
            exonic_file = input.exonic_file,
            nonsynonymous_file = input.nonsynonymous_file,
            deleterious_file = input.deleterious_file,
            output = output[0],
            var_filter_type = params.var_filter_type
        )

rule denovo_extract_deleterious:
    input:
        exonic = os.path.join(
            "results/annovar/output/variant_filtering/denovo",
            ("wgs_phased_{VAR_TYPE}.sample.{SAMPLE}.denovo.nosuperdup."
             "exonic.nonsyn.txt")
        )
    params:
        deleterious_count_cutoff = "2"
    output:
        exonic = os.path.join(
            "results/annovar/output/variant_filtering/denovo",
            ("wgs_phased_{VAR_TYPE}.sample.{SAMPLE}.denovo.nosuperdup."
             "exonic.deleterious.txt")
        )
    run:
        extract_deleterious_variants(input[0], output[0], 
                                     params.deleterious_count_cutoff)
    
rule denovo_extract_nonsynonymous:
    input:
        exonic = os.path.join(
            "results/annovar/output/variant_filtering/denovo",
            ("wgs_phased_{VAR_TYPE}.sample.{SAMPLE}.denovo.nosuperdup."
             "exonic.txt")
        ),
    output:
        exonic = os.path.join(
            "results/annovar/output/variant_filtering/denovo",
            ("wgs_phased_{VAR_TYPE}.sample.{SAMPLE}.denovo.nosuperdup."
             "exonic.nonsyn.txt")
        )
    run:
        extract_nonsynonymous_variants(input[0], output[0])

rule denovo_parse_coding_noncoding:
    input:
        os.path.join(
            "results/annovar/output/variant_filtering/denovo",
            "wgs_phased_{VAR_TYPE}.sample.{SAMPLE}.denovo.nosuperdup.txt"
        )
    output:
        exonic = os.path.join(
            "results/annovar/output/variant_filtering/denovo",
            ("wgs_phased_{VAR_TYPE}.sample.{SAMPLE}.denovo.nosuperdup."
             "exonic.txt")
        ),
        noncoding = os.path.join(
            "results/annovar/output/variant_filtering/denovo",
            ("wgs_phased_{VAR_TYPE}.sample.{SAMPLE}.denovo.nosuperdup."
             "noncoding.txt")
        )
    message:"""
    Split input file into coding and noncoding variants.
    
    Input: {input}
    Output: {output}
    """
    run:
        parse_coding_vs_noncoding(input[0], output.exonic, output.noncoding)

rule denovo_remove_superdup:
    input:
        os.path.join(
            "results/annovar/output/multiannotation/denovo",
            "wgs_phased_{VAR_TYPE}.sample.{SAMPLE}.denovo.hg19_multianno.txt"
        )
    output:
        os.path.join(
            "results/annovar/output/variant_filtering/denovo",
            "wgs_phased_{VAR_TYPE}.sample.{SAMPLE}.denovo.nosuperdup.txt"
        )
    run:
        remove_superdup(input[0], output[0])
    
rule denovo_table_annotation:
    input: 
        os.path.join(
            "results/annovar/input/denovo",
            "wgs_phased_{VAR_TYPE}.sample.{SAMPLE}.denovo.recode.avinput"
        ),
        "/cork/jchung/annovar/humandb/hg19_custom_dbsnp138.txt"
        
    params: 
        output_prefix = os.path.join(
            "results/annovar/output/multiannotation/denovo",
            "wgs_phased_{VAR_TYPE}.sample.{SAMPLE}.denovo"),
        protocol = ",".join(["refGene", "phastConsElements46way", 
                             "genomicSuperDups", "esp6500si_all", 
                             "1000g2012apr_all", "snp138", "generic",
                             "ljb23_all"]),
        operation = "g,r,r,f,f,f,f,f",
        build = "hg19",
        custom_db = "hg19_custom_dbsnp138.txt"
    output: 
        os.path.join(
            "results/annovar/output/multiannotation/denovo",
            "wgs_phased_{VAR_TYPE}.sample.{SAMPLE}.denovo.hg19_multianno.txt"
        )
    shell: """
    {annovar_dir}/table_annovar.pl \
    {input[0]} \
    {annovar_dir}/humandb/ \
    --protocol {params.protocol} \
    --operation {params.operation} \
    --buildver {params.build} \
    --genericdbfile {params.custom_db} \
    --outfile {params.output_prefix} \
    --otherinfo
    """

rule convert_denovo_to_avinput:
    input:
        os.path.join("results/annovar/input/denovo",
                     "wgs_phased_{VAR_TYPE}.sample.{SAMPLE}.denovo.recode.vcf")
    output:
        os.path.join(
            "results/annovar/input/denovo",
            "wgs_phased_{VAR_TYPE}.sample.{SAMPLE}.denovo.recode.avinput")
    shell: """
        {annovar_dir}/convert2annovar.pl \
        -format vcf4 {input} \
        -outfile {output} \
        -includeinfo \
        -filter PASS \
        -comment
    """

rule split_denovo_vcf_file:
    input:
        os.path.join("results/annovar/input/denovo",
                     "wgs_phased_{VAR_TYPE}.sample.{SAMPLE}.denovo.vcf")
    params:
        subject = "{SAMPLE}",
        output_prefix = os.path.join(
            "results/annovar/input/denovo",
            "wgs_phased_{VAR_TYPE}.sample.{SAMPLE}.denovo"
        )
    output: 
        os.path.join("results/annovar/input/denovo",
                     "wgs_phased_{VAR_TYPE}.sample.{SAMPLE}.denovo.recode.vcf")
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
        os.path.join("results/annovar/input/denovo",
                     "wgs_phased_{VAR_TYPE}.sample.BM1453.001.denovo.vcf")
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
        os.path.join("results/annovar/input/denovo",
                     "wgs_phased_{VAR_TYPE}.sample.BM1452.001.denovo.vcf")
    run:
        extract_denovo_variants(input[0], 
                                output[0], 
                                params.proband, 
                                params.father, 
                                params.mother)

################################################################################
# Remaining allele of chr22q11
################################################################################
rule check_22q11_genes:
    input:
        expand(
            os.path.join(
                "results/annovar/output/counts",
                ("wgs_phased_{VAR_TYPE}.sample.{SAMPLEID}."
                 "remaining22q11.filtering_summary.txt")
            ), 
            VAR_TYPE = var_type, SAMPLEID = probandid),
        expand(
            os.path.join(
                "results/annovar/output/counts",
                ("wgs_phased_{VAR_TYPE}.sample.{SAMPLEID}."
                 "remaining22q11.exonic.nonsyn.deleterious.var_function.txt")
            ), 
            VAR_TYPE = var_type, SAMPLEID = probandid),
        expand(
            os.path.join(
                "results/annovar/output/counts",
                ("wgs_phased_{VAR_TYPE}.sample.{SAMPLEID}.remaining22q11."
                 "exonic.nonsyn.deleterious.exonic_function.txt")
            ),
            VAR_TYPE = var_type, SAMPLEID = probandid)

rule chr22q11_summarize_variant_function:
    input:
        os.path.join(
            "results/annovar/output/variant_filtering/remaining_22q11",
            ("wgs_phased_{VAR_TYPE}.sample.{SAMPLEID}." 
             "remaining22q11.exonic.nonsyn.deleterious.txt")
        )
    output:
        os.path.join(
                "results/annovar/output/counts",
                ("wgs_phased_{VAR_TYPE}.sample.{SAMPLEID}."
                 "remaining22q11.exonic.nonsyn.deleterious.var_function.txt")
            )
    run:
        summarize_variant_function(input[0], output[0])

rule chr22q11_summarize_exonic_function:
    input:
        os.path.join(
            "results/annovar/output/variant_filtering/remaining_22q11",
            ("wgs_phased_{VAR_TYPE}.sample.{SAMPLEID}." 
             "remaining22q11.exonic.nonsyn.deleterious.txt")
        )
    output:
        os.path.join(
                "results/annovar/output/counts",
                ("wgs_phased_{VAR_TYPE}.sample.{SAMPLEID}.remaining22q11."
                 "exonic.nonsyn.deleterious.exonic_function.txt")
            )
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
    params:
        var_filter_type = "Variants on chr22q11.2"
    output:
        os.path.join(
            "results/annovar/output/counts",
            ("wgs_phased_{VAR_TYPE}.sample.{SAMPLEID}."
             "remaining22q11.filtering_summary.txt")
        )
    run:
        summarize_filtering_results(
            multianno_file = input.multianno_file,
            var_filter_file = input.chr22_file,
            superdup_file = "None",
            noncoding_file = input.noncoding_file,
            exonic_file = input.exonic_file,
            nonsynonymous_file = input.nonsynonymous_file,
            deleterious_file = input.deleterious_file,
            output = output[0],
            var_filter_type = params.var_filter_type
        )

rule chr22q11_extract_deleterious:
    input:
        os.path.join(
            "results/annovar/output/variant_filtering/remaining_22q11",
            ("wgs_phased_{VAR_TYPE}.sample.{SAMPLEID}." 
             "remaining22q11.exonic.nonsyn.txt")
        )
    params: 
        deleterious_cutoff = "2"
    output: 
        os.path.join(
            "results/annovar/output/variant_filtering/remaining_22q11",
            ("wgs_phased_{VAR_TYPE}.sample.{SAMPLEID}." 
             "remaining22q11.exonic.nonsyn.deleterious.txt")
        )
    run:
        extract_deleterious_variants(input[0], output[0], 
                                     params.deleterious_cutoff)

rule chr22q11_extract_nonsynonymous:
    input: 
        os.path.join(
            "results/annovar/output/variant_filtering/remaining_22q11",
            ("wgs_phased_{VAR_TYPE}.sample.{SAMPLEID}." 
             "remaining22q11.exonic.txt")
        )
    output: 
        os.path.join(
            "results/annovar/output/variant_filtering/remaining_22q11",
            ("wgs_phased_{VAR_TYPE}.sample.{SAMPLEID}." 
             "remaining22q11.exonic.nonsyn.txt")
        )
    run:
        extract_nonsynonymous_variants(input[0], output[0])

rule chr22q11_parse_coding_noncoding:
    input:
        os.path.join(
            "results/annovar/output/variant_filtering/remaining_22q11",
            "wgs_phased_{VAR_TYPE}.sample.{SAMPLEID}.remaining22q11.txt"
        )
    output:
        exonic = os.path.join(
            "results/annovar/output/variant_filtering/remaining_22q11",
            ("wgs_phased_{VAR_TYPE}.sample.{SAMPLEID}." 
             "remaining22q11.exonic.txt")
        ),
        noncoding = os.path.join(
            "results/annovar/output/variant_filtering/remaining_22q11",
            ("wgs_phased_{VAR_TYPE}.sample.{SAMPLEID}." 
             "remaining22q11.noncoding.txt")
        )
    message:"""
    Parse exonic and noncoding variants.
    
        Input: {input}
        Output: {output}
    """
    run:
        parse_coding_vs_noncoding(input[0], output.exonic, output.noncoding)

rule extract_22q11_genes:
    input:
        os.path.join(
            "results/annovar/output/multiannotation/indiv_samples/{SAMPLEID}",
            "wgs_phased_{VAR_TYPE}.sample.{SAMPLEID}.hg19_multianno.txt"
        )
    params:
        chromosome = "chr22",
        start = "18656000",
        stop = "21792000"
    output: 
        os.path.join(
            "results/annovar/output/variant_filtering/remaining_22q11",
            "wgs_phased_{VAR_TYPE}.sample.{SAMPLEID}.remaining22q11.txt"
        )
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
        expand(
            os.path.join(
                "results/annovar/output/counts",
                ("wgs_phased_{VAR_TYPE}.sample.{SAMPLEID}.freq.{FREQ}."
                 "filtering_summary.txt")
            ), 
            FREQ = var_freq, 
            VAR_TYPE = var_type, 
            SAMPLEID = sampleid
        ),
        expand(
            os.path.join(
                "results/annovar/output/counts",
                ("wgs_phased_{VAR_TYPE}.sample.{SAMPLEID}.freq.{FREQ}."
                 "nosuperdup.exonic.deleterious.var_function.txt")
            ), 
            FREQ = var_freq, 
            VAR_TYPE = var_type, 
            SAMPLEID = sampleid
        ),
        expand(
            os.path.join(
                "results/annovar/output/variant_filtering/{FREQ}",
                ("wgs_phased_{VAR_TYPE}.sample.{SAMPLEID}.freq.{FREQ}."
                 "nosuperdup.exonic.deleterious.txt")
            ), 
            FREQ = var_freq, 
            VAR_TYPE = var_type, 
            SAMPLEID = sampleid
        )

rule indiv_summarize_variant_function:
    input:
        os.path.join(
            "results/annovar/output/variant_filtering/{FREQ}",
            ("wgs_phased_{VAR_TYPE}.sample.{SAMPLEID}.freq.{FREQ}."
             "nosuperdup.exonic.deleterious.txt")
        )
    output:
        os.path.join(
            "results/annovar/output/counts",
            ("wgs_phased_{VAR_TYPE}.sample.{SAMPLEID}.freq.{FREQ}."
             "nosuperdup.exonic.deleterious.var_function.txt")
        )
    run:
        summarize_variant_function(input[0], output[0])

rule indiv_summarize_exonic_function:
    input:
        os.path.join(
            "results/annovar/output/variant_filtering/{FREQ}",
            ("wgs_phased_{VAR_TYPE}.sample.{SAMPLEID}.freq.{FREQ}."
             "nosuperdup.exonic.deleterious.txt")
        )
    output:
        os.path.join(
            "results/annovar/output/counts",
            ("wgs_phased_{VAR_TYPE}.sample.{SAMPLEID}.freq.{FREQ}."
             "nosuperdup.exonic.deleterious.exonic_function.txt")
        )
    run:
        summarize_exonic_functions(input[0], output[0])

rule indiv_summarize_filtering_steps:
    input:
        multianno_file = os.path.join(
            "results/annovar/output/multiannotation/indiv_samples/{SAMPLEID}",
            "wgs_phased_{VAR_TYPE}.sample.{SAMPLEID}.hg19_multianno.txt"
        ),
        var_filter_file = os.path.join(
            "results/annovar/output/variant_filtering/{FREQ}",
            "wgs_phased_{VAR_TYPE}.sample.{SAMPLEID}.freq.{FREQ}.txt"
        ),
        superdup_file = os.path.join(
            "results/annovar/output/variant_filtering/{FREQ}",
            ("wgs_phased_{VAR_TYPE}.sample.{SAMPLEID}.freq.{FREQ}."
             "nosuperdup.txt")
        ),
        noncoding_file = os.path.join(
            "results/annovar/output/variant_filtering/{FREQ}",
            ("wgs_phased_{VAR_TYPE}.sample.{SAMPLEID}.freq.{FREQ}."
             "nosuperdup.noncoding.txt")
        ),
        exonic_file = os.path.join(
            "results/annovar/output/variant_filtering/{FREQ}",
            ("wgs_phased_{VAR_TYPE}.sample.{SAMPLEID}.freq.{FREQ}."
             "nosuperdup.exonic.txt")
        ),
        nonsynonymous_file = os.path.join(
            "results/annovar/output/variant_filtering/{FREQ}",
            ("wgs_phased_{VAR_TYPE}.sample.{SAMPLEID}.freq.{FREQ}."
             "nosuperdup.exonic.nonsyn.txt")
        ),
        deleterious_file = os.path.join(
            "results/annovar/output/variant_filtering/{FREQ}",
            ("wgs_phased_{VAR_TYPE}.sample.{SAMPLEID}.freq.{FREQ}."
             "nosuperdup.exonic.deleterious.txt")
        )
    params:
        var_filter_type = "Allele frequency: {FREQ}"
    output:
        os.path.join(
            "results/annovar/output/counts",
            ("wgs_phased_{VAR_TYPE}.sample.{SAMPLEID}.freq.{FREQ}."
             "filtering_summary.txt")
        )
    run:
        summarize_filtering_results(
            multianno_file = input.multianno_file,
            var_filter_file = input.var_filter_file,
            superdup_file = input.superdup_file,
            noncoding_file = input.noncoding_file,
            exonic_file = input.exonic_file,
            nonsynonymous_file = input.nonsynonymous_file,
            deleterious_file = input.deleterious_file,
            output = output[0],
            var_filter_type = params.var_filter_type
        )

rule indiv_extract_deleterious:
    input:
        os.path.join(
            "results/annovar/output/variant_filtering/{FREQ}",
            ("wgs_phased_{VAR_TYPE}.sample.{SAMPLEID}.freq.{FREQ}."
             "nosuperdup.exonic.nonsyn.txt")
        )
    params:
        deleterious_count_cutoff = "2"
    output:
        os.path.join(
            "results/annovar/output/variant_filtering/{FREQ}",
            ("wgs_phased_{VAR_TYPE}.sample.{SAMPLEID}.freq.{FREQ}."
             "nosuperdup.exonic.deleterious.txt")
        )
    message: """
    Extract deleterious exonic variants.
    
    Input: {input}
    Output: {output}
    """
    run:
        extract_deleterious_variants(input[0], output[0], 
                                     params.deleterious_count_cutoff)

rule indiv_extract_nonsynonymous:
    input:
        os.path.join(
            "results/annovar/output/variant_filtering/{FREQ}",
            ("wgs_phased_{VAR_TYPE}.sample.{SAMPLEID}.freq.{FREQ}."
             "nosuperdup.exonic.txt")
        )
    output:
        os.path.join(
            "results/annovar/output/variant_filtering/{FREQ}",
            ("wgs_phased_{VAR_TYPE}.sample.{SAMPLEID}.freq.{FREQ}."
             "nosuperdup.exonic.nonsyn.txt")
        )
    message:"""
    Extract non-synonymous variants.
    
    Input: {input}
    Output: {output}
    """
    run:
        extract_nonsynonymous_variants(input[0], output[0])

rule indiv_parse_coding_noncoding:
    input:
        os.path.join(
            "results/annovar/output/variant_filtering/{FREQ}",
            ("wgs_phased_{VAR_TYPE}.sample.{SAMPLEID}.freq.{FREQ}."
             "nosuperdup.txt")
        )
    output:
        exonic = os.path.join(
            "results/annovar/output/variant_filtering/{FREQ}",
            ("wgs_phased_{VAR_TYPE}.sample.{SAMPLEID}.freq.{FREQ}."
             "nosuperdup.exonic.txt")
        ),
        noncoding = os.path.join(
            "results/annovar/output/variant_filtering/{FREQ}",
            ("wgs_phased_{VAR_TYPE}.sample.{SAMPLEID}.freq.{FREQ}."
             "nosuperdup.noncoding.txt")
        )
    message:"""
    Split input file into coding and noncoding variants.
    
    Input: {input}
    Output: {output}
    """
    run:
        parse_coding_vs_noncoding(input[0], output.exonic, output.noncoding)

rule indiv_remove_superdup:
    input:
        os.path.join(
            "results/annovar/output/variant_filtering/{FREQ}",
            "wgs_phased_{VAR_TYPE}.sample.{SAMPLEID}.freq.{FREQ}.txt"
        )
    output: 
        os.path.join(
            "results/annovar/output/variant_filtering/{FREQ}",
            "wgs_phased_{VAR_TYPE}.sample.{SAMPLEID}.freq.{FREQ}.nosuperdup.txt"
        )
    message:"""
    Removing variants found in genomic superdups.
    
    Input: {input}
    Output: {output}
    """
    run:
        remove_superdup(input[0], output[0])

rule indiv_filter_on_variant_frequency:
    input:
        os.path.join(
            "results/annovar/output/multiannotation/indiv_samples/{SAMPLEID}",
            "wgs_phased_{VAR_TYPE}.sample.{SAMPLEID}.hg19_multianno.txt"
        )
    params:
        frequency = "{FREQ}"
    output:
        os.path.join(
            "results/annovar/output/variant_filtering/{FREQ}",
            "wgs_phased_{VAR_TYPE}.sample.{SAMPLEID}.freq.{FREQ}.txt"
        )
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
        expand(
            os.path.join(
                "results/annovar/output/counts",
                "wgs_phased_{VAR_TYPE}.sample.{SAMPLEID}.pop_counts.txt"
            ), 
            VAR_TYPE = var_type, 
            SAMPLEID = sampleid
        ),
        expand(
            os.path.join(
                "results/annovar/output/counts",
                "wgs_phased_{VAR_TYPE}.all_sample.pop_counts.txt",
            ),
            VAR_TYPE = var_type)
    output:
        snp_output = os.path.join("results/annovar/output/counts",
                                  "wgs_phased_snp.all_samples.pop_counts.txt"),
        indel_ouptut = os.path.join(
                           "results/annovar/output/counts",
                           "wgs_phased_indel.all_samples.pop_counts.txt"
                       )
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
        os.path.join(
            "results/annovar/output/multiannotation/indiv_samples/{SAMPLEID}",
            "wgs_phased_{VAR_TYPE}.sample.{SAMPLEID}.hg19_multianno.txt"
        )
    output:
        os.path.join(
            "results/annovar/output/counts",
            "wgs_phased_{VAR_TYPE}.sample.{SAMPLEID}.pop_counts.txt"
        )
    message:"""
    Summarize counts of the different annotations.
    
    Get the counts of variants that have annotations for:
        
        * All variants
        * Found in dbSNP138
        * Found in custom dbsnp138
        * Found in 1kg 
        * Found in ESP6500
        * Found in any DB
        * Novel variants
    
    Input: {input}
    Output: {output}
    """
    run:
        count_population_database(input[0], output[0])

rule table_annotation_indiv:
    input: 
        os.path.join("results/annovar/input",
                     "qc_wgs_phased_{VAR_TYPE}.sample.{SAMPLEID}.avinput"),
        "/cork/jchung/annovar/humandb/hg19_custom_dbsnp138.txt"
    params: 
        output_prefix = os.path.join(
            "results/annovar/output/multiannotation/indiv_samples/{SAMPLEID}",
            "wgs_phased_{VAR_TYPE}.sample.{SAMPLEID}"),
        protocol = ",".join(["refGene", "phastConsElements46way", 
                             "genomicSuperDups", "esp6500si_all", 
                             "1000g2012apr_all", "snp138", "generic",
                             "ljb23_all"]),
        operation = "g,r,r,f,f,f,f,f",
        build = "hg19",
        custom_db = "hg19_custom_dbsnp138.txt"
    output: 
        os.path.join(
            "results/annovar/output/multiannotation/indiv_samples/{SAMPLEID}",
            "wgs_phased_{VAR_TYPE}.sample.{SAMPLEID}.hg19_multianno.txt"
        )
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
    {input[0]} \
    {annovar_dir}/humandb/ \
    --protocol {params.protocol} \
    --operation {params.operation} \
    --buildver {params.build} \
    --hgvs \
    --genericdbfile {params.custom_db} \
    --outfile {params.output_prefix} \
    --otherinfo
    """

################################################################################
# All sample CADD table annotations
################################################################################
# rule filter_cadd:
    # input:
        # os.path.join("results/annovar_output/cadd_annotation",
                     # "wgs_phased_snp.cadd.hg19_cadd_dropped")
    # params: 
        # scaled_cadd_cutoff = "20"
    # output: 
        # os.path.join("results/annovar_output/cadd_annotation",
                     # "wgs_phased_snp.cadd.hg19_cadd_dropped.filtered")
    # run:
        # in_file = open(input[0], "r")
        # out_file = open(output[0], "w")

        # for line in in_file:
            # split_line = line.split("\t")
            # cadd_score = split_line[1].split(",")
            # scaled_cadd = cadd_score[1]
            # if float(scaled_cadd) >= int(params.scaled_cadd_cutoff):
                # out_file.write(line)
        # in_file.close()
        # out_file.close()

# rule annotate_novel_cadd:
    # input: 
        # os.path.join("results/indiv_avinput",
                     # "qc_wgs_phased_snp.sample.BM1452.001.avinput")
    # params: 
        # output_prefix = os.path.join("results/annovar_output/cadd_annotation",
                                     # "wgs_phased_snp.cadd")
    # output: 
        # os.path.join("results/annovar_output/cadd_annotation",
                     # "wgs_phased_snp.cadd.hg19_cadd_dropped")
    # shell: """
        # $HOME/programs/annovar/annovar2014Jul14/annotate_variation.pl \
        # {input} \
        # $HOME/programs/annovar/annovar2014Jul14/humandb/ \
        # -filter \
        # -dbtype cadd \
        # -buildver hg19 \
        # -hgvs \
        # -out {params.output_prefix} \
        # -otherinfo
    # """

################################################################################
# Prepare vcf files for individual samples
################################################################################
rule convert_to_annovar_indiv:
    input: 
        os.path.join("results/annovar/input",
                     "qc_wgs_phased_{VAR_TYPE}.sample.{SAMPLEID}.recode.vcf")
    output: 
        os.path.join("results/annovar/input",
                     "qc_wgs_phased_{VAR_TYPE}.sample.{SAMPLEID}.avinput")
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
        os.path.join("results/annovar/input",
                     "wgs_phased_{VAR_TYPE}.sample.{SAMPLEID}.recode.vcf")
    params: 
        sample = "{SAMPLEID}"
    output: 
        os.path.join("results/annovar/input",
                     "qc_wgs_phased_{VAR_TYPE}.sample.{SAMPLEID}.recode.vcf")
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
        output_prefix = os.path.join("results/annovar/input",
                                     "wgs_phased_{VAR_TYPE}.sample.{SAMPLEID}")
    output: 
        os.path.join("results/annovar/input",
                     "wgs_phased_{VAR_TYPE}.sample.{SAMPLEID}.recode.vcf")
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
            os.path.join("results/annovar/output/counts/",
                     "wgs_phased_{VAR_TYPE}.all_sample.pop_counts.txt"), 
            VAR_TYPE = var_type
        )

rule all_sample_population_counts:
    input: 
        os.path.join("results/annovar/output/multiannotation/all_samples",
                     "wgs_phased_{VAR_TYPE}.hg19_multianno.txt")
    output: 
        os.path.join("results/annovar/output/counts/",
                     "wgs_phased_{VAR_TYPE}.all_sample.pop_counts.txt")
    message:"""
    Summarize counts of the different annotations.

    Get the counts of variants that were previously reported:
        
        * All variants
        * Found in dbSNP138
        * Found in generic dbsnp138
        * Found in 1kg 
        * Found in ESP6500
        * In any DB
        * Novel variants
    
    Input: {input}
    Output: {output}
    """
    run:
        count_population_database(input[0], output[0])

rule all_sample_summarize_variant_function:
    input:
        os.path.join(
            "results/annovar/output/variant_filtering/all_samples/{FREQ}",
            ("wgs_phased_{VAR_TYPE}.all_sample.freq.{FREQ}.nosuperdup."
             "exonic.deleterious.txt")
        )
    output:
        os.path.join(
            "results/annovar/output/counts",
            ("wgs_phased_{VAR_TYPE}.all_sample.freq.{FREQ}."
             "nosuperdup.exonic.deleterious.var_function.txt")
        )
    run:
        summarize_variant_function(input[0], output[0])

rule all_sample_summarize_exonic_function:
    input:
        os.path.join(
            "results/annovar/output/variant_filtering/all_samples/{FREQ}",
            ("wgs_phased_{VAR_TYPE}.all_sample.freq.{FREQ}.nosuperdup."
             "exonic.deleterious.txt")
        )
    output:
        os.path.join(
            "results/annovar/output/counts",
            ("wgs_phased_{VAR_TYPE}.all_sample.freq.{FREQ}."
             "nosuperdup.exonic.deleterious.exonic_function.txt")
        )
    run:
        summarize_exonic_functions(input[0], output[0])

rule all_sample_summarize_filtering_steps:
    input:
        multianno_file = os.path.join(
            "results/annovar/output/multiannotation/all_samples/",
            "wgs_phased_{VAR_TYPE}.hg19_multianno.txt"
        ),
        var_filter_file = os.path.join(
            "results/annovar/output/variant_filtering/all_samples/{FREQ}",
            "wgs_phased_{VAR_TYPE}.all_sample.freq.{FREQ}.txt"
        ),
        superdup_file = os.path.join(
            "results/annovar/output/variant_filtering/all_samples/",
            "{FREQ}/wgs_phased_{VAR_TYPE}.all_sample.freq.{FREQ}.nosuperdup.txt"
        ),
        noncoding_file = os.path.join(
            "results/annovar/output/variant_filtering/all_samples/{FREQ}",
            ("wgs_phased_{VAR_TYPE}.all_sample.freq.{FREQ}.nosuperdup."
             "noncoding.txt")
        ),
        exonic_file = os.path.join(
            "results/annovar/output/variant_filtering/all_samples/{FREQ}",
            "wgs_phased_{VAR_TYPE}.all_sample.freq.{FREQ}.nosuperdup.exonic.txt"
        ),
        nonsynonymous_file = os.path.join(
            "results/annovar/output/variant_filtering/all_samples/{FREQ}",
            ("wgs_phased_{VAR_TYPE}.all_sample.freq.{FREQ}.nosuperdup."
             "exonic.nonsyn.txt")
        ),
        deleterious_file = os.path.join(
            "results/annovar/output/variant_filtering/all_samples/{FREQ}",
            ("wgs_phased_{VAR_TYPE}.all_sample.freq.{FREQ}.nosuperdup."
             "exonic.deleterious.txt")
        )
    params:
        var_filter_type = "Allele frequency: {FREQ}"
    output:
        os.path.join(
            "results/annovar/output/counts",
            ("wgs_phased_{VAR_TYPE}.all_sample.freq.{FREQ}."
             "filtering_summary.txt")
        )
    run:
        summarize_filtering_results(
            multianno_file = input.multianno_file,
            var_filter_file = input.var_filter_file,
            superdup_file = input.superdup_file,
            noncoding_file = input.noncoding_file,
            exonic_file = input.exonic_file,
            nonsynonymous_file = input.nonsynonymous_file,
            deleterious_file = input.deleterious_file,
            output = output[0],
            var_filter_type = params.var_filter_type
        )

rule all_sample_extract_deleterious:
    input:
        os.path.join(
            "results/annovar/output/variant_filtering/all_samples/{FREQ}",
            ("wgs_phased_{VAR_TYPE}.freq.{FREQ}.nosuperdup."
             "exonic.nonsyn.txt")
        )
    params:
        deleterious_count_cutoff = "2"
    output:
        os.path.join(
            "results/annovar/output/variant_filtering/all_samples/{FREQ}",
            ("wgs_phased_{VAR_TYPE}.freq.{FREQ}.nosuperdup."
             "exonic.deleterious.txt")
        )
    run:
        extract_deleterious_variants(input[0], output[0], 
                                     params.deleterious_count_cutoff)

rule all_sample_extract_nonsynonymous:
    input:
        os.path.join(
            "results/annovar/output/variant_filtering/all_samples/{FREQ}",
            ("wgs_phased_{VAR_TYPE}.all_sample.freq.{FREQ}.nosuperdup."
             "exonic.txt")
        )
    output:
        os.path.join(
            "results/annovar/output/variant_filtering/all_samples/{FREQ}",
            ("wgs_phased_{VAR_TYPE}.all_sample.freq.{FREQ}.nosuperdup."
             "exonic.nonsyn.txt")
        )
    run:
        extract_nonsynonymous_variants(input[0], output[0])

rule all_sample_parse_coding_noncoding:
    input:
        os.path.join(
            "results/annovar/output/variant_filtering/all_samples/{FREQ}",
            "wgs_phased_{VAR_TYPE}.all_sample.freq.{FREQ}.nosuperdup.txt"
        )
    output:
        exonic = os.path.join(
            "results/annovar/output/variant_filtering/all_samples/{FREQ}",
            "wgs_phased_{VAR_TYPE}.all_sample.freq.{FREQ}.nosuperdup.exonic.txt"
        ),
        noncoding = os.path.join(
            "results/annovar/output/variant_filtering/all_samples/{FREQ}",
            ("wgs_phased_{VAR_TYPE}.all_sample.freq.{FREQ}.nosuperdup."
             "noncoding.txt")
        )
    run:
        parse_coding_vs_noncoding(input[0], output.exonic, output.noncoding)

rule all_sample_remove_superdup:
    input:
        os.path.join(
            "results/annovar/output/variant_filtering/all_samples/{FREQ}",
            "wgs_phased_{VAR_TYPE}.all_sample.freq.{FREQ}.txt")
    output:
        os.path.join(
            "results/annovar/output/variant_filtering/all_samples/{FREQ}",
            "wgs_phased_{VAR_TYPE}.all_sample.freq.{FREQ}.nosuperdup.txt")
    run:
        remove_superdup(input[0], output[0])

rule all_sample_filter_on_variant_frequency:
    input:
        os.path.join(
            "results/annovar/output/multiannotation/all_samples/",
            "wgs_phased_{VAR_TYPE}.hg19_multianno.txt")
    params:
        frequency = "{FREQ}"
    output:
        os.path.join(
            "results/annovar/output/variant_filtering/all_samples/{FREQ}",
            "wgs_phased_{VAR_TYPE}.all_sample.freq.{FREQ}.txt")
    run:
        if params.frequency == "all":
            shell("cp {input} {output}")
        elif params.frequency == "novel":
            extract_novel_variants(input[0], output[0])
        else:
            maf = float(params.frequency[4:])
            extract_rare_variants(input[0], output[0], maf)

rule table_annotation_all_samples:
    input: 
        "results/annovar/input/allsample/wgs_phased_{VAR_TYPE}.avinput",
        "/cork/jchung/annovar/humandb/hg19_custom_dbsnp138.txt"
    params: 
        output_prefix = os.path.join(
            "results/annovar/output/multiannotation/all_samples",
            "wgs_phased_{VAR_TYPE}"),
        protocol = ",".join(["refGene", "phastConsElements46way", 
                             "genomicSuperDups", "esp6500si_all", 
                             "1000g2012apr_all", "snp138", "generic",
                             "ljb23_all"]),
        operation = "g,r,r,f,f,f,f,f",
        build = "hg19",
        custom_db = "hg19_custom_dbsnp138.txt"
    output:
        os.path.join("results/annovar/output/multiannotation/all_samples/",
                     "wgs_phased_{VAR_TYPE}.hg19_multianno.txt")
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
    {input[0]} \
    {annovar_dir}/humandb/ \
    --protocol {params.protocol} \
    --operation {params.operation} \
    --buildver {params.build} \
    --genericdbfile {params.custom_db} \
    --outfile {params.output_prefix} \
    --otherinfo
    """

rule refgene_annotation_all_samples:
    input: 
        "results/annovar/input/allsample/wgs_phased_{VAR_TYPE}.avinput"
    params: 
        output_prefix = os.path.join(
            "results/annovar/output/refgene/all_samples",
            "wgs_phased_{VAR_TYPE}.refgene_hgvs"),
        build = "hg19"
    output:
        os.path.join("results/annovar/output/refgene/all_samples/",
                     "wgs_phased_{VAR_TYPE}.refgene_hgvs.variant_function"),
        os.path.join("results/annovar/output/refgene/all_samples/",
                     "wgs_phased_{VAR_TYPE}.refgene_hgvs.exonic_variant_function")
    message: """
    Perform gene-based annotation using annovar and output in HGVS format.
    
    Annotate variants using:
        
        * refGene

    Input: {input}
    Output: {output}
    """
    shell: """
    {annovar_dir}/annotate_variation.pl \
    --geneanno \
    --outfile {params.output_prefix} \
    --buildver {params.build} \
    --hgvs \
    --otherinfo \
    {input[0]} \
    {annovar_dir}/humandb/ \
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

########################################################################
# Create custom sn138 generic db.
# This is because there are some problems with the allele information
# found in the provided hg19_snp138.txt ANNOVAR database.
# Particularly: rs373880503 and rs376103961
########################################################################
rule dbsnp_avinput_to_genericdb:
    input: "/cork/jchung/annovar/dbsnp138/dbsnp138.avinput"
    output: "/cork/jchung/annovar/humandb/hg19_custom_dbsnp138.txt"
    shell: """
        awk -v OFS=$"\t" '{{print $1, $2, $3, $4, $5, $8}}' \
        {input} > {output}
    """

rule dbsnp_vcf_to_avinput:
    input: "/cork/jchung/annovar/dbsnp138/All.vcf"
    output: "/cork/jchung/annovar/dbsnp138/dbsnp138.avinput"
    shell: """
        {annovar_dir}/convert2annovar.pl \
        --format vcf4 {input} \
        --outfile {output} \
        --includeinfo \
    """
        
        
        