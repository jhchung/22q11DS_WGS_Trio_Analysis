# Analysis of whole genome sequencing data using VAAST
import os

def proband_id(family_id):
    return family_id + ".001"

def parent_id(family_id):
    return [family_id + s for s in [".100", ".200"]]

def combine_family(proband, parents):
    family = parents
    family.append(proband)
    return family

def filter_vaast_results(vaast_result_file, output_file, pvalue_cutoff):
    featureNames = []
    genomicPvalues = []
    with open(vaast_result_file, "r") as f:
        for a_line in f:
            if a_line.startswith(">"):
                featureNames.append(a_line)
            if a_line.startswith("genome_permutation_p"):
                genomicPvalues.append(a_line)
    
    # Find lines that have P-value less than cutoff
    genomicPvalues = extract_pvalue(genomicPvalues)
    significantPvalues = check_pvalue(genomicPvalues, pvalue_cutoff)
    significantFeatures = extract_significant_features(featureNames, significantPvalues)
    formattedSignificantFeatures = format_feature_name(significantFeatures)
        
    outputString = ",".join(formattedSignificantFeatures)
    with open(output_file, "w") as o:
        o.write(outputString)

def extract_pvalue(pvalue_strings):
    formatted_pvalue = [pval.replace("genome_permutation_p:", "") for pval in pvalue_strings]
    formatted_pvalue = [float(pval) for pval in formatted_pvalue]
    return formatted_pvalue

def check_pvalue(pvalues, pvalue_cutoff):
    significant_pvalues = []
    for pval in pvalues:
        if pval <= pvalue_cutoff:
            significant_pvalues.append(True)
        else:
            significant_pvalues.append(False)
    return significant_pvalues

def extract_significant_features(feature_names, significant_pvalues):
    significant_pvalues_index = [i for i, x in enumerate(significant_pvalues) if x == True]
    
    significant_features = []
    for sigIndex in significant_pvalues_index:
        significant_features.append(feature_names[sigIndex])
    return significant_features

def format_feature_name(feature_names):
    formatted_significant_features = []
    for feature in feature_names:
        formatted_feature = feature.split()[0]
        formatted_feature = formatted_feature.replace(">", "")
        formatted_significant_features.append(formatted_feature)
    return formatted_significant_features

workdir: "/home/jchung/projects/wgs/"
vcftools = "/apps1/vcftools/vcftools_0.1.11/cpp/vcftools"
vaast_sort_gff = "/home/jchung/programs/VAAST_Code_2.0.1/bin/vaast_tools/vaast_sort_gff"
vaast_converter = "/home/jchung/programs/VAAST_Code_2.0.1/bin/vaast_tools/vaast_converter"
quality_check = "/home/jchung/programs/VAAST_Code_2.0.1/bin/vaast_tools/quality-check.pl"

refseqFeatures = "/home/jchung/programs/VAAST_Code_2.0.1/data/hg19/Features/refGene_hg19.gff3"
refseqFeaturesSort = "/home/jchung/programs/VAAST_Code_2.0.1/data/hg19/Features/refGene_hg19.sorted.gff3"
hg19Fasta = "/home/jchung/programs/VAAST_Code_2.0.1/data/hg19/Fasta/vaast_hsap_chrs_hg19.fa"
variantMask = "/home/jchung/programs/VAAST_Code_2.0.1/data/hg19/Variant_Masking/40bp_pe.hg19.bed"
background1KG = "/home/jchung/programs/VAAST_Code_2.0.1/data/hg19/Background_CDR/1KG_refGene_Dec2011_CGDiv_NHLBI_NoCall.cdr"

famk1_id = "BM1452"
famk2_id = "BM1453"
family_ids = [famk1_id, famk2_id]

probandk1 = proband_id(famk1_id)
probandk2 = proband_id(famk2_id)

parentsk1 = parent_id(famk1_id)
parentsk2 = parent_id(famk2_id)

familyk1 = combine_family(probandk1, parentsk1)
familyk2 = combine_family(probandk2, parentsk2)

all_samples = familyk1 + familyk2
    
number_of_cores = 8

rule all:
    input: 
        expand("results/vaastOutput/VAAST/{FAM_ID}.AnalysisRecComp1KG.TopFeatures.vaast", FAM_ID = family_ids)

rule check_vaast_dominant_model:
    input: 
        expand("results/vaastOutput/VAAST/{FAM_ID}.AnalysisDomInc1KG.TopFeatures.vaast", FAM_ID = family_ids)

rule top_genes_vaast_dominant_model:
    input: 
        target = "results/vaastOutput/VST/{FAM_ID}.proband.cdr",
        background = background1KG,
        features = refseqFeaturesSort,
        sigFeatures = "results/vaastOutput/VAAST/{FAM_ID}.AnalysisDomInc1KG.sigFeatures.txt",
        variantMaskFile = variantMask
    params: 
        outputPrefix = "results/vaastOutput/VAAST/{FAM_ID}.AnalysisDomInc1KG.TopFeatures"
    output: 
        "results/vaastOutput/VAAST/{FAM_ID}.AnalysisDomInc1KG.TopFeatures.vaast"
    run: 
        with open(input.sigFeatures, "r") as sigFeaturesFile:
            sigFeatures = sigFeaturesFile.read()
        print(sigFeatures)
        shell("""
            VAAST \
            --mode lrt \
            -x {input.variantMaskFile} \
            --inheritance d \
            --penetrance i \
            --locus_heterogeneity y \
            -gp 1e10 \
            --features {sigFeatures} \
            --use_aas_info y \
            --codon_bias \
            --fast_gp \
            --mp2 {number_of_cores} \
            --outfile {params.outputPrefix} \
            -k \
            {input.features} \
            {input.background} \
            {input.target}
        """)
        
rule extract_significant_features_dominant_model:
    input: 
        "results/vaastOutput/VAAST/{FAM_ID}.AnalysisDomInc1KG.vaast"
    output: 
        "results/vaastOutput/VAAST/{FAM_ID}.AnalysisDomInc1KG.sigFeatures.txt"
    run:
        filter_vaast_results(input[0], output[0], 0.1)
        
rule vaast_dominant_model:
    input: 
        target = "results/vaastOutput/VST/{FAM_ID}.proband.cdr",
        background = background1KG,
        features = refseqFeaturesSort,
        variantMaskFile = variantMask
    params: 
        outputPrefix = "results/vaastOutput/VAAST/{FAM_ID}.AnalysisDomInc1KG"
    output: 
        "results/vaastOutput/VAAST/{FAM_ID}.AnalysisDomInc1KG.vaast"
    shell: """
    VAAST \
    --mode lrt \
    -x {input.variantMaskFile} \
    --inheritance d \
    --penetrance i \
    --locus_heterogeneity y \
    -gp 1e05 \
    --splice_site \
    --fast_gp \
    --use_aas_info y \
    --codon_bias \
    --mp1 {number_of_cores} \
    --outfile {params.outputPrefix} \
    -k \
    {input.features} \
    {input.background} \
    {input.target}"""
    
rule top_genes_vaast_recessive_model:
    input: 
        target = "results/vaastOutput/VST/{FAM_ID}.proband.cdr",
        trio = "results/vaastOutput/VST/{FAM_ID}.parents.union.cdr",
        background = background1KG,
        features = refseqFeaturesSort,
        sigFeatures = "results/vaastOutput/VAAST/{FAM_ID}.AnalysisRecComp1KG.SigFeatures.txt"
    params: 
        outputPrefix = "results/vaastOutput/VAAST/{FAM_ID}.AnalysisRecComp1KG.TopFeatures"
    output:
        "results/vaastOutput/VAAST/{FAM_ID}.AnalysisRecComp1KG.TopFeatures.vaast"
    run: 
        with open(input.sigFeatures, "r") as sigFeaturesFile:
            sigFeatures = sigFeaturesFile.read()
        print(sigFeatures)
        shell("""
        VAAST \
        --mode lrt \
        --inheritance r \
        --penetrance c \
        --locus_heterogeneity y \
        -gp 1e10 \
        --features {sigFeatures} \
        --use_aas_info y \
        --codon_bias \
        --fast_gp \
        --mp2 {number_of_cores} \
        --outfile {params.outputPrefix} \
        --trio {input.trio} \
        -k \
        {input.features} \
        {input.background} \
        {input.target}""")
    
rule extract_significant_features_recessive_model:
    input: 
        vaastResult = "results/vaastOutput/VAAST/{FAM_ID}.AnalysisRecComp1KG.vaast"
    output: 
        "results/vaastOutput/VAAST/{FAM_ID}.AnalysisRecComp1KG.SigFeatures.txt"
    run: 
        filter_vaast_results(input.vaastResult, output[0], 0.1)

# Run VAAST on BM1453 family. Assume that the mutations show a recessive 
# inheritance pattern and are incompletly penetrant.
rule vaast_recessive_model:
    input: 
        target = "results/vaastOutput/VST/{FAM_ID}.proband.cdr",
        trio = "results/vaastOutput/VST/{FAM_ID}.parents.union.cdr",
        background = background1KG,
        features = refseqFeaturesSort
    params: 
        outputPrefix = "results/vaastOutput/VAAST/{FAM_ID}.AnalysisRecComp1KG"
    output: 
        "results/vaastOutput/VAAST/{FAM_ID}.AnalysisRecComp1KG.vaast"
    shell: """
    VAAST \
    --mode lrt \
    --inheritance r \
    --penetrance c \
    --locus_heterogeneity y \
    -gp 1e05 \
    --fast_gp \
    --splice_site \
    --use_aas_info y \
    --codon_bias \
    --mp1 {number_of_cores} \
    --outfile {params.outputPrefix} \
    --trio {input.trio} \
    -k \
    {input.features} \
    {input.background} \
    {input.target}"""

rule run_VST_on_Parents:
    input: 
        "results/vaastOutput/VAT/refGene.{FAM_ID}.100.vat.gvf",
        "results/vaastOutput/VAT/refGene.{FAM_ID}.200.vat.gvf"
    output: 
        "results/vaastOutput/VST/{FAM_ID}.parents.union.cdr"
    shell: """
    VST \
    -o 'U(0..1)' \
    -chunk 100000000 \
    -b hg19 \
    {input[0]} {input[1]} \
    > {output} """

rule vst_on_proband:
    input: "results/vaastOutput/VAT/refGene.{FAM_ID}.001.vat.gvf"
    output: "results/vaastOutput/VST/{FAM_ID}.proband.cdr"
    shell: """
    VST \
    -o 'I(0)' \
    -chunk 100000000 \
    -b hg19 \
    {input} \
    > {output} """
   
rule check_vat_for_all_samples:
    input: 
        expand("results/vaastOutput/VAT/refGene.{ALL_SAMPLES}.vat.gvf", ALL_SAMPLES = all_samples)
    
rule run_vat:
    input: 
        "results/vaastOutput/gvfFiles/{ALL_SAMPLES}.gvf",
        refseqFeaturesSort,
        hg19Fasta
    output: 
        "results/vaastOutput/VAT/refGene.{ALL_SAMPLES}.vat.gvf"
    shell: """
    VAT \
    --features {input[1]} \
    --build hg19 \
    -c 100000000 \
    --fasta {input[2]} \
    {input[0]} \
    > {output}"""

rule convert_vcf_2_gvf:
    input: "data/forVAAST/wgs_phased_snp_pass_autosomal_genic.recode.vcf"
    params: outputDir = "results/vaastOutput/gvfFiles"
    output: expand("results/vaastOutput/gvfFiles/{ALL_SAMPLES}.gvf", ALL_SAMPLES = all_samples)
    shell: """
    {vaast_converter} \
    --build hg19 \
    --path {params.outputDir} \
    {input}"""
    
rule filter_wgs_variants:
    input: "data/wgs_phased_snp.vcf"
    params: outputPrefix = "data/forVAAST/wgs_phased_snp_pass_autosomal_genic"
    output: "data/forVAAST/wgs_phased_snp_pass_autosomal_genic.recode.vcf"
    shell: """
    {vcftools} \
    --remove-filtered-all \
    --vcf {input[0]} \
    --mac 1 \
    --max-mac 12 \
    --minQ 30 \
    --out {params.outputPrefix} \
    --recode \
    --recode-INFO-all """

rule sort_refseq_gff:
    input: refseqFeatures
    output: refseqFeaturesSort
    shell: """
    {vaast_sort_gff} \
    {input} > {output}"""

rule convert_refgene_bed_to_gff:
    input: refseqFeatures
    output: "data/forVAAST/refGene_hg19.bed"
    shell: """
        Rscript src/convertGffToBed.R {input} {output}
    """
