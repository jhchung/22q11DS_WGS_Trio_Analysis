# Analysis of whole genome sequencing data using VAAST
import os

workdir: "/home/jchung/projects/wgs/"
vcftools = "/apps1/vcftools/vcftools_0.1.11/cpp/vcftools"
vaast_sort_gff = "/home/jchung/programs/VAAST_Code_2.1.0/bin/vaast_tools/vaast_sort_gff"
vaast_converter = "/home/jchung/programs/VAAST_Code_2.1.0/bin/vaast_tools/vaast_converter"
quality_check = "/home/jchung/programs/VAAST_Code_2.1.0/bin/vaast_tools/quality-check.pl"

refseqFeatures = "/home/jchung/programs/VAAST_Code_2.1.0/data/hg19/Features/refGene_hg19.gff3"
refseqFeaturesSort = "/home/jchung/programs/VAAST_Code_2.1.0/data/hg19/Features/refGene_hg19.sorted.gff3"
hg19Fasta = "/home/jchung/programs/VAAST_Code_2.1.0/data/hg19/Fasta/vaast_hsap_chrs_hg19.fa"
background1KG = "/home/jchung/programs/VAAST_Code_2.1.0/data/hg19/Background_CDR/1KG_refGene_Dec2011_CGDiv_NHLBI_NoCall.cdr"

BM1452Target = ["BM1452.001"]
BM1453Target = ["BM1453.001"]
BM1452Background = "BM1452.100 BM1452.200".split()
BM1453Background = "BM1453.100 BM1453.200".split()
BM1452 = BM1452Target + BM1452Background
BM1453 = BM1453Target + BM1453Background
sampleIDs = BM1452 + BM1453

# subworkflow wesBackgroundWorkflow:
    # workdir: "/home/jchung/projects/wgs/"
    # snakefile: "/home/jchung/projects/wgs/src/vaastPrepareWESControlForBackground.snakefile"
    
rule all:
    input: "results/vaastOutput/VAAST/BM1452AnalysisRecCompWES.vaast",
           "results/vaastOutput/VAAST/BM1452AnalysisRecCompWESTopFeatures.vaast",
           "results/vaastOutput/VAAST/BM1452AnalysisRecComp1KG.vaast",
           "results/vaastOutput/VAAST/BM1452AnalysisRecComp1KGTopFeatures.vaast"
           
rule rerunVAASTAnalysisOnTopFeaturesWES:
    input: target = "results/vaastOutput/VST/BM1452Target.cdr",
           trio = "results/vaastOutput/VST/BM1452Background.union.cdr",
           background = "results/vaastOutput/wesBackground/VST/wesControls.cdr",
           features = refseqFeaturesSort,
           sigFeatures = "results/vaastOutput/VAAST/BM1452AnalysisRecCompWESSigFeatures.txt"
    params: outputPrefix = "results/vaastOutput/VAAST/BM1452AnalysisRecCompWESTopFeatures"
    output: "results/vaastOutput/VAAST/BM1452AnalysisRecCompWESTopFeatures.vaast"
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
        --significance 2.4e-6 \
        --features {sigFeatures} \
        --use_aas_info y \
        --codon_bias \
        --fast_gp \
        --mp2 10 \
        --outfile {params.outputPrefix} \
        --trio {input.trio} \
        -k \
        {input.features} \
        {input.background} \
        {input.target}""")
    
rule extractSignificantFeaturesWES:
    input: vaastResult = "results/vaastOutput/VAAST/BM1452AnalysisRecCompWES.vaast"
    output: sigFeatureFile = "results/vaastOutput/VAAST/BM1452AnalysisRecCompWESSigFeatures.txt"
    run: 
        featureNames = []
        genomicPvalues = []
        with open(input.vaastResult, "r") as f:
            for a_line in f:
                if a_line.startswith(">"):
                    featureNames.append(a_line)
                if a_line.startswith("genome_permutation_p"):
                    genomicPvalues.append(a_line)
        
        # Find lines that have P-value less than cutoff
        genomicPvalues = [pval.replace("genome_permutation_p:", "") for pval in genomicPvalues]
        genomicPvalues = [float(pval) for pval in genomicPvalues]
        
        significantPvalues = []
        for pval in genomicPvalues:
            if pval <= 0.1:
                significantPvalues.append(True)
            else:
                significantPvalues.append(False)
        
        significantPvaluesIndex = [i for i, x in enumerate(significantPvalues) if x == True]
        
        significantFeatures = []
        for sigIndex in significantPvaluesIndex:
            significantFeatures.append(featureNames[sigIndex])
            
        formattedSignificantFeatures = []
        for feature in significantFeatures:
            formattedFeature = feature.split()[0]
            formattedFeature = formattedFeature.replace(">", "")
            formattedSignificantFeatures.append(formattedFeature)
            
        outputString = ",".join(formattedSignificantFeatures)
        with open(output.sigFeatureFile, "w") as o:
            o.write(outputString)
    


rule runVAASTOnBM1452RecessiveCompleteWES:
    input: target = "results/vaastOutput/VST/BM1452Target.cdr",
           trio = "results/vaastOutput/VST/BM1452Background.union.cdr",
           background = "results/vaastOutput/wesBackground/VST/wesControls.cdr",
           features = refseqFeaturesSort
    params: outputPrefix = "results/vaastOutput/VAAST/BM1452AnalysisRecCompWES"
    output: "results/vaastOutput/VAAST/BM1452AnalysisRecCompWES.vaast"
    shell: """
    VAAST \
    --mode lrt \
    --inheritance r \
    --penetrance c \
    --locus_heterogeneity y \
    -gp 1e5 \
    --use_aas_info y \
    --codon_bias \
    --fast_gp \
    --mp1 10 \
    --outfile {params.outputPrefix} \
    --trio {input.trio} \
    -k \
    {input.features} \
    {input.background} \
    {input.target}"""
    
rule rerunVAASTAnalysisOnTopFeatures1KG:
    input: target = "results/vaastOutput/VST/BM1452Target.cdr",
           trio = "results/vaastOutput/VST/BM1452Background.union.cdr",
           background = background1KG,
           features = refseqFeaturesSort,
           sigFeatures = "results/vaastOutput/VAAST/BM1452AnalysisRecComp1KGSigFeatures.txt"
    params: outputPrefix = "results/vaastOutput/VAAST/BM1452AnalysisRecComp1KGTopFeatures"
    output: "results/vaastOutput/VAAST/BM1452AnalysisRecComp1KGTopFeatures.vaast"
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
        --significance 2.4e-6 \
        --features {sigFeatures} \
        --use_aas_info y \
        --codon_bias \
        --fast_gp \
        --mp2 10 \
        --outfile {params.outputPrefix} \
        --trio {input.trio} \
        -k \
        {input.features} \
        {input.background} \
        {input.target}""")
    
rule extractSignificantFeatures1KG:
    input: vaastResult = "results/vaastOutput/VAAST/BM1452AnalysisRecComp1KG.vaast"
    output: sigFeatureFile = "results/vaastOutput/VAAST/BM1452AnalysisRecComp1KGSigFeatures.txt"
    run: 
        featureNames = []
        genomicPvalues = []
        with open(input.vaastResult, "r") as f:
            for a_line in f:
                if a_line.startswith(">"):
                    featureNames.append(a_line)
                if a_line.startswith("genome_permutation_p"):
                    genomicPvalues.append(a_line)
        
        # Find lines that have P-value less than cutoff
        genomicPvalues = [pval.replace("genome_permutation_p:", "") for pval in genomicPvalues]
        genomicPvalues = [float(pval) for pval in genomicPvalues]
        
        significantPvalues = []
        for pval in genomicPvalues:
            if pval <= 0.1:
                significantPvalues.append(True)
            else:
                significantPvalues.append(False)
        
        significantPvaluesIndex = [i for i, x in enumerate(significantPvalues) if x == True]
        
        significantFeatures = []
        for sigIndex in significantPvaluesIndex:
            significantFeatures.append(featureNames[sigIndex])
            
        formattedSignificantFeatures = []
        for feature in significantFeatures:
            formattedFeature = feature.split()[0]
            formattedFeature = formattedFeature.replace(">", "")
            formattedSignificantFeatures.append(formattedFeature)
            
        outputString = ",".join(formattedSignificantFeatures)
        with open(output.sigFeatureFile, "w") as o:
            o.write(outputString)
rule runVAASTOnBM1452RecessiveComplete1KG:
    input: target = "results/vaastOutput/VST/BM1452Target.cdr",
           trio = "results/vaastOutput/VST/BM1452Background.union.cdr",
           background = background1KG,
           features = refseqFeaturesSort
    params: outputPrefix = "results/vaastOutput/VAAST/BM1452AnalysisRecComp1KG"
    output: "results/vaastOutput/VAAST/BM1452AnalysisRecComp1KG.vaast"
    shell: """
    VAAST \
    --mode lrt \
    --inheritance r \
    --penetrance c \
    --locus_heterogeneity y \
    -gp 1e05 \
    --fast_gp \
    --use_aas_info y \
    --codon_bias \
    --mp1 10 \
    --outfile {params.outputPrefix} \
    --trio {input.trio} \
    -k \
    {input.features} \
    {input.background} \
    {input.target}"""

# Run VAAST on BM1452 family. Assume that the mutations show a recessive 
# inheritance pattern and are incompletly penetrant.
rule runVAASTOnBM1452RecessiveComplete:
    input: target = "results/vaastOutput/VST/BM1452Target.cdr",
           trio = "results/vaastOutput/VST/BM1452Background.union.cdr",
           background = background1KG,
           features = refseqFeaturesSort
    params: outputPrefix = "results/vaastOutput/VAAST/BM1452Analysis"
    output: "results/vaastOutput/VAAST/BM1452Analysis.vaast"
    shell: """
    VAAST \
    --mode lrt \
    --inheritance r \
    --penetrance c \
    --locus_heterogeneity y \
    -gp 1e05 \
    --fast_gp \
    --mp1 4 \
    --outfile {params.outputPrefix} \
    --trio {input.trio} \
    {input.features} \
    {input.background} \
    {input.target}"""

rule BM1452QualityCheckBackground:
    input: "results/vaastOutput/VST/BM1452Target.cdr",
           "results/vaastOutput/wesBackground/VST/wesControls.cdr"
    output: "results/vaastOutput/VST/BM1452WesBackgroundQC.txt"
    shell: """
    echo {input[0]}
    echo {input[1]}
    {quality_check} \
    -bootstrap 1000 \
    {input[0]} \
    {input[1]} > {output}"""
    
rule runVSTOnBM1453Background:
    input: expand("results/vaastOutput/VAT/refGene.{BM1453BACKGROUND}.vat.gvf", BM1453BACKGROUND = BM1453Background)
    output: "results/vaastOutput/VST/BM1453Background.union.cdr"
    shell: """
    VST \
    -o 'U(0..1)' \
    -chunk 100000000 \
    -b hg19 \
    {input} \
    > {output} """

rule runVSTOnBM1453Target:
    input: expand("results/vaastOutput/VAT/refGene.{BM1453TARGET}.vat.gvf", BM1453TARGET = BM1453Target)
    output: "results/vaastOutput/VST/BM1453Target.cdr"
    shell: """
    VST \
    -o 'I(0)' \
    -chunk 100000000 \
    -b hg19 \
    {input} \
    > {output} """

rule runVSTOnBM1452Background:
    input: expand("results/vaastOutput/VAT/refGene.{BM1452BACKGROUND}.vat.gvf", BM1452BACKGROUND = BM1452Background)
    output: "results/vaastOutput/VST/BM1452Background.union.cdr"
    shell: """
    VST \
    -o 'U(0..1)' \
    -chunk 100000000 \
    -b hg19 \
    {input} \
    > {output} """

rule runVSTOnBM1452Target:
    input: expand("results/vaastOutput/VAT/refGene.{BM1452TARGET}.vat.gvf", BM1452TARGET = BM1452Target)
    output: "results/vaastOutput/VST/BM1452Target.cdr"
    shell: """
    VST \
    -o 'I(0)' \
    -chunk 100000000 \
    -b hg19 \
    {input} \
    > {output} """
    
rule annotateVariants:
    input: expand("results/vaastOutput/VAT/refGene.{SAMPLEIDS}.vat.gvf", SAMPLEIDS = sampleIDs)
    
rule runVAT:
    input: "results/vaastOutput/gvfFiles/{SAMPLEIDS}.gvf",
           refseqFeaturesSort,
           hg19Fasta
    output: "results/vaastOutput/VAT/refGene.{SAMPLEIDS}.vat.gvf"
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
    input: "data/forVAAST/wgs_phased_snp_pass_autosomal_genic.recode.vcf"
    params: outputDir = "results/vaastOutput/gvfFiles"
    output: expand("results/vaastOutput/gvfFiles/{SAMPLEIDS}.gvf", SAMPLEIDS = sampleIDs)
    shell: """
    {vaast_converter} \
    --build hg19 \
    --path {params.outputDir} \
    {input}"""

rule filterWgsFile:
    input: "data/wgs_phased_snp.vcf",
           "data/forVAAST/refGene_hg19.bed"
    params: outputPrefix = "data/forVAAST/wgs_phased_snp_pass_autosomal_genic"
    output: "data/forVAAST/wgs_phased_snp_pass_autosomal_genic.recode.vcf"
    shell: """
    {vcftools} \
    --remove-filtered-all \
    --vcf {input[0]} \
    --not-chr X \
    --not-chr Y \
    --mac 1 \
    --max-mac 12 \
    --minQ 30 \
    --bed {input[1]} \
    --out {params.outputPrefix} \
    --recode \
    --recode-INFO-all """

rule createRefGeneBedFromGff:
    input: refseqFeatures
    output: "data/forVAAST/refGene_hg19.bed"
    shell: """
        Rscript src/convertGffToBed.R {input} {output}
    """

