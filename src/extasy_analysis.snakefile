# Analysis of WGS data using eXtasy
def format_hpo_term(hpo_term):
    formatted_term = hpo_term.replace(":", "_")
    return formatted_term

def list_geneprio_files(hpo_terms, geneprio_dir =  "/home/jchung/programs/eXtasy/geneprios/res/"):
    hpo_terms = [format_hpo_term(t) for t in hpo_terms]
    hpo_files = [term + "_fgs.tsv" for term in hpo_terms]
    hpo_files = [geneprio_dir + file for file in hpo_files]
    hpo_files_string = ','.join(hpo_files)
    return hpo_files_string
    
workdir: "/home/jchung/projects/wgs/"
vcftools = "/apps1/vcftools/vcftools_0.1.11/cpp/vcftools"
SAMPLEID = "BM1452.001 BM1452.100 BM1452.200 BM1453.001 BM1453.100 BM1453.200".split()
bm1452_hpo_terms = ["HP:0001636", "HP:0000389", "HP:0007018", "HP:0000193", "HP:0000135", "HP:0001263", "HP:0001382", "HP:0002967", "HP:0000483", "HP:0000646", "HP:0000540"]
bm1453_hpo_terms = ["HP:0000220", "HP:0000389"]

rule all:
    input: 
        "results/extasy_output/BM1452/BM1452.001.PhasedSnp.recode.combined.extasy results/extasy_output/BM1453/BM1453.001.PhasedSnp.recode.combined.extasy".split()
   
rule run_extasy_bm1453:
    input:
        "data/filtered/BM1453.001.PhasedSnp.recode.vcf"
    params:
        geneprios = list_geneprio_files(bm1453_hpo_terms),
        output = "results/extasy_output/BM1453/BM1453.001.PhasedSnp.recode"
    output:
        "results/extasy_output/BM1453/BM1453.001.PhasedSnp.recode.combined.extasy"
    shell: """
    extasy.rb \
    -i {input} \
    -g {params.geneprios} \
    -p {params.output} \
    -c 
    """

rule run_extasy_bm1452:
    input: 
        "data/filtered/BM1452.001.PhasedSnp.recode.vcf"
    params: 
        geneprios = list_geneprio_files(bm1452_hpo_terms),
        output = "results/extasy_output/BM1452/BM1452.001.PhasedSnp.recode"
    output:
        "results/extasy_output/BM1452/BM1452.001.PhasedSnp.recode.combined.extasy"
    shell: """
    extasy.rb \
    -i {input} \
    -g {params.geneprios} \
    -p {params.output} \
    -c 
    """
    
rule check_recode:
    input: expand("data/filtered/{sampleId}.PhasedSnp.recode.vcf", sampleId = SAMPLEID)

rule export_individual_sample:
    input: "data/wgs_phased_snp.vcf"
    params: 
        sample = "{sampleId}",
        outPrefix = "data/filtered/{sampleId}.PhasedSnp"
    output: "data/filtered/{sampleId}.PhasedSnp.recode.vcf"
    shell: """
    {vcftools} \
    --vcf {input} \
    --indv {params.sample} \
    --remove-filtered-all \
    --non-ref-ac 1 \
    --max-non-ref-ac 2 \
    --out {params.outPrefix} \
    --recode
    """
