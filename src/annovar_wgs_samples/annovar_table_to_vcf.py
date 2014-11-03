import sys

def extract_vcf(input, output):
    """Extract VCF file from ANNOVAR table annotations.
    
    Args:
        input: Path to input file.
        output: Path to output file.
    """
    in_file = open(input, "r")
    out_file = open(output, "w")
    
    for line in in_file:
        if line.startswith("Chr\tStart"):
            annovar_header = line.split("\t")
            vcf_start_column = len(annovar_header) - 1
        elif line.startswith("#"):
            vcf_header = line.split()
            output_vcf_header = " ".join(vcf_header)
            
            out_file.write(output_vcf_header + "\n")
        elif not line.startswith("#"):
            data = line.split("\t")
            vcf_columns = data[vcf_start_column: ]
            vcf_columns = "\t".join(vcf_columns)
            out_file.write(vcf_columns)
    
    in_file.close()
    out_file.close()
    
if __name__ == '__main__':
    annovar_file = sys.argv[1]
    vcf_file = sys.argv[2]
    extract_vcf(annovar_file, vcf_file)