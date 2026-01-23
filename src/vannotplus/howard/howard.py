"""
0) load config
1) split the input file by chr
2) for each chr.vcf:
    for each line:
        get the variant id
        get the annotation from the db
        add the annotation to the INFO field, renaming it if required in config
        add the line to the output file
is there anything at this point that has to be computed? There could be if there were scores based on sample-specific annotations, but there isn't for now.
3) merge all chr.vcf into the output vcf

no need to call exomiser and barcode, vannot will do it with the output vcf
"""

def main_howard():
    pass

if __name__ == "__main__":
    main_howard()