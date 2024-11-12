"""
This is a proof of concept for bigwig file annotation
Not used by the rest of the code for now
"""

import cyvcf2
import pyBigWig


def main():
    VCF_PATH = "/home1/L_PROD/NGS/BAS/HOWARD/data/nicaises/configs/KLA2403985.naked.vcf"
    OUTPUT_PATH = "/home1/L_PROD/NGS/BAS/HOWARD/data/nicaises/KLA.fastphast.vcf"
    BW_PATHS = {
        "phastCons100ways": "/home1/DB/HOWARD/phastCons/100way/hg19/hg19.100way.phastCons.bw",
        "phyloP100way": "/home1/DB/HOWARD/phyloP/100way/hg19/hg19.100way.phyloP100way.bw",
        "GERP": "/home1/DB/HOWARD/GERP/All_hg19_RS.bw",
    }
    bw_loaded = {}
    input_vcf = cyvcf2.VCF(VCF_PATH)

    for bw_name, bw_path in BW_PATHS.items():
        input_vcf.add_info_to_header(
            {
                "ID": bw_name,
                "Number": 1,
                "Type": "Float",
                "Description": f"{bw_name} by vannotplus",
            }
        )
        bw_loaded[bw_name] = pyBigWig.open(bw_path)

    output_vcf = cyvcf2.Writer(OUTPUT_PATH, input_vcf)

    for variant in input_vcf:
        # bigwig are 0-based https://github.com/deeptools/pyBigWig?tab=readme-ov-file#a-note-on-coordinates
        for bw_name, bw_db in bw_loaded.items():
            res = bw_db.values(variant.CHROM, variant.POS - 1, variant.POS)
            variant.INFO[bw_name] = res[0]
        # phastcons_res = bw_phastcons.values(variant.CHROM, variant.POS - 1, variant.POS)
        # variant.INFO["phastCons100ways"] = phastcons_res[0]
        # phylop_res = bw_phylop.values(variant.CHROM, variant.POS - 1, variant.POS)
        # variant.INFO["phyloP100way"] = phylop_res[0]
        output_vcf.write_record(variant)

    output_vcf.close()
    print("done")


if __name__ == "__main__":
    main()
