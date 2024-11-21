from datetime import datetime
import glob
import os
import os.path as op
from os.path import join as osj
import subprocess
import tempfile

from cyvcf2 import cyvcf2
import pyBigWig


def get_output_vcf(output_dir, alfa_raw_file):
    REGION_DICT = {
        "ALFA_AFA": "ALFA African American",
        "ALFA_AFO": "ALFA African Others",
        "ALFA_AFR": "ALFA African",
        "ALFA_ASN": "ALFA Asian",
        "ALFA_EAS": "ALFA East Asian",
        "ALFA_EUR": "ALFA European",
        "ALFA_GLB": "ALFA total population",
        "ALFA_LAC": "ALFA Latin American 1",
        "ALFA_LEN": "ALFA Latin American 2",
        "ALFA_OAS": "ALFA Other Asian",
        "ALFA_OTR": "ALFA Other population",
        "ALFA_SAS": "ALFA South Asian",
    }

    world_region = op.basename(alfa_raw_file).split(".bb")[0]
    output_path = osj(output_dir, world_region + ".vcf")

    vcf_header = [
        "##fileformat=VCFv4.4",
        "##fileDate=" + datetime.today().strftime("%m/%d/%Y"),
        "##InputFile=" + alfa_raw_file,
        "\t".join(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"]),
    ]

    with tempfile.NamedTemporaryFile(mode="wt", delete=False) as tmp:
        for l in vcf_header:
            tmp.write(l + "\n")

    vcf = cyvcf2.Writer(output_path, cyvcf2.VCF(tmp.name))
    os.remove(tmp.name)

    vcf.add_info_to_header(
        {
            "ID": world_region,
            "Number": 1,
            "Type": "Float",
            "Description": f"Allele Frequency Aggregator, comming from chip array from almost 1M subject in dbGaP project, {REGION_DICT[world_region]}",
        }
    )
    vcf.close()

    return output_path


class Variant:
    def __init__(self, chrom, pos, rs, ref, alt, info):
        self.chr = chrom
        self.pos = pos
        self.id = rs
        self.ref = ref
        self.alt = alt
        self.info = info

    def __str__(self):
        return "\t".join(
            [self.chr, self.pos, self.id, self.ref, self.alt, ".", ".", self.info]
        )


def bigbedentry_to_variant(entry, world_region) -> list[Variant]:
    """
    Goal: convert entry such as
    (7004498, 7004498, 'rs2151341419\t208\t.\t7004498\t7004498\t156,101,151\tchr17:7004498,REF_AF(T)=0.7917;ALT_AF(C)=0.2083')
    into a VCF formatted variant
    """
    # only keep the last tuple and split it
    entry = entry[2].split("\t")
    rs = entry[0]
    if not (rs.startswith("rs") or rs == "."):
        raise ValueError(f"no rs found - bad bigbed entry: {entry}")

    # only keep the last part and split it
    entry = entry[-1].split(",", maxsplit=1)
    # entry should now be like ["chr17:7004498", "REF_AF(T)=0.7917;ALT_AF(C)=0.2083"]
    chrom = entry[0].split(":")[0]
    pos = entry[0].split(":")[1]

    entry = entry[1].split(";")
    # entry should now be like ["REF_AF(T)=0.7917", "ALT_AF(C)=0.2083"]
    ref = entry[0].split("REF_AF(")[1].split(")")[0]
    alts = entry[1].split("ALT_AF(")[1].split(")")[0]
    data = entry[1].split("=")[1]

    # this func returns a list containing one variant per existing alt
    res = []

    alts_list = alts.split(",")
    if len(alts_list) == 1:
        res = [Variant(chrom, pos, rs, ref, alts, f"{world_region}={data}")]
    elif len(alts_list) > 1:
        data_list = data.split(",")
        for i in range(len(alts_list)):
            variant = Variant(
                chrom, pos, rs, ref, alts_list[i], f"{world_region}={data_list[i]}"
            )
            res.append(variant)

    return res


def main():
    OUTPUT_DIR = "/home1/DB/HOWARD/ALFA/hg19/sam"
    BGZIP = "/home1/HUB/bin/htslib/1.14/bin/bgzip"
    ALFA_RAW_FILES = glob.glob("/home1/DB/HOWARD/ALFA/hg19/rawdata/*.bb")

    for raw in ALFA_RAW_FILES:
        world_region = op.basename(raw).split(".bb")[0]
        print("input:", raw)
        output_vcf = get_output_vcf(OUTPUT_DIR, raw)
        print("output:", output_vcf)

        with open(output_vcf, "a") as out:
            bb = pyBigWig.open(raw)
            for chrom, length in bb.chroms().items():
                for entry in bb.entries(chrom, 1, length):
                    # print(entry)
                    if not entry[-1].split("\t")[-1].endswith("REF_AF=0;ALT_AF=0"):
                        entries_as_variants = bigbedentry_to_variant(
                            entry, world_region
                        )
                        for var in entries_as_variants:
                            out.write(str(var) + "\n")

        subprocess.run(f"{BGZIP} {output_vcf}", shell=True, stdout=None, stderr=None)


if __name__ == "__main__":
    main()
