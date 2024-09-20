from datetime import datetime
import json
import os
from os.path import join as osj

from cyvcf2 import cyvcf2
import numpy as np

from vannotplus.config import APP_TO_PED, MOUNT
from vannotplus.family.barcode import run_shell
from vannotplus.family.ped9 import Ped

TEMPLATE = osj(os.path.dirname(__file__), "template.json")


def main_exomiser(input_vcf, output_dir, ped_dir, app):
    with open(TEMPLATE, "r") as f:
        template = json.load(f)
    print(template)

    ped_path = osj(ped_dir, APP_TO_PED[app])
    ped = Ped(ped_file=ped_path)

    vcf_in_container = input_vcf
    output_dir_in_container = output_dir
    for real_path, container_path in MOUNT.items():
        if input_vcf.startswith(real_path):
            vcf_in_container = input_vcf.replace(real_path, container_path)
        if output_dir.startswith(real_path):
            output_dir_in_container = output_dir.replace(real_path, container_path)

    vcf = cyvcf2.VCF(input_vcf)
    assembly = os.path.basename(vcf.get_header_type("reference")["reference"]).split(
        ".fa"
    )[0]

    sample_variant_dict = {}
    for s in vcf.samples:
        print("s", s)

        # write monosample VCF
        sample_vcf = cyvcf2.VCF(input_vcf, samples=s)
        writer = cyvcf2.Writer(osj(output_dir, s + "_work.vcf"), sample_vcf)
        writer.write_header()
        for variant in sample_vcf:
            writer.write_record(variant)
        writer.close()

        write_template(
            template,
            s,
            ped,
            osj(output_dir_in_container, s + "_work.vcf"),
            output_dir_in_container,
            output_dir,
            assembly,
        )
        template_file_in_container = osj(output_dir_in_container, s + "_template.json")

        cmd = "-XX:ParallelGCThreads=12  -XX:MaxHeapSize=75G  -jar /root/howard/tools/exomiser/14.0.0/bin/exomiser-cli-14.0.0.jar"
        cmd += f" --analysis={template_file_in_container}"
        cmd += " --spring.config.location=/root/howard/databases/exomiser/sam/application.properties"
        cmd += " --exomiser.data-directory=/root/howard/databases/exomiser/sam"
        cmd = docker_cmd(cmd)
        print(cmd)
        # run_shell(cmd)
        sample_variant_dict[s] = get_annotated_variants(osj(output_dir, s + ".vcf.gz"))

    annot_to_add = [
        "EXOMISER_P_VALUE",
        "EXOMISER_GENE_COMBINED_SCORE",
        "EXOMISER_GENE_PHENO_SCORE",
        "EXOMISER_GENE_VARIANT_SCORE",
        "EXOMISER_VARIANT_SCORE",
    ]
    for annot in annot_to_add:
        vcf.add_format_to_header(
            {
                "ID": annot,
                "Number": 1,
                "Type": "Float",
                "Description": "Computed by vannotplus 1.0.0",
            }
        )
    writer = cyvcf2.Writer(osj(output_dir, "merged.vcf"), vcf)
    writer.write_header()
    for variant in vcf:
        key = "_".join([variant.CHROM, str(variant.POS), variant.REF, str(variant.ALT)])

        for annot in annot_to_add:
            annot_list = []
            for s in vcf.samples:
                try:
                    annot_list.append(sample_variant_dict[s][key][annot])
                except KeyError:
                    annot_list.append(np.nan)

            variant.set_format(annot, np.array(annot_list, dtype=float))
        writer.write_record(variant)
    writer.close()


def get_annotated_variants(vcf_path: str) -> dict:
    vcf = cyvcf2.VCF(vcf_path)

    # verify description so indexes can be used safely later
    exomiser_header = vcf.get_header_type("Exomiser")
    if (
        "{RANK|ID|GENE_SYMBOL|ENTREZ_GENE_ID|MOI|P-VALUE|EXOMISER_GENE_COMBINED_SCORE|EXOMISER_GENE_PHENO_SCORE|EXOMISER_GENE_VARIANT_SCORE|EXOMISER_VARIANT_SCORE|CONTRIBUTING_VARIANT|WHITELIST_VARIANT|FUNCTIONAL_CLASS|HGVS|EXOMISER_ACMG_CLASSIFICATION|EXOMISER_ACMG_EVIDENCE|EXOMISER_ACMG_DISEASE_ID|EXOMISER_ACMG_DISEASE_NAME}"
        not in exomiser_header["Description"]
    ):
        raise ValueError(
            f"Unexpected Exomiser description in VCF header: {exomiser_header}"
        )

    res = {}
    for variant in vcf:
        key = "_".join([variant.CHROM, str(variant.POS), variant.REF, str(variant.ALT)])
        exomiser_data = variant.INFO["Exomiser"].split("|")
        res[key] = {
            "EXOMISER_P_VALUE": exomiser_data[5],
            "EXOMISER_GENE_COMBINED_SCORE": exomiser_data[6],
            "EXOMISER_GENE_PHENO_SCORE": exomiser_data[7],
            "EXOMISER_GENE_VARIANT_SCORE": exomiser_data[8],
            "EXOMISER_VARIANT_SCORE": exomiser_data[9],
        }
    return res


def write_template(
    template, s, ped, vcf_in_container, output_dir_in_container, output_dir, assembly
):
    template["phenopacket"]["subject"]["id"] = s
    template["phenopacket"]["subject"]["sex"] = ped[s].sex
    template["phenopacket"]["hpoIds"] = ped[s].HPO
    template["phenopacket"]["htsFiles"] = [
        {"uri": vcf_in_container, "htsFormat": "VCF", "genomeAssembly": assembly}
    ]
    # template["phenopacket"]["metaData"]["created"] = datetime.today().strftime(
    #     "%Y-%m-%dT%H:%M:%S%Z"
    # )
    template["outputOptions"]["outputDirectory"] = output_dir_in_container
    template["outputOptions"]["outputFileName"] = s

    template_file = osj(output_dir, s + "_template.json")
    with open(template_file, "w") as f:
        json.dump(template, f)


def docker_cmd(cmd):
    docker_args = [
        "docker",
        "run",
        "--rm",
        "--name",
        "sam_howard",
        "--env",
        f"http_proxy={os.environ['http_proxy']}",
        "--env",
        f"ftp_proxy={os.environ['ftp_proxy']}",
        "--env",
        f"https_proxy={os.environ['https_proxy']}",
        "-v",
        f"/home1/BAS/HOWARD/data:/data",
        "-v",
        f"/home1/BAS/HOWARD/databases:/databases",
        "--entrypoint",
        "/usr/bin/java",
        "howard:0.11.0",
    ]

    cmd = " ".join(docker_args) + " " + cmd
    return cmd


if __name__ == "__main__":
    INPUT_VCF = "/home1/L_PROD/NGS/BAS/HOWARD/data/nicaises/ped/jb_trio_output.vcf"
    # PED_DIR = "/home1/data/WORK_DIR_SAM/Ped_raw/Data"
    PED_DIR = "/home1/L_PROD/NGS/tmp/peds"
    APP = "WES_AGILENT"
    OUTPUT_DIR = "/home1/L_PROD/NGS/BAS/HOWARD/data/nicaises/test"
    # OUTPUT_DIR = "/home1/BAS/HOWARD/data/nicaises/test"

    main_exomiser(INPUT_VCF, OUTPUT_DIR, PED_DIR, APP)
