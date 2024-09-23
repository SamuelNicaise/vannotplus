from datetime import datetime
import json
import logging as log
import os
from os.path import join as osj

from cyvcf2 import cyvcf2
import numpy as np

from vannotplus.commons import load_ped, run_shell
from vannotplus.family.ped9 import Ped

TEMPLATE = osj(os.path.dirname(__file__), "template.json")


def main_exomiser(input_vcf, output_vcf, app, config):
    with open(TEMPLATE, "r") as f:
        template = json.load(f)
    output_dir = os.path.dirname(output_vcf)  # TODO: change output_dir to tmp dir
    output_dir_in_container = output_dir
    ped = load_ped(config, app)
    vcf_in_container = input_vcf
    for real_path, container_path in config["mount"].items():
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

        cmd = f"-XX:ParallelGCThreads={config['exomiser']['threads']}  -XX:MaxHeapSize={config['exomiser']['heap']}  -jar {config['exomiser']['jar']}"
        cmd += f" --analysis={template_file_in_container}"
        cmd += f" --spring.config.location={config['exomiser']['properties']}"
        cmd += f" --exomiser.data-directory={config['exomiser']['db']}"
        cmd = docker_cmd(config, cmd)
        run_shell(cmd)
        sample_variant_dict[s] = get_annotated_variants(osj(output_dir, s + ".vcf.gz"))

    annots_to_add = [
        "EXOMISER_P_VALUE",
        "EXOMISER_GENE_COMBINED_SCORE",
        "EXOMISER_GENE_PHENO_SCORE",
        "EXOMISER_GENE_VARIANT_SCORE",
        "EXOMISER_VARIANT_SCORE",
    ]
    for annot in annots_to_add:
        vcf.add_format_to_header(
            {
                "ID": annot,
                "Number": 1,
                "Type": "Float",
                "Description": "Computed by vannotplus 1.0.0",
            }
        )
    writer = cyvcf2.Writer(output_vcf, vcf)
    writer.write_header()
    for variant in vcf:
        key = "_".join([variant.CHROM, str(variant.POS), variant.REF, str(variant.ALT)])

        for annot in annots_to_add:
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


def docker_cmd(config, cmd):
    docker_cmd = "docker run --rm --name vannotplus_exomiser"
    for k, v in config["mount"].items():
        docker_cmd += f" -v {k}:{v}"
    docker_cmd += f" --entrypoint /usr/bin/java howard:{config['howard']['version']}"

    cmd = docker_cmd + " " + cmd
    return cmd
