import json
import logging as log
import os
from os.path import join as osj
import shutil
import tempfile
import random

from cyvcf2 import cyvcf2
import numpy as np

from vannotplus import __version__
from vannotplus.commons import get_variant_id, load_ped, run_shell
from vannotplus.family.ped9 import Ped

TEMPLATE = osj(os.path.dirname(__file__), "template.json")


def any_sample_has_HPOs(samples: list[str], ped: Ped) -> bool:
    for s in samples:
        if sample_has_HPOs(s, ped):
            return True
    return False


def sample_has_HPOs(sample: str, ped: Ped) -> bool:
    if sample in ped:
        if ped[sample].HPO not in (None, [], [""]):
            return True
    return False


def main_exomiser(input_vcf, output_vcf, app, config, remove_info_in_tmp=True):
    """
    Split input_vcf into one VCF per sample
    Write the JSON template with a phenopacket using the sample's HPOs (pedigree determined by app)
    Run Exomiser in a docker container for each sample
    Merge the Exomiser annotations back into the original VCF, as a sample-specific FORMAT field

    If remove_info_in_tmp is True, the INFO field is removed from the temporary monosample VCFs
    to avoid issues with Exomiser being limited to VCF version <= 4.2
    All of the original INFO fields are retained in the final output VCF no matter what.
    The only reason to use remove_info_in_tmp=False is improving performance, if the input VCF is already compatible with Exomiser.
    """
    with open(TEMPLATE, "r") as f:
        template = json.load(f)
    output_dir = os.path.dirname(output_vcf)
    tmp_dir = tempfile.mkdtemp(dir=output_dir)
    tmp_dir_in_container = tmp_dir
    ped = load_ped(config, app)
    # vcf_in_container = input_vcf
    for real_path, container_path in config["mount"].items():
        # if input_vcf.startswith(real_path):
        #   vcf_in_container = input_vcf.replace(real_path, container_path)
        if tmp_dir.startswith(real_path):
            tmp_dir_in_container = tmp_dir.replace(real_path, container_path)
    log.debug(f"tmp_dir: {tmp_dir}")
    log.debug(f"tmp_dir_in_container: {tmp_dir_in_container}")

    vcf = cyvcf2.VCF(input_vcf)
    assembly = os.path.basename(vcf.get_header_type("reference")["reference"]).split(
        ".fa"
    )[0]

    # make sure input and output dirs are mounted in docker no matter what the config is
    io_dirs = [os.path.dirname(input_vcf), os.path.dirname(output_vcf)]
    for io_dir in io_dirs:
        if io_dir == "":
            io_dir = os.path.abspath(".")
        if io_dir not in config["mount"]:
            config["mount"][io_dir] = io_dir
    log.debug(f"config[mount]:{config['mount']}")

    if not any_sample_has_HPOs(vcf.samples, ped):
        log.debug(
            f"No HPO found for samples in {input_vcf}, copying it to {output_vcf} without change"
        )
        vcf.close()
        shutil.copy(input_vcf, output_vcf)
        return

    sample_variant_dict = {}
    for s in vcf.samples:
        if not sample_has_HPOs(s, ped):
            log.info(f"No HPO found for sample {s}, skipping Exomiser for this sample")
            shutil.copy(input_vcf, osj(tmp_dir, s + ".vcf.gz"))
            sample_variant_dict[s] = get_annotated_variants(
                osj(tmp_dir, s + ".vcf.gz"), no_hpos=True
            )
            continue

        # write monosample VCF
        sample_vcf = cyvcf2.VCF(input_vcf, samples=s)
        writer = cyvcf2.Writer(osj(tmp_dir, s + "_exomiserinput.vcf"), sample_vcf)
        writer.write_header()
        for variant in sample_vcf:
            if not remove_info_in_tmp:
                writer.write_record(variant)
            else:
                # Exomiser 14.0.0 does not follow the 4.4 VCF spec allowing spaces in INFO fields
                l = str(variant).strip().split("\t")
                l[7] = "."
                infoless_variant = writer.variant_from_string("\t".join(l))
                writer.write_record(infoless_variant)
        writer.close()

        write_template(
            template,
            s,
            ped,
            osj(tmp_dir_in_container, s + "_exomiserinput.vcf"),
            tmp_dir_in_container,
            tmp_dir,
            assembly,
        )
        template_file_in_container = osj(tmp_dir_in_container, s + "_template.json")

        cmd = f"-XX:ParallelGCThreads={config['exomiser']['threads']}  -XX:MaxHeapSize={config['exomiser']['heap']}  -jar {config['exomiser']['jar']}"
        cmd += f" --analysis={template_file_in_container}"
        cmd += f" --spring.config.location={config['exomiser']['properties']}"
        cmd += f" --exomiser.data-directory={config['exomiser']['db']}"
        cmd = docker_cmd(config, cmd)
        run_shell(cmd)
        sample_variant_dict[s] = get_annotated_variants(osj(tmp_dir, s + ".vcf.gz"))

    try:
        annots_to_add = config["exomiser"]["annotations_to_add"]
    except KeyError:
        # Default annotation: only phenotype score
        annots_to_add = ["EXOMISER_GENE_PHENO_SCORE"]
    descriptions = {k: f"Exported from vannotplus {__version__}" for k in annots_to_add}
    descriptions["EXOMISER_GENE_PHENO_SCORE"] = (
        "Gene-specific phenotype relevance score computed by Exomiser. 1 means the gene is highly linked to the HPOs, 0 means no link. Based on semantic similarity of the patient's HPO terms to phenotypic annotations of genes in 1) OMIM, with inheritance consistency checked and adjusted if inconsistent (the score is halved if inheritance is incompatible) and 2) phenotypes of protein-protein associated neighboring genes derived from protein interaction networks (hiphive)."
    )

    for annot in annots_to_add:
        vcf.add_format_to_header(
            {
                "ID": annot,
                "Number": 1,
                "Type": "Float",
                "Description": descriptions[annot],
            }
        )
    writer = cyvcf2.Writer(output_vcf, vcf)
    writer.write_header()
    for variant in vcf:
        key = get_variant_id(variant)

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

    if log.root.level > 10:  # if log level > debug
        shutil.rmtree(tmp_dir)


def get_annotated_variants(vcf_path: str, no_hpos=False) -> dict:
    """
    Return a dict of dicts with Exomiser annotations for each variant in the VCF, such as:
    {
        "chr1_123456_A_T": {
            "EXOMISER_P_VALUE": "0.001",
            "EXOMISER_GENE_COMBINED_SCORE": "0.85",
            "EXOMISER_GENE_PHENO_SCORE": "0.9",
            "EXOMISER_GENE_VARIANT_SCORE": "0.8",
            "EXOMISER_VARIANT_SCORE": "0.75"
        },
        ...
    }
    If no_hpos is True, return an empty dict for each variant.
    """
    log.debug(f"get_annotated_variants::vcf_path:{vcf_path}")
    vcf = cyvcf2.VCF(vcf_path)

    if not no_hpos:
        # verify description so indexes can be used safely later
        exomiser_header = vcf.get_header_type("Exomiser")
        if (
            "{RANK|ID|GENE_SYMBOL|ENTREZ_GENE_ID|MOI|P-VALUE|EXOMISER_GENE_COMBINED_SCORE|EXOMISER_GENE_PHENO_SCORE|EXOMISER_GENE_VARIANT_SCORE|EXOMISER_VARIANT_SCORE|CONTRIBUTING_VARIANT|WHITELIST_VARIANT|FUNCTIONAL_CLASS|HGVS|EXOMISER_ACMG_CLASSIFICATION|EXOMISER_ACMG_EVIDENCE|EXOMISER_ACMG_DISEASE_ID|EXOMISER_ACMG_DISEASE_NAME}"
            not in exomiser_header["Description"]
        ):
            raise ValueError(
                f"Unexpected Exomiser description in VCF header: {exomiser_header} -- failing VCF: {vcf_path}"
            )

    res = {}
    for variant in vcf:
        key = get_variant_id(variant)
        if no_hpos:
            res[key] = {}
        else:
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
    template, s, ped, vcf_in_container, tmp_dir_in_container, tmp_dir, assembly
):
    template["analysis"]["proband"] = s
    if s in ped:
        template["analysis"]["hpoIds"] = ped[s].HPO
    else:
        template["analysis"]["hpoIds"] = []
    template["analysis"]["vcf"] = vcf_in_container
    template["analysis"]["genomeAssembly"] = assembly

    template["outputOptions"]["outputDirectory"] = tmp_dir_in_container
    template["outputOptions"]["outputFileName"] = s

    template_file = osj(tmp_dir, s + "_template.json")
    with open(template_file, "w") as f:
        json.dump(template, f)


def docker_cmd(config, cmd):
    random_tag = random.randint(1, 1000000)
    docker_name = f"VANNOTPLUS_exomiser_{random_tag}"
    docker_cmd = f"docker run --rm --name {docker_name}"
    for k, v in config["mount"].items():
        docker_cmd += f" -v {k}:{v}"
    docker_cmd += f" --entrypoint /usr/bin/java howard:{config['howard']['version']}"

    cmd = docker_cmd + " " + cmd
    return cmd
