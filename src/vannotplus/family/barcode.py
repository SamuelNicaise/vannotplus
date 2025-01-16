import logging as log
import os
from os.path import join as osj
import shutil
import tempfile

from cyvcf2 import cyvcf2
import numpy as np

from vannotplus.commons import load_ped, run_shell
from vannotplus.family.ped9 import Ped, Sample


def get_parental_aliases(sample: Sample, ped: Ped, is_mother: bool) -> list[str]:
    """
    For a given sample
        if is_mother = True:
            return maternal ID + its aliases
        else
            return paternal ID + its aliases
    """
    aliases = []
    if is_mother:
        if sample.maternal_id not in ("", None):
            aliases = [sample.maternal_id]
            if ped.get(sample.maternal_id, "") not in ("", None):
                aliases += ped[sample.maternal_id].alias
    else:
        if sample.paternal_id not in ("", None):
            aliases = [sample.paternal_id]
            if ped.get(sample.paternal_id, "") not in ("", None):
                aliases += ped[sample.paternal_id].alias
    return aliases


def main_barcode(input_vcf, output_vcf, app, config):
    vcf = cyvcf2.VCF(input_vcf)
    ped = load_ped(config, app)
    tmp_dir = tempfile.TemporaryDirectory()
    work_vcf = osj(tmp_dir.name, os.path.basename(input_vcf))
    shutil.copy2(input_vcf, work_vcf)

    # for each index,identify if one or more of its parents are in the VCF
    # if yes, run a howard familybarcode command on the VCF
    for s in vcf.samples:
        family_samples = [s]
        if s in ped:
            maternal_aliases = get_parental_aliases(ped[s], ped, True)
            for a in maternal_aliases:
                if a in vcf.samples:
                    family_samples.append(a)
            paternal_aliases = get_parental_aliases(ped[s], ped, False)
            for a in paternal_aliases:
                if a in vcf.samples:
                    family_samples.append(a)

        # howard command can be run iteratively with the same vcf in input and output, removing the need for a final merge
        if len(family_samples) > 1:
            log.info(f"Computing family: {family_samples}")
            cmd = config["howard"]["bin"]
            cmd += f" --input {work_vcf}"
            cmd += f" --output {work_vcf}"
            cmd += " --calculations='BARCODEFAMILY'"
            cmd += f" --family_pedigree='{','.join(family_samples)}'"
            log.debug(cmd)
            run_shell(cmd)

    shutil.move(work_vcf, output_vcf)
    tmp_dir.cleanup()


def get_sample_to_family_dict(input_vcf: cyvcf2.VCF, ped: Ped) -> dict[str, list[str]]:
    """
    for each sample in the VCF, identify if one or more of its parents are in the VCF
    for example
        if ped says:sample A has father M and mother F
        and vcf contains A and F
    this func returns:
    {
        A: [F],
        F: []
    }
    """
    sample_family_dict: dict[str, list[str]] = {}
    for s in input_vcf.samples:
        family_samples = [s]
        if s in ped:
            maternal_aliases = get_parental_aliases(ped[s], ped, True)
            for a in maternal_aliases:
                if a in input_vcf.samples:
                    family_samples.append(a)
            paternal_aliases = get_parental_aliases(ped[s], ped, False)
            for a in paternal_aliases:
                if a in input_vcf.samples:
                    family_samples.append(a)

        sample_family_dict[s] = family_samples

    return sample_family_dict


def get_families_indexes(input_vcf: cyvcf2.VCF, ped: Ped) -> list[list[int]]:
    """
    for each sample in the VCF, identify if one or more of its parents are in the VCF
    for example
        if ped says:
            sample A1 has father M1 and mother F1
            sample A2 has father M2 and mother F2
        and vcf contains (in order) the samples: A1, M1, F1, B, F2, A2
    this func returns:
    [[021], [1], [2], [3], [4], [54]]
    indexes are ordered and the mother is always before the father if both are available
    """
    families_indexes: list[list[int]] = []
    for i, s in enumerate(input_vcf.samples):
        family_samples = [i]
        if s in ped:
            maternal_aliases = get_parental_aliases(ped[s], ped, True)
            for a in maternal_aliases:
                if a in input_vcf.samples:
                    family_samples.append(input_vcf.samples.index(a))
            paternal_aliases = get_parental_aliases(ped[s], ped, False)
            for a in paternal_aliases:
                if a in input_vcf.samples:
                    family_samples.append(input_vcf.samples.index(a))

        families_indexes.append(family_samples)

    return families_indexes


def main_barcode_fast(
    input_vcf_path: str, output_vcf_path: str, app: str, config: dict
):
    # use a copied temporary vcf as input so its header can be modified for ease of development
    tmp_dir = tempfile.TemporaryDirectory()
    work_vcf_path = osj(tmp_dir.name, os.path.basename(input_vcf_path))
    shutil.copy2(input_vcf_path, work_vcf_path)
    # note: gts012=True is extremely important whenever using cyvcf2.Variant.gt_types, see cyvcf2 doc
    input_vcf = cyvcf2.VCF(work_vcf_path, gts012=True)

    ped = load_ped(config, app)

    # for each sample, get indexes corresponding to its parents in the input VCF if they exist
    families_indexes = get_families_indexes(input_vcf, ped)

    # change input header as variants will originate from input_vcf
    input_vcf.add_format_to_header(
        {
            "ID": "BCF",
            "Number": 1,
            "Type": "Integer",
            "Description": "Family barcode: for each sample in the family, assign 1 integer depending on genotype. 0 = wild type or unknown, 1 = heterozygous, 2 = homozygous. The family's sample list can be found in the BCFS tag.",
        }
    )
    input_vcf.add_format_to_header(
        {
            "ID": "BCFS",
            "Number": ".",
            "Type": "String",
            "Description": "Samples in the family barcode",
        }
    )
    output_vcf = cyvcf2.Writer(output_vcf_path, input_vcf)

    # then iterate over input_vcf and write variants with barcodes in output_vcf
    for var in input_vcf:
        if var.POS % 1000 == 0:
            print(var.CHROM, var.POS)
        bcf_list = []
        bcfs_list = []

        genotypes: np.ndarray = var.gt_types
        # consider that unknown (3) are wild type (0) in barcode
        genotypes[var.gt_types == 3] = 0

        for i, _ in enumerate(input_vcf.samples):
            if len(families_indexes[i]) == 1:
                # no family
                bcf_list.append(".")
                bcfs_list.append(".")
            else:
                barcode = "".join([str(v) for v in genotypes[families_indexes[i]]])
                barcode_samples = ",".join(
                    list(map(input_vcf.samples.__getitem__, families_indexes[i]))
                )
                bcf_list.append(barcode)
                bcfs_list.append(barcode_samples)

        var.set_format("BCF", np.asarray(bcf_list, dtype=np.bytes_))
        var.set_format("BCFS", np.asarray(bcfs_list, dtype=np.bytes_))
        output_vcf.write_record(var)

    output_vcf.close()
