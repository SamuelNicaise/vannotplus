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
    DEPRECATED

    for each sample in the VCF, identify if one or more of its parents are in the VCF
    for example
        if ped says:
            sample A1 has father M1 and mother F1
            sample A2 has father M2 and mother F2
        and vcf contains (in order) the samples: A1, M1, F1, B, F2, A2
    this func returns:
    [[0, 2, 1], [1], [2], [3], [4], [5, 4]]
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


def get_families_indexes_v2(input_vcf: cyvcf2.VCF, ped: Ped) -> list[list[int]]:
    """
    For each sample in the VCF, identify the indexes (VCF columns) of other samples in the family.
    Returns them in an ordered list of lists.

    Some implementation details :
    - If a sample has no family, its list will only contain itself
    - All samples of a given family have the same barcode, so the same indexes.
    - Index order is : all affected samples of the family, then parents (mother first), then their parents if any (mother first) (then their parents, etc), then remaining samples of the family.
    - The above implies that trios will have the usual barcode: index-mother-father
    - It is possible for mother/father samples to have a different family ID, in case of pooled parental samples. In this case the pooled parental samples will have no barcode. There should be a APP#POOL STARK tag to help detect this.

    TODO: does not actually work yet with grandparents. Only parents are looked up.

    for example
        if ped says:
            sample A1 has father M1 and mother F1
            sample A2 has father M2 (not in VCF) and mother F2
            sample A3 has father M3 and mother F3 and other members of the family X3 (affected) and Y3 (unaffected)
            ped does not contain sample B
        and vcf contains (in order) the samples: A1, M1, F1, B, F2, A2, A3, M3, F3, X3, Y3
                                         index:  0   1    2  3  4   5   6   7   8   9   10
    this func returns:
    [[0,2,1], [0,2,1], [0,2,1], [3], [5,4], [5,4], [6,9,7,8,10], [6,9,7,8,10], [6,9,7,8,10], [6,9,7,8,10], [6,9,7,8,10]]
    """
    families_indexes: list[list[int]] = []

    for s in input_vcf.samples:
        family = ped.get_family_from_sample(s)
        if family == None:
            indexes = samples_to_indexes(input_vcf.samples, [s])
        else:
            barcode_samples = get_samples_for_barcode(ped, family)
            indexes = samples_to_indexes(input_vcf.samples, barcode_samples)
        families_indexes.append(indexes)

    return families_indexes


def get_samples_for_barcode(ped: Ped, family: str) -> list[str]:
    """
    Get an ordered list of sample names. Order corresponds to specifications in get_families_indexes_v2 docstring.
    """
    samples = ped.get_samples_from_family(family)
    if len(samples) == 0:
        return []
    res = []

    # 1) Find ordered parents/grandparents/etc
    unique_parents = []
    for s in samples:
        parents = s.get_parents()
        for p in parents:
            if p not in unique_parents and p not in ("", None):
                unique_parents.append(p)

    # 2) Add affected samples except if they're going to be in the parents block
    for s in samples:
        if s.is_affected and s.individual_id not in unique_parents:
            res.append(s.individual_id)

    # 3) Add the parents
    res.extend(unique_parents)

    # 4) Add any remaining samples from the family
    for s in samples:
        if s.individual_id not in res:
            res.append(s.individual_id)

    # This process should not result in duplicates
    if len(res) != len(set(res)):
        raise ValueError(f"Duplicate samples detected in the result list: {res}")

    return res


# def recursive_get_parents(
#     ped: Ped, sample: Sample, generation: int = 0, parents: list[str] = []
# ) -> list[str]:
#     """
#     Recursively get all parents/grandparents/greatgrandparents/etc of a sample.
#     Parents are ordered by generation, with the closest parents first and the mother first in each pair of parents.

#     TODO: function works in standalone but not when repeatedly called from get_families_indexes_v2. Need to find the hidden side effect.
#     """
#     if sample.maternal_id not in ("", None):
#         if sample.maternal_id not in parents:
#             parents.append(sample.maternal_id)
#         recursive_get_parents(ped, ped[sample.maternal_id], generation + 1, parents)
#     if sample.paternal_id not in ("", None):
#         if sample.paternal_id not in parents:
#             parents.append(sample.paternal_id)
#         recursive_get_parents(ped, ped[sample.paternal_id], generation + 1, parents)
#     return parents


def samples_to_indexes(sample_list: list[str], samples_to_find: list[str]) -> list[int]:
    indexes = []
    for s in samples_to_find:
        if s in sample_list:
            indexes.append(sample_list.index(s))
    return indexes


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
    families_indexes = get_families_indexes_v2(input_vcf, ped)

    # change input header as variants will originate from input_vcf
    input_vcf.add_format_to_header(
        {
            "ID": "BCF",
            "Number": 1,
            "Type": "String",
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
        if var.POS % 100000 == 0:
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


if __name__ == "__main__":
    VCF = cyvcf2.VCF("/home1/HUB/bin/vannotplus/vannotplus/mort_subite.vcf")
    PED = Ped("/home1/L_PROD/NGS/PRODUCTION/ped_raw/MORT_SUBITE.json")
    # name = "DRFG73"
    # fam = ped.get_family_from_sample(name)
    # if fam:
    #     samples = ped.get_samples_from_family(fam)
    #     print([str(s) for s in samples])
    #     print(get_samples_for_barcode(ped, fam))
    #     print(samples_to_indexes(vcf.samples, get_samples_for_barcode(ped, fam)))
    # print(res)

    # l = vcf.samples
    # print(l.index("DRFG98"))
    # print(vcf.samples)
    # print("DRFG98" in l)
    # print(l.index("DRFG98"))
    # print(samples_to_indexes(l, ["DRFG98", "DRFG99"]))

    # indexes = get_families_indexes_v2(vcf, ped)
    # indexes_str = ";".join([",".join(map(str, sublist)) for sublist in indexes])
    # print(indexes_str)
    # print(indexes)

    # from vannotplus.commons import load_config

    # main_barcode_fast(
    #     "/home1/L_PROD/NGS/STARK18/data/VCF_MORT_SUBITE/adding_barcodes/res/MORT_SUBITE_merged_no_barcodes.vcf.gz",
    #     "/home1/L_PROD/NGS/STARK18/data/VCF_MORT_SUBITE/adding_barcodes/res/MORT_SUBITE_merged_with_barcodes.vcf.gz",
    #     "MORT_SUBITE",
    #     load_config("/home1/HUB/bin/vannotplus/vannotplus/src/vannotplus/config.yml"),
    # )
