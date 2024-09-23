import logging as log
import os
from os.path import join as osj
import shutil
import subprocess
import tempfile

from cyvcf2 import cyvcf2

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
