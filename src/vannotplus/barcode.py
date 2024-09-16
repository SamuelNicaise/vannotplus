import logging as log
import os
from os.path import join as osj
import shutil
import subprocess
import tempfile

from cyvcf2 import cyvcf2

from vannotplus.config import APP_TO_PED
from vannotplus.ped9 import Ped, Sample


def run_shell(cmd: str) -> None:
    """
    Run cmd in a separate shell
    Show stdout/stderr only if log level is log.DEBUG
    """
    log.debug(cmd)
    if log.root.level <= 10:
        redirect = None
    else:
        redirect = subprocess.DEVNULL
    subprocess.run(cmd, shell=True, stdout=redirect, stderr=redirect)


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


def main_barcode(input_vcf, output_vcf, ped_dir, app):
    # get samples list from input VCF
    vcf = cyvcf2.VCF(input_vcf)
    # identify ped file through a config file
    ped_path = osj(ped_dir, APP_TO_PED[app])
    # load ped file
    ped = Ped(ped_file=ped_path)
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

        if len(family_samples) > 1:
            print(family_samples)
            cmd = "/home1/data/conda/envs/howard_up_to_date/bin/howard calculation"
            cmd += f" --input {work_vcf}"
            cmd += f" --output {work_vcf}"
            cmd += " --calculations='BARCODEFAMILY'"
            cmd += f" --family_pedigree='{','.join(family_samples)}'"
            print(cmd)
            run_shell(cmd)

            # create howard command

    # merge back VCF for seamless VANNOT integration
    shutil.move(work_vcf, output_vcf)
    tmp_dir.cleanup()
